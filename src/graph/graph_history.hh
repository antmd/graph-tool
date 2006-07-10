// graph-tool -- a general graph modification and manipulation thingy
//
// Copyright (C) 2006  Tiago de Paula Peixoto <tiago@forked.de>
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

#include "../graph/node_graph.hh"
#include <string>
#include <iostream>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/line.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/exception.hpp>
#include <boost/regex.hpp>
#include <boost/algorithm/string.hpp>
#include <ext/hash_set>

#include <unistd.h>    /* standard unix functions, like getpid()         */
#include <sys/types.h> /* various type definitions, like pid_t           */
#include <signal.h>    /* signal name macros, and the signal() prototype */

#include <config.h>
#include "../sim/histogram.hh"

using namespace std;
using namespace earthquake;
using namespace boost;

class HistoryAction
{
public:
    HistoryAction(string file_name, bool compressed)
    :_file_name(file_name), _stream_initialized(false), _compressed(compressed)
    {
        // get the previous histories from file (if it exists already)
        // and recompress file if necessary
        ifstream input_file;
        iostreams::filtering_stream<boost::iostreams::input> input_stream;
        bool new_file = true;
        
        if (filesystem::exists(_file_name))
        {
            new_file = false;
            _temp_file = _file_name+"-graph-tool-temp";
            if (_compressed)
            {
                try
                {
                    filesystem::rename(_file_name, _temp_file);
                }
                catch (filesystem::filesystem_error &e)
                {
                    throw NodeGraphException("error renaming file " + _file_name + " to " + _temp_file + ": " + e.what());
                }
            }
            
            try 
            {    
                // old file or new file
                if (_compressed)
                    input_file.open(_temp_file.c_str(), ios_base::in | ios_base::binary);
                else
                    input_file.open(_file_name.c_str(), ios_base::in | ios_base::binary);
                if (_compressed)
                    input_stream.push(iostreams::gzip_decompressor());
                input_stream.push(input_file);
                input_stream.exceptions(ios_base::badbit);
            }
            catch (ios_base::failure &e)
            {
                throw NodeGraphException("error opening existing file " + _file_name + ": " + e.what());
            }
        }
        
        try 
        {                
            _file_stream.rdbuf()->pubsetbuf(0,0); // kill buffering
            _file_stream.open(file_name.c_str(), ios_base::out | ios_base::app | ios_base::binary);
            _file_stream.exceptions(ios_base::badbit | ios_base::failbit);
            _file_stream.setf(ios_base::unitbuf); // always flush output
            if (_compressed)
            {
                _stream.push(iostreams::gzip_compressor());
                _stream.push(_file_stream);
            }
            else
                _stream.push(_file_stream,0,0);
            _stream.exceptions(ios_base::badbit | ios_base::failbit);
            _stream.precision(10);
            _stream.setf(ios_base::scientific);
        }
        catch (ios_base::failure &e)
        {
            throw NodeGraphException("error opening file " + _file_name + ": " + e.what());
        }
        
        if (!new_file)
        {
            try
            {
                string line;
                while(input_stream)
                {
                    getline(input_stream, line);
                    if (line == "")
                        continue;
                    stringstream line_stream(line);
                    size_t index = numeric_limits<size_t>::max();
                    line_stream >> index;
                    if (index == numeric_limits<size_t>::max())
                        throw NodeGraphException("error reading existing file " + _file_name + ": invalid line \"" + line + "\"");
                    _previous_histories.insert(index-1);
                    if (_compressed)
                        _stream << line << endl;
                }
            }
            catch (ios_base::failure &e)
            {
                if (_compressed)
                    throw NodeGraphException("error recompressing existing file " + _file_name + ": " + e.what());
                else
                    throw NodeGraphException("error reading existing file " + _file_name + ": " + e.what());
            }
        }
    }
    virtual void operator()(ostream &s) = 0;
    virtual string FileName() {return _file_name;}
    virtual ostream& Stream() {return _stream;}
    virtual bool IsPreviousHistory(size_t index) {return _previous_histories.find(index) != _previous_histories.end();}
    virtual void SetPreviousHistory(size_t index) {_previous_histories.insert(index);}
    virtual ~HistoryAction()
    {
        _stream.pop();
        if (_compressed && filesystem::exists(_temp_file))
        {
            try 
            {
                filesystem::remove(_temp_file);
            }
            catch (filesystem::filesystem_error &e)
            {
                cerr << "Warning: error deleting file " + _temp_file + ": " + e.what();
            }
        }
    }
private:
    string _file_name;
    ofstream _file_stream;
    iostreams::filtering_stream<boost::iostreams::output> _stream;
    bool _stream_initialized;
    bool _compressed; 
    __gnu_cxx::hash_set<size_t> _previous_histories;
    string _temp_file;
};


class HistoryLineFilter : public iostreams::line_filter 
{
public:
    explicit HistoryLineFilter(size_t iter): _iter(iter) { }
private:
    string do_filter(const std::string& line)
    {
        stringstream s;
        s << _iter << "\t" << line;
        return s.str();
    }
    size_t _iter;
};

//==============================================================================
// Signal handling stuff
//==============================================================================
static bool _interrupt_signal = false;
void catch_sig_stop(int sig_num);


//==============================================================================
// WriteHistoryElements
//==============================================================================

void WriteHistoryElements(size_t history_index, list<HistoryAction*> &hist_actions)
{
    iostreams::filtering_stream<boost::iostreams::output> stream;
    stream.push(HistoryLineFilter(history_index+1));
    for (typeof(hist_actions.begin()) iter = hist_actions.begin(); iter != hist_actions.end(); ++iter)
    {
        if (_interrupt_signal)
            break;
        if ((*iter)->IsPreviousHistory(history_index))
            continue;
        else
            (*iter)->SetPreviousHistory(history_index);
        stream.push((*iter)->Stream());
        try 
        {
            (**iter)(stream);
            stream.strict_sync();
        }
        catch (std::ios_base::failure &e)
        {
            throw NodeGraphException("error writing to file " + (*iter)->FileName() + ": " + e.what());
        }
        stream.pop();
    }
}

//==============================================================================
// WriteHistory (vector)
//==============================================================================

void WriteHistory(list<HistoryAction*> &hist_actions, vector<size_t> hist_freq, pair<size_t,size_t> event_range, size_t max_edges, NodeGraph &ngraph)
{
    // setup signal handling
    _interrupt_signal = false;
    signal(SIGINT, catch_sig_stop);
    signal(SIGTERM, catch_sig_stop);
    signal(SIGQUIT, catch_sig_stop);
    signal(SIGHUP, catch_sig_stop);
    signal(SIGPWR, catch_sig_stop);
    
    pair<size_t,size_t> index_range;
    for(typeof(hist_freq.begin()) hist_iter = hist_freq.begin(); hist_iter != hist_freq.end(); ++hist_iter)
    {
        if (_interrupt_signal)
            break;
        if (*hist_iter < event_range.first || *hist_iter > event_range.second)
            continue;
        if (*hist_iter - 1 - event_range.first >= max_edges)
            index_range.first = *hist_iter - max_edges;
        else
            index_range.first = event_range.first;
        if (*hist_iter > 0)
            index_range.second = *hist_iter - 1;       
        else
            index_range.second = *hist_iter;
        ngraph.SetIndexRangeFilter(index_range);
	if (ngraph.GetNumberOfEdges() < index_range.second - index_range.first)
	    index_range.second = index_range.first + ngraph.GetNumberOfEdges() - 1;
        WriteHistoryElements(index_range.second, hist_actions);
    }
}

//==============================================================================
// WriteHistory (range)
//==============================================================================

void WriteHistory(list<HistoryAction*> &hist_actions, pair<size_t,size_t> hist_range, size_t hist_freq, pair<size_t,size_t> event_range, size_t max_edges, NodeGraph &ngraph)
{
    // setup signal handling
    _interrupt_signal = false;
    signal(SIGINT, catch_sig_stop);
    signal(SIGTERM, catch_sig_stop);
    signal(SIGQUIT, catch_sig_stop);
    signal(SIGHUP, catch_sig_stop);
    signal(SIGPWR, catch_sig_stop);
    
    pair<size_t,size_t> index_range;
    index_range.first = event_range.first;
    hist_range.second = min(hist_range.second, min(event_range.second, ngraph.GetNumberOfEdges() + hist_range.first));
    for(index_range.second = hist_range.first + hist_freq - 1; index_range.second < hist_range.second; index_range.second += hist_freq)
    {
        if (_interrupt_signal)
            break;
        if (index_range.second - index_range.first >= max_edges)
            index_range.first = index_range.second - max_edges + 1;
    
        ngraph.SetIndexRangeFilter(index_range);
        WriteHistoryElements(index_range.second, hist_actions);
    }
}

//==============================================================================
// WriteHistoryAvgInDeg (vector)
//==============================================================================
void WriteHistoryAvgInDeg(list<HistoryAction*> &hist_actions, vector<double> hist_freq, pair<size_t,size_t> event_range, NodeGraph &ngraph)
{
    // setup signal handling
    _interrupt_signal = false;
    signal(SIGINT, catch_sig_stop);
    signal(SIGTERM, catch_sig_stop);
    signal(SIGQUIT, catch_sig_stop);
    signal(SIGHUP, catch_sig_stop);
    signal(SIGPWR, catch_sig_stop);
    
    for(typeof(hist_freq.begin()) hist_iter = hist_freq.begin(); hist_iter != hist_freq.end(); ++hist_iter)
    {
        if (_interrupt_signal)
            break;
        ngraph.SetIndexRangeFilter(event_range);
        double avg_in_deg = ngraph.SetIndexRangeFilter(*hist_iter);
        if (avg_in_deg != 0.0)
            WriteHistoryElements(ngraph.GetIndexRangeFilter().second, hist_actions);
    }
}



template <class T>
class MemberWrap : public HistoryAction
{
public:
    MemberWrap(T member, string file_name, NodeGraph *ngraph, bool compressed=false)
        : HistoryAction(file_name, compressed), _member(member), _ngraph(ngraph) {}
    virtual ~MemberWrap(){}
    virtual void operator()(ostream &s) 
    { 
        typename T::result_type res = _member(_ngraph);
        s << res;
    }
    
private:
    T _member;
    NodeGraph *_ngraph;
};

template <class T>
MemberWrap<T> * make_member_wrap(T member, string file_name, NodeGraph *ngraph, bool compressed=false)
{
    return new MemberWrap<T>(member, file_name, ngraph, compressed);
}

template <class T>
class HistogramMeanWrap : public HistoryAction
{
public:
    HistogramMeanWrap(T member, string file_name, NodeGraph *ngraph, bool compressed = false)
        :HistoryAction(file_name, compressed), _member(member), _ngraph(ngraph){}
    virtual ~HistogramMeanWrap(){}
    
    virtual void operator()(ostream &s)
    {
        double mean = GetHistogramMean(_member(_ngraph));
        s <<  mean;
    }
    
private:
    T _member;
    NodeGraph *_ngraph;
};

template <class T>
HistogramMeanWrap<T> * 
make_histogram_mean_wrap(T member, string file_name, NodeGraph *ngraph, bool compressed=false)
{
    return new HistogramMeanWrap<T>(member, file_name, ngraph, compressed);
}


template <class T> 
class HistogramWrap : public HistoryAction
{
public:
    HistogramWrap (T member, string file_name, NodeGraph *ngraph, bool compressed = true, bool normalize = true) 
        : HistoryAction(file_name, compressed), _member(member), _ngraph(ngraph), _normalize(normalize){}
    virtual ~HistogramWrap(){}

    virtual void operator()(ostream &s)
    {
        double total = 0.0;
        typename T::result_type m = _member(_ngraph);
        if (_normalize)
            for (typeof(m.begin()) iter = m.begin(); iter != m.end(); iter++)
                total += iter->second;
        s.precision(10);
        s.setf(ios_base::scientific);
        for (typeof(m.begin()) iter = m.begin(); iter != m.end(); iter++)
        {
            s << iter->first << '\t';
            if (_normalize)
                s << double(iter->second)/total << std::endl;
            else
                s << double(iter->second) << std::endl;
        }
    }
private:
    T _member;
    NodeGraph *_ngraph;
    bool _normalize;
};

template <class T>
HistogramWrap<T> * make_histogram_wrap(T member, string file_name, NodeGraph *ngraph, bool compressed=true, bool normalize = true)
{
    return new HistogramWrap<T>(member, file_name, ngraph, compressed, normalize);
}

//==============================================================================
// Signal handling stuff
//==============================================================================
extern char *_argvz;
void catch_sig_stop(int sig_num)
{
    std::cerr << _argvz << ": received signal ";
    switch (sig_num)
    {
    case SIGINT:
        std::cerr << "SIGINT (Interrupt).";
        break;
    case SIGTERM:
        std::cerr << "SIGTERM (Termination).";
        break;
    case SIGQUIT:
        std::cerr << "SIGQUIT (Terminal quit).";
        break;
    case SIGHUP:
        std::cerr << "SIGHUP (Hangup).";
        break;
    case SIGPWR:
        std::cerr << "SIGPWR (Power failure restart).";
        break;
    }
    std::cerr << " Stopping the program and saving files. May take a long time." << std::endl;
    _interrupt_signal = true;
}
