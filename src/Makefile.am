## Process this file with automake to produce Makefile.in

AUTOMAKE_OPTIONS = subdir-objects

SUBDIRS = graph .

bin_SCRIPTS = graph-tool

EXTRA_DIST = $(bin_SCRIPTS)

# Also install graph-tool as a python module (graph_tool.py)
install-exec-hook:
	ln -s $(bindir)/graph-tool $(pythondir)/graph_tool.py
uninstall-hook:
	rm -f $(pythondir)/graph_tool.py

