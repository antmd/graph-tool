## Process this file with automake to produce Makefile.in

SUBDIRS = src

graphtooldocdir = ${PACKAGE_DOC_DIR}
nobase_dist_graphtooldoc_DATA = \
	README\
	COPYING\
	AUTHORS\
	INSTALL\
	NEWS\
	TODO\
	src/boost-workaround/LICENSE_1_0.txt\
	doc/centrality.rst \
	doc/correlations.rst \
	doc/price.py \
	doc/search_module.rst \
	doc/util.rst \
	doc/clustering.rst \
	doc/draw.rst \
	doc/index.rst \
	doc/spectral.rst \
	doc/community.rst  \
	doc/flow.rst \
	doc/Makefile \
	doc/quickstart.rst \
	doc/generation.rst \
	doc/mathjax.py \
	doc/run_action.rst \
	doc/stats.rst \
	doc/graph_tool.rst \
	doc/modules.rst \
	doc/topology.rst \
	doc/conf.py \
	doc/pyenv.py \
	doc/community.xml \
	doc/search_example.xml \
	doc/sphinxext/README.txt \
	doc/sphinxext/LICENSE.txt \
	doc/sphinxext/__init__.py \
	doc/sphinxext/setup.py \
	doc/sphinxext/comment_eater.py \
	doc/sphinxext/docscrape.py \
	doc/sphinxext/MANIFEST.in \
	doc/sphinxext/phantom_import.py \
	doc/sphinxext/compiler_unparse.py \
	doc/sphinxext/docscrape_sphinx.py \
	doc/sphinxext/numpydoc.py \
	doc/sphinxext/traitsdoc.py \
	doc/sphinxext/tests/test_docscrape.py \
	doc/gt_theme/static/flasky.css \
	doc/gt_theme/static/gt_style.css \
	doc/gt_theme/theme.conf \
	doc/.static/graph-icon.png \
	doc/.static/graph-icon.ico


EXTRA_DIST = $(graphtooldoc_DATA) autogen.sh

# Copy all the spec files. Of cource, only one is actually used.
dist-hook:
	for specfile in *.spec; do \
		if test -f $$specfile; then \
			cp -p $$specfile $(distdir); \
		fi \
	done