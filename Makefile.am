## Process this file with automake to produce Makefile.in

SUBDIRS = src

graphtooldocdir = ${PACKAGE_DOC_DIR}
graphtooldoc_DATA = \
	README\
	COPYING\
	AUTHORS\
	INSTALL\
	NEWS\
	TODO\
	src/boost-workaround/LICENSE_1_0.txt

EXTRA_DIST = $(graphtooldoc_DATA)

# Copy all the spec files. Of cource, only one is actually used.
dist-hook:
	for specfile in *.spec; do \
		if test -f $$specfile; then \
			cp -p $$specfile $(distdir); \
		fi \
	done

