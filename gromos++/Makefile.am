# This is the File that generates Makefile.in
# through automake.

docdir = $(prefix)/doc/$(PACKAGE)-$(VERSION)

SUBDIRS = src \
          programs \
	  data \
	  examples

EXTRA_DIST =	doc/doxygen.conf.in \
	     	doc/tail.html

doc:	doc/doxygen.conf.in programs/0index.doxy
	sed "s/<SRCDIR>/${top_srcdir}/g" ${top_srcdir}/doc/doxygen.conf.in > doxygen.conf
	mkdir -p doc
	doxygen doxygen.conf

doc-install:	doc
	install -d $(docdir)
	cp -r doc/html $(docdir)
