PDFIMAGES = $(shell ls *.pdf)
PNGIMAGES = ${PDFIMAGES:.pdf=.png}
HUGEIMAGES = $(shell find . -maxdepth 1 -name "*.pdf"  -size +1000k)
ORIGIMAGES = ${HUGEIMAGES:.pdf=.pdf-orig}

PNGCONV = gm convert -density 600 -resample 50 -trim -antialias -quality 9 -filter Cubic

all: $(ORIGIMAGES) $(PNGIMAGES)
	sphinx-build -E -b html . build

%.pdf-orig: %.pdf
	cp $< $@ 
	gm mogrify -resize 600 -trim -antialias -filter Cubic $<

%.png: %.pdf
	test -d $<-orig && $(PNGCONV) $<-orig $@  || true
	test -d $<-orig || $(PNGCONV) $< $@ || true

latex:
	sphinx-build -E -b latex . build

test:
	OMP_NUM_THREADS=1 sphinx-build -b doctest . build

push:
	rsync -rEvpLz build/* root@skewed.de:/var/www/graph-tool-doc/

push-dev:
	rsync -rEvpLz build/* root@skewed.de:/var/www/graph-tool-doc/dev/

clean:
	rm -rf build

