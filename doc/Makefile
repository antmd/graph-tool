all:
	sphinx-build -E -b html . build

test:
	OMP_NUM_THREADS=1 sphinx-build -b doctest . build

push:
	rsync -rEvpLz build/* root@forked.de:/var/www/graph-tool-doc/
