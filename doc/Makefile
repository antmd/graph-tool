all:
	sphinx-build -E -b html . build

latex:
	sphinx-build -E -b latex . build

test:
	OMP_NUM_THREADS=1 sphinx-build -b doctest . build

clean:
	rm -rf build

