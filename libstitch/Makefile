CC_SETUP=OMPI_CC=gcc
CXX_SETUP=OMPI_MPICXX=g++
CC=${CC_SETUP} /usr/local/bin/mpicc
CXX=${CXX_SETUP} /usr/local/bin/mpicxx
TARGET_DIR=../../stitch-master/local/

stitch_test: stitch_test.c stitch.lib
	${CC} -o stitch_test stitch_test.c -DSTITCH_PARALLEL -L. -lstitch -ldl -lpthread

stitch.lib: stitch.c stitch.h sqlite3.c sqlite3.h
	CC="${CC}" CXX="${CXX}" ${CC} -c stitch.c sqlite3.c -DSTITCH_PARALLEL -lstitch -ldl -lpthread
	ar rc libstitch.a stitch.o sqlite3.o

stitch_module: stitch.c stitch.h sqlite3.c sqlite3.h
	cd .. ; CC="gcc" CXX="g++" python3 setup.py build && CC="gcc" CXX="g++" python3 setup.py install

clean:
	cd .. ; python3 setup.py clean
	rm -f *.o libstitch.a *.so stitch_test
	rm -rf ${TARGET_DIR}
	rm -rf ../build
	find .. -name __pycache__ -exec rm -rf {} \; ; find .. -name "*.pyc" -exec rm -rf {} \; ; rm -f *.st

run:
	rm -f test.st ; mpirun -np 4 ./stitch_test
