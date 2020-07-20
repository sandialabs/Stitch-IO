CC=mpicc

stitch_test: stitch_test.c stitch.lib
	python3 setup.py build_ext --inplace
	${CC} -o stitch_test stitch_test.c -L./ -lstitch -ldl -lpthread

stitch.lib: stitch.c stitch.h sqlite3.c sqlite3.h
	${CC} -g -c stitch.c sqlite3.c
	ar rc libstitch.a stitch.o sqlite3.o

all: stitch_test stitch_lib
	${CC} -g -o stitch_test stitch.c sqlite3.c stitch_test.c
##	mpicc -g -o stitch_test stitch.c sqlite3.c stitch_test.c
##	python3 setup.py build
##	python3 setup.py install

clean:
	rm -f *.o libstitch.a *.so stitch_test
	rm -rf ./build
	rm -f *.st

run:
#	rm -f test.st && ./stitch_test
	rm -f test.st && mpirun -np 4 ./stitch_test
