all:
	g++ --std=c++14 -O2 -c edlib/src/edlib.cpp -o edlib.o -I edlib/include
	g++ --std=c++14 -O2 -c rr.cpp -o rr.o -I edlib/include
	g++ rr.o edlib.o -o rr -lz

clean :
	-rm *.o rr
