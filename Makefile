all: decompose

decompose: decompose.cc
	g++ decompose.cc -o decompose

test: decompose
	./decompose > latest.txt
	diff simple.txt latest.txt
