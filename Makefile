all: decompose

decompose: decompose.cc
	g++ decompose.cc -o decompose
