all: compile execute

compile:
	c++ -std=c++11 -o main.exe $(wildcard *.cpp) -larmadillo -O3 -march=native -Xpreprocessor -fopenmp -lomp

execute:
	./main.exe c 100000 2 2 1.0