all: compile execute

compile:
	c++ -o main.exe $(wildcard *.cpp) -larmadillo -O3 -march=native -fopenmp

execute:
	./main.exe 20 2 100000