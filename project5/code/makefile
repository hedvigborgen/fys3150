all: compile execute

compile:
	c++ -std=c++11 -o main.exe $(wildcard *.cpp) -larmadillo -O3 -march=native -Xpreprocessor -fopenmp -lomp

execute:
	./main.exe 1 0 1 100000 1000000 1 

execute1:
	./main.exe 20 0.6 0.2 1000 100000 1 0 "at the end"
