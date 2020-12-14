all: compile execute1 execute2 execute3

compile:
	c++ -std=c++11 -o main.exe $(wildcard *.cpp) -larmadillo -O3 -march=native -Xpreprocessor -fopenmp -lomp

execute1:
	./main.exe MCCs 1000000 0 8 0.2 0.1

execute2:
	./main.exe Alpha 1000000 0 8 0.2 0.1 1.0

execute3:
	./main.exe Parameters 1000000 1 1.0 0.85 0.5 1
