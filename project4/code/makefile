all: compile execute

compile:
	c++ -o main.exe $(wildcard *.cpp) -larmadillo -O3 -march=native

execute:
	./main.exe 20 1.0