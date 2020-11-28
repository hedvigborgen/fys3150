all: compile execute4c 
part4c: compile execute4c
part4f: compile execute4f
test: compile execute_test

compile:
	c++ -std=c++11 -o main.exe $(wildcard *.cpp) -larmadillo -O3 -march=native -Xpreprocessor -fopenmp -lomp

execute4c:
	./main.exe c 100000 2 2 1.0

execute4f:
	./main.exe f 100000 2 20

execute_test:
	./main.exe test 100000 2 2 1.0

