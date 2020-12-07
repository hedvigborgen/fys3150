all: compile execute

compile:
	c++ -std=c++11 -o main.exe $(wildcard *.cpp) -larmadillo -O3 -march=native -Xpreprocessor -fopenmp -lomp

execute:
	./main.exe 20 1000 100000 1 1 ../output/outputfile.dat
