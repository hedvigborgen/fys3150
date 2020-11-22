all: compile execute

compile:
	c++ -o main.exe $(wildcard *.cpp) -larmadillo

execute:
	./main.exe 2 100 1.0 >> ../output/data.dat