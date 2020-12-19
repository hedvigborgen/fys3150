all: compile execute1 execute2 execute3 execute4 execute5 execute6 execute7

compile:
	c++ -std=c++11 -o main.exe $(wildcard *.cpp) -larmadillo -O3 -march=native -Xpreprocessor -fopenmp -lomp

execute1:
	./main.exe MCCs 1000000 0 8 0.2 0.1

execute2:
	./main.exe StepLength 1000000 0 8 0.2 0.1 0.5

execute3:
	./main.exe Parameters 1000000 1 0.85 0.5 1.0

execute4:
	./main.exe Loop 1000000 1.0

execute5:
	./main.exe Virial0 1000000 0.01 0.01

execute6:
	./main.exe Virial1 1000000 0.01 0.01

execute7:
	./main.exe Virial2 1000000 0.01 0.01
