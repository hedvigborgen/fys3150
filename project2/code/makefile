all: compile execute

compile:
	c++ -o main.exe alt_main.cpp eigensolver.cpp -larmadillo -lblas -llapack

test_compile:
	c++ -o main.exe test.cpp eigensolver.cpp -larmadillo -lblas -llapack

execute:
	./main.exe
