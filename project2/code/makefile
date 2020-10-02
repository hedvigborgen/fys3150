all_test: test_compile execute

compile:
	c++ -o main.exe main.cpp eigensolver.cpp -larmadillo -lblas -llapack

test_compile:
	c++ -o main.exe test.cpp eigensolver.cpp -larmadillo -lblas -llapack


execute:
	./main.exe
