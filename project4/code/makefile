all: compile execute4c execute4d_T1_ordered execute4d_T2_ordered execute4d_T1_random execute4d_T2_random execute4e_T1 execute4e_T2 execute4f

part4c: compile execute4c
part4d: compile execute4d_T1_ordered execute4d_T2_ordered execute4d_T1_random execute4d_T2_random
part4e: compile execute4e_T1 execute4e_T2
part4f: compile execute4f


compile:
	c++ -std=c++11 -o main.exe $(wildcard *.cpp) -larmadillo -O3 -march=native -Xpreprocessor -fopenmp -lomp


execute4c:
	./main.exe c 100000 2 1.0 2


execute4d_T1_ordered:
	./main.exe d 100000 20 1.0 1
execute4d_T2_ordered:
	./main.exe d 100000 20 2.4 1
execute4d_T1_random:
	./main.exe d 100000 20 1.00 2
execute4d_T2_random:
	./main.exe d 100000 20 2.40 2


execute4e_T1:
	./main.exe e 100000 20 1.0
execute4e_T2:
	./main.exe e 100000 20 2.4


execute4f:
	./main.exe f 100000 20 2
