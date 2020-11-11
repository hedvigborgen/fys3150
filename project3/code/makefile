compile:
	c++ -o main.exe $(wildcard *.cpp) --std=c++11
execute:
	./main.exe 1 100000 0.001 ../input/ten_bodies.txt 10 2

example: compile execute

test_beta:
	./main.exe 2 1000 0.001 ../input/two_bodies_elliptical.txt 2 2

precession:
	./main.exe 3 1000000 0.0001 ../input/precession_mercury.txt 2 2

