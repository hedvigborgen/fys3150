all: compile execute

compile:
	c++ -o main.exe project1.cpp lib.cpp

compile_error:
	c++ -o error.exe error.cpp

execute:
	./main.exe

execute_error:
	./error.exe

all_alternative: compile_alternative execute

compile_alternative:
	c++ -o main.out $(wildcard *.cpp)
