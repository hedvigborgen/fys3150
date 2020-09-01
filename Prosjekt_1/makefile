all: compile execute

compile:
	c++ -o main.exe project1.cpp

execute:
	./main.exe

all_alternative: compile_alternative execute

compile_alternative:
	c++ -o main.out $(wildcard *.cpp)
