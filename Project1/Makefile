
program: main.cpp func.cpp
	c++ -std=c++11 -Wall main.cpp func.cpp -L/usr/local/lib -I/usr/local/include -larmadillo -llapack -lblas -o main.exe

test: test.cpp func.cpp
	c++ -std=c++11 -Wall test.cpp func.cpp -L/usr/local/lib -I/usr/local/include -larmadillo -llapack -lblas -o test.exe
