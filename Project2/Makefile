
program: main.cpp func.cpp utils.cpp
	c++ -std=c++11 -Wall main.cpp func.cpp utils.cpp -o main.exe

test: test.cpp func.cpp
	c++ -std=c++11 -Wall test.cpp  testFunc.cpp func.cpp utils.cpp bisect.cpp linalgUtils.cpp -L/usr/local/lib -I/usr/local/include -larmadillo -llapack -lblas -o test.exe
