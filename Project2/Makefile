
program: main.cpp func.cpp utils.cpp
	c++ -std=c++11 -Wall main.cpp jacobi.cpp utils.cpp -o main.exe

test: test.cpp
	c++ -std=c++11 -Wall test.cpp  testLinalg.cpp jacobi.cpp utils.cpp bisect.cpp linalgUtils.cpp -L/usr/local/lib -I/usr/local/include -larmadillo -llapack -lblas -o test.exe
