
program: main.cpp 
	c++ -std=c++11 -Wall main.cpp jacobi.cpp bisect.cpp compareArmadillo.cpp utils.cpp linalgUtils.cpp -L/usr/local/lib -I/usr/local/include -larmadillo -llapack -lblas -o main.exe

test: test.cpp
	c++ -std=c++11 -Wall test.cpp jacobi.cpp bisect.cpp utils.cpp linalgUtils.cpp compareArmadillo.cpp -L/usr/local/lib -I/usr/local/include -larmadillo -llapack -lblas -o test.exe
