
1D:
	clear
	c++ -std=c++11 -Wall main1D.cpp triDiaSolver.cpp PDE.cpp utils.cpp -o main1D.exe 
	./main1D.exe
	python plot_1D.py
	mv t1_dx* t2_dx* 1D_results/

2D:
	clear
	c++ -std=c++11 -Wall main2D.cpp triDiaSolver.cpp PDE_2D.cpp utils.cpp -o main2D.exe 
	./main2D.exe
	python plot_2D.py
	mv t1_h* t2_h* 2D_results/

parallel:
	clear
	mpic++ -std=c++11  -O3 parallel.cpp triDiaSolver.cpp PDE_2D.cpp utils.cpp -o parallel.exe
	mpiexec --hostfile hostfile -np 4 ./parallel.exe
	python plot_2D.py
	mv t1_h* t2_h* 2DP_results/

clean:
	rm *.txt *.png *.pdf