
1D:
	clear
	c++ -std=c++11 -Wall main1D.cpp triDiaSolver.cpp PDE.cpp utils.cpp -o main1D.exe 
	./main1D.exe
	python plot_1D.py
	mv t1_* t2_* 1D_results/

2D:
	clear
	c++ -std=c++11 -Wall main2D.cpp triDiaSolver.cpp PDE_2D.cpp utils.cpp -o main2D.exe 
	./main1D.exe


clean:
	rm *.txt *.png *.pdf