sim:
	g++ simulation.cpp -I ./include/optim ./functions/serafini2019.cpp ./functions/getIntermediateP.cpp ./functions/getSample.cpp ./functions/compute4thOrderMoments.cpp ./functions/counsell2015.cpp ./functions/kolmogorovD.cpp ./functions/kurtosis.cpp ./functions/skewness.cpp ./functions/fleishman1978.cpp -larmadillo -pthread -std=c++11 -o simulation




funcs:
	g++ -c ./functions/serafini2019.cpp ./functions/getIntermediateP.cpp ./functions/getSample.cpp ./functions/counsell2015.cpp ./functions/compute4thOrderMoments.cpp ./functions/kolmogorovD.cpp ./functions/skewness.cpp ./functions/kurtosis.cpp


sim_short:

	g++ simulation.cpp  -I ./functions serafini2019.o getIntermediateP.o getSample.o compute4thOrderMoments.o counsell2015.o kolmogorovD.o kurtosis.o skewness.o -larmadillo -pthread -std=c++11 -o simulation