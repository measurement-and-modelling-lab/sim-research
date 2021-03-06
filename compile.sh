g++ simulation.cpp				\
    -I ./functions				\
    ./functions/serafini2019.cpp		\
    ./functions/getIntermediateP.cpp		\
    ./functions/getSample.cpp		        \
    ./functions/compute4thOrderMoments.cpp	\
    ./functions/counsell2015.cpp		\
    ./functions/kolmogorovD.cpp			\
    ./functions/kurtosis.cpp			\
    ./functions/skewness.cpp			\
    -larmadillo -std=c++11 -o simulation

g++ seed_generator.cpp			\
    -std=c++11 -o seed_generator
