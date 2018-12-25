g++ main.cpp					\
    -I ./functions				\
    ./functions/ADF.cpp				\
    ./functions/ValeMaurelli.cpp		\
    ./functions/compute4thOrderMoments.cpp	\
    ./functions/counsell.cpp			\
    ./functions/ksD.cpp				\
    ./functions/kurtosis.cpp			\
    ./functions/skewness.cpp			\
    -larmadillo -std=c++11
