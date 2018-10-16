#include <iostream>
#include <array>
#include <algorithm>

int main() {

    std::array<float, 4> TMP{{20, -3.14/2, 5, -3.14/2}};

    // Add mu to every element in the matrix
    float mu = 1;
    transform(TMP.begin(), TMP.end(), TMP.begin(),
	      bind2nd(std::plus<double>(), mu));

    for(const auto& text : TMP) {
	std::cout << text << std::endl;
    }

}
