#include <iostream>
#include <string>
#include <vector>
#include <iostream>
using namespace std;
#include <sstream>
#include <fstream>
#include <unistd.h>


int main (int argc, char const *argv[]) {

std::ifstream csv("bivariate_data.csv");
std::string line;
std::vector <std::vector<std::string> > items;

if (csv.is_open()) {
        for (std::string row_line; std::getline(csv, row_line);)
        {
            items.emplace_back();
            std::istringstream row_stream(row_line);
            for(std::string column; std::getline(row_stream, column, ',');){
                items.back().push_back(column);
            }
        }
}
else {
    cout << "Unable to open file";
}
  for(int i=0;i<50; i++){
  	for(int j=0;j<2; j++){
  		cout << items[i][j] + ",";
  	}
  	cout << endl;
  	
  }

}

int KISS()
  {
  	unsigned int x = (getpid()*time(NULL)),y = rand(),z = time(NULL),c = time(NULL); /* Seed variables */
    unsigned long long t, a = 698769069ULL;
    x = 69069*x+12345;
    y ^= (y<<13); y ^= (y>>17); y ^= (y<<5); /* y must never be set to zero! */
    t = a*z+c; c = (t>>32); /* Also avoid setting z=c=0! */
    return (x+y+(z=t)) % 50000;
  }

 