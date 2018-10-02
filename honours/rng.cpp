#include <iostream>
#include <string>
#include <vector>
#include <iostream>
using namespace std;
#include <sstream>
#include <fstream>
#include <unistd.h>
#include <stdlib.h>
#include <fcntl.h>


unsigned int devrand(void){

	int fn;
	unsigned int r;
	fn = open("/dev/urandom", O_RDONLY);
	if (fn == -1)
	exit(-1); /* Failed! */
	if (read(fn, &r, 4) != 4)
	exit(-1); /* Failed! */
	close(fn);
	return r;

}
/* Initialise KISS generator using /dev/urandom */
unsigned int* init_KISS(){

	unsigned int x,y,c,z;
	x = devrand();
	while (!(y = devrand())); /* y must not be zero! */
	z = devrand();
	/* We don’t really need to set c as well but let's anyway… */
	/* NOTE: offset c by 1 to avoid z=c=0 */
	c = devrand() % 698769068 + 1; /* Should be less than 698769069 */
	//printf("%d\n %d\n %d\n %d\n",x,y,z,c );
	unsigned int seeds[4] = {x,y,z,c};
	return seeds;

}
int KISS(int size){

  	//unsigned int x = (getpid()*time(NULL)),y = rand(),z = time(NULL),c = time(NULL); /* Seed variables */
  	unsigned int x,y,z,c;
  	unsigned int* seed = init_KISS();
  	x=seed[0];
  	y=seed[1];
  	z=seed[2];
  	c=seed[3];
    unsigned long long t, a = 698769069ULL;
    x = 69069*x+12345;
    y ^= (y<<13); y ^= (y>>17); y ^= (y<<5); /* y must never be set to zero! */
    t = a*z+c; c = (t>>32); /* Also avoid setting z=c=0! */
    return (x+y+(z=t)) % size;

  }
std::vector <std::vector<std::string> > readCSVtoMatrix(string filename){

	std::ifstream csv(filename);
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
	 return items;
}

std::vector <std::vector<std::string> > randomSample(std::vector <std::vector<std::string> > population){
	unsigned int rand;
	std::vector <std::vector<std::string> > sample(
    population.size(),
    std::vector<string>('0'));
	for (int i=0;i<population.size();i++){
		rand = KISS(population.size());
		sample[i][0]=population[rand][0];	
		sample[i][1]=population[rand][1];	
	}
	return sample;

}

int main (int argc, char const *argv[]) {

	  std::vector <std::vector<std::string> > population = readCSVtoMatrix("bivariate_data.csv");
	  std::vector <std::vector<std::string> > sample = randomSample(population);
	  for(int i=0;i<population.size(); i++){
	  	for(int j=0;j<2; j++){
	  		cout << sample[i][j] + ",";
	  	}
	  	cout << endl;
	  	
	  }

}



 