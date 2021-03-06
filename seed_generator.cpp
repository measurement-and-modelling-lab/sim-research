#include <iostream>	// cout
#include <fcntl.h>	// open
#include <unistd.h>	// read

using namespace std;

unsigned int devrand(void) {

  int fn;
  unsigned int r;
  fn = open("/dev/urandom", O_RDONLY);
  if (fn == -1)
    exit(-1); // Failed!
  if (read(fn, &r, 4) != 4)
    exit(-1); // Failed!
  close(fn);

  return r;
}

int KISS() {

  int x, y, z, c;
  x = devrand();
  while (!(y = devrand()));
  z = devrand();
  c = devrand() % 698769068 + 1;
  unsigned long long t, a = 698769069ULL;
  x = 69069 * x + 12345;
  y ^= (y << 13);
  y ^= (y >> 17);
  y ^= (y << 5); // y must never be set to zero!
  t = a * z + c;
  c = (t >> 32); // Also avoid setting z=c=0!

  return x + y + (z = t);
}

int main(int argc, char const *argv[]) {

    int n = atof(argv[1]); // number of seeds to produce

  for (int i = 0; i < n; i++) {

    cout << KISS() << endl;

  }
}
