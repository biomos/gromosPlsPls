#include <iostream>
#include <fstream>
#include <sstream>
#include <string.h>
#include <stdio.h>
#include <vector>
#include <math.h>
#include <cstdlib>
#include <iomanip>
#include <cassert>



using namespace std;


vector <double> boxmul() {

//2 random vectors in a sphere, obtained by using the box-muller transform

//  Gives 2 random variables with a gaussian distribution, by using the polar Box-Muller method
//  Ref: Box,  G. E. P. and Muller, M. E. Ann. Math. Stat. 29 (1958)
//  610-611

//  The implementation is of the Marsaglia variant (ref??)

  vector <double> gaussrand;
  gaussrand.resize(2);

  double r1;
  double r2;

  double w = 10;

  while ((w>1)||(w==0)) {
    r1 =2*double(rand())/RAND_MAX -1;
    r2 = 2*double(rand())/RAND_MAX -1;

    w = r1*r1+r2*r2;
  }

  w =sqrt(-2*log(w)/w);

  gaussrand[0] = w*r1;
  gaussrand[1] = w*r2;

  return gaussrand;

}
