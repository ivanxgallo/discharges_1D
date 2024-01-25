#ifndef _AZAR_UTILS_H_
#define _AZAR_UTILS_H_

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <unistd.h>
#include <stdlib.h>

void randomize();
void randomize(unsigned int) ;
inline double dazar(double inf=0.0e0, double sup = 1.0e0){
  return (sup-inf)* double(rand())/double(RAND_MAX)+inf ;
}
inline int iazar(int inf , int sup){
  return int(floor(dazar(inf, sup+1))) ;
}

#endif
