//
// Azar
//
// 11Nov99 primera version
// 30Nov99 header externo y sobrecarga de randomize
//________________________
//JR//

#include "azar.h"

void randomize(unsigned int seed)
{
  srand(seed);	// Inicializa el random
}

void randomize()
{
  srand(time(0) * getpid());	// Inicializa el random
}



