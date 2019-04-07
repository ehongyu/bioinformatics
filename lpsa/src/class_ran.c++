//================================================================
//
//    RANDOM NUMBER GENERATOR    GGL
//    original reference for algorithm
//    P. Lewis, A. Goodman, J. Miller,IBM Sys. J., 2, 136 (1969)

#include "head.h"

GGL::GGL(TP initSeed)
{
  Seed=initSeed;
}
void GGL::SetSeed(TP newSeed)
{
  Seed=newSeed;
}

TP GGL::operator()()
{
  Seed=16807*Seed;
  Seed=Seed-int(Seed/TP(2147483647))*TP(2147483647);
  return Seed/2147483647;
}
