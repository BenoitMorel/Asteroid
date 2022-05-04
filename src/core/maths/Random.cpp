#include "Random.hpp"

Random::Distributions Random::_d;

void Random::setSeed(unsigned int seed)
{
  _d._rng.seed(seed);
}
int Random::getInt() 
{
  return _d._unii(_d._rng);
}

int Random::getInt(int max)
{
  return getInt() % max;
}

double Random::getProba() 
{
  return _d._uniproba(_d._rng);
}


