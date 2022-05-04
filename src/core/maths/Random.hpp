#pragma once
#include <random>


class Random {
public:
  Random() = delete;
  static void setSeed(unsigned int seed);
  static int getInt(); 
  static int getInt(int max); 
  static double getProba();
  static std::mt19937_64 &getRNG() {return _d._rng;}
private:
  struct Distributions {
    std::mt19937_64 _rng;
    std::uniform_int_distribution<int> _unii;
    std::uniform_real_distribution<double> _uniproba;
  };

  static Distributions _d;
  
};
