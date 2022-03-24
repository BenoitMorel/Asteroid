#pragma once

#include <string>

struct Arguments {
  Arguments(int argc, char * argv[]);
  void printHelp();
  void printCommand();
 
  int argc;
  char ** argv;

  std::string distanceMatrixListFile;
  std::string prefix;
  // if this argument is not set, we assume
  // single-copy families and that gene names
  // correspond to species names
  std::string geneSpeciesMapping;
  unsigned int seed;

};
