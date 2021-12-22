#pragma once

#include <string>

class Arguments {
public:

  Arguments(int argc, char * argv[]);
  void printHelp();
  void printCommand();
 
private:
  int argc;
  char ** argv;

  std::string inputGeneTreeFile;
  std::string startingSpeciesTree;
  std::string outputSpeciesTree;

  unsigned int seed;
public:
};
