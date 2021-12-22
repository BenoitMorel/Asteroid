#include <iostream>
#include <corax/corax.h>
#include <IO/Logger.hpp>
#include "Arguments.hpp"
//#include <trees/PLLUnrootedTree.hpp>

#include <memory>
#include <set>

/*
using GeneTrees = std::vector<std::shared_ptr<PLLUnrootedTree> >;

void parseTrees(const std::string &inputGeneTreeFile,
    GeneTrees &geneTrees)
{
  // Read all trees
  std::ifstream reader(family.startingGeneTree);
  std::string geneTreeStr;
  while (std::getline(reader, geneTreeStr)) {
    geneTrees.push_back(std::make_shared<PLLUnrootedTree>(
          geneTreeStr, false));
  }
}
*/


int main(int argc, char * argv[])
{
  Logger::timed << "Starting Asteroid..." << std::endl;
  Arguments arg(argc, argv);
  arg.printCommand();
  std::set<std::string> speciesLabels;
  
  Logger::info << "Todo: parallelization" << std::endl;
  Logger::info << "Todo: missing data pattern compression" << std::endl;

  /*
  GeneTrees geneTrees;
  parseTrees(arg.inputGeneTreeFile,
      geneTrees;
  */
  return 0;

}
