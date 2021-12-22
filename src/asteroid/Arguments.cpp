#include "Arguments.hpp"
#include <IO/Logger.hpp>



Arguments::Arguments(int argc, char * argv[]):
  argc(argc),
  argv(argv)
{
  if (argc == 1) {
    printHelp();
    ParallelContext::abort(0);
  }
  for (int i = 1; i < argc; ++i) {
    std::string arg(argv[i]);
    if (arg == "-h" || arg == "--help") {
      printHelp();
      ParallelContext::abort(0);
    } else if (arg == "-i" || arg == "--input-gene-trees") {
      inputGeneTreeFile = std::string(argv[++i]);
    } else if (arg == "-s" || arg == "--starting-species-tree") {
      outputSpeciesTree = std::string(argv[++i]);
    } else if (arg == "-o" || arg == "--output-species-tree") {
      outputSpeciesTree = std::string(argv[++i]);
    } else if (arg == "--seed") {
      seed = atoi(argv[++i]);
    } else {
      Logger::info << "Unrecognized argument " << arg << std::endl;
      Logger::info << "Aborting" << std::endl;
      ParallelContext::abort(1);
    }
  }
}

void Arguments::printHelp() 
{
  Logger::info << "Help message" << std::endl;
}

void Arguments::printCommand() {
  Logger::info << "Asteroid was called as follow:" << std::endl;
  for (int i = 0; i < argc; ++i) {
    Logger::info << argv[i] << " ";
  }
  Logger::info << std::endl << std::endl;
}

