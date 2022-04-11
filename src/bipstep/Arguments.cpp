#include "Arguments.hpp"
#include <IO/Logger.hpp>



Arguments::Arguments(int argc, char * argv[]):
  argc(argc),
  argv(argv),
  randomInsertion(false),
  seed(1)
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
    } else if (arg == "-p" || arg == "--prefix") {
      prefix = std::string(argv[++i]);
    } else if (arg == "-m" || arg == "--gene-species-mapping") {
      geneSpeciesMapping = std::string(argv[++i]);
    } else if (arg == "--random-insertion") {
      randomInsertion = true;
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
  Logger::info << "Arguments:" << std::endl;
  Logger::info << "-i | --input-gene-trees <STRING>      \t Path to file containing one gene tree per line, in newick format" << std::endl;
}

void Arguments::printCommand() {
  Logger::info << "Asteroid was called as follow:" << std::endl;
  for (int i = 0; i < argc; ++i) {
    Logger::info << argv[i] << " ";
  }
  Logger::info << std::endl << std::endl;
}

