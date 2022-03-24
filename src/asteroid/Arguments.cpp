#include "Arguments.hpp"
#include <IO/Logger.hpp>



Arguments::Arguments(int argc, char * argv[]):
  argc(argc),
  argv(argv),
  randomStartingTrees(1),
  bootstrapTrees(0),
  seed(1),
  noCorrection(false),
  minBL(-1.0)
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
    } else if (arg == "-r" || arg == "--random-starting-trees") {
      randomStartingTrees = atoi(argv[++i]);
    } else if (arg == "-b" || arg == "--bs-trees") {
      bootstrapTrees = atoi(argv[++i]);
    } else if (arg == "-p" || arg == "--prefix") {
      prefix = std::string(argv[++i]);
    } else if (arg == "-m" || arg == "--gene-species-mapping") {
      geneSpeciesMapping = std::string(argv[++i]);
    } else if (arg == "--seed") {
      seed = atoi(argv[++i]);
    } else if (arg == "-n" || arg == "--no-correction") {
      noCorrection = true;
    } else if (arg == "--min-bl") {
      minBL = atof(argv[++i]);
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
  Logger::info << "-r | --random-starting-trees <INTEGER>\t Number of starting random trees. 0 to start from an ASTRID tree. Default value is 1." << std::endl;
  Logger::info << "-b | --bs-trees <INTEGER>             \t Number of bootstrap trees to compute. Default is 0." << std::endl;
  Logger::info << "-p | --prefix <STRING>                \t Prefix for the output files." << std::endl;
  Logger::info << "-m | --gene-species-mapping <STRING>  \t Path to the mapping file. If unset, Asteroid assumes that the gene tree leaf labels correspond to the species names. See our wiki for the format." << std::endl;
  Logger::info << "--seed <INTEGER>                      \t Random seed. Default is 1." << std::endl;
  Logger::info << "-n | --no-correction                  \t When set, Asteroid disables its missing data correction, and the algorithm is then equivalent to ASTRID." << std::endl;
  Logger::info << "--min-bl <REAL>                       \t Branches with a length under this value will be contracted into polotomies before computing the distance matrices. DEfault value is -1.0 (no contraction)" << std::endl;
}

void Arguments::printCommand() {
  Logger::info << "Asteroid was called as follow:" << std::endl;
  for (int i = 0; i < argc; ++i) {
    Logger::info << argv[i] << " ";
  }
  Logger::info << std::endl << std::endl;
}

