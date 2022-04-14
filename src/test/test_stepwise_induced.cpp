#include <iostream>

#include <vector>
#include <maths/bitvector.hpp>
#include <maths/Random.hpp>
#include <trees/StepwiseTree.hpp>
#include <trees/InducedStepwiseTree.hpp>

int main(int argc, char ** argv)
{
  Random::setSeed(42);
  auto ITERATIONS = 1;
  auto LABELS = 100;
  for (unsigned int iter = 0; iter < ITERATIONS; ++iter) {
    std::cout << "NEW ITERATION " << iter << std::endl;
    std::vector<std::string> labels;
    for (unsigned int i = 0; i < LABELS; ++i) {
      labels.push_back(std::string("N") + std::to_string(i));
    }
    StepwiseTree superTree;
    BitVector coverage(labels.size(), false);
    for (unsigned int bit = 0; bit < coverage.size(); ++bit) {
      if (Random::getInt() % 2) {
        coverage.set(bit);
      }
    }
    if (coverage.count() < 2) {
      std::cout << "SKIPPING CURRENT ITERATION" << std::endl;
      continue;
    }
    InducedStepwiseTree inducedTree(superTree, coverage);
    superTree.addListener(&inducedTree); 
    for (unsigned int i = 0; i < labels.size(); ++i) {
      const auto &label = labels[i];
      std::cout << "Adding label " << label << " " << (coverage[i]) << std::endl;
      superTree.addLeaf(i, superTree.getRandomBranch());
      std::cout << superTree.getNewickString(labels) << std::endl;
      std::cout << inducedTree.getWrappedTree().getNewickString(labels) << std::endl;
      std::cout << std::endl;
    }
  }
  return 0;
}
