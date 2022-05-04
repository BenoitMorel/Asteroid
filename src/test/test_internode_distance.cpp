
#include <iostream>
#include <trees/PLLUnrootedTree.hpp>
#include <DistanceMethods/InternodeDistance.hpp>
#include <map>
#include <vector>

static void test_internode()
{
  std::cout << "Testing internode distance" << std::endl;
  PLLUnrootedTree tree("((a:1.0,b:1.0):1.0,(c:1.0,d:1.0):1.0,e:10.0);", false);

  std::vector<std::string> idToLabel;
  idToLabel.push_back("a");
  idToLabel.push_back("b");
  idToLabel.push_back("c");
  idToLabel.push_back("d");
  idToLabel.push_back("e");
  std::map<std::string, unsigned int> labelToId;
  for (unsigned int i = 0; i < idToLabel.size(); ++i) {
    labelToId.insert({idToLabel[i], i});
  }
  for (auto leaf: tree.getLeaves()) {
    leaf->data = reinterpret_cast<void *>(static_cast<intptr_t>(labelToId[leaf->label]));
  }
  /*
  DistanceMatrix matrix(idToLabel.size(),
      std::vector<double>(idToLabel.size(),
        std::numeric_limits<double>infinity());

    */  
  //InternodeDistance::computeFromGeneTree
  std::cout << "Success" << std::endl;
}

int main(int, char **)
{
  test_internode(); 
  return 0;
}

