
#include <iostream>
#include <trees/PLLUnrootedTree.hpp>
#include <DistanceMethods/InternodeDistance.hpp>
#include <map>
#include <vector>
#include <climits>


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
  UIntMatrix speciesMatrix(idToLabel.size(),
      std::vector<unsigned int>(idToLabel.size(), 0));

  InternodeDistance::computeFromSpeciesTree(tree, speciesMatrix);

  assert(speciesMatrix[0][0] == 0);
  assert(speciesMatrix[0][1] == 1);
  assert(speciesMatrix[1][0] == 1);
  assert(speciesMatrix[0][2] == 3);
  assert(speciesMatrix[2][0] == 3);
  assert(speciesMatrix[3][0] == 3);
  assert(speciesMatrix[0][3] == 3);
  assert(speciesMatrix[4][0] == 2);
  assert(speciesMatrix[0][4] == 2);
  assert(speciesMatrix[1][2] == 3);
  assert(speciesMatrix[2][1] == 3);
  assert(speciesMatrix[3][1] == 3);
  assert(speciesMatrix[1][3] == 3);
  assert(speciesMatrix[4][1] == 2);
  assert(speciesMatrix[1][4] == 2);
  assert(speciesMatrix[2][4] == 2);
  assert(speciesMatrix[4][2] == 2);
  assert(speciesMatrix[3][4] == 2);
  assert(speciesMatrix[4][3] == 2);


  IDParam params;
  params.useBL = true;
  DistanceMatrix blMatrix(idToLabel.size(),
      std::vector<double>(idToLabel.size(), std::numeric_limits<double>::infinity()));
  InternodeDistance::computeFromGeneTree(tree, blMatrix, params);
  assert(blMatrix[0][0] == 0.0);
  std::cerr << blMatrix[0][1] << std::endl;
  assert(blMatrix[0][1] == 2.0);
  assert(blMatrix[1][0] == 2.0);
  assert(blMatrix[0][2] == 4.0);
  assert(blMatrix[2][0] == 4.0);
  assert(blMatrix[3][0] == 4.0);
  
  assert(blMatrix[0][3] == 4.0);
  assert(blMatrix[4][0] == 12.0);
  assert(blMatrix[0][4] == 12.0);
  assert(blMatrix[1][2] == 4.0);
  assert(blMatrix[2][1] == 4.0);
  assert(blMatrix[3][1] == 4.0);
  assert(blMatrix[1][3] == 4.0);
  assert(blMatrix[4][1] == 12.0);
  assert(blMatrix[1][4] == 12.0);
  assert(blMatrix[2][4] == 12.0);
  assert(blMatrix[4][2] == 12.0);
  assert(blMatrix[3][4] == 12.0);
  assert(blMatrix[4][3] == 12.0);

  std::cout << "Success" << std::endl;
}

int main(int, char **)
{
  test_internode(); 
  return 0;
}

