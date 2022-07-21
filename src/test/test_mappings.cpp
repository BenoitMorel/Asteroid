#include <iostream>
#include <IO/GeneSpeciesMapping.hpp>
#include <trees/PLLUnrootedTree.hpp>
#include <set>

static void test_from_tree()
{
  std::cout << "Test mappings from tree" << std::endl;
  PLLUnrootedTree tree("((a,b),(c,d),(e,f));", false);
  GeneSpeciesMapping mapping;
  mapping.fillFromGeneTree(tree);
  assert(mapping.getCoveredSpecies().size() == 6);
  std::set<std::string> covered;
  covered.insert("a");
  covered.insert("b");
  covered.insert("c");
  covered.insert("d");
  covered.insert("e");
  covered.insert("f");
  for (auto c: covered) {
    assert(mapping.getSpecies(c) == c);
  }   
  std::cout << "Success!" << std::endl;
}

int main(int, char **)
{
  test_from_tree();
  return 0;

}
