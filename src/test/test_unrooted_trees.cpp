#include <iostream>
#include <trees/PLLUnrootedTree.hpp>
#include <set>

static void testParser() 
{
  std::vector<std::string> newicks;
  // binary tree
  newicks.push_back("((a,b),(c,d), (e,f));");
  // binary tree with branch lengths
  newicks.push_back("((a:0.5,b),(c,d:1.0):3.5555, (e:1.0,f));");
  // binary tree with labels
  newicks.push_back("((a,b)hey,(c,d)ho, (e,f)5);");
  // tree with polytomies
  newicks.push_back("((a,b,g),(c,d,y), (e,f,z,s,d,g), (nifnif, nafnaf, noufnouf));");

  for (auto newick: newicks) {
    std::cout << "Reading tree " << newick << std::endl;
    PLLUnrootedTree tree(newick, false);
    assert(tree.getLeavesNumber() > 0);
    std::cout << "Success!" << std::endl;
  }
}

static void testIsomorphic()
{
  PLLUnrootedTree tree1_1("((a,b),(c,d),(e,f));", false);
  PLLUnrootedTree tree1_2("((e,f),(d,c), (b,a));", false);
  PLLUnrootedTree tree2("((a,b),(c,e),(f,d));", false);
  std::cout << "Check isomorphic trees..." << std::endl;
  assert(PLLUnrootedTree::areIsomorphic(tree1_1, tree1_2));
  std::cout << "Success!" << std::endl;
  std::cout << "Check non-isomorphic trees..." << std::endl;
  assert(!PLLUnrootedTree::areIsomorphic(tree1_1, tree2));
  assert(!PLLUnrootedTree::areIsomorphic(tree1_2, tree2));
  std::cout << "Success!" << std::endl;
}


static void testInducedTree()
{
  std::cout << "Test induced trees..." << std::endl;
  PLLUnrootedTree tree1("((a,b),(c,d),(e,f));", false);
  PLLUnrootedTree tree2("((a,b),(c,f));", false);
  std::set<std::string> inducedLabels;
  inducedLabels.insert("a");
  inducedLabels.insert("b");
  inducedLabels.insert("c");
  inducedLabels.insert("f");
  BitVector okNodes(tree1.getLeavesNumber(), false);
  for (auto leaf: tree1.getLeaves()) {
    if (inducedLabels.end() != inducedLabels.find(leaf->label)) {
      okNodes.set(leaf->node_index);
    }
  }
  auto induced_1_2 = tree1.getInducedTree(okNodes);
  assert(PLLUnrootedTree::areIsomorphic(tree2, *induced_1_2));
  std::cout << "Success!" << std::endl;
}


int main(int, char **)
{
  testParser();
  testIsomorphic();
  testInducedTree(); 

  return 0;
}

