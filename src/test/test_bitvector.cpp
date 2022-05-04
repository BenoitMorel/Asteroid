#include <iostream>
#include <maths/bitvector.hpp>
#include <vector>
#include <set>


static void testCount()
{
  std::cout << "Testing bitvector count function..." << std::endl;
  unsigned int size = 50;
  BitVector v1(size, false);
  assert(v1.count() == 0);
  for (unsigned int i = 0; i < v1.size(); ++i) {
    if (i % 2) {
      v1.set(i);
    }
  }
  assert(v1.count() == size / 2);
  BitVector v2(50, true);
  assert(v2.count() == 50);
  std::cout << "Success!" << std::endl;
}


static void testSetOperations()
{
  std::cout << "Testing bitvector set operations..." << std::endl;
  unsigned int size = 50;
  BitVector v1(size, false);
  BitVector v2(size, false);
  BitVector vunion(size, false);
  BitVector vinter(size, false);
  v1.set(2);
  v1.set(10);
  v1.set(12);
  v2.set(3);
  v2.set(10);
  v2.set(13);
  vunion.set(2);
  vunion.set(3);
  vunion.set(10);
  vunion.set(12);
  vunion.set(13);
  vinter.set(10);
  assert((v1 & v2) == vinter);
  assert((v1 | v2) == vunion); 
  std::cout << "Success!" << std::endl; 
}


static void testHashtable()
{
  std::cout << "Testing bitvector hashtable" << std::endl;
  std::vector<BitVector> vec;
  vec.push_back(BitVector("00011100"));
  vec.push_back(BitVector("00011101"));
  vec.push_back(BitVector("00011110"));
  vec.push_back(BitVector("01011100"));
  vec.push_back(BitVector("11011100"));
  std::set<BitVector> hashtable;
  for (auto v: vec) {
    hashtable.insert(v);
  }
  assert(vec.size() == hashtable.size());
  for (auto v: vec) {
    assert(hashtable.find(v) != hashtable.end());
  }
  BitVector notIn("00110011");
  assert(hashtable.find(notIn) == hashtable.end());
  std::set<BitVector> hashtable2 = hashtable;
  hashtable2.insert(BitVector("00011100"));
  assert(hashtable == hashtable2);
  std::cout << "Success!" << std::endl; 
}


int main(int, char **)
{
  testCount();
  testSetOperations();
  testHashtable();
  return 0;
}




