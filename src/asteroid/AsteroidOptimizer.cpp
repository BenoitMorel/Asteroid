#include "AsteroidOptimizer.hpp" 
#include <IO/Logger.hpp>
#include <trees/PLLUnrootedTree.hpp>

double AsteroidOptimizer::eval(PLLUnrootedTree &tree)
{
  _lastScore = -_asteroid.computeLength(tree);
  Logger::timed << "score=" << _lastScore << std::endl;
  return _lastScore;
}

static bool sprYeldsSameTree(corax_unode_t *p, corax_unode_t *r)
{
  assert(p);
  assert(r);
  assert(p->next);
  return (r == p) || (r == p->next) || (r == p->next->next)
    || (r == p->back) || (r == p->next->back) || (r == p->next->next->back);
}

bool isSPRMoveValid(PLLUnrootedTree &tree,
    corax_unode_t *prune, 
    corax_unode_t *regraft)
{
  // regraft should not be a child of prune
  auto pruneChildren = tree.getPostOrderNodesFrom(prune->back);
  for (auto child: pruneChildren) {
    if (regraft == child) {
      return false;
    }
    if (regraft->back == child) {
      return false;
    }
  }
  return !sprYeldsSameTree(prune, regraft);
}

bool isInScores(const std::vector<double> &scores, double score)
{
  for (auto s: scores) {
    if (fabs(s-score) < 0.000001) {
      return true;
    }
  }
  return false;
}

void addInvolvedNode(corax_unode_t *node, 
    std::unordered_set<corax_unode_t *> &involved)
{
  involved.insert(node);
  if (node->next) {
    involved.insert(node->next);
    involved.insert(node->next->next);
  }
}

bool wasInvolved(corax_unode_t *node,
    std::unordered_set<corax_unode_t *> &involved)
{
  return involved.find(node) != involved.end();
}
  
bool AsteroidOptimizer::computeAndApplyBestSPR()
{
  unsigned int maxRadiusWithoutImprovement = 1;
  Logger::timed << "last score " << _lastScore << std::endl;
  std::vector<SPRMove> bestMoves;
  double epsilon = 0.00000001;
  _asteroid.getBestSPR(_speciesTree, 
      maxRadiusWithoutImprovement,
      bestMoves);
  if (bestMoves.size() == 0 || bestMoves[0].score < epsilon) {
    Logger::info << "Local search failed, trying with max radius..." << std::endl;
    _asteroid.getBestSPR(_speciesTree, 
        99999999,
        bestMoves);
    if (bestMoves.size() == 0 || bestMoves[0].score < epsilon) {
      Logger::info << "Global search failed, stopping" << std::endl;
      return false;
    }
  }
  bool better = false;
  unsigned int appliedMoves = 0;
  unsigned int maxAppliedMoves = 9999999; // no max
  corax_tree_rollback_t emptyRollback;
  std::vector<corax_tree_rollback_t> rollbacks;
  std::vector<double> hackScore;
  std::unordered_set<corax_unode_t *> involved;
  double expectedDiff = 0.0;
  for (auto move: bestMoves) {
    if (move.score > epsilon 
        && isSPRMoveValid(_speciesTree, move.pruneNode->back, move.regraftNode) 
        && !isInScores(hackScore, move.score)
        && !wasInvolved(move.pruneNode, involved)
        && !wasInvolved(move.regraftNode, involved)
      ) {
      addInvolvedNode(move.pruneNode->back, involved);
      addInvolvedNode(move.regraftNode, involved);
      rollbacks.push_back(emptyRollback);   
      hackScore.push_back(move.score);
      expectedDiff += move.score;
      auto ok = corax_utree_spr(move.pruneNode->back, 
          move.regraftNode, 
          &rollbacks.back());
      appliedMoves++;
      assert(ok);
      better = true;
    }
    if (appliedMoves >= maxAppliedMoves) {
      break;
    }
  }
  Logger::info << "Moves: " << appliedMoves << std::endl;
  double newScore = -_asteroid.computeLength(_speciesTree);
  double diff = newScore - _lastScore;
  if (appliedMoves) {
    Logger::info << "Expected: " << expectedDiff << " real:" << diff << " highest: " << bestMoves[0].score << std::endl;
  }
  if (newScore < _lastScore + epsilon) {
    Logger::info << "New score " << newScore << " worse than last score " << _lastScore << std::endl;
    for (int i = rollbacks.size() - 1; i >= 1; --i) {
      corax_tree_rollback(&rollbacks[i]);
    }
    newScore = -_asteroid.computeLength(_speciesTree);
    assert(newScore > _lastScore);
  }
  _lastScore = newScore;
  return better;
}


AsteroidOptimizer::AsteroidOptimizer(PLLUnrootedTree &speciesTree,
    const BoolMatrix &perFamilyCoverage,
    const UIntMatrix &gidToSpid,
    const std::vector<DistanceMatrix> &distanceMatrices):
  _speciesTree(speciesTree),
  _asteroid(speciesTree, perFamilyCoverage, gidToSpid, distanceMatrices),
  _lastScore(0.0)
{
  ParallelContext::barrier();
}

double AsteroidOptimizer::optimize()
{
  double startingScore = eval(_speciesTree); 
  Logger::info << "Starting score: " << startingScore << std::endl;
  bool ok = true;
  unsigned int it = 0;
  while (ok) {
    ok = computeAndApplyBestSPR();
    ++it;
  }
  Logger::info << "SPR search stopped after " << it << " iterations" << std::endl;
  return _lastScore;
}




