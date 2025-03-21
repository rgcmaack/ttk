/// \ingroup base
/// \class ttk:FTMTree
/// \author Charles Gueunet <charles.gueunet@lip6.fr>
/// \date Dec 2016.
///
///\brief TTK processing package that efficiently computes the
/// contour tree of scalar data and more
/// (data segmentation, topological simplification,
/// persistence diagrams, persistence curves, etc.).
///
///\param dataType Data type of the input scalar field (char, float,
/// etc.).

#include "FTMTree_MT.h"

#include <stack>

#define PRIOR(x)
// #define PRIOR(x) priority(x)

#ifdef __INTEL_COMPILER
#define HIGHER
#endif

#ifndef TTK_ENABLE_OPENMP
#define HIGHER
#endif

#ifdef TTK_ENABLE_OMP_PRIORITY
#define OPTIONAL_PRIORITY(value) priority(value)
#else
#define OPTIONAL_PRIORITY(value)
#endif

using namespace std;
using namespace ttk;
using namespace ftm;

FTMTree_MT::FTMTree_MT(Params *const params,
                       Scalars *const scalars,
                       TreeType type)
  : params_(params), scalars_(scalars) {

  this->setDebugMsgPrefix("FTMtree_MT");

  mt_data_.treeType = type;
}

FTMTree_MT::~FTMTree_MT() {

  // remove UF data structures
  if(mt_data_.ufs) {
    sort(mt_data_.ufs->begin(), mt_data_.ufs->end());
    auto it = unique(mt_data_.ufs->begin(), mt_data_.ufs->end());
    mt_data_.ufs->resize(std::distance(mt_data_.ufs->begin(), it));
    for(auto *addr : *mt_data_.ufs)
      if(addr)
        delete addr;
  }

  // if (mt_data_.propagation) {
  //    Already cleaned by ufs
  //    sort(mt_data_.propagation->begin(), mt_data_.propagation->end());
  //    auto it = unique(mt_data_.propagation->begin(),
  //    mt_data_.propagation->end());
  //    mt_data_.propagation->resize(std::distance(mt_data_.propagation->begin(),
  //    it)); for (auto* addr : *mt_data_.propagation) if(addr) delete addr;
  // }

  // remove containers
  if(mt_data_.superArcs) {
    delete mt_data_.superArcs;
    mt_data_.superArcs = nullptr;
  }
  if(mt_data_.nodes) {
    delete mt_data_.nodes;
    mt_data_.nodes = nullptr;
  }
  if(mt_data_.roots) {
    delete mt_data_.roots;
    mt_data_.roots = nullptr;
  }
  if(mt_data_.leaves) {
    delete mt_data_.leaves;
    mt_data_.leaves = nullptr;
  }
  if(mt_data_.vert2tree) {
    delete mt_data_.vert2tree;
    mt_data_.vert2tree = nullptr;
  }
  if(mt_data_.trunkSegments) {
    delete mt_data_.trunkSegments;
    mt_data_.trunkSegments = nullptr;
  }
  if(mt_data_.visitOrder) {
    delete mt_data_.visitOrder;
    mt_data_.visitOrder = nullptr;
  }
  if(mt_data_.ufs) {
    delete mt_data_.ufs;
    mt_data_.ufs = nullptr;
  }
  if(mt_data_.states) {
    delete mt_data_.states;
    mt_data_.states = nullptr;
  }
  if(mt_data_.propagation) {
    delete mt_data_.propagation;
    mt_data_.propagation = nullptr;
  }
  if(mt_data_.valences) {
    delete mt_data_.valences;
    mt_data_.valences = nullptr;
  }
  if(mt_data_.openedNodes) {
    delete mt_data_.openedNodes;
    mt_data_.openedNodes = nullptr;
  }

#ifdef TTK_ENABLE_FTM_TREE_STATS_TIME
  if(mt_data_.activeTasksStats) {
    delete mt_data_.activeTasksStats;
    mt_data_.activeTasksStats = nullptr;
  }
#endif
}

void FTMTree_MT::buildSegmentation() {

  const idSuperArc nbArcs = mt_data_.superArcs->size();

  // Make reserve

  // SuperArc i correspond to segment i,
  // one arc correspond to one segment
  vector<SimplexId> sizes(nbArcs);

  // get the size of each segment
  const idSuperArc arcChunkSize = getChunkSize(nbArcs);
  const idSuperArc arcChunkNb = getChunkCount(nbArcs);
  for(idSuperArc arcChunkId = 0; arcChunkId < arcChunkNb; ++arcChunkId) {
    // WHY shared(sizes) is needed ??
#ifdef TTK_ENABLE_OPENMP
#pragma omp task firstprivate(arcChunkId) shared(sizes) \
  OPTIONAL_PRIORITY(isPrior())
#endif
    {
      const idSuperArc lowerBound = arcChunkId * arcChunkSize;
      const idSuperArc upperBound
        = min(nbArcs, (arcChunkId + 1) * arcChunkSize);
      for(idSuperArc a = lowerBound; a < upperBound; ++a) {
        sizes[a]
          = max(SimplexId{0}, (*mt_data_.superArcs)[a].getNbVertSeen() - 1);
      }
    }
  }
#ifdef TTK_ENABLE_OPENMP
#pragma omp taskwait
#endif

  // change segments size using the created vector
  mt_data_.segments_.resize(sizes);

  Timer segmentsSet;

  // Fill segments using vert2tree

  // current status of the segmentation of this arc
  vector<SimplexId> posSegm(nbArcs, 0);

  // Segments are connex region of geometrie forming
  // the segmentation (sorted in ascending order)
  const SimplexId nbVert = scalars_->size;
  const SimplexId chunkSize = getChunkSize();
  const SimplexId chunkNb = getChunkCount();
  for(SimplexId chunkId = 0; chunkId < chunkNb; ++chunkId) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp task firstprivate(chunkId) shared(posSegm) \
  OPTIONAL_PRIORITY(isPrior())
#endif
    {
      const SimplexId lowerBound = chunkId * chunkSize;
      const SimplexId upperBound = min(nbVert, (chunkId + 1) * chunkSize);
      for(SimplexId i = lowerBound; i < upperBound; ++i) {
        const auto vert = scalars_->sortedVertices[i];
        if(isCorrespondingArc(vert)) {
          idSuperArc sa = getCorrespondingSuperArcId(vert);
          SimplexId vertToAdd;
          if((*mt_data_.visitOrder)[vert] != nullVertex) {
            // Opposite order for Split Tree
            vertToAdd = (*mt_data_.visitOrder)[vert];
            if(isST())
              vertToAdd = getSuperArc(sa)->getNbVertSeen() - vertToAdd - 2;
            mt_data_.segments_[sa][vertToAdd] = vert;
          } else if(mt_data_.trunkSegments->size() == 0) {
            // MT computation
#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic capture
#endif
            vertToAdd = posSegm[sa]++;
            mt_data_.segments_[sa][vertToAdd] = vert;
          }

        } // end is arc
      } // end for
    } // end task
  }
#ifdef TTK_ENABLE_OPENMP
#pragma omp taskwait
#endif

  printTime(segmentsSet, "segmentation set vertices", 4);

  if(mt_data_.trunkSegments->size() == 0) {
    // sort arc that have been filled by the trunk
    // only for MT
    Timer segmentsSortTime;
    for(idSuperArc a = 0; a < nbArcs; ++a) {
      if(posSegm[a]) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp task firstprivate(a) OPTIONAL_PRIORITY(isPrior())
#endif
        mt_data_.segments_[a].sort(scalars_);
      }
    }
#ifdef TTK_ENABLE_OPENMP
#pragma omp taskwait
#endif
    printTime(segmentsSortTime, "segmentation sort vertices", 4);
  } else {
    // Contour tree: we create the arc segmentation for arcs in the trunk
    Timer segmentsArcTime;
    for(idSuperArc a = 0; a < nbArcs; ++a) {
      // CT computation, we have already the vert list
      if((*mt_data_.trunkSegments)[a].size()) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp task firstprivate(a) OPTIONAL_PRIORITY(isPrior())
#endif
        mt_data_.segments_[a].createFromList(
          scalars_, (*mt_data_.trunkSegments)[a],
          mt_data_.treeType == TreeType::Split);
      }
    }
#ifdef TTK_ENABLE_OPENMP
#pragma omp taskwait
#endif

    printTime(segmentsArcTime, "segmentation arcs lists", 4);
  }

  // Update SuperArc region

  // ST have a segmentation wich is in the reverse-order of its build
  // ST have a segmentation sorted in ascending order as JT
  for(idSuperArc arcChunkId = 0; arcChunkId < arcChunkNb; ++arcChunkId) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp task firstprivate(arcChunkId) OPTIONAL_PRIORITY(isPrior())
#endif
    {
      const idSuperArc lowerBound = arcChunkId * arcChunkSize;
      const idSuperArc upperBound
        = min(nbArcs, (arcChunkId + 1) * arcChunkSize);
      for(idSuperArc a = lowerBound; a < upperBound; ++a) {
        // avoid empty region
        if(mt_data_.segments_[a].size()) {
          (*mt_data_.superArcs)[a].concat(
            mt_data_.segments_[a].begin(), mt_data_.segments_[a].end());
        }
      }
    }
  }
#ifdef TTK_ENABLE_OPENMP
#pragma omp taskwait
#endif
}

FTMTree_MT *FTMTree_MT::clone() const {
  FTMTree_MT *newMT = new FTMTree_MT(params_, scalars_, mt_data_.treeType);

  newMT->mt_data_.superArcs = mt_data_.superArcs;
  newMT->mt_data_.nodes = mt_data_.nodes;
  newMT->mt_data_.leaves = mt_data_.leaves;
  newMT->mt_data_.roots = mt_data_.roots;
  newMT->mt_data_.vert2tree = mt_data_.vert2tree;

  return newMT;
}

void FTMTree_MT::closeArcsUF(idNode closeNode, UF uf) {
  for(const auto &sa : uf->find()->getOpenedArcs()) {
    closeSuperArc(sa, closeNode);
  }
  uf->find()->clearOpenedArcs();
}

void FTMTree_MT::closeSuperArc(idSuperArc superArcId, idNode upNodeId) {
#ifndef TTK_ENABLE_KAMIKAZE

  if(superArcId >= getNumberOfSuperArcs()) {
    cout << "[Merge Tree] closeSuperArc on a inexisting arc !" << endl;
    return;
  }

  if(upNodeId >= getNumberOfNodes()) {
    cout << "[Merge Tree] closeOpenedArc on a inexisting node !" << endl;
    return;
  }

#endif
  (*mt_data_.superArcs)[superArcId].setUpNodeId(upNodeId);
  (*mt_data_.nodes)[upNodeId].addDownSuperArcId(superArcId);
}

void FTMTree_MT::delNode(idNode node) {
  Node *mainNode = getNode(node);

  if(mainNode->getNumberOfUpSuperArcs() == 0) {

    // Root: No Superarc
#ifndef TTK_ENABLE_KAMIKAZE
    if(mainNode->getNumberOfDownSuperArcs() != 1) {
      // Root with several childs: impossible /\ .
      cout << endl << "[FTMTree_MT]:delNode won't delete ";
      cout << mainNode->getVertexId() << " (root) with ";
      cout << static_cast<unsigned>(mainNode->getNumberOfDownSuperArcs())
           << " down ";
      cout << static_cast<unsigned>(mainNode->getNumberOfUpSuperArcs())
           << " up ";
      return;
    }
#endif

    idSuperArc downArc = mainNode->getDownSuperArcId(0);
    Node *downNode = getNode((*mt_data_.superArcs)[downArc].getDownNodeId());

    downNode->removeUpSuperArc(downArc);
    mainNode->clearDownSuperArcs();

  } else if(mainNode->getNumberOfDownSuperArcs() < 2) {
    // Have one up arc

    // We delete the upArc of this node,
    // if there is a down arc, we reattach it to the upNode

    idSuperArc upArc = mainNode->getUpSuperArcId(0);
    idNode upId = (*mt_data_.superArcs)[upArc].getUpNodeId();
    Node *upNode = getNode(upId);

    upNode->removeDownSuperArc(upArc);
    mainNode->clearUpSuperArcs();

    if(mainNode->getNumberOfDownSuperArcs()) {
      // Have one down arc

      // Reconnect
      idSuperArc downArc = mainNode->getDownSuperArcId(0);
      (*mt_data_.superArcs)[downArc].setUpNodeId(upId);
      upNode->addDownSuperArcId(downArc);
      mainNode->clearDownSuperArcs();

      // Segmentation
      (*mt_data_.superArcs)[downArc].concat((*mt_data_.superArcs)[upArc]);
    }
  }
#ifndef TTK_ENABLE_KAMIKAZE
  else
    cerr << "delete node with multiple childrens " << endl;
#endif
}

void FTMTree_MT::finalizeSegmentation(void) {
  for(auto &arc : *mt_data_.superArcs) {
    arc.createSegmentation(scalars_);
  }
}

tuple<SimplexId, SimplexId>
  FTMTree_MT::getBoundsFromVerts(const vector<SimplexId> &trunkVerts) const {
  SimplexId begin, stop;

  if(isST()) {
    begin = 0;
    stop = scalars_->offsets[trunkVerts[0]];
  } else {
    begin = scalars_->offsets[trunkVerts[0]];
    stop = scalars_->size;
  }

  return make_tuple(begin, stop);
}

Node *FTMTree_MT::getDownNode(const SuperArc *a) {
  return &((*mt_data_.nodes)[a->getDownNodeId()]);
}

idNode FTMTree_MT::getDownNodeId(const SuperArc *a) {
  return a->getDownNodeId();
}

Node *FTMTree_MT::getLowerNode(const SuperArc *a) {
  if(isST())
    return getUpNode(a);

  return getDownNode(a);
}

idNode FTMTree_MT::getLowerNodeId(const SuperArc *a) {
  if(isST())
    return getUpNodeId(a);

  return getDownNodeId(a);
}

Node *FTMTree_MT::getUpNode(const SuperArc *a) {
  return &((*mt_data_.nodes)[a->getUpNodeId()]);
}

idNode FTMTree_MT::getUpNodeId(const SuperArc *a) {
  return a->getUpNodeId();
}

Node *FTMTree_MT::getUpperNode(const SuperArc *a) {
  if(isST())
    return getDownNode(a);

  return getUpNode(a);
}

idNode FTMTree_MT::getUpperNodeId(const SuperArc *a) {
  if(isST())
    return getDownNodeId(a);

  return getUpNodeId(a);
}

idNode FTMTree_MT::getVertInRange(const vector<SimplexId> &range,
                                  const SimplexId v,
                                  const idNode last) const {
  idNode idRes = last;
  const idNode rangeSize = range.size();
  while(idRes + 1 < rangeSize && comp_.vertLower(range[idRes + 1], v)) {
    ++idRes;
  }
  return idRes;
}

idSuperArc FTMTree_MT::insertNode(Node *node, const bool segm) {
  // Normal insert : existing arc stay below inserted (JT example)
  //  *   - <- upNodeId
  //  | \ |   <- newSA
  //  |   * <- newNodeId
  //  |   |   <- currentSA
  //  - - -
  // already present
  if(isCorrespondingNode(node->getVertexId())) {
    Node *myNode = vertex2Node(node->getVertexId());
    // If it has been hidden / replaced we need to re-make it
    idSuperArc correspondingArcId = myNode->getUpSuperArcId(0);
    updateCorrespondingArc(myNode->getVertexId(), correspondingArcId);
  }

  idNode upNodeId, newNodeId;
  idSuperArc currentSA, newSA;
  SimplexId origin;

  // Create new node
  currentSA = getCorrespondingSuperArcId(node->getVertexId());
  upNodeId = (*mt_data_.superArcs)[currentSA].getUpNodeId();
  origin = (*mt_data_.nodes)[(*mt_data_.superArcs)[currentSA].getDownNodeId()]
             .getOrigin();
  newNodeId = makeNode(node, origin);

  // Connectivity
  // Insert only node inside the partition : created arc don t cross
  newSA = makeSuperArc(newNodeId, upNodeId);

  (*mt_data_.superArcs)[currentSA].setUpNodeId(newNodeId);
  (*mt_data_.nodes)[upNodeId].removeDownSuperArc(currentSA);
  (*mt_data_.nodes)[newNodeId].addDownSuperArcId(currentSA);

  // cut the vertex list at the node position and
  // give each arc its part.
  if(segm) {
    if(mt_data_.treeType == TreeType::Split) {
      (*mt_data_.superArcs)[newSA].concat(
        get<1>((*mt_data_.superArcs)[currentSA].splitBack(
          node->getVertexId(), scalars_)));
    } else {
      (*mt_data_.superArcs)[newSA].concat(
        get<1>((*mt_data_.superArcs)[currentSA].splitFront(
          node->getVertexId(), scalars_)));
    }
  }

  return newSA;
}

idNode FTMTree_MT::makeNode(SimplexId vertexId, SimplexId term) {
#ifndef TTK_ENABLE_KAMIKAZE
  if(vertexId < 0 || vertexId >= scalars_->size) {
    this->printMsg({{"make node, wrong vertex :", std::to_string(vertexId)},
                    {" on ", std::to_string(scalars_->size)}},
                   debug::Priority::ERROR);
    return -1;
  }
#endif

  if(isCorrespondingNode(vertexId)) {
    return getCorrespondingNodeId(vertexId);
  }

  idNode newNodeId = mt_data_.nodes->getNext();
  (*mt_data_.nodes)[newNodeId].setVertexId(vertexId);
  (*mt_data_.nodes)[newNodeId].setTerminaison(term);
  updateCorrespondingNode(vertexId, newNodeId);

  return newNodeId;
}

idNode FTMTree_MT::makeNode(const Node *const n, SimplexId ttkNotUsed(term)) {
  return makeNode(n->getVertexId());
}

idSuperArc FTMTree_MT::makeSuperArc(idNode downNodeId, idNode upNodeId)

{
  idSuperArc newSuperArcId = mt_data_.superArcs->getNext();
  (*mt_data_.superArcs)[newSuperArcId].setDownNodeId(downNodeId);
  (*mt_data_.superArcs)[newSuperArcId].setUpNodeId(upNodeId);

  (*mt_data_.nodes)[downNodeId].addUpSuperArcId(newSuperArcId);
  (*mt_data_.nodes)[upNodeId].addDownSuperArcId(newSuperArcId);

  return newSuperArcId;
}

void FTMTree_MT::move(FTMTree_MT *mt) {
  // we already have common data
  mt_data_.superArcs = mt->mt_data_.superArcs;
  mt->mt_data_.superArcs = nullptr;
  mt_data_.nodes = mt->mt_data_.nodes;
  mt->mt_data_.nodes = nullptr;
  mt_data_.leaves = mt->mt_data_.leaves;
  mt->mt_data_.leaves = nullptr;
  mt_data_.roots = mt->mt_data_.roots;
  mt->mt_data_.roots = nullptr;
  mt_data_.vert2tree = mt->mt_data_.vert2tree;
  mt->mt_data_.vert2tree = nullptr;
}

void FTMTree_MT::normalizeIds(void) {
  Timer normTime;
  sortLeaves(true);

  auto getNodeParentArcNb
    = [&](const idNode curNode, const bool goUp) -> idSuperArc {
    if(goUp) {
      return getNode(curNode)->getNumberOfUpSuperArcs();
    }

    return getNode(curNode)->getNumberOfDownSuperArcs();
  };

  auto getNodeParentArc
    = [&](const idNode curNode, const bool goUp, idSuperArc i) -> idSuperArc {
    if(goUp) {
      return getNode(curNode)->getUpSuperArcId(i);
    }

    return getNode(curNode)->getDownSuperArcId(i);
  };

  auto getArcParentNode
    = [&](const idSuperArc curArc, const bool goUp) -> idNode {
    if(goUp) {
      return getSuperArc(curArc)->getUpNodeId();
    }

    return getSuperArc(curArc)->getDownNodeId();
  };

  std::queue<tuple<idNode, bool>> q;
  std::stack<tuple<idNode, bool>> qr;
  for(const idNode n : *mt_data_.leaves) {
    bool goUp = isJT() || isST() || getNode(n)->getNumberOfUpSuperArcs();
    if(goUp)
      q.emplace(make_tuple(n, goUp));
    else
      qr.emplace(make_tuple(n, goUp));
  }

  while(!qr.empty()) {
    q.emplace(qr.top());
    qr.pop();
  }

  // Normalized id
  idSuperArc nIdMin = 0;
  idSuperArc nIdMax = getNumberOfSuperArcs() - 1;

  vector<bool> seenUp(getNumberOfSuperArcs(), false);
  vector<bool> seenDown(getNumberOfSuperArcs(), false);

  while(!q.empty()) {
    bool goUp;
    idNode curNodeId;
    tie(curNodeId, goUp) = q.front();
    q.pop();

    if(goUp)
      sortUpArcs(curNodeId);
    else
      sortDownArcs(curNodeId);

    // Assign arc above
    const idSuperArc nbArcParent = getNodeParentArcNb(curNodeId, goUp);
    for(idSuperArc pid = 0; pid < nbArcParent; pid++) {
      const idSuperArc currentArcId = getNodeParentArc(curNodeId, goUp, pid);
      if(goUp) {
        if(getSuperArc(currentArcId)->getNormalizedId() == nullSuperArc) {
          getSuperArc(currentArcId)->setNormalizeIds(nIdMin++);
        }
        if(!seenUp[currentArcId]) {
          q.emplace(make_tuple(getArcParentNode(currentArcId, goUp), goUp));
          seenUp[currentArcId] = true;
        }
      } else {
        if(getSuperArc(currentArcId)->getNormalizedId() == nullSuperArc) {
          getSuperArc(currentArcId)->setNormalizeIds(nIdMax--);
        }
        if(!seenDown[currentArcId]) {
          q.emplace(make_tuple(getArcParentNode(currentArcId, goUp), goUp));
          seenDown[currentArcId] = true;
        }
      }
    }
  }

#ifndef TTK_ENABLE_KAMIKAZE
  if(std::abs((long)nIdMax - (long)nIdMin) > 1) {
    this->printMsg({"error during normalize, tree compromized: ",
                    std::to_string(nIdMin), " ", std::to_string(nIdMax)},
                   debug::Priority::ERROR);
  }
#endif

  printTime(normTime, "normalize ids", 4);
}

idSuperArc FTMTree_MT::openSuperArc(idNode downNodeId) {
#ifndef TTK_ENABLE_KAMIKAZE
  if(downNodeId >= getNumberOfNodes()) {
    this->printErr("openSuperArc on a inexisting node !");
    return -2;
  }
#endif

  idSuperArc newSuperArcId = mt_data_.superArcs->getNext();
  (*mt_data_.superArcs)[newSuperArcId].setDownNodeId(downNodeId);
  (*mt_data_.nodes)[downNodeId].addUpSuperArcId(newSuperArcId);

  return newSuperArcId;
}

string FTMTree_MT::printArc(idSuperArc a) {
  const SuperArc *sa = getSuperArc(a);
  stringstream res;
  const SimplexId dv = getNode(sa->getDownNodeId())->getVertexId();
  const SimplexId uv = getNode(sa->getUpNodeId())->getVertexId();
  res << a;
  res << " : ";
  if(dv != nullVertex) {
    res << dv << " -- ";
  } else {
    res << "XX -- ";
  }
  if(uv != nullVertex) {
    res << uv;
  } else {
    res << "XX";
  }

  res.seekg(0, ios::end);
  while(res.tellg() < 25) {
    res << " ";
    res.seekg(0, ios::end);
  }
  res.seekg(0, ios::beg);

  res << "segm #" << sa->regionSize() << " / " << scalars_->size; // << " -> ";

  res.seekg(0, ios::end);

  while(res.tellg() < 45) {
    res << " ";
    res.seekg(0, ios::end);
  }
  res.seekg(0, ios::beg);

  res << sa->printReg();
  return res.str();
}

string FTMTree_MT::printNode(idNode n) {
  const Node *node = getNode(n);
  stringstream res;
  res << n;
  res << " : (";
  res << node->getVertexId() << ") \\ ";

  for(idSuperArc i = 0; i < node->getNumberOfDownSuperArcs(); ++i) {
    res << "+";
    res << node->getDownSuperArcId(i) << " ";
  }

  res << " / ";

  for(idSuperArc i = 0; i < node->getNumberOfUpSuperArcs(); ++i) {
    res << "+";
    res << node->getUpSuperArcId(i) << " ";
  }

  return res.str();
}

void FTMTree_MT::printParams(void) const {
  if(debugLevel_ > 1) {
    if(debugLevel_ > 2) {
      this->printMsg(ttk::debug::Separator::L1);
    }
    this->printMsg("number of threads : " + std::to_string(threadNumber_));
    if(debugLevel_ > 2) {
      this->printMsg("* debug lvl  : " + std::to_string(debugLevel_));
      string tt;
      if(params_->treeType == TreeType::Contour) {
        tt = "Contour";
      } else if(params_->treeType == TreeType::Join) {
        tt = "Join";
      } else if(params_->treeType == TreeType::Split) {
        tt = "Split";
      } else if(params_->treeType == TreeType::Join_Split) {
        tt = "Join + Split";
      }
      this->printMsg("* tree type  : " + tt);
      this->printMsg(ttk::debug::Separator::L1);
    }
  }
}

int FTMTree_MT::printTime(Timer &t,
                          const string &s,
                          const int debugLevel) const {

  if(this->debugLevel_ >= debugLevel) {
    stringstream st;

    for(int i = 3; i < debugLevel; i++)
      st << "-";
    st << s;

#ifdef TTK_ENABLE_FTM_TREE_PROCESS_SPEED
    if(nbScalars == -1) {
      nbScalars = scalars_->size;
    }
    int speed = nbScalars / t.getElapsedTime();
    st << " at " << speed << " vert/s";
#endif

    this->printMsg(st.str(), 1, t.getElapsedTime(), this->threadNumber_);
  }
  return 1;
}

void FTMTree_MT::printTree2() {
#ifdef TTK_ENABLE_OPENMP
#pragma omp critical
#endif
  {
    cout << "Nodes----------" << endl;
    for(idNode nid = 0; nid < getNumberOfNodes(); nid++) {
      cout << printNode(nid) << endl;
    }

    cout << "Arcs-----------" << endl;
    for(idSuperArc said = 0; said < getNumberOfSuperArcs(); ++said) {
      cout << printArc(said) << endl;
    }

    cout << "Leaves" << endl;
    for(const auto &l : *mt_data_.leaves)
      cout << " " << (*mt_data_.nodes)[l].getVertexId();
    cout << endl;

    cout << "Roots" << endl;
    for(const auto &r : *mt_data_.roots)
      cout << " " << (*mt_data_.nodes)[r].getVertexId();
    cout << endl;
  }
}

void FTMTree_MT::sortLeaves(const bool para) {
  auto indirect_sort = [&](const idNode a, const idNode b) {
    return comp_.vertLower(
      getNode(a)->getVertexId(), getNode(b)->getVertexId());
  };

  if(para) {
#ifdef __clang__
    std::sort(mt_data_.leaves->begin(), mt_data_.leaves->end(), indirect_sort);
#else
#ifndef _MSC_VER
#ifdef TTK_ENABLE_OPENMP
    __gnu_parallel::sort(
      mt_data_.leaves->begin(), mt_data_.leaves->end(), indirect_sort);
#else
    std::sort(mt_data_.leaves->begin(), mt_data_.leaves->end(), indirect_sort);
#endif
#else
    std::sort(mt_data_.leaves->begin(), mt_data_.leaves->end(), indirect_sort);
#endif
#endif
  } else {
    std::sort(mt_data_.leaves->begin(), mt_data_.leaves->end(), indirect_sort);
  }
}

vector<idNode> FTMTree_MT::sortedNodes(const bool para) {
  vector<idNode> sortedNodes(mt_data_.nodes->size());
  std::iota(sortedNodes.begin(), sortedNodes.end(), 0);

  auto indirect_sort = [&](const idNode a, const idNode b) {
    return comp_.vertLower(
      getNode(a)->getVertexId(), getNode(b)->getVertexId());
  };

  if(para) {
#ifdef __clang__
    std::sort(sortedNodes.begin(), sortedNodes.end(), indirect_sort);
#else
#ifndef _MSC_VER
#ifdef TTK_ENABLE_OPENMP
    __gnu_parallel::sort(sortedNodes.begin(), sortedNodes.end(), indirect_sort);
#else
    std::sort(sortedNodes.begin(), sortedNodes.end(), indirect_sort);
#endif
#else
    std::sort(sortedNodes.begin(), sortedNodes.end(), indirect_sort);
#endif
#endif
  } else {
#ifdef TTK_ENABLE_OPENMP
#pragma omp single
#endif
    { std::sort(sortedNodes.begin(), sortedNodes.end(), indirect_sort); }
  }

  return sortedNodes;
}

SimplexId FTMTree_MT::trunkCTSegmentation(const vector<SimplexId> &trunkVerts,
                                          const SimplexId begin,
                                          const SimplexId stop) {
  const int nbTasksThreads = 40;
  const auto sizeBackBone = abs(stop - begin);
  const auto chunkSize = getChunkSize(sizeBackBone, nbTasksThreads);
  const auto chunkNb = getChunkCount(sizeBackBone, nbTasksThreads);
  // si pas efficace vecteur de la taille de node ici a la place de acc
  idNode lastVertInRange = 0;
  mt_data_.trunkSegments->resize(getNumberOfSuperArcs());
  for(SimplexId chunkId = 0; chunkId < chunkNb; ++chunkId) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp task firstprivate(chunkId, lastVertInRange) shared(trunkVerts) \
  OPTIONAL_PRIORITY(isPrior())
#endif
    {
      vector<SimplexId> regularList;
      if(params_->segm) {
        regularList.reserve(25);
      }
      const SimplexId lowerBound = begin + chunkId * chunkSize;
      const SimplexId upperBound
        = min(stop, (begin + (chunkId + 1) * chunkSize));
      if(lowerBound != upperBound) {
        const SimplexId pos = isST() ? upperBound - 1 : lowerBound;
        lastVertInRange
          = getVertInRange(trunkVerts, scalars_->sortedVertices[pos], 0);
      }
      for(SimplexId v = lowerBound; v < upperBound; ++v) {
        const SimplexId s
          = isST() ? scalars_->sortedVertices[lowerBound + upperBound - 1 - v]
                   : scalars_->sortedVertices[v];
        if(isCorrespondingNull(s)) {
          const idNode oldVertInRange = lastVertInRange;
          lastVertInRange = getVertInRange(trunkVerts, s, lastVertInRange);
          const idSuperArc thisArc = upArcFromVert(trunkVerts[lastVertInRange]);
          updateCorrespondingArc(s, thisArc);

          if(params_->segm) {
            if(oldVertInRange == lastVertInRange) {
              regularList.emplace_back(s);
            } else {
              // accumulated to have only one atomic update when needed
              const idSuperArc oldArc
                = upArcFromVert(trunkVerts[oldVertInRange]);
              if(regularList.size()) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp critical
#endif
                {
                  (*mt_data_.trunkSegments)[oldArc].emplace_back(regularList);
                  regularList.clear();
                }
              }
              // hand.vtu, sequential: 28554
              regularList.emplace_back(s);
            }
          }
        }
      }
      // force increment last arc
      const idNode baseNode
        = getCorrespondingNodeId(trunkVerts[lastVertInRange]);
      const idSuperArc upArc = getNode(baseNode)->getUpSuperArcId(0);
      if(regularList.size()) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp critical
#endif
        {
          (*mt_data_.trunkSegments)[upArc].emplace_back(regularList);
          regularList.clear();
        }
      }
    }
  }
#ifdef TTK_ENABLE_OPENMP
#pragma omp taskwait
#endif
  // count added
  SimplexId tot = 0;
#ifdef TTK_ENABLE_FTM_TREE_PROCESS_SPEED
  for(const auto &l : *mt_data_.trunkSegments) {
    SimplexId arcSize = 0;
    for(const auto &v : l) {
      arcSize += v.size();
    }
    tot += arcSize;
  }
#endif
  return tot;
}

SimplexId FTMTree_MT::trunkSegmentation(const vector<SimplexId> &trunkVerts,
                                        const SimplexId begin,
                                        const SimplexId stop) {
  // Assign missing vert to the good arc
  // and also add the corresponding number for
  // futur arc reserve
  const int nbTasksThreads = 40;
  const auto sizeBackBone = abs(stop - begin);
  const auto chunkSize = getChunkSize(sizeBackBone, nbTasksThreads);
  const auto chunkNb = getChunkCount(sizeBackBone, nbTasksThreads);
  // si pas efficace vecteur de la taille de node ici a la place de acc
  SimplexId tot = 0;
  for(SimplexId chunkId = 0; chunkId < chunkNb; ++chunkId) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp task firstprivate(chunkId) shared(trunkVerts, tot) \
  OPTIONAL_PRIORITY(isPrior())
#endif
    {
      idNode lastVertInRange = 0;
      SimplexId acc = 0;

      const SimplexId lowerBound = begin + chunkId * chunkSize;
      const SimplexId upperBound
        = min(stop, (begin + (chunkId + 1) * chunkSize));
      for(SimplexId v = lowerBound; v < upperBound; ++v) {
        const SimplexId s
          = isST() ? scalars_->sortedVertices[lowerBound + upperBound - 1 - v]
                   : scalars_->sortedVertices[v];
        if(isCorrespondingNull(s)) {
          const idNode oldVertInRange = lastVertInRange;
          lastVertInRange = getVertInRange(trunkVerts, s, lastVertInRange);
          const idSuperArc thisArc = upArcFromVert(trunkVerts[lastVertInRange]);
          updateCorrespondingArc(s, thisArc);

          if(params_->segm) {
            if(oldVertInRange == lastVertInRange) {
              ++acc;
            } else {
              // accumulated to have only one atomic update when needed
              const idSuperArc oldArc
                = upArcFromVert(trunkVerts[oldVertInRange]);
              getSuperArc(oldArc)->atomicIncVisited(acc);
#ifdef TTK_ENABLE_FTM_TREE_PROCESS_SPEED
#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic update
#endif
              tot += acc;
#endif
              acc = 1;
            }
          }
        }
      }
      // force increment last arc
      const idNode baseNode
        = getCorrespondingNodeId(trunkVerts[lastVertInRange]);
      const idSuperArc upArc = getNode(baseNode)->getUpSuperArcId(0);
      getSuperArc(upArc)->atomicIncVisited(acc);
#ifdef TTK_ENABLE_FTM_TREE_PROCESS_SPEED
#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic update
#endif
      tot += acc;
#endif
    } // end task
  }
#ifdef TTK_ENABLE_OPENMP
#pragma omp taskwait
#endif
  return tot;
}

ostream &ttk::ftm::operator<<(ostream &o, SuperArc const &a) {
  o << a.getDownNodeId() << " <>> " << a.getUpNodeId();
  return o;
}

ostream &ttk::ftm::operator<<(ostream &o, Node const &n) {
  o << n.getNumberOfDownSuperArcs() << " .-. " << n.getNumberOfUpSuperArcs();
  return o;
}
