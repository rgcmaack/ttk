/// \ingroup base
/// \class ttk::ScalarFieldCriticalPointsPL
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date June 2015.
///
/// \brief TTK processing package for the computation of critical points in PL
/// scalar fields defined on PL manifolds.
///
/// This class computes the list of critical points of the input scalar field
/// and classify them according to their type.
///
/// \param dataType Data type of the input scalar field (char, float,
/// etc.).
///
/// \b Related \b publication \n
/// "Critical points and curvature for embedded polyhedral surfaces" \n
/// Thomas Banchoff \n
/// American Mathematical Monthly, 1970.
///
///  Progressive Approach used by default
/// \b Related \b publication \n
/// "A Progressive Approach to Scalar Field Topology" \n
/// Jules Vidal, Pierre Guillou, Julien Tierny\n
/// IEEE Transactions on Visualization and Computer Graphics, 2021
///
/// \sa ttkScalarFieldCriticalPointsPL.cpp %for a usage example.

#pragma once

#include <map>

// base code includes
#include <ScalarFieldCriticalPoints.h>
#include <Triangulation.h>

const bool saddleCandidateLookup[28] = {
  false, // (1) 0,0,0 - 0
  false, // (-) 0,0,1 - 1
  false, // (-) 0,0,2 - 2 
  true,  // (2) 0,0,3 - 3
  false, // (-) 0,1,0 - 4
  false, // (-) 0,1,1 - 5 
  false, // (-) 0,1,2 - 6
  false, // (-) 0,1,3 - 7
  true,  // (2) 0,2,0 - 8
  false, // (-) 0,2,1 - 9
  true,  // (2) 0,2,2 -10
  true,  // (3) 0,2,3 -11
  false, // (-) 0,3,0 -12
  false, // (-) 0,3,1 -13
  false, // (-) 0,3,2 -14
  false, // (-) 0,3,3 -15
  true,  // (2) 1,0,0 -16
  true,  // (2) 1,0,1 -17
  false, // (-) 1,0,2 -18
  true,  // (3) 1,0,3 -19
  true,  // (2) 1,1,0 -20
  true,  // (2) 1,1,1 -21
  false, // (-) 1,1,2 -22
  true,  // (3) 1,1,3 -23
  true,  // (3) 1,2,0 -24
  true,  // (3) 1,2,1 -25
  true,  // (3) 1,2,2 -26
  true   // (4) 1,2,3 -27
};

namespace ttk {

  class ScalarFieldCriticalPointsPL : public virtual Debug {

  public:
    ScalarFieldCriticalPointsPL();

    /**
     * Execute the package.
     * \param offsets Pointer to order field on vertices
     * \param triangulation Triangulation
     * \return Returns 0 upon success, negative values otherwise.
     *
     * @pre For this function to behave correctly in the absence of
     * the VTK wrapper, ttk::preconditionOrderArray() needs to be
     * called to fill the @p offsets buffer prior to any
     * computation (the VTK wrapper already includes a mecanism to
     * automatically generate such a preconditioned buffer).
     * @see examples/c++/main.cpp for an example use.
     */
    template <class triangulationType = AbstractTriangulation>
    int execute(const SimplexId *const offsets,
                const triangulationType *triangulation);

    template <typename triangulationType>
    int computeExtrema(
      const SimplexId *const offsets,
      bool *const saddleCandidates,
      SimplexId *const ascManifold,
      SimplexId *const dscManifold,
      int &numberOfMinima,
      int &numberOfMaxima,
      const triangulationType &triangulation) const;

    template <typename triangulationType>
    int computeFinalSegmentation(
      const SimplexId numMaxima,
      const SimplexId *const ascendingManifold,
      const SimplexId *const descendingManifold,
      SimplexId *const morseSmaleManifold,
      const triangulationType &triangulation) const;

    template <typename triangulationType>
    int computeSaddleCandidates_3D(
      bool *const saddleCandidates,
      const SimplexId *const morseSmaleManifold,
      const triangulationType &triangulation) const;

    int findSaddlesFromCadidates(
      const bool *const saddleCandidates,
      const SimplexId *const orderArr,
      const AbstractTriangulation *triangulation) const;

    inline void setDomainDimension(const int &dimension) {
      dimension_ = dimension;
    }

    inline void
      setOutput(std::vector<std::pair<SimplexId, char>> *criticalPoints) {
      criticalPoints_ = criticalPoints;
    }

    inline void
      preconditionTriangulation(AbstractTriangulation *triangulation) {
      // pre-condition functions
      if(triangulation) {
        triangulation->preconditionVertexNeighbors();
        triangulation->preconditionVertexStars();
        triangulation->preconditionBoundaryVertices();
      }

      setDomainDimension(triangulation->getDimensionality());
      setVertexNumber(triangulation->getNumberOfVertices());
    }

    inline void setVertexLinkEdgeLists(
      const std::vector<std::vector<std::pair<SimplexId, SimplexId>>>
        *edgeList) {
      vertexLinkEdgeLists_ = edgeList;
    }

    /// Set the number of vertices in the scalar field.
    /// \param vertexNumber Number of vertices in the data-set.
    /// \return Returns 0 upon success, negative values otherwise.
    int setVertexNumber(const SimplexId &vertexNumber) {
      vertexNumber_ = vertexNumber;
      return 0;
    }

    void setNonManifold(const bool b) {
      forceNonManifoldCheck = b;
    }

    void displayStats();

  protected:
    int dimension_{};
    SimplexId vertexNumber_{};
    const std::vector<std::vector<std::pair<SimplexId, SimplexId>>>
      *vertexLinkEdgeLists_{};
    std::vector<std::pair<SimplexId, char>> *criticalPoints_{};

    bool forceNonManifoldCheck{false};
    bool db_omp{true};
  };
} // namespace ttk

// template functions
template <class triangulationType>
int ttk::ScalarFieldCriticalPointsPL::execute(
  const SimplexId *const offsets, const triangulationType *triangulation) {

  ttk::Timer globalTimer;
  printMsg(ttk::debug::Separator::L1);

  SimplexId numMaxima, numMinima;

  const SimplexId nV = triangulation->getNumberOfVertices();

  SimplexId *ascManifold = new SimplexId[nV];
  SimplexId *decManifold = new SimplexId[nV];
  SimplexId *mscManifold = new SimplexId[nV];
  bool *saddleCandidates = new bool[nV];

  computeExtrema(offsets, saddleCandidates, ascManifold, decManifold, numMinima, numMaxima, triangulation);
  computeFinalSegmentation(numMaxima, ascManifold, decManifold, mscManifold, triangulation);
  computeSaddleCandidates_3D(saddleCandidates, mscManifold, triangulation);
  findSaddlesFromCadidates(saddleCandidates, offsets, triangulation);

  this->printMsg("Processed " + std::to_string(nV) + " vertices",
                 0, // progress form 0-1
                 globalTimer.getElapsedTime(), // elapsed time so far
                 this->threadNumber_);
  printMsg(ttk::debug::Separator::L1);
  return 0;
}

template <typename triangulationType>
int ttk::ScalarFieldCriticalPointsPL::computeExtrema(
  const SimplexId *const orderArr,
  bool *const saddleCandidates,
  SimplexId *const ascManifold,
  SimplexId *const dscManifold,
  int &numberOfMinima,
  int &numberOfMaxima,
  const triangulationType &triangulation) const {

  ttk::Timer localTimer;
  numberOfMinima = 0;
  numberOfMaxima = 0;
  
  this->printMsg(ttk::debug::Separator::L1);
  // print the progress of the current subprocedure (currently 0%)
  this->printMsg("Computing Scalar Field Extrema",
                 0, // progress form 0-1
                 0, // elapsed time so far
                 this->threadNumber_, ttk::debug::LineMode::REPLACE);

  /* compute the Descending Maifold iterating over each vertex, searching
   * the biggest neighbor and compressing its path to its maximum */
  const SimplexId nVertices = triangulation->getNumberOfVertices();

  // vertices that may still be compressed
  std::vector<SimplexId>* activeVertices
    = new std::vector<SimplexId>();
  SimplexId nActiveVertices;
  std::map<SimplexId, int> minMap, maxMap;

//#ifndef TTK_ENABLE_OPENMP
if(!db_omp) {
  // find maxima and intialize vector of not fully compressed vertices
  for(SimplexId i = 0; i < nVertices; i++) {
    saddleCandidates[i] = false;
    SimplexId neighborId;
    const SimplexId numNeighbors =
      triangulation->getVertexNeighborNumber(i);      

    bool hasBiggerNeighbor = false;
    SimplexId &dmi = dscManifold[i];
    dmi = i;

    bool hasSmallerNeighbor = false;
    SimplexId &ami = ascManifold[i];
    ami = i;

    // check all neighbors
    for(SimplexId n = 0; n < numNeighbors; n++) {
      triangulation->getVertexNeighbor(i, n, neighborId);

      if(orderArr[neighborId] < orderArr[ami]) {
        ami = neighborId;
        hasSmallerNeighbor = true;
      } else if(orderArr[neighborId] > orderArr[dmi]) {
        dmi = neighborId;
        hasBiggerNeighbor = true;
      }
    }

    if(hasBiggerNeighbor || hasSmallerNeighbor) {
      activeVertices->push_back(i);
    } 

    if(!hasBiggerNeighbor) {
      criticalPoints_->push_back(
        std::make_pair(i, (char)(CriticalType::Local_maximum)));
      maxMap.insert({i, numberOfMaxima});
      numberOfMaxima++;
    }

    if(!hasSmallerNeighbor) {
      criticalPoints_->push_back(
        std::make_pair(i, (char)(CriticalType::Local_minimum)));
      minMap.insert({i, numberOfMinima});
      numberOfMinima++;
    }
  }

  nActiveVertices = activeVertices->size();
  std::vector<SimplexId>* newActiveVert;

  // compress paths until no changes occur
  while(nActiveVertices > 0) {
    newActiveVert = new std::vector<SimplexId>();
      
    for(SimplexId i = 0; i < nActiveVertices; i++) {
      SimplexId &v = (*activeVertices)[i];
      SimplexId &vDsc = dscManifold[v];
      SimplexId &vAsc = ascManifold[v];

      // compress paths
      vDsc = dscManifold[vDsc];
      vAsc = ascManifold[vAsc];

      // check if not fully compressed
      if(vDsc != dscManifold[vDsc] || vAsc != ascManifold[vAsc]) {
        newActiveVert->push_back(v);
      }
    }

    delete activeVertices;
    activeVertices = newActiveVert;

    nActiveVertices = activeVertices->size();
  }

  for(SimplexId i = 0; i < nVertices; i++) {
    dscManifold[i] = maxMap[dscManifold[i]];
    ascManifold[i] = minMap[ascManifold[i]];
  }
}
//#else // TTK_ENABLE_OPENMP
else {
  #pragma omp parallel num_threads(this->threadNumber_)
  {
    std::vector<SimplexId> lActiveVertices; //active verticies per thread

    // find the biggest neighbor for each vertex
    #pragma omp for schedule(static)
    for(SimplexId i = 0; i < nVertices; i++) {
      saddleCandidates[i] = false;
      SimplexId neighborId;
      SimplexId numNeighbors =
        triangulation->getVertexNeighborNumber(i);

      bool hasBiggerNeighbor = false;
      SimplexId &dmi = dscManifold[i];
      dmi = i;

      bool hasSmallerNeighbor = false;
      SimplexId &ami = ascManifold[i];
      ami = i;

      // check all neighbors
      for(SimplexId n = 0; n < numNeighbors; n++) {
        triangulation->getVertexNeighbor(i, n, neighborId);

        if(orderArr[neighborId] < orderArr[ami]) {
          ami = neighborId;
          hasSmallerNeighbor = true;
        } else if(orderArr[neighborId] > orderArr[dmi]) {
          dmi = neighborId;
          hasBiggerNeighbor = true;
        }
      }

      if(hasBiggerNeighbor || hasSmallerNeighbor) {
        lActiveVertices.push_back(i);
      } 

      if(!hasBiggerNeighbor || !hasSmallerNeighbor) {
        #pragma omp critical
        {
          if(!hasBiggerNeighbor) {
            criticalPoints_->push_back(
              std::make_pair(i, (char)(CriticalType::Local_maximum)));
            maxMap.insert({i, numberOfMaxima});
            numberOfMaxima++;
          }

          if(!hasSmallerNeighbor) {
            criticalPoints_->push_back(
              std::make_pair(i, (char)(CriticalType::Local_minimum)));
            minMap.insert({i, numberOfMinima});
            numberOfMinima++;
          }
        }
      }
    }

    #pragma omp critical
    activeVertices->insert(activeVertices->end(),
      lActiveVertices.begin(), lActiveVertices.end());

    lActiveVertices.clear();

    #pragma omp barrier

    #pragma omp single
    {
      nActiveVertices = activeVertices->size();
    }

    // compress paths until no changes occur
    while(nActiveVertices > 0) {
      #pragma omp barrier

      #pragma omp for schedule(static)
      for(SimplexId i = 0; i < nActiveVertices; i++) {
        SimplexId &v = (*activeVertices)[i];
        SimplexId &vDsc = dscManifold[v];
        SimplexId &vAsc = ascManifold[v];

        // compress paths
        vDsc = dscManifold[vDsc];
        vAsc = ascManifold[vAsc];

        // check if fully compressed
        if(vDsc != dscManifold[vDsc] || vAsc != ascManifold[vAsc]) {
          lActiveVertices.push_back(v);
        }
      }

      #pragma omp single
      {
        activeVertices->clear();
      }

      #pragma omp critical
      activeVertices->insert(activeVertices->end(),
        lActiveVertices.begin(), lActiveVertices.end());

      lActiveVertices.clear();

      #pragma omp barrier

      #pragma omp single
      {
        nActiveVertices = activeVertices->size();
      }
    }

    #pragma omp for schedule(static)
    for(SimplexId i = 0; i < nVertices; i++) {
      dscManifold[i] = maxMap[dscManifold[i]];
      ascManifold[i] = minMap[ascManifold[i]];
    }
  }
} //#endif // TTK_ENABLE_OPENMP

  delete activeVertices;
  this->printMsg("Computing Scalar Field Extrema", 0.95,
                 localTimer.getElapsedTime(), this->threadNumber_);

  return 1; // return success
}

template <typename triangulationType>
int ttk::ScalarFieldCriticalPointsPL::computeFinalSegmentation(
  const SimplexId numberOfMinima,
  const SimplexId *const ascendingManifold,
  const SimplexId *const descendingManifold,
  SimplexId *const morseSmaleManifold,
  const triangulationType &triangulation) const {
  ttk::Timer localTimer;

  this->printMsg("Computing MSC Manifold",
                 0, localTimer.getElapsedTime(), this->threadNumber_,
                 ttk::debug::LineMode::REPLACE);
  
  const size_t nVerts = triangulation->getNumberOfVertices();

//#ifndef TTK_ENABLE_OPENMP
if(!db_omp) {
  for(size_t i = 0; i < nVerts; ++i) {
    morseSmaleManifold[i] = ascendingManifold[i] * numberOfMinima +
      descendingManifold[i];
  }
//#else // TTK_ENABLE_OPENMP
} else {
  #pragma omp parallel for schedule(static) num_threads(threadNumber_)
  for(size_t i = 0; i < nVerts; ++i) {
    morseSmaleManifold[i] = ascendingManifold[i] * numberOfMinima +
      descendingManifold[i];
  }
} //#endif // TTK_ENABLE_OPENMP

  this->printMsg("Computed MSC Manifold",
                 1, localTimer.getElapsedTime(), this->threadNumber_);

  return 0;
}

template <typename triangulationType>
int ttk::ScalarFieldCriticalPointsPL::computeSaddleCandidates_3D(
  bool *const saddleCandidates,
  const SimplexId *const morseSmaleManifold,
  const triangulationType &triangulation) const {

  ttk::Timer localTimer;

  // print the progress of the current subprocedure (currently 0%)
  this->printMsg("Computing 2-Separatrices 3D",
                 0, // progress form 0-1
                 0, // elapsed time so far
                 this->threadNumber_, ttk::debug::LineMode::REPLACE);

  const SimplexId numTetra = triangulation->getNumberOfCells();

//#ifndef TTK_ENABLE_OPENMP
if(!db_omp) {
  for(SimplexId tet = 0; tet < numTetra; tet++) {
    SimplexId vertices[4];
    triangulation->getCellVertex(tet, 0, vertices[0]);
    triangulation->getCellVertex(tet, 1, vertices[1]);
    triangulation->getCellVertex(tet, 2, vertices[2]);
    triangulation->getCellVertex(tet, 3, vertices[3]);

    const SimplexId msm[4] = {
      morseSmaleManifold[vertices[0]], morseSmaleManifold[vertices[1]],
      morseSmaleManifold[vertices[2]], morseSmaleManifold[vertices[3]]};

    const unsigned char index1 = (msm[0] == msm[1]) ? 0x00 : 0x10; // 0 : 1
    const unsigned char index2 = (msm[0] == msm[2]) ? 0x00 :       // 0
                                 (msm[1] == msm[2]) ? 0x04 : 0x08; // 1 : 2
    const unsigned char index3 = (msm[0] == msm[3]) ? 0x00 :       // 0
                                 (msm[1] == msm[3]) ? 0x01 :       // 1
                                 (msm[2] == msm[3]) ? 0x02 : 0x03; // 2 : 3

    const unsigned char lookupIndex = index1 | index2 | index3;

    if(saddleCandidateLookup[lookupIndex]) {
      saddleCandidates[vertices[0]] = true;
      saddleCandidates[vertices[1]] = true;
      saddleCandidates[vertices[2]] = true;
      saddleCandidates[vertices[3]] = true;      
    }
  }
//#else // TTK_ENABLE_OPENMP
} else {
  #pragma omp parallel num_threads(this->threadNumber_)
  {
    #pragma omp for schedule(static)
    for(SimplexId tet = 0; tet < numTetra; tet++) {
      SimplexId vertices[4];
      triangulation->getCellVertex(tet, 0, vertices[0]);
      triangulation->getCellVertex(tet, 1, vertices[1]);
      triangulation->getCellVertex(tet, 2, vertices[2]);
      triangulation->getCellVertex(tet, 3, vertices[3]);

      const SimplexId msm[4] = {
        morseSmaleManifold[vertices[0]], morseSmaleManifold[vertices[1]],
        morseSmaleManifold[vertices[2]], morseSmaleManifold[vertices[3]]};

      const unsigned char index1 = (msm[0] == msm[1]) ? 0x00 : 0x10; // 0 : 1
      const unsigned char index2 = (msm[0] == msm[2]) ? 0x00 :       // 0
                                   (msm[1] == msm[2]) ? 0x04 : 0x08; // 1 : 2
      const unsigned char index3 = (msm[0] == msm[3]) ? 0x00 :       // 0
                                   (msm[1] == msm[3]) ? 0x01 :       // 1
                                   (msm[2] == msm[3]) ? 0x02 : 0x03; // 2 : 3

      const unsigned char lookupIndex = index1 | index2 | index3;

      if(saddleCandidateLookup[lookupIndex]) {
        saddleCandidates[vertices[0]] = true;
        saddleCandidates[vertices[1]] = true;
        saddleCandidates[vertices[2]] = true;
        saddleCandidates[vertices[3]] = true;
      }
    }
  }
} // #if TTK_OMPENMP
  this->printMsg("Computed 2-Separatrices 3D",
                 1, localTimer.getElapsedTime(), this->threadNumber_);

  return 0;
}

int ttk::ScalarFieldCriticalPointsPL::findSaddlesFromCadidates(
  const bool *const saddleCandidates,
  const SimplexId *const orderArr,
  const AbstractTriangulation *triangulation) const {

  ttk::Timer localTimer;

  this->printMsg("Computing for saddles",
                 0, localTimer.getElapsedTime(), this->threadNumber_,
                 ttk::debug::LineMode::REPLACE);

  const int dim = triangulation->getDimensionality();
  const SimplexId nVertices = triangulation->getNumberOfVertices();

  ScalarFieldCriticalPoints sfcp;
  sfcp.setDomainDimension(dim);
  
//#ifndef TTK_ENABLE_OPENMP
if(!db_omp) {
  for(SimplexId i = 0; i < nVertices; i++) {
    if(saddleCandidates[i]) {
      const CriticalType critC = (CriticalType)sfcp.getCriticalType(
        i, orderArr, triangulation);

      if(critC == CriticalType::Saddle1 || critC == CriticalType::Saddle2 ||
        critC == CriticalType::Degenerate) {
        criticalPoints_->push_back(std::make_pair(i, (char)(critC)));
      }
    }
  }
} else {
  #pragma omp parallel for schedule(static)
  for(SimplexId i = 0; i < nVertices; i++) {
    if(saddleCandidates[i]) {
      const CriticalType critC = (CriticalType)sfcp.getCriticalType(
        i, orderArr, triangulation);

      if(critC == CriticalType::Saddle1 || critC == CriticalType::Saddle2 ||
        critC == CriticalType::Degenerate) {
        #pragma omp critical
        criticalPoints_->push_back(std::make_pair(i, (char)(critC)));
      }
    }
  }
}

  this->printMsg("Computed saddles" ,
                 1, localTimer.getElapsedTime(), this->threadNumber_);
                 
  return 1;
}