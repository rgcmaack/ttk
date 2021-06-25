/// \ingroup base
/// \class ttk::morseSmaleComplexOrder
/// \author Robin G. C. Maack <maack@rhrk.uni-kl.de>
/// \date June 2021.

#pragma once

#include <Debug.h>
#include <Triangulation.h>
#include <unordered_map>
#include <map>
#include <list>
#include <vtk-9.0/vtkDataArray.h>

namespace ttk {

  /**
   * Basic concept of critical point representing its vertexid and index of 
   * criticality
   */
  struct CriticalPoint {
    explicit CriticalPoint() = default;

    explicit CriticalPoint(const int index, const SimplexId id)
      : index_{index}, id_{id} {
    }

    int index_{-1};
    SimplexId id_{-1};
  };

  /**
   * The morseSmaleComplexOrder class provides methods to compute a
   * MSC variant that only uses the input field and its order field.
   */
  class morseSmaleComplexOrder : virtual public Debug {

  public:
    morseSmaleComplexOrder();

    int preconditionTriangulation(
      ttk::AbstractTriangulation *triangulation) const {
      return triangulation->preconditionVertexNeighbors();
    };

    /**
     * Set the input order field.
     */
    inline int setInputOrderField(void *const data) {
      inputOrderField_ = data;
      return 0;
    };

    /**
     * Set the data pointers to the output segmentation scalar fields.
     */
    inline int setOutputMorseComplexes(void *const ascendingManifold,
                                       void *const descendingManifold,
                                       void *const morseSmaleManifold) {
      outputAscendingManifold_ = ascendingManifold;
      outputDescendingManifold_ = descendingManifold;
      outputMorseSmaleManifold_ = morseSmaleManifold;
      return 0;
    };

    /**
     * Set the output critical points data pointers.
     */
    inline int setOutputCriticalPoints(
      std::vector<std::array<float, 3>> *const criticalPoints_points,
      std::vector<char> *const criticalPoints_points_cellDimensons,
      std::vector<SimplexId> *const criticalPoints_points_cellIds) {
      outputCriticalPoints_points_ = criticalPoints_points;
      outputCriticalPoints_points_cellDimensions_
        = criticalPoints_points_cellDimensons;
      outputCriticalPoints_points_cellIds_ = criticalPoints_points_cellIds;
      return 0;
    };

    template <typename dataType, typename triangulationType>
    int execute(const triangulationType &triangulation);

    template <class dataType, class triangulationType>
    int computeDescendingManifold(
      SimplexId *const descendingManifold,
      std::vector<CriticalPoint> &criticalPoints,
      SimplexId &numberOfMaxima,
      const triangulationType &triangulation) const;

    template <class dataType, class triangulationType>
    int computeAscendingManifold(
      SimplexId *const ascendingManifold,
      std::vector<CriticalPoint> &criticalPoints,
      SimplexId &numberOfMinima,
      const triangulationType &triangulation) const;

    template <typename triangulationType>
    int computeFinalSegmentation(
      const SimplexId numMaxima,
      const SimplexId numMinima,
      const SimplexId *const ascendingManifold,
      const SimplexId *const descendingManifold,
      SimplexId *const morseSmaleManifold,
      const triangulationType &triangulation) const;
      
  protected:
    bool ComputeAscendingSeparatrices1{false};
    bool ComputeDescendingSeparatrices1{false};
    bool ComputeSaddleConnectors{false};
    bool ComputeAscendingSeparatrices2{false};
    bool ComputeDescendingSeparatrices2{false};
    bool ReturnSaddleConnectors{false};

    // input
    const void *inputOrderField_{};

    // critical points
    std::vector<std::array<float, 3>> *outputCriticalPoints_points_{};
    std::vector<char> *outputCriticalPoints_points_cellDimensions_{};
    std::vector<SimplexId> *outputCriticalPoints_points_cellIds_{};

    // segmentation
    void *outputAscendingManifold_{};
    void *outputDescendingManifold_{};
    void *outputMorseSmaleManifold_{};

  }; // morseSmaleComplexOrder class

} // namespace ttk

template <typename dataType, typename triangulationType>
int ttk::morseSmaleComplexOrder::execute(const triangulationType &triangulation)
{

#ifndef TTK_ENABLE_KAMIKAZE
  if(!inputOrderField_) {
    this->printErr("Input scalar field pointer is null.");
    return -1;
  }
#endif

  Timer t;

  // nullptr_t is implicitly convertible and comparable to any pointer type
  // or pointer-to-member type.
  SimplexId *ascendingManifold
    = static_cast<SimplexId *>(outputAscendingManifold_);
  SimplexId *descendingManifold
    = static_cast<SimplexId *>(outputDescendingManifold_);
  SimplexId *morseSmaleManifold
    = static_cast<SimplexId *>(outputMorseSmaleManifold_);

  std::vector<CriticalPoint> criticalPoints;
  {
    Timer tmp;

    SimplexId numberOfMaxima{};
    SimplexId numberOfMinima{};

    if(ascendingManifold)
      computeAscendingManifold<dataType, triangulationType>(
                                ascendingManifold, criticalPoints,
                                numberOfMinima, triangulation);

    if(descendingManifold)
      computeDescendingManifold<dataType, triangulationType>(
                                descendingManifold, criticalPoints,
                                numberOfMaxima, triangulation);

    if(ascendingManifold and descendingManifold and morseSmaleManifold)
      computeFinalSegmentation<triangulationType>(
                                numberOfMinima, numberOfMaxima,
                                ascendingManifold, descendingManifold,
                                morseSmaleManifold, triangulation);

    if(ascendingManifold or descendingManifold) {
      this->printMsg("Segmentation computed", 1.0, tmp.getElapsedTime(),
                     this->threadNumber_);
    }
  }

  this->printMsg("Data-set ("
                   + std::to_string(triangulation.getNumberOfVertices())
                   + " points) processed",
                 1.0, t.getElapsedTime(), this->threadNumber_);

  return 0;
}

template <class dataType, class triangulationType>
int ttk::morseSmaleComplexOrder::computeDescendingManifold(
  SimplexId *const descendingManifold,
  std::vector<CriticalPoint> &criticalPoints,
  SimplexId &numberOfMaxima,
  const triangulationType &triangulation) const {

  const dataType *inputField = static_cast<const dataType *>(inputOrderField_);

  // start global timer
  ttk::Timer globalTimer;
  int iterations = 1;
  numberOfMaxima = 0;

  // print horizontal separator
  this->printMsg(ttk::debug::Separator::L1); // L1 is the '=' separator

  // print input parameters in table format
  this->printMsg({
    {"#Threads", std::to_string(this->threadNumber_)},
    {"#Vertices", std::to_string(triangulation.getNumberOfVertices())},
  });
  this->printMsg(ttk::debug::Separator::L1);

  // -----------------------------------------------------------------------
  // Compute MSC variant
  // -----------------------------------------------------------------------
  {
    // start a local timer for this subprocedure
    ttk::Timer localTimer;

    // print the progress of the current subprocedure (currently 0%)
    this->printMsg("Computing Descending Manifold",
                   0, // progress form 0-1
                   0, // elapsed time so far
                   this->threadNumber_);

    /* compute the Descending Maifold iterating over each vertex, searching
     * the biggest neighbor and compressing its path to its maximum */
    SimplexId nVertices = triangulation.getNumberOfVertices();
    std::unordered_map<SimplexId, SimplexId> maxima;

    // vertices that may still be compressed
    std::vector<SimplexId>* activeVertices
      = new std::vector<SimplexId>();
    SimplexId nActiveVertices;

#ifndef TTK_ENABLE_OPENMP
    // find maxima and intialize vector of not fully compressed vertices
    for(SimplexId i = 0; i < nVertices; i++) {
      SimplexId neighborId;
      SimplexId numNeighbors =
        triangulation.getVertexNeighborNumber(i);
      bool hasBiggerNeighbor = false;

      descendingManifold[i] = i;

      // check all neighbors
      for(SimplexId n = 0; n < numNeighbors; n++) {
        triangulation.getVertexNeighbor(i, n, neighborId);
        if(inputField[neighborId] > inputField[descendingManifold[i]]) {
          descendingManifold[i] = neighborId;
          hasBiggerNeighbor = true;
        }
      }
      if(hasBiggerNeighbor) {
        activeVertices->push_back(i);
      }
      else {
        maxima.insert( {i, numberOfMaxima++} );
      }
    }

    nActiveVertices = activeVertices->size();
    std::vector<SimplexId>* newActiveVert;

    // compress paths until no changes occur
    while(nActiveVertices > 0) {
      newActiveVert = new std::vector<SimplexId>();

      std::string msgit = "Iteration: " + std::to_string(iterations) +
        "(" + std::to_string(nActiveVertices) + "/" +
        std::to_string(nVertices) + ")";
      this->printMsg(msgit, 1.0 - ((double)nActiveVertices / nVertices),
        localTimer.getElapsedTime(), this->threadNumber_);
        
      for(SimplexId i = 0; i < nActiveVertices; i++) {
        SimplexId v = activeVertices->at(i);
        dataType &vo = descendingManifold[v];
        // compress path
        vo = descendingManifold[vo];

        if(vo != descendingManifold[vo]) {
          newActiveVert->push_back(v);
        }
      }

      iterations += 1;
      delete activeVertices;
      activeVertices = &*newActiveVert;

      nActiveVertices = activeVertices->size();
    }

    // make indices dense
    for(SimplexId i = 0; i < nVertices; i++) {
      descendingManifold[i] = maxima[descendingManifold[i]];
    }
#else // TTK_ENABLE_OPENMP
    #pragma omp parallel num_threads(this->threadNumber_) \
            shared(iterations, maxima, numberOfMaxima, \
            activeVertices, nActiveVertices)
    {
      std::vector<SimplexId> lActiveVertices; //active verticies per thread

      #pragma omp for schedule(dynamic, 8)
      // find the biggest neighbor for each vertex
      for(SimplexId i = 0; i < nVertices; i++) {
        SimplexId neighborId;
        SimplexId numNeighbors =
          triangulation.getVertexNeighborNumber(i);
        bool hasBiggerNeighbor = false;

        descendingManifold[i] = i;

        // check all neighbors
        for(SimplexId n = 0; n < numNeighbors; n++) {
          triangulation.getVertexNeighbor(i, n, neighborId);
          if(inputField[neighborId] > inputField[descendingManifold[i]]) {
            descendingManifold[i] = neighborId;
            hasBiggerNeighbor = true;
          }
        }

        if(hasBiggerNeighbor) {
            lActiveVertices.push_back(i);
        }
        else {
          #pragma omp critical
          maxima.insert( {i, numberOfMaxima} );

          #pragma omp atomic update
          numberOfMaxima += 1;
        }
      }

      #pragma omp single nowait
      {
        std::string msgm =
          "Computed " + std::to_string(numberOfMaxima) + " Maxima";
        this->printMsg(msgm,
          0.1, localTimer.getElapsedTime(), this->threadNumber_);
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

        #pragma omp single nowait
        {
          std::string msgit = "Iteration: " + std::to_string(iterations) +
            "(" + std::to_string(nActiveVertices) + "/" +
            std::to_string(nVertices) + ")";
          double prog = 0.9 - ((double)nActiveVertices / nVertices) * 0.8;
          this->printMsg(msgit, prog,
            localTimer.getElapsedTime(), this->threadNumber_);
        }

        /* std::list< std::tuple<SimplexId, SimplexId> > lChanges; */

        #pragma omp for schedule(dynamic)
        for(SimplexId i = 0; i < nActiveVertices; i++) {
          SimplexId v = activeVertices->at(i);
          SimplexId &vo = descendingManifold[v];

          /* // save changes
          lChanges.push_front(
            std::make_tuple(v, descendingManifold[descendingManifold[v]]));*/
          vo = descendingManifold[vo];

          // check if fully compressed
          if(vo != descendingManifold[vo]) {
            lActiveVertices.push_back(v);
          }
        }
        #pragma omp barrier

        /* // apply changes
        for(auto it = lChanges.begin(); it != lChanges.end(); ++it) {
          descendingManifold[std::get<0>(*it)] = std::get<1>(*it);
        }*/

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
          iterations += 1;
        }
      }

      #pragma omp single
      {
        this->printMsg("Compressed Paths",
          0.95, // progress
          localTimer.getElapsedTime(), this->threadNumber_);
      }

      // set critical point indices
      #pragma omp for schedule(dynamic, 8)
      for(SimplexId i = 0; i < nVertices; i++) {
        descendingManifold[i] = maxima.at(descendingManifold[i]);
      }
    }
#endif // TTK_ENABLE_OPENMP

    const int dim = triangulation.getDimensionality();

    float x, y, z;
    for (auto& it: maxima) {
      criticalPoints.push_back(CriticalPoint(dim, it.first));
      triangulation.getVertexPoint(it.first, x, y, z);
      outputCriticalPoints_points_->push_back({x, y, z});
      outputCriticalPoints_points_cellDimensions_->push_back(dim);
      outputCriticalPoints_points_cellIds_->push_back(it.first);
    }

    delete activeVertices;

    // print the progress of the current subprocedure with elapsed time
    this->printMsg("Computed Descending Manifold",
                   1, // progress
                   localTimer.getElapsedTime(), this->threadNumber_);
  }

  // ---------------------------------------------------------------------
  // print global performance
  // ---------------------------------------------------------------------
  {
    this->printMsg(ttk::debug::Separator::L2); // horizontal '-' separator

    const std::string maxMsg = "#Maxima: " + std::to_string(numberOfMaxima);
    this->printMsg(maxMsg, 1, globalTimer.getElapsedTime() );

    const std::string itMsg = "#Iterations: " + std::to_string(iterations);
    this->printMsg(itMsg, 1, globalTimer.getElapsedTime() );

    this->printMsg("Completed", 1, globalTimer.getElapsedTime() );
    this->printMsg(ttk::debug::Separator::L1); // horizontal '=' separator
  }

  return 1; // return success
}


template <class dataType, class triangulationType>
int ttk::morseSmaleComplexOrder::computeAscendingManifold(
  SimplexId *const ascendingManifold,
  std::vector<CriticalPoint> &criticalPoints,
  SimplexId &numberOfMinima,
  const triangulationType &triangulation) const {

  const dataType *inputField = static_cast<const dataType *>(inputOrderField_);
  
  // start global timer
  ttk::Timer globalTimer;
  int iterations = 1;

  // print horizontal separator
  this->printMsg(ttk::debug::Separator::L1); // L1 is the '=' separator

  // print input parameters in table format
  this->printMsg({
    {"#Threads", std::to_string(this->threadNumber_)},
    {"#Vertices", std::to_string(triangulation.getNumberOfVertices())},
  });
  this->printMsg(ttk::debug::Separator::L1);

  // -----------------------------------------------------------------------
  // Compute MSC variant
  // -----------------------------------------------------------------------
  {
    // start a local timer for this subprocedure
    ttk::Timer localTimer;

    // print the progress of the current subprocedure (currently 0%)
    this->printMsg("Computing Ascending Manifold",
                   0, // progress form 0-1
                   0, // elapsed time so far
                   this->threadNumber_);

    /* compute the Descending Maifold iterating over each vertex, searching
     * the biggest neighbor and compressing its path to its maximum */

    SimplexId nVertices = triangulation.getNumberOfVertices();
    std::unordered_map<SimplexId, SimplexId> minima;

    // vertices that may still be compressed
    std::vector<SimplexId>* activeVertices
      = new std::vector<SimplexId>();
    SimplexId nActiveVertices;

#ifndef TTK_ENABLE_OPENMP
    // find minima and intialize vector of not fully compressed vertices
    for(SimplexId i = 0; i < nVertices; i++) {
      SimplexId neighborId;
      SimplexId numNeighbors =
        triangulation.getVertexNeighborNumber(i);
      bool hasBiggerNeighbor = false;

      ascendingManifold[i] = i;

      // check all neighbors
      for(SimplexId n = 0; n < numNeighbors; n++) {
        triangulation.getVertexNeighbor(i, n, neighborId);
        if(inputField[neighborId] < inputField[ascendingManifold[i]]) {
          ascendingManifold[i] = neighborId;
          hasBiggerNeighbor = true;
        }
      }
      if(hasBiggerNeighbor) {
        activeVertices->push_back(i);
      }
      else {
        minima.insert( {i, numberOfMinima++} );
      }
    }

    nActiveVertices = activeVertices->size();
    std::vector<SimplexId>* newActiveVert;

    // compress paths until no changes occur
    while(nActiveVertices > 0) {
      newActiveVert = new std::vector<SimplexId>();

      std::string msgit = "Iteration: " + std::to_string(iterations) +
        "(" + std::to_string(nActiveVertices) + "/" +
        std::to_string(nVertices) + ")";
      this->printMsg(msgit, 1.0 - ((double)nActiveVertices / nVertices),
        localTimer.getElapsedTime(), this->threadNumber_);
        
      for(SimplexId i = 0; i < nActiveVertices; i++) {
        SimplexId v = activeVertices->at(i);
        dataType &vo = ascendingManifold[v];
        // compress path
        vo = ascendingManifold[vo];

        if(vo != ascendingManifold[vo]) {
          newActiveVert->push_back(v);
        }
      }

      iterations += 1;
      delete activeVertices;
      activeVertices = &*newActiveVert;

      nActiveVertices = activeVertices->size();
    }

    for(SimplexId i = 0; i < nVertices; i++) {
      ascendingManifold[i] = minima[ascendingManifold[i]];
    }
#else // TTK_ENABLE_OPENMP
    #pragma omp parallel num_threads(this->threadNumber_) \
            shared(iterations, minima, numberOfMinima, \
            activeVertices, nActiveVertices)
    {
      std::vector<SimplexId> lActiveVertices; //active verticies per thread

      #pragma omp for schedule(dynamic, 8)
      // find the biggest neighbor for each vertex
      for(SimplexId i = 0; i < nVertices; i++) {
        SimplexId neighborId;
        SimplexId numNeighbors =
          triangulation.getVertexNeighborNumber(i);
        bool hasBiggerNeighbor = false;

        ascendingManifold[i] = i;

        // check all neighbors
        for(SimplexId n = 0; n < numNeighbors; n++) {
          triangulation.getVertexNeighbor(i, n, neighborId);
          if(inputField[neighborId] < inputField[ascendingManifold[i]]) {
            ascendingManifold[i] = neighborId;
            hasBiggerNeighbor = true;
          }
        }

        if(hasBiggerNeighbor) {
            lActiveVertices.push_back(i);
        }
        else {
          #pragma omp critical
          minima.insert( {i, numberOfMinima} );

          #pragma omp atomic update
          numberOfMinima += 1;
        }
      }

      #pragma omp single nowait
      {
        std::string msgm =
          "Computed " + std::to_string(numberOfMinima) + " Minima";
        this->printMsg(msgm,
          0.1, localTimer.getElapsedTime(), this->threadNumber_);
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

        #pragma omp single nowait
        {
          std::string msgit = "Iteration: " + std::to_string(iterations) +
            "(" + std::to_string(nActiveVertices) + "/" +
            std::to_string(nVertices) + ")";
          double prog = 0.9 - ((double)nActiveVertices / nVertices) * 0.8;
          this->printMsg(msgit, prog,
            localTimer.getElapsedTime(), this->threadNumber_);
        }

        /* std::list< std::tuple<SimplexId, SimplexId> > lChanges; */

        #pragma omp for schedule(dynamic)
        for(SimplexId i = 0; i < nActiveVertices; i++) {
          SimplexId v = activeVertices->at(i);
          SimplexId &vo = ascendingManifold[v];

          /* // save changes
          lChanges.push_front(
            std::make_tuple(v, ascendingManifold[ascendingManifold[v]]));*/
          vo = ascendingManifold[vo];

          // check if fully compressed
          if(vo != ascendingManifold[vo]) {
            lActiveVertices.push_back(v);
          }
        }
        #pragma omp barrier

        /* // apply changes
        for(auto it = lChanges.begin(); it != lChanges.end(); ++it) {
          ascendingManifold[std::get<0>(*it)] = std::get<1>(*it);
        }*/

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
          iterations += 1;
        }
      }

      #pragma omp single
      {
        this->printMsg("Compressed Paths",
          0.95, // progress
          localTimer.getElapsedTime(), this->threadNumber_);
      }

      // set critical point indices
      #pragma omp for schedule(dynamic, 8)
      for(SimplexId i = 0; i < nVertices; i++) {
        ascendingManifold[i] = minima.at(ascendingManifold[i]);
      }
    }
#endif // TTK_ENABLE_OPENMP

    float x, y, z;
    for (auto& it: minima) {
      criticalPoints.push_back(CriticalPoint(0, it.first));
      triangulation.getVertexPoint(it.first, x, y, z);
      outputCriticalPoints_points_->push_back({x, y, z});
      outputCriticalPoints_points_cellDimensions_->push_back(0);
      outputCriticalPoints_points_cellIds_->push_back(it.first);
    }

    delete activeVertices;

    // print the progress of the current subprocedure with elapsed time
    this->printMsg("Computed Ascending Manifold",
                   1, // progress
                   localTimer.getElapsedTime(), this->threadNumber_);
  }

  // ---------------------------------------------------------------------
  // print global performance
  // ---------------------------------------------------------------------
  {
    this->printMsg(ttk::debug::Separator::L2); // horizontal '-' separator

    const std::string maxMsg = "#Minima: " + std::to_string(numberOfMinima);
    this->printMsg(maxMsg, 1, globalTimer.getElapsedTime() );

    const std::string itMsg = "#Iterations: " + std::to_string(iterations);
    this->printMsg(itMsg, 1, globalTimer.getElapsedTime() );

    this->printMsg("Completed", 1, globalTimer.getElapsedTime() );
    this->printMsg(ttk::debug::Separator::L1); // horizontal '=' separator
  }

  return 1; // return success
}

template <typename triangulationType>
int ttk::morseSmaleComplexOrder::computeFinalSegmentation(
  const SimplexId numberOfMaxima,
  const SimplexId numberOfMinima,
  const SimplexId *const ascendingManifold,
  const SimplexId *const descendingManifold,
  SimplexId *const morseSmaleManifold,
  const triangulationType &triangulation) const {
  
  const size_t nVerts = triangulation.getNumberOfVertices();

  // associate a unique "sparse region id" to each (ascending, descending) pair
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < nVerts; ++i) {
    const auto a = ascendingManifold[i];
    const auto d = descendingManifold[i];
    if(a == -1 || d == -1) {
      morseSmaleManifold[i] = -1;
    } else {
      morseSmaleManifold[i] = a * numberOfMaxima + d;
    }
  }

  // store the "sparse region ids" by copying the morseSmaleManifold output
  std::vector<SimplexId> sparseRegionIds(
    morseSmaleManifold, morseSmaleManifold + nVerts);

  // get unique "sparse region ids"
  PSORT(this->threadNumber_)(sparseRegionIds.begin(), sparseRegionIds.end());
  const auto last = std::unique(sparseRegionIds.begin(), sparseRegionIds.end());
  sparseRegionIds.erase(last, sparseRegionIds.end());

  // "sparse region id" -> "dense region id"
  std::map<SimplexId, size_t> sparseToDenseRegionId{};

  for(size_t i = 0; i < sparseRegionIds.size(); ++i) {
    sparseToDenseRegionId[sparseRegionIds[i]] = i;
  }

  // update region id on all vertices: "sparse id" -> "dense id"

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < nVerts; ++i) {
    morseSmaleManifold[i] = sparseToDenseRegionId[morseSmaleManifold[i]];
  }

  return 0;
}