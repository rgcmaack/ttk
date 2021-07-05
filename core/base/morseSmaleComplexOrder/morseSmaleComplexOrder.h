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
#include <stack>

#include <cstddef>

using ttk::SimplexId;

#ifdef __cpp_lib_hardware_interference_size
    using std::hardware_constructive_interference_size;
    using std::hardware_destructive_interference_size;
#else
    // 64 bytes on x86-64 │ L1_CACHE_BYTES │ L1_CACHE_SHIFT │ __cacheline_aligned │ ...
    constexpr std::size_t hardware_constructive_interference_size
        = 2 * sizeof(std::max_align_t);
    constexpr std::size_t hardware_destructive_interference_size
        = 2 * sizeof(std::max_align_t);
#endif

namespace ttk {

  namespace omsc {
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

    struct Separatrix {
      // default :
      explicit Separatrix()
        : isValid_{}, source_{}, destination_{}, geometry_{} {
      }

      // initialization with only the extrema :
      explicit Separatrix(const SimplexId &extremum,
                          const SimplexId segmentGeometry)
        : isValid_{true}, source_{extremum} {
          geometry_.push_back(segmentGeometry);
        }

      // initialization with one segment :
      explicit Separatrix(const bool isValid,
                          const SimplexId &saddle,
                          const SimplexId &extremum,
                          const SimplexId segmentGeometry)
        : isValid_{isValid}, source_{saddle}, destination_{extremum} {
        geometry_.push_back(segmentGeometry);
      }

      // initialization with multiple segments :
      explicit Separatrix(const bool isValid,
                          const SimplexId &saddle,
                          const SimplexId &extremum,
                          const std::vector<SimplexId> &geometry)
        : isValid_{isValid}, source_{saddle}, destination_{extremum},
          geometry_{geometry} {
      }

      explicit Separatrix(const Separatrix &separatrix)
        : isValid_{separatrix.isValid_}, source_{separatrix.source_},
          destination_{separatrix.destination_},
          geometry_{separatrix.geometry_} {
      }

      explicit Separatrix(Separatrix &&separatrix) noexcept
        : isValid_{separatrix.isValid_}, source_{separatrix.source_},
          destination_{separatrix.destination_},
          geometry_{std::move(separatrix.geometry_)} {
      }

      Separatrix &operator=(Separatrix &&separatrix) noexcept {
        isValid_ = separatrix.isValid_;
        source_ = separatrix.source_;
        destination_ = separatrix.destination_;
        geometry_ = std::move(separatrix.geometry_);

        return *this;
      }

      /**
       * Flag indicating if this separatrix can be processed.
       */
      bool isValid_;

      /**
       * Source vertex of the separatrix.
       */
      SimplexId source_;

      /**
       * Destination vertex of the separatrix.
       */
      SimplexId destination_;

      /**
       * Container of ids. Each id addresses a separate
       * container corresponding to a dense representation
       * of the geometry (i.e. separatricesGeometry).
       */
      std::vector<SimplexId> geometry_;
    };
  }

  /**
   * The morseSmaleComplexOrder class provides methods to compute a
   * MSC variant that only uses the input field and its order field.
   */
  class morseSmaleComplexOrder : virtual public Debug {

  public:
    morseSmaleComplexOrder();

    int preconditionTriangulation(
      ttk::AbstractTriangulation *triangulation) const {
        int success = triangulation->preconditionVertexNeighbors();
        success += triangulation->preconditionVertexStars();
        success += triangulation->preconditionEdgeTriangles();
        success += triangulation->preconditionTriangleEdges();
      return success;
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
      std::vector<omsc::CriticalPoint> &criticalPoints,
      SimplexId &numberOfMaxima,
      const triangulationType &triangulation) const;

    template <class dataType, class triangulationType>
    int computeAscendingManifold(
      SimplexId *const ascendingManifold,
      std::vector<omsc::CriticalPoint> &criticalPoints,
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

    template <typename triangulationType>
    int computeSeparatrices1_2D(
      const SimplexId numMaxima,
      const SimplexId numMinima,
      std::vector<omsc::CriticalPoint> &criticalPoints,    
      std::vector<omsc::Separatrix> &separatrices,
      const SimplexId *const morseSmaleManifold,
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

    bool db_omp = true;

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

  // print horizontal separator
  this->printMsg(ttk::debug::Separator::L1); // L1 is the '=' separator

  // print input parameters in table format
  this->printMsg({
    {"#Threads", std::to_string(this->threadNumber_)},
    {"#Vertices", std::to_string(triangulation.getNumberOfVertices())},
  });
  this->printMsg(ttk::debug::Separator::L1);

  SimplexId *ascendingManifold
    = static_cast<SimplexId *>(outputAscendingManifold_);
  SimplexId *descendingManifold
    = static_cast<SimplexId *>(outputDescendingManifold_);
  SimplexId *morseSmaleManifold
    = static_cast<SimplexId *>(outputMorseSmaleManifold_);

  std::vector<omsc::CriticalPoint> criticalPoints;
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

  this->printMsg( std::to_string(triangulation.getNumberOfVertices())
                   + " verticies processed",
                 1.0, t.getElapsedTime(), this->threadNumber_);

  return 0;
}

template <class dataType, class triangulationType>
int ttk::morseSmaleComplexOrder::computeDescendingManifold(
  SimplexId *const descendingManifold,
  std::vector<omsc::CriticalPoint> &criticalPoints,
  SimplexId &numberOfMaxima,
  const triangulationType &triangulation) const {

  const dataType *inputField = static_cast<const dataType *>(inputOrderField_);

  // start global timer
  ttk::Timer globalTimer;
  int iterations = 1;
  numberOfMaxima = 0;

  // -----------------------------------------------------------------------
  // Compute MSC variant
  // -----------------------------------------------------------------------
  {
    // start a local timer for this subprocedure
    ttk::Timer localTimer;

    this->printMsg(ttk::debug::Separator::L1);
    // print the progress of the current subprocedure (currently 0%)
    this->printMsg("Computing Descending Manifold",
                   0, // progress form 0-1
                   0, // elapsed time so far
                   this->threadNumber_, ttk::debug::LineMode::REPLACE);

    /* compute the Descending Maifold iterating over each vertex, searching
     * the biggest neighbor and compressing its path to its maximum */
    const SimplexId nVertices = triangulation.getNumberOfVertices();
    std::unordered_map<SimplexId, SimplexId> maxima;

    // vertices that may still be compressed
    std::vector<SimplexId>* activeVertices
      = new std::vector<SimplexId>();
    SimplexId nActiveVertices;

//#ifndef TTK_ENABLE_OPENMP
if(!db_omp) {
    // find maxima and intialize vector of not fully compressed vertices
    for(SimplexId i = nVertices - 1; i >= 0; i--) {
      const SimplexId numNeighbors =
        triangulation.getVertexNeighborNumber(i);

      SimplexId neighborId;
      bool hasBiggerNeighbor = false;

      SimplexId &dmi = descendingManifold[i];
      dmi = i;

      // check all neighbors
      for(SimplexId n = 0; n < numNeighbors; n++) {
        triangulation.getVertexNeighbor(i, n, neighborId);
        if(inputField[neighborId] > inputField[dmi]) {
          dmi = neighborId;
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
        localTimer.getElapsedTime(), this->threadNumber_,
        ttk::debug::LineMode::REPLACE);
        
      for(SimplexId i = 0; i < nActiveVertices; i++) {
        SimplexId v = activeVertices->at(i);
        SimplexId &vo = descendingManifold[v];

        // compress path
        vo = descendingManifold[vo];

        // check if not fully compressed
        if(vo != descendingManifold[vo]) {
          newActiveVert->push_back(v);
        }
      }

      iterations += 1;
      delete activeVertices;
      activeVertices = newActiveVert;

      nActiveVertices = activeVertices->size();
    }

    // make indices dense
    for(SimplexId i = 0; i < nVertices; i++) {
      descendingManifold[i] = maxima[descendingManifold[i]];
    }
}
//#else // TTK_ENABLE_OPENMP
else {
    const size_t numCacheLineEntries =
      hardware_destructive_interference_size / sizeof(SimplexId);

    #pragma omp parallel num_threads(this->threadNumber_) \
            shared(iterations, maxima, numberOfMaxima, \
            activeVertices, nActiveVertices)
    {
      std::vector<SimplexId> lActiveVertices; //active verticies per thread

      // find the biggest neighbor for each vertex
      #pragma omp for schedule(dynamic, numCacheLineEntries)
      for(SimplexId i = nVertices - 1; i >= 0; i--) {
        SimplexId neighborId;
        SimplexId numNeighbors =
          triangulation.getVertexNeighborNumber(i);
        bool hasBiggerNeighbor = false;
        SimplexId &dmi = descendingManifold[i];

        dmi = i;

        // check all neighbors
        for(SimplexId n = 0; n < numNeighbors; n++) {
          triangulation.getVertexNeighbor(i, n, neighborId);
          if(inputField[neighborId] > inputField[dmi]) {
            dmi = neighborId;
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
            localTimer.getElapsedTime(), this->threadNumber_,
            ttk::debug::LineMode::REPLACE);
        }

        /* std::list< std::tuple<SimplexId, SimplexId> > lChanges; */

        #pragma omp for schedule(guided)
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

      #pragma omp single nowait
      {
        this->printMsg("Compressed Paths",
          0.95, // progress
          localTimer.getElapsedTime(), this->threadNumber_,
          ttk::debug::LineMode::REPLACE);
      }

      // set critical point indices
      #pragma omp for schedule(dynamic, numCacheLineEntries)
      for(SimplexId i = 0; i < nVertices; i++) {
        descendingManifold[i] = maxima.at(descendingManifold[i]);
      }
    }
}
//#endif // TTK_ENABLE_OPENMP

    const int dim = triangulation.getDimensionality();

    float x, y, z;
    for (auto& it: maxima) {
      criticalPoints.push_back(omsc::CriticalPoint(dim, it.first));
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

    const std::string itMsg = "#Iterations: " + std::to_string(iterations - 1);
    this->printMsg(itMsg, 1, globalTimer.getElapsedTime() );

    this->printMsg("Completed", 1, globalTimer.getElapsedTime() );
    this->printMsg(ttk::debug::Separator::L1); // horizontal '=' separator
  }

  return 1; // return success
}


template <class dataType, class triangulationType>
int ttk::morseSmaleComplexOrder::computeAscendingManifold(
  SimplexId *const ascendingManifold,
  std::vector<omsc::CriticalPoint> &criticalPoints,
  SimplexId &numberOfMinima,
  const triangulationType &triangulation) const {

  const dataType *inputField = static_cast<const dataType *>(inputOrderField_);

  // start global timer
  ttk::Timer globalTimer;
  int iterations = 1;

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
                   this->threadNumber_, ttk::debug::LineMode::REPLACE);

    /* compute the Descending Maifold iterating over each vertex, searching
     * the biggest neighbor and compressing its path to its maximum */

    SimplexId nVertices = triangulation.getNumberOfVertices();
    std::unordered_map<SimplexId, SimplexId> minima;

    numberOfMinima = 0;

    // vertices that may still be compressed
    std::vector<SimplexId>* activeVertices
      = new std::vector<SimplexId>();
    SimplexId nActiveVertices;

//#ifndef TTK_ENABLE_OPENMP
if(!db_omp) {
    // find minima and intialize vector of not fully compressed vertices
    for(SimplexId i = 0; i < nVertices; i++) {
      const SimplexId numNeighbors =
        triangulation.getVertexNeighborNumber(i);

      SimplexId neighborId;
      bool hasSmallerNeighbor = false;

      SimplexId &ami = ascendingManifold[i];
      ami = i;

      // check all neighbors
      for(SimplexId n = 0; n < numNeighbors; n++) {
        triangulation.getVertexNeighbor(i, n, neighborId);
        if(inputField[neighborId] < inputField[ami]) {
          ami = neighborId;
          hasSmallerNeighbor = true;
        }
      }
      if(hasSmallerNeighbor) {
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
        localTimer.getElapsedTime(), this->threadNumber_,
        ttk::debug::LineMode::REPLACE);
        
      for(SimplexId i = 0; i < nActiveVertices; i++) {
        SimplexId v = activeVertices->at(i);
        SimplexId &vo = ascendingManifold[v];

        // compress path
        vo = ascendingManifold[vo];

        if(vo != ascendingManifold[vo]) {
          newActiveVert->push_back(v);
        }
      }

      iterations += 1;
      delete activeVertices;
      activeVertices = newActiveVert;

      nActiveVertices = activeVertices->size();
    }

    for(SimplexId i = 0; i < nVertices; i++) {
      ascendingManifold[i] = minima[ascendingManifold[i]];
    }
}
//#else // TTK_ENABLE_OPENMP
else {
    const size_t numCacheLineEntries =
      hardware_destructive_interference_size / sizeof(SimplexId);

    #pragma omp parallel num_threads(this->threadNumber_) \
            shared(iterations, minima, numberOfMinima, \
            activeVertices, nActiveVertices)
    {
      std::vector<SimplexId> lActiveVertices; //active verticies per thread

      #pragma omp for schedule(dynamic, numCacheLineEntries)
      // find the biggest neighbor for each vertex
      for(SimplexId i = 0; i < nVertices; i++) {
        SimplexId neighborId;
        SimplexId numNeighbors =
          triangulation.getVertexNeighborNumber(i);
        bool hasSmallerNeighbor = false;

        ascendingManifold[i] = i;

        // check all neighbors
        for(SimplexId n = 0; n < numNeighbors; n++) {
          triangulation.getVertexNeighbor(i, n, neighborId);
          if(inputField[neighborId] < inputField[ascendingManifold[i]]) {
            ascendingManifold[i] = neighborId;
            hasSmallerNeighbor = true;
          }
        }

        if(hasSmallerNeighbor) {
            lActiveVertices.push_back(i);
        }
        else {
          #pragma omp critical
          minima.insert( {i, numberOfMinima} );

          #pragma omp atomic update
          numberOfMinima += 1;
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

        #pragma omp single nowait
        {
          std::string msgit = "Iteration: " + std::to_string(iterations) +
            "(" + std::to_string(nActiveVertices) + "/" +
            std::to_string(nVertices) + ")";
          double prog = 0.9 - ((double)nActiveVertices / nVertices) * 0.8;
          this->printMsg(msgit, prog,
            localTimer.getElapsedTime(), this->threadNumber_,
            ttk::debug::LineMode::REPLACE);
        }

        /* std::list< std::tuple<SimplexId, SimplexId> > lChanges; */

        #pragma omp for schedule(guided)
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

      #pragma omp single nowait
      {
        this->printMsg("Compressed Paths",
          0.95, // progress
          localTimer.getElapsedTime(), this->threadNumber_,
          ttk::debug::LineMode::REPLACE);
      }

      // set critical point indices
      #pragma omp for schedule(dynamic, numCacheLineEntries)
      for(SimplexId i = 0; i < nVertices; i++) {
        ascendingManifold[i] = minima.at(ascendingManifold[i]);
      }
    }
}
//#endif // TTK_ENABLE_OPENMP

    float x, y, z;
    for (auto& it: minima) {
      criticalPoints.push_back(omsc::CriticalPoint(0, it.first));
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

template <typename triangulationType>
int ttk::morseSmaleComplexOrder::computeSeparatrices1_2D(
  const SimplexId numMaxima,
  const SimplexId numMinima,
  std::vector<omsc::CriticalPoint> &criticalPoints, 
  std::vector<omsc::Separatrix> &separatrices,
  const SimplexId *const morseSmaleManifold,
  const triangulationType &triangulation) const {
    
  std::vector<omsc::CriticalPoint> maxima;
  std::vector<omsc::CriticalPoint> minima;
  const int dim = triangulation.getDimensionality();

  // retrieve maxima and minima
  for(const auto& c: criticalPoints) {
    if(c.index_ == 0) {
      maxima.push_back(c);
    }
    else if(c.index_ == dim) {
      minima.push_back(c);
    }
  }

  SimplexId numSeparatrices = 0;

/*#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_) schedule(dynamic) \
        shared(numSeparatrices, separatrices)
#endif*/
  // get edges splitting msc labels on their vertices
  // and the containing triangle
  for(SimplexId i = 0; i < numMaxima; ++i) {
    SimplexId &maxRef = maxima[i];

    // separatrixId, triangleId, edgeId
    std::stack<std::tuple<int, SimplexId, SimplexId>> s;
  
    SimplexId numStar = triangulation.getVertexStarNumber(maxRef);
    for(SimplexId starId = 0; starId < numStar; ++starId) {

      SimplexId triangleId, edgeId, vertexId0, vertexId1;
      triangulation.getVertexStar(maxRef, starId, triangleId);

      // get the edge that does not contain the maximum
      for(int e = 0; e < 3; e++) {
        triangulation.getTriangleEdge(triangleId, e, edgeId);
        triangulation.getEdgeVertex(edgeId, 0, vertexId0);
        triangulation.getEdgeVertex(edgeId, 1, vertexId1);

        if(vertexId0 != maxRef && vertexId1 != maxRef)
          break;
      }

      if(morseSmaleManifold[vertexId0] != morseSmaleManifold[vertexId1]) {
        s.push(std::make_tuple(numSeparatrices++, triangleId, edgeId));
        separatrices.push_back(omsc::Separatrix(maxRef, edgeId));
      }
    }

    while(!s.empty()) {
      std::tuple<int, SimplexId, SimplexId> sepTup = s.top(); s.pop();
      int sepId = std::get<0>(sepTup);
      SimplexId currTriangleId = std::get<1>(sepTup);
      SimplexId currEdgeId = std::get<2>(sepTup);

      // get next triangle by searching the edge triangles
      SimplexId nextTriangleId;
      SimplexId numTriangles = triangulation.getEdgeTriangleNumber(currEdgeId);
      triangulation.getEdgeTriangle(currEdgeId, 0, nextTriangleId);
      if(numTriangles > 1 && nextTriangleId == currTriangleId) {
        triangulation.getEdgeTriangle(currEdgeId, 1, nextTriangleId);
      }
      else { // border of the geometry
      }

      // get triangle edges excluding the current one
      SimplexId nextEdgeIds[2];
      int foundEdges = 0;
      for(int e = 0; e < 3; ++e) {
        SimplexId edgeId;
        triangulation.getTriangleEdge(nextTriangleId, e, edgeId);

        if(edgeId != currEdgeId) {
          nextEdgeIds[foundEdges++] = edgeId;
        }
      }

      SimplexId vertexId00, vertexId01, vertexId10, vertexId11;
      triangulation.getEgdeVertex(nextEdgeIds[0], 0, vertexId00);
      triangulation.getEgdeVertex(nextEdgeIds[0], 1, vertexId01);
      triangulation.getEgdeVertex(nextEdgeIds[1], 0, vertexId10);
      triangulation.getEgdeVertex(nextEdgeIds[1], 1, vertexId11);

      bool isSplittingMSC0 =
        morseSmaleManifold[vertexId00] == morseSmaleManifold[vertexId01]
        ? 0 : 1;
      bool isSplittingMSC1 =
        morseSmaleManifold[vertexId10] == morseSmaleManifold[vertexId11]
        ? 0 : 1;
      
      if(isSplittingMSC0 != isSplittingMSC1) { // path is extended
        SimplexId nextEdgeId =
          isSplittingMSC0 ? nextEdgeIds[0] : nextEdgeIds[1];
        separatrices[sepId].geometry_.push_back(nextEdgeId);
        s.push(std::make_tuple(sepId, nextTriangleId, nextEdgeId));
      }
      else { // path is splitting
        
      }
    }

    
  }

  return 0;
  
}