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
#include <set>

#include <cstddef>

using ttk::SimplexId;

#ifdef __cpp_lib_hardware_interference_size
    using std::hardware_constructive_interference_size;
    using std::hardware_destructive_interference_size;
#else
    constexpr std::size_t hardware_constructive_interference_size
        = 2 * sizeof(std::max_align_t);
    constexpr std::size_t hardware_destructive_interference_size
        = 2 * sizeof(std::max_align_t);
#endif

namespace ttk {

  namespace omsc {
    /**
     * Basic concept of critical point representing its id and index of 
     * criticality and simplex dimension - default 0 = vertex
     */
    struct CriticalPoint {
      explicit CriticalPoint() = default;

      explicit CriticalPoint(const SimplexId id, const int index)
        : id_{id}, index_{index} {
      }

      explicit CriticalPoint(const SimplexId id, const int index, const int dim)
        : id_{id}, index_{index}, dim_{dim} {
      }

      SimplexId id_{-1};
      int index_{-1};
      int dim_{0};
    };

    struct Separatrix {
      // default :
      explicit Separatrix()
        : source_{}, destination_{}, geometry_{}, destIsVertex_{},
          onBoundary_{} {
      }

      // initialization with only the extrema :
      explicit Separatrix(const SimplexId &extremum,
                          const SimplexId segmentGeometry)
        : source_{extremum} {
          geometry_.push_back(segmentGeometry);
        }

      // initialization with one segment :
      explicit Separatrix(const bool isValid,
                          const SimplexId &saddle,
                          const SimplexId &extremum,
                          const SimplexId segmentGeometry)
        : source_{saddle}, destination_{extremum},
          destIsVertex_{destIsVertex_} {
        geometry_.push_back(segmentGeometry);
      }

      // initialization with multiple segments :
      explicit Separatrix(const bool isValid,
                          const SimplexId &saddle,
                          const SimplexId &extremum,
                          const std::vector<SimplexId> &geometry)
        : source_{saddle}, destination_{extremum}, geometry_{geometry},
          destIsVertex_{destIsVertex_} {
      }

      explicit Separatrix(const Separatrix &separatrix)
        : source_{separatrix.source_},
          destination_{separatrix.destination_},
          geometry_{separatrix.geometry_},
          destIsVertex_{separatrix.destIsVertex_} {
      }

      explicit Separatrix(Separatrix &&separatrix) noexcept
        : source_{separatrix.source_},
          destination_{separatrix.destination_},
          geometry_{std::move(separatrix.geometry_)},
          destIsVertex_{separatrix.destIsVertex_} {
      }

      Separatrix &operator=(Separatrix &&separatrix) noexcept {
        source_ = separatrix.source_;
        destination_ = separatrix.destination_;
        geometry_ = std::move(separatrix.geometry_);
        destIsVertex_ = separatrix.destIsVertex_;

        return *this;
      }

      /**
       * Source vertex of the separatrix.
       */
      SimplexId source_;

      /**
       * Destination vertex of the separatrix.
       */
      SimplexId destination_{-1};

      /**
       * Container of ids. Each id addresses a separate
       * container corresponding to a dense representation
       * of the geometry (i.e. separatricesGeometry).
       */
      std::vector<SimplexId> geometry_;

      /**
       * Flag indicating that saddle is a vertex.
       */
      bool destIsVertex_;

      /**
       * Flag indicating that the separatix (partially) uses the boundary
       */
      bool onBoundary_;

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
        success += triangulation->preconditionVertexEdges();
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

    inline int setOutputSeparatrices1(
      SimplexId *const separatrices1_numberOfPoints,
      std::vector<float> *const separatrices1_points,
      std::vector<char> *const separatrices1_points_smoothingMask,
      std::vector<char> *const separatrices1_points_cellDimensions,
      std::vector<SimplexId> *const separatrices1_points_cellIds,
      SimplexId *const separatrices1_numberOfCells,
      std::vector<SimplexId> *const separatrices1_cells_connectivity,
      std::vector<SimplexId> *const separatrices1_cells_sourceIds,
      std::vector<SimplexId> *const separatrices1_cells_destinationIds,
      std::vector<SimplexId> *const separatrices1_cells_separatrixIds,
      std::vector<char> *const separatrices1_cells_separatrixTypes,
      std::vector<SimplexId> *const s1_cells_separatrixFunctionMaximaId,
      std::vector<SimplexId> *const s1_cells_separatrixFunctionMinimaId,
      std::vector<char> *const separatrices1_cells_isOnBoundary) {
      outputSeparatrices1_numberOfPoints_ = separatrices1_numberOfPoints;
      outputSeparatrices1_points_ = separatrices1_points;
      outputSeparatrices1_points_smoothingMask_
        = separatrices1_points_smoothingMask;
      outputSeparatrices1_points_cellDimensions_
        = separatrices1_points_cellDimensions;
      outputSeparatrices1_points_cellIds_ = separatrices1_points_cellIds;
      outputSeparatrices1_numberOfCells_ = separatrices1_numberOfCells;
      outputSeparatrices1_cells_connectivity_
        = separatrices1_cells_connectivity;
      outputSeparatrices1_cells_sourceIds_ = separatrices1_cells_sourceIds;
      outputSeparatrices1_cells_destinationIds_
        = separatrices1_cells_destinationIds;
      outputSeparatrices1_cells_separatrixIds_
        = separatrices1_cells_separatrixIds;
      outputSeparatrices1_cells_separatrixTypes_
        = separatrices1_cells_separatrixTypes;
      outputS1_cells_separatrixFunctionMaximaId_
        = s1_cells_separatrixFunctionMaximaId;
      outputS1_cells_separatrixFunctionMinimaId_
        = s1_cells_separatrixFunctionMinimaId;
      outputSeparatrices1_cells_isOnBoundary_
        = separatrices1_cells_isOnBoundary;
      return 0;
    }

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
      std::vector<omsc::CriticalPoint> &criticalPoints,    
      std::vector<omsc::Separatrix> &separatrices,
      const SimplexId *const morseSmaleManifold,
      const triangulationType &triangulation) const;

    template <typename triangulationType>
    void separatrixInnerPath(
      const SimplexId separatrixId,
      const SimplexId startTriangleId,
      const SimplexId startEdgeId,
      const omsc::CriticalPoint critPoint,
      const bool wasSplitted,
      SimplexId &globalSepId,
      omsc::Separatrix **separatrices,
      omsc::CriticalPoint **saddles,
      const SimplexId *const morseSmaleManifold,
      const triangulationType &triangulation) const;

    template <typename triangulationType>
    void separatrixBorderPath(
      const SimplexId separatrixId,
      const SimplexId startVertexId,
      const SimplexId startEdgeId,
      const omsc::CriticalPoint critPoint,
      const bool wasSplitted,
      SimplexId &globalSepId,
      omsc::Separatrix **separatrices,
      omsc::CriticalPoint **saddles,
      const SimplexId *const morseSmaleManifold,
      const triangulationType &triangulation) const;

    template <typename triangulationType>
    bool getNextTriangle(
      const SimplexId currentTriangle,
      const SimplexId currentEdge,
      SimplexId &nextTriangle, 
      SimplexId *nextEdges,
      const triangulationType &triangulation) const;

    template <typename triangulationType>
    bool hasVertexSingleLabelNeighborhood(
      const SimplexId currentVertex,
      const SimplexId *const morseSmaleManifold,
      const triangulationType &triangulation) const;

    template <typename triangulationType>
    bool isVertexOnBorder(
      const SimplexId currentVertex,
      const triangulationType &triangulation) const;

    template <typename triangulationType>
    bool nextBorderEdge(
      const SimplexId currentVertex,
      const SimplexId currentEdge,
      SimplexId &nextVertex,
      SimplexId &nextEdge,
      const triangulationType &triangulation) const;

    template <typename triangulationType>
    int setSeparatrices1(
      const std::vector<omsc::Separatrix> &separatrices,
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

    // segmentation
    void *outputAscendingManifold_{};
    void *outputDescendingManifold_{};
    void *outputMorseSmaleManifold_{};

    // critical points
    std::vector<std::array<float, 3>> *outputCriticalPoints_points_{};
    std::vector<char> *outputCriticalPoints_points_cellDimensions_{};
    std::vector<SimplexId> *outputCriticalPoints_points_cellIds_{};

    //1-separatrices
    SimplexId *outputSeparatrices1_numberOfPoints_{};
    std::vector<float> *outputSeparatrices1_points_{};
    std::vector<char> *outputSeparatrices1_points_smoothingMask_{};
    std::vector<char> *outputSeparatrices1_points_cellDimensions_{};
    std::vector<SimplexId> *outputSeparatrices1_points_cellIds_{};
    SimplexId *outputSeparatrices1_numberOfCells_{};
    std::vector<SimplexId> *outputSeparatrices1_cells_connectivity_{};
    std::vector<SimplexId> *outputSeparatrices1_cells_sourceIds_{};
    std::vector<SimplexId> *outputSeparatrices1_cells_destinationIds_{};
    std::vector<SimplexId> *outputSeparatrices1_cells_separatrixIds_{};
    std::vector<char> *outputSeparatrices1_cells_separatrixTypes_{};
    std::vector<SimplexId> *outputS1_cells_separatrixFunctionMaximaId_{};
    std::vector<SimplexId> *outputS1_cells_separatrixFunctionMinimaId_{};
    std::vector<char> *outputSeparatrices1_cells_isOnBoundary_{};

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
  std::vector<omsc::Separatrix> separatrices1;
  {
    Timer tmp;

    const int dim = triangulation.getDimensionality();

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

    this->printMsg("!!!Dim = " + std::to_string(dim));

    if(dim == 2) {
      computeSeparatrices1_2D<triangulationType>(
                                criticalPoints, separatrices1,
                                morseSmaleManifold, triangulation);

      this->printMsg("Write 1-seps");
      setSeparatrices1<triangulationType>(separatrices1, triangulation);
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
      criticalPoints.push_back(omsc::CriticalPoint(it.first, dim));
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
      criticalPoints.push_back(omsc::CriticalPoint(it.first, 0));
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
  std::vector<omsc::CriticalPoint> &criticalPoints, 
  std::vector<omsc::Separatrix> &separatrixVector,
  const SimplexId *const morseSmaleManifold,
  const triangulationType &triangulation) const {

  // start a local timer for this subprocedure
  ttk::Timer localTimer;

  this->printMsg(ttk::debug::Separator::L1);
  // print the progress of the current subprocedure (currently 0%)
  this->printMsg("Computing 1-separatrices by split",
                 0, // progress form 0-1
                 0, // elapsed time so far
                 this->threadNumber_);//, ttk::debug::LineMode::REPLACE);

  const size_t maxSepNumber = criticalPoints.size() * 4;
  omsc::Separatrix **separatrices = new omsc::Separatrix*[maxSepNumber]();
  omsc::CriticalPoint **saddles = new omsc::CriticalPoint*[maxSepNumber]();

  SimplexId globalSepId = 0;

  this->printMsg("#CritPoints: " + std::to_string(criticalPoints.size()));

  // retrieve maxima and minima
  for(const auto& c: criticalPoints) {
    const bool singleLabel = hasVertexSingleLabelNeighborhood(
      c.id_, morseSmaleManifold, triangulation);

    const bool onBorder = isVertexOnBorder(c.id_, triangulation);

    if(!singleLabel && !onBorder) { // regular inner critical point
      this->printMsg("[0]REG: " + std::to_string(c.id_));
      SimplexId numStar = triangulation.getVertexStarNumber(c.id_);
      for(SimplexId starId = 0; starId < numStar; ++starId) {

        SimplexId triangleId, edgeId, vertexId0, vertexId1;
        triangulation.getVertexStar(c.id_, starId, triangleId);

        for(int e = 0; e < 3; e++) {
          triangulation.getTriangleEdge(triangleId, e, edgeId);
          triangulation.getEdgeVertex(edgeId, 0, vertexId0);
          triangulation.getEdgeVertex(edgeId, 1, vertexId1);

          if(vertexId0 != c.id_ && vertexId1 != c.id_) {
            if(morseSmaleManifold[vertexId0] != morseSmaleManifold[vertexId1]) {
            SimplexId newSepId = globalSepId++;
            auto newSep = new omsc::Separatrix(c.id_, edgeId);
            newSep->geometry_.push_back(edgeId);
            separatrices[newSepId] = newSep;

            separatrixInnerPath(newSepId, triangleId, edgeId,
              c, false, globalSepId, separatrices, saddles,
              morseSmaleManifold, triangulation);
            }
            break;
          }     
        }
      }
    }
    else if(!singleLabel && onBorder) { // regular border critical point
    this->printMsg("[1]ROB: " + std::to_string(c.id_));
      SimplexId numStar = triangulation.getVertexStarNumber(c.id_);
      for(SimplexId starId = 0; starId < numStar; ++starId) {

        SimplexId triangleId, edgeId, vertexId0, vertexId1;
        triangulation.getVertexStar(c.id_, starId, triangleId);

        for(int e = 0; e < 3; e++) {
          triangulation.getTriangleEdge(triangleId, e, edgeId);
          triangulation.getEdgeVertex(edgeId, 0, vertexId0);
          triangulation.getEdgeVertex(edgeId, 1, vertexId1);

          if(vertexId0 != c.id_ && vertexId1 != c.id_) { // inner path
            if(morseSmaleManifold[vertexId0] != morseSmaleManifold[vertexId1]) {
              SimplexId newSepId = globalSepId++;
              auto newSep = new omsc::Separatrix(c.id_, edgeId);
              newSep->geometry_.push_back(edgeId);
              separatrices[newSepId] = newSep;

              this->printMsg("-Inner Path " + std::to_string(newSepId) +
              "- T:" + std::to_string(triangleId) +
              " E: " + std::to_string(edgeId));
              separatrixInnerPath(newSepId, triangleId, edgeId,
                c, false, globalSepId, separatrices, saddles,
                morseSmaleManifold, triangulation);
            }
          }
          else { // border paths
            if(triangulation.getEdgeTriangleNumber(edgeId) < 2) {
              SimplexId newSepId0 = globalSepId++;
              auto newSep0 = new omsc::Separatrix(c.id_, edgeId);
              separatrices[newSepId0] = newSep0;

              this->printMsg("-Border Path " + std::to_string(newSepId0) +
            "- V:" + std::to_string(vertexId0) +
            " E: " + std::to_string(edgeId));
              separatrixBorderPath(newSepId0, vertexId0, edgeId,
                c, false, globalSepId, separatrices, saddles,
                morseSmaleManifold, triangulation);
              
              SimplexId newSepId1 = globalSepId++;
              auto newSep = new omsc::Separatrix(c.id_, edgeId);
              separatrices[newSepId1] = newSep;

              this->printMsg("-Border Path " + std::to_string(newSepId1) +
            "- V:" + std::to_string(vertexId1) +
            " E: " + std::to_string(edgeId));
              separatrixBorderPath(newSepId1, vertexId1, edgeId,
                c, false, globalSepId, separatrices, saddles,
                morseSmaleManifold, triangulation);
            }
          } 
        }
      }
    }
    else if(singleLabel && onBorder) { // single label on border
      this->printMsg("[2]SOB: " + std::to_string(c.id_));
      SimplexId numStar = triangulation.getVertexStarNumber(c.id_);
      for(SimplexId starId = 0; starId < numStar; ++starId) {

        SimplexId triangleId, edgeId, vertexId0, vertexId1;
        triangulation.getVertexStar(c.id_, starId, triangleId);
       
        for(int e = 0; e < 3; e++) {
          triangulation.getTriangleEdge(triangleId, e, edgeId);
          
          if(triangulation.getEdgeTriangleNumber(edgeId) < 2) {
            triangulation.getEdgeVertex(edgeId, 0, vertexId0);
            triangulation.getEdgeVertex(edgeId, 1, vertexId1);

            SimplexId newSepId = globalSepId++;
              auto newSep = new omsc::Separatrix(c.id_, edgeId);
            separatrices[newSepId] = newSep;

            separatrixBorderPath(newSepId, c.id_, edgeId,
              c, false, globalSepId, separatrices, saddles,
              morseSmaleManifold, triangulation);
            
            break;
          }
        }
      }
    }
    else { //single label not on border
      this->printMsg("[3]SNB: " + std::to_string(c.id_));
    }
  }

  this->printMsg("write existing separatrices");


  // write existing separatrices
  for(int s = 0; s < globalSepId; ++s) {
    if(separatrices[s]) {
      separatrixVector.push_back(*separatrices[s]);
    }
  }

  this->printMsg("get unique critical points");

  // get unique critical points
  std::set<std::tuple<SimplexId, int>> critPointSet;
  {
    for(int c = 0; c < maxSepNumber; ++c) {
      if(saddles[c]) {
        critPointSet.insert(
          std::make_tuple(saddles[c]->id_, saddles[c]->dim_));
      }
    }
  }

  this->printMsg("write critical points");

  // write critical points
  std::array<float, 3UL> pos;
  for (auto it: critPointSet) {
    this->printMsg("ID: " + std::to_string(std::get<0>(it)) +
     "Dim: " + std::to_string(std::get<1>(it)));
    criticalPoints.push_back(
      omsc::CriticalPoint(std::get<0>(it), 1, std::get<1>(it)));

    if(std::get<1>(it) == 0) { // vertex
      triangulation.getVertexPoint(std::get<0>(it), pos[0], pos[1], pos[2]);
    } else { // edge
      triangulation.getEdgeIncenter(std::get<0>(it), pos.data());
    }
    
    outputCriticalPoints_points_->push_back(pos);
    outputCriticalPoints_points_cellDimensions_->push_back(1);
    outputCriticalPoints_points_cellIds_->push_back(std::get<0>(it));
  }
  this->printMsg("Test2");
  return 0;
}

template <typename triangulationType>
void ttk::morseSmaleComplexOrder::separatrixInnerPath(
  const SimplexId separatrixId,
  const SimplexId startTriangleId,
  const SimplexId startEdgeId,
  const omsc::CriticalPoint critPoint,
  const bool wasSplitted,
  SimplexId &globalSepId,
  omsc::Separatrix **separatrices,
  omsc::CriticalPoint **saddles,
  const SimplexId *const morseSmaleManifold,
  const triangulationType &triangulation) const {
  
  SimplexId currTriangleId = startTriangleId;
  SimplexId currEdgeId = startEdgeId;

  while(true) {

    SimplexId nextTriangleId, nextEdgeIds[2];
    bool foundNextTriangle = this->getNextTriangle(
      currTriangleId, currEdgeId, nextTriangleId, nextEdgeIds, triangulation);

    if(!foundNextTriangle) { // on border of geometry
      this->printMsg("--geometry border");
      if(critPoint.index_ == 0) { // minima - split path
        SimplexId verts[2];
        triangulation.getEdgeVertex(currEdgeId, 0, verts[0]);
        triangulation.getEdgeVertex(currEdgeId, 1, verts[1]);

        SimplexId newSepId = globalSepId++;
        auto newSep = new omsc::Separatrix(*separatrices[separatrixId]);
        separatrices[newSepId] = newSep;

        separatrixBorderPath(separatrixId, verts[0], currEdgeId, critPoint,
          true, globalSepId, separatrices, saddles,
          morseSmaleManifold, triangulation);
        
        separatrixBorderPath(newSepId, verts[1], currEdgeId, critPoint,
          true, globalSepId, separatrices, saddles,
          morseSmaleManifold, triangulation);
      }
      else { // maxima - create saddle
        separatrices[separatrixId]->geometry_.pop_back();
        separatrices[separatrixId]->destination_ = currEdgeId;
        separatrices[separatrixId]->destIsVertex_ = false;

        auto newSaddle = new omsc::CriticalPoint(currEdgeId, 1, 1);
        saddles[separatrixId] = newSaddle;
      }
      return;
    }

    SimplexId vertexId00, vertexId01, vertexId10, vertexId11;
    triangulation.getEdgeVertex(nextEdgeIds[0], 0, vertexId00);
    triangulation.getEdgeVertex(nextEdgeIds[0], 1, vertexId01);
    triangulation.getEdgeVertex(nextEdgeIds[1], 0, vertexId10);
    triangulation.getEdgeVertex(nextEdgeIds[1], 1, vertexId11);

    bool isSplittingMSC0 =
      morseSmaleManifold[vertexId00] == morseSmaleManifold[vertexId01]
      ? 0 : 1;
    bool isSplittingMSC1 =
      morseSmaleManifold[vertexId10] == morseSmaleManifold[vertexId11]
      ? 0 : 1;

    if(isSplittingMSC0 != isSplittingMSC1) { // path is extended
      //this->printMsg("--extend path");
      SimplexId nextEdgeId =
        isSplittingMSC0 ? nextEdgeIds[0] : nextEdgeIds[1];
      separatrices[separatrixId]->geometry_.push_back(nextEdgeId);
      currTriangleId = nextTriangleId;
      currEdgeId = nextEdgeId;
      continue;
    } else { // path is splitting or saddle (either edge or vertex saddle)

      SimplexId nnTriId0, nnTriId1, nnEdgeIds0[2], nnEdgeIds1[2];
      bool foundNextTriangle0 = this->getNextTriangle(
        nextTriangleId, nextEdgeIds[0],
        nnTriId0, nnEdgeIds0, triangulation);
      bool foundNextTriangle1 = this->getNextTriangle(
        nextTriangleId, nextEdgeIds[1],
        nnTriId1, nnEdgeIds1, triangulation);

      std::set<SimplexId> triangleSet0;
      triangleSet0.insert(morseSmaleManifold[vertexId00]);
      triangleSet0.insert(morseSmaleManifold[vertexId01]);
      triangleSet0.insert(morseSmaleManifold[vertexId10]);
      triangleSet0.insert(morseSmaleManifold[vertexId11]);
      std::set<SimplexId> triangleSet1(triangleSet0);

      if(foundNextTriangle0) {
        SimplexId nnVertId00, nnVertId01, nnVertId02;
        triangulation.getTriangleVertex(nnTriId0, 0, nnVertId00);
        triangulation.getTriangleVertex(nnTriId0, 1, nnVertId01);
        triangulation.getTriangleVertex(nnTriId0, 2, nnVertId02);
        triangleSet0.insert(morseSmaleManifold[nnVertId00]);
        triangleSet0.insert(morseSmaleManifold[nnVertId01]);
        triangleSet0.insert(morseSmaleManifold[nnVertId02]);
      }

      if(foundNextTriangle1) {
        SimplexId nnVertId10, nnVertId11, nnVertId12;
        triangulation.getTriangleVertex(nnTriId1, 0, nnVertId10);
        triangulation.getTriangleVertex(nnTriId1, 1, nnVertId11);
        triangulation.getTriangleVertex(nnTriId1, 2, nnVertId12);
        triangleSet1.insert(morseSmaleManifold[nnVertId10]);
        triangleSet1.insert(morseSmaleManifold[nnVertId11]);
        triangleSet1.insert(morseSmaleManifold[nnVertId12]);
      }

      bool edge0IsSaddle = triangleSet0.size() > 3;
      bool edge1IsSaddle = triangleSet1.size() > 3;

      if(edge0IsSaddle && edge1IsSaddle) { // saddle at nextEdgeIds[0] and [1]
        separatrices[separatrixId]->destination_ = nextEdgeIds[0];
        separatrices[separatrixId]->destIsVertex_ = false;
        auto newSaddle0 = new omsc::CriticalPoint(nextEdgeIds[0], 1, 1);
        saddles[separatrixId] = newSaddle0;

        const SimplexId newSepId = globalSepId++;
        auto newSeperatrix = new omsc::Separatrix(*separatrices[separatrixId]);
        separatrices[newSepId] = newSeperatrix;
        separatrices[newSepId]->destination_ = nextEdgeIds[1];
        separatrices[newSepId]->destIsVertex_ = false;
        auto newSaddle1 = new omsc::CriticalPoint(nextEdgeIds[1], 1, 1);
        saddles[newSepId] = newSaddle1;
        this->printMsg("--double edge saddle");
      } else if(edge0IsSaddle) { // nextEdgeIds[0] is a saddle
        separatrices[separatrixId]->destination_ = nextEdgeIds[0];
        separatrices[separatrixId]->destIsVertex_ = false;

        auto newSaddle = new omsc::CriticalPoint(nextEdgeIds[0], 1, 1);
        saddles[separatrixId] = newSaddle;
        this->printMsg("--0 single edge saddle");
      } else if(edge1IsSaddle) { // nextEdgeIds[1] is a saddle
        separatrices[separatrixId]->destination_ = nextEdgeIds[1];
        separatrices[separatrixId]->destIsVertex_ = false;

        auto newSaddle = new omsc::CriticalPoint(nextEdgeIds[1], 1, 1);
        saddles[separatrixId] = newSaddle;
        this->printMsg("--1 single edge saddle");
      } else { // check for vertex saddle 
        bool isVertexSaddle = false;
        SimplexId vertexSaddleId;
        for(int triVert = 0; triVert < 3; ++triVert) {

          if(isVertexSaddle) {
            break;
          }    
          triangulation.getTriangleVertex(
            nextTriangleId, triVert, vertexSaddleId);

          SimplexId triVertNN = triangulation.getVertexNeighborNumber(
            vertexSaddleId);

          std::set<SimplexId> commonVertexSet;
          for(int v = 0; v < triVertNN; ++v) {
            SimplexId vertNeigh;
            triangulation.getVertexNeighbor(vertexSaddleId, v, vertNeigh);
            if(vertNeigh == critPoint.id_) {
              commonVertexSet.clear();
              break;
            }
            commonVertexSet.insert(morseSmaleManifold[vertNeigh]);
          }
          if(commonVertexSet.size() > 3) {
            isVertexSaddle = true;
            break;
          }
        }

        if(isVertexSaddle) { // vertexSaddleId is a saddle         
          separatrices[separatrixId]->destination_ = vertexSaddleId;
          separatrices[separatrixId]->destIsVertex_ = true;

          auto newSaddle = new omsc::CriticalPoint(vertexSaddleId, 1, 0);
          saddles[separatrixId] = newSaddle;
          this->printMsg("--vertex saddle");
        } else { // no saddle - split the separatrix
          this->printMsg("--split");
          separatrices[separatrixId]->geometry_.push_back(nextEdgeIds[0]);

          SimplexId newSepId = globalSepId++;
          auto newSep = new omsc::Separatrix(*separatrices[separatrixId]);
          newSep->geometry_.push_back(nextEdgeIds[0]);
          separatrices[newSepId] = newSep;

          separatrixInnerPath(newSepId, nextTriangleId, nextEdgeIds[0],
            critPoint, false, globalSepId, separatrices, saddles,
            morseSmaleManifold, triangulation);

          separatrixInnerPath(separatrixId, nextTriangleId, nextEdgeIds[1],
            critPoint, false, globalSepId, separatrices, saddles,
            morseSmaleManifold, triangulation);
        }
      }
    }
    return;
  }
}

template <typename triangulationType>
void ttk::morseSmaleComplexOrder::separatrixBorderPath(
  const SimplexId separatrixId,
  const SimplexId startVertexId,
  const SimplexId startEdgeId,
  const omsc::CriticalPoint critPoint,
  const bool wasSplitted,
  SimplexId &globalSepId,
  omsc::Separatrix **separatrices,
  omsc::CriticalPoint **saddles,
  const SimplexId *const morseSmaleManifold,
  const triangulationType &triangulation) const {

  SimplexId currentEdge = startEdgeId, currentVertex = startVertexId;
  SimplexId nextVertex, nextEdge;

  nextBorderEdge(
    currentVertex, currentEdge, nextVertex, nextEdge, triangulation);
  
  separatrices[separatrixId]->geometry_.push_back(nextEdge);
  currentVertex = nextVertex;
  currentEdge = nextEdge;
  nextBorderEdge(
    currentVertex, currentEdge, nextVertex, nextEdge, triangulation);

  while(
    morseSmaleManifold[currentVertex] == morseSmaleManifold[nextVertex]) {
    separatrices[separatrixId]->geometry_.push_back(nextEdge);
    currentVertex = nextVertex;
    currentEdge = nextEdge;
    nextBorderEdge(
      currentVertex, currentEdge, nextVertex, nextEdge, triangulation);
  }

  separatrices[separatrixId]->onBoundary_ = true;
  separatrices[separatrixId]->destination_ = nextEdge;
  separatrices[separatrixId]->destIsVertex_ = false;
}

template <typename triangulationType>
bool ttk::morseSmaleComplexOrder::getNextTriangle(
  const SimplexId currentTriangle,
  const SimplexId currentEdge,
  SimplexId &nextTriangle, 
  SimplexId *nextEdges,
  const triangulationType &triangulation) const {
  
  SimplexId numTriangles = triangulation.getEdgeTriangleNumber(currentEdge);

  if(numTriangles <= 1) {
    return false;
  }
  triangulation.getEdgeTriangle(currentEdge, 0, nextTriangle);

  if(nextTriangle == currentTriangle) {
    triangulation.getEdgeTriangle(currentEdge, 1, nextTriangle);
  }

  // get triangle edges excluding the current one
  int foundEdges = 0;
  for(int e = 0; e < 3; ++e) {
    SimplexId edgeId;
    triangulation.getTriangleEdge(nextTriangle, e, edgeId);

    if(edgeId != currentEdge) {
      nextEdges[foundEdges++] = edgeId;
    }
  }
  
  return true;
}

template <typename triangulationType>
bool ttk::morseSmaleComplexOrder::hasVertexSingleLabelNeighborhood(
  const SimplexId currentVertex,
  const SimplexId *const morseSmaleManifold,
  const triangulationType &triangulation) const {
  
  SimplexId mSM = morseSmaleManifold[currentVertex];
  const SimplexId neighborNumber =
    triangulation.getVertexNeighborNumber(currentVertex);

  for(SimplexId n = 0; n < neighborNumber; ++n) {
    SimplexId neighbor;
    triangulation.getVertexNeighbor(currentVertex, n, neighbor);

    if(morseSmaleManifold[neighbor] != mSM) {
      return false;
    }
  }
  return true;
}

template <typename triangulationType>
bool ttk::morseSmaleComplexOrder::isVertexOnBorder(
  const SimplexId currentVertex,
  const triangulationType &triangulation) const {

  const SimplexId edgeNumber =
  triangulation.getVertexEdgeNumber(currentVertex);
  for(SimplexId e = 0; e < edgeNumber; ++e) {
    SimplexId edge;
    triangulation.getVertexEdge(currentVertex, e, edge);
    
    if(triangulation.getEdgeTriangleNumber(edge) < 2) {
      return true;
    }
  }

  return false;
}

template <typename triangulationType>
bool ttk::morseSmaleComplexOrder::nextBorderEdge(
  const SimplexId currentVertex,
  const SimplexId currentEdge,
  SimplexId &nextVertex,
  SimplexId &nextEdge,
  const triangulationType &triangulation) const {

  SimplexId edgeNumber =
    triangulation.getVertexEdgeNumber(currentVertex);

  for(SimplexId e = 0; e < edgeNumber; ++e) {
    triangulation.getVertexEdge(currentVertex, e, nextEdge);

    SimplexId edgeTriangleNumber =
      triangulation.getEdgeTriangleNumber(nextEdge);
    
    if(currentEdge != nextEdge && edgeTriangleNumber < 2) {
      triangulation.getEdgeVertex(nextEdge, 0, nextVertex);
      if(nextVertex == currentVertex) {
        triangulation.getEdgeVertex(nextEdge, 1, nextVertex);
      }
      return true;
    }
  }
  return false;
}

template <typename triangulationType>
int ttk::morseSmaleComplexOrder::setSeparatrices1(
  const std::vector<omsc::Separatrix> &separatrices,
  const triangulationType &triangulation) const {

#ifndef TTK_ENABLE_KAMIKAZE
  if(outputSeparatrices1_numberOfPoints_ == nullptr) {
    this->printErr("1-separatrices pointer to numberOfPoints is null.");
    return -1;
  }
  if(outputSeparatrices1_points_ == nullptr) {
    this->printErr("1-separatrices pointer to points is null.");
    return -1;
  }
  if(outputSeparatrices1_numberOfCells_ == nullptr) {
    this->printErr("1-separatrices pointer to numberOfCells is null.");
    return -1;
  }
  if(outputSeparatrices1_cells_connectivity_ == nullptr) {
    this->printErr("1-separatrices pointer to cells is null.");
    return -1;
  }
  if(inputOrderField_ == nullptr) {
    this->printErr(
      " 1-separatrices pointer to the input scalar field is null.");
    return -1;
  }
#endif

  this->printMsg("T0");

  //auto separatrixFunctionMaxima = outputS1_cells_separatrixFunctionMaximaId_;
  //auto separatrixFunctionMinima = outputS1_cells_separatrixFunctionMinimaId_;  

  // total number of separatrices points
  auto npoints{static_cast<size_t>(*outputSeparatrices1_numberOfPoints_)};
  // list of valid geometryId to flatten loops
  std::vector<SimplexId> validGeomIds{};
  // number of separatrices
  size_t numSep = separatrices.size();
  size_t numCells = 0;

  this->printMsg("T1");

  // count total number of points and cells, flatten geometryId loops
  for(size_t i = 0; i < numSep; ++i) {
    const auto &sep = separatrices[i];
    npoints += sep.geometry_.size() + 2;
    numCells += sep.geometry_.size() + 1;
  }

  this->printMsg("T2");

  const int dimensionality = triangulation.getDimensionality();

  // resize arrays
  outputSeparatrices1_points_->resize(3 * npoints);
  auto &points = *outputSeparatrices1_points_;
  outputSeparatrices1_cells_connectivity_->resize(2 * numCells);
  auto &cellsConn = *outputSeparatrices1_cells_connectivity_;
  if(outputSeparatrices1_cells_sourceIds_ != nullptr)
    outputSeparatrices1_cells_sourceIds_->resize(numSep);
  if(outputSeparatrices1_cells_destinationIds_ != nullptr)
    outputSeparatrices1_cells_destinationIds_->resize(numSep);
  if(outputSeparatrices1_cells_separatrixIds_ != nullptr)
    outputSeparatrices1_cells_separatrixIds_->resize(numSep);
  if(outputSeparatrices1_cells_isOnBoundary_ != nullptr)
    outputSeparatrices1_cells_isOnBoundary_->resize(numSep);

  SimplexId currPId = 0, currCId = 0;
  this->printMsg("T3");
  for(size_t i = 0; i < numSep; ++i) {

    const auto &sep = separatrices[i];
    const auto &src = sep.source_;
    const auto &dst = sep.destination_;
    const bool destIsVertex = sep.destIsVertex_;
    const auto onBoundary = sep.onBoundary_;

    std::array<float, 3> pt{};
    triangulation.getVertexPoint(src, pt[0], pt[1], pt[2]);
    points[3 * currPId + 0] = pt[0];
    points[3 * currPId + 1] = pt[1];
    points[3 * currPId + 2] = pt[2];

    //if(outputSeparatrices1_points_cellDimensions_ != nullptr)
    //  (*outputSeparatrices1_points_cellDimensions_)[currPId] = 0;
    //if(outputSeparatrices1_points_cellIds_ != nullptr)
    //  (*outputSeparatrices1_points_cellIds_)[currPId] = src;

    currPId += 1;

    for(const auto geom : sep.geometry_) {
      triangulation.getEdgeIncenter(geom, pt.data());

      points[3 * currPId + 0] = pt[0];
      points[3 * currPId + 1] = pt[1];
      points[3 * currPId + 2] = pt[2];

      //if(outputSeparatrices1_points_cellDimensions_ != nullptr)
      //  (*outputSeparatrices1_points_cellDimensions_)[currPId] = 1;
      //if(outputSeparatrices1_points_cellIds_ != nullptr)
      //  (*outputSeparatrices1_points_cellIds_)[currPId] = geom;

      cellsConn[2 * currCId + 0] = currPId - 1;
      cellsConn[2 * currCId + 1] = currPId;

      currPId += 1;
      currCId += 1;
    }

    if(destIsVertex) {
      triangulation.getVertexPoint(dst, pt[0], pt[1], pt[2]);
    }
    else {
      triangulation.getEdgeIncenter(dst, pt.data());
    }
    
    points[3 * currPId + 0] = pt[0];
    points[3 * currPId + 1] = pt[1];
    points[3 * currPId + 2] = pt[2];

    //if(outputSeparatrices1_points_cellDimensions_ != nullptr)
    //  (*outputSeparatrices1_points_cellDimensions_)[currPId] = 1;
    //if(outputSeparatrices1_points_cellIds_ != nullptr)
    //  (*outputSeparatrices1_points_cellIds_)[currPId] = src;

    cellsConn[2 * currCId + 0] = currPId - 1;
    cellsConn[2 * currCId + 1] = currPId;

    currPId += 1;
    currCId += 1;

    if(outputSeparatrices1_cells_sourceIds_ != nullptr)
      (*outputSeparatrices1_cells_sourceIds_)[i] = src;
    if(outputSeparatrices1_cells_destinationIds_ != nullptr)
      (*outputSeparatrices1_cells_destinationIds_)[i] = dst;
    if(outputSeparatrices1_cells_isOnBoundary_ != nullptr)
      (*outputSeparatrices1_cells_isOnBoundary_)[i] = onBoundary;
  }

  this->printMsg("T4");

  // update pointers
  *outputSeparatrices1_numberOfPoints_ = npoints;
  *outputSeparatrices1_numberOfCells_ = numCells;

  return 0;
}