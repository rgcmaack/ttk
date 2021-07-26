/// \ingroup base
/// \class ttk::MorseSmaleComplexQuasi
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
#include <unordered_set>
#include <bitset>

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

int tetraederLookup[28][12] = {
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 0,0,0
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 0,0,1
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 0,0,2
  { 2,  4,  5, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 0,0,3
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 0,1,0
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 0,1,1
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 0,1,2
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 0,1,3
  { 1,  3,  5, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 0,2,0
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 0,2,1
  { 1,  2,  4,  1,  3,  4, -1, -1, -1, -1, -1, -1}, // 0,2,2
  { 1,  3,  5,  2,  4,  5,  3,  4,  5,  1,  2,  5}, // 0,2,3
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 0,3,0
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 0,3,1
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 0,3,2
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 0,3,3
  { 0,  3,  4, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 1,0,0
  { 0,  2,  3,  2,  3,  5, -1, -1, -1, -1, -1, -1}, // 1,0,1
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 1,0,2
  { 0,  3,  4,  2,  4,  5,  0,  2,  4,  3,  4,  5}, // 1,0,3
  { 0,  1,  4,  1,  4,  5, -1, -1, -1, -1, -1, -1}, // 1,1,0
  { 0,  1,  2, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 1,1,1
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 1,1,2
  { 0,  1,  2,  2,  4,  5,  1,  2,  5,  0,  2,  4}, // 1,1,3
  { 0,  3,  4,  1,  3,  5,  0,  1,  3,  3,  4,  5}, // 1,2,0
  { 0,  1,  2,  1,  3,  5, -1, -1, -1, -1, -1, -1}, // 1,2,1
  { 0,  1,  2,  0,  3,  4,  0,  2,  4,  0,  1,  3}, // 1,2,2
  { 6,  6,  6,  6,  6,  6, -1, -1, -1, -1, -1, -1}  // 1,2,3
};

namespace ttk {

  namespace mscq {
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
        : source_{}, destination_{}, geometry_{}, onBoundary_{} {
      }

      // initialization with only the extrema :
      explicit Separatrix(const SimplexId segmentGeometry) {
          geometry_.push_back(segmentGeometry);
        }

      // initialization with one segment :
      explicit Separatrix(const bool isValid,
                          const SimplexId &saddle,
                          const SimplexId &extremum,
                          const SimplexId segmentGeometry)
        : source_{saddle}, destination_{extremum} {
        geometry_.push_back(segmentGeometry);
      }

      // initialization with multiple segments :
      explicit Separatrix(const bool isValid,
                          const SimplexId &saddle,
                          const SimplexId &extremum,
                          const std::vector<SimplexId> &geometry)
        : source_{saddle}, destination_{extremum}, geometry_{geometry} {
      }

      explicit Separatrix(const Separatrix &separatrix)
        : source_{separatrix.source_},
          destination_{separatrix.destination_},
          geometry_{separatrix.geometry_},
          onBoundary_{separatrix.onBoundary_} {
      }

      explicit Separatrix(Separatrix &&separatrix) noexcept
        : source_{separatrix.source_},
          destination_{separatrix.destination_},
          geometry_{std::move(separatrix.geometry_)},
          onBoundary_{separatrix.onBoundary_} {
      }

      Separatrix &operator=(Separatrix &&separatrix) noexcept {
        source_ = separatrix.source_;
        destination_ = separatrix.destination_;
        geometry_ = std::move(separatrix.geometry_);
        onBoundary_ = separatrix.onBoundary_;

        return *this;
      }

      /**
       * Source triangle of the separatrix.
       */
      SimplexId source_{-1};

      /**
       * Destination triangle of the separatrix.
       */
      SimplexId destination_{-1};

      /**
       * Container of ids. Each id addresses a separate
       * container corresponding to a dense representation
       * of the geometry (i.e. separatricesGeometry).
       */
      std::vector<SimplexId> geometry_;

      /**
       * Flag indicating that the separatix (partially) uses the boundary
       */
      bool onBoundary_{false};

    };

    struct DoublyLinkedElement {
      explicit DoublyLinkedElement() {
      }

      // default :
      explicit DoublyLinkedElement(SimplexId newEdge) {
        neighbors_[0] = newEdge;
      }

      void insert(SimplexId newEdge) {
        neighbors_[1] = newEdge;
        hasTwoNeighbors = true;
      }

      SimplexId next(SimplexId edgeFrom) {
        if(!hasTwoNeighbors) {
          return -1;
        }

        if(neighbors_[0] == edgeFrom) {
          return neighbors_[1];
        }

        return neighbors_[0];
      }

      SimplexId neighbors_[2];
      bool hasTwoNeighbors{false};
    };
  }

  /**
   * The MorseSmaleComplexQuasi class provides methods to compute a
   * MSC variant that only uses the input field and its order field.
   */
  class MorseSmaleComplexQuasi : virtual public Debug {

  public:
    MorseSmaleComplexQuasi();

    int preconditionTriangulation(
      ttk::AbstractTriangulation *triangulation) const {
        int success = triangulation->preconditionVertexNeighbors();
        success += triangulation->preconditionVertexStars();
        success += triangulation->preconditionEdgeTriangles();
        success += triangulation->preconditionTriangleEdges();
        success += triangulation->preconditionVertexEdges();
        success += triangulation->preconditionCellEdges();
        success += triangulation->preconditionCellTriangles();
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

    inline int setOutputSeparatrices2(
      SimplexId *const separatrices2_numberOfPoints,
      std::vector<float> *const separatrices2_points,
      SimplexId *const separatrices2_numberOfCells,
      std::vector<SimplexId> *const separatrices2_cells_connectivity) {
      outputSeparatrices2_numberOfPoints_ = separatrices2_numberOfPoints;
      outputSeparatrices2_points_ = separatrices2_points;
      outputSeparatrices2_numberOfCells_ = separatrices2_numberOfCells;
      outputSeparatrices2_cells_connectivity_
        = separatrices2_cells_connectivity;
      return 0;
    }

    inline long long getSparseId(
      const SimplexId edgeId0, const SimplexId edgeId1) const {
      if(edgeId0 < edgeId1) {
        return edgeId0 * num_MSC_regions_ + edgeId1;
      }

      return edgeId1 * num_MSC_regions_ + edgeId0;
    }

    template <typename dataType, typename triangulationType>
    int execute(const triangulationType &triangulation);

    template <class dataType, class triangulationType>
    int computeDescendingManifold(
      SimplexId *const descendingManifold,
      std::vector<mscq::CriticalPoint> &criticalPoints,
      SimplexId &numberOfMaxima,
      const triangulationType &triangulation) const;

    template <class dataType, class triangulationType>
    int computeAscendingManifold(
      SimplexId *const ascendingManifold,
      std::vector<mscq::CriticalPoint> &criticalPoints,
      SimplexId &numberOfMinima,
      const triangulationType &triangulation) const;

    template <typename triangulationType>
    int computeFinalSegmentation(
      const SimplexId numMaxima,
      const SimplexId numMinima,
      const SimplexId *const ascendingManifold,
      const SimplexId *const descendingManifold,
      SimplexId &numManifolds,
      SimplexId *const morseSmaleManifold,
      const triangulationType &triangulation) const;

    template <class dataType, class triangulationType>
    int computeSeparatrices1_2D(
      std::vector<mscq::CriticalPoint> &criticalPoints,    
      std::vector<mscq::Separatrix> &separatrices,
      const SimplexId *const morseSmaleManifold,
      const triangulationType &triangulation) const;

    template <typename triangulationType>
    bool nextBorderEdge(
      const SimplexId currentVertex,
      const SimplexId currentEdge,
      SimplexId &nextVertex,
      SimplexId &nextEdge,
      const triangulationType &triangulation) const;

    template <class dataType, class triangulationType>
    bool isEdgeSaddle(
      const SimplexId edgeID,
      const triangulationType &triangulation) const;

    template <typename triangulationType>
    int setSeparatrices1(
      const std::vector<mscq::Separatrix> &separatrices,
      const triangulationType &triangulation) const;

    template <class dataType, class triangulationType>
    int computeSeparatrices1_3D(
      std::vector<std::array<float, 9>> &trianglePos,    
      std::vector<mscq::Separatrix> &separatrices,
      const SimplexId *const morseSmaleManifold,
      const triangulationType &triangulation) const;

    int setSeparatrices2(
      std::vector<std::array<float, 9>> trianglePos) const;

    template <class triangulationType>
    inline int getEdgeIncenter(
      SimplexId vertexId0, SimplexId vertexId1, float incenter[3],
      const triangulationType &triangulation) const {
      float p[6];
      triangulation.getVertexPoint(vertexId0, p[0], p[1], p[2]);
      triangulation.getVertexPoint(vertexId1, p[3], p[4], p[5]);

      incenter[0] = 0.5 * p[0] + 0.5 * p[3];
      incenter[1] = 0.5 * p[1] + 0.5 * p[4];
      incenter[2] = 0.5 * p[2] + 0.5 * p[5];

      return 0;
    }
      
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

    // 1-separatrices
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

    // 2-separatrices
    SimplexId *outputSeparatrices2_numberOfPoints_{};
    std::vector<float> *outputSeparatrices2_points_{};
    SimplexId *outputSeparatrices2_numberOfCells_{};
    std::vector<SimplexId> *outputSeparatrices2_cells_connectivity_{};

    // Misc
    SimplexId num_MSC_regions_{};

    // Debug
    bool db_omp = true;

  }; // MorseSmaleComplexQuasi class

} // namespace ttk

template <typename dataType, typename triangulationType>
int ttk::MorseSmaleComplexQuasi::execute(const triangulationType &triangulation)
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

  std::vector<mscq::CriticalPoint> criticalPoints;
  std::vector<mscq::Separatrix> separatrices1;
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
                                num_MSC_regions_,
                                morseSmaleManifold, triangulation);

    if(ascendingManifold or descendingManifold) {
      this->printMsg("Segmentation computed", 1.0, tmp.getElapsedTime(),
                     this->threadNumber_);
    }

    if(dim == 2) {
      computeSeparatrices1_2D<dataType, triangulationType>(
                                criticalPoints, separatrices1,
                                morseSmaleManifold, triangulation);

      this->printMsg("Write 1-seps");
      setSeparatrices1<triangulationType>(separatrices1, triangulation);
    } else if (dim == 3) {
    std::vector<std::array<float, 9>> trianglePos;

      computeSeparatrices1_3D<dataType, triangulationType>(
                              trianglePos, separatrices1,
                              morseSmaleManifold, triangulation);

      setSeparatrices2(trianglePos);
    }
  }

  this->printMsg( std::to_string(triangulation.getNumberOfVertices())
                   + " verticies processed",
                 1.0, t.getElapsedTime(), this->threadNumber_);

  return 0;
}

template <class dataType, class triangulationType>
int ttk::MorseSmaleComplexQuasi::computeDescendingManifold(
  SimplexId *const descendingManifold,
  std::vector<mscq::CriticalPoint> &criticalPoints,
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
      criticalPoints.push_back(mscq::CriticalPoint(it.first, dim));
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
int ttk::MorseSmaleComplexQuasi::computeAscendingManifold(
  SimplexId *const ascendingManifold,
  std::vector<mscq::CriticalPoint> &criticalPoints,
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
      criticalPoints.push_back(mscq::CriticalPoint(it.first, 0));
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
int ttk::MorseSmaleComplexQuasi::computeFinalSegmentation(
  const SimplexId numberOfMaxima,
  const SimplexId numberOfMinima,
  const SimplexId *const ascendingManifold,
  const SimplexId *const descendingManifold,
  SimplexId &numManifolds,
  SimplexId *const morseSmaleManifold,
  const triangulationType &triangulation) const{
  
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

  numManifolds = sparseRegionIds.size();

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

template <class dataType, class triangulationType>
int ttk::MorseSmaleComplexQuasi::computeSeparatrices1_2D(
  std::vector<mscq::CriticalPoint> &criticalPoints, 
  std::vector<mscq::Separatrix> &separatrixVector,
  const SimplexId *const morseSmaleManifold,
  const triangulationType &triangulation) const {

  // start a local timer for this subprocedure
  ttk::Timer localTimer;

  this->printMsg(ttk::debug::Separator::L1);
  this->printMsg("Computing 1-separatrices",
                 0, // progress form 0-1
                 0, // elapsed time so far
                 this->threadNumber_);//, ttk::debug::LineMode::REPLACE);
  // sparseId, vector of edgeID tuples
  std::unordered_map<
    long long, std::vector<std::tuple<SimplexId, SimplexId>>*> sepPieces;

  // sparseId, triangle0, edge0, triangle1, edge1
  std::unordered_map<long long, std::tuple<
    SimplexId, SimplexId, SimplexId, SimplexId>> triConnectors;
  std::unordered_set<long long> unhandledTriConnectors;

  // morseSmaleManifoldId, edgeId and corresponding vertexID
  std::unordered_map<
    SimplexId, std::tuple<SimplexId, SimplexId>> borderSepSeeds;

  const SimplexId numTriangles = triangulation.getNumberOfTriangles();

  for(SimplexId tri = 0; tri < numTriangles; tri++) {
    SimplexId edges[3];
    triangulation.getTriangleEdge(tri, 0, edges[0]);
    triangulation.getTriangleEdge(tri, 1, edges[1]);
    triangulation.getTriangleEdge(tri, 2, edges[2]);

    bool sameLabel[3];
    SimplexId msmEdge[3][2];

    for(SimplexId e = 0; e < 3; ++e) {
      SimplexId v0, v1;
      triangulation.getEdgeVertex(edges[e], 0, v0);
      triangulation.getEdgeVertex(edges[e], 1, v1);

      msmEdge[e][0] = morseSmaleManifold[v0];
      msmEdge[e][1] = morseSmaleManifold[v1];

      sameLabel[e] =
        morseSmaleManifold[v0] == morseSmaleManifold[v1];
    }

    int sameLabelCount =
      (int)sameLabel[0] + (int)sameLabel[1] + (int)sameLabel[2];

    if(sameLabelCount == 0) { // 3 labels on triangle
      for(SimplexId e = 0; e < 3; ++e) {
        const long long sparseID = getSparseId(msmEdge[e][0], msmEdge[e][1]);
        unhandledTriConnectors.insert(sparseID);

        if(triConnectors.find(sparseID) == triConnectors.end()) {
          triConnectors.insert(
            {sparseID, std::make_tuple(tri, edges[e], -1, -1)});
        } else {
          std::get<2>(triConnectors[sparseID]) = tri;
          std::get<3>(triConnectors[sparseID]) = edges[e];
        }
      }
    } else if (sameLabelCount == 1) { // 2 labels on triangle

      SimplexId e0, e1;
      if(sameLabel[0]) {
        e0 = 1; e1 = 2;
      } else if(sameLabel[1]) {
        e0 = 0; e1 = 2;
      } else {
        e0 = 0; e1 = 1;
      }

      const long long sparseID = getSparseId(msmEdge[e0][0], msmEdge[e0][1]);

      if(sepPieces.find(sparseID) == sepPieces.end()) { // new sparseID
        auto sepPiece = new std::vector<std::tuple<SimplexId, SimplexId>>;
        sepPiece->push_back(std::make_tuple(edges[e0], edges[e1]));
        sepPieces.insert({sparseID, sepPiece});
      } else { // existing sparseID
        sepPieces[sparseID]->push_back(
          std::make_tuple(edges[e0], edges[e1]));
      }
    }
  }

  for (auto& it_seps: sepPieces) { // put the sepPieces in sequence
    auto sepVectorPtr = it_seps.second;
    const SimplexId vecSize = sepVectorPtr->size();

    SimplexId maxIndex = 0;

    // write all Edges into a doubly Linked List
    mscq::DoublyLinkedElement linkedEdges[vecSize+1];
    std::unordered_map<SimplexId, SimplexId> edgeIdToArrIndex;
    {
      for (auto& it_vec: *sepVectorPtr) {
        const SimplexId firstEdge = std::get<0>(it_vec);
        const SimplexId secondEdge = std::get<1>(it_vec);

        if(edgeIdToArrIndex.find(firstEdge) == edgeIdToArrIndex.end()) {
          SimplexId newIndex = maxIndex++;
          edgeIdToArrIndex.insert({firstEdge, newIndex});
          linkedEdges[newIndex] = mscq::DoublyLinkedElement(secondEdge);
        } else {
          linkedEdges[edgeIdToArrIndex[firstEdge]].insert(secondEdge);
        }

        if(edgeIdToArrIndex.find(secondEdge) == edgeIdToArrIndex.end()) {
          SimplexId newIndex = maxIndex++;
          edgeIdToArrIndex.insert({secondEdge, newIndex});
          linkedEdges[newIndex] = mscq::DoublyLinkedElement(firstEdge);
        } else {
          linkedEdges[edgeIdToArrIndex[secondEdge]].insert(firstEdge);
        }
      }
    }

    // find a starting edge
    SimplexId startEdge{-1};
    {
      for (auto& it: edgeIdToArrIndex) {
        if(!linkedEdges[it.second].hasTwoNeighbors) {
          startEdge = it.first;
          break;
        }
      }
    }

    // write separatrix
    mscq::Separatrix newSep(startEdge);

    SimplexId v0, v1;
    triangulation.getEdgeVertex(startEdge, 0, v0);
    triangulation.getEdgeVertex(startEdge, 1, v1);

    const long long sparseID = getSparseId(
      morseSmaleManifold[v0], morseSmaleManifold[v1]);

    unhandledTriConnectors.erase(sparseID); // remove handled IDs

    //bool endsOnBorder = false;
    bool splitAtStart = false;
 
    if(triConnectors.find(sparseID) != triConnectors.end()) {
      SimplexId splitTriangle = std::get<0>(triConnectors[sparseID]);
      SimplexId splitEdge = std::get<1>(triConnectors[sparseID]);

      splitAtStart = splitEdge == startEdge;

      if(splitAtStart) {
        newSep.source_ = splitTriangle;
      } else {
        newSep.destination_ = splitTriangle;
      }

      splitTriangle = std::get<2>(triConnectors[sparseID]);
      splitEdge = std::get<3>(triConnectors[sparseID]);

      if(splitTriangle != -1) {
        if(splitAtStart) {
          newSep.destination_ = splitTriangle;
        } else {
         newSep.source_ = splitTriangle;
        }
      } else { //ends at border edge
        //endsOnBorder = true;
      }
    }
    
    SimplexId previousId = startEdge;
    SimplexId nextId = linkedEdges[edgeIdToArrIndex[startEdge]].neighbors_[0];
    newSep.geometry_.push_back(nextId);

    while(linkedEdges[edgeIdToArrIndex[nextId]].hasTwoNeighbors) {
      SimplexId tempPreviousId = nextId;
      nextId = linkedEdges[edgeIdToArrIndex[nextId]].next(previousId);
      previousId = tempPreviousId;
      newSep.geometry_.push_back(nextId);
    }

    separatrixVector.push_back(newSep);

    // add border seedpoints for border separatrices
    /*if(endsOnBorder) {
      SimplexId borderEdge = newSep.geometry_.front();
      if(splitAtStart) {
        borderEdge = newSep.geometry_.back();
      }

      /*triangulation.getEdgeVertex(borderEdge, 0, v0);
      triangulation.getEdgeVertex(borderEdge, 1, v1);

      if(isEdgeSaddle<dataType, triangulationType>(borderEdge, triangulation)) {
        std::array<float, 3> pt{};
        triangulation.getEdgeIncenter(borderEdge, pt.data());
        outputCriticalPoints_points_->push_back(pt);
        outputCriticalPoints_points_cellDimensions_->push_back(1);
        outputCriticalPoints_points_cellIds_->push_back(borderEdge);
        criticalPoints.push_back(mscq::CriticalPoint(borderEdge, 1, 1));
      }*/

      /*if(borderSepSeeds.find(morseSmaleManifold[v0]) == borderSepSeeds.end()) {
        borderSepSeeds.insert({
          morseSmaleManifold[v0], std::make_tuple(borderEdge, v0)});
      }
      if(borderSepSeeds.find(morseSmaleManifold[v1]) == borderSepSeeds.end()) {
        borderSepSeeds.insert({
          morseSmaleManifold[v1], std::make_tuple(borderEdge, v1)});
      }
    }*/
  }

  // merge two splitting triangles into a saddle
  for(auto u : unhandledTriConnectors) {
    const auto& conn = triConnectors[u];
    SimplexId splitTriangle0 = std::get<0>(conn);
    SimplexId splitEdge0 = std::get<1>(conn);
    SimplexId splitTriangle1 = std::get<2>(conn);

    mscq::Separatrix newSep(splitEdge0);
    newSep.source_ = splitTriangle0;
    newSep.destination_ = splitTriangle1;
    separatrixVector.push_back(newSep);

    std::array<float, 3> pt{};
    triangulation.getEdgeIncenter(splitEdge0, pt.data());
    outputCriticalPoints_points_->push_back(pt);
    outputCriticalPoints_points_cellDimensions_->push_back(1);
    outputCriticalPoints_points_cellIds_->push_back(splitEdge0);
    criticalPoints.push_back(mscq::CriticalPoint(splitEdge0, 1, 1));
  }

  /*for(auto b : borderSepSeeds) {
    SimplexId currentEdge = std::get<0>(b.second);
    SimplexId currentVertex = std::get<1>(b.second);

    mscq::Separatrix newSep(currentVertex);
    newSep.source_ = currentEdge;

    SimplexId nextEdge, nextVertex;
    nextBorderEdge(
      currentVertex, currentEdge, nextVertex, nextEdge, triangulation);
    
    while(
    morseSmaleManifold[currentVertex] == morseSmaleManifold[nextVertex]) {
    newSep.geometry_.push_back(currentVertex);
    currentVertex = nextVertex;
    currentEdge = nextEdge;
    nextBorderEdge(
      currentVertex, currentEdge, nextVertex, nextEdge, triangulation);
    }

    newSep.destination_ = nextEdge;
    newSep.onBoundary_ = true;
    separatrixVector.push_back(newSep);
  }*/

  return 0;
}

template <typename triangulationType>
bool ttk::MorseSmaleComplexQuasi::nextBorderEdge(
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

template <class dataType, class triangulationType>
bool ttk::MorseSmaleComplexQuasi::isEdgeSaddle(
  const SimplexId edgeID,
  const triangulationType &triangulation) const {
  
  const dataType *inputField = static_cast<const dataType *>(inputOrderField_);

  SimplexId v0, v1;
  triangulation.getEdgeVertex(edgeID, 0, v0);
  triangulation.getEdgeVertex(edgeID, 1, v1);

  bool hasBiggerNeighbor0 = false;
  bool hasSmallerNeigbor0 = false;

  SimplexId numNeigborsV0 = triangulation.getVertexNeighborNumber(v0);
  for(int i = 0; i < numNeigborsV0; ++i) {
    SimplexId neighbor;
    triangulation.getVertexNeighbor(v0, i, neighbor);

    if(neighbor != v1) {
      hasBiggerNeighbor0 |= (inputField[v0] < inputField[neighbor]);
      hasSmallerNeigbor0 |= (inputField[v0] > inputField[neighbor]);
    }
  }

  bool hasBiggerNeighbor1 = false;
  bool hasSmallerNeigbor1 = false;

  SimplexId numNeigborsV1 = triangulation.getVertexNeighborNumber(v1);
  for(int i = 0; i < numNeigborsV1; ++i) {
    SimplexId neighbor;
    triangulation.getVertexNeighbor(v1, i, neighbor);

    if(neighbor != v0) {
      hasBiggerNeighbor1 |= (inputField[v1] < inputField[neighbor]);
      hasSmallerNeigbor1 |= (inputField[v1] > inputField[neighbor]);
    }
  }
  
  return hasBiggerNeighbor0 && hasSmallerNeigbor0 &&
         hasBiggerNeighbor1 && hasSmallerNeigbor1;
}

template <typename triangulationType>
int ttk::MorseSmaleComplexQuasi::setSeparatrices1(
  const std::vector<mscq::Separatrix> &separatrices,
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

  // total number of separatrices points
  auto npoints{static_cast<size_t>(*outputSeparatrices1_numberOfPoints_)};
  // list of valid geometryId to flatten loops
  std::vector<SimplexId> validGeomIds{};
  // number of separatrices
  size_t numSep = separatrices.size();
  size_t numCells = 0;

  // count total number of points and cells, flatten geometryId loops
  for(size_t i = 0; i < numSep; ++i) {
    const auto &sep = separatrices[i];

    npoints += sep.geometry_.size() + 2;
    numCells += sep.geometry_.size() + 1;
  }

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
  for(size_t i = 0; i < numSep; ++i) {

    const auto &sep = separatrices[i];
    const auto &src = sep.source_;
    const auto &dst = sep.destination_;

    std::array<float, 3> pt{};

    if(src != -1) {
      triangulation.getTriangleIncenter(src, pt.data());   
    } else {
      triangulation.getEdgeIncenter(sep.geometry_.front(), pt.data());
    }

    points[3 * currPId + 0] = pt[0];
    points[3 * currPId + 1] = pt[1];
    points[3 * currPId + 2] = pt[2];

    currPId += 1;

    for(const auto geom : sep.geometry_) {
      triangulation.getEdgeIncenter(geom, pt.data());

      points[3 * currPId + 0] = pt[0];
      points[3 * currPId + 1] = pt[1];
      points[3 * currPId + 2] = pt[2];

      cellsConn[2 * currCId + 0] = currPId - 1;
      cellsConn[2 * currCId + 1] = currPId;

      currPId += 1;
      currCId += 1;
    }

    if(dst != -1) {
      triangulation.getTriangleIncenter(dst, pt.data());   
    } else {
      triangulation.getEdgeIncenter(sep.geometry_.back(), pt.data());
    }

    points[3 * currPId + 0] = pt[0];
    points[3 * currPId + 1] = pt[1];
    points[3 * currPId + 2] = pt[2];

    cellsConn[2 * currCId + 0] = currPId - 1;
    cellsConn[2 * currCId + 1] = currPId;

    currPId += 1;
    currCId += 1;

    if(outputSeparatrices1_cells_sourceIds_ != nullptr)
      (*outputSeparatrices1_cells_sourceIds_)[i] = src;
    if(outputSeparatrices1_cells_destinationIds_ != nullptr)
      (*outputSeparatrices1_cells_destinationIds_)[i] = dst;
    if(outputSeparatrices1_cells_isOnBoundary_ != nullptr)
      (*outputSeparatrices1_cells_isOnBoundary_)[i] = false;
  }

  // update pointers
  *outputSeparatrices1_numberOfPoints_ = npoints;
  *outputSeparatrices1_numberOfCells_ = numCells;

  return 0;
}

template <class dataType, class triangulationType>
int ttk::MorseSmaleComplexQuasi::computeSeparatrices1_3D(
  std::vector<std::array<float, 9>> &trianglePos,    
  std::vector<mscq::Separatrix> &separatrices,
  const SimplexId *const morseSmaleManifold,
  const triangulationType &triangulation) const {

  if(outputSeparatrices2_numberOfPoints_ == nullptr) {
    this->printErr("2-separatrices pointer to numberOfPoints is null.");
    return -1;
  }
  if(outputSeparatrices2_points_ == nullptr) {
    this->printErr("2-separatrices pointer to points is null.");
    return -1;
  }
  if(outputSeparatrices2_numberOfCells_ == nullptr) {
    this->printErr("2-separatrices pointer to numberOfCells is null.");
    return -1;
  }
  if(outputSeparatrices2_cells_connectivity_ == nullptr) {
    this->printErr("2-separatrices pointer to cells is null.");
    return -1;
  }
  
  int caseStatistics[28] = {0};
  const SimplexId numTetra = triangulation.getNumberOfCells();

  this->printMsg("Start tri calculation");

  for(SimplexId tet = 0; tet < numTetra; tet++) {
    SimplexId vertices[4];
    triangulation.getCellVertex(tet, 0, vertices[0]);
    triangulation.getCellVertex(tet, 1, vertices[1]);
    triangulation.getCellVertex(tet, 2, vertices[2]);
    triangulation.getCellVertex(tet, 3, vertices[3]);

    const SimplexId &msmV0 = morseSmaleManifold[vertices[0]];
    const SimplexId &msmV1 = morseSmaleManifold[vertices[1]];
    const SimplexId &msmV2 = morseSmaleManifold[vertices[2]];
    const SimplexId &msmV3 = morseSmaleManifold[vertices[3]];

    unsigned char index1 = (msmV0 == msmV1) ? 0x00 : 0x10; // 0 : 1
    unsigned char index2 = (msmV0 == msmV2) ? 0x00 :       // 0 :
                           (msmV1 == msmV2) ? 0x04 : 0x08; // 1 : 2
    unsigned char index3 = (msmV0 == msmV3) ? 0x00 :       // 0
                           (msmV1 == msmV3) ? 0x01 :       // 1
                           (msmV2 == msmV3) ? 0x02 : 0x03; // 2 : 3

    unsigned char lookupIndex = index1 | index2 | index3;
    int *tetEdgeIndices = tetraederLookup[lookupIndex];

    caseStatistics[lookupIndex] += 1;

    //this->printMsg("Index computed");

    //std::bitset<8> x(lookupIndex);
    //this->printMsg(std::to_string(lookupIndex));
    //std::cout << x;

    if(tetEdgeIndices[0] == -1) { // 1 label on tetraeder
      //this->printMsg("Case -1");
      continue;
    } else {
      //this->printMsg("Case 0");
      float edgeCenters[6][3];
      getEdgeIncenter(vertices[0], vertices[1], edgeCenters[0], triangulation);
      getEdgeIncenter(vertices[0], vertices[2], edgeCenters[1], triangulation);
      getEdgeIncenter(vertices[0], vertices[3], edgeCenters[2], triangulation);
      getEdgeIncenter(vertices[1], vertices[2], edgeCenters[3], triangulation);
      getEdgeIncenter(vertices[1], vertices[3], edgeCenters[4], triangulation);
      getEdgeIncenter(vertices[2], vertices[3], edgeCenters[5], triangulation);    

      if(tetEdgeIndices[0] == 6) { // 4 labels on tetraeder
      //this->printMsg("Case 1");
        float tetCenter[3];
        triangulation.getCellIncenter(tet, 3, tetCenter);

        // vertex 0
        trianglePos.push_back({
          edgeCenters[0][0], edgeCenters[0][1], edgeCenters[0][2], 
          edgeCenters[1][0], edgeCenters[1][1], edgeCenters[1][2], 
          tetCenter[0], tetCenter[1], tetCenter[2]
        });
        trianglePos.push_back({
          edgeCenters[0][0], edgeCenters[0][1], edgeCenters[0][2], 
          edgeCenters[2][0], edgeCenters[2][1], edgeCenters[2][2], 
          tetCenter[0], tetCenter[1], tetCenter[2]
        });
        trianglePos.push_back({
          edgeCenters[1][0], edgeCenters[1][1], edgeCenters[1][2], 
          edgeCenters[2][0], edgeCenters[2][1], edgeCenters[2][2], 
          tetCenter[0], tetCenter[1], tetCenter[2]
        });

        // vertex 1
        trianglePos.push_back({
          edgeCenters[1][0], edgeCenters[1][1], edgeCenters[1][2], 
          edgeCenters[3][0], edgeCenters[3][1], edgeCenters[3][2], 
          tetCenter[0], tetCenter[1], tetCenter[2]
        });
        trianglePos.push_back({
          edgeCenters[3][0], edgeCenters[3][1], edgeCenters[3][2], 
          edgeCenters[5][0], edgeCenters[5][1], edgeCenters[5][2], 
          tetCenter[0], tetCenter[1], tetCenter[2]
        });
        trianglePos.push_back({
          edgeCenters[1][0], edgeCenters[1][1], edgeCenters[1][2], 
          edgeCenters[5][0], edgeCenters[5][1], edgeCenters[5][2], 
          tetCenter[0], tetCenter[1], tetCenter[2]
        });

        // vertex 2
        trianglePos.push_back({
          edgeCenters[2][0], edgeCenters[2][1], edgeCenters[2][2], 
          edgeCenters[4][0], edgeCenters[4][1], edgeCenters[4][2], 
          tetCenter[0], tetCenter[1], tetCenter[2]
        });
        trianglePos.push_back({
          edgeCenters[4][0], edgeCenters[4][1], edgeCenters[4][2], 
          edgeCenters[5][0], edgeCenters[5][1], edgeCenters[5][2], 
          tetCenter[0], tetCenter[1], tetCenter[2]
        });
        trianglePos.push_back({
          edgeCenters[2][0], edgeCenters[2][1], edgeCenters[2][2], 
          edgeCenters[5][0], edgeCenters[5][1], edgeCenters[5][2], 
          tetCenter[0], tetCenter[1], tetCenter[2]
        });

        // vertex 3 - 0,2 and 2,4 already exist
        trianglePos.push_back({
          edgeCenters[0][0], edgeCenters[0][1], edgeCenters[0][2], 
          edgeCenters[4][0], edgeCenters[4][1], edgeCenters[4][2], 
          tetCenter[0], tetCenter[1], tetCenter[2]
        });
        trianglePos.push_back({
          edgeCenters[3][0], edgeCenters[3][1], edgeCenters[3][2], 
          edgeCenters[4][0], edgeCenters[4][1], edgeCenters[4][2], 
          tetCenter[0], tetCenter[1], tetCenter[2]
        });
        trianglePos.push_back({
          edgeCenters[0][0], edgeCenters[0][1], edgeCenters[0][2], 
          edgeCenters[3][0], edgeCenters[3][1], edgeCenters[3][2], 
          tetCenter[0], tetCenter[1], tetCenter[2]
        });
      } else { // 2 or 3 labels on tetraeder
        //this->printMsg("Case 2");
        trianglePos.push_back({
          edgeCenters[tetEdgeIndices[0]][0],
          edgeCenters[tetEdgeIndices[0]][1],
          edgeCenters[tetEdgeIndices[0]][2], 
          edgeCenters[tetEdgeIndices[1]][0],
          edgeCenters[tetEdgeIndices[1]][1],
          edgeCenters[tetEdgeIndices[1]][2], 
          edgeCenters[tetEdgeIndices[2]][0],
          edgeCenters[tetEdgeIndices[2]][1],
          edgeCenters[tetEdgeIndices[2]][2]
        });

        if(tetEdgeIndices[3] != -1) {
          trianglePos.push_back({
            edgeCenters[tetEdgeIndices[3]][0],
            edgeCenters[tetEdgeIndices[3]][1],
            edgeCenters[tetEdgeIndices[3]][2], 
            edgeCenters[tetEdgeIndices[4]][0],
            edgeCenters[tetEdgeIndices[4]][1],
            edgeCenters[tetEdgeIndices[4]][2], 
            edgeCenters[tetEdgeIndices[5]][0],
            edgeCenters[tetEdgeIndices[5]][1],
            edgeCenters[tetEdgeIndices[5]][2]
          });
          if(tetEdgeIndices[6] != -1) {
            trianglePos.push_back({
              edgeCenters[tetEdgeIndices[6]][0],
              edgeCenters[tetEdgeIndices[6]][1],
              edgeCenters[tetEdgeIndices[6]][2], 
              edgeCenters[tetEdgeIndices[7]][0],
              edgeCenters[tetEdgeIndices[7]][1],
              edgeCenters[tetEdgeIndices[7]][2], 
              edgeCenters[tetEdgeIndices[8]][0],
              edgeCenters[tetEdgeIndices[8]][1],
              edgeCenters[tetEdgeIndices[8]][2]
            });
            trianglePos.push_back({
              edgeCenters[tetEdgeIndices[9]][0],
              edgeCenters[tetEdgeIndices[9]][1],
              edgeCenters[tetEdgeIndices[9]][2], 
              edgeCenters[tetEdgeIndices[10]][0],
              edgeCenters[tetEdgeIndices[10]][1],
              edgeCenters[tetEdgeIndices[10]][2], 
              edgeCenters[tetEdgeIndices[11]][0],
              edgeCenters[tetEdgeIndices[11]][1],
              edgeCenters[tetEdgeIndices[11]][2]
            });
          }
        }
      }
    }
  }

  this->printMsg("Case AAAA: " + std::to_string(caseStatistics[0]));
  this->printMsg("Case AAAD: " + std::to_string(caseStatistics[3]));
  this->printMsg("Case AACA: " + std::to_string(caseStatistics[8]));
  this->printMsg("Case AACC: " + std::to_string(caseStatistics[10]));
  this->printMsg("Case AACD: " + std::to_string(caseStatistics[11]));
  this->printMsg("Case ABAA: " + std::to_string(caseStatistics[16]));
  this->printMsg("Case ABAB: " + std::to_string(caseStatistics[17]));
  this->printMsg("Case ABAD: " + std::to_string(caseStatistics[19]));
  this->printMsg("Case ABBA: " + std::to_string(caseStatistics[20]));
  this->printMsg("Case ABBB: " + std::to_string(caseStatistics[21]));
  this->printMsg("Case ABBD: " + std::to_string(caseStatistics[23]));
  this->printMsg("Case ABCA: " + std::to_string(caseStatistics[24]));
  this->printMsg("Case ABCB: " + std::to_string(caseStatistics[25]));
  this->printMsg("Case ABCC: " + std::to_string(caseStatistics[26]));
  this->printMsg("Case ABCD: " + std::to_string(caseStatistics[27]));
  


  return 0;
}

int ttk::MorseSmaleComplexQuasi::setSeparatrices2(
  std::vector<std::array<float, 9>> trianglePos) const {

  const int numTriangles = trianglePos.size();

  size_t npoints = numTriangles * 3;
  size_t numCells = numTriangles;

  // resize arrays
  outputSeparatrices2_points_->resize(3 * npoints);
  auto &points = *outputSeparatrices2_points_;
  outputSeparatrices2_cells_connectivity_->resize(3 * numCells);
  auto &cellsConn = *outputSeparatrices2_cells_connectivity_;

  this->printMsg("Start writing tris" + std::to_string(numTriangles));

  int lastPtId = 0;
  int lastCellId = 0;

  for(int tri = 0; tri < numTriangles; ++tri) {
    points[9 * tri + 0] = trianglePos[tri][0];
    points[9 * tri + 1] = trianglePos[tri][1];
    points[9 * tri + 2] = trianglePos[tri][2];

    points[9 * tri + 3] = trianglePos[tri][3];
    points[9 * tri + 4] = trianglePos[tri][4];
    points[9 * tri + 5] = trianglePos[tri][5];

    points[9 * tri + 6] = trianglePos[tri][6];
    points[9 * tri + 7] = trianglePos[tri][7];
    points[9 * tri + 8] = trianglePos[tri][8];

    cellsConn[3 * tri]     = 3 * tri;
    cellsConn[3 * tri + 1] = 3 * tri + 1;
    cellsConn[3 * tri + 2] = 3 * tri + 2;

    lastPtId = 9 * tri + 8;
    lastCellId = 3 * tri + 2;
  }

  this->printMsg(std::to_string(lastPtId) + " - " + std::to_string(3 * npoints));
  this->printMsg(std::to_string(lastCellId) + " - " + std::to_string(3 * numCells));
  this->printMsg("Done writing" + std::to_string(numTriangles));

  // update pointers
  *outputSeparatrices2_numberOfPoints_ = npoints;
  *outputSeparatrices2_numberOfCells_ = numCells;
  this->printMsg("End writing 2-seps");

  return 1;
}