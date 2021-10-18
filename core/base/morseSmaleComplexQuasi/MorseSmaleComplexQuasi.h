/// \ingroup base
/// \class ttk::MorseSmaleComplexQuasi
/// \author Robin G. C. Maack <maack@rhrk.uni-kl.de>
/// \date June 2021.

#pragma once

#include <Debug.h>
#include <Triangulation.h>
#include <UnionFind.h>

#include <unordered_map>
#include <map>
#include <list>
#include <stack>
#include <set>
#include <unordered_set>
#include <bitset>
#include <vector>
#include <algorithm>

#include <cstddef>
#include <bits/stdc++.h> 

using ttk::SimplexId;

// For cacheline width calculation
#ifdef __cpp_lib_hardware_interference_size
    using std::hardware_destructive_interference_size;
#else
    constexpr std::size_t hardware_destructive_interference_size
        = 2 * sizeof(std::max_align_t);
#endif

// (#Labels) 2nd label, 3rd label, 4th label - index
// first label is always 0 
int tetraederLookup[28][15] = {
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (1) 0,0,0 - 0
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (-) 0,0,1 - 1
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (-) 0,0,2 - 2 
  { 2,  4,  5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (2) 0,0,3 - 3
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (-) 0,1,0 - 4
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (-) 0,1,1 - 5 
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (-) 0,1,2 - 6
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (-) 0,1,3 - 7
  { 1,  3,  5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (2) 0,2,0 - 8
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (-) 0,2,1 - 9
  { 1,  2,  4,  1,  3,  4, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (2) 0,2,2 -10
  { 5,  8,  9,  1,  9,  8,  1,  9,  3,  2,  9,  4,  2,  9,  8}, // (3) 0,2,3 -11
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (-) 0,3,0 -12
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (-) 0,3,1 -13
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (-) 0,3,2 -14
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (-) 0,3,3 -15
  { 0,  3,  4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (2) 1,0,0 -16
  { 0,  2,  3,  2,  3,  5, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (2) 1,0,1 -17
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (-) 1,0,2 -18
  { 4,  7,  9,  0,  7,  9,  0,  3,  9,  2,  7,  9,  2,  5,  9}, // (3) 1,0,3 -19
  { 0,  1,  4,  1,  4,  5, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (2) 1,1,0 -20
  { 0,  1,  2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (2) 1,1,1 -21
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (-) 1,1,2 -22
  { 2,  7,  8,  1,  7,  8,  0,  1,  7,  4,  7,  8,  4,  5,  8}, // (3) 1,1,3 -23
  { 3,  6,  9,  4,  6,  9,  0,  4,  6,  1,  6,  9,  1,  5,  9}, // (3) 1,2,0 -24
  { 1,  6,  8,  2,  6,  8,  2,  6,  0,  3,  8,  5,  3,  8,  6}, // (3) 1,2,1 -25
  { 0,  6,  7,  2,  6,  1,  2,  6,  7,  3,  7,  4,  3,  7,  6}, // (3) 1,2,2 -26
  {10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10}  // (4) 1,2,3 -27
};

bool tetraederLookupFast[28] = {
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

      SimplexId id_{-1};
      int index_{-1};
    };

    struct Separatrix {
      // default :
      explicit Separatrix()
        : geometry_{}, endsOnBoundary_{false}, onBoundary_{false} {
      }

      // initialization with one segment :
      explicit Separatrix(const SimplexId segmentGeometry) {
        geometry_.push_back(segmentGeometry);
      }

      explicit Separatrix(
        const SimplexId segmentGeometry,
        const bool endsOnBoundary,
        const bool onBoundary) {
        geometry_.push_back(segmentGeometry);
        endsOnBoundary_  = endsOnBoundary;
        onBoundary_ = onBoundary;
      }

      // initialization with multiple segments :
      //explicit Separatrix(const std::vector<SimplexId> &geometry)
      //  : geometry_{geometry} {
      //}

      explicit Separatrix(const Separatrix &separatrix)
        : geometry_{separatrix.geometry_},
          endsOnBoundary_{separatrix.endsOnBoundary_},
          onBoundary_{separatrix.onBoundary_} {
      }

      explicit Separatrix(Separatrix &&separatrix) noexcept
          : geometry_{std::move(separatrix.geometry_)},
          endsOnBoundary_{separatrix.endsOnBoundary_},
          onBoundary_{separatrix.onBoundary_} {
      }

      Separatrix &operator=(Separatrix &&separatrix) noexcept {
        geometry_ = std::move(separatrix.geometry_);
        endsOnBoundary_ = separatrix.endsOnBoundary_;
        onBoundary_ = separatrix.onBoundary_;

        return *this;
      }

      Separatrix &operator=(const Separatrix &separatrix) noexcept {
        geometry_ = std::move(separatrix.geometry_);
        endsOnBoundary_ = separatrix.endsOnBoundary_;
        onBoundary_ = separatrix.onBoundary_;

        return *this;
      }

      /**
       * Container of ids. Each id addresses a separate
       * container corresponding to a dense representation
       * of the geometry (i.e. separatricesGeometry).
       */
      std::vector<SimplexId> geometry_;

      /**
       * Flag indicating that the separatix lies on the boundary
       */
      bool endsOnBoundary_{false};

      /**
       * Flag indicating that the separatix lies on the boundary
       */
      bool onBoundary_{false};
    };

    struct DoublyLinkedElement {
      explicit DoublyLinkedElement() {
      }

      explicit DoublyLinkedElement(SimplexId neighbor) {
        neighbors_[0] = neighbor;
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
      bool visited{false};
    };
  }

  /**
   * The MorseSmaleComplexQuasi class provides methods to compute a
   * MSC variant that only uses the input field and its order field.
   */
  class MorseSmaleComplexQuasi : virtual public Debug {

  public:
    MorseSmaleComplexQuasi();

    enum class SEPARATRICES_MANIFOLD {
      MORSESMALE = 0,
      ASCENDING = 1,
      DESCENDING = 2
    };

    int preconditionTriangulation(
      ttk::AbstractTriangulation *triangulation) const {
        int success = triangulation->preconditionVertexNeighbors();
        success += triangulation->preconditionTriangleEdges();
        success += triangulation->preconditionCellTriangles();
        success += triangulation->preconditionVertexStars();
        success += triangulation->preconditionEdgeTriangles();
        success += triangulation->preconditionTriangleStars();
        success += triangulation->preconditionBoundaryVertices();
        success += triangulation->preconditionBoundaryTriangles();
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
      std::vector<SimplexId> *const separatrices2_cells_connectivity,
      std::vector<SimplexId> *const separatrices2_cells_cases) {
      outputSeparatrices2_numberOfPoints_ = separatrices2_numberOfPoints;
      outputSeparatrices2_points_ = separatrices2_points;
      outputSeparatrices2_numberOfCells_ = separatrices2_numberOfCells;
      outputSeparatrices2_cells_connectivity_
        = separatrices2_cells_connectivity;
      outputSeparatrices2_cells_cases = separatrices2_cells_cases;
      return 0;
    }

    inline long long getSparseId(
      const long long simID0, const long long simID1) const {
      if(simID0 < simID1) {
        return simID0 * num_MSC_regions_ + simID1;
      }
      return simID1 * num_MSC_regions_ + simID0;
    }

    inline long long getSparseId(
      const long long simID0, const long long simID1,
      const long long simID2) const {
      if(simID0 < simID1) {
        if(simID1 < simID2) {
          return simID0 * num_MSC_regions_ * num_MSC_regions_ +
          simID1 * num_MSC_regions_ + simID2;
        } else {
          if(simID0 < simID2) {
            return simID0 * num_MSC_regions_ * num_MSC_regions_ +
            simID2 * num_MSC_regions_ + simID1;
          } else {
            return simID2 * num_MSC_regions_ * num_MSC_regions_ +
            simID0 * num_MSC_regions_ + simID1;
          }
        }
      } else {
        if(simID0 < simID2) {
          return simID1 * num_MSC_regions_ * num_MSC_regions_ +
          simID0 * num_MSC_regions_ + simID2;
        } else {
          if(simID1 < simID2) {
            return simID1 * num_MSC_regions_ * num_MSC_regions_ +
            simID2 * num_MSC_regions_ + simID0;
          } else {
            return simID2 * num_MSC_regions_ * num_MSC_regions_ +
            simID1 * num_MSC_regions_ + simID0;
          }
        }
      }
    }

    template <typename dataType, typename triangulationType>
    int execute(
      const triangulationType &triangulation, 
      SEPARATRICES_MANIFOLD sepManifoldType,
      const bool computeSeparatrices,
      const bool fastSeparatrices);

    template <typename dataType, typename triangulationType>
    int computeManifold(
      SimplexId *const manifold,
      std::vector<mscq::CriticalPoint> &criticalPoints,
      SimplexId &numberOfMaxima,
      const triangulationType &triangulation,
      const bool ascending) const;

    template <typename triangulationType>
    int computeFinalSegmentation(
      const SimplexId numMaxima,
      const std::vector<mscq::CriticalPoint> &criticalPoints,
      const SimplexId *const ascendingManifold,
      const SimplexId *const descendingManifold,
      SimplexId &numManifolds,
      SimplexId *const morseSmaleManifold,
      const triangulationType &triangulation) const;

    template <typename dataType, typename triangulationType>
    int computeSeparatrices1Pieces_2D(
      std::unordered_set<SimplexId> &saddleCandidates,
      std::unordered_map<long long, std::vector<
        std::tuple<SimplexId, SimplexId>>*> &sepPieces,
      std::unordered_map<long long, std::vector<
        std::tuple<SimplexId, SimplexId> >> &triConnectors,
      const SimplexId *const morseSmaleManifold,
      const triangulationType &triangulation) const;

    template <typename dataType, typename triangulationType>
    int computeSeparatrices1Inner_2D(
      std::unordered_set<SimplexId> &saddleCandidates,   
      std::vector<mscq::Separatrix> &separatrices,
      std::unordered_map<long long, std::vector<
        std::tuple<SimplexId, SimplexId>>*> sepPieces,
      std::unordered_map<long long, std::vector<
        std::tuple<SimplexId, SimplexId> >> triConnectors,
      const SimplexId *const morseSmaleManifold,
      const triangulationType &triangulation) const;

    template <typename dataType, typename triangulationType>
    int computeSeparatrices1Border_2D(
      std::unordered_set<SimplexId> &saddleCandidates,   
      std::vector<mscq::Separatrix> &separatrices,
      const SimplexId *const morseSmaleManifold,
      const triangulationType &triangulation) const;

    template <typename triangulationType>
    int setSeparatrices1_2D(
      const std::vector<mscq::Separatrix> &separatrices,
      const triangulationType &triangulation) const;

    template <typename dataType, typename triangulationType>
    int compute2Separatrices_3D(
      std::vector<std::array<float, 9>> &trianglePos,    
      std::vector<SimplexId> &caseData,
      std::unordered_set<SimplexId> *saddleCandidates,
      std::unordered_map<long long, std::vector<
        std::tuple<SimplexId, SimplexId, bool>>> &pieces1separatrices,
      const SimplexId *const morseSmaleManifold,
      const triangulationType &triangulation) const;

    template <typename dataType, typename triangulationType>
    int compute1SeparatricesInner_3D(
      std::unordered_map<long long, std::vector<
      std::tuple<SimplexId, SimplexId, bool>>> &pieces1separatrices,
      std::vector<mscq::Separatrix> &separatrices1,
      std::unordered_set<SimplexId> &borderSeeds,
      const SimplexId *const morseSmaleManifold,
      const triangulationType &triangulation) const;

    template <typename triangulationType>
    int compute1SeparatricesBorder_3D(
      std::unordered_set<SimplexId> &borderSeeds,
      std::vector<mscq::Separatrix> &separatrices1,
      const SimplexId *const morseSmaleManifold,
      const triangulationType &triangulation) const;

    template <typename dataType, typename triangulationType>
    int computeSeparatrices_3D_fast(
      std::vector<std::array<float, 9>> &trianglePos,    
      std::vector<SimplexId> &caseData,
      const SimplexId *const morseSmaleManifold,
      const triangulationType &triangulation) const;

    template <typename triangulationType>
    int createBorderPath_3D(
      mscq::Separatrix &separatrix,
      SimplexId startTriangle,
      SimplexId startEdge,
      const SimplexId *const morseSmaleManifold,
      const triangulationType &triangulation
    ) const;

    int extractMSCPaths(
      std::vector<mscq::CriticalPoint> &criticalPoints,    
      std::vector<mscq::Separatrix> &separatrices) const;

    template <typename triangulationType>
    int setSeparatrices1_3D(
      const std::vector<mscq::Separatrix> &separatrices,
      const triangulationType &triangulation) const;

    int setSeparatrices2_3D(
      std::vector<std::array<float, 9>> &trianglePos,
      std::vector<SimplexId> &caseData) const;

    template <typename dataType, typename triangulationType>
    int findSaddlesFromCadidates(
      std::vector<mscq::CriticalPoint> &criticalPoints, 
      const std::unordered_set<SimplexId> &saddleCandidates,
      const triangulationType &triangulation) const;

    template <typename dataType, typename triangulationType>
    int setCriticalPoints(
      std::vector<mscq::CriticalPoint> &criticalPoints,
      const triangulationType &triangulation) const;

    template <typename dataType, typename triangulationType>
    std::pair<ttk::SimplexId, ttk::SimplexId> getNumberOfLowerUpperComponents(
      const SimplexId vertexId,
      const dataType * const offsets,
      const triangulationType &triangulation) const;

    template <typename triangulationType>
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

    inline int getCenter(
      float pos0[3], float pos1[3], float pos2[3], float incenter[3]) const {
      float d[3];
      d[0] = Geometry::distance(pos1, pos2);
      d[1] = Geometry::distance(pos0, pos2);
      d[2] = Geometry::distance(pos0, pos1);
      const float sum = d[0] + d[1] + d[2];

      d[0] = d[0] / sum;
      d[1] = d[1] / sum;
      d[2] = d[2] / sum;

      incenter[0] = d[0] * pos0[0] + d[1] * pos1[0] + d[2] * pos2[0];
      incenter[1] = d[0] * pos0[1] + d[1] * pos1[1] + d[2] * pos2[1];
      incenter[2] = d[0] * pos0[2] + d[1] * pos1[2] + d[2] * pos2[2];

      return 0;
    }
      
  protected:
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
    std::vector<SimplexId> *outputSeparatrices2_cells_cases{};

    // Misc
    SimplexId num_MSC_regions_{};

    // Debug
    bool db_omp = true;

  }; // MorseSmaleComplexQuasi class

} // namespace ttk

template <typename dataType, typename triangulationType>
int ttk::MorseSmaleComplexQuasi::execute(
  const triangulationType &triangulation, 
  SEPARATRICES_MANIFOLD sepManifoldType,
  const bool computeSeparatrices,
  const bool fastSeparatrices)
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
  {
    Timer tmp;

    const int dim = triangulation.getDimensionality();

    SimplexId numberOfMaxima{0};
    SimplexId numberOfMinima{0};

    if(ascendingManifold)
      computeManifold<dataType, triangulationType>(
                                ascendingManifold, criticalPoints,
                                numberOfMinima, triangulation, true);

    if(descendingManifold)                            
      computeManifold<dataType, triangulationType>(
                                descendingManifold, criticalPoints,
                                numberOfMaxima, triangulation, false);

    if(ascendingManifold && descendingManifold && morseSmaleManifold) {
      computeFinalSegmentation<triangulationType>(
                                numberOfMinima, criticalPoints,
                                ascendingManifold, descendingManifold,
                                num_MSC_regions_, 
                                morseSmaleManifold, triangulation);
    }

    if(computeSeparatrices) {
      SimplexId *sepManifold;
      bool SepManifoldIsValid = false;

      switch(sepManifoldType) {
        case SEPARATRICES_MANIFOLD::MORSESMALE :
          if(ascendingManifold && descendingManifold && morseSmaleManifold) {
            sepManifold = static_cast<SimplexId *>(outputMorseSmaleManifold_);
            SepManifoldIsValid = true;
          }
          break;
        case SEPARATRICES_MANIFOLD::ASCENDING :
          if(ascendingManifold) {
            sepManifold = static_cast<SimplexId *>(outputAscendingManifold_);
            num_MSC_regions_ = numberOfMinima;
            SepManifoldIsValid = true;
          }
          break;
        case SEPARATRICES_MANIFOLD::DESCENDING :
          if(descendingManifold) {
            sepManifold = static_cast<SimplexId *>(outputDescendingManifold_);
            num_MSC_regions_ = numberOfMaxima;
            SepManifoldIsValid = true;
          }
          break;
        default:
          break;
      }

      if(SepManifoldIsValid) {
        std::unordered_set<SimplexId> saddleCandidates;
        std::vector<mscq::Separatrix> separatrices1;

        if(dim == 2) {
          // sparseId, vector of edgeID tuples
          std::unordered_map<long long, std::vector<
            std::tuple<SimplexId, SimplexId>>*> sepPieces;

          // 3 label triangles
          std::unordered_map<long long, std::vector<
            std::tuple<SimplexId, SimplexId> >> triConnectors;
          computeSeparatrices1Pieces_2D<dataType, triangulationType>(
            saddleCandidates, sepPieces, triConnectors,
            sepManifold, triangulation
          );

          this->printMsg("YXX_" + std::to_string(sepPieces.size()));

          computeSeparatrices1Inner_2D<dataType, triangulationType>(
            saddleCandidates, separatrices1, sepPieces, triConnectors,
            sepManifold, triangulation
          );

          this->printMsg("XXX_" + std::to_string(separatrices1.size()));

          computeSeparatrices1Border_2D<dataType, triangulationType>(
            saddleCandidates, separatrices1, sepManifold, triangulation
          );

          this->printMsg("XXX_" + std::to_string(separatrices1.size()));

          findSaddlesFromCadidates<dataType, triangulationType>(
            criticalPoints, saddleCandidates, triangulation
          );

          this->printMsg("XXX_" + std::to_string(separatrices1.size()));

          setCriticalPoints<dataType, triangulationType>(
            criticalPoints, triangulation
          );

          setSeparatrices1_2D<triangulationType>(
            separatrices1, triangulation
          );
        } else if(dim == 3) {
          std::vector<std::array<float, 9>> trianglePos; // 2-separatrices
          std::vector<SimplexId> caseData; // trianglePos array
          std::unordered_map<long long, std::vector<
            std::tuple<SimplexId, SimplexId, bool>>> pieces1separatrices;
          std::unordered_set<SimplexId> borderSeeds;

          if(fastSeparatrices) {
            computeSeparatrices_3D_fast<dataType, triangulationType>(
                                      trianglePos, caseData,
                                      sepManifold, triangulation);

            setSeparatrices2_3D(trianglePos, caseData);
          } else {
            

            compute2Separatrices_3D<dataType, triangulationType>(
              trianglePos, caseData, &saddleCandidates, pieces1separatrices,
              sepManifold, triangulation
            );

            findSaddlesFromCadidates<dataType, triangulationType>(
              criticalPoints, saddleCandidates, triangulation
            );

            compute1SeparatricesInner_3D<dataType, triangulationType>(
              pieces1separatrices, separatrices1, borderSeeds,
              sepManifold, triangulation
            );

            compute1SeparatricesBorder_3D<triangulationType>(
              borderSeeds, separatrices1, sepManifold, triangulation
            );

            setCriticalPoints<dataType, triangulationType>(
              criticalPoints, triangulation
            );

            setSeparatrices1_3D<triangulationType>(
              separatrices1, triangulation
            );
            setSeparatrices2_3D(trianglePos, caseData);
          }
        }
      }
    }
  }

  this->printMsg(ttk::debug::Separator::L1);
  this->printMsg( std::to_string(triangulation.getNumberOfVertices())
                   + " verticies processed",
                 1.0, t.getElapsedTime(), this->threadNumber_);
  this->printMsg(ttk::debug::Separator::L0);

  return 0;
}


template <typename dataType, typename triangulationType>
int ttk::MorseSmaleComplexQuasi::computeManifold(
  SimplexId *const descendingManifold,
  std::vector<mscq::CriticalPoint> &criticalPoints,
  SimplexId &numberOfMaxima,
  const triangulationType &triangulation,
  const bool ascending) const {

  const dataType *inputField = static_cast<const dataType *>(inputOrderField_);

  // start global timer
  ttk::Timer localTimer;
  int iterations = 1;
  numberOfMaxima = 0;

  const std::string manifoldStr = ascending ? "Ascending" : "Descending";
  const std::string minMaxStr = ascending ? "Minima" : "Maxima";

  // -----------------------------------------------------------------------
  // Compute MSC variant
  // -----------------------------------------------------------------------
  {
    this->printMsg(ttk::debug::Separator::L1);
    // print the progress of the current subprocedure (currently 0%)
    this->printMsg("Computing " + manifoldStr + " Manifold",
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
    for(SimplexId i = 0; i < nVertices; i++) {
      const SimplexId numNeighbors =
        triangulation.getVertexNeighborNumber(i);

      SimplexId neighborId;
      bool hasBiggerNeighbor = false;

      SimplexId &dmi = descendingManifold[i];
      dmi = i;

      // check all neighbors
      for(SimplexId n = 0; n < numNeighbors; n++) {
        triangulation.getVertexNeighbor(i, n, neighborId);
        if(ascending) {
          if(inputField[neighborId] < inputField[dmi]) {
            dmi = neighborId;
            hasBiggerNeighbor = true;
          }
        } else {
          if(inputField[neighborId] > inputField[dmi]) {
            dmi = neighborId;
            hasBiggerNeighbor = true;
          }
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
      for(SimplexId i = 0; i < nVertices; i++) {
        SimplexId neighborId;
        SimplexId numNeighbors =
          triangulation.getVertexNeighborNumber(i);
        bool hasBiggerNeighbor = false;
        SimplexId &dmi = descendingManifold[i];

        dmi = i;

        // check all neighbors
        for(SimplexId n = 0; n < numNeighbors; n++) {
          triangulation.getVertexNeighbor(i, n, neighborId);
          if(ascending) {
            if(inputField[neighborId] < inputField[dmi]) {
              dmi = neighborId;
              hasBiggerNeighbor = true;
            }
          } else {
            if(inputField[neighborId] > inputField[dmi]) {
              dmi = neighborId;
              hasBiggerNeighbor = true;
            }
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
} //#endif // TTK_ENABLE_OPENMP

    const int index =
      ascending ? 0 : triangulation.getDimensionality();
    //const int dim = triangulation.getDimensionality();

    for (auto& it: maxima) {
      criticalPoints.push_back(mscq::CriticalPoint(it.first, index));
    }

    delete activeVertices;
  }

  {
    this->printMsg("Computed " + manifoldStr + " Manifold",
                   1, localTimer.getElapsedTime(), this->threadNumber_);

    this->printMsg("#" + minMaxStr +  ": " + std::to_string(numberOfMaxima),
                   1, localTimer.getElapsedTime(), this->threadNumber_);

    this->printMsg("#Iterations: " + std::to_string(iterations - 1),
                   1, localTimer.getElapsedTime(), this->threadNumber_);

    this->printMsg(ttk::debug::Separator::L1); // horizontal '=' separator
  }

  return 1; // return success
}


template <typename triangulationType>
int ttk::MorseSmaleComplexQuasi::computeFinalSegmentation(
  const SimplexId numberOfMinima,
  const std::vector<mscq::CriticalPoint> &criticalPoints,
  const SimplexId *const ascendingManifold,
  const SimplexId *const descendingManifold,
  SimplexId &numManifolds,
  SimplexId *const morseSmaleManifold,
  const triangulationType &triangulation) const{
  
  const size_t nVerts = triangulation.getNumberOfVertices();

//#ifndef TTK_ENABLE_OPENMP
if(!db_omp) {
  for(size_t i = 0; i < nVerts; ++i) {
    morseSmaleManifold[i] = ascendingManifold[i] * numberOfMinima +
      descendingManifold[i];
  }
//#else // TTK_ENABLE_OPENMP
} else {
  const size_t numCacheLineEntries =
  hardware_destructive_interference_size / sizeof(SimplexId);
  #pragma omp parallel num_threads(threadNumber_)
  {
    #pragma omp for schedule(static, numCacheLineEntries) nowait
    for(size_t i = 0; i < nVerts; ++i) {
      morseSmaleManifold[i] = ascendingManifold[i] * numberOfMinima +
        descendingManifold[i];
    }
  }
}

  // associate a unique "sparse region id" to each (ascending, descending) pair

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

template <typename dataType, typename triangulationType>
int ttk::MorseSmaleComplexQuasi::computeSeparatrices1Pieces_2D(
  std::unordered_set<SimplexId> &saddleCandidates,
  std::unordered_map<long long, std::vector<
    std::tuple<SimplexId, SimplexId>>*> &sepPieces,
  std::unordered_map<long long, std::vector<
    std::tuple<SimplexId, SimplexId> >> &triConnectors,
  const SimplexId *const morseSmaleManifold,
  const triangulationType &triangulation) const {
  // start a local timer for this subprocedure
  ttk::Timer localTimer;

  this->printMsg("Computing 1-separatrices pieces",
                 0, // progress form 0-1
                 0, // elapsed time so far
                 this->threadNumber_, ttk::debug::LineMode::REPLACE);

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

        if(triConnectors.find(sparseID) == triConnectors.end()) {
          triConnectors.insert(
            {sparseID, std::vector<std::tuple<SimplexId, SimplexId>>()} );
        }
        triConnectors[sparseID].push_back(std::make_tuple(edges[e], tri));

        SimplexId saddleCandidateVertex[3];
        triangulation.getTriangleVertex(tri, 0, saddleCandidateVertex[0]);
        triangulation.getTriangleVertex(tri, 1, saddleCandidateVertex[1]);
        triangulation.getTriangleVertex(tri, 2, saddleCandidateVertex[2]);

        saddleCandidates.insert(saddleCandidateVertex[0]);
        saddleCandidates.insert(saddleCandidateVertex[1]);
        saddleCandidates.insert(saddleCandidateVertex[2]);
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

  this->printMsg("Computed 1-separatrices pieces",
                 0, // progress form 0-1
                 localTimer.getElapsedTime(), // elapsed time so far
                 this->threadNumber_);

  return 1;
}

template <typename dataType, typename triangulationType>
int ttk::MorseSmaleComplexQuasi::computeSeparatrices1Inner_2D(
  std::unordered_set<SimplexId> &saddleCandidates,   
  std::vector<mscq::Separatrix> &separatrices,
  std::unordered_map<long long, std::vector<
    std::tuple<SimplexId, SimplexId>>*> sepPieces,
  std::unordered_map<long long, std::vector<
    std::tuple<SimplexId, SimplexId> >> triConnectors,
  const SimplexId *const morseSmaleManifold,
  const triangulationType &triangulation) const {
  // start a local timer for this subprocedure
  ttk::Timer localTimer;

  this->printMsg("Computing inner 1-separatrices",
                 0, // progress form 0-1
                 0, // elapsed time so far
                 this->threadNumber_, ttk::debug::LineMode::REPLACE);

  // morseSmaleManifoldId, edgeId and corresponding vertexID
  std::unordered_map<
    SimplexId, std::tuple<SimplexId, SimplexId>> borderSepSeeds;

  for (auto& it_seps: sepPieces) { // put the sepPieces in sequence
    
    auto sepVectorPtr = it_seps.second;
    const SimplexId vecSize = sepVectorPtr->size();

    SimplexId maxIndex = 0;

    // write all Edges into a doubly Linked List
    std::vector<mscq::DoublyLinkedElement> linkedEdges;
    linkedEdges.reserve(vecSize + 1);
    std::unordered_map<SimplexId, SimplexId> edgeIdToArrIndex;
    {
      for (auto& it_vec: *sepVectorPtr) {
        const SimplexId firstEdge = std::get<0>(it_vec);
        const SimplexId secondEdge = std::get<1>(it_vec);

        if(edgeIdToArrIndex.find(firstEdge) == edgeIdToArrIndex.end()) {
          SimplexId newIndex = maxIndex++;
          edgeIdToArrIndex.insert({firstEdge, newIndex});
          linkedEdges.push_back(mscq::DoublyLinkedElement(secondEdge));
        } else {
          linkedEdges[edgeIdToArrIndex[firstEdge]].insert(secondEdge);
        }

        if(edgeIdToArrIndex.find(secondEdge) == edgeIdToArrIndex.end()) {
          SimplexId newIndex = maxIndex++;
          edgeIdToArrIndex.insert({secondEdge, newIndex});
          linkedEdges.push_back(mscq::DoublyLinkedElement(firstEdge));
        } else {
          linkedEdges[edgeIdToArrIndex[secondEdge]].insert(firstEdge);
        }
      }
    }

    // find a starting edge
    std::unordered_map<SimplexId, bool> startIds;
    for(auto& it: edgeIdToArrIndex) {
      if(!linkedEdges[it.second].hasTwoNeighbors) {
        startIds.insert({it.first, false});
      }
    }

    for (const auto& it_starts: startIds) { // detect paths
      if(it_starts.second) { // Is simplex already part of a separatrix?
        continue;
      }
      // write separatrix
      const SimplexId startEdge = it_starts.first;
      mscq::Separatrix newSep;
      newSep.onBoundary_ = false;

      SimplexId startEdgeVerts[2];
      triangulation.getEdgeVertex(startEdge, 0, startEdgeVerts[0]);
      triangulation.getEdgeVertex(startEdge, 1, startEdgeVerts[1]);

      const long long sparseID = getSparseId(
        morseSmaleManifold[startEdgeVerts[0]],
        morseSmaleManifold[startEdgeVerts[1]]
      );

      auto &triConnVect = triConnectors[sparseID];
      const size_t triConnVectSize = triConnVect.size();

      bool startsOnBoundary = true;
      bool endsOnBoundary = true;

      // Find the starting Triangle, if it exists
      for(size_t ec = 0; ec < triConnVectSize; ++ec) {
        if(std::get<0>(triConnVect[ec]) == startEdge) {
          newSep.geometry_.push_back(std::get<1>(triConnVect[ec]));
          triConnectors[sparseID].erase( triConnectors[sparseID].begin() + ec);
          startsOnBoundary = false;
          break;
        }
      }
      
      SimplexId previousId = startEdge;
      newSep.geometry_.push_back(startEdge);
      SimplexId nextId = linkedEdges[edgeIdToArrIndex[startEdge]].neighbors_[0];
      linkedEdges[edgeIdToArrIndex[previousId]].visited = true;
      newSep.geometry_.push_back(nextId);

      while(linkedEdges[edgeIdToArrIndex[nextId]].hasTwoNeighbors &&
        !linkedEdges[edgeIdToArrIndex[nextId]].visited) {
        SimplexId tempPreviousId = nextId;
        nextId = linkedEdges[edgeIdToArrIndex[nextId]].next(previousId);
        previousId = tempPreviousId;
        newSep.geometry_.push_back(nextId);
        linkedEdges[edgeIdToArrIndex[previousId]].visited = true;
      }

      // Find the ending Triangle, if it exists
      for(size_t ec = 0; ec < triConnVectSize; ++ec) {
        if(std::get<0>(triConnVect[ec]) == nextId) {
          newSep.geometry_.push_back(std::get<1>(triConnVect[ec]));
          triConnectors[sparseID].erase( triConnectors[sparseID].begin() + ec);
          endsOnBoundary = false;
          break;
        }
      }
      
      startIds[nextId] = true;

      if(endsOnBoundary) {
        SimplexId &dEdge = newSep.geometry_.back();

        SimplexId vert0, vert1;
        triangulation.getEdgeVertex(dEdge, 0, vert0);
        triangulation.getEdgeVertex(dEdge, 1, vert1);

        saddleCandidates.insert(vert0);
        saddleCandidates.insert(vert1);

        newSep.endsOnBoundary_ = true;
      }
      if(startsOnBoundary) {
        SimplexId &sEdge = newSep.geometry_.front();

        SimplexId vert0, vert1;
        triangulation.getEdgeVertex(sEdge, 0, vert0);
        triangulation.getEdgeVertex(sEdge, 1, vert1);

        saddleCandidates.insert(vert0);
        saddleCandidates.insert(vert1);

        newSep.endsOnBoundary_ = true;

        std::reverse(newSep.geometry_.begin(), newSep.geometry_.end());
      }

      separatrices.push_back(newSep);
    }
  }

  // merge two splitting triangles into a saddle
  for(auto u : triConnectors) {
    if(u.second.size() > 0) {
      for(size_t i = 0; i < u.second.size(); ++i) {
        for (size_t j = i+1; j < u.second.size(); ++j) {
          if(std::get<0>(u.second[i]) == std::get<0>(u.second[j])) {
            SimplexId &splitEdge0 = std::get<0>(u.second[i]);
            SimplexId &splitTriangle0 = std::get<1>(u.second[i]);
            SimplexId &splitTriangle1 = std::get<1>(u.second[j]);

            mscq::Separatrix newSep;
            newSep.geometry_.push_back(splitTriangle0);
            newSep.geometry_.push_back(splitEdge0);
            newSep.geometry_.push_back(splitTriangle1);
            separatrices.push_back(newSep);

            u.second.erase( u.second.begin() + j);
          }
        }
      }
    }
  }
  
  this->printMsg("Computed inner 1-separatrices",
                 0, // progress form 0-1
                 localTimer.getElapsedTime(), // elapsed time so far
                 this->threadNumber_);

  return 1;
}

template <typename dataType, typename triangulationType>
int ttk::MorseSmaleComplexQuasi::computeSeparatrices1Border_2D(
  std::unordered_set<SimplexId> &saddleCandidates,   
  std::vector<mscq::Separatrix> &separatrices,
  const SimplexId *const morseSmaleManifold,
  const triangulationType &triangulation) const {
  // start a local timer for this subprocedure
  ttk::Timer localTimer;

  this->printMsg("Computing border 1-separatrices",
                 0, // progress form 0-1
                 0, // elapsed time so far
                 this->threadNumber_, ttk::debug::LineMode::REPLACE);



  this->printMsg("Computed border 1-separatrices",
                 0, // progress form 0-1
                 localTimer.getElapsedTime(), // elapsed time so far
                 this->threadNumber_);
  
  return 1;
}

template <typename triangulationType>
int ttk::MorseSmaleComplexQuasi::setSeparatrices1_2D(
  const std::vector<mscq::Separatrix> &separatrices,
  const triangulationType &triangulation) const {

  ttk::Timer localTimer;

  this->printMsg("Writing 1-separatrices:" + std::to_string(separatrices.size()),
                 1, 0, this->threadNumber_,
                 ttk::debug::LineMode::REPLACE);

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

  // number of separatrices
  size_t numSep = separatrices.size();
  size_t npoints = 0;
  size_t numCells = 0;

  // count total number of points and cells, flatten geometryId loops
  for(size_t i = 0; i < numSep; ++i) {
    const auto &sep = separatrices[i];

    npoints += sep.geometry_.size(); // +0?
    numCells += sep.geometry_.size() - 1; // -1?
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

    std::array<float, 3> pt{};

    if(!sep.onBoundary_) { // inner separatrix
      triangulation.getTriangleIncenter(
        (sep.geometry_[0]), pt.data());

      points[3 * currPId + 0] = pt[0];
      points[3 * currPId + 1] = pt[1];
      points[3 * currPId + 2] = pt[2];

      currPId += 1;

      size_t geomSize = sep.geometry_.size() - 1;

      for(size_t j = 1; j < geomSize; ++j) {
        triangulation.getEdgeIncenter(sep.geometry_[j], pt.data());

        points[3 * currPId + 0] = pt[0];
        points[3 * currPId + 1] = pt[1];
        points[3 * currPId + 2] = pt[2];

        cellsConn[2 * currCId + 0] = currPId - 1;
        cellsConn[2 * currCId + 1] = currPId;

        currPId += 1;
        currCId += 1;
      }

      if(sep.endsOnBoundary_){
        triangulation.getEdgeIncenter(sep.geometry_[geomSize], pt.data());
      } else {
        triangulation.getTriangleIncenter(
          sep.geometry_[geomSize], pt.data()
        );
      }

      points[3 * currPId + 0] = pt[0];
      points[3 * currPId + 1] = pt[1];
      points[3 * currPId + 2] = pt[2];

      cellsConn[2 * currCId + 0] = currPId - 1;
      cellsConn[2 * currCId + 1] = currPId;

      currPId += 1;
      currCId += 1;

      if(outputSeparatrices1_cells_isOnBoundary_ != nullptr)
        (*outputSeparatrices1_cells_isOnBoundary_)[i] = false;

    } else { // boundary separatrix

      /*triangulation.getTriangleIncenter(sep.geometry_.front(), pt.data());

      points[3 * currPId + 0] = pt[0];
      points[3 * currPId + 1] = pt[1];
      points[3 * currPId + 2] = pt[2];

      currPId += 1;

      size_t geomSize = sep.geometry_.size() - 1;

      for(size_t j = 1; j < geomSize; ++j) {
        triangulation.getEdgeIncenter(sep.geometry_[j], pt.data());

        points[3 * currPId + 0] = pt[0];
        points[3 * currPId + 1] = pt[1];
        points[3 * currPId + 2] = pt[2];

        cellsConn[2 * currCId + 0] = currPId - 1;
        cellsConn[2 * currCId + 1] = currPId;

        currPId += 1;
        currCId += 1;
      }

      triangulation.getTriangleIncenter(sep.geometry_.back(), pt.data());

      points[3 * currPId + 0] = pt[0];
      points[3 * currPId + 1] = pt[1];
      points[3 * currPId + 2] = pt[2];

      cellsConn[2 * currCId + 0] = currPId - 1;
      cellsConn[2 * currCId + 1] = currPId;

      currPId += 1;
      currCId += 1;

      if(outputSeparatrices1_cells_isOnBoundary_ != nullptr)
        (*outputSeparatrices1_cells_isOnBoundary_)[i] = true;*/
    }
  }

  // update pointers
  *outputSeparatrices1_numberOfPoints_ = npoints;
  *outputSeparatrices1_numberOfCells_ = numCells;

  this->printMsg("Wrote 1-separatrices",
                 1, localTimer.getElapsedTime(), this->threadNumber_);

  return 0;
}

template <typename dataType, typename triangulationType>
int ttk::MorseSmaleComplexQuasi::compute2Separatrices_3D(
  std::vector<std::array<float, 9>> &trianglePos,    
  std::vector<SimplexId> &caseData,
  std::unordered_set<SimplexId> *saddleCandidates,
  std::unordered_map<long long, std::vector<
    std::tuple<SimplexId, SimplexId, bool>>> &pieces1separatrices,
  const SimplexId *const morseSmaleManifold,
  const triangulationType &triangulation) const {
  
  ttk::Timer localTimer;

  // print the progress of the current subprocedure (currently 0%)
  this->printMsg("Computing 2-Separatrices 3D",
                 0, // progress form 0-1
                 0, // elapsed time so far
                 this->threadNumber_, ttk::debug::LineMode::REPLACE);

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

  const SimplexId numTetra = triangulation.getNumberOfCells();
  const bool testFindSaddles = saddleCandidates != nullptr;

//#ifndef TTK_ENABLE_OPENMP
if(!db_omp) {
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
    unsigned char index2 = (msmV0 == msmV2) ? 0x00 :       // 0
                           (msmV1 == msmV2) ? 0x04 : 0x08; // 1 : 2
    unsigned char index3 = (msmV0 == msmV3) ? 0x00 :       // 0
                           (msmV1 == msmV3) ? 0x01 :       // 1
                           (msmV2 == msmV3) ? 0x02 : 0x03; // 2 : 3

    unsigned char lookupIndex = index1 | index2 | index3;
    int *tetEdgeIndices = tetraederLookup[lookupIndex];

    if(tetEdgeIndices[0] == -1) { // 1 label on tetraeder
      continue;
    } else {
      float edgeCenters[10][3];
      getEdgeIncenter(vertices[0], vertices[1], edgeCenters[0], triangulation);
      getEdgeIncenter(vertices[0], vertices[2], edgeCenters[1], triangulation);
      getEdgeIncenter(vertices[0], vertices[3], edgeCenters[2], triangulation);
      getEdgeIncenter(vertices[1], vertices[2], edgeCenters[3], triangulation);
      getEdgeIncenter(vertices[1], vertices[3], edgeCenters[4], triangulation);
      getEdgeIncenter(vertices[2], vertices[3], edgeCenters[5], triangulation);

      float vertPos[4][3];
      triangulation.getVertexPoint(
        vertices[0], vertPos[0][0], vertPos[0][1], vertPos[0][2]);
      triangulation.getVertexPoint(
        vertices[1], vertPos[1][0], vertPos[1][1], vertPos[1][2]);
      triangulation.getVertexPoint(
        vertices[2], vertPos[2][0], vertPos[2][1], vertPos[2][2]);
      triangulation.getVertexPoint(
        vertices[3], vertPos[3][0], vertPos[3][1], vertPos[3][2]);

      getCenter(vertPos[0], vertPos[1], vertPos[2], edgeCenters[6]);
      getCenter(vertPos[0], vertPos[1], vertPos[3], edgeCenters[7]);
      getCenter(vertPos[0], vertPos[2], vertPos[3], edgeCenters[8]);
      getCenter(vertPos[1], vertPos[2], vertPos[3], edgeCenters[9]);

      if(tetEdgeIndices[0] == 10) { // 4 labels on tetraeder
        if(testFindSaddles) {
          saddleCandidates->insert(vertices[0]);
          saddleCandidates->insert(vertices[1]);
          saddleCandidates->insert(vertices[2]);
          saddleCandidates->insert(vertices[3]);
        }

        float tetCenter[3];
        triangulation.getCellIncenter(tet, 3, tetCenter);

        // vertex 0
        trianglePos.push_back({
          edgeCenters[5][0], edgeCenters[5][1], edgeCenters[5][2], 
          edgeCenters[8][0], edgeCenters[8][1], edgeCenters[8][2], 
          tetCenter[0], tetCenter[1], tetCenter[2]
        });
        trianglePos.push_back({
          edgeCenters[8][0], edgeCenters[8][1], edgeCenters[8][2], 
          edgeCenters[1][0], edgeCenters[1][1], edgeCenters[1][2], 
          tetCenter[0], tetCenter[1], tetCenter[2]
        });
        trianglePos.push_back({
          edgeCenters[1][0], edgeCenters[1][1], edgeCenters[1][2], 
          edgeCenters[6][0], edgeCenters[6][1], edgeCenters[6][2], 
          tetCenter[0], tetCenter[1], tetCenter[2]
        });
        trianglePos.push_back({
          edgeCenters[6][0], edgeCenters[6][1], edgeCenters[6][2], 
          edgeCenters[3][0], edgeCenters[3][1], edgeCenters[3][2], 
          tetCenter[0], tetCenter[1], tetCenter[2]
        });
        trianglePos.push_back({
          edgeCenters[3][0], edgeCenters[3][1], edgeCenters[3][2], 
          edgeCenters[9][0], edgeCenters[9][1], edgeCenters[9][2], 
          tetCenter[0], tetCenter[1], tetCenter[2]
        });
        trianglePos.push_back({
          edgeCenters[9][0], edgeCenters[9][1], edgeCenters[9][2], 
          edgeCenters[5][0], edgeCenters[5][1], edgeCenters[5][2], 
          tetCenter[0], tetCenter[1], tetCenter[2]
        });

        // vertex 1
        trianglePos.push_back({
          edgeCenters[8][0], edgeCenters[8][1], edgeCenters[8][2], 
          edgeCenters[2][0], edgeCenters[2][1], edgeCenters[2][2], 
          tetCenter[0], tetCenter[1], tetCenter[2]
        });
        trianglePos.push_back({
          edgeCenters[2][0], edgeCenters[2][1], edgeCenters[2][2], 
          edgeCenters[7][0], edgeCenters[7][1], edgeCenters[7][2], 
          tetCenter[0], tetCenter[1], tetCenter[2]
        });
        trianglePos.push_back({
          edgeCenters[7][0], edgeCenters[7][1], edgeCenters[7][2], 
          edgeCenters[0][0], edgeCenters[0][1], edgeCenters[0][2], 
          tetCenter[0], tetCenter[1], tetCenter[2]
        });
        trianglePos.push_back({
          edgeCenters[0][0], edgeCenters[0][1], edgeCenters[0][2], 
          edgeCenters[6][0], edgeCenters[6][1], edgeCenters[6][2], 
          tetCenter[0], tetCenter[1], tetCenter[2]
        });

        // vertex 2
        trianglePos.push_back({
          edgeCenters[7][0], edgeCenters[7][1], edgeCenters[7][2], 
          edgeCenters[4][0], edgeCenters[4][1], edgeCenters[4][2], 
          tetCenter[0], tetCenter[1], tetCenter[2]
        });
        trianglePos.push_back({
          edgeCenters[4][0], edgeCenters[4][1], edgeCenters[4][2], 
          edgeCenters[9][0], edgeCenters[9][1], edgeCenters[9][2], 
          tetCenter[0], tetCenter[1], tetCenter[2]
        });

        caseData.insert(caseData.end(),{
          lookupIndex, lookupIndex, lookupIndex,
          lookupIndex, lookupIndex, lookupIndex,
          lookupIndex, lookupIndex, lookupIndex,
          lookupIndex, lookupIndex, lookupIndex}
        );

        for(int tri = 0; tri < 4; ++tri) {
          SimplexId triID;
          triangulation.getCellTriangle(tet, tri, triID);

          SimplexId triVertIds[3];
          triangulation.getTriangleVertex(triID, 0, triVertIds[0]);
          triangulation.getTriangleVertex(triID, 1, triVertIds[1]);
          triangulation.getTriangleVertex(triID, 2, triVertIds[2]);

          const long long sparseID = getSparseId(
            morseSmaleManifold[triVertIds[0]],
            morseSmaleManifold[triVertIds[1]],
            morseSmaleManifold[triVertIds[2]]
          );

          if(pieces1separatrices.find(sparseID) == pieces1separatrices.end()) {
            pieces1separatrices.insert({sparseID,
              std::vector<std::tuple<SimplexId, SimplexId, bool>>()}
            );
          }
          pieces1separatrices[sparseID].push_back(
            std::make_tuple(triID, tet, true));
        }
      } else { // 2 or 3 labels on tetraeder
        if(testFindSaddles) {
          if(tetEdgeIndices[6] == -1) {

            if(triangulation.isVertexOnBoundary(vertices[0]))
              saddleCandidates->insert(vertices[0]);

            if(triangulation.isVertexOnBoundary(vertices[1]))
              saddleCandidates->insert(vertices[1]);

            if(triangulation.isVertexOnBoundary(vertices[2]))
              saddleCandidates->insert(vertices[2]);

            if(triangulation.isVertexOnBoundary(vertices[3]))
              saddleCandidates->insert(vertices[3]);
          }
        }

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

        if(tetEdgeIndices[3] != -1) { // 2 or 3 labels on tetraeder
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

          if(tetEdgeIndices[6] != -1) { // 3 labels on tetraeder
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
            trianglePos.push_back({
              edgeCenters[tetEdgeIndices[12]][0],
              edgeCenters[tetEdgeIndices[12]][1],
              edgeCenters[tetEdgeIndices[12]][2], 
              edgeCenters[tetEdgeIndices[13]][0],
              edgeCenters[tetEdgeIndices[13]][1],
              edgeCenters[tetEdgeIndices[13]][2], 
              edgeCenters[tetEdgeIndices[14]][0],
              edgeCenters[tetEdgeIndices[14]][1],
              edgeCenters[tetEdgeIndices[14]][2]
            });

            caseData.insert(caseData.end(), {
              lookupIndex, lookupIndex, lookupIndex, lookupIndex, lookupIndex
            });

            SimplexId sepTriangles[2];
            bool foundFirst = false;
            long long sparseID;

            // find triangles that have 3 unique labels
            for(int tri = 0; tri < 4; ++tri) {
              SimplexId TriIDBuffer;
              triangulation.getCellTriangle(tet, tri, TriIDBuffer);

              SimplexId vert[3];
              triangulation.getTriangleVertex(TriIDBuffer, 0, vert[0]);
              triangulation.getTriangleVertex(TriIDBuffer, 1, vert[1]);
              triangulation.getTriangleVertex(TriIDBuffer, 2, vert[2]);

              if( // 3 unique labels on triangle
                morseSmaleManifold[vert[0]] != morseSmaleManifold[vert[1]] &&
                morseSmaleManifold[vert[0]] != morseSmaleManifold[vert[2]] &&
                morseSmaleManifold[vert[1]] != morseSmaleManifold[vert[2]]) {
                if(foundFirst) {
                  sepTriangles[1] = TriIDBuffer;
                  break;
                } else {
                  sepTriangles[0] = TriIDBuffer;
                  sparseID = getSparseId(morseSmaleManifold[vert[0]],
                    morseSmaleManifold[vert[1]], morseSmaleManifold[vert[2]]
                  );
                  foundFirst = true;
                }
              }
            }

            if(pieces1separatrices.find(sparseID) ==
              pieces1separatrices.end()) {
              pieces1separatrices.insert({
                sparseID,
                std::vector<std::tuple<SimplexId, SimplexId, bool>>()});
            }
            pieces1separatrices[sparseID].push_back(
              std::make_tuple(sepTriangles[0], sepTriangles[1], false)
            );
          } else {
            caseData.insert(caseData.end(), {lookupIndex, lookupIndex});
          }
        } else {
          caseData.push_back(lookupIndex);
        }
      }
    }
  }
//#else // TTK_ENABLE_OPENMP
} else {
  const size_t numCacheLineEntries =
      hardware_destructive_interference_size / sizeof(SimplexId);

  #pragma omp parallel num_threads(this->threadNumber_)
  {
    std::vector<std::array<float, 9>> trianglePosLocal;
    std::vector<SimplexId> caseDataLocal;
    std::vector<std::tuple<long long, SimplexId, SimplexId, bool>> sepPcsLocal;
    std::vector<SimplexId> saddleCandidatesLocal;

    #pragma omp for schedule(dynamic, numCacheLineEntries) nowait
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
      unsigned char index2 = (msmV0 == msmV2) ? 0x00 :       // 0
                             (msmV1 == msmV2) ? 0x04 : 0x08; // 1 : 2
      unsigned char index3 = (msmV0 == msmV3) ? 0x00 :       // 0
                             (msmV1 == msmV3) ? 0x01 :       // 1
                             (msmV2 == msmV3) ? 0x02 : 0x03; // 2 : 3

      unsigned char lookupIndex = index1 | index2 | index3;
      int *tetEdgeIndices = tetraederLookup[lookupIndex];

      if(tetEdgeIndices[0] == -1) { // 1 label on tetraeder
        continue;
      } else {
        float edgeCenters[10][3];
        getEdgeIncenter(
          vertices[0], vertices[1], edgeCenters[0], triangulation);
        getEdgeIncenter(
          vertices[0], vertices[2], edgeCenters[1], triangulation);
        getEdgeIncenter(
          vertices[0], vertices[3], edgeCenters[2], triangulation);
        getEdgeIncenter(
          vertices[1], vertices[2], edgeCenters[3], triangulation);
        getEdgeIncenter(
          vertices[1], vertices[3], edgeCenters[4], triangulation);
        getEdgeIncenter(
          vertices[2], vertices[3], edgeCenters[5], triangulation);

        float vertPos[4][3];
        triangulation.getVertexPoint(
          vertices[0], vertPos[0][0], vertPos[0][1], vertPos[0][2]);
        triangulation.getVertexPoint(
          vertices[1], vertPos[1][0], vertPos[1][1], vertPos[1][2]);
        triangulation.getVertexPoint(
          vertices[2], vertPos[2][0], vertPos[2][1], vertPos[2][2]);
        triangulation.getVertexPoint(
          vertices[3], vertPos[3][0], vertPos[3][1], vertPos[3][2]);

        getCenter(vertPos[0], vertPos[1], vertPos[2], edgeCenters[6]);
        getCenter(vertPos[0], vertPos[1], vertPos[3], edgeCenters[7]);
        getCenter(vertPos[0], vertPos[2], vertPos[3], edgeCenters[8]);
        getCenter(vertPos[1], vertPos[2], vertPos[3], edgeCenters[9]);

        if(tetEdgeIndices[0] == 10) { // 4 labels on tetraeder
          if(testFindSaddles) {
            saddleCandidatesLocal.push_back(vertices[0]);
            saddleCandidatesLocal.push_back(vertices[1]);
            saddleCandidatesLocal.push_back(vertices[2]);
            saddleCandidatesLocal.push_back(vertices[3]);
          }

          float tetCenter[3];
          triangulation.getCellIncenter(tet, 3, tetCenter);

          // vertex 0
          trianglePosLocal.push_back({
            edgeCenters[5][0], edgeCenters[5][1], edgeCenters[5][2], 
            edgeCenters[8][0], edgeCenters[8][1], edgeCenters[8][2], 
            tetCenter[0], tetCenter[1], tetCenter[2]
          });
          trianglePosLocal.push_back({
            edgeCenters[8][0], edgeCenters[8][1], edgeCenters[8][2], 
            edgeCenters[1][0], edgeCenters[1][1], edgeCenters[1][2], 
            tetCenter[0], tetCenter[1], tetCenter[2]
          });
          trianglePosLocal.push_back({
            edgeCenters[1][0], edgeCenters[1][1], edgeCenters[1][2], 
            edgeCenters[6][0], edgeCenters[6][1], edgeCenters[6][2], 
            tetCenter[0], tetCenter[1], tetCenter[2]
          });
          trianglePosLocal.push_back({
            edgeCenters[6][0], edgeCenters[6][1], edgeCenters[6][2], 
            edgeCenters[3][0], edgeCenters[3][1], edgeCenters[3][2], 
            tetCenter[0], tetCenter[1], tetCenter[2]
          });
          trianglePosLocal.push_back({
            edgeCenters[3][0], edgeCenters[3][1], edgeCenters[3][2], 
            edgeCenters[9][0], edgeCenters[9][1], edgeCenters[9][2], 
            tetCenter[0], tetCenter[1], tetCenter[2]
          });
          trianglePosLocal.push_back({
            edgeCenters[9][0], edgeCenters[9][1], edgeCenters[9][2], 
            edgeCenters[5][0], edgeCenters[5][1], edgeCenters[5][2], 
            tetCenter[0], tetCenter[1], tetCenter[2]
          });

          // vertex 1
          trianglePosLocal.push_back({
            edgeCenters[8][0], edgeCenters[8][1], edgeCenters[8][2], 
            edgeCenters[2][0], edgeCenters[2][1], edgeCenters[2][2], 
            tetCenter[0], tetCenter[1], tetCenter[2]
          });
          trianglePosLocal.push_back({
            edgeCenters[2][0], edgeCenters[2][1], edgeCenters[2][2], 
            edgeCenters[7][0], edgeCenters[7][1], edgeCenters[7][2], 
            tetCenter[0], tetCenter[1], tetCenter[2]
          });
          trianglePosLocal.push_back({
            edgeCenters[7][0], edgeCenters[7][1], edgeCenters[7][2], 
            edgeCenters[0][0], edgeCenters[0][1], edgeCenters[0][2], 
            tetCenter[0], tetCenter[1], tetCenter[2]
          });
          trianglePosLocal.push_back({
            edgeCenters[0][0], edgeCenters[0][1], edgeCenters[0][2], 
            edgeCenters[6][0], edgeCenters[6][1], edgeCenters[6][2], 
            tetCenter[0], tetCenter[1], tetCenter[2]
          });

          // vertex 2
          trianglePosLocal.push_back({
            edgeCenters[7][0], edgeCenters[7][1], edgeCenters[7][2], 
            edgeCenters[4][0], edgeCenters[4][1], edgeCenters[4][2], 
            tetCenter[0], tetCenter[1], tetCenter[2]
          });
          trianglePosLocal.push_back({
            edgeCenters[4][0], edgeCenters[4][1], edgeCenters[4][2], 
            edgeCenters[9][0], edgeCenters[9][1], edgeCenters[9][2], 
            tetCenter[0], tetCenter[1], tetCenter[2]
          });

          caseDataLocal.insert(caseDataLocal.end(),{
            lookupIndex, lookupIndex, lookupIndex,
            lookupIndex, lookupIndex, lookupIndex,
            lookupIndex, lookupIndex, lookupIndex,
            lookupIndex, lookupIndex, lookupIndex}
          );

          for(int tri = 0; tri < 4; ++tri) {
            SimplexId triID;
            triangulation.getCellTriangle(tet, tri, triID);

            SimplexId triVertIds[3];
            triangulation.getTriangleVertex(triID, 0, triVertIds[0]);
            triangulation.getTriangleVertex(triID, 1, triVertIds[1]);
            triangulation.getTriangleVertex(triID, 2, triVertIds[2]);

            const long long sparseID = getSparseId(
              morseSmaleManifold[triVertIds[0]],
              morseSmaleManifold[triVertIds[1]],
              morseSmaleManifold[triVertIds[2]]
            );

            sepPcsLocal.push_back(std::make_tuple(sparseID, triID, tet, true));
          }
        } else { // 2 or 3 labels on tetraeder
          if(testFindSaddles) {
            if(tetEdgeIndices[6] == -1) {
              if(triangulation.isVertexOnBoundary(vertices[0]))
                saddleCandidatesLocal.push_back(vertices[0]);

              if(triangulation.isVertexOnBoundary(vertices[1]))
                saddleCandidatesLocal.push_back(vertices[1]);

              if(triangulation.isVertexOnBoundary(vertices[2]))
                saddleCandidatesLocal.push_back(vertices[2]);

              if(triangulation.isVertexOnBoundary(vertices[3]))
                saddleCandidatesLocal.push_back(vertices[3]);
            }
          }

          trianglePosLocal.push_back({
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

          if(tetEdgeIndices[3] != -1) { // 2 or 3 labels on tetraeder
            trianglePosLocal.push_back({
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

            if(tetEdgeIndices[6] != -1) { // 3 labels on tetraeder
              trianglePosLocal.push_back({
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
              trianglePosLocal.push_back({
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
              trianglePosLocal.push_back({
                edgeCenters[tetEdgeIndices[12]][0],
                edgeCenters[tetEdgeIndices[12]][1],
                edgeCenters[tetEdgeIndices[12]][2], 
                edgeCenters[tetEdgeIndices[13]][0],
                edgeCenters[tetEdgeIndices[13]][1],
                edgeCenters[tetEdgeIndices[13]][2], 
                edgeCenters[tetEdgeIndices[14]][0],
                edgeCenters[tetEdgeIndices[14]][1],
                edgeCenters[tetEdgeIndices[14]][2]
              });

              caseDataLocal.insert(caseDataLocal.end(), {
                lookupIndex, lookupIndex, lookupIndex, lookupIndex, lookupIndex
              });

              SimplexId sepTriangles[2];
              bool foundFirst = false;
              long long sparseID;

              // find triangles that have 3 unique labels
              for(int tri = 0; tri < 4; ++tri) {
                SimplexId TriIDBuffer;
                triangulation.getCellTriangle(tet, tri, TriIDBuffer);

                SimplexId vert[3];
                triangulation.getTriangleVertex(TriIDBuffer, 0, vert[0]);
                triangulation.getTriangleVertex(TriIDBuffer, 1, vert[1]);
                triangulation.getTriangleVertex(TriIDBuffer, 2, vert[2]);

                if( // 3 unique labels on triangle
                morseSmaleManifold[vert[0]] != morseSmaleManifold[vert[1]] &&
                morseSmaleManifold[vert[0]] != morseSmaleManifold[vert[2]] &&
                morseSmaleManifold[vert[1]] != morseSmaleManifold[vert[2]]) {
                  if(foundFirst) {
                    sepTriangles[1] = TriIDBuffer;
                    break;
                  } else {
                    sepTriangles[0] = TriIDBuffer;
                    sparseID = getSparseId(morseSmaleManifold[vert[0]],
                      morseSmaleManifold[vert[1]], morseSmaleManifold[vert[2]]
                    );
                    foundFirst = true;
                  }
                }
              }

              sepPcsLocal.push_back(std::make_tuple(
                sparseID, sepTriangles[0], sepTriangles[1], false));
            } else {
              caseDataLocal.insert(caseDataLocal.end(),
                {lookupIndex, lookupIndex});
            }
          } else {
            caseDataLocal.push_back(lookupIndex);
          }
        }
      }
    }

    #pragma omp critical
    {
      trianglePos.insert(trianglePos.end(),
        trianglePosLocal.begin(), trianglePosLocal.end());
      caseData.insert(caseData.end(),
        caseDataLocal.begin(), caseDataLocal.end());
    }

    #pragma omp critical
    {
      for(auto& it_sepPcs1 : sepPcsLocal){
        const long long &sparseID = std::get<0>(it_sepPcs1);
        if(pieces1separatrices.find(sparseID) == pieces1separatrices.end()) { // new sparseID
          pieces1separatrices.insert({sparseID,
            std::vector<std::tuple<SimplexId, SimplexId, bool>>()}
          );
        }
        pieces1separatrices[sparseID].push_back(std::make_tuple(
          std::get<1>(it_sepPcs1),
          std::get<2>(it_sepPcs1),
          std::get<3>(it_sepPcs1))
        );
      }
    }

    #pragma omp critical
    {
      std::copy(saddleCandidatesLocal.begin(), saddleCandidatesLocal.end(),
        std::inserter(*saddleCandidates, saddleCandidates->end()));
    }
  }
} // #if TTK_OMPENMP
  this->printMsg("Computed 2-Separatrices 3D",
                 1, localTimer.getElapsedTime(), this->threadNumber_);

  return 0;
}

template <typename dataType, typename triangulationType>
int ttk::MorseSmaleComplexQuasi::compute1SeparatricesInner_3D(
  std::unordered_map<long long, std::vector<
    std::tuple<SimplexId, SimplexId, bool>>> &pieces1separatrices,
  std::vector<mscq::Separatrix> &separatrices1,
  std::unordered_set<SimplexId> &borderSeeds,
  const SimplexId *const morseSmaleManifold,
  const triangulationType &triangulation) const {

  ttk::Timer localTimer;

  // print the progress of the current subprocedure (currently 0%)
  this->printMsg("Computing inner 1-Separatrices 3D",
                 0, // progress form 0-1
                 0, // elapsed time so far
                 this->threadNumber_, ttk::debug::LineMode::REPLACE);
  
  const auto numTris = triangulation.getNumberOfTriangles();

if(!db_omp) {
  for (auto& it_seps: pieces1separatrices) { // put the sepPieces in sequence
    auto sepVectorPtr = it_seps.second;
    const size_t vecSize = sepVectorPtr.size();
    SimplexId maxIndex = 0;

    // write all Triangles into a doubly Linked List
    mscq::DoublyLinkedElement linkedTris[vecSize * 2];
    std::unordered_map<SimplexId, SimplexId> triIdToArrIndex;
    {
      for(std::size_t sepI = 0; sepI < vecSize; ++sepI) {
        const SimplexId firstTri = std::get<0>(sepVectorPtr[sepI]);
        const SimplexId secondTri = std::get<1>(sepVectorPtr[sepI]);
        const bool is4LabelSplit = std::get<2>(sepVectorPtr[sepI]);

        if(is4LabelSplit) {
          if(triIdToArrIndex.find(firstTri) == triIdToArrIndex.end()) {
            SimplexId newIndex = maxIndex++;
            triIdToArrIndex.insert({firstTri, newIndex});
            linkedTris[newIndex] =
              mscq::DoublyLinkedElement((secondTri + numTris));

            newIndex = maxIndex++;
            triIdToArrIndex.insert({secondTri + numTris, newIndex});
            linkedTris[newIndex] = mscq::DoublyLinkedElement(firstTri);

          } else {
            linkedTris[triIdToArrIndex[firstTri]].insert(secondTri  + numTris);

            SimplexId newIndex = maxIndex++;
            triIdToArrIndex.insert({(secondTri + numTris), newIndex});
            linkedTris[newIndex] = mscq::DoublyLinkedElement(firstTri);
          }
        } else {
          if(triIdToArrIndex.find(firstTri) == triIdToArrIndex.end()) {
            SimplexId newIndex = maxIndex++;
            triIdToArrIndex.insert({firstTri, newIndex});
            linkedTris[newIndex] = mscq::DoublyLinkedElement(secondTri);

          } else {
            linkedTris[triIdToArrIndex[firstTri]].insert(secondTri);
          }

          if(triIdToArrIndex.find(secondTri) == triIdToArrIndex.end()) {
            SimplexId newIndex = maxIndex++;
            triIdToArrIndex.insert({secondTri, newIndex});
            linkedTris[newIndex] = mscq::DoublyLinkedElement(firstTri);

          } else {
            linkedTris[triIdToArrIndex[secondTri]].insert(firstTri);
          }
        }
      }
    }

    std::unordered_map<SimplexId, bool> startIds;
    for(auto& it: triIdToArrIndex) {
      if(!linkedTris[it.second].hasTwoNeighbors) {
        startIds.insert({it.first, false});
      }
    }

    for(const auto& it_starts: startIds) { // detect paths
      if(it_starts.second) { // Is simplex already part of a separatrix?
        continue;
      }
      // write separatrix
      const SimplexId startTri = it_starts.first;
      mscq::Separatrix newSep(startTri);
      
      SimplexId previousId = startTri;
      SimplexId nextId = linkedTris[triIdToArrIndex[startTri]].neighbors_[0];
      linkedTris[triIdToArrIndex[previousId]].visited = true;
      newSep.geometry_.push_back(nextId);

      while(linkedTris[triIdToArrIndex[nextId]].hasTwoNeighbors &&
        !linkedTris[triIdToArrIndex[nextId]].visited) {
        SimplexId tempPreviousId = nextId;
        nextId = linkedTris[triIdToArrIndex[nextId]].next(previousId);
        previousId = tempPreviousId;
        newSep.geometry_.push_back(nextId);
        linkedTris[triIdToArrIndex[previousId]].visited = true;
      }
      
      startIds[nextId] = true;

      if(newSep.geometry_.front() >= numTris ||
        newSep.geometry_.back() >= numTris) {
        
        if(newSep.geometry_.front() < newSep.geometry_.back()) {
          std::reverse(newSep.geometry_.begin(), newSep.geometry_.end());
        }

        if(newSep.geometry_.back() < numTris) { // ends on boundary
          borderSeeds.insert(newSep.geometry_.back());
          newSep.endsOnBoundary_ = true;
        }

        separatrices1.push_back(newSep);
      } else {
        borderSeeds.insert(newSep.geometry_.front());
        borderSeeds.insert(newSep.geometry_.back());
      }
    }
  }
} else {
  #pragma omp parallel num_threads(this->threadNumber_)
  {
    std::vector<mscq::Separatrix> separatrices1Local;
    std::vector<SimplexId> borderSeedsLocal;

    #pragma omp for schedule(dynamic, 4) nowait
    for(size_t b = 0; b < pieces1separatrices.bucket_count(); b++) {
      for(auto it_s = pieces1separatrices.begin(b);
        it_s != pieces1separatrices.end(b); it_s++){
        auto sepVectorPtr = it_s->second;
        const size_t vecSize = sepVectorPtr.size();
        SimplexId maxIndex = 0;

        // write all Triangles into a doubly Linked List
        mscq::DoublyLinkedElement linkedTris[vecSize * 2];
        std::unordered_map<SimplexId, SimplexId> triIdToArrIndex;
        {
          for(std::size_t sepI = 0; sepI < vecSize; ++sepI) {
            const SimplexId firstTri = std::get<0>(sepVectorPtr[sepI]);
            const SimplexId secondTri = std::get<1>(sepVectorPtr[sepI]);
            const bool is4LabelSplit = std::get<2>(sepVectorPtr[sepI]);

            if(is4LabelSplit) {
              if(triIdToArrIndex.find(firstTri) == triIdToArrIndex.end()) {
                SimplexId newIndex = maxIndex++;
                triIdToArrIndex.insert({firstTri, newIndex});
                linkedTris[newIndex] =
                  mscq::DoublyLinkedElement((secondTri + numTris));

                newIndex = maxIndex++;
                triIdToArrIndex.insert({secondTri + numTris, newIndex});
                linkedTris[newIndex] = mscq::DoublyLinkedElement(firstTri);

              } else {
                linkedTris[triIdToArrIndex[firstTri]].insert(
                  secondTri  + numTris
                );

                SimplexId newIndex = maxIndex++;
                triIdToArrIndex.insert({(secondTri + numTris), newIndex});
                linkedTris[newIndex] = mscq::DoublyLinkedElement(firstTri);
              }
            } else {
              if(triIdToArrIndex.find(firstTri) == triIdToArrIndex.end()) {
                SimplexId newIndex = maxIndex++;
                triIdToArrIndex.insert({firstTri, newIndex});
                linkedTris[newIndex] = mscq::DoublyLinkedElement(secondTri);

              } else {
                linkedTris[triIdToArrIndex[firstTri]].insert(secondTri);
              }

              if(triIdToArrIndex.find(secondTri) == triIdToArrIndex.end()) {
                SimplexId newIndex = maxIndex++;
                triIdToArrIndex.insert({secondTri, newIndex});
                linkedTris[newIndex] = mscq::DoublyLinkedElement(firstTri);

              } else {
                linkedTris[triIdToArrIndex[secondTri]].insert(firstTri);
              }
            }
          }
        }

        std::unordered_map<SimplexId, bool> startIds;
        for(auto& it: triIdToArrIndex) {
          if(!linkedTris[it.second].hasTwoNeighbors) {
            startIds.insert({it.first, false});
          }
        }

        for (const auto& it_starts: startIds) { // detect paths
          if(it_starts.second) { // Is simplex already part of a separatrix?
            continue;
          }
          // write separatrix
          const SimplexId startTri = it_starts.first;
          mscq::Separatrix newSep(startTri);

          SimplexId previousId = startTri;
          SimplexId nextId =
            linkedTris[triIdToArrIndex[startTri]].neighbors_[0];
          linkedTris[triIdToArrIndex[previousId]].visited = true;
          newSep.geometry_.push_back(nextId);

          while(linkedTris[triIdToArrIndex[nextId]].hasTwoNeighbors &&
            !linkedTris[triIdToArrIndex[nextId]].visited) {
            SimplexId tempPreviousId = nextId;
            nextId = linkedTris[triIdToArrIndex[nextId]].next(previousId);
            previousId = tempPreviousId;
            newSep.geometry_.push_back(nextId);
            linkedTris[triIdToArrIndex[previousId]].visited = true;
          }

          startIds[nextId] = true;

          if(newSep.geometry_.front() >= numTris ||
            newSep.geometry_.back() >= numTris) {
            
            if(newSep.geometry_.front() < newSep.geometry_.back()) {
              std::reverse(newSep.geometry_.begin(), newSep.geometry_.end());
            }

            if(newSep.geometry_.back() < numTris) { // ends on boundary
              borderSeedsLocal.push_back(newSep.geometry_.back());
              newSep.endsOnBoundary_ = true;
            }

            separatrices1Local.push_back(newSep);
          } else {
            borderSeedsLocal.push_back(newSep.geometry_.front());
            borderSeedsLocal.push_back(newSep.geometry_.back());
          }
        }

        for (auto& it: triIdToArrIndex) { // detect cycles
          if(!linkedTris[it.second].visited) {
            mscq::Separatrix newSep(std::get<0>(it));

            SimplexId previousId = std::get<0>(it);
            SimplexId nextId = linkedTris[std::get<1>(it)].neighbors_[0];
            linkedTris[triIdToArrIndex[previousId]].visited = true;
            newSep.geometry_.push_back(nextId);

            while(!linkedTris[triIdToArrIndex[nextId]].visited) {
              SimplexId tempPreviousId = nextId;
              nextId = linkedTris[triIdToArrIndex[nextId]].next(previousId);
              previousId = tempPreviousId;
              newSep.geometry_.push_back(nextId);
              linkedTris[triIdToArrIndex[previousId]].visited = true;
            }

            separatrices1Local.push_back(newSep);
          }
        }
      }
    }
    #pragma omp critical
    {
      separatrices1.insert(separatrices1.end(),
        separatrices1Local.begin(), separatrices1Local.end());
    }
    #pragma omp critical
    {
      std::copy(borderSeedsLocal.begin(), borderSeedsLocal.end(),
        std::inserter(borderSeeds, borderSeeds.end()));
    }
  }
} // TTK_ENABLE_OPENMP

  this->printMsg("Computed inner 1-Separatrices 3D",
                 1, localTimer.getElapsedTime(), this->threadNumber_);

  return 0;
}

template <typename triangulationType>
int ttk::MorseSmaleComplexQuasi::compute1SeparatricesBorder_3D(
  std::unordered_set<SimplexId> &borderSeeds,
  std::vector<mscq::Separatrix> &separatrices1,
  const SimplexId *const morseSmaleManifold,
  const triangulationType &triangulation) const {

  ttk::Timer localTimer;

  // print the progress of the current subprocedure (currently 0%)
  this->printMsg("Computing border 1-Separatrices 3D",
                 0, // progress form 0-1
                 localTimer.getElapsedTime(), // elapsed time so far
                 this->threadNumber_, ttk::debug::LineMode::REPLACE);

  std::unordered_map<long long,
    std::vector<std::tuple<SimplexId, SimplexId>>> borderPathSeeds;

  for (const auto& it: borderSeeds) {
    SimplexId triEdges[3];
    triangulation.getTriangleEdge(it, 0, triEdges[0]);
    triangulation.getTriangleEdge(it, 1, triEdges[1]);
    triangulation.getTriangleEdge(it, 2, triEdges[2]);

    SimplexId triEdgeVertices[3][2];
    triangulation.getEdgeVertex(triEdges[0], 0, triEdgeVertices[0][0]);
    triangulation.getEdgeVertex(triEdges[0], 1, triEdgeVertices[0][1]);
    triangulation.getEdgeVertex(triEdges[1], 0, triEdgeVertices[1][0]);
    triangulation.getEdgeVertex(triEdges[1], 1, triEdgeVertices[1][1]);
    triangulation.getEdgeVertex(triEdges[2], 0, triEdgeVertices[2][0]);
    triangulation.getEdgeVertex(triEdges[2], 1, triEdgeVertices[2][1]);

    long long sparseIds[3];
    sparseIds[0] = getSparseId(
      morseSmaleManifold[triEdgeVertices[0][0]],
      morseSmaleManifold[triEdgeVertices[0][1]]
    );
    sparseIds[1] = getSparseId(
      morseSmaleManifold[triEdgeVertices[1][0]],
      morseSmaleManifold[triEdgeVertices[1][1]]
    );
    sparseIds[2] = getSparseId(
      morseSmaleManifold[triEdgeVertices[2][0]],
      morseSmaleManifold[triEdgeVertices[2][1]]
    );

    if(borderPathSeeds.find(sparseIds[0]) == borderPathSeeds.end()) {
      borderPathSeeds.insert({
        sparseIds[0], std::vector<std::tuple<SimplexId, SimplexId>>()}
      );
    }
    if(borderPathSeeds.find(sparseIds[1]) == borderPathSeeds.end()) {
      borderPathSeeds.insert({
        sparseIds[1], std::vector<std::tuple<SimplexId, SimplexId>>()}
      );
    }
    if(borderPathSeeds.find(sparseIds[2]) == borderPathSeeds.end()) {
      borderPathSeeds.insert({
        sparseIds[2], std::vector<std::tuple<SimplexId, SimplexId>>()}
      );
    }

    borderPathSeeds[sparseIds[0]].push_back(
      std::make_tuple(it, triEdges[0])
    );
    borderPathSeeds[sparseIds[1]].push_back(
      std::make_tuple(it, triEdges[1])
    );
    borderPathSeeds[sparseIds[2]].push_back(
      std::make_tuple(it, triEdges[2])
    );
  }

  for (auto& it: borderPathSeeds) {
    if(it.second.size() <= 2) { // one border path with this sparseID
      mscq::Separatrix newSep;
      newSep.onBoundary_ = true;
      createBorderPath_3D(
        newSep,
        std::get<0>(it.second.front()),
        std::get<1>(it.second.front()),
        morseSmaleManifold,
        triangulation
      );
      separatrices1.push_back(newSep);
    } else { // multiple border paths with this sparseID
      for(size_t bps = 0; bps < it.second.size(); bps++) {

        // create the path first
        mscq::Separatrix newSep;
        newSep.onBoundary_ = true;
        createBorderPath_3D(
          newSep,
          std::get<0>(it.second[bps]),
          std::get<1>(it.second[bps]),
          morseSmaleManifold,
          triangulation
        );
        separatrices1.push_back(newSep);

        //remove the path from the list
        for(size_t bps2 = 0; bps2 < it.second.size(); bps2++) {
          if(newSep.geometry_.back() == std::get<0>(it.second[bps2]) &&
            newSep.geometry_[newSep.geometry_.size()- 2] == 
            std::get<1>(it.second[bps2])) {
            it.second.erase(it.second.begin() + bps2);
            continue;
          }
        }
      }
    }
  }

  this->printMsg("Computed border 1-Separatrices 3D",
                 1, localTimer.getElapsedTime(), this->threadNumber_);

  return 1;
}

template <typename triangulationType>
int ttk::MorseSmaleComplexQuasi::createBorderPath_3D(
  mscq::Separatrix &separatrix,
  SimplexId startTriangle,
  SimplexId startEdge,
  const SimplexId *const morseSmaleManifold,
  const triangulationType &triangulation) const {
  
  separatrix.geometry_.push_back(startTriangle);
  separatrix.geometry_.push_back(startEdge);

  SimplexId currTriangle;
  SimplexId currEdge;
  SimplexId prevTriangle = startTriangle;
  SimplexId prevEdge = startEdge;

  while(true) {
    int triangleNumber = triangulation.getEdgeTriangleNumber(prevEdge);
    for(int i = 0; i < triangleNumber; ++i) {
      triangulation.getEdgeTriangle(prevEdge, i, currTriangle);
      if(triangulation.isTriangleOnBoundary(currTriangle) &&
        currTriangle != prevTriangle) {
        break;
      }
    }

    SimplexId triVerts[3];
    triangulation.getTriangleVertex(currTriangle, 0, triVerts[0]);
    triangulation.getTriangleVertex(currTriangle, 1, triVerts[1]);
    triangulation.getTriangleVertex(currTriangle, 2, triVerts[2]);

    // 3 labels on triangle -> end of path reached
    if(morseSmaleManifold[triVerts[0]] != morseSmaleManifold[triVerts[1]] &&
      morseSmaleManifold[triVerts[0]] != morseSmaleManifold[triVerts[2]] &&
      morseSmaleManifold[triVerts[1]] != morseSmaleManifold[triVerts[2]]) {
      break;
    }

    // find next edge
    for(int i = 0; i < 3; ++i) {
      triangulation.getTriangleEdge(currTriangle, i, currEdge);

      if(currEdge != prevEdge) {
        SimplexId edgeVerts[2];
        triangulation.getEdgeVertex(currEdge, 0, edgeVerts[0]);
        triangulation.getEdgeVertex(currEdge, 1, edgeVerts[1]);

        if(morseSmaleManifold[edgeVerts[0]] !=
          morseSmaleManifold[edgeVerts[1]]) {
          break;
        }
      }
    }

    separatrix.geometry_.push_back(currEdge);

    prevEdge = currEdge;
    prevTriangle = currTriangle;
  }

  separatrix.geometry_.push_back(currTriangle);

  return 1;
}

template <typename dataType, typename triangulationType>
int ttk::MorseSmaleComplexQuasi::computeSeparatrices_3D_fast(
  std::vector<std::array<float, 9>> &trianglePos,    
  std::vector<SimplexId> &caseData,
  const SimplexId *const morseSmaleManifold,
  const triangulationType &triangulation) const {
  ttk::Timer localTimer;

  this->printMsg(ttk::debug::Separator::L1);
  this->printMsg("Computing 2-Separatrices 3D Fast",
                 0, // progress form 0-1
                 0, // elapsed time so far
                 this->threadNumber_, ttk::debug::LineMode::REPLACE);

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

  const SimplexId numTetra = triangulation.getNumberOfCells();

//#ifndef TTK_ENABLE_OPENMP
if(!db_omp) {
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
    unsigned char index2 = (msmV0 == msmV2) ? 0x00 :       // 0
                           (msmV1 == msmV2) ? 0x04 : 0x08; // 1 : 2
    unsigned char index3 = (msmV0 == msmV3) ? 0x00 :       // 0
                           (msmV1 == msmV3) ? 0x01 :       // 1
                           (msmV2 == msmV3) ? 0x02 : 0x03; // 2 : 3

    unsigned char lookupIndex = index1 | index2 | index3;

    if(tetraederLookupFast[lookupIndex]) { // <= 2 label on tetraeder
      float vertPos[4][3];
      triangulation.getVertexPoint(
        vertices[0], vertPos[0][0], vertPos[0][1], vertPos[0][2]);
      triangulation.getVertexPoint(
        vertices[1], vertPos[1][0], vertPos[1][1], vertPos[1][2]);
      triangulation.getVertexPoint(
        vertices[2], vertPos[2][0], vertPos[2][1], vertPos[2][2]);
      triangulation.getVertexPoint(
        vertices[3], vertPos[3][0], vertPos[3][1], vertPos[3][2]);
      
      trianglePos.push_back({
        vertPos[0][0], vertPos[0][1], vertPos[0][2], 
        vertPos[1][0], vertPos[1][1], vertPos[1][2], 
        vertPos[2][0], vertPos[2][1], vertPos[2][2]
      });
      trianglePos.push_back({
        vertPos[0][0], vertPos[0][1], vertPos[0][2], 
        vertPos[1][0], vertPos[1][1], vertPos[1][2], 
        vertPos[3][0], vertPos[3][1], vertPos[3][2]
      });
      trianglePos.push_back({
        vertPos[0][0], vertPos[0][1], vertPos[0][2], 
        vertPos[2][0], vertPos[2][1], vertPos[2][2], 
        vertPos[3][0], vertPos[3][1], vertPos[3][2]
      });
      trianglePos.push_back({
        vertPos[1][0], vertPos[1][1], vertPos[1][2], 
        vertPos[2][0], vertPos[2][1], vertPos[2][2], 
        vertPos[3][0], vertPos[3][1], vertPos[3][2]
      });

      caseData.insert(caseData.end(),{
        lookupIndex, lookupIndex, lookupIndex, lookupIndex
      });
    }
  }
//#else // TTK_ENABLE_OPENMP
} else {
  const size_t numCacheLineEntries =
    hardware_destructive_interference_size / sizeof(SimplexId);

  #pragma omp parallel num_threads(this->threadNumber_)
  {
    std::vector<std::array<float, 9>> trianglePosLocal;
    std::vector<SimplexId> caseDataLocal;

    #pragma omp for schedule(dynamic, numCacheLineEntries) nowait
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
      unsigned char index2 = (msmV0 == msmV2) ? 0x00 :       // 0
                             (msmV1 == msmV2) ? 0x04 : 0x08; // 1 : 2
      unsigned char index3 = (msmV0 == msmV3) ? 0x00 :       // 0
                             (msmV1 == msmV3) ? 0x01 :       // 1
                             (msmV2 == msmV3) ? 0x02 : 0x03; // 2 : 3

      unsigned char lookupIndex = index1 | index2 | index3;

      if(tetraederLookupFast[lookupIndex]) { // <= 2 label on tetraeder
        float vertPos[4][3];
        triangulation.getVertexPoint(
          vertices[0], vertPos[0][0], vertPos[0][1], vertPos[0][2]);
        triangulation.getVertexPoint(
          vertices[1], vertPos[1][0], vertPos[1][1], vertPos[1][2]);
        triangulation.getVertexPoint(
          vertices[2], vertPos[2][0], vertPos[2][1], vertPos[2][2]);
        triangulation.getVertexPoint(
          vertices[3], vertPos[3][0], vertPos[3][1], vertPos[3][2]);

        trianglePosLocal.push_back({
          vertPos[0][0], vertPos[0][1], vertPos[0][2], 
          vertPos[1][0], vertPos[1][1], vertPos[1][2], 
          vertPos[2][0], vertPos[2][1], vertPos[2][2]
        });
        trianglePosLocal.push_back({
          vertPos[0][0], vertPos[0][1], vertPos[0][2], 
          vertPos[1][0], vertPos[1][1], vertPos[1][2], 
          vertPos[3][0], vertPos[3][1], vertPos[3][2]
        });
        trianglePosLocal.push_back({
          vertPos[0][0], vertPos[0][1], vertPos[0][2], 
          vertPos[2][0], vertPos[2][1], vertPos[2][2], 
          vertPos[3][0], vertPos[3][1], vertPos[3][2]
        });
        trianglePosLocal.push_back({
          vertPos[1][0], vertPos[1][1], vertPos[1][2], 
          vertPos[2][0], vertPos[2][1], vertPos[2][2], 
          vertPos[3][0], vertPos[3][1], vertPos[3][2]
        });

        caseDataLocal.insert(caseDataLocal.end(),{
          lookupIndex, lookupIndex, lookupIndex, lookupIndex
        });
      }
    }

    #pragma omp critical
    {
      trianglePos.insert(trianglePos.end(),
        trianglePosLocal.begin(), trianglePosLocal.end());
      caseData.insert(caseData.end(),
        caseDataLocal.begin(), caseDataLocal.end());
    }
  }
} // TTK_ENABLE_OPENMP
  this->printMsg("Computed Separatrices 3D FAST",
                 1, localTimer.getElapsedTime(), this->threadNumber_);

  this->printMsg("#2-sep triangles: " + std::to_string(trianglePos.size()),
                 1, localTimer.getElapsedTime(), this->threadNumber_);

  this->printMsg("Completed",
                 1, localTimer.getElapsedTime(), this->threadNumber_);
  this->printMsg(ttk::debug::Separator::L1); // horizontal '=' separator

  return 0;
}

template <typename triangulationType>
int ttk::MorseSmaleComplexQuasi::setSeparatrices1_3D(
  const std::vector<mscq::Separatrix> &separatrices,
  const triangulationType &triangulation) const {

  ttk::Timer localTimer;

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

  this->printMsg("Writing 1-separatrices",
                  1,
                  localTimer.getElapsedTime(),
                  this->threadNumber_,
                  ttk::debug::LineMode::REPLACE);

  const auto numTris = triangulation.getNumberOfTriangles();
  // number of separatrices
  size_t numSep = separatrices.size();
  size_t npoints = 0;
  size_t numCells = 0;

  // count total number of points and cells, flatten geometryId loops
  for(size_t i = 0; i < numSep; ++i) {
    const auto &sep = separatrices[i];

    npoints += sep.geometry_.size() + 1; // +0?
    numCells += sep.geometry_.size(); // -1?
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

    std::array<float, 3> pt{};

    if(!sep.onBoundary_) { // inner separatrix
      triangulation.getTetraIncenter(
        (sep.geometry_.front() - numTris), pt.data());

      points[3 * currPId + 0] = pt[0];
      points[3 * currPId + 1] = pt[1];
      points[3 * currPId + 2] = pt[2];

      currPId += 1;

      size_t geomSize = sep.geometry_.size() - 1;

      for(size_t j = 1; j < geomSize; ++j) {
        triangulation.getTriangleIncenter(sep.geometry_[j], pt.data());

        points[3 * currPId + 0] = pt[0];
        points[3 * currPId + 1] = pt[1];
        points[3 * currPId + 2] = pt[2];

        cellsConn[2 * currCId + 0] = currPId - 1;
        cellsConn[2 * currCId + 1] = currPId;

        currPId += 1;
        currCId += 1;
      }

      if(sep.endsOnBoundary_){
        triangulation.getTriangleIncenter(sep.geometry_[geomSize], pt.data());
      } else {
        triangulation.getTetraIncenter(
          sep.geometry_[geomSize] - numTris, pt.data()
        );
      }

      points[3 * currPId + 0] = pt[0];
      points[3 * currPId + 1] = pt[1];
      points[3 * currPId + 2] = pt[2];

      cellsConn[2 * currCId + 0] = currPId - 1;
      cellsConn[2 * currCId + 1] = currPId;

      currPId += 1;
      currCId += 1;

      if(outputSeparatrices1_cells_isOnBoundary_ != nullptr)
        (*outputSeparatrices1_cells_isOnBoundary_)[i] = false;

    } else { // boundary separatrix

      triangulation.getTriangleIncenter(sep.geometry_.front(), pt.data());

      points[3 * currPId + 0] = pt[0];
      points[3 * currPId + 1] = pt[1];
      points[3 * currPId + 2] = pt[2];

      currPId += 1;

      size_t geomSize = sep.geometry_.size() - 1;

      for(size_t j = 1; j < geomSize; ++j) {
        triangulation.getEdgeIncenter(sep.geometry_[j], pt.data());

        points[3 * currPId + 0] = pt[0];
        points[3 * currPId + 1] = pt[1];
        points[3 * currPId + 2] = pt[2];

        cellsConn[2 * currCId + 0] = currPId - 1;
        cellsConn[2 * currCId + 1] = currPId;

        currPId += 1;
        currCId += 1;
      }

      triangulation.getTriangleIncenter(sep.geometry_.back(), pt.data());

      points[3 * currPId + 0] = pt[0];
      points[3 * currPId + 1] = pt[1];
      points[3 * currPId + 2] = pt[2];

      cellsConn[2 * currCId + 0] = currPId - 1;
      cellsConn[2 * currCId + 1] = currPId;

      currPId += 1;
      currCId += 1;

      if(outputSeparatrices1_cells_isOnBoundary_ != nullptr)
        (*outputSeparatrices1_cells_isOnBoundary_)[i] = true;
    }
  }

  // update pointers
  *outputSeparatrices1_numberOfPoints_ = npoints;
  *outputSeparatrices1_numberOfCells_ = numCells;

  this->printMsg("Wrote 1-separatrices",
                 1, localTimer.getElapsedTime(), this->threadNumber_);

  return 0;
}

int ttk::MorseSmaleComplexQuasi::extractMSCPaths(
  std::vector<mscq::CriticalPoint> &criticalPoints,    
  std::vector<mscq::Separatrix> &separatrices) const {

  // "vertices" of the "graph" with a list of all its neighbors
  std::unordered_map<SimplexId, std::set<SimplexId>> vertexNeighbors;
  
  // "edges" of the "graph", all edges that are connected to a "vertex"
  std::unordered_map<SimplexId, std::vector<size_t>> vertexEdges;

  const size_t sepsSize = separatrices.size();
  
  for(size_t sep = 0; sep < sepsSize; sep++) {
    const SimplexId& start = separatrices[sep].geometry_.front();
    const SimplexId& end = separatrices[sep].geometry_.back();

    // Check if id already exists
    if(vertexNeighbors.find(start) == vertexNeighbors.end()) {
      vertexNeighbors.insert({start, std::set<SimplexId>()});
    }
    if(vertexNeighbors.find(end) == vertexNeighbors.end()) {
      vertexNeighbors.insert({end, std::set<SimplexId>()});
    }
    if(vertexEdges.find(start) == vertexEdges.end()) {
      vertexEdges.insert({start, std::vector<size_t>()});
    }
    if(vertexEdges.find(end) == vertexEdges.end()) {
      vertexEdges.insert({end, std::vector<size_t>()});
    }

    // insert verticies and edges
    vertexNeighbors[start].insert(end);
    vertexNeighbors[end].insert(start);
    vertexEdges[start].push_back(sep);
    vertexEdges[end].push_back(sep);
  }

  return 1;
}

int ttk::MorseSmaleComplexQuasi::setSeparatrices2_3D(
  std::vector<std::array<float, 9>> &trianglePos,
  std::vector<SimplexId> &caseData) const {

  ttk::Timer localTimer;

  this->printMsg("Writing 2-separatrices",
                  1,
                  localTimer.getElapsedTime(),
                  this->threadNumber_,
                  ttk::debug::LineMode::REPLACE);

  const int numTriangles = trianglePos.size();

  size_t npoints = numTriangles * 3;
  size_t numCells = numTriangles;

  // resize arrays
  outputSeparatrices2_points_->resize(3 * npoints);
  auto &points = *outputSeparatrices2_points_;
  outputSeparatrices2_cells_connectivity_->resize(3 * numCells);
  auto &cellsConn = *outputSeparatrices2_cells_connectivity_;
  outputSeparatrices2_cells_cases->resize(numTriangles);
  auto &cellCase = *outputSeparatrices2_cells_cases;
  cellCase = caseData;

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for schedule(static) if(numTriangles > 1000000)
#endif
  for(int tri = 0; tri < numTriangles; ++tri) {
    auto &triPos = trianglePos[tri];
    points[9 * tri + 0] = triPos[0];
    points[9 * tri + 1] = triPos[1];
    points[9 * tri + 2] = triPos[2];

    points[9 * tri + 3] = triPos[3];
    points[9 * tri + 4] = triPos[4];
    points[9 * tri + 5] = triPos[5];

    points[9 * tri + 6] = triPos[6];
    points[9 * tri + 7] = triPos[7];
    points[9 * tri + 8] = triPos[8];

    cellsConn[3 * tri]     = 3 * tri;
    cellsConn[3 * tri + 1] = 3 * tri + 1;
    cellsConn[3 * tri + 2] = 3 * tri + 2;
  }

  // update pointers
  *outputSeparatrices2_numberOfPoints_ = npoints;
  *outputSeparatrices2_numberOfCells_ = numCells;

  this->printMsg("Wrote 2-separatrices",
                 1, localTimer.getElapsedTime(), this->threadNumber_);

  return 1;
}

template <typename dataType, typename triangulationType>
int ttk::MorseSmaleComplexQuasi::findSaddlesFromCadidates(
  std::vector<mscq::CriticalPoint> &criticalPoints, 
  const std::unordered_set<SimplexId> &saddleCandidates,
  const triangulationType &triangulation) const {

  ttk::Timer localTimer;

  const int dim = triangulation.getDimensionality();

  const dataType *inputField = static_cast<const dataType *>(inputOrderField_);

  int numSaddles = 0;

  this->printMsg("Searching for saddles",
                 1, localTimer.getElapsedTime(), this->threadNumber_,
                 ttk::debug::LineMode::REPLACE);
  
if(!db_omp) {
  for(const auto& candidate: saddleCandidates) {
    SimplexId lowerComponents, upperComponents;
    std::tie(lowerComponents, upperComponents)
      = getNumberOfLowerUpperComponents(
      candidate, inputField, triangulation
    );

    if(dim == 3) {
      if(lowerComponents > 1 && upperComponents == 1) {
        criticalPoints.push_back(mscq::CriticalPoint(candidate, 1));
        numSaddles += 1;
      } else if(lowerComponents == 1 && upperComponents > 1) {
        criticalPoints.push_back(mscq::CriticalPoint(candidate, 2));
        numSaddles += 1;
      }
    } else {
      if((lowerComponents > 1 && upperComponents == 1) ||
        (lowerComponents == 1 && upperComponents > 1) ||
        (lowerComponents > 1 && upperComponents > 1)) {
        criticalPoints.push_back(mscq::CriticalPoint(candidate, 1));
        numSaddles += 1;
      }
    }
    
  }
} else {
  #pragma omp parallel
  {
    #pragma omp single
    {
      for(auto it = saddleCandidates.begin();
        it != saddleCandidates.end(); it++) {
        #pragma omp task
        {
          SimplexId lowerComponents, upperComponents;
          bool foundSaddle = false;
          int criticalityIndex = -1;

          std::tie(lowerComponents, upperComponents)
            = getNumberOfLowerUpperComponents(
            *it, inputField, triangulation
          );

          if(dim == 3) {
            if(lowerComponents > 1 && upperComponents == 1) {
              foundSaddle = true;
              criticalityIndex = 1;
            } else if(lowerComponents == 1 && upperComponents > 1) {
              foundSaddle = true;
              criticalityIndex = 2;
            }
          } else {
            if((lowerComponents > 1 && upperComponents == 1) ||
              (lowerComponents == 1 && upperComponents > 1) ||
              (lowerComponents > 1 && upperComponents > 1)) {
              foundSaddle = true;
              criticalityIndex = 1;
            }
          }

          if(foundSaddle) {
            #pragma omp atomic update
            numSaddles += 1;

            #pragma omp critical
            {
              criticalPoints.push_back(
                mscq::CriticalPoint(*it, criticalityIndex)
              );
            }
          }
        }
      }
    }
  }
}

  this->printMsg("Found saddles" ,
                 1, localTimer.getElapsedTime(), this->threadNumber_);

  return 1;
}

template <typename dataType, typename triangulationType>
int ttk::MorseSmaleComplexQuasi::setCriticalPoints(
  std::vector<mscq::CriticalPoint> &criticalPoints,
  const triangulationType &triangulation) const {

  ttk::Timer localTimer;

  this->printMsg("Writing saddles",
                 1, localTimer.getElapsedTime(), this->threadNumber_,
                 ttk::debug::LineMode::REPLACE);
  
  float x, y, z;

  for(const auto crit : criticalPoints) {
    triangulation.getVertexPoint(crit.id_, x, y, z);
    outputCriticalPoints_points_->push_back({x, y, z});
    outputCriticalPoints_points_cellDimensions_->push_back(crit.index_);
    outputCriticalPoints_points_cellIds_->push_back(crit.id_);
  }

  this->printMsg("Wrote saddles",
                 1, localTimer.getElapsedTime(), this->threadNumber_);

  return 1;
}

template <typename dataType, typename triangulationType>
std::pair<ttk::SimplexId, ttk::SimplexId>
  ttk::MorseSmaleComplexQuasi::getNumberOfLowerUpperComponents(
    const SimplexId vertexId,
    const dataType * const offsets,
    const triangulationType &triangulation) const {

  SimplexId neighborNumber = triangulation.getVertexNeighborNumber(vertexId);
  std::vector<SimplexId> lowerNeighbors, upperNeighbors;

  for(SimplexId i = 0; i < neighborNumber; i++) {
    SimplexId neighborId = 0;
    triangulation.getVertexNeighbor(vertexId, i, neighborId);

    if(offsets[neighborId] < offsets[vertexId]) {
      lowerNeighbors.push_back(neighborId);
    }

    // upper link
    if(offsets[neighborId] > offsets[vertexId]) {
      upperNeighbors.push_back(neighborId);
    }
  }

  // shortcut, if min or max do not construct the complete star
  if(lowerNeighbors.empty()) {
    // minimum
    return std::make_pair(0, 1);
  }

  if(upperNeighbors.empty()) {
    // maximum
    return std::make_pair(1, 0);
  }

  // now do the actual work
  std::vector<UnionFind> lowerSeeds(lowerNeighbors.size());
  std::vector<UnionFind *> lowerList(lowerNeighbors.size());
  std::vector<UnionFind> upperSeeds(upperNeighbors.size());
  std::vector<UnionFind *> upperList(upperNeighbors.size());

  for(SimplexId i = 0; i < (SimplexId)lowerSeeds.size(); i++) {
    lowerList[i] = &(lowerSeeds[i]);
  }
  for(SimplexId i = 0; i < (SimplexId)upperSeeds.size(); i++) {
    upperList[i] = &(upperSeeds[i]);
  }

  SimplexId vertexStarSize = triangulation.getVertexStarNumber(vertexId);

  for(SimplexId i = 0; i < vertexStarSize; i++) {
    SimplexId cellId = 0;
    triangulation.getVertexStar(vertexId, i, cellId);

    SimplexId cellSize = triangulation.getCellVertexNumber(cellId);
    for(SimplexId j = 0; j < cellSize; j++) {
      SimplexId neighborId0 = -1;
      triangulation.getCellVertex(cellId, j, neighborId0);

      if(neighborId0 != vertexId) {
        // we are on the link

        bool lower0 = offsets[neighborId0] < offsets[vertexId];

        // connect it to everybody except himself and vertexId
        for(SimplexId k = j + 1; k < cellSize; k++) {

          SimplexId neighborId1 = -1;
          triangulation.getCellVertex(cellId, k, neighborId1);

          if((neighborId1 != neighborId0) && (neighborId1 != vertexId)) {

            bool lower1 = offsets[neighborId1] < offsets[vertexId];

            std::vector<SimplexId> *neighbors = &lowerNeighbors;
            std::vector<UnionFind *> *seeds = &lowerList;

            if(!lower0) {
              neighbors = &upperNeighbors;
              seeds = &upperList;
            }

            if(lower0 == lower1) {
              // connect their union-find sets!
              SimplexId lowerId0 = -1, lowerId1 = -1;
              for(SimplexId l = 0; l < (SimplexId)neighbors->size(); l++) {
                if((*neighbors)[l] == neighborId0) {
                  lowerId0 = l;
                }
                if((*neighbors)[l] == neighborId1) {
                  lowerId1 = l;
                }
              }
              if((lowerId0 != -1) && (lowerId1 != -1)) {
                (*seeds)[lowerId0] = UnionFind::makeUnion(
                  (*seeds)[lowerId0], (*seeds)[lowerId1]);
                (*seeds)[lowerId1] = (*seeds)[lowerId0];
              }
            }
          }
        }
      }
    }
  }

  // let's remove duplicates now

  // update the UF if necessary
  for(SimplexId i = 0; i < (SimplexId)lowerList.size(); i++)
    lowerList[i] = lowerList[i]->find();
  for(SimplexId i = 0; i < (SimplexId)upperList.size(); i++)
    upperList[i] = upperList[i]->find();

  std::vector<UnionFind *>::iterator it;
  std::sort(lowerList.begin(), lowerList.end());
  it = unique(lowerList.begin(), lowerList.end());
  lowerList.resize(distance(lowerList.begin(), it));

  std::sort(upperList.begin(), upperList.end());
  it = unique(upperList.begin(), upperList.end());
  upperList.resize(distance(upperList.begin(), it));

  if(debugLevel_ >= (int)(debug::Priority::VERBOSE)) {
    printMsg("Vertex #" + std::to_string(vertexId)
               + ": lowerLink-#CC=" + std::to_string(lowerList.size())
               + " upperLink-#CC=" + std::to_string(upperList.size()),
             debug::Priority::VERBOSE);
  }

  return std::make_pair(lowerList.size(), upperList.size());
}