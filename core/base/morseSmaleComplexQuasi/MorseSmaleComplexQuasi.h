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
    struct Simplex {
      explicit Simplex() = default;

      explicit Simplex(const SimplexId id, const int index)
        : id_{id}, index_{index} {
      }

      bool operator==(const  Simplex &other) const {
        return id_ == other.id_ && index_ == other.index_;
      }

      bool operator<(const Simplex &other) const {
        return index_ < other.index_   ? true
               : other.index_ < index_ ? false
               : id_ < other.id_       ? true
                                       : false;
      }

      bool operator>(const Simplex &other) const {
        return index_ < other.index_   ? false
               : other.index_ < index_ ? true
               : id_ < other.id_       ? false
                                       : true;
      }

      SimplexId id_{-1};
      int index_{-1};
    };

    struct Separatrix {
      // default :
      explicit Separatrix()
        : geometry_{}, dims_{}, onBoundary_{false},
        startAtHighDimSimplex_{true}, endAtHighDimSimplex_{true} {
      }

      // initialization with one segment :
      explicit Separatrix(const SimplexId segmentGeometry, const int dim) {
        geometry_.push_back(segmentGeometry);
        dims_.push_back(dim);
      }

      explicit Separatrix(const Separatrix &separatrix)
        : geometry_{separatrix.geometry_}, dims_{std::move(separatrix.dims_)},
          incline_{separatrix.incline_},
          critTypeStart_{separatrix.critTypeStart_},
          critTypeEnd_{separatrix.critTypeEnd_},
          onBoundary_{separatrix.onBoundary_},
          startAtHighDimSimplex_{separatrix.startAtHighDimSimplex_},
          endAtHighDimSimplex_{separatrix.endAtHighDimSimplex_},
          critIdStart_{separatrix.critIdStart_},
          critIdEnd_{separatrix.critIdEnd_}, lowerNeighbors_{std::move(
                                               separatrix.lowerNeighbors_)},
          regionIds_{std::move(separatrix.regionIds_)} {
      }

      explicit Separatrix(Separatrix &&separatrix) noexcept
        : geometry_{std::move(separatrix.geometry_)},
          dims_{std::move(separatrix.dims_)}, incline_{separatrix.incline_},
          critTypeStart_{separatrix.critTypeStart_},
          critTypeEnd_{separatrix.critTypeEnd_},
          onBoundary_{separatrix.onBoundary_},
          startAtHighDimSimplex_{separatrix.startAtHighDimSimplex_},
          endAtHighDimSimplex_{separatrix.endAtHighDimSimplex_},
          critIdStart_{separatrix.critIdStart_},
          critIdEnd_{separatrix.critIdEnd_},
          lowerNeighbors_{std::move(separatrix.lowerNeighbors_)},
          regionIds_{std::move(separatrix.regionIds_)} {
      }

      Separatrix &operator=(Separatrix &&separatrix) noexcept {
        geometry_ = std::move(separatrix.geometry_);
        dims_ = std::move(separatrix.dims_);
        incline_ = separatrix.incline_;
        critTypeStart_ = separatrix.critTypeStart_;
        critTypeEnd_ = separatrix.critTypeEnd_;
        onBoundary_ = separatrix.onBoundary_;
        startAtHighDimSimplex_ = separatrix.startAtHighDimSimplex_;
        endAtHighDimSimplex_ = separatrix.endAtHighDimSimplex_;
        critIdStart_ = separatrix.critIdStart_;
        critIdEnd_ = separatrix.critIdEnd_;
        lowerNeighbors_ = std::move(separatrix.lowerNeighbors_);
        regionIds_ = std::move(separatrix.regionIds_);

        return *this;
      }

      Separatrix &operator=(const Separatrix &separatrix) noexcept {
        geometry_ = std::move(separatrix.geometry_);
        dims_ = std::move(separatrix.dims_);
        incline_ = separatrix.incline_;
        critTypeStart_ = separatrix.critTypeStart_;
        critTypeEnd_ = separatrix.critTypeEnd_;
        onBoundary_ = separatrix.onBoundary_;
        startAtHighDimSimplex_ = separatrix.startAtHighDimSimplex_;
        endAtHighDimSimplex_ = separatrix.endAtHighDimSimplex_;
        critIdStart_ = separatrix.critIdStart_;
        critIdEnd_ = separatrix.critIdEnd_;
        lowerNeighbors_ = std::move(separatrix.lowerNeighbors_);
        regionIds_ = std::move(separatrix.regionIds_);

        return *this;
      }

      void insertGeometry(SimplexId s, int d) {
        geometry_.push_back(s);
        dims_.push_back(d);
      }

      void reverse() {
        std::reverse(geometry_.begin(), geometry_.end());
        std::reverse(dims_.begin(), dims_.end());

        int temp = critTypeStart_;
        critTypeStart_ = critTypeEnd_;
        critTypeEnd_ = temp;
      }

      size_t length(){
        return geometry_.size();
      }

      void getSimplexAt(size_t pos, mscq::Simplex &simp) {
        if(pos < geometry_.size()) {
          simp.id_ = geometry_[pos];
          simp.index_ = dims_[pos];
        }
      }

      SimplexId Id_{-1};

      /**
       * Container of ids. Each id addresses a separate
       * container corresponding to a dense representation
       * of the geometry (i.e. separatricesGeometry).
       */
      std::vector<SimplexId> geometry_;

      std::vector<int> dims_;

      int incline_{0};

      int critTypeStart_{-1};

      int critTypeEnd_{-1};

      /**
       * Flag indicating that the separatix is located on the boundary
       */
      bool onBoundary_{false};

      bool startAtHighDimSimplex_{true};
      bool endAtHighDimSimplex_{true};

      SimplexId critIdStart_{-1};
      SimplexId critIdEnd_{-1};

      SimplexId lowIdStart_{-1};
      SimplexId lowIdEnd_{-1};

      std::vector<mscq::Separatrix> lowerNeighbors_;
      std::set<SimplexId> regionIds_;
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

    struct SaddleSaddlePath {
      SaddleSaddlePath() {
      }

      SaddleSaddlePath(mscq::Simplex svi,
                       mscq::Simplex li,
                       mscq::Separatrix *sep,
                       bool reversed) {
        saddleVertId_ = svi;
        lastId_ = li;
        path_.push_back(sep);
        reversed_.push_back(reversed);
        visitedIds_.insert(sep->Id_);
        length_ = sep->length();
      }

      SaddleSaddlePath(const mscq::Simplex &svi,
                       const mscq::Simplex &li,
                       ttk::mscq::Separatrix &sep,
                       bool reversed) {
        saddleVertId_ = svi;
        lastId_ = li;
        path_.push_back(&sep);
        reversed_.push_back(reversed);
        visitedIds_.insert(sep.Id_);
        length_ = sep.length();
      }

      explicit SaddleSaddlePath(const SaddleSaddlePath &ssp)
        : saddleVertId_{ssp.saddleVertId_}, lastId_{ssp.lastId_},
          path_{ssp.path_}, reversed_{std::move(ssp.reversed_)},
          visitedIds_{ssp.visitedIds_}, length_{ssp.length_} {
      }

      explicit SaddleSaddlePath(SaddleSaddlePath &&ssp) noexcept
        : saddleVertId_{ssp.saddleVertId_}, lastId_{ssp.lastId_},
          path_{std::move(ssp.path_)}, reversed_{std::move(ssp.reversed_)},
          visitedIds_{std::move(ssp.visitedIds_)}, length_{ssp.length_} {
      }

      SaddleSaddlePath &operator=(SaddleSaddlePath &&ssp) noexcept {
        saddleVertId_ = ssp.saddleVertId_;
        lastId_ = ssp.lastId_;
        path_ = std::move(ssp.path_);
        reversed_ = std::move(ssp.reversed_);
        visitedIds_ = std::move(ssp.visitedIds_);
        length_ = ssp.length_;

        return *this;
      }

      SaddleSaddlePath &operator=(const SaddleSaddlePath &ssp) noexcept {
        saddleVertId_ = ssp.saddleVertId_;
        lastId_ = ssp.lastId_;
        path_ = std::move(ssp.path_);
        reversed_ = std::move(ssp.reversed_);
        visitedIds_ = std::move(ssp.visitedIds_);
        length_ = ssp.length_;

        return *this;
      }

      void insertPathPiece(mscq::Separatrix *sep,
                           mscq::Simplex newLastId,
                           bool reversed) {
        lastId_ = newLastId;
        path_.push_back(sep);
        reversed_.push_back(reversed);
        visitedIds_.insert(sep->Id_);
        length_ += sep->length() - 1;
      }

      mscq::Simplex saddleVertId_;
      mscq::Simplex lastId_;
      std::vector<mscq::Separatrix *> path_;
      std::vector<bool> reversed_;
      std::set<SimplexId> visitedIds_;
      SimplexId length_{0};
    };
  } // namespace mscq

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
      int success = 0;
      success += triangulation->preconditionVertexNeighbors();
      success += triangulation->preconditionVertexEdges();
      success += triangulation->preconditionVertexTriangles();
      success += triangulation->preconditionVertexStars();

      success += triangulation->preconditionEdgeTriangles();

      success += triangulation->preconditionTriangleEdges();
      success += triangulation->preconditionTriangleStars();

      success += triangulation->preconditionCellTriangles();

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
      SimplexId *const neighbor,
      std::vector<mscq::Simplex> &criticalPoints,
      SimplexId &numberOfMaxima,
      const triangulationType &triangulation,
      const bool ascending) const;

    template <typename triangulationType>
    int computeFinalSegmentation(
      const SimplexId numMaxima,
      const std::vector<mscq::Simplex> &criticalPoints,
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

    template <typename triangulationType>
    int split1SeparatricesAtCrit_3D(
      std::vector<mscq::Simplex> &criticalPoints,
      std::vector<mscq::Separatrix> &separatrices1,
      const SimplexId *const morseSmaleManifold,
      const triangulationType &triangulation) const;

    template <typename dataType, typename triangulationType>
    int findSaddleSaddleConnectors(
      std::vector<mscq::Simplex> &criticalPoints,
      std::vector<mscq::Separatrix> &separatrices,
      std::vector<mscq::SaddleSaddlePath> &ssPaths,
      const triangulationType &triangulation) const;

    template <typename dataType, typename triangulationType>
    int computeSeparatrices_3D_fast(
      std::vector<std::array<float, 9>> &trianglePos,
      std::vector<SimplexId> &caseData,
      const SimplexId *const morseSmaleManifold,
      const triangulationType &triangulation) const;

    template <typename triangulationType>
    int createBorderPath_3D(mscq::Separatrix &separatrix,
                            SimplexId startTriangle,
                            SimplexId startEdge,
                            const SimplexId *const morseSmaleManifold,
                            const triangulationType &triangulation) const;

    template <typename dataType, typename triangulationType>
    int computeIntegralLines(
      std::vector<mscq::Separatrix> &plIntLines,
      const SimplexId *const ascendingNeighbor,
      const SimplexId *const descendingNeighbor,
      const SimplexId *const ascendingManifold,
      const SimplexId *const descendingManifold,
      const std::vector<mscq::Simplex> &criticalPoints,
      const triangulationType &triangulation) const;

    template <typename triangulationType>
    int setSeparatrices1_3D(const std::vector<mscq::Separatrix> &separatrices,
                            const triangulationType &triangulation) const;

    template<typename triangulationType>
    int setSaddleSaddleConnectors_3D(
        const std::vector<mscq::SaddleSaddlePath> &ssp,
        const triangulationType &triangulation) const;

    int setSeparatrices2_3D(std::vector<std::array<float, 9>> &trianglePos,
                            std::vector<SimplexId> &caseData) const;

    template <typename dataType, typename triangulationType>
    int findSaddlesFromCadidates(
      std::vector<mscq::Simplex> &criticalPoints,
      const std::unordered_set<SimplexId> &saddleCandidates,
      const triangulationType &triangulation) const;

    template <typename dataType, typename triangulationType>
    int setCriticalPoints(std::vector<mscq::Simplex> &criticalPoints,
                          const triangulationType &triangulation) const;

    template <typename dataType, typename triangulationType>
    std::pair<ttk::SimplexId, ttk::SimplexId> getNumberOfLowerUpperComponents(
      const SimplexId vertexId,
      const dataType *const offsets,
      const triangulationType &triangulation) const;

    template <typename dataType, typename triangulationType>
    inline double
      getAvgSimplexOffset(const SimplexId simplexId,
                         const int dim,
                         const dataType *const offsets,
                         const triangulationType &triangulation) const {
      switch(dim) {
        case 0 :
          return getAvgVertexOffset(simplexId, offsets, triangulation);
        case 1:
          return getAvgEdgeOffset(simplexId, offsets, triangulation);
        case 2:
          return getAvgTriOffset(simplexId, offsets, triangulation);
        case 3:
          return getAvgTetOffset(simplexId, offsets, triangulation);
        default:
          return -1.0;
      }
    }

    template <typename dataType, typename triangulationType>
    inline double
      getAvgVertexOffset(const SimplexId vertexId,
                         const dataType *const offsets,
                         const triangulationType &triangulation) const {
      return offsets[vertexId];
    }

    template <typename dataType, typename triangulationType>
    inline double
      getAvgEdgeOffset(const SimplexId edgeId,
                       const dataType *const offsets,
                       const triangulationType &triangulation) const {
      SimplexId verts[2];
      triangulation.getEdgeVertex(edgeId, 0, verts[0]);
      triangulation.getEdgeVertex(edgeId, 1, verts[1]);

      return (offsets[verts[0]] + offsets[verts[1]]) / 2.0;
    }

    template <typename dataType, typename triangulationType>
    inline double
      getAvgTriOffset(const SimplexId triId,
                      const dataType *const offsets,
                      const triangulationType &triangulation) const {
      SimplexId verts[3];
      triangulation.getEdgeVertex(triId, 0, verts[0]);
      triangulation.getEdgeVertex(triId, 1, verts[1]);
      triangulation.getEdgeVertex(triId, 2, verts[2]);

      return (offsets[verts[0]] + offsets[verts[1]] + offsets[verts[2]]) / 3.0;
    }

    template <typename dataType, typename triangulationType>
    inline double
      getAvgTetOffset(const SimplexId tetId,
                      const dataType *const offsets,
                      const triangulationType &triangulation) const {
      SimplexId verts[4];
      triangulation.getEdgeVertex(tetId, 0, verts[0]);
      triangulation.getEdgeVertex(tetId, 1, verts[1]);
      triangulation.getEdgeVertex(tetId, 2, verts[2]);
      triangulation.getEdgeVertex(tetId, 3, verts[3]);

      return (offsets[verts[0]] + offsets[verts[1]] + offsets[verts[2]]
              + offsets[verts[3]])
             / 4.0;
    }

    template <typename dataType, typename triangulationType>
    SimplexId getSmallesVertexOfVertex(
      const SimplexId vertexId,
      const dataType *const offsets,
      const triangulationType &triangulation) const;

    template <typename dataType, typename triangulationType>
    SimplexId getSmallesVertexOfEdge(
      const SimplexId edgeId,
      const dataType *const offsets,
      const triangulationType &triangulation) const;

    template <typename dataType, typename triangulationType>
    SimplexId getSmallesVertexOfTri(
      const SimplexId triId,
      const dataType *const offsets,
      const triangulationType &triangulation) const;

    template <typename dataType, typename triangulationType>
    SimplexId getSmallesVertexOfTet(
      const SimplexId tetId,
      const dataType *const offsets,
      const triangulationType &triangulation) const;

    template <typename triangulationType>
    inline int getEdgeIncenter(SimplexId vertexId0,
                               SimplexId vertexId1,
                               float incenter[3],
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

  const SimplexId nV = triangulation.getNumberOfVertices();

  SimplexId *ascendingNeighbor  = (SimplexId*) malloc(nV * sizeof(SimplexId));
  SimplexId *descendingNeighbor = (SimplexId*) malloc(nV * sizeof(SimplexId));

  std::vector<mscq::Simplex> criticalPoints;
  {
    Timer tmp;

    const int dim = triangulation.getDimensionality();

    SimplexId numberOfMaxima{0};
    SimplexId numberOfMinima{0};

    if(ascendingManifold)
      computeManifold<dataType, triangulationType>(
                                ascendingManifold, ascendingNeighbor,
                                criticalPoints, numberOfMinima,
                                triangulation, true);

    if(descendingManifold)                            
      computeManifold<dataType, triangulationType>(
                                descendingManifold, descendingNeighbor,
                                criticalPoints, numberOfMaxima,
                                triangulation, false);

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

          computeSeparatrices1Inner_2D<dataType, triangulationType>(
            saddleCandidates, separatrices1, sepPieces, triConnectors,
            sepManifold, triangulation
          );

          computeSeparatrices1Border_2D<dataType, triangulationType>(
            saddleCandidates, separatrices1, sepManifold, triangulation
          );

          findSaddlesFromCadidates<dataType, triangulationType>(
            criticalPoints, saddleCandidates, triangulation
          );

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
          std::vector<mscq::Separatrix> plIntLines;
      
          if(fastSeparatrices) {
            computeSeparatrices_3D_fast<dataType, triangulationType>(
                                      trianglePos, caseData,
                                      sepManifold, triangulation);

            setSeparatrices2_3D(trianglePos, caseData);
          } else {
            std::vector<mscq::SaddleSaddlePath> ssp;

            compute2Separatrices_3D<dataType, triangulationType>(
              trianglePos, caseData, &saddleCandidates, pieces1separatrices,
              sepManifold, triangulation);

            findSaddlesFromCadidates<dataType, triangulationType>(
              criticalPoints, saddleCandidates, triangulation);

            compute1SeparatricesInner_3D<dataType, triangulationType>(
              pieces1separatrices, separatrices1, borderSeeds,
              sepManifold, triangulation);

            compute1SeparatricesBorder_3D<triangulationType>(
              borderSeeds, separatrices1, sepManifold, triangulation);

            computeIntegralLines<dataType, triangulationType>(
              plIntLines, ascendingNeighbor, descendingNeighbor,
              ascendingManifold, descendingManifold, criticalPoints,
              triangulation);

            split1SeparatricesAtCrit_3D<triangulationType>(
              criticalPoints, separatrices1, sepManifold, triangulation);

            findSaddleSaddleConnectors<dataType, triangulationType>(
              criticalPoints, separatrices1, ssp, triangulation);

            setSaddleSaddleConnectors_3D<triangulationType>(ssp, triangulation);

            setCriticalPoints<dataType, triangulationType>(
              criticalPoints, triangulation);

            //setIntegralLines_3D<triangulationType>(
            //  plIntLines, triangulation);

            //setSeparatrices1_3D<triangulationType>(plIntLines, triangulation);

            //setSeparatrices1_3D<triangulationType>(
            //  separatrices1, triangulation);

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
  SimplexId *const manifold,
  SimplexId *const neighbor,
  std::vector<mscq::Simplex> &criticalPoints,
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

      SimplexId &dmi = manifold[i];
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

    memcpy(neighbor, manifold, nVertices * sizeof(SimplexId));

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
        SimplexId &vo = manifold[v];

        // compress path
        vo = manifold[vo];

        // check if not fully compressed
        if(vo != manifold[vo]) {
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
      manifold[i] = maxima[manifold[i]];
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
      #pragma omp for schedule(static)
      for(SimplexId i = 0; i < nVertices; i++) {
        SimplexId neighborId;
        SimplexId numNeighbors =
          triangulation.getVertexNeighborNumber(i);
        bool hasBiggerNeighbor = false;
        SimplexId &dmi = manifold[i];

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
        memcpy(neighbor, manifold, nVertices * sizeof(SimplexId));
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

        #pragma omp for schedule(guided)
        for(SimplexId i = 0; i < nActiveVertices; i++) {
          SimplexId v = activeVertices->at(i);
          SimplexId &vo = manifold[v];

          vo = manifold[vo];

          // check if fully compressed
          if(vo != manifold[vo]) {
            lActiveVertices.push_back(v);
          }
        }
        #pragma omp barrier

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
        manifold[i] = maxima.at(manifold[i]);
      }
    }
} //#endif // TTK_ENABLE_OPENMP

    const int index =
      ascending ? 0 : triangulation.getDimensionality();
    //const int dim = triangulation.getDimensionality();

    for (auto& it: maxima) {
      criticalPoints.push_back(mscq::Simplex(it.first, index));
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
  const std::vector<mscq::Simplex> &criticalPoints,
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

  this->printMsg("Computing 1-separatrix pieces",
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

  this->printMsg("Computed 1-separatrix pieces",
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
          newSep.insertGeometry(std::get<1>(triConnVect[ec]), 2);
          triConnectors[sparseID].erase( triConnectors[sparseID].begin() + ec);
          startsOnBoundary = false;
          break;
        }
      }
      
      SimplexId previousId = startEdge;
      newSep.insertGeometry(startEdge, 1);
      SimplexId nextId = linkedEdges[edgeIdToArrIndex[startEdge]].neighbors_[0];
      linkedEdges[edgeIdToArrIndex[previousId]].visited = true;
      newSep.insertGeometry(nextId, 1);

      while(linkedEdges[edgeIdToArrIndex[nextId]].hasTwoNeighbors &&
        !linkedEdges[edgeIdToArrIndex[nextId]].visited) {
        SimplexId tempPreviousId = nextId;
        nextId = linkedEdges[edgeIdToArrIndex[nextId]].next(previousId);
        previousId = tempPreviousId;
        newSep.insertGeometry(nextId, 1);
        linkedEdges[edgeIdToArrIndex[previousId]].visited = true;
      }

      // Find the ending Triangle, if it exists
      for(size_t ec = 0; ec < triConnVectSize; ++ec) {
        if(std::get<0>(triConnVect[ec]) == nextId) {
          newSep.insertGeometry(std::get<1>(triConnVect[ec]), 2);
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

        newSep.endAtHighDimSimplex_ = false;
      }
      if(startsOnBoundary) {
        SimplexId &sEdge = newSep.geometry_.front();

        SimplexId vert0, vert1;
        triangulation.getEdgeVertex(sEdge, 0, vert0);
        triangulation.getEdgeVertex(sEdge, 1, vert1);

        saddleCandidates.insert(vert0);
        saddleCandidates.insert(vert1);

        newSep.endAtHighDimSimplex_ = false;

        newSep.reverse();
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
            newSep.insertGeometry(splitTriangle0, 2);
            newSep.insertGeometry(splitEdge0, 1);
            newSep.insertGeometry(splitTriangle1, 2);
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

  int (triangulationType::*getPosFuncs[3])(SimplexId, float *) const = {
    &triangulationType::getVertexCoordinates,
    &triangulationType::getEdgeIncenter,
    &triangulationType::getTriangleIncenter
  };

  SimplexId currPId = 0, currCId = 0;
  for(size_t i = 0; i < numSep; ++i) {
    const auto &sep = separatrices[i];

    std::array<float, 3> pt{};

    (triangulation.*getPosFuncs[sep.dims_[0]])(sep.geometry_[0], pt.data());

    points[3 * currPId + 0] = pt[0];
    points[3 * currPId + 1] = pt[1];
    points[3 * currPId + 2] = pt[2];

    currPId += 1;

    for(size_t j = 1; j < sep.geometry_.size(); ++j) {
      (triangulation.*getPosFuncs[sep.dims_[j]])(sep.geometry_[j], pt.data());

      points[3 * currPId + 0] = pt[0];
      points[3 * currPId + 1] = pt[1];
      points[3 * currPId + 2] = pt[2];

      cellsConn[2 * currCId + 0] = currPId - 1;
      cellsConn[2 * currCId + 1] = currPId;

      currPId += 1;
      currCId += 1;
    }

    if(outputSeparatrices1_cells_isOnBoundary_ != nullptr)
      (*outputSeparatrices1_cells_isOnBoundary_)[i] = false;
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
  const bool findSaddles = saddleCandidates != nullptr;

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
      if(findSaddles) {
        saddleCandidates->insert(vertices[0]);
        saddleCandidates->insert(vertices[1]);
        saddleCandidates->insert(vertices[2]);
        saddleCandidates->insert(vertices[3]);
      }

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
  #pragma omp parallel num_threads(this->threadNumber_)
  {
    std::vector<std::array<float, 9>> trianglePosLocal;
    std::vector<SimplexId> caseDataLocal;
    std::vector<std::tuple<long long, SimplexId, SimplexId, bool>> sepPcsLocal;
    std::vector<SimplexId> saddleCandidatesLocal;

    #pragma omp for schedule(static) nowait
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
        if(findSaddles) {
          saddleCandidatesLocal.push_back(vertices[0]);
          saddleCandidatesLocal.push_back(vertices[1]);
          saddleCandidatesLocal.push_back(vertices[2]);
          saddleCandidatesLocal.push_back(vertices[3]);
        }

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
  #pragma omp parallel num_threads(this->threadNumber_)
  {
    std::vector<std::array<float, 9>> trianglePosLocal;
    std::vector<SimplexId> caseDataLocal;

    #pragma omp for schedule(static) nowait
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
      mscq::Separatrix newSep(startTri, 2);
      
      SimplexId previousId = startTri;
      SimplexId nextId = linkedTris[triIdToArrIndex[startTri]].neighbors_[0];
      linkedTris[triIdToArrIndex[previousId]].visited = true;
      newSep.insertGeometry(nextId, 2);

      while(linkedTris[triIdToArrIndex[nextId]].hasTwoNeighbors &&
        !linkedTris[triIdToArrIndex[nextId]].visited) {
        SimplexId tempPreviousId = nextId;
        nextId = linkedTris[triIdToArrIndex[nextId]].next(previousId);
        previousId = tempPreviousId;
        newSep.insertGeometry(nextId, 2);
        linkedTris[triIdToArrIndex[previousId]].visited = true;
      }
      
      startIds[nextId] = true;

      if(newSep.geometry_.front() >= numTris
         || newSep.geometry_.back() >= numTris) {

        if(newSep.geometry_.front() < newSep.geometry_.back()) {
          newSep.reverse();
        }

        if(newSep.geometry_.back() < numTris) { // ends on boundary
          borderSeeds.insert(newSep.geometry_.back());
          newSep.endAtHighDimSimplex_ = false;
        } else {
          // reverse hack for identifying tetrahedra
          newSep.geometry_.back() = newSep.geometry_.back() - numTris;
          newSep.dims_.back() += 1;
          newSep.endAtHighDimSimplex_ = true;
        }

        // reverse hack for identifying tetrahedra
        newSep.geometry_.front() = newSep.geometry_.front() - numTris;
        newSep.dims_.front() += 1;
        newSep.startAtHighDimSimplex_ = true;

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
          mscq::Separatrix newSep(startTri, 2);

          SimplexId previousId = startTri;
          SimplexId nextId =
            linkedTris[triIdToArrIndex[startTri]].neighbors_[0];
          linkedTris[triIdToArrIndex[previousId]].visited = true;
          newSep.insertGeometry(nextId, 2);

          while(linkedTris[triIdToArrIndex[nextId]].hasTwoNeighbors &&
            !linkedTris[triIdToArrIndex[nextId]].visited) {
            SimplexId tempPreviousId = nextId;
            nextId = linkedTris[triIdToArrIndex[nextId]].next(previousId);
            previousId = tempPreviousId;
            newSep.insertGeometry(nextId, 2);
            linkedTris[triIdToArrIndex[previousId]].visited = true;
          }

          startIds[nextId] = true;

          if(newSep.geometry_.front() >= numTris ||
            newSep.geometry_.back() >= numTris) {
            
            if(newSep.geometry_.front() < newSep.geometry_.back()) {
              newSep.reverse();
            }

            if(newSep.geometry_.back() < numTris) { // ends on boundary
              borderSeedsLocal.push_back(newSep.geometry_.back());
              newSep.endAtHighDimSimplex_ = false;
            } else {
              // reverse hack for identifying tetrahedra
              newSep.geometry_.back() = newSep.geometry_.back() - numTris;
              newSep.dims_.back() += 1;
              newSep.endAtHighDimSimplex_ = true;
            }
            
            // reverse hack for identifying tetrahedra
            newSep.geometry_.front() = newSep.geometry_.front() - numTris;
            newSep.dims_.front() += 1;
            newSep.startAtHighDimSimplex_ = true;

            separatrices1Local.push_back(newSep);
          } else {
            borderSeedsLocal.push_back(newSep.geometry_.front());
            borderSeedsLocal.push_back(newSep.geometry_.back());
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
            newSep.geometry_[newSep.geometry_.size() - 2] == 
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
int ttk::MorseSmaleComplexQuasi::split1SeparatricesAtCrit_3D(
  std::vector<mscq::Simplex> &criticalPoints,
  std::vector<mscq::Separatrix> &separatrices1,
  const SimplexId *const morseSmaleManifold,
  const triangulationType &triangulation) const {

  ttk::Timer localTimer;
  this->printMsg("Splitting 1-Separatrices at crits",
                 0, // progress form 0-1
                 localTimer.getElapsedTime(), // elapsed time so far
                 this->threadNumber_, ttk::debug::LineMode::REPLACE);

  std::unordered_map<SimplexId, int> critIdToIndex;
  critIdToIndex.reserve(criticalPoints.size());

  std::unordered_map<SimplexId, std::set<SimplexId>> vertCritMap;
  std::unordered_map<SimplexId, std::set<SimplexId>> edgeCritMap;
  std::unordered_map<SimplexId, std::set<SimplexId>> triCritMap;
  std::unordered_map<SimplexId, std::set<SimplexId>> tetCritMap;

  // Create lookup maps for tets, tris and edges that contain critical points
  for(const auto& c: criticalPoints) {
    critIdToIndex.insert({c.id_, c.index_});

    vertCritMap.insert({c.id_, std::set<SimplexId>()});
    vertCritMap[c.id_].insert(c.id_);

    const SimplexId numTris = triangulation.getVertexTriangleNumber(c.id_);
    SimplexId triId;

    for (SimplexId tri = 0; tri < numTris; ++tri) {
      triangulation.getVertexTriangle(c.id_, tri, triId);

      if(triCritMap.find(triId) == triCritMap.end()) {
        triCritMap.insert({triId, std::set<SimplexId>()});
      }
      triCritMap[triId].insert(c.id_);
    }

    if(triangulation.isVertexOnBoundary(c.id_)) {
      const SimplexId numEdges = triangulation.getVertexEdgeNumber(c.id_);
      SimplexId edgeId;

      for (SimplexId edge = 0; edge < numEdges; ++edge) {
        triangulation.getVertexEdge(c.id_, edge, edgeId);

        if(edgeCritMap.find(edgeId) == edgeCritMap.end()) {
          edgeCritMap.insert({edgeId, std::set<SimplexId>()});
        }
        edgeCritMap[edgeId].insert(c.id_);
      }
    } else {
      const SimplexId numTets = triangulation.getVertexStarNumber(c.id_);
      SimplexId tetId;

      for (SimplexId tet = 0; tet < numTets; ++tet) {
        triangulation.getVertexEdge(c.id_, tet, tetId);

        if(tetCritMap.find(tetId) == tetCritMap.end()) {
          tetCritMap.insert({tetId, std::set<SimplexId>()});
        }
        tetCritMap[tetId].insert(c.id_);
      }
    }
  }

  std::unordered_map<SimplexId, std::set<SimplexId>>* critMaps[4] = {
    &vertCritMap, &edgeCritMap, &triCritMap, &tetCritMap
  };

  std::vector<mscq::Separatrix> oldSeparatrices(separatrices1);
  separatrices1.clear();
  separatrices1.reserve(2 * oldSeparatrices.size());

  const size_t numSep = oldSeparatrices.size();

  for(size_t i = 0; i < numSep; ++i) {
    const auto &sep = oldSeparatrices[i];

    std::vector<size_t> splitPoints;
    std::vector<std::set<SimplexId>> splitPointCriticalPoints;
    std::set<SimplexId> foundCriticalPoints;

    // first simplex contains a critical point?
    auto critIdFront = critMaps[sep.dims_[0]]->find(sep.geometry_[0]);
    if(critIdFront != critMaps[sep.dims_[0]]->end()) {
      splitPoints.push_back(0);
      splitPointCriticalPoints.push_back(critIdFront->second);
      for(auto critE : critIdFront->second) {
        foundCriticalPoints.insert(critE);
      }
    }

    // find geometry parts containing a critical point
    const size_t geoSize = sep.geometry_.size() - 2;
    for(size_t geo = 2; geo < geoSize; ++geo) {
      auto critIdsMid = critMaps[sep.dims_[geo]]->find(sep.geometry_[geo]);
      if(critIdsMid != critMaps[sep.dims_[geo]]->end()) {
        std::set<SimplexId> newCritSet;
        for(auto c : critIdsMid->second) {
          if(foundCriticalPoints.find(c) ==
            foundCriticalPoints.end()) {
            foundCriticalPoints.insert(c);
            newCritSet.insert(c);
          }
        }

        if(!newCritSet.empty()) {
          splitPointCriticalPoints.push_back(newCritSet);
          splitPoints.push_back(geo);
        }
      }
    }

    // last simplex contains a critical point?
    auto critIdsBack = critMaps[sep.dims_.back()]->find(sep.geometry_.back());
    if(critIdsBack != critMaps[sep.dims_.back()]->end()) {
      std::set<SimplexId> newCritSet;
      for(auto c : critIdsBack->second) {
        if(foundCriticalPoints.find(sep.geometry_.back()) ==
          foundCriticalPoints.end()) {
          foundCriticalPoints.insert(c);
          newCritSet.insert(c);
        }
      }

      if(!newCritSet.empty()) {
        splitPointCriticalPoints.push_back(newCritSet);
        splitPoints.push_back(sep.geometry_.size() - 1);
      }
    }

    // create new Separatrices if critical points were found
    if(splitPoints.empty()) {
      separatrices1.push_back(sep);
    }
    else {
      std::vector<mscq::Separatrix> startSeps;
      size_t nextInsertion = 0;
      size_t splitPointIdx = 0;

      if(splitPoints[splitPointIdx] == 0) { // critical point at start simplex
        for(const auto &spcp : splitPointCriticalPoints[splitPointIdx]) {
          mscq::Separatrix newSep(spcp, 0);
          newSep.endAtHighDimSimplex_ = sep.endAtHighDimSimplex_;
          newSep.startAtHighDimSimplex_ = false;
          newSep.critIdStart_ = spcp;
          newSep.critTypeStart_ = critIdToIndex[spcp];
          startSeps.push_back(newSep);
        }
        nextInsertion = 1;
        splitPointIdx = 1;
      } else { // No critical point at start simplex
        mscq::Separatrix newSep;
        newSep.endAtHighDimSimplex_ = sep.endAtHighDimSimplex_;
        newSep.startAtHighDimSimplex_ = true;
        startSeps.push_back(newSep);
      }

      for(; splitPointIdx < splitPoints.size(); ++splitPointIdx) {
        for(size_t ss = 0; ss < startSeps.size(); ++ss) {
          for(const auto &spcp : splitPointCriticalPoints[splitPointIdx]) {
            mscq::Separatrix newSep(startSeps[ss]);
            if(splitPoints[splitPointIdx] < sep.geometry_.size() - 1) {
              newSep.endAtHighDimSimplex_ = false;
            }
            newSep.critIdEnd_ = spcp;
            newSep.critTypeEnd_ = critIdToIndex[spcp];
            if(splitPoints[splitPointIdx] > nextInsertion) {
              newSep.geometry_.insert(
                newSep.geometry_.end(), sep.geometry_.begin() + nextInsertion,
                sep.geometry_.begin() + (splitPoints[splitPointIdx]));
              newSep.dims_.insert(
                newSep.dims_.end(), sep.dims_.begin() + nextInsertion,
                sep.dims_.begin() + (splitPoints[splitPointIdx]));
            }
            newSep.insertGeometry(spcp, 0);

            newSep.onBoundary_ = sep.onBoundary_;
            separatrices1.push_back(newSep);
          }
        }

        startSeps.clear();

        for(const auto &spcp : splitPointCriticalPoints[splitPointIdx]) {
          mscq::Separatrix newSep(spcp, 0);
          newSep.startAtHighDimSimplex_ = false;
          newSep.onBoundary_ = sep.onBoundary_;
          newSep.critIdStart_ = spcp;
          newSep.critTypeStart_ = critIdToIndex[spcp];
          startSeps.push_back(newSep);
        }

        nextInsertion = splitPoints[splitPointIdx] + 1;
      }

      // last simplex is no critical point
      if(nextInsertion < sep.geometry_.size()) {
        for(size_t ss = 0; ss < startSeps.size(); ++ss) {
          if(nextInsertion < sep.geometry_.size()) {
            startSeps[ss].geometry_.insert(
              startSeps[ss].geometry_.end(),
              sep.geometry_.begin() + nextInsertion, sep.geometry_.end());
            startSeps[ss].dims_.insert(startSeps[ss].dims_.end(),
                                       sep.dims_.begin() + nextInsertion,
                                       sep.dims_.end());
          }
          startSeps[ss].onBoundary_ = sep.onBoundary_;
          startSeps[ss].endAtHighDimSimplex_ = sep.endAtHighDimSimplex_;
          separatrices1.push_back(startSeps[ss]);
        }
      }
    }
  }

  this->printMsg("Splitted 1-Separatrices at crits",
                 1, // progress form 0-1
                 localTimer.getElapsedTime(), // elapsed time so far
                 this->threadNumber_);

  return 1;
}

template <typename dataType, typename triangulationType>
SimplexId ttk::MorseSmaleComplexQuasi::getSmallesVertexOfVertex(
  const SimplexId vertexId,
  const dataType *const offsets,
  const triangulationType &triangulation) const {
  return vertexId;
}

template <typename dataType, typename triangulationType>
SimplexId ttk::MorseSmaleComplexQuasi::getSmallesVertexOfEdge(
  const SimplexId edgeId,
  const dataType *const offsets,
  const triangulationType &triangulation) const {
  SimplexId verts[2];
  triangulation.getEdgeVertex(edgeId, 0, verts[0]);
  triangulation.getEdgeVertex(edgeId, 1, verts[1]);

  return offsets[verts[0]] < offsets[verts[1]] ? verts[0] : verts[1];
}

template <typename dataType, typename triangulationType>
SimplexId ttk::MorseSmaleComplexQuasi::getSmallesVertexOfTri(
  const SimplexId triId,
  const dataType *const offsets,
  const triangulationType &triangulation) const {
  SimplexId verts[3];
  triangulation.getTriangleVertex(triId, 0, verts[0]);
  triangulation.getTriangleVertex(triId, 1, verts[1]);
  triangulation.getTriangleVertex(triId, 2, verts[2]);

  const dataType offs[3]
    = {offsets[verts[0]], offsets[verts[1]], offsets[verts[2]]};

  size_t index = std::min_element(offs, offs + 3) - offs;

  return verts[index];
}

template <typename dataType, typename triangulationType>
SimplexId ttk::MorseSmaleComplexQuasi::getSmallesVertexOfTet(
  const SimplexId tetId,
  const dataType *const offsets,
  const triangulationType &triangulation) const {
  SimplexId verts[4];
  triangulation.getCellVertex(tetId, 0, verts[0]);
  triangulation.getCellVertex(tetId, 1, verts[1]);
  triangulation.getCellVertex(tetId, 2, verts[2]);
  triangulation.getCellVertex(tetId, 3, verts[3]);

  const dataType offs[4] = {
    offsets[verts[0]], offsets[verts[1]], offsets[verts[2]], offsets[verts[3]]};

  size_t index = std::min_element(offs, offs + 4) - offs;

  return verts[index];
}

template <typename dataType, typename triangulationType>
int ttk::MorseSmaleComplexQuasi::findSaddleSaddleConnectors(
  std::vector<mscq::Simplex> &criticalPoints,
  std::vector<mscq::Separatrix> &separatrices,
  std::vector<mscq::SaddleSaddlePath> &ssPaths,
  const triangulationType &triangulation) const {

  // start a local timer for this subprocedure
  ttk::Timer localTimer;

  this->printMsg("Computing saddle saddle connectors",
                 0, // progress form 0-1
                 0, // elapsed time so far
                 this->threadNumber_);//, ttk::debug::LineMode::REPLACE);

  const dataType *offsets = static_cast<const dataType *>(inputOrderField_);

  // map vertex Ids to their correspoding critical point
  std::map<SimplexId, size_t> vertexToCritId;
  for(size_t i = 0; i < criticalPoints.size(); ++i) {
    vertexToCritId.insert({criticalPoints[i].id_, i});
  }

  // queue to store unfinished saddle saddle paths
  std::queue<mscq::SaddleSaddlePath> sspQueue;

  // each simplex maps to a vector of all connected Separatrices
  std::map<mscq::Simplex, std::vector<mscq::Separatrix*>> simpSepMap;

  const size_t sepsSize = separatrices.size();

  for(size_t sep = 0; sep < sepsSize; sep++) {
    mscq::Separatrix &sepRef = separatrices[sep];
    sepRef.Id_ = sep;

    // ascending/descending/neutral separatrix
    int incline = 0;

    double prevOffset = getAvgSimplexOffset(
      sepRef.geometry_[0], sepRef.dims_[0], offsets, triangulation);
    double currOffset;

    for(size_t g = 1; g < sepRef.geometry_.size(); ++g) {
      currOffset = getAvgSimplexOffset(
        sepRef.geometry_[g], sepRef.dims_[g], offsets, triangulation);

      incline
        += currOffset < prevOffset ? -1 : currOffset > prevOffset ? 1 : 0;

      prevOffset = currOffset;
    }

    if(incline == 0) {
      prevOffset = getAvgSimplexOffset(
        sepRef.geometry_.front(), sepRef.dims_.front(), offsets, triangulation);
      currOffset = getAvgSimplexOffset(
        sepRef.geometry_.back(), sepRef.dims_.back(), offsets, triangulation);

      incline += currOffset < prevOffset   ? -1
                        : currOffset > prevOffset ? 1
                                                  : 0;
    }

    sepRef.incline_ = incline < 0 ? -1 : incline > 0 ? 1 : 0;

    mscq::Simplex start, end;
    sepRef.getSimplexAt(0, start);
    sepRef.getSimplexAt(sepRef.geometry_.size() - 1, end);

    if(sepRef.critTypeStart_ != -1) {
      if(sepRef.critTypeStart_ == 2 && sepRef.incline_ <= 0) {
        sspQueue.push(mscq::SaddleSaddlePath(start, end, sepRef, false));
        }
      else if(sepRef.critTypeStart_ == 0 || sepRef.critTypeStart_ == 3) {
        continue;
      }
    }

    if(sepRef.critTypeEnd_ != -1) {
      if(sepRef.critTypeEnd_ == 2 && sepRef.incline_ >= 0) {
        sspQueue.push(mscq::SaddleSaddlePath(end, start, sepRef, true));
      } else if(sepRef.critTypeEnd_ == 0 || sepRef.critTypeEnd_ == 3) {
        continue;
      }
    }

    //this->printMsg("Incline: " + std::to_string(sepRef.incline_));
    {
      if(simpSepMap.find(start) == simpSepMap.end()) {
        simpSepMap.insert({start, std::vector<mscq::Separatrix *>()});
      }

    }

    // Check if id already exists and insert function value decreasing direction
    if(sepRef.incline_ <= 0) {
      if(simpSepMap.find(start) == simpSepMap.end()) {
        simpSepMap.insert({start, std::vector<mscq::Separatrix *>()});
      }
      simpSepMap[start].push_back(&sepRef);
    }
    if(sepRef.incline_ >= 0) {
      if(simpSepMap.find(end) == simpSepMap.end()) {
        simpSepMap.insert({end, std::vector<mscq::Separatrix *>()});
      }
      simpSepMap[end].push_back(&sepRef);
    }
  }

  /*for(auto &simsep : simpSepMap) {
    this->printMsg("SSsize of " + std::to_string(simsep.first.id_) + "("
                   + std::to_string(simsep.first.index_)
                   + ") : " + std::to_string(simsep.second.size()));
  }*/

  // filter separatrices with the same start and end / take the shortest one
  for(auto & simsep : simpSepMap) {
    std::map<mscq::Simplex, std::vector<mscq::Separatrix *>> sepGroups;
    //size_t simsepsize = simsep.second.size();
    for(auto sep : simsep.second) {

      mscq::Simplex start, end;
      sep->getSimplexAt(0, start);
      sep->getSimplexAt(sep->geometry_.size() - 1, end);

      if(simsep.first == start) {
        if(sepGroups.find(end) == sepGroups.end()) {
          sepGroups.insert({end, std::vector<mscq::Separatrix *>()});
        }
        sepGroups[end].push_back(sep);
      } else {
        if(sepGroups.find(start) == sepGroups.end()) {
          sepGroups.insert({start, std::vector<mscq::Separatrix *>()});
        }
        sepGroups[start].push_back(sep);
      }
    }

    simsep.second.clear();

    for(auto &sepgrp : sepGroups) {
      mscq::Separatrix *shortestSep = sepgrp.second.front();
      if(sepgrp.second.size() > 1) {
        for(size_t gsep = 1; gsep < sepgrp.second.size(); ++gsep) {
          if(shortestSep->length() > sepgrp.second[gsep]->length()) {
            shortestSep = sepgrp.second[gsep];
          }
        }
      }
      simsep.second.push_back(shortestSep);
    }

    //this->printMsg("Sepnum: " + std::to_string(simsep.second.size()) + " / "
    //               + std::to_string(simsepsize));
  }

  this->printMsg("Queue size: " + std::to_string(sspQueue.size()));
  size_t maxQueueSize = sspQueue.size();

  std::map<SimplexId, size_t> critPairToSSP;

  while(!sspQueue.empty()) {

    if(sspQueue.size() > maxQueueSize) {
      maxQueueSize = sspQueue.size();
      if(sspQueue.size() > 10000) {
        return 1;
      }
    }

    mscq::SaddleSaddlePath &refSSP = sspQueue.front();

    // path complete?
    if(refSSP.lastId_.index_ == 0
      && vertexToCritId.find(refSSP.lastId_.id_) != vertexToCritId.end()) {
      if(criticalPoints[vertexToCritId[refSSP.lastId_.id_]].index_ == 1) {
        const SimplexId uniqueCritId
          = (vertexToCritId[refSSP.saddleVertId_.id_] * criticalPoints.size())
            + vertexToCritId[refSSP.lastId_.id_];

        if(critPairToSSP.find(uniqueCritId) == critPairToSSP.end()) {
          ssPaths.push_back(refSSP);
          critPairToSSP.insert({uniqueCritId, ssPaths.size() - 1});
        } else {
          if(ssPaths[critPairToSSP[uniqueCritId]].length_ > refSSP.length_) {
            ssPaths[critPairToSSP[uniqueCritId]] = refSSP;
          }
        }
      }
    }

    const size_t numSeps = simpSepMap[refSSP.lastId_].size();

    //this->printMsg("Queuesize: " + std::to_string(sspQueue.size())
    //               + " : numSeps " + std::to_string(numSeps));

    // enqueue extended paths
    for(size_t i = 0; i < numSeps; ++i) {
      mscq::Separatrix *sepPtr = simpSepMap[refSSP.lastId_][i];

      mscq::Simplex start, end, newLastId;
      sepPtr->getSimplexAt(0, start);
      sepPtr->getSimplexAt(sepPtr->geometry_.size() - 1, end);

      bool reversed = false;

      if(refSSP.lastId_ == start) {
        newLastId = end;
      } else {
        reversed = true;
        newLastId = start;
      }

      if(refSSP.visitedIds_.find(simpSepMap[refSSP.lastId_][i]->Id_)
         == refSSP.visitedIds_.end()) {

        mscq::SaddleSaddlePath newPath(refSSP);
        newPath.insertPathPiece(sepPtr, newLastId, reversed);
        sspQueue.push(newPath);
      }
    }

    sspQueue.pop();
  }

  this->printMsg("Max Queue size: " + std::to_string(maxQueueSize));
  this->printMsg("#Paths: " + std::to_string(ssPaths.size()));

  this->printMsg("Computed saddle saddle connectors",
                 0, // progress form 0-1
                 localTimer.getElapsedTime(), // elapsed time so far
                 this->threadNumber_);

  return 1;
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

  this->printMsg("Writing 1-separatrices", 1, localTimer.getElapsedTime(),
                 this->threadNumber_, ttk::debug::LineMode::REPLACE);

  // previous points and cell sizes
  const size_t prevNpoints = *outputSeparatrices1_numberOfPoints_,
               prevNcells = *outputSeparatrices1_numberOfCells_;

  // number of separatrices
  const size_t numSep = separatrices.size();
  size_t npoints = 0;
  size_t numCells = 0;

  // count total number of points and cells
  for(size_t i = 0; i < numSep; ++i) {
    const auto &sep = separatrices[i];

    const int numSepPts = sep.geometry_.size();

    npoints += numSepPts;
    numCells += numSepPts - 1;
  }

  // resize arrays
  outputSeparatrices1_points_->resize((prevNpoints + npoints) * 3);
  auto &points = *outputSeparatrices1_points_;
  outputSeparatrices1_cells_connectivity_->resize((prevNcells + numCells) * 2);
  auto &cellsConn = *outputSeparatrices1_cells_connectivity_;
  if(outputSeparatrices1_cells_sourceIds_ != nullptr)
    outputSeparatrices1_cells_sourceIds_->resize(numSep);
  if(outputSeparatrices1_cells_destinationIds_ != nullptr)
    outputSeparatrices1_cells_destinationIds_->resize(numSep);
  if(outputSeparatrices1_cells_separatrixIds_ != nullptr)
    outputSeparatrices1_cells_separatrixIds_->resize(numSep);
  if(outputSeparatrices1_cells_isOnBoundary_ != nullptr)
    outputSeparatrices1_cells_isOnBoundary_->resize(numSep);

  int (triangulationType::*getPosFuncs[4])(SimplexId, float *) const
    = {&triangulationType::getVertexCoordinates,
       &triangulationType::getEdgeIncenter,
       &triangulationType::getTriangleIncenter,
       &triangulationType::getTetraIncenter};

  SimplexId currPId = prevNpoints, currCId = prevNcells;
  for(size_t i = 0; i < numSep; ++i) {
    const auto &sep = separatrices[i];

    std::array<float, 3> pt{};

    (triangulation.*getPosFuncs[sep.dims_[0]])(sep.geometry_[0], pt.data());

    points[3 * currPId + 0] = pt[0];
    points[3 * currPId + 1] = pt[1];
    points[3 * currPId + 2] = pt[2];

    currPId += 1;

    for(size_t j = 1; j < sep.geometry_.size(); ++j) {
      (triangulation.*getPosFuncs[sep.dims_[j]])(sep.geometry_[j], pt.data());

      points[3 * currPId + 0] = pt[0];
      points[3 * currPId + 1] = pt[1];
      points[3 * currPId + 2] = pt[2];

      cellsConn[2 * currCId + 0] = currPId - 1;
      cellsConn[2 * currCId + 1] = currPId;

      currPId += 1;
      currCId += 1;
    }

    if(outputSeparatrices1_cells_isOnBoundary_ != nullptr)
      (*outputSeparatrices1_cells_isOnBoundary_)[i] = false;
  }

  // update pointers
  *outputSeparatrices1_numberOfPoints_ = prevNpoints + npoints;
  *outputSeparatrices1_numberOfCells_ = prevNcells + numCells;

  this->printMsg("Wrote " + std::to_string(numSep) + " 1-separatrices",
    1, localTimer.getElapsedTime(), this->threadNumber_);

  return 0;
}

template <typename triangulationType>
int ttk::MorseSmaleComplexQuasi::setSaddleSaddleConnectors_3D(
  const std::vector<mscq::SaddleSaddlePath> &ssp,
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

  this->printMsg("Writing saddle saddle connectors", 1, localTimer.getElapsedTime(),
                 this->threadNumber_); //,ttk::debug::LineMode::REPLACE);

  // previous points and cell sizes
  const size_t prevNpoints = *outputSeparatrices1_numberOfPoints_,
               prevNcells = *outputSeparatrices1_numberOfCells_;

  this->printMsg("PRevpts: " + std::to_string(prevNpoints) + " : "
                 + std::to_string(prevNcells));

  // number of separatrices
  const size_t numSSPaths = ssp.size(); // separatrices.size();
  size_t npoints = 0;
  size_t numCells = 0;

  // count total number of points and cells
  for(size_t i = 0; i < numSSPaths; ++i) {
    const size_t numSeps = ssp[i].path_.size();
    for(size_t j = 0; j < numSeps; ++j) {
      mscq::Separatrix *sep = ssp[i].path_[j];

      int numSepPts = sep->geometry_.size();

      npoints += numSepPts;
      numCells += numSepPts - 1;
    }

    npoints -= (numSeps - 1);
  }

  this->printMsg("Npts: " + std::to_string(npoints) + " : "
                 + std::to_string(numCells));

  // resize arrays
  outputSeparatrices1_points_->resize((prevNpoints + npoints) * 3);
  auto &points = *outputSeparatrices1_points_;
  outputSeparatrices1_cells_connectivity_->resize((prevNcells + numCells) * 2);
  auto &cellsConn = *outputSeparatrices1_cells_connectivity_;
  /*if(outputSeparatrices1_cells_sourceIds_ != nullptr)
    outputSeparatrices1_cells_sourceIds_->resize(numSep);
  if(outputSeparatrices1_cells_destinationIds_ != nullptr)
    outputSeparatrices1_cells_destinationIds_->resize(numSep);
  if(outputSeparatrices1_cells_separatrixIds_ != nullptr)
    outputSeparatrices1_cells_separatrixIds_->resize(numSep);
  if(outputSeparatrices1_cells_isOnBoundary_ != nullptr)
    outputSeparatrices1_cells_isOnBoundary_->resize(numSep);*/

  int (triangulationType::*getPosFuncs[4])(SimplexId, float *) const
    = {&triangulationType::getVertexCoordinates,
       &triangulationType::getEdgeIncenter,
       &triangulationType::getTriangleIncenter,
       &triangulationType::getTetraIncenter};

  std::array<float, 3> pt{};

  SimplexId currPId = prevNpoints, currCId = prevNcells;
  for(size_t i = 0; i < numSSPaths; ++i) {

    const size_t numSeps = ssp[i].path_.size();

    (triangulation.*getPosFuncs[ssp[i].path_[0]->dims_.back()])(
      ssp[i].path_[0]->geometry_.back(), pt.data());

    points[3 * currPId + 0] = pt[0];
    points[3 * currPId + 1] = pt[1];
    points[3 * currPId + 2] = pt[2];

    currPId += 1;

    for(size_t j = 0; j < numSeps; ++j) {
      mscq::Separatrix *sep = ssp[i].path_[j];

      if(sep->geometry_.size() == 0) {
        continue;
      }

      const bool reversed = ssp[i].reversed_[j];

      if(reversed) {
        for(size_t k = 1; k < sep->geometry_.size(); ++k) {
          const int index = sep->geometry_.size() - k - 1;
          (triangulation.*getPosFuncs[sep->dims_[index]])(
            sep->geometry_[index], pt.data());

          points[3 * currPId + 0] = pt[0];
          points[3 * currPId + 1] = pt[1];
          points[3 * currPId + 2] = pt[2];

          cellsConn[2 * currCId + 0] = currPId - 1;
          cellsConn[2 * currCId + 1] = currPId;

          currPId += 1;
          currCId += 1;
        }
      } else {
        for(size_t k = 1; k < sep->geometry_.size(); ++k) {
          (triangulation.*getPosFuncs[sep->dims_[k]])(
            sep->geometry_[k], pt.data());

          points[3 * currPId + 0] = pt[0];
          points[3 * currPId + 1] = pt[1];
          points[3 * currPId + 2] = pt[2];

          cellsConn[2 * currCId + 0] = currPId - 1;
          cellsConn[2 * currCId + 1] = currPId;

          currPId += 1;
          currCId += 1;
        }
      }
    }
  }

  this->printMsg("Pid: " + std::to_string(currPId) + " : "
                 + std::to_string(currCId));

  // update pointers
  *outputSeparatrices1_numberOfPoints_ = prevNpoints + npoints;
  *outputSeparatrices1_numberOfCells_ = prevNcells + numCells;

  this->printMsg("Wrote saddle saddle connectors", 1, localTimer.getElapsedTime(),
                 this->threadNumber_);

  return 0;
}

template <typename triangulationType>
int ttk::MorseSmaleComplexQuasi::createBorderPath_3D(
  mscq::Separatrix &separatrix,
  SimplexId startTriangle,
  SimplexId startEdge,
  const SimplexId *const morseSmaleManifold,
  const triangulationType &triangulation) const {

  separatrix.insertGeometry(startTriangle, 2);
  separatrix.insertGeometry(startEdge, 1);

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

    separatrix.insertGeometry(currEdge, 1);

    prevEdge = currEdge;
    prevTriangle = currTriangle;
  }

  separatrix.insertGeometry(currTriangle, 2);

  return 1;
}

template <typename dataType, typename triangulationType>
int ttk::MorseSmaleComplexQuasi::computeIntegralLines(
  std::vector<mscq::Separatrix> &plIntLines,
  const SimplexId *const ascendingNeighbor,
  const SimplexId *const descendingNeighbor,
  const SimplexId *const ascendingManifold,
  const SimplexId *const descendingManifold,
  const std::vector<mscq::Simplex> &criticalPoints,
  const triangulationType &triangulation) const {
    
  ttk::Timer localTimer;

  // print the progress of the current subprocedure (currently 0%)
  this->printMsg("Computing integral lines",
                 0, // progress form 0-1
                 localTimer.getElapsedTime(), // elapsed time so far
                 this->threadNumber_, ttk::debug::LineMode::REPLACE);

  std::vector<std::tuple<SimplexId, SimplexId, int>> critsByIndex[4];
  std::unordered_map<SimplexId, int> critMap;

  for(const auto& crit : criticalPoints) {
    std::map<SimplexId, SimplexId> reachableExtrema;
    const bool isSaddle1 = crit.index_ == 1;
    const bool isSaddle2 = crit.index_ == 2;

    if(isSaddle1 || isSaddle2) {
      SimplexId numNeighbors = triangulation.getVertexNeighborNumber(crit.id_);
      std::set<SimplexId> neigborSet;

      for(SimplexId i = 0; i < numNeighbors; ++i) {
        SimplexId neighbor;
        triangulation.getVertexNeighbor(crit.id_, i, neighbor);

        neigborSet.insert(neighbor);
      }

      for(SimplexId i = 0; i < numNeighbors; ++i) {
        SimplexId neighbor;
        triangulation.getVertexNeighbor(crit.id_, i, neighbor);

        if(isSaddle1) { // negative integration
          if(ascendingNeighbor[neighbor] != crit.id_
             && neigborSet.find(ascendingNeighbor[neighbor]) == neigborSet.end()
             && neigborSet.find(ascendingNeighbor[ascendingNeighbor[neighbor]])
                  == neigborSet.end()) {
            const SimplexId extremum = ascendingManifold[neighbor];
            if(reachableExtrema.find(extremum) == reachableExtrema.end()) {
              reachableExtrema.insert({extremum, neighbor});
            }
          }
        } else { // positive integration
          if(descendingNeighbor[neighbor] != crit.id_
             && neigborSet.find(descendingNeighbor[neighbor])
                == neigborSet.end()
             && neigborSet.find(
                descendingNeighbor[descendingNeighbor[neighbor]])
                == neigborSet.end()) {
            const SimplexId extremum = descendingManifold[neighbor];
            if(reachableExtrema.find(extremum) == reachableExtrema.end()) {
              reachableExtrema.insert({extremum, neighbor});
            }
          }
        }

        const int index = isSaddle1 ? 0 : 1;

        for(const auto& re : reachableExtrema) {
          critsByIndex[index].push_back(
            std::make_tuple(crit.id_, std::get<1>(re), crit.index_));
        }
      }
    }

    critMap.insert({crit.id_, crit.index_});
  }

  for(const auto& c : critsByIndex[0]) { // 1-saddle -> minimum
    mscq::Separatrix plIntLine;

    SimplexId currentVert = std::get<1>(c);
    plIntLine.insertGeometry(std::get<0>(c), 0);

    while(critMap.find(currentVert) == critMap.end()) {
      plIntLine.insertGeometry(currentVert, 0);

      currentVert = ascendingNeighbor[currentVert];
    }

    plIntLine.insertGeometry(currentVert, 0);

    plIntLines.push_back(plIntLine);
  }

  for(const auto& c : critsByIndex[1]) { // 2-saddle maximum
    mscq::Separatrix plIntLine;

    SimplexId currentVert = std::get<1>(c);
    plIntLine.insertGeometry(std::get<0>(c), 0);

    while(critMap.find(currentVert) == critMap.end()) {
      plIntLine.insertGeometry(currentVert, 0);

      currentVert = descendingNeighbor[currentVert];
    }

    plIntLine.insertGeometry(currentVert, 0);

    plIntLines.push_back(plIntLine);
  }

  this->printMsg("Computed integral lines",
                 1, localTimer.getElapsedTime(), this->threadNumber_);

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
  std::vector<mscq::Simplex> &criticalPoints, 
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
        criticalPoints.push_back(mscq::Simplex(candidate, 1));
        numSaddles += 1;
      } else if(lowerComponents == 1 && upperComponents > 1) {
        criticalPoints.push_back(mscq::Simplex(candidate, 2));
        numSaddles += 1;
      }
    } else {
      if((lowerComponents > 1 && upperComponents == 1) ||
        (lowerComponents == 1 && upperComponents > 1) ||
        (lowerComponents > 1 && upperComponents > 1)) {
        criticalPoints.push_back(mscq::Simplex(candidate, 1));
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
                mscq::Simplex(*it, criticalityIndex)
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
  std::vector<mscq::Simplex> &criticalPoints,
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