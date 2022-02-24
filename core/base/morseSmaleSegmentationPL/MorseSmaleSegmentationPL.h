/// \ingroup base
/// \class ttk::MorseSmaleSegmentationPL
/// \author Robin G. C. Maack <maack@rhrk.uni-kl.de>
/// \date June 2021.

#pragma once

#include <Debug.h>
#include <Triangulation.h>
#include <UnionFind.h>
#include <ScalarFieldCriticalPoints.h>

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
const int tetraederLookup[28][15] = {
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

const int tetraederLabelLookup[27][10] = {
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (1) 0,0,0 - 0
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (-) 0,0,1 - 1
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (-) 0,0,2 - 2 
  { 0,  3, -1, -1, -1, -1, -1, -1, -1, -1}, // (2) 0,0,3 - 3
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (-) 0,1,0 - 4
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (-) 0,1,1 - 5 
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (-) 0,1,2 - 6
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (-) 0,1,3 - 7
  { 0,  2, -1, -1, -1, -1, -1, -1, -1, -1}, // (2) 0,2,0 - 8
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (-) 0,2,1 - 9
  { 0,  2,  0,  2, -1, -1, -1, -1, -1, -1}, // (2) 0,2,2 -10
  { 2,  3,  0,  2,  0,  2,  0,  3,  0,  3}, // (3) 0,2,3 -11
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (-) 0,3,0 -12
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (-) 0,3,1 -13
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (-) 0,3,2 -14
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (-) 0,3,3 -15
  { 0,  1, -1, -1, -1, -1, -1, -1, -1, -1}, // (2) 1,0,0 -16
  { 0,  1,  0,  1, -1, -1, -1, -1, -1, -1}, // (2) 1,0,1 -17
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (-) 1,0,2 -18
  { 1,  3,  0,  1,  0,  1,  0,  3,  0,  3}, // (3) 1,0,3 -19
  { 0,  1,  0,  1, -1, -1, -1, -1, -1, -1}, // (2) 1,1,0 -20
  { 0,  1, -1, -1, -1, -1, -1, -1, -1, -1}, // (2) 1,1,1 -21
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (-) 1,1,2 -22
  { 0,  3,  0,  1,  0,  1,  1,  3,  1,  3}, // (3) 1,1,3 -23
  { 1,  2,  0,  1,  0,  1,  0,  2,  0,  2}, // (3) 1,2,0 -24
  { 0,  2,  0,  1,  0,  1,  1,  2,  1,  2}, // (3) 1,2,1 -25
  { 0,  1,  0,  2,  0,  2,  1,  2,  1,  2}  // (3) 1,2,2 -26
};

const size_t tetraederNumTriangles[28] = {
  {0}, // (1) 0,0,0 - 0
  {0}, // (-) 0,0,1 - 1
  {0}, // (-) 0,0,2 - 2 
  {1}, // (2) 0,0,3 - 3
  {0}, // (-) 0,1,0 - 4
  {0}, // (-) 0,1,1 - 5 
  {0}, // (-) 0,1,2 - 6
  {0}, // (-) 0,1,3 - 7
  {1}, // (2) 0,2,0 - 8
  {0}, // (-) 0,2,1 - 9
  {2}, // (2) 0,2,2 -10
  {5}, // (3) 0,2,3 -11
  {0}, // (-) 0,3,0 -12
  {0}, // (-) 0,3,1 -13
  {0}, // (-) 0,3,2 -14
  {0}, // (-) 0,3,3 -15
  {1}, // (2) 1,0,0 -16
  {2}, // (2) 1,0,1 -17
  {0}, // (-) 1,0,2 -18
  {5}, // (3) 1,0,3 -19
  {2}, // (2) 1,1,0 -20
  {1}, // (2) 1,1,1 -21
  {0}, // (-) 1,1,2 -22
  {5}, // (3) 1,1,3 -23
  {5}, // (3) 1,2,0 -24
  {5}, // (3) 1,2,1 -25
  {5}, // (3) 1,2,2 -26
  {12} // (4) 1,2,3 -27
};

const bool tetraederLookupIsMultiLabel[28] = {
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

const bool tetraederLookupIs2Label[28] = {
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
  false, // (3) 0,2,3 -11
  false, // (-) 0,3,0 -12
  false, // (-) 0,3,1 -13
  false, // (-) 0,3,2 -14
  false, // (-) 0,3,3 -15
  true,  // (2) 1,0,0 -16
  true,  // (2) 1,0,1 -17
  false, // (-) 1,0,2 -18
  false, // (3) 1,0,3 -19
  true,  // (2) 1,1,0 -20
  true,  // (2) 1,1,1 -21
  false, // (-) 1,1,2 -22
  false, // (3) 1,1,3 -23
  false, // (3) 1,2,0 -24
  false, // (3) 1,2,1 -25
  false, // (3) 1,2,2 -26
  false  // (4) 1,2,3 -27
};

const bool tetraederLookupIs3Label[28] = {
  false, // (1) 0,0,0 - 0
  false, // (-) 0,0,1 - 1
  false, // (-) 0,0,2 - 2 
  false,  // (2) 0,0,3 - 3
  false, // (-) 0,1,0 - 4
  false, // (-) 0,1,1 - 5 
  false, // (-) 0,1,2 - 6
  false, // (-) 0,1,3 - 7
  false, // (2) 0,2,0 - 8
  false, // (-) 0,2,1 - 9
  false, // (2) 0,2,2 -10
  true,  // (3) 0,2,3 -11
  false, // (-) 0,3,0 -12
  false, // (-) 0,3,1 -13
  false, // (-) 0,3,2 -14
  false, // (-) 0,3,3 -15
  false, // (2) 1,0,0 -16
  false, // (2) 1,0,1 -17
  false, // (-) 1,0,2 -18
  true,  // (3) 1,0,3 -19
  false, // (2) 1,1,0 -20
  false, // (2) 1,1,1 -21
  false, // (-) 1,1,2 -22
  true,  // (3) 1,1,3 -23
  true,  // (3) 1,2,0 -24
  true,  // (3) 1,2,1 -25
  true,  // (3) 1,2,2 -26
  false  // (4) 1,2,3 -27
};


const int tetraederLookupFast[28] = {
  {-1}, // (1) 0,0,0 - 0
  {-1}, // (-) 0,0,1 - 1
  {-1}, // (-) 0,0,2 - 2 
  { 0}, // (2) 0,0,3 - 3
  {-1}, // (-) 0,1,0 - 4
  {-1}, // (-) 0,1,1 - 5 
  {-1}, // (-) 0,1,2 - 6
  {-1}, // (-) 0,1,3 - 7
  { 1}, // (2) 0,2,0 - 8
  {-1}, // (-) 0,2,1 - 9
  {-1}, // (2) 0,2,2 -10
  {-1}, // (3) 0,2,3 -11
  {-1}, // (-) 0,3,0 -12
  {-1}, // (-) 0,3,1 -13
  {-1}, // (-) 0,3,2 -14
  {-1}, // (-) 0,3,3 -15
  { 2}, // (2) 1,0,0 -16
  {-1}, // (2) 1,0,1 -17
  {-1}, // (-) 1,0,2 -18
  {-1}, // (3) 1,0,3 -19
  {-1}, // (2) 1,1,0 -20
  { 3}, // (2) 1,1,1 -21
  {-1}, // (-) 1,1,2 -22
  {-1}, // (3) 1,1,3 -23
  {-1}, // (3) 1,2,0 -24
  {-1}, // (3) 1,2,1 -25
  {-1}, // (3) 1,2,2 -26
  {-1}  // (4) 1,2,3 -27
};

const int tetraederLookupSplitBasisns2Label[22][8] = {
  {-1, -1, -1, -1, -1, -1, -1, -1}, // (1) 0,0,0 - 0
  {-1, -1, -1, -1, -1, -1, -1, -1}, // (-) 0,0,1 - 1
  {-1, -1, -1, -1, -1, -1, -1, -1}, // (-) 0,0,2 - 2 
  { 0,  3,  1,  3,  2,  3, -1, -1}, // (2) 0,0,3 - 3
  {-1, -1, -1, -1, -1, -1, -1, -1}, // (-) 0,1,0 - 4
  {-1, -1, -1, -1, -1, -1, -1, -1}, // (-) 0,1,1 - 5 
  {-1, -1, -1, -1, -1, -1, -1, -1}, // (-) 0,1,2 - 6
  {-1, -1, -1, -1, -1, -1, -1, -1}, // (-) 0,1,3 - 7
  { 0,  2,  1,  2,  3,  2, -1, -1}, // (2) 0,2,0 - 8
  {-1, -1, -1, -1, -1, -1, -1, -1}, // (-) 0,2,1 - 9
  { 0,  2,  1,  3,  0,  3,  1,  2}, // (2) 0,2,2 -10
  {-1, -1, -1, -1, -1, -1, -1, -1}, // (3) 0,2,3 -11
  {-1, -1, -1, -1, -1, -1, -1, -1}, // (-) 0,3,0 -12
  {-1, -1, -1, -1, -1, -1, -1, -1}, // (-) 0,3,1 -13
  {-1, -1, -1, -1, -1, -1, -1, -1}, // (-) 0,3,2 -14
  {-1, -1, -1, -1, -1, -1, -1, -1}, // (-) 0,3,3 -15
  { 0,  1,  2,  1,  3,  1, -1, -1}, // (2) 1,0,0 -16
  { 0,  1,  2,  3,  0,  3,  2,  1}, // (2) 1,0,1 -17
  {-1, -1, -1, -1, -1, -1, -1, -1}, // (-) 1,0,2 -18
  {-1, -1, -1, -1, -1, -1, -1, -1}, // (3) 1,0,3 -19
  { 0,  1,  3,  2,  0,  2,  3,  1}, // (2) 1,1,0 -20
  { 1,  0,  2,  0,  3,  0, -1, -1}, // (2) 1,1,1 -21
};

const int tetraederLookupSplitBasisns3Label[27][4] = {
  {-1, -1, -1, -1}, // (1) 0,0,0 - 0
  {-1, -1, -1, -1}, // (-) 0,0,1 - 1
  {-1, -1, -1, -1}, // (-) 0,0,2 - 2 
  {-1, -1, -1, -1}, // (2) 0,0,3 - 3
  {-1, -1, -1, -1}, // (-) 0,1,0 - 4
  {-1, -1, -1, -1}, // (-) 0,1,1 - 5 
  {-1, -1, -1, -1}, // (-) 0,1,2 - 6
  {-1, -1, -1, -1}, // (-) 0,1,3 - 7
  {-1, -1, -1, -1}, // (2) 0,2,0 - 8
  {-1, -1, -1, -1}, // (-) 0,2,1 - 9
  {-1, -1, -1, -1}, // (2) 0,2,2 -10
  { 0,  2,  1,  3}, // (3) 0,2,3 -11
  {-1, -1, -1, -1}, // (-) 0,3,0 -12
  {-1, -1, -1, -1}, // (-) 0,3,1 -13
  {-1, -1, -1, -1}, // (-) 0,3,2 -14
  {-1, -1, -1, -1}, // (-) 0,3,3 -15
  {-1, -1, -1, -1}, // (2) 1,0,0 -16
  {-1, -1, -1, -1}, // (2) 1,0,1 -17
  {-1, -1, -1, -1}, // (-) 1,0,2 -18
  { 0,  1,  2,  3}, // (3) 1,0,3 -19
  {-1, -1, -1, -1}, // (2) 1,1,0 -20
  {-1, -1, -1, -1}, // (2) 1,1,1 -21
  {-1, -1, -1, -1}, // (-) 1,1,2 -22
  { 1,  0,  2,  3}, // (3) 1,1,3 -23
  { 0,  1,  3,  2}, // (3) 1,2,0 -24
  { 1,  0,  3,  2}, // (3) 1,2,1 -25
  { 2,  0,  3,  1}, // (3) 1,2,2 -26
};

const int lookupOtherLabels[4][3] = {
  {1,2,3},
  {0,2,3},
  {0,1,3},
  {0,1,2},
};

namespace ttk {

  namespace mscq {
    struct Separatrix {
      // default :
      explicit Separatrix()
        : geometry_{} {
      }

      // initialization with one segment :
      explicit Separatrix(const SimplexId segmentGeometry) {
        geometry_.push_back(segmentGeometry);
      }

      explicit Separatrix(const Separatrix &separatrix)
        : geometry_{separatrix.geometry_},
          type_{separatrix.type_} {
      }

      explicit Separatrix(Separatrix &&separatrix) noexcept
        : geometry_{std::move(separatrix.geometry_)},
          type_{separatrix.type_} {
      }

      Separatrix &operator=(Separatrix &&separatrix) noexcept {
        geometry_ = std::move(separatrix.geometry_);
        type_ = separatrix.type_;

        return *this;
      }

      Separatrix &operator=(const Separatrix &separatrix) noexcept {
        geometry_ = std::move(separatrix.geometry_);
        type_ = separatrix.type_;

        return *this;
      }

      void reverse() {
        std::reverse(geometry_.begin(), geometry_.end());
      }

      size_t length() {
        return geometry_.size();
      }

      SimplexId Id_{-1};

      /**
       * Container of ids. Each id addresses a separate
       * container corresponding to a dense representation
       * of the geometry (i.e. separatricesGeometry).
       */
      std::vector<SimplexId> geometry_;

      int type_{-1};
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
  } // namespace mscq

  /**
   * The MorseSmaleSegmentationPL class provides methods to compute a
   * MSC variant that only uses the input field and its order field.
   */
  class MorseSmaleSegmentationPL : virtual public Debug {

  public:
    MorseSmaleSegmentationPL();

    enum class SEPARATRICES_MANIFOLD {
      MORSESMALE = 0,
      ASCENDING = 1,
      DESCENDING = 2
    };
    
    enum class SEPARATRICES1_MODE {
      NONE = 0,
      DETAILED = 1,
      SIMPLE = 2
    };

    enum class SEPARATRICES2_MODE {
      NONE = 0,
      WALLS = 1,
      SEPARATEBASINSFINE = 2,
      SEPARATEBASINSFAST = 3
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
      std::vector<SimplexId> *const outputSeparatrices2_cells_connectivity,
      std::vector<SimplexId> *const outputSeparatrices2_cells_mscIds,
      std::vector<SimplexId> *const separatrices2_cells_cases) {
      outputSeparatrices2_numberOfPoints_ = separatrices2_numberOfPoints;
      outputSeparatrices2_points_ = separatrices2_points;
      outputSeparatrices2_numberOfCells_ = separatrices2_numberOfCells;
      outputSeparatrices2_cells_connectivity_ =
        outputSeparatrices2_cells_connectivity;
      outputSeparatrices2_cells_mscIds_ = outputSeparatrices2_cells_mscIds;
      outputSeparatrices2_cells_cases = separatrices2_cells_cases;
      return 0;
    }

    inline long long getSparseId(
      const long long simID0, const long long simID1,
      const SimplexId &numMSCRegions) const {
      if(simID0 < simID1) {
        return simID0 * numMSCRegions + simID1;
      }
      return simID1 * numMSCRegions + simID0;
    }

    inline long long getSparseId(const long long simID0,
                                 const long long simID1,
                                 const long long simID2,
                                 const SimplexId &numMSCRegions) const {
      if(simID0 < simID1) {
        if(simID1 < simID2) {
          return simID0 * numMSCRegions * numMSCRegions +
          simID1 * numMSCRegions + simID2;
        } else {
          if(simID0 < simID2) {
            return simID0 * numMSCRegions * numMSCRegions +
            simID2 * numMSCRegions + simID1;
          } else {
            return simID2 * numMSCRegions * numMSCRegions +
            simID0 * numMSCRegions + simID1;
          }
        }
      } else {
        if(simID0 < simID2) {
          return simID1 * numMSCRegions * numMSCRegions +
          simID0 * numMSCRegions + simID2;
        } else {
          if(simID1 < simID2) {
            return simID1 * numMSCRegions * numMSCRegions +
            simID2 * numMSCRegions + simID0;
          } else {
            return simID2 * numMSCRegions * numMSCRegions +
            simID1 * numMSCRegions + simID0;
          }
        }
      }
    }

    template <typename triangulationType>
    int execute(
      const SEPARATRICES_MANIFOLD sepManifoldType,
      const SEPARATRICES1_MODE sep1Mode,
      const SEPARATRICES2_MODE sep2Mode, 
      const bool computeSaddles,
      const triangulationType &triangulation);

    template <typename triangulationType>
    int computeManifold(
      std::vector<std::pair<SimplexId, char>> *criticalPoints,
      SimplexId *const manifold,
      SimplexId *const neighbor,
      std::map<SimplexId, SimplexId> &critMap,
      SimplexId &numberOfExtrema,
      const SimplexId *const orderArr,
      const triangulationType &triangulation,
      const bool ascending) const;

    template <typename triangulationType>
    int computeManifolds(
      std::vector<std::pair<SimplexId, char>> *criticalPoints,
      SimplexId *const ascManifold,
      SimplexId *const dscManifold,
      SimplexId *const ascNeighbor,
      SimplexId *const dscNeighbor,
      std::map<SimplexId, SimplexId> &ascCritMap,
      std::map<SimplexId, SimplexId> &dscCritMap,
      int &numberOfMinima,
      int &numberOfMaxima,
      const SimplexId *const orderArr,
      const triangulationType &triangulation) const;

    template <typename triangulationType>
    int computeFinalSegmentation(
      const SimplexId numMaxima,
      const SimplexId *const ascendingManifold,
      const SimplexId *const descendingManifold,
      SimplexId &numManifolds,
      SimplexId *const morseSmaleManifold,
      const triangulationType &triangulation) const;

    template <typename triangulationType>
    int computeSeparatrices1Pieces_2D(
      std::unordered_set<SimplexId> &sep2VertSet,
      std::unordered_map<long long, std::vector<
        std::tuple<SimplexId, SimplexId>>*> &sepPieces,
      std::unordered_map<long long, std::vector<
        std::tuple<SimplexId, SimplexId> >> &triConnectors,
      const SimplexId &numMSCRegions,
      const SimplexId *const morseSmaleManifold,
      const triangulationType &triangulation) const;

    template <typename triangulationType>
    int computeSeparatrices1Inner_2D(
      std::unordered_set<SimplexId> &sep2VertSet,
      std::vector<mscq::Separatrix> &separatrices,
      std::unordered_map<long long,
                         std::vector<std::tuple<SimplexId, SimplexId>> *>
        sepPieces,
      std::unordered_map<long long,
                         std::vector<std::tuple<SimplexId, SimplexId>>>
        triConnectors,
      const SimplexId &numMSCRegions,
      const SimplexId *const morseSmaleManifold,
      const triangulationType &triangulation) const;

    /*template <typename triangulationType>
    int computeSeparatrices1Border_2D(
      std::unordered_set<SimplexId> &sep2VertSet,   
      std::vector<mscq::Separatrix> &separatrices,
      const SimplexId *const morseSmaleManifold,
      const triangulationType &triangulation) const;*/

    template <typename triangulationType>
    int setSeparatrices1Inner_2D(
      const std::vector<mscq::Separatrix> &separatrices,
      const triangulationType &triangulation) const;

    int getSaddles(
      const SimplexId * const offsets,
      std::vector<std::pair<SimplexId, char>> &criticalPoints,
      const AbstractTriangulation &triangulation) const;

    template <typename triangulationType>
    int computeSeparatrices2_3D(
      std::vector<std::array<float, 9>> &trianglePos,
      std::vector<SimplexId> &caseData,
      std::vector<long long> &msLabels,
      const SimplexId &numMSCRegions,
      const SimplexId *const morseSmaleManifold,
      const triangulationType &triangulation) const;

    int computemsLabelMap(
      const std::vector<long long> &msLabels,
      std::map<long long, SimplexId> &msLabelMap) const;

    template <typename triangulationType>
    int computeBasinSeparation_3D_fine(
      std::vector<std::array<float, 9>> &trianglePos,
      std::vector<SimplexId> &caseData,
      std::vector<long long> &msLabels,
      const SimplexId *const morseSmaleManifold,
      const triangulationType &triangulation) const;

    template <typename triangulationType>
    int computeBasinSeparation_3D_fast(
      std::vector<std::array<float, 9>> &trianglePos,
      std::vector<SimplexId> &caseData,
      std::vector<long long> &msLabels,
      const SimplexId *const morseSmaleManifold,
      const triangulationType &triangulation) const;

    template <typename triangulationType>
    int createBorderPath_3D(mscq::Separatrix &separatrix,
                            SimplexId startTriangle,
                            SimplexId startEdge,
                            const SimplexId *const morseSmaleManifold,
                            const triangulationType &triangulation) const;

    template <typename triangulationType>
    int computeSaddleExtremaConnectors_3D(
      std::vector<mscq::Separatrix> &sadExtrConns,
      const SimplexId *const ascendingNeighbor,
      const SimplexId *const descendingNeighbor,
      const SimplexId *const ascendingManifold,
      const SimplexId *const descendingManifold,
      const SimplexId * const orderArr,
      const std::vector<std::pair<SimplexId, char>> &criticalPoints,
      const triangulationType &triangulation) const;

    template <typename triangulationType>
    int computeSaddleExtremaConnectors_3D_fast(
      std::vector<mscq::Separatrix> &sadExtrConns,
      const SimplexId *const ascendingManifold,
      const SimplexId *const descendingManifold,
      const std::map<SimplexId, SimplexId> &ascCritMap,
      const std::map<SimplexId, SimplexId> &dscCritMap,
      const std::vector<std::pair<SimplexId, char>> &criticalPoints,
      const triangulationType &triangulation) const;

    template <typename scalarType, typename triangulationType>
    int computeSaddleSaddleConnectors_3D(
      std::vector<mscq::Separatrix> &sadSadConns,
      const SimplexId *const morseSmaleManifold,
      const SimplexId *const descManifold,
      const std::unordered_set<SimplexId> *sep2VertSet,
      const std::vector<std::pair<SimplexId, char>> &criticalPoints,
      const triangulationType &triangulation) const;

    template <typename triangulationType>
    int setSeparatrices1_3D(const std::vector<mscq::Separatrix> &separatrices,
                            const triangulationType &triangulation) const;

    int setSeparatrices2_3D(
      const std::vector<std::array<float, 9>> &trianglePos,
      const std::vector<SimplexId> &caseData,
      const std::vector<long long> &msLabels,
      std::map<long long, SimplexId> &msLabelMap) const;

    template <typename triangulationType>
    int setCriticalPoints(std::vector<std::pair<SimplexId, char>> &criticalPoints,
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

    inline int interpolatePoints(
      float pos0[3], float pos1[3], float lambda, float result[3]) const {

      result[0] = lambda * pos0[0] + (1-lambda) * pos1[0];
      result[1] = lambda * pos0[1] + (1-lambda) * pos1[1];
      result[2] = lambda * pos0[2] + (1-lambda) * pos1[2];

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
      SimplexId *outputSeparatrices1_numberOfCells_{};
      std::vector<float> *outputSeparatrices1_points_{};
      std::vector<char> *outputSeparatrices1_cells_separatrixTypes_{};


      std::vector<char> *outputSeparatrices1_points_smoothingMask_{};
      std::vector<char> *outputSeparatrices1_points_cellDimensions_{};
      std::vector<SimplexId> *outputSeparatrices1_points_cellIds_{};
      std::vector<SimplexId> *outputSeparatrices1_cells_connectivity_{};
      std::vector<SimplexId> *outputSeparatrices1_cells_sourceIds_{};
      std::vector<SimplexId> *outputSeparatrices1_cells_destinationIds_{};
      std::vector<SimplexId> *outputSeparatrices1_cells_separatrixIds_{};
      std::vector<SimplexId> *outputS1_cells_separatrixFunctionMaximaId_{};
      std::vector<SimplexId> *outputS1_cells_separatrixFunctionMinimaId_{};
      std::vector<char> *outputSeparatrices1_cells_isOnBoundary_{};

      // 2-separatrices
      SimplexId *outputSeparatrices2_numberOfPoints_{};
      std::vector<float> *outputSeparatrices2_points_{};
      SimplexId *outputSeparatrices2_numberOfCells_{};
      std::vector<SimplexId> *outputSeparatrices2_cells_cases{};
      std::vector<SimplexId> *outputSeparatrices2_cells_mscIds_{};
      std::vector<SimplexId> *outputSeparatrices2_cells_connectivity_{};

      // Debug
      bool db_omp = true;

    }; // MorseSmaleSegmentationPL class

} // namespace ttk

template <typename triangulationType>
int ttk::MorseSmaleSegmentationPL::execute(
  const SEPARATRICES_MANIFOLD sepManifoldType,
  const SEPARATRICES1_MODE sep1Mode,
  const SEPARATRICES2_MODE sep2Mode,
  const bool computeSaddles,
  const triangulationType &triangulation)
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

  const SimplexId *orderArr = static_cast<const SimplexId *>(inputOrderField_);
  const SimplexId nV = triangulation.getNumberOfVertices();
  const int dim = triangulation.getDimensionality();

  SimplexId *ascendingManifold
    = static_cast<SimplexId *>(outputAscendingManifold_);
  SimplexId *descendingManifold
    = static_cast<SimplexId *>(outputDescendingManifold_);
  SimplexId *morseSmaleManifold
    = static_cast<SimplexId *>(outputMorseSmaleManifold_);

  SimplexId *ascendingNeighbor  = (SimplexId*) malloc(nV * sizeof(SimplexId));
  SimplexId *descendingNeighbor = (SimplexId*) malloc(nV * sizeof(SimplexId));

  std::vector<std::pair<SimplexId, char>> criticalPoints;
  std::vector<std::pair<SimplexId, char>>* critPointPtr;

  std::map<SimplexId, SimplexId> ascCritMap, dscCritMap;

  if(computeSaddles) {
    critPointPtr = nullptr;
    getSaddles(orderArr, criticalPoints, triangulation);
  } else {
    critPointPtr = &criticalPoints;
  }

  SimplexId numberOfMaxima{0};
  SimplexId numberOfMinima{0};
  SimplexId numberOfMSCRegions{0};
  
  SimplexId *sepManifold;
  bool SepManifoldIsValid = false;

  switch(sepManifoldType) {
    case SEPARATRICES_MANIFOLD::MORSESMALE :
      if(ascendingManifold && descendingManifold && morseSmaleManifold) {
        computeManifolds(
          critPointPtr, ascendingManifold, descendingManifold,
          ascendingNeighbor, descendingNeighbor, ascCritMap, dscCritMap,
          numberOfMinima, numberOfMaxima, orderArr, triangulation);

        computeFinalSegmentation(numberOfMinima, ascendingManifold,
          descendingManifold, numberOfMSCRegions, morseSmaleManifold,
          triangulation);

        sepManifold = static_cast<SimplexId *>(outputMorseSmaleManifold_);
        SepManifoldIsValid = true;
      }
      break;
    case SEPARATRICES_MANIFOLD::ASCENDING :
      if(ascendingManifold) {
        computeManifold<triangulationType>(
          critPointPtr, ascendingManifold, ascendingNeighbor, ascCritMap,
          numberOfMinima, orderArr, triangulation, true);
        
        sepManifold = static_cast<SimplexId *>(outputAscendingManifold_);
        numberOfMSCRegions = numberOfMinima;
        SepManifoldIsValid = true;
      }
      break;
    case SEPARATRICES_MANIFOLD::DESCENDING :
      if(descendingManifold) {
        computeManifold<triangulationType>(
          critPointPtr, descendingManifold, descendingNeighbor, dscCritMap,
          numberOfMaxima, orderArr, triangulation, false);

        sepManifold = static_cast<SimplexId *>(outputDescendingManifold_);
        numberOfMSCRegions = numberOfMaxima;
        SepManifoldIsValid = true;
      }
      break;
    default:
      break;
  }

  if(dim == 2) {
    std::vector<mscq::Separatrix> separatrices1;
    std::unordered_set<SimplexId> sep2VertSet;

    // sparseId, vector of edgeID tuples
    std::unordered_map<long long, std::vector<
      std::tuple<SimplexId, SimplexId>>*> sepPieces;

    // 3 label triangles
    std::unordered_map<long long, std::vector<
      std::tuple<SimplexId, SimplexId> >> triConnectors;
    computeSeparatrices1Pieces_2D<triangulationType>(
      sep2VertSet, sepPieces, triConnectors,
      numberOfMSCRegions, sepManifold, triangulation);

    computeSeparatrices1Inner_2D<triangulationType>(
      sep2VertSet, separatrices1, sepPieces, triConnectors,
      numberOfMSCRegions, sepManifold, triangulation);

    /*computeSeparatrices1Border_2D<triangulationType>(
      sep2VertSet, separatrices1, sepManifold, triangulation
    );*/

    /*findSaddlesFromCadidates<triangulationType>(
      criticalPoints, sep2VertSet, triangulation
    );*/

    setCriticalPoints<triangulationType>(
      criticalPoints, triangulation
    );

    setSeparatrices1Inner_2D<triangulationType>(
      separatrices1, triangulation
    );
  } else if(dim == 3) {
    std::vector<std::array<float, 9>> trianglePos;
    std::vector<long long> msLabels;
    std::map<long long, SimplexId> msLabelMap;
    std::vector<SimplexId> caseData;
    std::vector<mscq::Separatrix> sadExtrConns;
    
    if(SepManifoldIsValid) {
      if(sep2Mode == SEPARATRICES2_MODE::WALLS) {
        computeSeparatrices2_3D<triangulationType>(
          trianglePos, caseData, msLabels, numberOfMSCRegions, sepManifold,
          triangulation);

        computemsLabelMap(msLabels, msLabelMap);

      } else if(sep2Mode == SEPARATRICES2_MODE::SEPARATEBASINSFINE) {
        computeBasinSeparation_3D_fine<triangulationType>(
                                  trianglePos, caseData, msLabels,
                                  sepManifold, triangulation);

        for(int i = -1; i < numberOfMSCRegions; ++i) {
          msLabelMap.insert({i,i});
        }

      } else if(sep2Mode == SEPARATRICES2_MODE::SEPARATEBASINSFAST) {
        computeBasinSeparation_3D_fast<triangulationType>(
                                  trianglePos, caseData, msLabels,
                                  sepManifold, triangulation);

        for(int i = -1; i < numberOfMSCRegions; ++i) {
          msLabelMap.insert({i,i});
        }
      }
    }

    if(ascendingManifold && descendingManifold && computeSaddles) {
      if(sep1Mode == SEPARATRICES1_MODE::SIMPLE) {
        computeSaddleExtremaConnectors_3D_fast<triangulationType>(
        sadExtrConns, ascendingManifold, descendingManifold,
        dscCritMap, ascCritMap, criticalPoints, triangulation);

      } else if(sep1Mode == SEPARATRICES1_MODE::DETAILED) {
        computeSaddleExtremaConnectors_3D<triangulationType>(
        sadExtrConns, ascendingNeighbor, descendingNeighbor,
        ascendingManifold, descendingManifold,
        orderArr, criticalPoints, triangulation);
      }
    }

    setSeparatrices1_3D<triangulationType>(sadExtrConns, triangulation);            

    setSeparatrices2_3D(trianglePos, caseData, msLabels, msLabelMap);
  }

  setCriticalPoints<triangulationType>(criticalPoints, triangulation);

  free(ascendingNeighbor);
  free(descendingNeighbor);

  this->printMsg(ttk::debug::Separator::L1);
  this->printMsg( std::to_string(triangulation.getNumberOfVertices())
                   + " verticies processed",
                 1.0, t.getElapsedTime(), this->threadNumber_);
  this->printMsg(ttk::debug::Separator::L0);

  return 0;
}


template <typename triangulationType>
int ttk::MorseSmaleSegmentationPL::computeManifold(
  std::vector<std::pair<SimplexId, char>> *criticalPoints,
  SimplexId *const manifold,
  SimplexId *const neighbor,
  std::map<SimplexId, SimplexId> &critMap,
  SimplexId &numberOfExtrema,
  const SimplexId * const orderArr,
  const triangulationType &triangulation,
  const bool ascending) const {

  // start global timer
  ttk::Timer localTimer;
  int iterations = 1;
  numberOfExtrema = 0;

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
    std::map<SimplexId, SimplexId> extremas;

    // vertices that may still be compressed
    std::vector<SimplexId>* activeVertices = new std::vector<SimplexId>();
    SimplexId nActiveVertices;

if(!db_omp) { //#ifndef TTK_ENABLE_OPENMP
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
          if(orderArr[neighborId] < orderArr[dmi]) {
            dmi = neighborId;
            hasBiggerNeighbor = true;
          }
        } else {
          if(orderArr[neighborId] > orderArr[dmi]) {
            dmi = neighborId;
            hasBiggerNeighbor = true;
          }
        }
      }
      if(hasBiggerNeighbor) {
        activeVertices->push_back(i);
      }
      else {
        critMap.insert({numberOfExtrema, i});
        extremas.insert({i, numberOfExtrema++});
      }

      neighbor[i] = dmi;
    }

    nActiveVertices = activeVertices->size();
    std::vector<SimplexId>* newActiveVert;

    // compress paths until no changes occur
    while(nActiveVertices > 0) {
      newActiveVert = new std::vector<SimplexId>();
        
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
      manifold[i] = extremas[manifold[i]];
    }
} else { //#else // TTK_ENABLE_OPENMP

    #pragma omp parallel num_threads(this->threadNumber_) \
            shared(iterations, extremas, numberOfExtrema, \
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
        SimplexId &mi = manifold[i];

        mi = i;

        // check all neighbors
        for(SimplexId n = 0; n < numNeighbors; n++) {
          triangulation.getVertexNeighbor(i, n, neighborId);
          if(ascending) {
            if(orderArr[neighborId] < orderArr[mi]) {
              mi = neighborId;
              hasBiggerNeighbor = true;
            }
          } else {
            if(orderArr[neighborId] > orderArr[mi]) {
              mi = neighborId;
              hasBiggerNeighbor = true;
            }
          }
        }

        if(hasBiggerNeighbor) {
            lActiveVertices.push_back(i);
        }
        else {
          #pragma omp critical
          {
            critMap.insert({numberOfExtrema, i});
            extremas.insert({i, numberOfExtrema});
            numberOfExtrema += 1;
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
      #pragma omp for schedule(static)
      for(SimplexId i = 0; i < nVertices; i++) {
        manifold[i] = extremas.at(manifold[i]);
      }
    }
} //#endif // TTK_ENABLE_OPENMP

    const int index =
      ascending ? 0 : triangulation.getDimensionality();

    if(criticalPoints) {
      for (auto& it: extremas) {
        criticalPoints->push_back(std::make_pair(it.first, index));
      }
    }
  
    delete activeVertices;
  }

  {
    this->printMsg("Computed " + manifoldStr + " Manifold",
                   1, localTimer.getElapsedTime(), this->threadNumber_);
  }

  return 1; // return success
}

template <typename triangulationType>
int ttk::MorseSmaleSegmentationPL::computeManifolds(
  std::vector<std::pair<SimplexId, char>> *criticalPoints,
  SimplexId *const ascManifold,
  SimplexId *const dscManifold,
  SimplexId *const ascNeighbor,
  SimplexId *const dscNeighbor,
  std::map<SimplexId, SimplexId> &ascCritMap,
  std::map<SimplexId, SimplexId> &dscCritMap,
  int &numberOfMinima,
  int &numberOfMaxima,
  const SimplexId *const orderArr,
  const triangulationType &triangulation) const {

  ttk::Timer localTimer;
  numberOfMinima = 0;
  numberOfMaxima = 0;

  const bool computeCriticalPoints = criticalPoints != nullptr;
  
  this->printMsg(ttk::debug::Separator::L1);
  // print the progress of the current subprocedure (currently 0%)
  this->printMsg("Computing Manifolds",
                 0, // progress form 0-1
                 0, // elapsed time so far
                 this->threadNumber_, ttk::debug::LineMode::REPLACE);

  /* compute the Descending Maifold iterating over each vertex, searching
   * the biggest neighbor and compressing its path to its maximum */
  const SimplexId nVertices = triangulation.getNumberOfVertices();

  // vertices that may still be compressed
  std::vector<SimplexId>* activeVertices
    = new std::vector<SimplexId>();
  SimplexId nActiveVertices;
  std::map<SimplexId, int> minMap, maxMap;


if(!db_omp) { //#ifndef TTK_ENABLE_OPENMP
  // find maxima and intialize vector of not fully compressed vertices
  for(SimplexId i = 0; i < nVertices; i++) {
    SimplexId neighborId;
    const SimplexId numNeighbors =
      triangulation.getVertexNeighborNumber(i);      

    bool hasBiggerNeighbor = false;
    SimplexId &dmi = dscManifold[i];
    dmi = i;

    bool hasSmallerNeighbor = false;
    SimplexId &ami = ascManifold[i];
    ami = i;

    // check all neighbors
    for(SimplexId n = 0; n < numNeighbors; n++) {
      triangulation.getVertexNeighbor(i, n, neighborId);

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
      if(computeCriticalPoints) {
        criticalPoints->push_back(
          std::make_pair(i, (char)(CriticalType::Local_maximum)));
      }

      dscCritMap.insert({numberOfMaxima, i});
      maxMap.insert({i, numberOfMaxima});
      numberOfMaxima++;
    }

    if(!hasSmallerNeighbor) {
      if(computeCriticalPoints) {
        criticalPoints->push_back(
          std::make_pair(i, (char)(CriticalType::Local_minimum)));
      }

      ascCritMap.insert({numberOfMinima, i});
      minMap.insert({i, numberOfMinima});
      numberOfMinima++;
    }

    ascNeighbor[i] = ami;
    dscNeighbor[i] = dmi;
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

} else { //#else // TTK_ENABLE_OPENMP

  #pragma omp parallel num_threads(this->threadNumber_)
  {
    std::vector<SimplexId> lActiveVertices; //active verticies per thread

    // find the biggest neighbor for each vertex
    #pragma omp for schedule(static)
    for(SimplexId i = 0; i < nVertices; i++) {
      SimplexId neighborId;
      SimplexId numNeighbors =
        triangulation.getVertexNeighborNumber(i);

      bool hasBiggerNeighbor = false;
      SimplexId &dmi = dscManifold[i];
      dmi = i;

      bool hasSmallerNeighbor = false;
      SimplexId &ami = ascManifold[i];
      ami = i;

      // check all neighbors
      for(SimplexId n = 0; n < numNeighbors; n++) {
        triangulation.getVertexNeighbor(i, n, neighborId);

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
            if(computeCriticalPoints) {
              criticalPoints->push_back(
              std::make_pair(i, (char)(CriticalType::Local_maximum)));
            }
            
            dscCritMap.insert({numberOfMaxima, i});
            maxMap.insert({i, numberOfMaxima});
            numberOfMaxima++;
          }

          if(!hasSmallerNeighbor) {
            if(computeCriticalPoints) {
              criticalPoints->push_back(
                std::make_pair(i, (char)(CriticalType::Local_minimum)));
            }

            ascCritMap.insert({numberOfMinima, i});
            minMap.insert({i, numberOfMinima});
            numberOfMinima++;
          }
        }
      }

      ascNeighbor[i] = ami;
      dscNeighbor[i] = dmi;
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

    #pragma omp barrier

    #pragma omp for schedule(static)
    for(SimplexId i = 0; i < nVertices; i++) {
      dscManifold[i] = maxMap[dscManifold[i]];
      ascManifold[i] = minMap[ascManifold[i]];
    }
  }

} //#endif // TTK_ENABLE_OPENMP

  delete activeVertices;
  this->printMsg("Computed Manifolds", 1.0,
                 localTimer.getElapsedTime(), this->threadNumber_);

  return 1; // return success
}

template <typename triangulationType>
int ttk::MorseSmaleSegmentationPL::computeFinalSegmentation(
  const SimplexId numberOfMinima,
  const SimplexId *const ascendingManifold,
  const SimplexId *const descendingManifold,
  SimplexId &numManifolds,
  SimplexId *const morseSmaleManifold,
  const triangulationType &triangulation) const {
  ttk::Timer localTimer;

  this->printMsg("Computing MSC Manifold",
                 0, localTimer.getElapsedTime(), this->threadNumber_,
                 ttk::debug::LineMode::REPLACE);
  
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
  TTK_PSORT(this->threadNumber_, sparseRegionIds.begin(), sparseRegionIds.end());
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

  this->printMsg("Computed MSC Manifold",
                 1, localTimer.getElapsedTime(), this->threadNumber_);

  return 0;
}

template <typename triangulationType>
int ttk::MorseSmaleSegmentationPL::computeSeparatrices1Pieces_2D(
  std::unordered_set<SimplexId> &sep2VertSet,
  std::unordered_map<long long, std::vector<std::tuple<SimplexId, SimplexId>> *>
    &sepPieces,
  std::unordered_map<long long, std::vector<std::tuple<SimplexId, SimplexId>>>
    &triConnectors,
  const SimplexId &numMSCRegions,
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
        const long long sparseID
          = getSparseId(msmEdge[e][0], msmEdge[e][1], numMSCRegions);

        if(triConnectors.find(sparseID) == triConnectors.end()) {
          triConnectors.insert(
            {sparseID, std::vector<std::tuple<SimplexId, SimplexId>>()} );
        }
        triConnectors[sparseID].push_back(std::make_tuple(edges[e], tri));

        SimplexId saddleCandidateVertex[3];
        triangulation.getTriangleVertex(tri, 0, saddleCandidateVertex[0]);
        triangulation.getTriangleVertex(tri, 1, saddleCandidateVertex[1]);
        triangulation.getTriangleVertex(tri, 2, saddleCandidateVertex[2]);

        sep2VertSet.insert(saddleCandidateVertex[0]);
        sep2VertSet.insert(saddleCandidateVertex[1]);
        sep2VertSet.insert(saddleCandidateVertex[2]);
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

      const long long sparseID
        = getSparseId(msmEdge[e0][0], msmEdge[e0][1], numMSCRegions);

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

template <typename triangulationType>
int ttk::MorseSmaleSegmentationPL::computeSeparatrices1Inner_2D(
  std::unordered_set<SimplexId> &sep2VertSet,
  std::vector<mscq::Separatrix> &separatrices,
  std::unordered_map<long long, std::vector<std::tuple<SimplexId, SimplexId>> *>
    sepPieces,
  std::unordered_map<long long, std::vector<std::tuple<SimplexId, SimplexId>>>
    triConnectors,
  const SimplexId &numMSCRegions,
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

      SimplexId startEdgeVerts[2];
      triangulation.getEdgeVertex(startEdge, 0, startEdgeVerts[0]);
      triangulation.getEdgeVertex(startEdge, 1, startEdgeVerts[1]);

      const long long sparseID
        = getSparseId(morseSmaleManifold[startEdgeVerts[0]],
                      morseSmaleManifold[startEdgeVerts[1]], numMSCRegions);

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

        sep2VertSet.insert(vert0);
        sep2VertSet.insert(vert1);
      }
      if(startsOnBoundary) {
        SimplexId &sEdge = newSep.geometry_.front();

        SimplexId vert0, vert1;
        triangulation.getEdgeVertex(sEdge, 0, vert0);
        triangulation.getEdgeVertex(sEdge, 1, vert1);

        sep2VertSet.insert(vert0);
        sep2VertSet.insert(vert1);

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
                 1, // progress form 0-1
                 localTimer.getElapsedTime(), // elapsed time so far
                 this->threadNumber_);

  return 1;
}

/*template <typename triangulationType>
int ttk::MorseSmaleSegmentationPL::computeSeparatrices1Border_2D(
  std::unordered_set<SimplexId> &sep2VertSet,   
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
                 1, // progress form 0-1
                 localTimer.getElapsedTime(), // elapsed time so far
                 this->threadNumber_);
  
  return 1;
}*/

template <typename triangulationType>
int ttk::MorseSmaleSegmentationPL::setSeparatrices1Inner_2D(
  const std::vector<mscq::Separatrix> &separatrices,
  const triangulationType &triangulation) const {

  ttk::Timer localTimer;

  this->printMsg("Writing 1-separatrices",
                 0, 0, this->threadNumber_,
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

    triangulation.getTriangleIncenter(sep.geometry_[0], pt.data());

    points[3 * currPId + 0] = pt[0];
    points[3 * currPId + 1] = pt[1];
    points[3 * currPId + 2] = pt[2];

    currPId += 1;

    for(size_t j = 1; j < sep.geometry_.size() - 1; ++j) {
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
      (*outputSeparatrices1_cells_isOnBoundary_)[i] = false;
  }

  // update pointers
  *outputSeparatrices1_numberOfPoints_ = npoints;
  *outputSeparatrices1_numberOfCells_ = numCells;

  this->printMsg("Wrote 1-separatrices",
                 1, localTimer.getElapsedTime(), this->threadNumber_);

  return 0;
}

int ttk::MorseSmaleSegmentationPL::getSaddles(
  const SimplexId * const orderArr,
  std::vector<std::pair<SimplexId, char>> &criticalPoints,
  const AbstractTriangulation &triangulation) const {

  ttk::Timer localTimer;

  this->printMsg("Computing Critical Points", 0, localTimer.getElapsedTime(),
                 this->threadNumber_, ttk::debug::LineMode::REPLACE);

  const int dim = triangulation.getDimensionality();

  ScalarFieldCriticalPoints sfcp;
  sfcp.setDomainDimension(dim);
  sfcp.setOutput(&criticalPoints);
  sfcp.setThreadNumber(this->threadNumber_);
  sfcp.execute(orderArr, &triangulation);

  this->printMsg("Computed Critical Points", 1, localTimer.getElapsedTime(),
                 this->threadNumber_);

  return 1;
  
}

template <typename triangulationType>
int ttk::MorseSmaleSegmentationPL::computeSeparatrices2_3D(
  std::vector<std::array<float, 9>> &trianglePos,
  std::vector<SimplexId> &caseData,
  std::vector<long long> &msLabels,
  const SimplexId &numMSCRegions,
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

//#ifndef TTK_ENABLE_OPENMP
if(!db_omp) {
  for(SimplexId tet = 0; tet < numTetra; tet++) {
    SimplexId vertices[4];
    triangulation.getCellVertex(tet, 0, vertices[0]);
    triangulation.getCellVertex(tet, 1, vertices[1]);
    triangulation.getCellVertex(tet, 2, vertices[2]);
    triangulation.getCellVertex(tet, 3, vertices[3]);

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
    const int *tetEdgeIndices = tetraederLookup[lookupIndex];
    const int *tetVertLabel = tetraederLabelLookup[lookupIndex];

    if(tetraederLookupIsMultiLabel[lookupIndex]) {
      float edgeCenters[10][3];
      // the 6 edge centers
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

      // the 4 triangle centers
      getCenter(vertPos[0], vertPos[1], vertPos[2], edgeCenters[6]);
      getCenter(vertPos[0], vertPos[1], vertPos[3], edgeCenters[7]);
      getCenter(vertPos[0], vertPos[2], vertPos[3], edgeCenters[8]);
      getCenter(vertPos[1], vertPos[2], vertPos[3], edgeCenters[9]);

      if(tetEdgeIndices[0] == 10) { // 4 labels on tetraeder
        float tetCenter[3];
        triangulation.getCellIncenter(tet, 3, tetCenter);

        // vertex 0
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
          edgeCenters[7][0], edgeCenters[7][1], edgeCenters[7][2],
          edgeCenters[4][0], edgeCenters[4][1], edgeCenters[4][2],
          tetCenter[0], tetCenter[1], tetCenter[2]
        });
        trianglePos.push_back({
          edgeCenters[4][0], edgeCenters[4][1], edgeCenters[4][2],
          edgeCenters[9][0], edgeCenters[9][1], edgeCenters[9][2],
          tetCenter[0], tetCenter[1], tetCenter[2]
        });
        trianglePos.push_back({
          edgeCenters[9][0], edgeCenters[9][1], edgeCenters[9][2], 
          edgeCenters[5][0], edgeCenters[5][1], edgeCenters[5][2], 
          tetCenter[0], tetCenter[1], tetCenter[2]
        });
        trianglePos.push_back({
          edgeCenters[5][0], edgeCenters[5][1], edgeCenters[5][2],
          edgeCenters[8][0], edgeCenters[8][1], edgeCenters[8][2],
          tetCenter[0], tetCenter[1], tetCenter[2]
        });

        long long sparseMSIds[] = {
          getSparseId(msm[0], msm[1], numMSCRegions),
          getSparseId(msm[0], msm[2], numMSCRegions),
          getSparseId(msm[0], msm[3], numMSCRegions),
          getSparseId(msm[1], msm[2], numMSCRegions),
          getSparseId(msm[1], msm[3], numMSCRegions),
          getSparseId(msm[2], msm[3], numMSCRegions)
        };
        msLabels.insert(msLabels.end(), {
          sparseMSIds[0], sparseMSIds[0],
          sparseMSIds[1], sparseMSIds[1],
          sparseMSIds[2], sparseMSIds[2],
          sparseMSIds[3], sparseMSIds[3],
          sparseMSIds[4], sparseMSIds[4],
          sparseMSIds[5], sparseMSIds[5]
        });

        caseData.insert(caseData.end(), {
          lookupIndex, lookupIndex, lookupIndex,
          lookupIndex, lookupIndex, lookupIndex,
          lookupIndex, lookupIndex, lookupIndex,
          lookupIndex, lookupIndex, lookupIndex}
        );
      } else { // 2 or 3 labels on tetraeder
        const size_t numTris = tetraederNumTriangles[lookupIndex];
        long long sparseIds[numTris];
        for(size_t t = 0; t < numTris; ++t) {
          trianglePos.push_back({
            edgeCenters[tetEdgeIndices[(t * 3)    ]][0],
            edgeCenters[tetEdgeIndices[(t * 3)    ]][1],
            edgeCenters[tetEdgeIndices[(t * 3)    ]][2],
            edgeCenters[tetEdgeIndices[(t * 3) + 1]][0],
            edgeCenters[tetEdgeIndices[(t * 3) + 1]][1],
            edgeCenters[tetEdgeIndices[(t * 3) + 1]][2],
            edgeCenters[tetEdgeIndices[(t * 3) + 2]][0],
            edgeCenters[tetEdgeIndices[(t * 3) + 2]][1],
            edgeCenters[tetEdgeIndices[(t * 3) + 2]][2]
          });

          sparseIds[t] = getSparseId(msm[tetVertLabel[t * 2]],
                                     msm[tetVertLabel[(t * 2) + 1]],
                                     numMSCRegions);
          msLabels.push_back(sparseIds[t]);

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
    std::vector<long long> msLabelsLocal;
    std::vector<std::tuple<long long, SimplexId, SimplexId, bool>> sepPcsLocal;
    std::vector<SimplexId> sep2VertSetLocal;
    std::vector<SimplexId> sep2ForkSetLocal;

    #pragma omp for schedule(static) nowait
    for(SimplexId tet = 0; tet < numTetra; tet++) {
      SimplexId vertices[4];
      triangulation.getCellVertex(tet, 0, vertices[0]);
      triangulation.getCellVertex(tet, 1, vertices[1]);
      triangulation.getCellVertex(tet, 2, vertices[2]);
      triangulation.getCellVertex(tet, 3, vertices[3]);

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
      const int *tetEdgeIndices = tetraederLookup[lookupIndex];
      const int *tetVertLabel = tetraederLabelLookup[lookupIndex];

      if(tetraederLookupIsMultiLabel[lookupIndex]) {
        float edgeCenters[10][3];
        // the 6 edge centers
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

        // the 4 triangle centers
        getCenter(vertPos[0], vertPos[1], vertPos[2], edgeCenters[6]);
        getCenter(vertPos[0], vertPos[1], vertPos[3], edgeCenters[7]);
        getCenter(vertPos[0], vertPos[2], vertPos[3], edgeCenters[8]);
        getCenter(vertPos[1], vertPos[2], vertPos[3], edgeCenters[9]);

        if(tetEdgeIndices[0] == 10) { // 4 labels on tetraeder
          float tetCenter[3];
          triangulation.getCellIncenter(tet, 3, tetCenter);

          long long sparseMSIds[6] = {
            getSparseId(msm[0], msm[1], numMSCRegions),
            getSparseId(msm[0], msm[2], numMSCRegions),
            getSparseId(msm[0], msm[3], numMSCRegions),
            getSparseId(msm[1], msm[2], numMSCRegions),
            getSparseId(msm[1], msm[3], numMSCRegions),
            getSparseId(msm[2], msm[3], numMSCRegions)
          };

          msLabelsLocal.insert(msLabelsLocal.end(), {
            sparseMSIds[0], sparseMSIds[0],
            sparseMSIds[1], sparseMSIds[1],
            sparseMSIds[2], sparseMSIds[2],
            sparseMSIds[3], sparseMSIds[3],
            sparseMSIds[4], sparseMSIds[4],
            sparseMSIds[5], sparseMSIds[5]
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
            edgeCenters[7][0], edgeCenters[7][1], edgeCenters[7][2],
            edgeCenters[4][0], edgeCenters[4][1], edgeCenters[4][2],
            tetCenter[0], tetCenter[1], tetCenter[2]
          });
          trianglePosLocal.push_back({
            edgeCenters[4][0], edgeCenters[4][1], edgeCenters[4][2],
            edgeCenters[9][0], edgeCenters[9][1], edgeCenters[9][2],
            tetCenter[0], tetCenter[1], tetCenter[2]
          });
          trianglePosLocal.push_back({
            edgeCenters[9][0], edgeCenters[9][1], edgeCenters[9][2], 
            edgeCenters[5][0], edgeCenters[5][1], edgeCenters[5][2], 
            tetCenter[0], tetCenter[1], tetCenter[2]
          });
          trianglePosLocal.push_back({
            edgeCenters[5][0], edgeCenters[5][1], edgeCenters[5][2],
            edgeCenters[8][0], edgeCenters[8][1], edgeCenters[8][2],
            tetCenter[0], tetCenter[1], tetCenter[2]
          });        

          caseDataLocal.insert(caseDataLocal.end(), {
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
          const size_t numTris = tetraederNumTriangles[lookupIndex];
          long long sparseIds[numTris];
          for(size_t t = 0; t < numTris; ++t) {
            sparseIds[t] = getSparseId(msm[tetVertLabel[t * 2]],
                                       msm[tetVertLabel[(t * 2) + 1]],
                                       numMSCRegions);

            msLabelsLocal.push_back(sparseIds[t]);

            trianglePosLocal.push_back({
              edgeCenters[tetEdgeIndices[(t * 3)    ]][0],
              edgeCenters[tetEdgeIndices[(t * 3)    ]][1],
              edgeCenters[tetEdgeIndices[(t * 3)    ]][2],
              edgeCenters[tetEdgeIndices[(t * 3) + 1]][0],
              edgeCenters[tetEdgeIndices[(t * 3) + 1]][1],
              edgeCenters[tetEdgeIndices[(t * 3) + 1]][2],
              edgeCenters[tetEdgeIndices[(t * 3) + 2]][0],
              edgeCenters[tetEdgeIndices[(t * 3) + 2]][1],
              edgeCenters[tetEdgeIndices[(t * 3) + 2]][2]
            });

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
      msLabels.insert(msLabels.end(),
        msLabelsLocal.begin(), msLabelsLocal.end());
    }
  }
} // #if TTK_OMPENMP
  this->printMsg("Computed 2-Separatrices 3D",
                 1, localTimer.getElapsedTime(), this->threadNumber_);

  return 0;
}

int ttk::MorseSmaleSegmentationPL::computemsLabelMap(
  const std::vector<long long> &msLabels,
  std::map<long long, SimplexId> &msLabelMap) const {

  std::vector<long long> mscLabelCopy(msLabels);

  // get unique "sparse region ids"
  TTK_PSORT(this->threadNumber_, mscLabelCopy.begin(), mscLabelCopy.end());
  const auto last = std::unique(mscLabelCopy.begin(), mscLabelCopy.end());
  mscLabelCopy.erase(last, mscLabelCopy.end());

  for(size_t i = 0; i < mscLabelCopy.size(); ++i) {
    msLabelMap.insert({mscLabelCopy[i], i});
  }

  return 0;
}

template <typename triangulationType>
int ttk::MorseSmaleSegmentationPL::computeBasinSeparation_3D_fine(
  std::vector<std::array<float, 9>> &trianglePos,    
  std::vector<SimplexId> &caseData,
  std::vector<long long> &msLabels,
  const SimplexId *const morseSmaleManifold,
  const triangulationType &triangulation) const {
  ttk::Timer localTimer;

  this->printMsg("Computing 2-Separatrices 3D[F]",
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

  const float d0 = 0.52;
  const float d1 = 1 - d0;

//#ifndef TTK_ENABLE_OPENMP
if(!db_omp) {
  for(SimplexId tet = 0; tet < numTetra; tet++) {
    SimplexId vertices[4];
    triangulation.getCellVertex(tet, 0, vertices[0]);
    triangulation.getCellVertex(tet, 1, vertices[1]);
    triangulation.getCellVertex(tet, 2, vertices[2]);
    triangulation.getCellVertex(tet, 3, vertices[3]);

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

    if(tetraederLookupIsMultiLabel[lookupIndex]) {
      float vPos[4][3];
      triangulation.getVertexPoint(
        vertices[0], vPos[0][0], vPos[0][1], vPos[0][2]);
      triangulation.getVertexPoint(
        vertices[1], vPos[1][0], vPos[1][1], vPos[1][2]);
      triangulation.getVertexPoint(
        vertices[2], vPos[2][0], vPos[2][1], vPos[2][2]);
      triangulation.getVertexPoint(
        vertices[3], vPos[3][0], vPos[3][1], vPos[3][2]);

      if(tetraederLookupIs2Label[lookupIndex]) { // 2 labels (eg. AAAB / AABB)
        const int *vIds = tetraederLookupSplitBasisns2Label[lookupIndex];

        float vert00[3], vert01[3], vert02[3],
              vert10[3], vert11[3], vert12[3];

        interpolatePoints(vPos[vIds[0]], vPos[vIds[1]], d0, vert00);
        interpolatePoints(vPos[vIds[2]], vPos[vIds[3]], d0, vert01);
        interpolatePoints(vPos[vIds[4]], vPos[vIds[5]], d0, vert02);
        interpolatePoints(vPos[vIds[0]], vPos[vIds[1]], d1, vert10);
        interpolatePoints(vPos[vIds[2]], vPos[vIds[3]], d1, vert11);
        interpolatePoints(vPos[vIds[4]], vPos[vIds[5]], d1, vert12);

        trianglePos.push_back({
          vert00[0], vert00[1], vert00[2], 
          vert01[0], vert01[1], vert01[2], 
          vert02[0], vert02[1], vert02[2]
        });
        trianglePos.push_back({
          vert10[0], vert10[1], vert10[2], 
          vert11[0], vert11[1], vert11[2], 
          vert12[0], vert12[1], vert12[2]
        });

        msLabels.push_back(msm[vIds[0]]);
        msLabels.push_back(msm[vIds[1]]);

        caseData.push_back(lookupIndex);
        caseData.push_back(lookupIndex);

        if(vIds[6] != -1) { // 2 vertices per label (e.g. AABB)
          float vert03[3], vert13[3];

          interpolatePoints(vPos[vIds[6]], vPos[vIds[7]], d0, vert03);
          interpolatePoints(vPos[vIds[6]], vPos[vIds[7]], d1, vert13);

          trianglePos.push_back({
            vert00[0], vert00[1], vert00[2], 
            vert01[0], vert01[1], vert01[2], 
            vert03[0], vert03[1], vert03[2]
          });
          trianglePos.push_back({
            vert10[0], vert10[1], vert10[2], 
            vert11[0], vert11[1], vert11[2], 
            vert13[0], vert13[1], vert13[2]
          });

          msLabels.push_back(msm[vIds[0]]);
          msLabels.push_back(msm[vIds[1]]);

          caseData.push_back(lookupIndex);
          caseData.push_back(lookupIndex);
        }
      } else if(tetraederLookupIs3Label[lookupIndex]) {
        const int *vIds = tetraederLookupSplitBasisns3Label[lookupIndex];

        float vert00[3], vert01[3], vert02[3], vert03[3],
              vert10[3], vert11[3], vert12[3],
              vert20[3], vert21[3], vert22[3];

        interpolatePoints(vPos[vIds[0]], vPos[vIds[1]], d0, vert00);
        interpolatePoints(vPos[vIds[2]], vPos[vIds[3]], d0, vert01);
        interpolatePoints(vPos[vIds[0]], vPos[vIds[3]], d0, vert02);
        interpolatePoints(vPos[vIds[2]], vPos[vIds[1]], d0, vert03);

        interpolatePoints(
          vPos[vIds[1]], vPos[lookupOtherLabels[vIds[1]][0]], d0, vert10);
        interpolatePoints(
          vPos[vIds[1]], vPos[lookupOtherLabels[vIds[1]][1]], d0, vert11);
        interpolatePoints(
          vPos[vIds[1]], vPos[lookupOtherLabels[vIds[1]][2]], d0, vert12);

        interpolatePoints(
          vPos[vIds[3]], vPos[lookupOtherLabels[vIds[3]][0]], d0, vert20);
        interpolatePoints(
          vPos[vIds[3]], vPos[lookupOtherLabels[vIds[3]][1]], d0, vert21);
        interpolatePoints(
          vPos[vIds[3]], vPos[lookupOtherLabels[vIds[3]][2]], d0, vert22);

        trianglePos.push_back({
          vert00[0], vert00[1], vert00[2], 
          vert01[0], vert01[1], vert01[2], 
          vert02[0], vert02[1], vert02[2]
        });
        trianglePos.push_back({
          vert00[0], vert00[1], vert00[2], 
          vert01[0], vert01[1], vert01[2], 
          vert03[0], vert03[1], vert03[2]
        });
        trianglePos.push_back({
          vert10[0], vert10[1], vert10[2], 
          vert11[0], vert11[1], vert11[2], 
          vert12[0], vert12[1], vert12[2]
        });
        trianglePos.push_back({
          vert20[0], vert20[1], vert20[2], 
          vert21[0], vert21[1], vert21[2], 
          vert22[0], vert22[1], vert22[2]
        });

        msLabels.push_back(msm[vIds[0]]);
        msLabels.push_back(msm[vIds[0]]);
        msLabels.push_back(msm[vIds[1]]);
        msLabels.push_back(msm[vIds[3]]);

        caseData.push_back(lookupIndex);
        caseData.push_back(lookupIndex);
        caseData.push_back(lookupIndex);
        caseData.push_back(lookupIndex);
      } else { // 4 labels
        float vert00[3], vert01[3], vert02[3],
              vert10[3], vert11[3], vert12[3],
              vert20[3], vert21[3], vert22[3],
              vert30[3], vert31[3], vert32[3];

        interpolatePoints(vPos[0], vPos[1], d0, vert00);
        interpolatePoints(vPos[0], vPos[2], d0, vert01);
        interpolatePoints(vPos[0], vPos[3], d0, vert02);

        interpolatePoints(vPos[1], vPos[0], d0, vert10);
        interpolatePoints(vPos[1], vPos[2], d0, vert11);
        interpolatePoints(vPos[1], vPos[3], d0, vert12);

        interpolatePoints(vPos[2], vPos[0], d0, vert20);
        interpolatePoints(vPos[2], vPos[1], d0, vert21);
        interpolatePoints(vPos[2], vPos[3], d0, vert22);

        interpolatePoints(vPos[3], vPos[0], d0, vert30);
        interpolatePoints(vPos[3], vPos[1], d0, vert31);
        interpolatePoints(vPos[3], vPos[2], d0, vert32);

        trianglePos.push_back({
          vert00[0], vert00[1], vert00[2], 
          vert01[0], vert01[1], vert01[2], 
          vert02[0], vert02[1], vert02[2]
        });
        trianglePos.push_back({
          vert10[0], vert10[1], vert10[2], 
          vert11[0], vert11[1], vert11[2], 
          vert12[0], vert12[1], vert12[2]
        });
        trianglePos.push_back({
          vert20[0], vert20[1], vert20[2], 
          vert21[0], vert21[1], vert21[2], 
          vert22[0], vert22[1], vert22[2]
        });
        trianglePos.push_back({
          vert30[0], vert30[1], vert30[2], 
          vert31[0], vert31[1], vert31[2], 
          vert32[0], vert32[1], vert32[2]
        });

        msLabels.push_back(msm[0]);
        msLabels.push_back(msm[1]);
        msLabels.push_back(msm[2]);
        msLabels.push_back(msm[3]);

        caseData.push_back(lookupIndex);
        caseData.push_back(lookupIndex);
        caseData.push_back(lookupIndex);
        caseData.push_back(lookupIndex);
      }
    }
  }
//#else // TTK_ENABLE_OPENMP
} else {
  #pragma omp parallel num_threads(this->threadNumber_)
  {
    std::vector<std::array<float, 9>> trianglePosLocal;
    std::vector<SimplexId> caseDataLocal;
    std::vector<long long> msLabelsLocal;
    std::vector<SimplexId> sep2VertSetLocal;

    #pragma omp for schedule(static) nowait
    for(SimplexId tet = 0; tet < numTetra; tet++) {
      SimplexId vertices[4];
      triangulation.getCellVertex(tet, 0, vertices[0]);
      triangulation.getCellVertex(tet, 1, vertices[1]);
      triangulation.getCellVertex(tet, 2, vertices[2]);
      triangulation.getCellVertex(tet, 3, vertices[3]);

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

      if(tetraederLookupIsMultiLabel[lookupIndex]) {
        float vPos[4][3];
        triangulation.getVertexPoint(
          vertices[0], vPos[0][0], vPos[0][1], vPos[0][2]);
        triangulation.getVertexPoint(
          vertices[1], vPos[1][0], vPos[1][1], vPos[1][2]);
        triangulation.getVertexPoint(
          vertices[2], vPos[2][0], vPos[2][1], vPos[2][2]);
        triangulation.getVertexPoint(
          vertices[3], vPos[3][0], vPos[3][1], vPos[3][2]);

        if(tetraederLookupIs2Label[lookupIndex]) { // 2 labels (eg. AAAB / AABB)
          const int *vIds = tetraederLookupSplitBasisns2Label[lookupIndex];

          float vert00[3], vert01[3], vert02[3],
                vert10[3], vert11[3], vert12[3];

          interpolatePoints(vPos[vIds[0]], vPos[vIds[1]], d0, vert00);
          interpolatePoints(vPos[vIds[2]], vPos[vIds[3]], d0, vert01);
          interpolatePoints(vPos[vIds[4]], vPos[vIds[5]], d0, vert02);
          interpolatePoints(vPos[vIds[0]], vPos[vIds[1]], d1, vert10);
          interpolatePoints(vPos[vIds[2]], vPos[vIds[3]], d1, vert11);
          interpolatePoints(vPos[vIds[4]], vPos[vIds[5]], d1, vert12);

          trianglePosLocal.push_back({
            vert00[0], vert00[1], vert00[2], 
            vert01[0], vert01[1], vert01[2], 
            vert02[0], vert02[1], vert02[2]
          });
          trianglePosLocal.push_back({
            vert10[0], vert10[1], vert10[2], 
            vert11[0], vert11[1], vert11[2], 
            vert12[0], vert12[1], vert12[2]
          });

          msLabelsLocal.push_back(msm[vIds[0]]);
          msLabelsLocal.push_back(msm[vIds[1]]);

          caseDataLocal.push_back(lookupIndex);
          caseDataLocal.push_back(lookupIndex);

          if(vIds[6] != -1) { // 2 vertices per label (e.g. AABB)
            float vert03[3], vert13[3];

            interpolatePoints(vPos[vIds[6]], vPos[vIds[7]], d0, vert03);
            interpolatePoints(vPos[vIds[6]], vPos[vIds[7]], d1, vert13);

            trianglePosLocal.push_back({
              vert00[0], vert00[1], vert00[2], 
              vert01[0], vert01[1], vert01[2], 
              vert03[0], vert03[1], vert03[2]
            });
            trianglePosLocal.push_back({
              vert10[0], vert10[1], vert10[2], 
              vert11[0], vert11[1], vert11[2], 
              vert13[0], vert13[1], vert13[2]
            });

            msLabelsLocal.push_back(msm[vIds[0]]);
            msLabelsLocal.push_back(msm[vIds[1]]);

            caseDataLocal.push_back(lookupIndex);
            caseDataLocal.push_back(lookupIndex);
          }
        } else if(tetraederLookupIs3Label[lookupIndex]) {
          const int *vIds = tetraederLookupSplitBasisns3Label[lookupIndex];

          float vert00[3], vert01[3], vert02[3], vert03[3],
                vert10[3], vert11[3], vert12[3],
                vert20[3], vert21[3], vert22[3];

          interpolatePoints(vPos[vIds[0]], vPos[vIds[1]], d0, vert00);
          interpolatePoints(vPos[vIds[2]], vPos[vIds[3]], d0, vert01);
          interpolatePoints(vPos[vIds[0]], vPos[vIds[3]], d0, vert02);
          interpolatePoints(vPos[vIds[2]], vPos[vIds[1]], d0, vert03);

          interpolatePoints(
            vPos[vIds[1]], vPos[lookupOtherLabels[vIds[1]][0]], d0, vert10);
          interpolatePoints(
            vPos[vIds[1]], vPos[lookupOtherLabels[vIds[1]][1]], d0, vert11);
          interpolatePoints(
            vPos[vIds[1]], vPos[lookupOtherLabels[vIds[1]][2]], d0, vert12);

          interpolatePoints(
            vPos[vIds[3]], vPos[lookupOtherLabels[vIds[3]][0]], d0, vert20);
          interpolatePoints(
            vPos[vIds[3]], vPos[lookupOtherLabels[vIds[3]][1]], d0, vert21);
          interpolatePoints(
            vPos[vIds[3]], vPos[lookupOtherLabels[vIds[3]][2]], d0, vert22);

          trianglePosLocal.push_back({
            vert00[0], vert00[1], vert00[2], 
            vert01[0], vert01[1], vert01[2], 
            vert02[0], vert02[1], vert02[2]
          });
          trianglePosLocal.push_back({
            vert00[0], vert00[1], vert00[2], 
            vert01[0], vert01[1], vert01[2], 
            vert03[0], vert03[1], vert03[2]
          });
          trianglePosLocal.push_back({
            vert10[0], vert10[1], vert10[2], 
            vert11[0], vert11[1], vert11[2], 
            vert12[0], vert12[1], vert12[2]
          });
          trianglePosLocal.push_back({
            vert20[0], vert20[1], vert20[2], 
            vert21[0], vert21[1], vert21[2], 
            vert22[0], vert22[1], vert22[2]
          });

          msLabelsLocal.push_back(msm[vIds[0]]);
          msLabelsLocal.push_back(msm[vIds[0]]);
          msLabelsLocal.push_back(msm[vIds[1]]);
          msLabelsLocal.push_back(msm[vIds[3]]);

          caseDataLocal.push_back(lookupIndex);
          caseDataLocal.push_back(lookupIndex);
          caseDataLocal.push_back(lookupIndex);
          caseDataLocal.push_back(lookupIndex);
        } else { // 4 labels
          float vert00[3], vert01[3], vert02[3],
                vert10[3], vert11[3], vert12[3],
                vert20[3], vert21[3], vert22[3],
                vert30[3], vert31[3], vert32[3];

          interpolatePoints(vPos[0], vPos[1], d0, vert00);
          interpolatePoints(vPos[0], vPos[2], d0, vert01);
          interpolatePoints(vPos[0], vPos[3], d0, vert02);

          interpolatePoints(vPos[1], vPos[0], d0, vert10);
          interpolatePoints(vPos[1], vPos[2], d0, vert11);
          interpolatePoints(vPos[1], vPos[3], d0, vert12);

          interpolatePoints(vPos[2], vPos[0], d0, vert20);
          interpolatePoints(vPos[2], vPos[1], d0, vert21);
          interpolatePoints(vPos[2], vPos[3], d0, vert22);

          interpolatePoints(vPos[3], vPos[0], d0, vert30);
          interpolatePoints(vPos[3], vPos[1], d0, vert31);
          interpolatePoints(vPos[3], vPos[2], d0, vert32);

          trianglePosLocal.push_back({
            vert00[0], vert00[1], vert00[2], 
            vert01[0], vert01[1], vert01[2], 
            vert02[0], vert02[1], vert02[2]
          });
          trianglePosLocal.push_back({
            vert10[0], vert10[1], vert10[2], 
            vert11[0], vert11[1], vert11[2], 
            vert12[0], vert12[1], vert12[2]
          });
          trianglePosLocal.push_back({
            vert20[0], vert20[1], vert20[2], 
            vert21[0], vert21[1], vert21[2], 
            vert22[0], vert22[1], vert22[2]
          });
          trianglePosLocal.push_back({
            vert30[0], vert30[1], vert30[2], 
            vert31[0], vert31[1], vert31[2], 
            vert32[0], vert32[1], vert32[2]
          });

          msLabelsLocal.push_back(msm[0]);
          msLabelsLocal.push_back(msm[1]);
          msLabelsLocal.push_back(msm[2]);
          msLabelsLocal.push_back(msm[3]);

          caseDataLocal.push_back(lookupIndex);
          caseDataLocal.push_back(lookupIndex);
          caseDataLocal.push_back(lookupIndex);
          caseDataLocal.push_back(lookupIndex);
        }
      }
    }

    #pragma omp critical
    {
      trianglePos.insert(trianglePos.end(),
        trianglePosLocal.begin(), trianglePosLocal.end());
      caseData.insert(caseData.end(),
        caseDataLocal.begin(), caseDataLocal.end());
      msLabels.insert(msLabels.end(),
        msLabelsLocal.begin(), msLabelsLocal.end());
    }
  }
} // TTK_ENABLE_OPENMP
  this->printMsg("Computed 2-Separatrices 3D[F]",
                 1, localTimer.getElapsedTime(), this->threadNumber_);

  return 0;
}

template <typename triangulationType>
int ttk::MorseSmaleSegmentationPL::computeBasinSeparation_3D_fast(
  std::vector<std::array<float, 9>> &trianglePos,    
  std::vector<SimplexId> &caseData,
  std::vector<long long> &mscLabels,
  const SimplexId *const morseSmaleManifold,
  const triangulationType &triangulation) const {
  ttk::Timer localTimer;

  this->printMsg("Computing 2-Separatrices 3D[F]",
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

    if(tetraederLookupIsMultiLabel[lookupIndex]) {
      if(tetraederLookupFast[lookupIndex] != -1) { // <= 2 labels on tetraeder
        float vertPos[4][3];
        triangulation.getVertexPoint(
          vertices[0], vertPos[0][0], vertPos[0][1], vertPos[0][2]);
        triangulation.getVertexPoint(
          vertices[1], vertPos[1][0], vertPos[1][1], vertPos[1][2]);
        triangulation.getVertexPoint(
          vertices[2], vertPos[2][0], vertPos[2][1], vertPos[2][2]);
        triangulation.getVertexPoint(
          vertices[3], vertPos[3][0], vertPos[3][1], vertPos[3][2]);

        switch(tetraederLookupFast[lookupIndex]) {
          case 0:
            trianglePos.push_back({
              vertPos[0][0], vertPos[0][1], vertPos[0][2], 
              vertPos[1][0], vertPos[1][1], vertPos[1][2], 
              vertPos[2][0], vertPos[2][1], vertPos[2][2]
            });
            mscLabels.push_back(msm[0]);
            break;
          case 1:
            trianglePos.push_back({
              vertPos[0][0], vertPos[0][1], vertPos[0][2], 
              vertPos[1][0], vertPos[1][1], vertPos[1][2], 
              vertPos[3][0], vertPos[3][1], vertPos[3][2]
            });
            mscLabels.push_back(msm[0]);
            break;
          case 2:
            trianglePos.push_back({
              vertPos[0][0], vertPos[0][1], vertPos[0][2], 
              vertPos[2][0], vertPos[2][1], vertPos[2][2], 
              vertPos[3][0], vertPos[3][1], vertPos[3][2]
            });
            mscLabels.push_back(msm[0]);
            break;
          case 3:
            trianglePos.push_back({
              vertPos[1][0], vertPos[1][1], vertPos[1][2], 
              vertPos[2][0], vertPos[2][1], vertPos[2][2], 
              vertPos[3][0], vertPos[3][1], vertPos[3][2]
            });
            mscLabels.push_back(msm[1]);
            break;
        }

        caseData.push_back(lookupIndex);
      }
    }
  }
//#else // TTK_ENABLE_OPENMP
} else {
  #pragma omp parallel num_threads(this->threadNumber_)
  {
    std::vector<std::array<float, 9>> trianglePosLocal;
    std::vector<SimplexId> caseDataLocal;
    std::vector<long long> mscLabelsLocal;

    #pragma omp for schedule(static) nowait
    for(SimplexId tet = 0; tet < numTetra; tet++) {
      SimplexId vertices[4];
      triangulation.getCellVertex(tet, 0, vertices[0]);
      triangulation.getCellVertex(tet, 1, vertices[1]);
      triangulation.getCellVertex(tet, 2, vertices[2]);
      triangulation.getCellVertex(tet, 3, vertices[3]);

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

      if(tetraederLookupIsMultiLabel[lookupIndex]) {
        if(tetraederLookupFast[lookupIndex] != -1) { // <= 2 label on tetraeder
          float vertPos[4][3];
          triangulation.getVertexPoint(
            vertices[0], vertPos[0][0], vertPos[0][1], vertPos[0][2]);
          triangulation.getVertexPoint(
            vertices[1], vertPos[1][0], vertPos[1][1], vertPos[1][2]);
          triangulation.getVertexPoint(
            vertices[2], vertPos[2][0], vertPos[2][1], vertPos[2][2]);
          triangulation.getVertexPoint(
            vertices[3], vertPos[3][0], vertPos[3][1], vertPos[3][2]);

          switch(tetraederLookupFast[lookupIndex]) {
            case 0:
              trianglePosLocal.push_back({
                vertPos[0][0], vertPos[0][1], vertPos[0][2], 
                vertPos[1][0], vertPos[1][1], vertPos[1][2], 
                vertPos[2][0], vertPos[2][1], vertPos[2][2]
              });
              mscLabelsLocal.push_back(msm[0]);
              break;
            case 1:
              trianglePosLocal.push_back({
                vertPos[0][0], vertPos[0][1], vertPos[0][2], 
                vertPos[1][0], vertPos[1][1], vertPos[1][2], 
                vertPos[3][0], vertPos[3][1], vertPos[3][2]
              });
              mscLabelsLocal.push_back(msm[0]);
              break;
            case 2:
              trianglePosLocal.push_back({
                vertPos[0][0], vertPos[0][1], vertPos[0][2], 
                vertPos[2][0], vertPos[2][1], vertPos[2][2], 
                vertPos[3][0], vertPos[3][1], vertPos[3][2]
              });
              mscLabelsLocal.push_back(msm[0]);
              break;
            case 3:
              trianglePosLocal.push_back({
                vertPos[1][0], vertPos[1][1], vertPos[1][2], 
                vertPos[2][0], vertPos[2][1], vertPos[2][2], 
                vertPos[3][0], vertPos[3][1], vertPos[3][2]
              });
              mscLabelsLocal.push_back(msm[1]);
              break;
          }

          caseDataLocal.push_back(lookupIndex);
        }
      }
    }

    #pragma omp critical
    {
      trianglePos.insert(trianglePos.end(),
        trianglePosLocal.begin(), trianglePosLocal.end());
      caseData.insert(caseData.end(),
        caseDataLocal.begin(), caseDataLocal.end());
      mscLabels.insert(mscLabels.end(),
        mscLabelsLocal.begin(), mscLabelsLocal.end());
    }
  }
} // TTK_ENABLE_OPENMP
  this->printMsg("Computed 2-Separatrices 3D[F]",
                 1, localTimer.getElapsedTime(), this->threadNumber_);

  return 0;
}

template <typename triangulationType>
int ttk::MorseSmaleSegmentationPL::computeSaddleExtremaConnectors_3D(
  std::vector<mscq::Separatrix> &sadExtrConns,
  const SimplexId *const ascNeighbor,
  const SimplexId *const descNeighbor,
  const SimplexId *const ascManifold,
  const SimplexId *const descManifold,
  const SimplexId *const orderArr,
  const std::vector<std::pair<SimplexId, char>> &criticalPoints,
  const triangulationType &triangulation) const {
    
  ttk::Timer localTimer;
  this->printMsg("Computing crit point connectors",
                 0, localTimer.getElapsedTime(), 
                 this->threadNumber_, ttk::debug::LineMode::REPLACE);

  std::vector<std::pair<SimplexId, SimplexId>> reachableExtrema[2];
  std::queue<mscq::Separatrix> saddleSaddleQueue[2];
  std::unordered_map<SimplexId, int> critMap;
  std::set<SimplexId> saddle1Set, saddle2Set;

  for(size_t c = 0; c < criticalPoints.size(); ++c) {
    const std::pair<SimplexId, char> &crit = criticalPoints[c];
    critMap.insert({crit.first, c});
    
    const bool isSaddle1 = crit.second == (char)CriticalType::Saddle1;
    const bool isSaddle2 = crit.second == (char)CriticalType::Saddle2;
    const int type = isSaddle1 ? 0 : 1;

    if(isSaddle1 || isSaddle2) {
      std::map<SimplexId, SimplexId> neigborSeeds;

      if(isSaddle1) {
        saddle1Set.insert(crit.first);
      } else {
        saddle2Set.insert(crit.first);
      }

      const SimplexId numNeighbors
        = triangulation.getVertexNeighborNumber(crit.first);

      for(SimplexId i = 0; i < numNeighbors; ++i) {
        SimplexId neighbor;
        triangulation.getVertexNeighbor(crit.first, i, neighbor);

        if(isSaddle1) {
          if(orderArr[neighbor] < orderArr[crit.first]) { // s1-min
            const SimplexId& ascManNeigh = ascManifold[neighbor];
            if(neigborSeeds.find(ascManNeigh) == neigborSeeds.end()) {
              neigborSeeds.insert({ascManNeigh, neighbor});
            } else if(
                orderArr[neigborSeeds[ascManNeigh]] > orderArr[neighbor]) {
              neigborSeeds[ascManNeigh] = neighbor;
            }

          }

        } else { // 2-saddle
          if(orderArr[neighbor] > orderArr[crit.first]) { // s2-max
            const SimplexId& dscManNeigh = ascManifold[neighbor];
            if(neigborSeeds.find(dscManNeigh)  == neigborSeeds.end()) {
              neigborSeeds.insert({dscManNeigh, neighbor});
            } else {
              if(orderArr[neigborSeeds[dscManNeigh]] < orderArr[neighbor]) {
                neigborSeeds[descManifold[neighbor]] = neighbor;
              }
            }
          }
        }
      }

      for(const auto& s : neigborSeeds) {
        reachableExtrema[type].push_back(std::make_pair(crit.first, s.second));
      }
    }
  }

  this->printMsg("Computing crit point connectors",
                 0.25, localTimer.getElapsedTime(), 
                 this->threadNumber_, ttk::debug::LineMode::REPLACE);

  std::map<long long, size_t> s1ToMin;
  for(const auto& c : reachableExtrema[0]) { // 1-saddle -> minimum
    mscq::Separatrix sadExtrConn(c.first);
    sadExtrConn.type_ = 0;

    SimplexId currentVert = c.second;

    while(ascNeighbor[currentVert] != currentVert) {
      sadExtrConn.geometry_.push_back(currentVert);
      currentVert = ascNeighbor[currentVert];
    }

    sadExtrConn.geometry_.push_back(currentVert);
    

    long long uniqueId = getSparseId(
      critMap[c.first], critMap[currentVert], criticalPoints.size());

    if(s1ToMin.find(uniqueId) == s1ToMin.end()) {
      s1ToMin.insert({uniqueId, sadExtrConns.size()});
      sadExtrConns.push_back(sadExtrConn);
    } else if(sadExtrConns[s1ToMin[uniqueId]].length() > sadExtrConn.length()) {
      sadExtrConns[s1ToMin[uniqueId]] = sadExtrConn;
    }
  }

  std::map<long long, size_t> s2ToMax;
  for(const auto& c : reachableExtrema[1]) { // 2-saddle maximum
    mscq::Separatrix sadExtrConn(c.first);
    sadExtrConn.type_ = 1;

    SimplexId currentVert = c.second;

    while(descNeighbor[currentVert] != currentVert) {
      sadExtrConn.geometry_.push_back(currentVert);
      currentVert = descNeighbor[currentVert];
    }

    sadExtrConn.geometry_.push_back(currentVert);

    long long uniqueId = getSparseId(
      critMap[c.first], critMap[currentVert], criticalPoints.size());

    if(s2ToMax.find(uniqueId) == s2ToMax.end()) {
      s2ToMax.insert({uniqueId, sadExtrConns.size()});
      sadExtrConns.push_back(sadExtrConn);
    } else if(sadExtrConns[s2ToMax[uniqueId]].length() > sadExtrConn.length()) {
      sadExtrConns[s2ToMax[uniqueId]] = sadExtrConn;
    }
  }

  this->printMsg("Computed crit point connectors",
                 1, localTimer.getElapsedTime(), this->threadNumber_);

  return 1;
}

template <typename triangulationType>
int ttk::MorseSmaleSegmentationPL::computeSaddleExtremaConnectors_3D_fast(
  std::vector<mscq::Separatrix> &sadExtrConns,
  const SimplexId *const ascendingManifold,
  const SimplexId *const descendingManifold,
  const std::map<SimplexId, SimplexId> &ascCritMap,
  const std::map<SimplexId, SimplexId> &dscCritMap,
  const std::vector<std::pair<SimplexId, char>> &criticalPoints,
  const triangulationType &triangulation) const {
    
  ttk::Timer localTimer;

  // print the progress of the current subprocedure (currently 0%)
  this->printMsg("Computing crit point connectors[F]",
                 0, // progress form 0-1
                 localTimer.getElapsedTime(), // elapsed time so far
                 this->threadNumber_, ttk::debug::LineMode::REPLACE);

  for(const auto& crit : criticalPoints) {
    const bool isSaddle1 = crit.second == (char)CriticalType::Saddle1;
    const bool isSaddle2 = crit.second == (char)CriticalType::Saddle2;

    if(isSaddle1 || isSaddle2) {
      SimplexId numNeighbors
        = triangulation.getVertexNeighborNumber(crit.first);
      std::set<SimplexId> neigborSet;

      for(SimplexId i = 0; i < numNeighbors; ++i) {
        SimplexId neighbor;
        triangulation.getVertexNeighbor(crit.first, i, neighbor);

        if(isSaddle1) {
          neighbor = ascCritMap.at(ascendingManifold[neighbor]);
        } else {
          neighbor = dscCritMap.at(descendingManifold[neighbor]);
        }

        if(neigborSet.find(neighbor) == neigborSet.end()) {
          mscq::Separatrix sadExtrConn(crit.first);
          if(isSaddle1) {
            sadExtrConn.type_ = 0;
          } else {
            sadExtrConn.type_ = 2;
          }

          sadExtrConn.geometry_.push_back(neighbor);
          sadExtrConns.push_back(sadExtrConn);

          neigborSet.insert(neighbor);
        }
      }
    }
  }

  this->printMsg("Computed crit point connectors[F]",
                 1, // progress form 0-1
                 localTimer.getElapsedTime(), // elapsed time so far
                 this->threadNumber_);

  return 1;
}

template <typename triangulationType>
int ttk::MorseSmaleSegmentationPL::setCriticalPoints(
  std::vector<std::pair<SimplexId, char>> &criticalPoints,
  const triangulationType &triangulation) const {

  ttk::Timer localTimer;

  this->printMsg("Writing critical points",
                 0, localTimer.getElapsedTime(), this->threadNumber_,
                 ttk::debug::LineMode::REPLACE);
  
  float x, y, z;

  for(const auto crit : criticalPoints) {
    triangulation.getVertexPoint(crit.first, x, y, z);
    outputCriticalPoints_points_->push_back({x, y, z});
    outputCriticalPoints_points_cellDimensions_->push_back(crit.second);
    outputCriticalPoints_points_cellIds_->push_back(crit.first);
  }

  this->printMsg("Wrote " + std::to_string(criticalPoints.size()) +
                 " critical points",
                 1, localTimer.getElapsedTime(), this->threadNumber_);

  return 1;
}

template <typename triangulationType>
int ttk::MorseSmaleSegmentationPL::setSeparatrices1_3D(
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
  if(outputSeparatrices1_cells_separatrixTypes_ != nullptr)
    outputSeparatrices1_cells_separatrixTypes_->resize((prevNcells + numCells));

  SimplexId currPId = prevNpoints, currCId = prevNcells;
  for(size_t i = 0; i < numSep; ++i) {
    const auto &sep = separatrices[i];

    std::array<float, 3> pt{};
    triangulation.getVertexPoint(sep.geometry_[0], pt[0], pt[1], pt[2]);

    points[3 * currPId + 0] = pt[0];
    points[3 * currPId + 1] = pt[1];
    points[3 * currPId + 2] = pt[2];

    currPId += 1;

    for(size_t j = 1; j < sep.geometry_.size(); ++j) {
      triangulation.getVertexPoint(sep.geometry_[j], pt[0], pt[1], pt[2]);

      points[3 * currPId + 0] = pt[0];
      points[3 * currPId + 1] = pt[1];
      points[3 * currPId + 2] = pt[2];

      cellsConn[2 * currCId + 0] = currPId - 1;
      cellsConn[2 * currCId + 1] = currPId;

      if(outputSeparatrices1_cells_separatrixTypes_ != nullptr)
        (*outputSeparatrices1_cells_separatrixTypes_)[currCId] = sep.type_;

      currPId += 1;
      currCId += 1;
    }
  }

  // update pointers
  *outputSeparatrices1_numberOfPoints_ = prevNpoints + npoints;
  *outputSeparatrices1_numberOfCells_ = prevNcells + numCells;

  this->printMsg("Wrote " + std::to_string(numSep) + " 1-separatrices",
    1, localTimer.getElapsedTime(), this->threadNumber_);

  return 0;
}

int ttk::MorseSmaleSegmentationPL::setSeparatrices2_3D(
  const std::vector<std::array<float, 9>> &trianglePos,
  const std::vector<SimplexId> &caseData,
  const std::vector<long long> &msLabels,
  std::map<long long, SimplexId> &msLabelMap) const {

  ttk::Timer localTimer;

  this->printMsg("Writing 2-separatrices",
                  0,
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
  auto &cellsMSCIds = *outputSeparatrices2_cells_mscIds_;
  outputSeparatrices2_cells_mscIds_->resize(numTriangles);
  auto &cellsCase = *outputSeparatrices2_cells_cases;
  cellsCase = caseData;

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for schedule(static) if(numTriangles > 100000)
#endif
  for(int tri = 0; tri < numTriangles; ++tri) {
    auto &triPos = trianglePos[tri];
    points[9 * tri    ] = triPos[0];
    points[9 * tri + 1] = triPos[1];
    points[9 * tri + 2] = triPos[2];

    points[9 * tri + 3] = triPos[3];
    points[9 * tri + 4] = triPos[4];
    points[9 * tri + 5] = triPos[5];

    points[9 * tri + 6] = triPos[6];
    points[9 * tri + 7] = triPos[7];
    points[9 * tri + 8] = triPos[8];

    cellsConn[3 * tri    ] = 3 * tri;
    cellsConn[3 * tri + 1] = 3 * tri + 1;
    cellsConn[3 * tri + 2] = 3 * tri + 2;

    cellsMSCIds[tri] = msLabelMap[msLabels[tri]];
  }

  // update pointers
  *outputSeparatrices2_numberOfPoints_ = npoints;
  *outputSeparatrices2_numberOfCells_ = numCells;

  this->printMsg("Wrote " + std::to_string(numTriangles) + " 2-separatrices",
                 1, localTimer.getElapsedTime(), this->threadNumber_);

  return 1;
}
