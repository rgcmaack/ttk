/// \ingroup base
/// \class ttk::MorseSmaleSegmentationPL
/// \author Robin G. C. Maack <maack@rhrk.uni-kl.de>
/// \date June 2021.
///
/// \brief TTK processing package for the computation of basin segmentations.
/// 
/// \details This TTK processing package allows to create a segmentation of 
/// basins (Areas flowing towards Minima, Maxima). It first finds the maximum
/// and minimum for each vertex, then surrounds areas with the same minimum and
/// maximum with a wall. Saddles can be computed and their integral lines, 
/// towards their corresponding extrema are generated.
/// 
///
/// \sa ttk::Triangulation
/// \sa ttk::ScalarFieldCriticalPoints

#pragma once

#include <Debug.h>
#include <Triangulation.h>
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

const bool tetraederLookupIs4Label[28] = {
  false, // (1) 0,0,0 - 0
  false, // (-) 0,0,1 - 1
  false, // (-) 0,0,2 - 2 
  false, // (2) 0,0,3 - 3
  false, // (-) 0,1,0 - 4
  false, // (-) 0,1,1 - 5 
  false, // (-) 0,1,2 - 6
  false, // (-) 0,1,3 - 7
  false, // (2) 0,2,0 - 8
  false, // (-) 0,2,1 - 9
  false, // (2) 0,2,2 -10
  false, // (3) 0,2,3 -11
  false, // (-) 0,3,0 -12
  false, // (-) 0,3,1 -13
  false, // (-) 0,3,2 -14
  false, // (-) 0,3,3 -15
  false, // (2) 1,0,0 -16
  false, // (2) 1,0,1 -17
  false, // (-) 1,0,2 -18
  false, // (3) 1,0,3 -19
  false, // (2) 1,1,0 -20
  false, // (2) 1,1,1 -21
  false, // (-) 1,1,2 -22
  false, // (3) 1,1,3 -23
  false, // (3) 1,2,0 -24
  false, // (3) 1,2,1 -25
  false, // (3) 1,2,2 -26
  true   // (4) 1,2,3 -27
};


const int tetraederLookupFastCase[28] = {
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

const bool tetraederLookupFast[28] = {
  {false}, // (1) 0,0,0 - 0
  {false}, // (-) 0,0,1 - 1
  {false}, // (-) 0,0,2 - 2 
  {true}, // (2) 0,0,3 - 3
  {false}, // (-) 0,1,0 - 4
  {false}, // (-) 0,1,1 - 5 
  {false}, // (-) 0,1,2 - 6
  {false}, // (-) 0,1,3 - 7
  {true}, // (2) 0,2,0 - 8
  {false}, // (-) 0,2,1 - 9
  {false}, // (2) 0,2,2 -10
  {false}, // (3) 0,2,3 -11
  {false}, // (-) 0,3,0 -12
  {false}, // (-) 0,3,1 -13
  {false}, // (-) 0,3,2 -14
  {false}, // (-) 0,3,3 -15
  {true}, // (2) 1,0,0 -16
  {false}, // (2) 1,0,1 -17
  {false}, // (-) 1,0,2 -18
  {false}, // (3) 1,0,3 -19
  {false}, // (2) 1,1,0 -20
  {true}, // (2) 1,1,1 -21
  {false}, // (-) 1,1,2 -22
  {false}, // (3) 1,1,3 -23
  {false}, // (3) 1,2,0 -24
  {false}, // (3) 1,2,1 -25
  {false}, // (3) 1,2,2 -26
  {false}  // (4) 1,2,3 -27
};

const int tetraederLookupFastTri[4][3] = {
  {0, 1, 2},
  {0, 1, 3},
  {0, 2, 3},
  {1, 2, 3}
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

const int tetraederLookupSplitBasisns3Label[27][11] = {
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (1) 0,0,0 - 0
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (-) 0,0,1 - 1
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (-) 0,0,2 - 2 
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (2) 0,0,3 - 3
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (-) 0,1,0 - 4
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (-) 0,1,1 - 5 
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (-) 0,1,2 - 6
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (-) 0,1,3 - 7
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (2) 0,2,0 - 8
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (-) 0,2,1 - 9
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (2) 0,2,2 -10
  { 0,  1,  2,  2,  1,  3,  4,  3,  2,  3,  5}, // (3) 0,2,3 -11
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (-) 0,3,0 -12
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (-) 0,3,1 -13
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (-) 0,3,2 -14
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (-) 0,3,3 -15
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (2) 1,0,0 -16
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (2) 1,0,1 -17
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (-) 1,0,2 -18
  { 0,  2,  0,  1,  2,  5,  3,  3,  3,  1,  4}, // (3) 1,0,3 -19
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (2) 1,1,0 -20
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (2) 1,1,1 -21
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // (-) 1,1,2 -22
  { 1,  0,  4,  1,  2,  1,  5,  2,  0,  3,  2}, // (3) 1,1,3 -23
  { 0,  0,  1,  0,  3,  4,  5,  3,  1,  2,  3}, // (3) 1,2,0 -24
  { 1,  0,  3,  0,  3,  2,  5,  2,  0,  2,  1}, // (3) 1,2,1 -25
  { 2,  1,  3,  0,  3,  2,  4,  1,  0,  1,  0}, // (3) 1,2,2 -26
};

const int lookupOtherLabels[4][3] = {
  {1,2,3},
  {0,2,3},
  {0,1,3},
  {0,1,2},
};

const bool triangleLookupIsMultiLabel[7] = {
  false, // (1) 0,0 - 0
  false, // (-) 0,1 - 1
  true,  // (2) 0,2 - 2 
  false, // (-) 0,3 - 3
  true,  // (2) 1,0 - 4
  true,  // (2) 1,1 - 5 
  true   // (3) 1,2 - 6
};

const bool triangleLookupIs2Label[7] = {
  false, // (1) 0,0 - 0
  false, // (-) 0,1 - 1
  true,  // (2) 0,2 - 2 
  false, // (-) 0,3 - 3
  true,  // (2) 1,0 - 4
  true,  // (2) 1,1 - 5 
  false  // (3) 1,2 - 6
};

const int triangleLookupEdgeVerts[7][4] = {
  {-1, -1, -1, -1}, // (1) 0,0 - 0
  {-1, -1, -1, -1}, // (-) 0,1 - 1
  { 1,  2,  0,  2}, // (2) 0,2 - 2 
  {-1, -1, -1, -1}, // (-) 0,3 - 3
  { 0,  1,  1,  2}, // (2) 1,0 - 4
  { 0,  1,  0,  2}, // (2) 1,1 - 5
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
      SEPARATEBASINSFAST = 3,
      EXPERIMENT = 4,
    };    

    int preconditionTriangulation(
      ttk::AbstractTriangulation *triangulation) const {
      int success = 0;
      success += triangulation->preconditionVertexNeighbors();
      success += triangulation->preconditionBoundaryVertices();
      success += triangulation->preconditionVertexStars();

      if(triangulation->getDimensionality() ==  2) {
        
      }

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
      SimplexId *const separatrices1_numberOfCells,
      std::vector<SimplexId> *const separatrices1_cells_connectivity,
      std::vector<char> *const separatrices1_cells_separatrixTypes,
      std::vector<float> *const separatrices1_cells_extremaDistance,
      std::vector<SimplexId> *const separatrices1_cells_extremaDistanceAbs) {
      outputSeparatrices1_numberOfPoints_ = separatrices1_numberOfPoints;
      outputSeparatrices1_points_ = separatrices1_points;
      outputSeparatrices1_numberOfCells_ = separatrices1_numberOfCells;
      outputSeparatrices1_cells_connectivity_
        = separatrices1_cells_connectivity;
      outputSeparatrices1_cells_separatrixTypes_
        = separatrices1_cells_separatrixTypes;
      outputseparatrices1_cells_extremaDistance_ =
        separatrices1_cells_extremaDistance;
      outputseparatrices1_cells_extremaDistanceAbs_ =
        separatrices1_cells_extremaDistanceAbs;
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
      outputSeparatrices2_cells_cases_ = separatrices2_cells_cases;
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
      SimplexId &numberOfMinima,
      SimplexId &numberOfMaxima,
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
      std::vector<std::array<float, 6>> &edgePos,
      std::vector<long long> &msLabels,
      std::unordered_set<SimplexId> &sep2VertSet,
      const SimplexId &numMSCRegions,
      const SimplexId *const morseSmaleManifold,
      const triangulationType &triangulation) const;

    template <typename triangulationType>
    int findSaddlesFromCadidates(
      std::vector<std::pair<SimplexId, char>> &criticalPoints,
      const std::unordered_set<SimplexId> &sep2VertSet,
      const triangulationType &triangulation) const;

    int setSeparatrices2_2D(
      const std::vector<std::array<float, 6>> &edgePos,
      const std::vector<long long> &msLabels,
      std::map<long long, SimplexId> &msLabelMap) const;

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

    template <typename triangulationType>
    int computeSeparatrices2_3D_split(
      std::vector<std::array<float, 9>> &trianglePos,
      std::vector<SimplexId> &caseData,
      std::vector<long long> &msLabels,
      const SimplexId &numAscRegions,
      const SimplexId &numDscRegions,
      const SimplexId *const ascManifold,
      const SimplexId *const dscManifold,
      const triangulationType &triangulation) const;

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

    int computeMSLabelMap(
      const std::vector<long long> &msLabels,
      std::map<long long, SimplexId> &msLabelMap) const;

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

      incenter[0] = 0.5 * (p[0] + p[3]);
      incenter[1] = 0.5 * (p[1] + p[4]);
      incenter[2] = 0.5 * (p[2] + p[5]);

      return 0;
    }

    inline int getCenter(
      float pos0[3], float pos1[3], float incenter[3]) const {
      incenter[0] = 0.5 * (pos0[0] + pos1[0]);
      incenter[1] = 0.5 * (pos0[1] + pos1[1]);
      incenter[2] = 0.5 * (pos0[2] + pos1[2]);

      return 0;
    }

    inline int getCenter(
      float pos0[3], float pos1[3], float pos2[3], float incenter[3]) const {
      incenter[0] = 0.3333 * (pos0[0] + pos1[0] + pos2[0]);
      incenter[1] = 0.3333 * (pos0[1] + pos1[1] + pos2[1]);
      incenter[2] = 0.3333 * (pos0[2] + pos1[2] + pos2[2]);

      return 0;
    }

    inline int getCenter(
      float pos0[3], float pos1[3], float pos2[3], float pos3[3],
      float incenter[3]) const {
      incenter[0] = 0.25 * (pos0[0] + pos1[0] + pos2[0] + pos3[0]);
      incenter[1] = 0.25 * (pos0[1] + pos1[1] + pos2[1] + pos3[1]);
      incenter[2] = 0.25 * (pos0[2] + pos1[2] + pos2[2] + pos3[2]);

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
      std::vector<SimplexId> *outputSeparatrices1_cells_connectivity_{};
      std::vector<float> *outputseparatrices1_cells_extremaDistance_{};
      std::vector<SimplexId> *outputseparatrices1_cells_extremaDistanceAbs_{};

      // 2-separatrices
      SimplexId *outputSeparatrices2_numberOfPoints_{};
      std::vector<float> *outputSeparatrices2_points_{};
      SimplexId *outputSeparatrices2_numberOfCells_{};
      std::vector<SimplexId> *outputSeparatrices2_cells_cases_{};
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

  SimplexId *ascendingNeighbor = nullptr;
  SimplexId *descendingNeighbor = nullptr;

  if(computeSaddles && sep1Mode == SEPARATRICES1_MODE::DETAILED && dim == 3) {
    ascendingNeighbor = (SimplexId*) malloc(nV * sizeof(SimplexId));
    descendingNeighbor = (SimplexId*) malloc(nV * sizeof(SimplexId));
  }

  std::vector<std::pair<SimplexId, char>> criticalPoints;
  std::vector<std::pair<SimplexId, char>>* critPointPtr;

  std::map<SimplexId, SimplexId> ascCritMap, dscCritMap;

  if(computeSaddles && dim == 3) {
    critPointPtr = nullptr;
    getSaddles(orderArr, criticalPoints, triangulation);
  } else {
    critPointPtr = &criticalPoints;
  }

  SimplexId numberOfMaxima{0};
  SimplexId numberOfMinima{0};
  SimplexId numberOfMSCRegions{0};
  
  SimplexId *sepManifold = nullptr;
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

  if(dim == 3) {
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

        computeMSLabelMap(msLabels, msLabelMap);

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

      } else if(sep2Mode == SEPARATRICES2_MODE::EXPERIMENT) {
        computeSeparatrices2_3D_split<triangulationType>(
          trianglePos, caseData, msLabels, numberOfMinima, numberOfMinima,
          ascendingManifold, descendingManifold, triangulation);

        computeMSLabelMap(msLabels, msLabelMap);
      }
    }

    if(ascendingManifold && descendingManifold && computeSaddles) {
      if(sep1Mode == SEPARATRICES1_MODE::SIMPLE) {
        computeSaddleExtremaConnectors_3D_fast<triangulationType>(
        sadExtrConns, ascendingManifold, descendingManifold,
        ascCritMap, dscCritMap, criticalPoints, triangulation);

      } else if(sep1Mode == SEPARATRICES1_MODE::DETAILED) {
        computeSaddleExtremaConnectors_3D<triangulationType>(
        sadExtrConns, ascendingNeighbor, descendingNeighbor,
        ascendingManifold, descendingManifold,
        orderArr, criticalPoints, triangulation);
      }

      if(computeSaddles && sep1Mode == SEPARATRICES1_MODE::DETAILED) {
        free(ascendingNeighbor);
        free(descendingNeighbor);
      }
    }

    setSeparatrices1_3D<triangulationType>(sadExtrConns, triangulation);            

    setSeparatrices2_3D(trianglePos, caseData, msLabels, msLabelMap);
  } else if(dim == 2) {
    std::vector<std::array<float, 6>> edgePos;
    std::vector<long long> msLabels;
    std::unordered_set<SimplexId> sep2VertSet;
    std::map<long long, SimplexId> msLabelMap;

    computeSeparatrices1Pieces_2D(
      edgePos, msLabels, sep2VertSet, numberOfMSCRegions, 
      morseSmaleManifold, triangulation);

    if(computeSaddles) {
      findSaddlesFromCadidates(criticalPoints, sep2VertSet, triangulation);
    }

    computeMSLabelMap(msLabels, msLabelMap);

    setSeparatrices2_2D(edgePos, msLabels, msLabelMap);
  }

  setCriticalPoints<triangulationType>(criticalPoints, triangulation);

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
  this->printMsg(ttk::debug::Separator::L1);

  const std::string manifoldStr = ascending ? "Ascending" : "Descending";
  const std::string minMaxStr = ascending ? "Minima" : "Maxima";
  // print the progress of the current subprocedure (currently 0%)
  this->printMsg("Computing " + manifoldStr + " Manifold",
                  0, // progress form 0-1
                  0, // elapsed time so far
                  this->threadNumber_, ttk::debug::LineMode::REPLACE);

  int iterations = 1;
  numberOfExtrema = 0;

  
  const bool computeNeighbor = (neighbor != nullptr);

  // -----------------------------------------------------------------------
  // Compute MSC variant
  // -----------------------------------------------------------------------
  {
    /* compute the Descending Maifold iterating over each vertex, searching
     * the biggest neighbor and compressing its path to its maximum */
    const SimplexId nVertices = triangulation.getNumberOfVertices();
    std::map<SimplexId, SimplexId> extremas;

if(!db_omp) { //#ifndef TTK_ENABLE_OPENMP
    // vertices that may still be compressed
    std::vector<SimplexId>* activeVertices = new std::vector<SimplexId>();
    SimplexId nActiveVertices;

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

    delete activeVertices;
} else { //#else // TTK_ENABLE_OPENMP

    #pragma omp parallel num_threads(threadNumber_)
    {
      std::vector<SimplexId>* lActiveVertices = new std::vector<SimplexId>();

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
            lActiveVertices->push_back(i);
        } else {
          #pragma omp critical
          {
            critMap.insert({numberOfExtrema, i});
            extremas.insert({i, numberOfExtrema});
            numberOfExtrema += 1;
          }
        }

        if(computeNeighbor) {
          neighbor[i] = mi;
        }
      }

      size_t lnActiveVertices = lActiveVertices->size();
      std::vector<SimplexId>* newActiveVert = new std::vector<SimplexId>(); 

      // compress paths until no changes occur
      while(lnActiveVertices > 0) {

        for(size_t i = 0; i < lnActiveVertices; i++) {
          SimplexId &v = (*lActiveVertices)[i];
          SimplexId &vo = manifold[v];

          vo = manifold[vo];

          // check if fully compressed
          if(vo != manifold[vo]) {
            lActiveVertices->push_back(v);
          }
        }

        delete lActiveVertices;
        lActiveVertices = newActiveVert;
        lnActiveVertices = lActiveVertices->size();

        newActiveVert = new std::vector<SimplexId>();
      }

      delete newActiveVert;

      #pragma omp single
      {
        int index = 0;
        for (auto & e: extremas) {
          e.second = index++;
        }
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
  SimplexId &numberOfMinima,
  SimplexId &numberOfMaxima,
  const SimplexId *const orderArr,
  const triangulationType &triangulation) const {

  ttk::Timer localTimer;

  this->printMsg(ttk::debug::Separator::L1);
  // print the progress of the current subprocedure (currently 0%)
  this->printMsg("Computing Manifolds",
                 0, // progress form 0-1
                 0, // elapsed time so far
                 this->threadNumber_, ttk::debug::LineMode::REPLACE);

  numberOfMinima = 0;
  numberOfMaxima = 0;

  const bool computeCriticalPoints = criticalPoints != nullptr;
  const bool computeNeighbor =
    (ascNeighbor != nullptr || dscNeighbor != nullptr);

  /* compute the Descending Maifold iterating over each vertex, searching
   * the biggest neighbor and compressing its path to its maximum */
  const SimplexId nVertices = triangulation.getNumberOfVertices();

  // vertices that may still be compressed
  
  
  std::map<SimplexId, int> minMap, maxMap;

if(!db_omp) { //#ifndef TTK_ENABLE_OPENMP
  SimplexId nActiveVertices;
  std::vector<SimplexId>* activeVertices = new std::vector<SimplexId>();
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

    if(computeNeighbor) {
      ascNeighbor[i] = ami;
      dscNeighbor[i] = dmi;
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

  delete activeVertices;

} else { //#else // TTK_ENABLE_OPENMP

  #pragma omp parallel num_threads(threadNumber_)
  {
    std::vector<SimplexId>* lActiveVertices = new std::vector<SimplexId>();

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
        lActiveVertices->push_back(i);
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

      if(computeNeighbor) {
        ascNeighbor[i] = ami;
        dscNeighbor[i] = dmi;
      }
    }

    size_t lnActiveVertices = lActiveVertices->size();
    std::vector<SimplexId>* newActiveVert = new std::vector<SimplexId>();

    // compress paths until no changes occur
    while(lnActiveVertices > 0) {
      for(size_t i = 0; i < lnActiveVertices; i++) {
        SimplexId &v = (*lActiveVertices)[i];
        SimplexId &vDsc = dscManifold[v];
        SimplexId &vAsc = ascManifold[v];

        // compress paths
        vDsc = dscManifold[vDsc];
        vAsc = ascManifold[vAsc];

        // check if fully compressed
        if(vDsc != dscManifold[vDsc] || vAsc != ascManifold[vAsc]) {
          newActiveVert->push_back(v);
        }
      }

      delete lActiveVertices;
      lActiveVertices = newActiveVert;
      lnActiveVertices = lActiveVertices->size();

      newActiveVert = new std::vector<SimplexId>();
    }

    delete newActiveVert;

    #pragma omp single nowait
    {
      int index = 0;
      for (auto & e: maxMap) {
        e.second = index++;
      }
    }

    #pragma omp single
    {
      int index = 0;
      for (auto & e: minMap) {
        e.second = index++;
      }
    }

    #pragma omp for schedule(static)
    for(SimplexId i = 0; i < nVertices; i++) {
      dscManifold[i] = maxMap[dscManifold[i]];
      ascManifold[i] = minMap[ascManifold[i]];
    }
  }

} //#endif // TTK_ENABLE_OPENMP

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
  std::vector<std::array<float, 6>> &edgePos,
  std::vector<long long> &msLabels,
  std::unordered_set<SimplexId> &sep2VertSet,
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
    SimplexId vertices[3];
    triangulation.getTriangleVertex(tri, 0, vertices[0]);
    triangulation.getTriangleVertex(tri, 1, vertices[1]);
    triangulation.getTriangleVertex(tri, 2, vertices[2]);

    const SimplexId msm[3] = {
      morseSmaleManifold[vertices[0]], morseSmaleManifold[vertices[1]],
      morseSmaleManifold[vertices[2]]};

    const unsigned char index0 = (msm[0] == msm[1]) ? 0x00 : 0x04; // 0 : 1
    const unsigned char index1 = (msm[0] == msm[2]) ? 0x00 :       // 0
                                 (msm[1] == msm[2]) ? 0x01 : 0x02; // 1 : 2

    const unsigned char lookupIndex = index0 | index1;

    if(triangleLookupIsMultiLabel[lookupIndex]) {
      if(triangleLookupIs2Label[lookupIndex]) {

        const int * edgeVerts = triangleLookupEdgeVerts[lookupIndex];

        float edgeCenters[2][3];
        getEdgeIncenter(vertices[edgeVerts[0]], vertices[edgeVerts[1]],
          edgeCenters[0], triangulation);
        getEdgeIncenter(vertices[edgeVerts[2]], vertices[edgeVerts[3]],
          edgeCenters[1], triangulation);

        edgePos.push_back({
          edgeCenters[0][0], edgeCenters[0][1], edgeCenters[0][2],
          edgeCenters[1][0], edgeCenters[1][1], edgeCenters[1][2]});

        const long long sparseID = getSparseId(
          msm[edgeVerts[0]], msm[edgeVerts[1]], numMSCRegions);

        msLabels.push_back(sparseID);

        if(triangulation.isVertexOnBoundary(vertices[0]))
          sep2VertSet.insert(vertices[0]);
        if(triangulation.isVertexOnBoundary(vertices[1]))
          sep2VertSet.insert(vertices[1]);
        if(triangulation.isVertexOnBoundary(vertices[2]))
          sep2VertSet.insert(vertices[2]);
      } else {
        float edgeCenters[4][3];
        getEdgeIncenter(
          vertices[0], vertices[1], edgeCenters[0], triangulation);
        getEdgeIncenter(
          vertices[1], vertices[2], edgeCenters[1], triangulation);
        getEdgeIncenter(
          vertices[2], vertices[0], edgeCenters[2], triangulation);
        triangulation.getTriangleIncenter(tri, edgeCenters[3]);

        edgePos.push_back({
          edgeCenters[0][0], edgeCenters[0][1], edgeCenters[0][2],
          edgeCenters[3][0], edgeCenters[3][1], edgeCenters[3][2]});

        edgePos.push_back({
          edgeCenters[1][0], edgeCenters[1][1], edgeCenters[1][2],
          edgeCenters[3][0], edgeCenters[3][1], edgeCenters[3][2]});

        edgePos.push_back({
          edgeCenters[2][0], edgeCenters[2][1], edgeCenters[2][2],
          edgeCenters[3][0], edgeCenters[3][1], edgeCenters[3][2]});

        const long long sparseID[3] = {
        getSparseId(msm[0], msm[1], numMSCRegions),
        getSparseId(msm[1], msm[2], numMSCRegions),
        getSparseId(msm[2], msm[0], numMSCRegions)};

        msLabels.push_back(sparseID[0]);
        msLabels.push_back(sparseID[1]);
        msLabels.push_back(sparseID[2]);

        sep2VertSet.insert(vertices[0]);
        sep2VertSet.insert(vertices[1]);
        sep2VertSet.insert(vertices[2]);
      }
    }
  }

  this->printMsg("Computed 1-separatrix pieces", 1, localTimer.getElapsedTime(),
                 this->threadNumber_);

  return 1;
}

template <typename triangulationType>
int ttk::MorseSmaleSegmentationPL::findSaddlesFromCadidates(
  std::vector<std::pair<SimplexId, char>> &criticalPoints,
  const std::unordered_set<SimplexId> &sep2VertSet,
  const triangulationType &triangulation) const {

  ttk::Timer localTimer;

  const int dim = triangulation.getDimensionality();

  const SimplexId * orderArr
    = static_cast<const SimplexId *>(inputOrderField_);

  //int numSaddles = 0;

  this->printMsg("Computing for saddles",
                 0, localTimer.getElapsedTime(), this->threadNumber_,
                 ttk::debug::LineMode::REPLACE);

  ScalarFieldCriticalPoints sfcp;
  sfcp.setDomainDimension(dim);
  
if(!db_omp) {
  for(SimplexId candidate: sep2VertSet) {
    const CriticalType critC = (CriticalType)sfcp.getCriticalType(
      candidate, orderArr, (AbstractTriangulation*)&triangulation);

    if(critC == CriticalType::Saddle1) {
      criticalPoints.push_back(std::make_pair(candidate, 1));
    } else if(critC == CriticalType::Saddle2) {
      criticalPoints.push_back(std::make_pair(candidate, 2));
    }
  }
} else {
  #pragma omp parallel num_threads(threadNumber_)
  {
    #pragma omp single
    {
      for(auto it = sep2VertSet.begin();
        it != sep2VertSet.end(); it++) {
        #pragma omp task
        {
          const CriticalType critC = (CriticalType)sfcp.getCriticalType(
            *it, orderArr, (AbstractTriangulation*)&triangulation);
          
          if(critC == CriticalType::Saddle1) {
            criticalPoints.push_back(std::make_pair(*it, 1));
          } else if(critC == CriticalType::Saddle2) {
            criticalPoints.push_back(std::make_pair(*it, 2));
          }
        }
      }
    }
  }
}

  this->printMsg("Computed saddles", 1, localTimer.getElapsedTime(),
    this->threadNumber_);

  return 1;
}

int ttk::MorseSmaleSegmentationPL::setSeparatrices2_2D(
  const std::vector<std::array<float, 6>> &edgePos,
  const std::vector<long long> &msLabels,
  std::map<long long, SimplexId> &msLabelMap) const {

  ttk::Timer localTimer;

  this->printMsg("Writing 2-separatrices", 0, localTimer.getElapsedTime(),
                  this->threadNumber_, ttk::debug::LineMode::REPLACE);

  const int numEdges = edgePos.size();

  size_t npoints = numEdges * 2;
  size_t numCells = numEdges;

  // resize arrays
  outputSeparatrices2_points_->resize(3 * npoints);
  auto &points = *outputSeparatrices2_points_;
  outputSeparatrices2_cells_connectivity_->resize(2 * numCells);
  auto &cellsConn = *outputSeparatrices2_cells_connectivity_;
  outputSeparatrices2_cells_cases_->resize(numEdges);
  auto &cellsMSCIds = *outputSeparatrices2_cells_mscIds_;
  outputSeparatrices2_cells_mscIds_->resize(numEdges);
  auto &cellsCase = *outputSeparatrices2_cells_cases_;
  outputSeparatrices2_cells_cases_->resize(numEdges);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_) schedule(static)
#endif
  for(int edge = 0; edge < numEdges; ++edge) {
    auto &edgeP = edgePos[edge];
    points[6 * edge    ] = edgeP[0];
    points[6 * edge + 1] = edgeP[1];
    points[6 * edge + 2] = edgeP[2];

    points[6 * edge + 3] = edgeP[3];
    points[6 * edge + 4] = edgeP[4];
    points[6 * edge + 5] = edgeP[5];

    cellsConn[2 * edge    ] = 2 * edge;
    cellsConn[2 * edge + 1] = 2 * edge + 1;

    cellsMSCIds[edge] = msLabelMap[msLabels[edge]];
    cellsCase[edge] = msLabelMap[msLabels[edge]];
  }

  // update pointers
  *outputSeparatrices2_numberOfPoints_ = npoints;
  *outputSeparatrices2_numberOfCells_ = numCells;

  this->printMsg("WroteX " + std::to_string(numEdges) + " 2-separatrices",
                 1, localTimer.getElapsedTime(), this->threadNumber_);

  return 1;
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
  this->printMsg("Computing 2-Separatrices 3D[Walls]",
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

      float vertPos[4][3];
      triangulation.getVertexPoint(
        vertices[0], vertPos[0][0], vertPos[0][1], vertPos[0][2]);
      triangulation.getVertexPoint(
        vertices[1], vertPos[1][0], vertPos[1][1], vertPos[1][2]);
      triangulation.getVertexPoint(
        vertices[2], vertPos[2][0], vertPos[2][1], vertPos[2][2]);
      triangulation.getVertexPoint(
        vertices[3], vertPos[3][0], vertPos[3][1], vertPos[3][2]);

      float edgeCenters[10][3];
      // 6 edge centers
      getCenter(vertPos[0], vertPos[1], edgeCenters[0]);
      getCenter(vertPos[0], vertPos[2], edgeCenters[1]);
      getCenter(vertPos[0], vertPos[3], edgeCenters[2]);
      getCenter(vertPos[1], vertPos[2], edgeCenters[3]);
      getCenter(vertPos[1], vertPos[3], edgeCenters[4]);
      getCenter(vertPos[2], vertPos[3], edgeCenters[5]);

      // 4 triangle centers
      getCenter(vertPos[0], vertPos[1], vertPos[2], edgeCenters[6]);
      getCenter(vertPos[0], vertPos[1], vertPos[3], edgeCenters[7]);
      getCenter(vertPos[0], vertPos[2], vertPos[3], edgeCenters[8]);
      getCenter(vertPos[1], vertPos[2], vertPos[3], edgeCenters[9]);

      if(tetEdgeIndices[0] == 10) { // 4 labels on tetraeder
        float tetCenter[3];
        getCenter(vertPos[0], vertPos[1], vertPos[2], vertPos[3], tetCenter);

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
  size_t triangleStartIndex[threadNumber_ + 1];
  #pragma omp parallel num_threads(threadNumber_)
  {
    const int tid = omp_get_thread_num();
    std::vector<std::array<float, 9>> trianglePosLocal;
    std::vector<SimplexId> caseDataLocal;
    std::vector<long long> msLabelsLocal;

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
        const int *tetEdgeIndices = tetraederLookup[lookupIndex];
        const int *tetVertLabel = tetraederLabelLookup[lookupIndex];

        float vertPos[4][3];
        triangulation.getVertexPoint(
          vertices[0], vertPos[0][0], vertPos[0][1], vertPos[0][2]);
        triangulation.getVertexPoint(
          vertices[1], vertPos[1][0], vertPos[1][1], vertPos[1][2]);
        triangulation.getVertexPoint(
          vertices[2], vertPos[2][0], vertPos[2][1], vertPos[2][2]);
        triangulation.getVertexPoint(
          vertices[3], vertPos[3][0], vertPos[3][1], vertPos[3][2]);

        float edgeCenters[10][3];
        // 6 edge centers
        getCenter(vertPos[0], vertPos[1], edgeCenters[0]);
        getCenter(vertPos[0], vertPos[2], edgeCenters[1]);
        getCenter(vertPos[0], vertPos[3], edgeCenters[2]);
        getCenter(vertPos[1], vertPos[2], edgeCenters[3]);
        getCenter(vertPos[1], vertPos[3], edgeCenters[4]);
        getCenter(vertPos[2], vertPos[3], edgeCenters[5]);

        // 4 triangle centers
        getCenter(vertPos[0], vertPos[1], vertPos[2], edgeCenters[6]);
        getCenter(vertPos[0], vertPos[1], vertPos[3], edgeCenters[7]);
        getCenter(vertPos[0], vertPos[2], vertPos[3], edgeCenters[8]);
        getCenter(vertPos[1], vertPos[2], vertPos[3], edgeCenters[9]);

        if(tetEdgeIndices[0] == 10) { // 4 labels on tetraeder
          float tetCenter[3];
          getCenter(vertPos[0], vertPos[1], vertPos[2], vertPos[3], tetCenter);

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

          trianglePosLocal.insert(trianglePosLocal.end(), {{
            edgeCenters[7][0], edgeCenters[7][1], edgeCenters[7][2], 
            edgeCenters[0][0], edgeCenters[0][1], edgeCenters[0][2], 
            tetCenter[0], tetCenter[1], tetCenter[2]
          },{
            edgeCenters[0][0], edgeCenters[0][1], edgeCenters[0][2], 
            edgeCenters[6][0], edgeCenters[6][1], edgeCenters[6][2], 
            tetCenter[0], tetCenter[1], tetCenter[2]
          },{
            edgeCenters[8][0], edgeCenters[8][1], edgeCenters[8][2],
            edgeCenters[1][0], edgeCenters[1][1], edgeCenters[1][2],
            tetCenter[0], tetCenter[1], tetCenter[2]
          },{
            edgeCenters[1][0], edgeCenters[1][1], edgeCenters[1][2],
            edgeCenters[6][0], edgeCenters[6][1], edgeCenters[6][2],
            tetCenter[0], tetCenter[1], tetCenter[2]
          },{
            edgeCenters[8][0], edgeCenters[8][1], edgeCenters[8][2], 
            edgeCenters[2][0], edgeCenters[2][1], edgeCenters[2][2], 
            tetCenter[0], tetCenter[1], tetCenter[2]
          },{
            edgeCenters[2][0], edgeCenters[2][1], edgeCenters[2][2], 
            edgeCenters[7][0], edgeCenters[7][1], edgeCenters[7][2], 
            tetCenter[0], tetCenter[1], tetCenter[2]
          },{
            edgeCenters[6][0], edgeCenters[6][1], edgeCenters[6][2],
            edgeCenters[3][0], edgeCenters[3][1], edgeCenters[3][2],
            tetCenter[0], tetCenter[1], tetCenter[2]
          },{
            edgeCenters[3][0], edgeCenters[3][1], edgeCenters[3][2],
            edgeCenters[9][0], edgeCenters[9][1], edgeCenters[9][2],
            tetCenter[0], tetCenter[1], tetCenter[2]
          },{
            edgeCenters[7][0], edgeCenters[7][1], edgeCenters[7][2],
            edgeCenters[4][0], edgeCenters[4][1], edgeCenters[4][2],
            tetCenter[0], tetCenter[1], tetCenter[2]
          },{
            edgeCenters[4][0], edgeCenters[4][1], edgeCenters[4][2],
            edgeCenters[9][0], edgeCenters[9][1], edgeCenters[9][2],
            tetCenter[0], tetCenter[1], tetCenter[2]
          },{
            edgeCenters[9][0], edgeCenters[9][1], edgeCenters[9][2], 
            edgeCenters[5][0], edgeCenters[5][1], edgeCenters[5][2], 
            tetCenter[0], tetCenter[1], tetCenter[2]
          },{
            edgeCenters[5][0], edgeCenters[5][1], edgeCenters[5][2],
            edgeCenters[8][0], edgeCenters[8][1], edgeCenters[8][2],
            tetCenter[0], tetCenter[1], tetCenter[2]
          }});        

          caseDataLocal.insert(caseDataLocal.end(), {
            lookupIndex, lookupIndex, lookupIndex,
            lookupIndex, lookupIndex, lookupIndex,
            lookupIndex, lookupIndex, lookupIndex,
            lookupIndex, lookupIndex, lookupIndex}
          );
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

    triangleStartIndex[tid + 1] = trianglePosLocal.size();

    #pragma omp barrier

    #pragma omp single
    {
      triangleStartIndex[0] = 0;

      // Count triangle number and create iterator start indices
      for(int t = 1; t < threadNumber_ + 1; ++t) {
        triangleStartIndex[t] += triangleStartIndex[t-1];
      }

      trianglePos.resize(triangleStartIndex[threadNumber_]);
    }

    std::copy(trianglePosLocal.begin(), trianglePosLocal.end(),
      trianglePos.begin() + triangleStartIndex[tid]);

    #pragma omp single
    {
      caseData.resize(triangleStartIndex[threadNumber_]);
    }

    std::copy(caseDataLocal.begin(), caseDataLocal.end(),
      caseData.begin() + triangleStartIndex[tid]);

    #pragma omp single
    {
      msLabels.resize(triangleStartIndex[threadNumber_]);
    }

    std::copy(msLabelsLocal.begin(), msLabelsLocal.end(),
      msLabels.begin() + triangleStartIndex[tid]);
  }
} // #if TTK_OMPENMP
  this->printMsg("Computed 2-Separatrices 3D[Walls]",
                 1, localTimer.getElapsedTime(), this->threadNumber_);

  return 0;
}

template <typename triangulationType>
int ttk::MorseSmaleSegmentationPL::computeSeparatrices2_3D_split(
  std::vector<std::array<float, 9>> &trianglePos,
  std::vector<SimplexId> &caseData,
  std::vector<long long> &msLabels,
  const SimplexId &numAscRegions,
  const SimplexId &numDscRegions,
  const SimplexId *const ascManifold,
  const SimplexId *const dscManifold,
  const triangulationType &triangulation) const {

  ttk::Timer localTimer;

  // print the progress of the current subprocedure (currently 0%)
  this->printMsg("Computing 2-Separatrices 3D[Split]",
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
  const long long ascSparseIDOffset = (numAscRegions * numAscRegions);

//#ifndef TTK_ENABLE_OPENMP
if(!db_omp) {
  for(SimplexId tet = 0; tet < numTetra; tet++) {
    SimplexId vertices[4];
    triangulation.getCellVertex(tet, 0, vertices[0]);
    triangulation.getCellVertex(tet, 1, vertices[1]);
    triangulation.getCellVertex(tet, 2, vertices[2]);
    triangulation.getCellVertex(tet, 3, vertices[3]);

    const SimplexId asc[4] = {
      ascManifold[vertices[0]], ascManifold[vertices[1]],
      ascManifold[vertices[2]], ascManifold[vertices[3]]};
    
    const SimplexId dsc[4] = {
      dscManifold[vertices[0]], dscManifold[vertices[1]],
      dscManifold[vertices[2]], dscManifold[vertices[3]]};

    const unsigned char index1A = (asc[0] == asc[1]) ? 0x00 : 0x10; // 0 : 1
    const unsigned char index2A = (asc[0] == asc[2]) ? 0x00 :       // 0
                                  (asc[1] == asc[2]) ? 0x04 : 0x08; // 1 : 2
    const unsigned char index3A = (asc[0] == asc[3]) ? 0x00 :       // 0
                                  (asc[1] == asc[3]) ? 0x01 :       // 1
                                  (asc[2] == asc[3]) ? 0x02 : 0x03; // 2 : 3

    const unsigned char index1D = (dsc[0] == dsc[1]) ? 0x00 : 0x10; // 0 : 1
    const unsigned char index2D = (dsc[0] == dsc[2]) ? 0x00 :       // 0
                                  (dsc[1] == dsc[2]) ? 0x04 : 0x08; // 1 : 2
    const unsigned char index3D = (dsc[0] == dsc[3]) ? 0x00 :       // 0
                                  (dsc[1] == dsc[3]) ? 0x01 :       // 1
                                  (dsc[2] == dsc[3]) ? 0x02 : 0x03; // 2 : 3

    const unsigned char lookupIndexA = index1A | index2A | index3A;
    const int *tetEdgeIndicesA = tetraederLookup[lookupIndexA];
    const int *tetVertLabelA = tetraederLabelLookup[lookupIndexA];

    const unsigned char lookupIndexD = index1D | index2D | index3D;
    const int *tetEdgeIndicesD = tetraederLookup[lookupIndexD];
    const int *tetVertLabelD = tetraederLabelLookup[lookupIndexD];

    if(tetraederLookupIsMultiLabel[lookupIndexA] 
    || tetraederLookupIsMultiLabel[lookupIndexD]) {

      float vertPos[4][3];
      triangulation.getVertexPoint(
        vertices[0], vertPos[0][0], vertPos[0][1], vertPos[0][2]);
      triangulation.getVertexPoint(
        vertices[1], vertPos[1][0], vertPos[1][1], vertPos[1][2]);
      triangulation.getVertexPoint(
        vertices[2], vertPos[2][0], vertPos[2][1], vertPos[2][2]);
      triangulation.getVertexPoint(
        vertices[3], vertPos[3][0], vertPos[3][1], vertPos[3][2]);

      float edgeCenters[10][3];
      // 6 edge centers
      getCenter(vertPos[0], vertPos[1], edgeCenters[0]);
      getCenter(vertPos[0], vertPos[2], edgeCenters[1]);
      getCenter(vertPos[0], vertPos[3], edgeCenters[2]);
      getCenter(vertPos[1], vertPos[2], edgeCenters[3]);
      getCenter(vertPos[1], vertPos[3], edgeCenters[4]);
      getCenter(vertPos[2], vertPos[3], edgeCenters[5]);

      // 4 triangle centers
      getCenter(vertPos[0], vertPos[1], vertPos[2], edgeCenters[6]);
      getCenter(vertPos[0], vertPos[1], vertPos[3], edgeCenters[7]);
      getCenter(vertPos[0], vertPos[2], vertPos[3], edgeCenters[8]);
      getCenter(vertPos[1], vertPos[2], vertPos[3], edgeCenters[9]);

      if(tetEdgeIndicesA[0] != -1) {
        if(tetEdgeIndicesA[0] == 10) { // 4 labels on tetraeder
          float tetCenter[3];
          getCenter(vertPos[0], vertPos[1], vertPos[2], vertPos[3], tetCenter);

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

          long long sparseAscIds[] = {
            getSparseId(asc[0], asc[1], numAscRegions),
            getSparseId(asc[0], asc[2], numAscRegions),
            getSparseId(asc[0], asc[3], numAscRegions),
            getSparseId(asc[1], asc[2], numAscRegions),
            getSparseId(asc[1], asc[3], numAscRegions),
            getSparseId(asc[2], asc[3], numAscRegions)
          };
          msLabels.insert(msLabels.end(), {
            sparseAscIds[0], sparseAscIds[0],
            sparseAscIds[1], sparseAscIds[1],
            sparseAscIds[2], sparseAscIds[2],
            sparseAscIds[3], sparseAscIds[3],
            sparseAscIds[4], sparseAscIds[4],
            sparseAscIds[5], sparseAscIds[5]
          });

          caseData.insert(caseData.end(), {
            lookupIndexA, lookupIndexA, lookupIndexA,
            lookupIndexA, lookupIndexA, lookupIndexA,
            lookupIndexA, lookupIndexA, lookupIndexA,
            lookupIndexA, lookupIndexA, lookupIndexA}
          );
        } else { // 2 or 3 labels on tetraeder
          const size_t numTris = tetraederNumTriangles[lookupIndexA];
          long long sparseIds[numTris];
          for(size_t t = 0; t < numTris; ++t) {
            trianglePos.push_back({
              edgeCenters[tetEdgeIndicesA[(t * 3)    ]][0],
              edgeCenters[tetEdgeIndicesA[(t * 3)    ]][1],
              edgeCenters[tetEdgeIndicesA[(t * 3)    ]][2],
              edgeCenters[tetEdgeIndicesA[(t * 3) + 1]][0],
              edgeCenters[tetEdgeIndicesA[(t * 3) + 1]][1],
              edgeCenters[tetEdgeIndicesA[(t * 3) + 1]][2],
              edgeCenters[tetEdgeIndicesA[(t * 3) + 2]][0],
              edgeCenters[tetEdgeIndicesA[(t * 3) + 2]][1],
              edgeCenters[tetEdgeIndicesA[(t * 3) + 2]][2]
            });

            sparseIds[t] = getSparseId(asc[tetVertLabelA[t * 2]],
                                       asc[tetVertLabelA[(t * 2) + 1]],
                                       numAscRegions);
            msLabels.push_back(sparseIds[t]);

            caseData.push_back(lookupIndexA);
          }
        }
      }
      if(tetEdgeIndicesD[0] != -1) {
        if(tetEdgeIndicesD[0] == 10) { // 4 labels on tetraeder
          float tetCenter[3];
          getCenter(vertPos[0], vertPos[1], vertPos[2], vertPos[3], tetCenter);

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

          long long sparseDscIds[] = {
            getSparseId(asc[0], asc[1], numDscRegions),
            getSparseId(asc[0], asc[2], numDscRegions),
            getSparseId(asc[0], asc[3], numDscRegions),
            getSparseId(asc[1], asc[2], numDscRegions),
            getSparseId(asc[1], asc[3], numDscRegions),
            getSparseId(asc[2], asc[3], numDscRegions)
          };
          msLabels.insert(msLabels.end(), {
            sparseDscIds[0], sparseDscIds[0],
            sparseDscIds[1], sparseDscIds[1],
            sparseDscIds[2], sparseDscIds[2],
            sparseDscIds[3], sparseDscIds[3],
            sparseDscIds[4], sparseDscIds[4],
            sparseDscIds[5], sparseDscIds[5]
          });

          caseData.insert(caseData.end(), {
            lookupIndexD, lookupIndexD, lookupIndexD,
            lookupIndexD, lookupIndexD, lookupIndexD,
            lookupIndexD, lookupIndexD, lookupIndexD,
            lookupIndexD, lookupIndexD, lookupIndexD}
          );
        } else { // 2 or 3 labels on tetraeder
          const size_t numTris = tetraederNumTriangles[lookupIndexD];
          long long sparseIds[numTris];
          for(size_t t = 0; t < numTris; ++t) {
            trianglePos.push_back({
              edgeCenters[tetEdgeIndicesD[(t * 3)    ]][0],
              edgeCenters[tetEdgeIndicesD[(t * 3)    ]][1],
              edgeCenters[tetEdgeIndicesD[(t * 3)    ]][2],
              edgeCenters[tetEdgeIndicesD[(t * 3) + 1]][0],
              edgeCenters[tetEdgeIndicesD[(t * 3) + 1]][1],
              edgeCenters[tetEdgeIndicesD[(t * 3) + 1]][2],
              edgeCenters[tetEdgeIndicesD[(t * 3) + 2]][0],
              edgeCenters[tetEdgeIndicesD[(t * 3) + 2]][1],
              edgeCenters[tetEdgeIndicesD[(t * 3) + 2]][2]
            });

            sparseIds[t] = getSparseId(dsc[tetVertLabelD[t * 2]],
                                       dsc[tetVertLabelD[(t * 2) + 1]],
                                       numDscRegions);
            msLabels.push_back(sparseIds[t] + ascSparseIDOffset);

            caseData.push_back(lookupIndexD);
          }
        }
      }
    }
  }
//#else // TTK_ENABLE_OPENMP
} else {
  size_t triangleStartIndex[threadNumber_ + 1];
  #pragma omp parallel num_threads(threadNumber_)
  {
    const int tid = omp_get_thread_num();
    std::vector<std::array<float, 9>> trianglePosLocal;
    std::vector<SimplexId> caseDataLocal;
    std::vector<long long> msLabelsLocal;

    #pragma omp for schedule(static) nowait
    for(SimplexId tet = 0; tet < numTetra; tet++) {
      SimplexId vertices[4];
      triangulation.getCellVertex(tet, 0, vertices[0]);
      triangulation.getCellVertex(tet, 1, vertices[1]);
      triangulation.getCellVertex(tet, 2, vertices[2]);
      triangulation.getCellVertex(tet, 3, vertices[3]);

      const SimplexId asc[4] = {
        ascManifold[vertices[0]], ascManifold[vertices[1]],
        ascManifold[vertices[2]], ascManifold[vertices[3]]};
      
      const SimplexId dsc[4] = {
        dscManifold[vertices[0]], dscManifold[vertices[1]],
        dscManifold[vertices[2]], dscManifold[vertices[3]]};

      const unsigned char index1A = (asc[0] == asc[1]) ? 0x00 : 0x10; // 0 : 1
      const unsigned char index2A = (asc[0] == asc[2]) ? 0x00 :       // 0
                                    (asc[1] == asc[2]) ? 0x04 : 0x08; // 1 : 2
      const unsigned char index3A = (asc[0] == asc[3]) ? 0x00 :       // 0
                                    (asc[1] == asc[3]) ? 0x01 :       // 1
                                    (asc[2] == asc[3]) ? 0x02 : 0x03; // 2 : 3

      const unsigned char index1D = (dsc[0] == dsc[1]) ? 0x00 : 0x10; // 0 : 1
      const unsigned char index2D = (dsc[0] == dsc[2]) ? 0x00 :       // 0
                                    (dsc[1] == dsc[2]) ? 0x04 : 0x08; // 1 : 2
      const unsigned char index3D = (dsc[0] == dsc[3]) ? 0x00 :       // 0
                                    (dsc[1] == dsc[3]) ? 0x01 :       // 1
                                    (dsc[2] == dsc[3]) ? 0x02 : 0x03; // 2 : 3

      const unsigned char lookupIndexA = index1A | index2A | index3A;
      const int *tetEdgeIndicesA = tetraederLookup[lookupIndexA];
      const int *tetVertLabelA = tetraederLabelLookup[lookupIndexA];

      const unsigned char lookupIndexD = index1D | index2D | index3D;
      const int *tetEdgeIndicesD = tetraederLookup[lookupIndexD];
      const int *tetVertLabelD = tetraederLabelLookup[lookupIndexD];

      if(tetraederLookupIsMultiLabel[lookupIndexA] 
      || tetraederLookupIsMultiLabel[lookupIndexD]) {
        float vertPos[4][3];
        triangulation.getVertexPoint(
          vertices[0], vertPos[0][0], vertPos[0][1], vertPos[0][2]);
        triangulation.getVertexPoint(
          vertices[1], vertPos[1][0], vertPos[1][1], vertPos[1][2]);
        triangulation.getVertexPoint(
          vertices[2], vertPos[2][0], vertPos[2][1], vertPos[2][2]);
        triangulation.getVertexPoint(
          vertices[3], vertPos[3][0], vertPos[3][1], vertPos[3][2]);

        float edgeCenters[10][3];
        // 6 edge centers
        getCenter(vertPos[0], vertPos[1], edgeCenters[0]);
        getCenter(vertPos[0], vertPos[2], edgeCenters[1]);
        getCenter(vertPos[0], vertPos[3], edgeCenters[2]);
        getCenter(vertPos[1], vertPos[2], edgeCenters[3]);
        getCenter(vertPos[1], vertPos[3], edgeCenters[4]);
        getCenter(vertPos[2], vertPos[3], edgeCenters[5]);

        // 4 triangle centers
        getCenter(vertPos[0], vertPos[1], vertPos[2], edgeCenters[6]);
        getCenter(vertPos[0], vertPos[1], vertPos[3], edgeCenters[7]);
        getCenter(vertPos[0], vertPos[2], vertPos[3], edgeCenters[8]);
        getCenter(vertPos[1], vertPos[2], vertPos[3], edgeCenters[9]);

        if(tetEdgeIndicesA[0] != -1) {
          if(tetEdgeIndicesA[0] == 10) { // 4 labels on tetraeder
            float tetCenter[3];
            getCenter(
              vertPos[0], vertPos[1], vertPos[2], vertPos[3], tetCenter);

            // vertex 0
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

            long long sparseAscIds[] = {
              getSparseId(asc[0], asc[1], numAscRegions),
              getSparseId(asc[0], asc[2], numAscRegions),
              getSparseId(asc[0], asc[3], numAscRegions),
              getSparseId(asc[1], asc[2], numAscRegions),
              getSparseId(asc[1], asc[3], numAscRegions),
              getSparseId(asc[2], asc[3], numAscRegions)
            };
            msLabelsLocal.insert(msLabelsLocal.end(), {
              sparseAscIds[0], sparseAscIds[0],
              sparseAscIds[1], sparseAscIds[1],
              sparseAscIds[2], sparseAscIds[2],
              sparseAscIds[3], sparseAscIds[3],
              sparseAscIds[4], sparseAscIds[4],
              sparseAscIds[5], sparseAscIds[5]
            });

            caseDataLocal.insert(caseDataLocal.end(), {
              lookupIndexA, lookupIndexA, lookupIndexA,
              lookupIndexA, lookupIndexA, lookupIndexA,
              lookupIndexA, lookupIndexA, lookupIndexA,
              lookupIndexA, lookupIndexA, lookupIndexA}
            );
          } else { // 2 or 3 labels on tetraeder
            const size_t numTris = tetraederNumTriangles[lookupIndexA];
            long long sparseIds[numTris];
            for(size_t t = 0; t < numTris; ++t) {
              trianglePosLocal.push_back({
                edgeCenters[tetEdgeIndicesA[(t * 3)    ]][0],
                edgeCenters[tetEdgeIndicesA[(t * 3)    ]][1],
                edgeCenters[tetEdgeIndicesA[(t * 3)    ]][2],
                edgeCenters[tetEdgeIndicesA[(t * 3) + 1]][0],
                edgeCenters[tetEdgeIndicesA[(t * 3) + 1]][1],
                edgeCenters[tetEdgeIndicesA[(t * 3) + 1]][2],
                edgeCenters[tetEdgeIndicesA[(t * 3) + 2]][0],
                edgeCenters[tetEdgeIndicesA[(t * 3) + 2]][1],
                edgeCenters[tetEdgeIndicesA[(t * 3) + 2]][2]
              });

              sparseIds[t] = getSparseId(asc[tetVertLabelA[t * 2]],
                                         asc[tetVertLabelA[(t * 2) + 1]],
                                         numAscRegions);
              msLabelsLocal.push_back(sparseIds[t]);

              //caseDataLocal.push_back(lookupIndexA);
              caseDataLocal.push_back(0);
            }
          }
        }
        if(tetEdgeIndicesD[0] != -1) {
          if(tetEdgeIndicesD[0] == 10) { // 4 labels on tetraeder
            float tetCenter[3];
            getCenter(
              vertPos[0], vertPos[1], vertPos[2], vertPos[3], tetCenter);

            // vertex 0
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

            long long sparseDscIds[] = {
              getSparseId(asc[0], asc[1], numDscRegions),
              getSparseId(asc[0], asc[2], numDscRegions),
              getSparseId(asc[0], asc[3], numDscRegions),
              getSparseId(asc[1], asc[2], numDscRegions),
              getSparseId(asc[1], asc[3], numDscRegions),
              getSparseId(asc[2], asc[3], numDscRegions)
            };
            msLabelsLocal.insert(msLabelsLocal.end(), {
              sparseDscIds[0], sparseDscIds[0],
              sparseDscIds[1], sparseDscIds[1],
              sparseDscIds[2], sparseDscIds[2],
              sparseDscIds[3], sparseDscIds[3],
              sparseDscIds[4], sparseDscIds[4],
              sparseDscIds[5], sparseDscIds[5]
            });

            caseDataLocal.insert(caseDataLocal.end(), {
              lookupIndexD, lookupIndexD, lookupIndexD,
              lookupIndexD, lookupIndexD, lookupIndexD,
              lookupIndexD, lookupIndexD, lookupIndexD,
              lookupIndexD, lookupIndexD, lookupIndexD}
            );
          } else { // 2 or 3 labels on tetraeder
            const size_t numTris = tetraederNumTriangles[lookupIndexD];
            long long sparseIds[numTris];
            for(size_t t = 0; t < numTris; ++t) {
              trianglePosLocal.push_back({
                edgeCenters[tetEdgeIndicesD[(t * 3)    ]][0],
                edgeCenters[tetEdgeIndicesD[(t * 3)    ]][1],
                edgeCenters[tetEdgeIndicesD[(t * 3)    ]][2],
                edgeCenters[tetEdgeIndicesD[(t * 3) + 1]][0],
                edgeCenters[tetEdgeIndicesD[(t * 3) + 1]][1],
                edgeCenters[tetEdgeIndicesD[(t * 3) + 1]][2],
                edgeCenters[tetEdgeIndicesD[(t * 3) + 2]][0],
                edgeCenters[tetEdgeIndicesD[(t * 3) + 2]][1],
                edgeCenters[tetEdgeIndicesD[(t * 3) + 2]][2]
              });

              sparseIds[t] = getSparseId(dsc[tetVertLabelD[t * 2]],
                                         dsc[tetVertLabelD[(t * 2) + 1]],
                                         numDscRegions);
              msLabelsLocal.push_back(sparseIds[t] + ascSparseIDOffset);

              //caseDataLocal.push_back(lookupIndexD);
              caseDataLocal.push_back(1);
            }
          }
        }
      }
    }

    triangleStartIndex[tid + 1] = trianglePosLocal.size();

    #pragma omp barrier

    #pragma omp single
    {
      triangleStartIndex[0] = 0;

      // Count triangle number and create iterator start indices
      for(int t = 1; t < threadNumber_ + 1; ++t) {
        triangleStartIndex[t] += triangleStartIndex[t-1];
      }

      trianglePos.resize(triangleStartIndex[threadNumber_]);
    }

    std::copy(trianglePosLocal.begin(), trianglePosLocal.end(),
      trianglePos.begin() + triangleStartIndex[tid]);

    #pragma omp single
    {
      caseData.resize(triangleStartIndex[threadNumber_]);
    }

    std::copy(caseDataLocal.begin(), caseDataLocal.end(),
      caseData.begin() + triangleStartIndex[tid]);

    #pragma omp single
    {
      msLabels.resize(triangleStartIndex[threadNumber_]);
    }

    std::copy(msLabelsLocal.begin(), msLabelsLocal.end(),
      msLabels.begin() + triangleStartIndex[tid]);
  }
} // #if TTK_OMPENMP
  this->printMsg("Computed 2-Separatrices 3D[Split]", 1,
    localTimer.getElapsedTime(), this->threadNumber_);

  return 0;
}

int ttk::MorseSmaleSegmentationPL::computeMSLabelMap(
  const std::vector<long long> &msLabels,
  std::map<long long, SimplexId> &msLabelMap) const {

  std::vector<long long> msLabelCopy(msLabels);

  // get unique "sparse region ids"
  TTK_PSORT(this->threadNumber_, msLabelCopy.begin(), msLabelCopy.end());
  const auto last = std::unique(msLabelCopy.begin(), msLabelCopy.end());
  msLabelCopy.erase(last, msLabelCopy.end());

  for(size_t i = 0; i < msLabelCopy.size(); ++i) {
    msLabelMap.insert({msLabelCopy[i], i});
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

  this->printMsg("Computing 2-Separatrices 3D[Fine]",
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

  const float diff = 0.02;
  const float d0 = 0.5 + diff;
  const float d1 = 0.5 - diff;
  const float dc = diff * 2 ;

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
            vert03[0], vert03[1], vert03[2]});
          trianglePos.push_back({
            vert10[0], vert10[1], vert10[2], 
            vert11[0], vert11[1], vert11[2], 
            vert13[0], vert13[1], vert13[2]});

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
          vert02[0], vert02[1], vert02[2]});
        trianglePos.push_back({
          vert00[0], vert00[1], vert00[2], 
          vert01[0], vert01[1], vert01[2], 
          vert03[0], vert03[1], vert03[2]});
        trianglePos.push_back({
          vert10[0], vert10[1], vert10[2], 
          vert11[0], vert11[1], vert11[2], 
          vert12[0], vert12[1], vert12[2]});
        trianglePos.push_back({
          vert20[0], vert20[1], vert20[2], 
          vert21[0], vert21[1], vert21[2], 
          vert22[0], vert22[1], vert22[2]});

        msLabels.push_back(msm[vIds[0]]);
        msLabels.push_back(msm[vIds[0]]);
        msLabels.push_back(msm[vIds[1]]);
        msLabels.push_back(msm[vIds[3]]);

        caseData.push_back(lookupIndex);
        caseData.push_back(lookupIndex);
        caseData.push_back(lookupIndex);
        caseData.push_back(lookupIndex);
      } else { // 4 labels
        float tetCenter[3];
        getCenter(vPos[0], vPos[1], vPos[2], vPos[3], tetCenter);

        // the 4 triangle centers
        float triCenter[4][3];
        getCenter(vPos[0], vPos[1], vPos[2], triCenter[0]);
        getCenter(vPos[0], vPos[1], vPos[3], triCenter[1]);
        getCenter(vPos[0], vPos[2], vPos[3], triCenter[2]);
        getCenter(vPos[1], vPos[2], vPos[3], triCenter[3]);

        float vert00[3], vert01[3], vert02[3], vert0tet[3],
              vert0t0[3], vert0t1[3], vert0t2[3], 
              vert10[3], vert11[3], vert12[3], vert1tet[3],
              vert1t0[3], vert1t1[3], vert1t2[3], 
              vert20[3], vert21[3], vert22[3], vert2tet[3],
              vert2t0[3], vert2t1[3], vert2t2[3], 
              vert30[3], vert31[3], vert32[3], vert3tet[3],
              vert3t0[3], vert3t1[3], vert3t2[3]; 

        interpolatePoints(vPos[0], vPos[1], dc, vert00);
        interpolatePoints(vPos[0], vPos[2], dc, vert01);
        interpolatePoints(vPos[0], vPos[3], dc, vert02);
        interpolatePoints(vPos[0], tetCenter, dc, vert0tet);
        interpolatePoints(vPos[0], triCenter[0], dc, vert0t0);
        interpolatePoints(vPos[0], triCenter[1], dc, vert0t1);
        interpolatePoints(vPos[0], triCenter[2], dc, vert0t2);

        interpolatePoints(vPos[1], vPos[0], dc, vert10);
        interpolatePoints(vPos[1], vPos[2], dc, vert11);
        interpolatePoints(vPos[1], vPos[3], dc, vert12);
        interpolatePoints(vPos[1], tetCenter, dc, vert1tet);
        interpolatePoints(vPos[1], triCenter[0], dc, vert1t0);
        interpolatePoints(vPos[1], triCenter[1], dc, vert1t1);
        interpolatePoints(vPos[1], triCenter[2], dc, vert1t2);

        interpolatePoints(vPos[2], vPos[0], dc, vert20);
        interpolatePoints(vPos[2], vPos[1], dc, vert21);
        interpolatePoints(vPos[2], vPos[3], dc, vert22);
        interpolatePoints(vPos[2], tetCenter, dc, vert2tet);
        interpolatePoints(vPos[2], triCenter[0], dc, vert2t0);
        interpolatePoints(vPos[2], triCenter[1], dc, vert2t1);
        interpolatePoints(vPos[2], triCenter[2], dc, vert2t2);

        interpolatePoints(vPos[3], vPos[0], dc, vert30);
        interpolatePoints(vPos[3], vPos[1], dc, vert31);
        interpolatePoints(vPos[3], vPos[2], dc, vert32);
        interpolatePoints(vPos[3], tetCenter, dc, vert3tet);
        interpolatePoints(vPos[3], triCenter[0], dc, vert3t0);
        interpolatePoints(vPos[3], triCenter[1], dc, vert3t1);
        interpolatePoints(vPos[3], triCenter[2], dc, vert3t2);

        // Label Vert 0
        trianglePos.push_back({
          vert00[0],   vert00[1],   vert00[2], 
          vert0t0[0],  vert0t0[1],  vert0t0[2], 
          vert0tet[0], vert0tet[1], vert0tet[2]});
        trianglePos.push_back({
          vert00[0],   vert00[1],   vert00[2], 
          vert0t1[0],  vert0t1[1],  vert0t1[2], 
          vert0tet[0], vert0tet[1], vert0tet[2]});
        trianglePos.push_back({
          vert01[0],   vert01[1],   vert01[2], 
          vert0t1[0],  vert0t1[1],  vert0t1[2], 
          vert0tet[0], vert0tet[1], vert0tet[2]});
        trianglePos.push_back({
          vert01[0],   vert01[1],   vert01[2], 
          vert0t2[0],  vert0t2[1],  vert0t2[2], 
          vert0tet[0], vert0tet[1], vert0tet[2]});
        trianglePos.push_back({
          vert02[0],   vert02[1],   vert02[2], 
          vert0t2[0],  vert0t2[1],  vert0t2[2], 
          vert0tet[0], vert0tet[1], vert0tet[2]});
        trianglePos.push_back({
          vert02[0],   vert02[1],   vert02[2], 
          vert0t0[0],  vert0t0[1],  vert0t0[2], 
          vert0tet[0], vert0tet[1], vert0tet[2]});

        
        // Label Vert 0
        trianglePos.push_back({
          vert10[0],   vert10[1],   vert10[2], 
          vert1t0[0],  vert1t0[1],  vert1t0[2], 
          vert1tet[0], vert1tet[1], vert1tet[2]});
        trianglePos.push_back({
          vert10[0],   vert10[1],   vert10[2], 
          vert1t1[0],  vert1t1[1],  vert1t1[2], 
          vert1tet[0], vert1tet[1], vert1tet[2]});
        trianglePos.push_back({
          vert11[0],   vert11[1],   vert11[2], 
          vert1t1[0],  vert1t1[1],  vert1t1[2], 
          vert1tet[0], vert1tet[1], vert1tet[2]});
        trianglePos.push_back({
          vert11[0],   vert11[1],   vert11[2], 
          vert1t2[0],  vert1t2[1],  vert1t2[2], 
          vert1tet[0], vert1tet[1], vert1tet[2]});
        trianglePos.push_back({
          vert12[0],   vert12[1],   vert12[2], 
          vert1t2[0],  vert1t2[1],  vert1t2[2], 
          vert1tet[0], vert1tet[1], vert1tet[2]});
        trianglePos.push_back({
          vert12[0],   vert12[1],   vert12[2], 
          vert1t0[0],  vert1t0[1],  vert1t0[2], 
          vert1tet[0], vert1tet[1], vert1tet[2]});


        // Label Vert 0
        trianglePos.push_back({
          vert20[0],   vert20[1],   vert20[2], 
          vert2t0[0],  vert2t0[1],  vert2t0[2], 
          vert2tet[0], vert2tet[1], vert2tet[2]});
        trianglePos.push_back({
          vert20[0],   vert20[1],   vert20[2], 
          vert2t1[0],  vert2t1[1],  vert2t1[2], 
          vert2tet[0], vert2tet[1], vert2tet[2]});
        trianglePos.push_back({
          vert21[0],   vert21[1],   vert21[2], 
          vert2t1[0],  vert2t1[1],  vert2t1[2], 
          vert2tet[0], vert2tet[1], vert2tet[2]});
        trianglePos.push_back({
          vert21[0],   vert21[1],   vert21[2], 
          vert2t2[0],  vert2t2[1],  vert2t2[2], 
          vert2tet[0], vert2tet[1], vert2tet[2]});
        trianglePos.push_back({
          vert22[0],   vert22[1],   vert22[2], 
          vert2t2[0],  vert2t2[1],  vert2t2[2], 
          vert2tet[0], vert2tet[1], vert2tet[2]});
        trianglePos.push_back({
          vert22[0],   vert22[1],   vert22[2], 
          vert2t0[0],  vert2t0[1],  vert2t0[2], 
          vert2tet[0], vert2tet[1], vert2tet[2]});


        // Label Vert 0
        trianglePos.push_back({
          vert30[0],   vert30[1],   vert30[2], 
          vert3t0[0],  vert3t0[1],  vert3t0[2], 
          vert3tet[0], vert3tet[1], vert3tet[2]});
        trianglePos.push_back({
          vert30[0],   vert30[1],   vert30[2], 
          vert3t1[0],  vert3t1[1],  vert3t1[2], 
          vert3tet[0], vert3tet[1], vert3tet[2]});
        trianglePos.push_back({
          vert31[0],   vert31[1],   vert31[2], 
          vert3t1[0],  vert3t1[1],  vert3t1[2], 
          vert3tet[0], vert3tet[1], vert3tet[2]});
        trianglePos.push_back({
          vert31[0],   vert31[1],   vert31[2], 
          vert3t2[0],  vert3t2[1],  vert3t2[2], 
          vert3tet[0], vert3tet[1], vert3tet[2]});
        trianglePos.push_back({
          vert32[0],   vert32[1],   vert32[2], 
          vert3t2[0],  vert3t2[1],  vert3t2[2], 
          vert3tet[0], vert3tet[1], vert3tet[2]});
        trianglePos.push_back({
          vert32[0],   vert32[1],   vert32[2], 
          vert3t0[0],  vert3t0[1],  vert3t0[2], 
          vert3tet[0], vert3tet[1], vert3tet[2]});

        msLabels.insert(msLabels.end(), {
          msm[0], msm[0], msm[0], msm[0], msm[0], msm[0],
          msm[1], msm[1], msm[1], msm[1], msm[1], msm[1],
          msm[2], msm[2], msm[2], msm[2], msm[2], msm[2],
          msm[3], msm[3], msm[3], msm[3], msm[3], msm[3]});

        caseData.insert(caseData.end(), {
          lookupIndex, lookupIndex, lookupIndex,
          lookupIndex, lookupIndex, lookupIndex,
          lookupIndex, lookupIndex, lookupIndex,
          lookupIndex, lookupIndex, lookupIndex,
          lookupIndex, lookupIndex, lookupIndex,
          lookupIndex, lookupIndex, lookupIndex,
          lookupIndex, lookupIndex, lookupIndex,
          lookupIndex, lookupIndex, lookupIndex});
      }
    }
  }
//#else // TTK_ENABLE_OPENMP
} else {
  size_t triangleStartIndex[threadNumber_ + 1];
  #pragma omp parallel num_threads(threadNumber_)
  {
    const int tid = omp_get_thread_num();
    std::vector<std::array<float, 9>> trianglePosLocal;
    std::vector<SimplexId> caseDataLocal;
    std::vector<long long> msLabelsLocal;
    std::vector<SimplexId> sep2VertSetLocal;

    #pragma omp for schedule(guided,8) nowait
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

          float triCenter[4][3];
          getCenter(vPos[0], vPos[1], vPos[2], triCenter[0]);
          getCenter(vPos[0], vPos[1], vPos[3], triCenter[1]);
          getCenter(vPos[0], vPos[2], vPos[3], triCenter[2]);
          getCenter(vPos[1], vPos[2], vPos[3], triCenter[3]);

          float edgeCenters[10][3];

          getCenter(vPos[0], vPos[1], edgeCenters[0]);
          getCenter(vPos[0], vPos[2], edgeCenters[1]);
          getCenter(vPos[0], vPos[3], edgeCenters[2]);
          getCenter(vPos[1], vPos[2], edgeCenters[3]);
          getCenter(vPos[1], vPos[3], edgeCenters[4]);
          getCenter(vPos[2], vPos[3], edgeCenters[5]);


          float edge00[3], edge01[3], edge02[3], edge03[3], tri00[3], tri01[3],
                edge10[3], edge11[3], edge12[3], tri10[3], tri11[3],
                edge20[3], edge21[3], edge22[3], tri20[3], tri21[3];

          interpolatePoints(vPos[vIds[0]], edgeCenters[vIds[1]],  dc, edge00);
          interpolatePoints(vPos[vIds[0]], edgeCenters[vIds[2]],  dc, edge01);
          interpolatePoints(vPos[vIds[0]], triCenter[vIds[3]],    dc, tri00);
          interpolatePoints(vPos[vIds[4]], edgeCenters[vIds[5]],  dc, edge02);
          interpolatePoints(vPos[vIds[4]], edgeCenters[vIds[6]],  dc, edge03);
          interpolatePoints(vPos[vIds[4]], triCenter[vIds[7]],    dc, tri01);
 
          interpolatePoints(vPos[vIds[8]], edgeCenters[vIds[1]],  dc, edge10);
          interpolatePoints(vPos[vIds[8]], edgeCenters[vIds[5]],  dc, edge11);
          interpolatePoints(vPos[vIds[8]], edgeCenters[vIds[10]], dc, edge12);
          interpolatePoints(vPos[vIds[8]], triCenter[vIds[3]],    dc, tri10);
          interpolatePoints(vPos[vIds[8]], triCenter[vIds[7]],    dc, tri11);
 
          interpolatePoints(vPos[vIds[9]], edgeCenters[vIds[2]],  dc, edge20);
          interpolatePoints(vPos[vIds[9]], edgeCenters[vIds[6]],  dc, edge21);
          interpolatePoints(vPos[vIds[9]], edgeCenters[vIds[10]], dc, edge22);
          interpolatePoints(vPos[vIds[9]], triCenter[vIds[3]],    dc, tri20);
          interpolatePoints(vPos[vIds[9]], triCenter[vIds[7]],    dc, tri21);

          trianglePosLocal.push_back({
            edge00[0], edge00[1], edge00[2], 
            edge02[0], edge02[1], edge02[2], 
            tri00[0], tri00[1], tri00[2]});
          trianglePosLocal.push_back({
            edge02[0], edge02[1], edge02[2], 
            tri00[0], tri00[1], tri00[2], 
            tri01[0], tri01[1], tri01[2]});
          trianglePosLocal.push_back({
            edge01[0], edge01[1], edge01[2], 
            edge03[0], edge03[1], edge03[2], 
            tri00[0], tri00[1], tri00[2]});
          trianglePosLocal.push_back({
            edge03[0], edge03[1], edge03[2], 
            tri00[0], tri00[1], tri00[2], 
            tri01[0], tri01[1], tri01[2]});

          trianglePosLocal.push_back({
            edge10[0], edge10[1], edge10[2], 
            edge11[0], edge11[1], edge11[2], 
            tri10[0], tri10[1], tri10[2]});
          trianglePosLocal.push_back({
            edge11[0], edge11[1], edge11[2], 
            tri10[0], tri10[1], tri10[2], 
            tri11[0], tri11[1], tri11[2]});
          trianglePosLocal.push_back({
            edge12[0], edge12[1], edge12[2], 
            tri10[0], tri10[1], tri10[2], 
            tri11[0], tri11[1], tri11[2]});

          trianglePosLocal.push_back({
            edge20[0], edge20[1], edge20[2], 
            edge21[0], edge21[1], edge21[2], 
            tri20[0], tri20[1], tri20[2]});
          trianglePosLocal.push_back({
            edge21[0], edge21[1], edge21[2], 
            tri20[0], tri20[1], tri20[2], 
            tri21[0], tri21[1], tri21[2]});
          trianglePosLocal.push_back({
            edge22[0], edge22[1], edge22[2], 
            tri20[0], tri20[1], tri20[2], 
            tri21[0], tri21[1], tri21[2]});


          msLabelsLocal.insert(msLabelsLocal.end(), {
            msm[vIds[0]], msm[vIds[0]], msm[vIds[0]], msm[vIds[0]],
            msm[vIds[8]], msm[vIds[8]], msm[vIds[8]],
            msm[vIds[9]], msm[vIds[9]], msm[vIds[9]]});

          caseDataLocal.insert(caseDataLocal.end(), {
            lookupIndex, lookupIndex, lookupIndex,
            lookupIndex, lookupIndex, lookupIndex,
            lookupIndex, lookupIndex, lookupIndex, lookupIndex});

        } else { // 4 labels
          float tetCenter[3];
          getCenter(vPos[0], vPos[1], vPos[2], vPos[3], tetCenter);

          // the 4 triangle centers
          float triCenter[4][3];
          getCenter(vPos[0], vPos[1], vPos[2], triCenter[0]);
          getCenter(vPos[0], vPos[1], vPos[3], triCenter[1]);
          getCenter(vPos[0], vPos[2], vPos[3], triCenter[2]);
          getCenter(vPos[1], vPos[2], vPos[3], triCenter[3]);

          float vert00[3], vert01[3], vert02[3], vert0tet[3],
                vert0t0[3], vert0t1[3], vert0t2[3], 
                vert10[3], vert11[3], vert12[3], vert1tet[3],
                vert1t0[3], vert1t1[3], vert1t2[3], 
                vert20[3], vert21[3], vert22[3], vert2tet[3],
                vert2t0[3], vert2t1[3], vert2t2[3], 
                vert30[3], vert31[3], vert32[3], vert3tet[3],
                vert3t0[3], vert3t1[3], vert3t2[3]; 

          interpolatePoints(vPos[0], vPos[1], d0, vert00);
          interpolatePoints(vPos[0], vPos[2], d0, vert01);
          interpolatePoints(vPos[0], vPos[3], d0, vert02);
          interpolatePoints(vPos[0], tetCenter, dc, vert0tet);
          interpolatePoints(vPos[0], triCenter[0], dc, vert0t0);
          interpolatePoints(vPos[0], triCenter[1], dc, vert0t1);
          interpolatePoints(vPos[0], triCenter[2], dc, vert0t2);

          interpolatePoints(vPos[1], vPos[0], d0, vert10);
          interpolatePoints(vPos[1], vPos[2], d0, vert11);
          interpolatePoints(vPos[1], vPos[3], d0, vert12);
          interpolatePoints(vPos[1], tetCenter, dc, vert1tet);
          interpolatePoints(vPos[1], triCenter[0], dc, vert1t0);
          interpolatePoints(vPos[1], triCenter[1], dc, vert1t1);
          interpolatePoints(vPos[1], triCenter[3], dc, vert1t2);

          interpolatePoints(vPos[2], vPos[0], d0, vert20);
          interpolatePoints(vPos[2], vPos[1], d0, vert21);
          interpolatePoints(vPos[2], vPos[3], d0, vert22);
          interpolatePoints(vPos[2], tetCenter, dc, vert2tet);
          interpolatePoints(vPos[2], triCenter[0], dc, vert2t0);
          interpolatePoints(vPos[2], triCenter[2], dc, vert2t1);
          interpolatePoints(vPos[2], triCenter[3], dc, vert2t2);

          interpolatePoints(vPos[3], vPos[0], d0, vert30);
          interpolatePoints(vPos[3], vPos[1], d0, vert31);
          interpolatePoints(vPos[3], vPos[2], d0, vert32);
          interpolatePoints(vPos[3], tetCenter, dc, vert3tet);
          interpolatePoints(vPos[3], triCenter[1], dc, vert3t0);
          interpolatePoints(vPos[3], triCenter[2], dc, vert3t1);
          interpolatePoints(vPos[3], triCenter[3], dc, vert3t2);

          // Label Vert 0
          trianglePosLocal.push_back({
            vert00[0],   vert00[1],   vert00[2], 
            vert0t0[0],  vert0t0[1],  vert0t0[2], 
            vert0tet[0], vert0tet[1], vert0tet[2]});
          trianglePosLocal.push_back({
            vert00[0],   vert00[1],   vert00[2], 
            vert0t1[0],  vert0t1[1],  vert0t1[2], 
            vert0tet[0], vert0tet[1], vert0tet[2]});
          trianglePosLocal.push_back({
            vert01[0],   vert01[1],   vert01[2], 
            vert0t0[0],  vert0t0[1],  vert0t0[2], 
            vert0tet[0], vert0tet[1], vert0tet[2]});
          trianglePosLocal.push_back({
            vert01[0],   vert01[1],   vert01[2], 
            vert0t2[0],  vert0t2[1],  vert0t2[2], 
            vert0tet[0], vert0tet[1], vert0tet[2]});
          trianglePosLocal.push_back({
            vert02[0],   vert02[1],   vert02[2], 
            vert0t2[0],  vert0t2[1],  vert0t2[2], 
            vert0tet[0], vert0tet[1], vert0tet[2]});
          trianglePosLocal.push_back({
            vert02[0],   vert02[1],   vert02[2], 
            vert0t1[0],  vert0t1[1],  vert0t1[2], 
            vert0tet[0], vert0tet[1], vert0tet[2]});


          // Label Vert 1
          trianglePosLocal.push_back({
            vert10[0],   vert10[1],   vert10[2], 
            vert1t0[0],  vert1t0[1],  vert1t0[2], 
            vert1tet[0], vert1tet[1], vert1tet[2]});
          trianglePosLocal.push_back({
            vert10[0],   vert10[1],   vert10[2], 
            vert1t1[0],  vert1t1[1],  vert1t1[2], 
            vert1tet[0], vert1tet[1], vert1tet[2]});
          trianglePosLocal.push_back({
            vert11[0],   vert11[1],   vert11[2], 
            vert1t0[0],  vert1t0[1],  vert1t0[2], 
            vert1tet[0], vert1tet[1], vert1tet[2]});
          trianglePosLocal.push_back({
            vert11[0],   vert11[1],   vert11[2], 
            vert1t2[0],  vert1t2[1],  vert1t2[2], 
            vert1tet[0], vert1tet[1], vert1tet[2]});
          trianglePosLocal.push_back({
            vert12[0],   vert12[1],   vert12[2], 
            vert1t2[0],  vert1t2[1],  vert1t2[2], 
            vert1tet[0], vert1tet[1], vert1tet[2]});
          trianglePosLocal.push_back({
            vert12[0],   vert12[1],   vert12[2], 
            vert1t1[0],  vert1t1[1],  vert1t1[2], 
            vert1tet[0], vert1tet[1], vert1tet[2]});


          // Label Vert 2
          trianglePosLocal.push_back({
            vert20[0],   vert20[1],   vert20[2], 
            vert2t0[0],  vert2t0[1],  vert2t0[2], 
            vert2tet[0], vert2tet[1], vert2tet[2]});
          trianglePosLocal.push_back({
            vert20[0],   vert20[1],   vert20[2], 
            vert2t1[0],  vert2t1[1],  vert2t1[2], 
            vert2tet[0], vert2tet[1], vert2tet[2]});
          trianglePosLocal.push_back({
            vert21[0],   vert21[1],   vert21[2], 
            vert2t0[0],  vert2t0[1],  vert2t0[2], 
            vert2tet[0], vert2tet[1], vert2tet[2]});
          trianglePosLocal.push_back({
            vert21[0],   vert21[1],   vert21[2], 
            vert2t2[0],  vert2t2[1],  vert2t2[2], 
            vert2tet[0], vert2tet[1], vert2tet[2]});
          trianglePosLocal.push_back({
            vert22[0],   vert22[1],   vert22[2], 
            vert2t2[0],  vert2t2[1],  vert2t2[2], 
            vert2tet[0], vert2tet[1], vert2tet[2]});
          trianglePosLocal.push_back({
            vert22[0],   vert22[1],   vert22[2], 
            vert2t1[0],  vert2t1[1],  vert2t1[2], 
            vert2tet[0], vert2tet[1], vert2tet[2]});


          // Label Vert 3
          trianglePosLocal.push_back({
            vert30[0],   vert30[1],   vert30[2], 
            vert3t0[0],  vert3t0[1],  vert3t0[2], 
            vert3tet[0], vert3tet[1], vert3tet[2]});
          trianglePosLocal.push_back({
            vert30[0],   vert30[1],   vert30[2], 
            vert3t1[0],  vert3t1[1],  vert3t1[2], 
            vert3tet[0], vert3tet[1], vert3tet[2]});
          trianglePosLocal.push_back({
            vert31[0],   vert31[1],   vert31[2], 
            vert3t0[0],  vert3t0[1],  vert3t0[2], 
            vert3tet[0], vert3tet[1], vert3tet[2]});
          trianglePosLocal.push_back({
            vert31[0],   vert31[1],   vert31[2], 
            vert3t2[0],  vert3t2[1],  vert3t2[2], 
            vert3tet[0], vert3tet[1], vert3tet[2]});
          trianglePosLocal.push_back({
            vert32[0],   vert32[1],   vert32[2], 
            vert3t2[0],  vert3t2[1],  vert3t2[2], 
            vert3tet[0], vert3tet[1], vert3tet[2]});
          trianglePosLocal.push_back({
            vert32[0],   vert32[1],   vert32[2], 
            vert3t1[0],  vert3t1[1],  vert3t1[2], 
            vert3tet[0], vert3tet[1], vert3tet[2]});

          msLabelsLocal.insert(msLabelsLocal.end(), {
            msm[0], msm[0], msm[0], msm[0], msm[0], msm[0],
            msm[1], msm[1], msm[1], msm[1], msm[1], msm[1],
            msm[2], msm[2], msm[2], msm[2], msm[2], msm[2],
            msm[3], msm[3], msm[3], msm[3], msm[3], msm[3]});

          caseDataLocal.insert(caseDataLocal.end(), {
            lookupIndex, lookupIndex, lookupIndex,
            lookupIndex, lookupIndex, lookupIndex,
            lookupIndex, lookupIndex, lookupIndex,
            lookupIndex, lookupIndex, lookupIndex,
            lookupIndex, lookupIndex, lookupIndex,
            lookupIndex, lookupIndex, lookupIndex,
            lookupIndex, lookupIndex, lookupIndex,
            lookupIndex, lookupIndex, lookupIndex});
        }
      }
    }

    triangleStartIndex[tid + 1] = trianglePosLocal.size();

    #pragma omp barrier

    #pragma omp single
    {
      triangleStartIndex[0] = 0;

      // Count triangle number and create iterator start indices
      for(int t = 1; t < threadNumber_ + 1; ++t) {
        triangleStartIndex[t] += triangleStartIndex[t-1];
      }

      trianglePos.resize(triangleStartIndex[threadNumber_]);
    }

    std::copy(trianglePosLocal.begin(), trianglePosLocal.end(),
      trianglePos.begin() + triangleStartIndex[tid]);

    #pragma omp single
    {
      caseData.resize(triangleStartIndex[threadNumber_]);
    }

    std::copy(caseDataLocal.begin(), caseDataLocal.end(),
      caseData.begin() + triangleStartIndex[tid]);

    #pragma omp single
    {
      msLabels.resize(triangleStartIndex[threadNumber_]);
    }

    std::copy(msLabelsLocal.begin(), msLabelsLocal.end(),
      msLabels.begin() + triangleStartIndex[tid]);
  }
} // TTK_ENABLE_OPENMP
  this->printMsg("Computed 2-Separatrices 3D[Fine]",
                 1, localTimer.getElapsedTime(), this->threadNumber_);

  return 0;
}

template <typename triangulationType>
int ttk::MorseSmaleSegmentationPL::computeBasinSeparation_3D_fast(
  std::vector<std::array<float, 9>> &trianglePos,    
  std::vector<SimplexId> &caseData,
  std::vector<long long> &msLabels,
  const SimplexId *const morseSmaleManifold,
  const triangulationType &triangulation) const {
  ttk::Timer localTimer;

  this->printMsg("Computing 2-Separatrices 3D[Fast]",
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
      if(tetraederLookupFast[lookupIndex]) { // <= 2 labels on tetraeder
        float vertPos[4][3];
        triangulation.getVertexPoint(
          vertices[0], vertPos[0][0], vertPos[0][1], vertPos[0][2]);
        triangulation.getVertexPoint(
          vertices[1], vertPos[1][0], vertPos[1][1], vertPos[1][2]);
        triangulation.getVertexPoint(
          vertices[2], vertPos[2][0], vertPos[2][1], vertPos[2][2]);
        triangulation.getVertexPoint(
          vertices[3], vertPos[3][0], vertPos[3][1], vertPos[3][2]);

        switch(tetraederLookupFastCase[lookupIndex]) {
          case 0:
            trianglePos.push_back({
              vertPos[0][0], vertPos[0][1], vertPos[0][2], 
              vertPos[1][0], vertPos[1][1], vertPos[1][2], 
              vertPos[2][0], vertPos[2][1], vertPos[2][2]
            });
            msLabels.push_back(msm[0]);
            break;
          case 1:
            trianglePos.push_back({
              vertPos[0][0], vertPos[0][1], vertPos[0][2], 
              vertPos[1][0], vertPos[1][1], vertPos[1][2], 
              vertPos[3][0], vertPos[3][1], vertPos[3][2]
            });
            msLabels.push_back(msm[0]);
            break;
          case 2:
            trianglePos.push_back({
              vertPos[0][0], vertPos[0][1], vertPos[0][2], 
              vertPos[2][0], vertPos[2][1], vertPos[2][2], 
              vertPos[3][0], vertPos[3][1], vertPos[3][2]
            });
            msLabels.push_back(msm[0]);
            break;
          case 3:
            trianglePos.push_back({
              vertPos[1][0], vertPos[1][1], vertPos[1][2], 
              vertPos[2][0], vertPos[2][1], vertPos[2][2], 
              vertPos[3][0], vertPos[3][1], vertPos[3][2]
            });
            msLabels.push_back(msm[1]);
            break;
        }

        caseData.push_back(lookupIndex);
      }
    }
  }
//#else // TTK_ENABLE_OPENMP
} else {
  size_t triangleStartIndex[threadNumber_ + 1];
  #pragma omp parallel num_threads(threadNumber_)
  {
    const int tid = omp_get_thread_num();
    std::vector<std::array<float, 9>> trianglePosLocal;
    std::vector<SimplexId> caseDataLocal;
    std::vector<long long> msLabelsLocal;

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

      if(tetraederLookupFast[lookupIndex]) {
        const int id0 =
          tetraederLookupFastTri[tetraederLookupFastCase[lookupIndex]][0];
        const int id1 =
          tetraederLookupFastTri[tetraederLookupFastCase[lookupIndex]][1];
        const int id2 =
          tetraederLookupFastTri[tetraederLookupFastCase[lookupIndex]][2];
        
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
          vertPos[id0][0], vertPos[id0][1], vertPos[id0][2], 
          vertPos[id1][0], vertPos[id1][1], vertPos[id1][2], 
          vertPos[id2][0], vertPos[id2][1], vertPos[id2][2]
        });
        msLabelsLocal.push_back(msm[id0]);

        caseDataLocal.push_back(lookupIndex);
      }
    }

    triangleStartIndex[tid + 1] = trianglePosLocal.size();

    #pragma omp barrier

    #pragma omp single
    {
      triangleStartIndex[0] = 0;

      // Count triangle number and create iterator start indices
      for(int t = 1; t < threadNumber_ + 1; ++t) {
        triangleStartIndex[t] += triangleStartIndex[t-1];
      }

      trianglePos.resize(triangleStartIndex[threadNumber_]);
    }

    std::copy(trianglePosLocal.begin(), trianglePosLocal.end(),
      trianglePos.begin() + triangleStartIndex[tid]);

    #pragma omp single
    {
      caseData.resize(triangleStartIndex[threadNumber_]);
    }

    std::copy(caseDataLocal.begin(), caseDataLocal.end(),
      caseData.begin() + triangleStartIndex[tid]);

    #pragma omp single
    {
      msLabels.resize(triangleStartIndex[threadNumber_]);
    }

    std::copy(msLabelsLocal.begin(), msLabelsLocal.end(),
      msLabels.begin() + triangleStartIndex[tid]);
  }
} // TTK_ENABLE_OPENMP
  this->printMsg("Computed 2-Separatrices 3D[Fast]",
                 1, localTimer.getElapsedTime(), this->threadNumber_);

  return 0;
}

template <typename triangulationType>
int ttk::MorseSmaleSegmentationPL::computeSaddleExtremaConnectors_3D(
  std::vector<mscq::Separatrix> &sadExtrConns,
  const SimplexId *const ascNeighbor,
  const SimplexId *const dscNeighbor,
  const SimplexId *const ascManifold,
  const SimplexId *const dscManifold,
  const SimplexId *const orderArr,
  const std::vector<std::pair<SimplexId, char>> &criticalPoints,
  const triangulationType &triangulation) const {
    
  ttk::Timer localTimer;
  this->printMsg("Computing crit point connectors",
                 0, localTimer.getElapsedTime(), 
                 this->threadNumber_, ttk::debug::LineMode::REPLACE);

size_t sadExtrStartIndex[threadNumber_ + 1];
#pragma omp parallel num_threads(threadNumber_)
{
  const int tid = omp_get_thread_num();
  std::vector<mscq::Separatrix> sadExtrConnsLocal;

  std::vector<std::pair<SimplexId, SimplexId>> reachableExtrema[2];

  #pragma omp for schedule(static) nowait
  for(size_t c = 0; c < criticalPoints.size(); ++c) {
    const std::pair<SimplexId, char> &crit = criticalPoints[c];
    
    const bool isSaddle1 = crit.second == (char)CriticalType::Saddle1;
    const bool isSaddle2 = crit.second == (char)CriticalType::Saddle2;

    if(isSaddle1 || isSaddle2) {
      std::map<SimplexId, SimplexId> neigborSeedsAsc;
      std::map<SimplexId, SimplexId> neigborSeedsDsc;

      const SimplexId numNeighbors
        = triangulation.getVertexNeighborNumber(crit.first);

      for(SimplexId i = 0; i < numNeighbors; ++i) {
        SimplexId neighbor;
        triangulation.getVertexNeighbor(crit.first, i, neighbor);

        if(isSaddle1) {
          const SimplexId& ascManNeigh = ascManifold[neighbor];

          if(neigborSeedsAsc.find(ascManNeigh) == neigborSeedsAsc.end()) {
            neigborSeedsAsc.insert({ascManNeigh, neighbor});
          } else if(
            orderArr[neigborSeedsAsc[ascManNeigh]] > orderArr[neighbor]) {
            neigborSeedsAsc[ascManNeigh] = neighbor;
          }
        } else { // 2-saddle
          const SimplexId& dscManNeigh = dscManifold[neighbor];

          if(neigborSeedsDsc.find(dscManNeigh)  == neigborSeedsDsc.end()) {
            neigborSeedsDsc.insert({dscManNeigh, neighbor});
          } else
            if(orderArr[neigborSeedsDsc[dscManNeigh]] < orderArr[neighbor]) {
            neigborSeedsDsc[dscManNeigh] = neighbor;
          }
        }
      }

      for(const auto& s : neigborSeedsAsc) {
        reachableExtrema[0].push_back(std::make_pair(crit.first, s.second));
      }
      for(const auto& s : neigborSeedsDsc) {
        reachableExtrema[1].push_back(std::make_pair(crit.first, s.second));
      }
    }
  }

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

    sadExtrConnsLocal.push_back(sadExtrConn);
  }

  std::map<long long, size_t> s2ToMax;
  for(const auto& c : reachableExtrema[1]) { // 2-saddle maximum
    mscq::Separatrix sadExtrConn(c.first);
    sadExtrConn.type_ = 1;

    SimplexId currentVert = c.second;

    while(dscNeighbor[currentVert] != currentVert) {
      sadExtrConn.geometry_.push_back(currentVert);
      currentVert = dscNeighbor[currentVert];
    }

    sadExtrConn.geometry_.push_back(currentVert);

    sadExtrConnsLocal.push_back(sadExtrConn);
  }

  sadExtrStartIndex[tid + 1] = sadExtrConnsLocal.size();

  #pragma omp barrier

  #pragma omp single
  {
    sadExtrStartIndex[0] = 0;

    // Count triangle number and create iterator start indices
    for(int t = 1; t < threadNumber_ + 1; ++t) {
      sadExtrStartIndex[t] += sadExtrStartIndex[t-1];
    }

    sadExtrConns.resize(sadExtrStartIndex[threadNumber_]);
  }

  std::copy(sadExtrConnsLocal.begin(), sadExtrConnsLocal.end(),
    sadExtrConns.begin() + sadExtrStartIndex[tid]);
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

size_t sadExtrStartIndex[threadNumber_ + 1];
#pragma omp parallel num_threads(threadNumber_)
{
  const int tid = omp_get_thread_num();
  std::vector<mscq::Separatrix> sadExtrConnsLocal;

  #pragma omp for schedule(static) nowait
  for(size_t c = 0; c < criticalPoints.size(); ++c) {
    const std::pair<SimplexId, char> &crit = criticalPoints[c];
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
            sadExtrConn.type_ = 1;
          }

          sadExtrConn.geometry_.push_back(neighbor);
          sadExtrConnsLocal.push_back(sadExtrConn);

          neigborSet.insert(neighbor);
        }
      }
    }
  }

  sadExtrStartIndex[tid + 1] = sadExtrConnsLocal.size();

  #pragma omp barrier

  #pragma omp single
  {
    sadExtrStartIndex[0] = 0;

    // Count triangle number and create iterator start indices
    for(int t = 1; t < threadNumber_ + 1; ++t) {
      sadExtrStartIndex[t] += sadExtrStartIndex[t-1];
    }

    sadExtrConns.resize(sadExtrStartIndex[threadNumber_]);
  }

  std::copy(sadExtrConnsLocal.begin(), sadExtrConnsLocal.end(),
    sadExtrConns.begin() + sadExtrStartIndex[tid]);
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
  auto &extremaDistance = *outputseparatrices1_cells_extremaDistance_;
  outputseparatrices1_cells_extremaDistance_->resize(prevNcells + numCells);
  auto &extremaDistanceAbs = *outputseparatrices1_cells_extremaDistanceAbs_;
  outputseparatrices1_cells_extremaDistanceAbs_->resize(prevNcells + numCells);

  SimplexId currPId = prevNpoints, currCId = prevNcells;
  for(size_t i = 0; i < numSep; ++i) {
    const auto &sep = separatrices[i];

    const size_t geoSize = sep.geometry_.size() - 2 > 0 ? 
      sep.geometry_.size() - 2 : 1;

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

      const int jInv = geoSize - (j - 1);
      extremaDistance[currCId] = (float)(jInv) / (float)geoSize;
      extremaDistanceAbs[currCId] = jInv;

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

  const size_t numTriangles = trianglePos.size();

  size_t npoints = numTriangles * 3;
  size_t numCells = numTriangles;

  // resize arrays
  outputSeparatrices2_points_->resize(3 * npoints);
  auto &points = *outputSeparatrices2_points_;
  outputSeparatrices2_cells_connectivity_->resize(3 * numCells);
  auto &cellsConn = *outputSeparatrices2_cells_connectivity_;
  outputSeparatrices2_cells_cases_->resize(numTriangles);
  auto &cellsMSCIds = *outputSeparatrices2_cells_mscIds_;
  outputSeparatrices2_cells_mscIds_->resize(numTriangles);
  auto &cellsCase = *outputSeparatrices2_cells_cases_;
  cellsCase = caseData;

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_) schedule(static)
#endif
  for(size_t tri = 0; tri < numTriangles; ++tri) {
    const size_t ninetri = 9 * tri;
    const size_t threetri = 3 * tri;
    auto &triPos = trianglePos[tri];
    points[ninetri    ] = triPos[0];
    points[ninetri + 1] = triPos[1];
    points[ninetri + 2] = triPos[2];

    points[ninetri + 3] = triPos[3];
    points[ninetri + 4] = triPos[4];
    points[ninetri + 5] = triPos[5];

    points[ninetri + 6] = triPos[6];
    points[ninetri + 7] = triPos[7];
    points[ninetri + 8] = triPos[8];

    cellsConn[threetri    ] = threetri;
    cellsConn[threetri + 1] = threetri + 1;
    cellsConn[threetri + 2] = threetri + 2;

    cellsMSCIds[tri] = msLabelMap[msLabels[tri]];
  }

  // update pointers
  *outputSeparatrices2_numberOfPoints_ = npoints;
  *outputSeparatrices2_numberOfCells_ = numCells;

  this->printMsg("Wrote " + std::to_string(numTriangles) + " 2-separatrices",
                 1, localTimer.getElapsedTime(), this->threadNumber_);

  return 1;
}
