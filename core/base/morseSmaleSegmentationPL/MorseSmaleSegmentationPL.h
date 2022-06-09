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

/**
 * Lookup table to retrieve the triangle coordinates for each tetrahedron case.
 * The comments next to each line first indicate the number of unique vertex 
 * labels on the tetrahedron, the 2nd label, 3rd label, 4th label, and its index
 * (the first label is alway considered to be 0). Every 3 entries in each row
 * represent the locations on the tetrahedron to connect to a triangle, where 
 * 0-4 are edge centers, 5-8 triangle centers, and 9 is the tetraherdron center.
 * This allows to directly retrieve the triangulation from this lookup table
 */
const int tetLookupWall[28][15] = {
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

/**
 * Lookup table providing the labels of each triangle created by tetLookupWall.
 * Every two entries describe the two labels of both vertices that are separated
 * by each triangle, allowing to calculate a sparse Id for each triangle.
 */
const int tetLookupWallLabel[27][10] = {
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

/**
 * Lookup table providing number of triangles created by tetLookupWall.
 */
const size_t tetLookupNumWallTriangles[28] = {
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

/**
 * Lookup table providing the number of triangles created by tetLookupWall in
 * each case, using the "Fine" Mode that splits every triangle into two slightly
 * shifted triangles. This gives each basin its own geometry.
 */
const size_t tetLookupNumFineTriangles[28] = {
  {0},  // (1) 0,0,0 - 0
  {0},  // (-) 0,0,1 - 1
  {0},  // (-) 0,0,2 - 2 
  {2},  // (2) 0,0,3 - 3
  {0},  // (-) 0,1,0 - 4
  {0},  // (-) 0,1,1 - 5 
  {0},  // (-) 0,1,2 - 6
  {0},  // (-) 0,1,3 - 7
  {2},  // (2) 0,2,0 - 8
  {0},  // (-) 0,2,1 - 9
  {4},  // (2) 0,2,2 -10
  {10}, // (3) 0,2,3 -11
  {0},  // (-) 0,3,0 -12
  {0},  // (-) 0,3,1 -13
  {0},  // (-) 0,3,2 -14
  {0},  // (-) 0,3,3 -15
  {2},  // (2) 1,0,0 -16
  {4},  // (2) 1,0,1 -17
  {0},  // (-) 1,0,2 -18
  {10}, // (3) 1,0,3 -19
  {4},  // (2) 1,1,0 -20
  {2},  // (2) 1,1,1 -21
  {0},  // (-) 1,1,2 -22
  {10}, // (3) 1,1,3 -23
  {10}, // (3) 1,2,0 -24
  {10}, // (3) 1,2,1 -25
  {10}, // (3) 1,2,2 -26
  {24}  // (4) 1,2,3 -27
};

/**
 * Lookup table providing the vertex with the 2nd and 3rd(if available) label of 
 * a tetrahedron. They are given as the index of the vertex in the tetrahedron.
 * The first label always comes from the 0th vertex of the tetrahedron.
 */
const int tetraederUniqueLabels[27][2] = {
  {-1, -1}, // (1) 0,0,0 - 0
  {-1, -1}, // (-) 0,0,1 - 1
  {-1, -1}, // (-) 0,0,2 - 2 
  { 3, -1}, // (2) 0,0,3 - 3
  {-1, -1}, // (-) 0,1,0 - 4
  {-1, -1}, // (-) 0,1,1 - 5 
  {-1, -1}, // (-) 0,1,2 - 6
  {-1, -1}, // (-) 0,1,3 - 7
  { 2, -1}, // (2) 0,2,0 - 8
  {-1, -1}, // (-) 0,2,1 - 9
  { 2, -1}, // (2) 0,2,2 -10
  { 2,  3}, // (3) 0,2,3 -11
  {-1, -1}, // (-) 0,3,0 -12
  {-1, -1}, // (-) 0,3,1 -13
  {-1, -1}, // (-) 0,3,2 -14
  {-1, -1}, // (-) 0,3,3 -15
  { 1, -1}, // (2) 1,0,0 -16
  { 1, -1}, // (2) 1,0,1 -17
  {-1, -1}, // (-) 1,0,2 -18
  { 1,  3}, // (3) 1,0,3 -19
  { 1, -1}, // (2) 1,1,0 -20
  { 1, -1}, // (2) 1,1,1 -21
  {-1, -1}, // (-) 1,1,2 -22
  { 1,  3}, // (3) 1,1,3 -23
  { 1,  2}, // (3) 1,2,0 -24
  { 1,  2}, // (3) 1,2,1 -25
  { 1, -1}  // (3) 1,2,2 -26
};

/** Lookup table providing if a lookupIndex is a multi label case. */
const bool tetLookupIsMultiLabel[28] = {
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


/**  Lookup table providing if a lookupIndex is a two label case. */
const bool tetLookupIs2Label[28] = {
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

/**  Lookup table providing if a lookupIndex is a three label case. */
const bool tetLookupIs3Label[28] = {
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

/**
 * Lookup table providing the vertex index(local to the tetrahedron) to specify
 * the label for the triangles using the "Fast" Mode.
 */
const int tetLookupFastCase[28] = {
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

/** Does this tetrahedron have 3 vertices with the same label? */
const bool tetLookupFast[28] = {
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

/** 
 * Provide the vertex indices, local to the tetrahedron, for the cases in 
 * lookup table tetLookupFastCase. These three vertices have the same label.
 */
const int tetLookupFastTri[4][3] = {
  {0, 1, 2},
  {0, 1, 3},
  {0, 2, 3},
  {1, 2, 3}
};

/**
 * Vertex indices for interpolation in the "Fine" Mode where two labels are
 * present at one tetrahedron. Every pair of twovertex indices is used to
 * interpolate along their connecting edge. Therefore, three edges are given
 * when three vertices have the same label, as onetriangle is enough to separate
 * the remaining vertex. In the case of two vertices with the same label, two
 * triangles are needed, so 4 edge centers are given by the 8 vertex indices.
 */
const int tetLookupSplitBasins2Label[22][8] = {
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

/**
 * Vertex indices for interpolation in the "Fine" Mode where three labels are
 * present at one tetrahedron. In the case of three labels, two vertices have
 * the same label. For those two vertices the vertex index local to the
 * tetrahedron, the two edges connected to the vertices that are needed for
 * the triangulation, and the triangle index are given at positions 0-3 and 4-7
 * respectively. Position 8 and 9 give the vertex index of the two remaining
 * vertices and Position 10 gives the edge index conntecting them.
 */
const int tetLookupSplitBasisns3Label[27][11] = {
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

/** Retrieve all labels in the range [0,3] ignoring the index */
const int lookupOtherLabels[4][3] = {
  {1,2,3}, // 0
  {0,2,3}, // 1
  {0,1,3}, // 2
  {0,1,2}, // 3
};

/** Multi label triangles cases. */
const bool triangleLookupIsMultiLabel[7] = {
  false, // (1) 0,0 - 0
  false, // (-) 0,1 - 1
  true,  // (2) 0,2 - 2 
  false, // (-) 0,3 - 3
  true,  // (2) 1,0 - 4
  true,  // (2) 1,1 - 5 
  true   // (3) 1,2 - 6
};

/** Two label triangle cases. */
const bool triangleLookupIs2Label[7] = {
  false, // (1) 0,0 - 0
  false, // (-) 0,1 - 1
  true,  // (2) 0,2 - 2 
  false, // (-) 0,3 - 3
  true,  // (2) 1,0 - 4
  true,  // (2) 1,1 - 5 
  false  // (3) 1,2 - 6
};

/** 
 * Two label triangle cases. Every 2 entries represent the edge idices that
 * are used to create the splitting edge.
 */
const int triangleLookupEdgeVerts[7][4] = {
  {-1, -1, -1, -1}, // (1) 0,0 - 0
  {-1, -1, -1, -1}, // (-) 0,1 - 1
  { 1,  2,  0,  2}, // (2) 0,2 - 2 
  {-1, -1, -1, -1}, // (-) 0,3 - 3
  { 0,  1,  1,  2}, // (2) 1,0 - 4
  { 0,  1,  0,  2}, // (2) 1,1 - 5
};

/** Two label triangle cases. */
const int triangleLookupEdgeCount[7] = {
  0, // (1) 0,0 - 0
  0, // (-) 0,1 - 1
  1,  // (2) 0,2 - 2 
  0, // (-) 0,3 - 3
  1,  // (2) 1,0 - 4
  1,  // (2) 1,1 - 5 
  3  // (3) 1,2 - 6
};


namespace ttk {

  namespace mscq {

    /**
     * Utility class representing Piecewise Linear Integral Lines.
     */
    struct IntegralLine {


      explicit IntegralLine()
        : geometry_{}, type_{-1} {
      }

      /** Initialization with one segment. */
      explicit IntegralLine(const SimplexId segmentGeometry) {
        geometry_.push_back(segmentGeometry);
      }

      /** Assignment by moving the content to the new Integral Line */
      IntegralLine &operator=(const IntegralLine &separatrix) noexcept {
        geometry_ = std::move(separatrix.geometry_);
        type_ = separatrix.type_;

        return *this;
      }

      /** Container of vertexIds of all points of the Integral Line. */
      std::vector<SimplexId> geometry_;

      /** Type of IntegralLine, e.g. ascending or descending. */
      int type_{-1};
    };
  } // namespace mscq

  /**
   * The MorseSmaleSegmentationPL class provides methods to compute a
   * MSC variant that only uses the input field and its order field.
   */
  class MorseSmaleSegmentationPL : virtual public Debug {

  public:
    MorseSmaleSegmentationPL();

    /** @brief Type of segmentation */
    enum class SEPARATRICES_MANIFOLD {
      MORSESMALE = 0,
      ASCENDING = 1,
      DESCENDING = 2
    };
    
    /** @brief Type of 1-separatrix output */
    enum class SEPARATRICES1_MODE {
      NONE = 0,
      DETAILED = 1,
      SIMPLE = 2
    };

    /** @brief Type of 2-separatrix output */
    enum class SEPARATRICES2_MODE {
      NONE = 0,
      WALLS = 1,
      SEPARATEBASINSFINE = 2,
      SEPARATEBASINSFAST = 3,
      EXPERIMENT = 4,
    };    

    /**
     * Set the input triangulation and preprocess the needed
     * mesh traversal queries.
     * 
     * @param calculateSaddles Should Saddles be calculated?
     * @param triangulation Triangulation
     */
    int preconditionTriangulation(
      const bool calculateSaddles,
      ttk::AbstractTriangulation *triangulation) const {
      int success = 0;
      success += triangulation->preconditionVertexNeighbors();
      
      if(calculateSaddles) {
        success += triangulation->preconditionBoundaryVertices();
        success += triangulation->preconditionVertexStars();
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
     * Set the data pointers to the output segmentation fields.
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
     * Set all necessary critical point output variables
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

    /**
     * Set all necessary 1-separatrix output variables
     */
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

    /**
     * Set all necessary 2-separatrix output variables
     */
    inline int setOutputSeparatrices2(
      SimplexId *const separatrices2_numberOfPoints,
      std::vector<float> *const separatrices2_points,
      SimplexId *const separatrices2_numberOfCells,
      std::vector<SimplexId> *const outputSeparatrices2_cells_connectivity,
      std::vector<SimplexId> *const outputSeparatrices2_cells_mscIds) {
      outputSeparatrices2_numberOfPoints_ = separatrices2_numberOfPoints;
      outputSeparatrices2_points_ = separatrices2_points;
      outputSeparatrices2_numberOfCells_ = separatrices2_numberOfCells;
      outputSeparatrices2_cells_connectivity_ =
        outputSeparatrices2_cells_connectivity;
      outputSeparatrices2_cells_mscIds_ = outputSeparatrices2_cells_mscIds;
      return 0;
    }

    /**
     * Execute the package.
     * @param[in] sepManifoldType Type of segmentation
     * @param[in] sep1Mode Type of 1-separatrices
     * @param[in] sep2Mode Type of 2-separatrices
     * @param[in] computeSaddles Should Saddles be calculated?
     * @param[in] triangulation Triangulation
     * @return Returns 0 upon success, negative values otherwise.
     *
     * @pre For this function to behave correctly setInputOrderField
     * needs to be called to fill the inputOrderField_ buffer prior to any
     * computation. MorseSmaleSegmentationPL::preconditionTriangulation must be 
     * called before execution to enable Vertex based operations. 
     * setOutputMorseComplexes, setOutputCriticalPoints, setOutputSeparatrices1,
     * and setOutputSeparatrices2 need to be called to set all necessary ouput
     * buffers.
     * @see examples/c++/main.cpp for an example use.
     */
    template <typename triangulationType>
    int execute(
      const SEPARATRICES_MANIFOLD sepManifoldType,
      const SEPARATRICES1_MODE sep1Mode,
      const SEPARATRICES2_MODE sep2Mode, 
      const bool computeSaddles,
      const triangulationType &triangulation);

    /**
     * @brief Compute the ascending or descending segmentation.
     * 
     * @param[out] criticalPoints List of critical Points with Id and Type
     * @param[out] manifold Scalar field to save segementation to
     * @param[out] neighbor Scalar field saving the largest/smallest neighbor Id
     * @param[out] critMap Map from critical point index to vertex index.
     * @param[out] numberOfExtrema number of extrema after computation
     * @param[in] orderArr Input order array
     * @param[in] triangulation Triangulation
     * @param[in] ascending Should the Ascending segmentation be computed
     * @return int 0 on success.
     */
    template <typename triangulationType>
    int computeSingleSegmentation(
      std::vector<std::pair<SimplexId, char>> *criticalPoints,
      SimplexId *const manifold,
      SimplexId *const neighbor,
      std::map<SimplexId, SimplexId> &critMap,
      SimplexId &numberOfExtrema,
      const SimplexId *const orderArr,
      const triangulationType &triangulation,
      const bool ascending) const;

    /**
     * @brief Compute the ascending and descending segmentation.
     * 
     * @param[out] criticalPoints List of critical Points with Id and Type
     * @param[out] ascManifold Scalar field to write the ascending segementation to
     * @param[out] dscManifold Scalar field to write the descending segementation to
     * @param[out] ascNeighbor Scalar field to write the smallest neighbor to
     * @param[out] dscNeighbor Scalar field to write the largest neighbor to
     * @param[out] ascCritMap Map from critical point index to vertex index
     * @param[out] dscCritMap Map from critical point index to vertex index
     * @param[out] numberOfMinima number of minima after computation
     * @param[out] numberOfMaxima number of maxima after computation
     * @param[in] orderArr Input order array
     * @param[in] triangulation Triangulation
     * @return int 0 on success. 
     */
    template <typename triangulationType>
    int computeAscDscSegmentations(
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

    /**
     * @brief Intersect ascending and descending Segmentation
     * 
     * @param[out] numManifolds 
     * @param[out] morseSmaleManifold 
     * @param[in] numMaxima maxima count
     * @param[in] ascendingManifold Ascending segmentation scalar field
     * @param[in] descendingManifold Descending segmentation scalar field
     * @param[in] triangulation 
     * @return int 0 on success.
     */
    template <typename triangulationType>
    int computeFinalSegmentation(
      SimplexId &numManifolds,
      SimplexId *const morseSmaleManifold,
      const SimplexId numMaxima,
      const SimplexId *const ascendingManifold,
      const SimplexId *const descendingManifold,
      const triangulationType &triangulation) const;

    /**
     * @brief Retrieve 1-separatrices in 2D domains
     * 
     * @param[out] sep2VertSet Set of all vertices belonging to 3-label tris
     * @param[in] numMSCRegions Count of regions given by the segmentation
     * @param[in] morseSmaleManifold segmentation of the domain
     * @param[in] computeSaddles Should saddles be computed
     * @param[in] triangulation Triangulation
     * @return int 0 on success.
     */
    template <typename triangulationType>
    int computeSeparatrices1Pieces_2D(
      std::vector<SimplexId> &sep2VertSet,
      const SimplexId &numMSCRegions,
      const SimplexId *const morseSmaleManifold,
      const bool computeSaddles,
      const triangulationType &triangulation) const;

    /**
     * @brief Check the Saddle criterion on a list of vertices
     * 
     * @param[out] criticalPoints List of critical points
     * @param[in] sep2VertSet List of Vertex Ids
     * @param[in] triangulation Triangulation
     * @return int 0 on success.
     */
    template <typename triangulationType>
    int findSaddlesFromCadidates(
      std::vector<std::pair<SimplexId, char>> &criticalPoints,
      const std::vector<SimplexId> &sep2VertSet,
      const triangulationType &triangulation) const;

    /**
     * @brief Compute all critical points
     * 
     * @param[out] criticalPoints Critical points list
     * @param[in] orderArr orderArray
     * @param[in] triangulation Triangulation
     * @return int 0 on success
     */
    int getCriticalPoints(
      std::vector<std::pair<SimplexId, char>> &criticalPoints,
      const SimplexId * const offsets,
      const AbstractTriangulation &triangulation) const;

    /**
     * @brief Compute Separating Walls [Wall Mode]
     * 
     * @param[in] numMSCRegions MS regions count
     * @param[in] morseSmaleManifold MS segmentation
     * @param[in] triangulation Triangulation
     * @return int 0 on success
     */
    template <typename triangulationType>
    int computeSeparatrices2_3D_wall(
      const SimplexId &numMSCRegions,
      const SimplexId *const morseSmaleManifold,
      const triangulationType &triangulation) const;

    /**
     * @brief Compute Separating Walls for Asc and Desc separately [Wall Mode]
     * 
     * @param[in] numAscRegions Ascending regions count
     * @param[in] numDscRegions Descending regions count
     * @param[in] ascManifold Ascending segmentation
     * @param[in] dscManifold Descending segmentation
     * @param[in] triangulation Triangulation
     * @return int 0 on success
     */
    template <typename triangulationType>
    int computeSeparatrices2_3D_split(
      const SimplexId &numAscRegions,
      const SimplexId &numDscRegions,
      const SimplexId *const ascManifold,
      const SimplexId *const dscManifold,
      const triangulationType &triangulation) const;

    /**
     * @brief Compute Separating Walls [Basin Separation Mode]
     * 
     * @param[in] morseSmaleManifold MS segmentation
     * @param[in] triangulation Triangulation
     * @return int 0 on success
     */
    template <typename triangulationType>
    int computeBasinSeparation_3D_fine(
      const SimplexId *const morseSmaleManifold,
      const triangulationType &triangulation) const;

    /**
     * @brief Compute Separating Walls [Fast Basin Separation Mode]
     * 
     * @param[in] morseSmaleManifold MS Segmentation
     * @param[in] triangulation Triangulation
     * @return int 0 on success
     */
    template <typename triangulationType>
    int computeBasinSeparation_3D_fast(
      const SimplexId *const morseSmaleManifold,
      const triangulationType &triangulation) const;

    /**
     * @brief Computes Extrema saddle connectors as PL integral lines
     * 
     * @param[out] sadExtrConns Resulting Saddle Extrema connectors
     * @param[in] ascendingNeighbor Scalar field representing smallest neighbors
     * @param[in] descendingNeighbor Scalar field representing largest neighbors
     * @param[in] ascendingManifold Ascending segmentation
     * @param[in] descendingManifold Descensing sementation
     * @param[in] orderArr Order Array
     * @param[in] criticalPoints Critical points list
     * @param[in] triangulation Triangulation
     * @return int 0 on success
     */
    template <typename triangulationType>
    int computeSaddleExtremaConnectors_3D(
      std::vector<mscq::IntegralLine> &sadExtrConns,
      const SimplexId *const ascendingNeighbor,
      const SimplexId *const descendingNeighbor,
      const SimplexId *const ascendingManifold,
      const SimplexId *const descendingManifold,
      const SimplexId * const orderArr,
      const std::vector<std::pair<SimplexId, char>> &criticalPoints,
      const triangulationType &triangulation) const;

    /**
     * @brief Computes Extrema saddle connectors as straight lines
     * 
     * @param[out] sadExtrConns Resulting Saddle Extrema connectors
     * @param[in] ascendingManifold Ascending sementation
     * @param[in] descendingManifold Descending segmentation
     * @param[in] ascCritMap Map: Minima index to its vertex Id
     * @param[in] dscCritMap Map: Maximum index to its vertex Id
     * @param[in] criticalPoints Critical points list
     * @param[in] triangulation Triangulation
     * @return int 0 on success
     */
    template <typename triangulationType>
    int computeSaddleExtremaConnectors_3D_fast(
      std::vector<mscq::IntegralLine> &sadExtrConns,
      const SimplexId *const ascendingManifold,
      const SimplexId *const descendingManifold,
      const std::map<SimplexId, SimplexId> &ascCritMap,
      const std::map<SimplexId, SimplexId> &dscCritMap,
      const std::vector<std::pair<SimplexId, char>> &criticalPoints,
      const triangulationType &triangulation) const;

    /**
     * @brief Write the 1-separatrices to memory
     * 
     * @param[in] separatrices 1-separatrix list
     * @param[in] triangulation Triangulation
     * @return int 0 on success
     */
    template <typename triangulationType>
    int setSeparatrices1_3D(const std::vector<mscq::IntegralLine> &separatrices,
                            const triangulationType &triangulation) const;

    /**
     * @brief Write the critical points to memory
     * 
     * @param[in] criticalPoints Critical points list
     * @param[in] triangulation Triangulation
     * @return int 0 on success
     */
    template <typename triangulationType>
    int setCriticalPoints(std::vector<std::pair<SimplexId, char>> &criticalPoints,
                          const triangulationType &triangulation) const;

    /**
     * Get a consistent Sparse Id of three indicies. Used to retrieve a
     * unique Id for a combination of three IDs in the same range [0, numIds-1]
     * that is as small as possible and consistent by considering the smallest
     * ID first, and the largest second.
     * 
     * @param simID0 First Id
     * @param simID1 Second Id
     * @param numIds Number of available Ids
     * @return consistent sparse ID
     */
    inline long long getSparseId(
      const long long simID0, const long long simID1,
      const SimplexId &numMSCRegions) const {
      if(simID0 < simID1) {
        return simID0 * numMSCRegions + simID1;
      }
      return simID1 * numMSCRegions + simID0;
    }

    /**
     * Get a consistent Sparse Id of three indicies. Used to retrieve a
     * unique Id for a combination of three IDs in the same range [0, numIds-1]
     * that is as small as possible and consistent by considering the smallest
     * ID first, the second smallest second and the largest third.
     * 
     * @param simID0 First Id
     * @param simID1 Second Id
     * @param simID2 Third Id
     * @param numIds Number of available Ids
     * @return consistent sparse ID
     */
    inline long long getSparseId(const long long simID0,
                                 const long long simID1,
                                 const long long simID2,
                                 const SimplexId &numIds) const {
      if(simID0 < simID1) {
        if(simID1 < simID2) {
          return simID0 * numIds * numIds + simID1 * numIds + simID2;
        } else {
          if(simID0 < simID2) {
            return simID0 * numIds * numIds + simID2 * numIds + simID1;
          } else {
            return simID2 * numIds * numIds + simID0 * numIds + simID1;
          }
        }
      } else {
        if(simID0 < simID2) {
          return simID1 * numIds * numIds + simID0 * numIds + simID2;
        } else {
          if(simID1 < simID2) {
            return simID1 * numIds * numIds + simID2 * numIds + simID0;
          } else {
            return simID2 * numIds * numIds + simID1 * numIds + simID0;
          }
        }
      }
    }

    /**
     * @brief Get the center of two points
     * 
     * @param[in] pos0 Position 0
     * @param[in] pos1 Position 1
     * @param[out] incenter Resulting position
     * @return int 0 on success
     */
    inline int getCenter(
      float pos0[3], float pos1[3], float incenter[3]) const {
      incenter[0] = 0.5 * (pos0[0] + pos1[0]);
      incenter[1] = 0.5 * (pos0[1] + pos1[1]);
      incenter[2] = 0.5 * (pos0[2] + pos1[2]);

      return 0;
    }

    /**
     * @brief Get the center of three points
     * 
     * @param[in] pos0 Position 0
     * @param[in] pos1 Position 1
     * @param[in] pos2 Position 2
     * @param[out] incenter Resulting position
     * @return int 0 on success
     */
    inline int getCenter(
      float pos0[3], float pos1[3], float pos2[3], float incenter[3]) const {
      incenter[0] = 0.3333 * (pos0[0] + pos1[0] + pos2[0]);
      incenter[1] = 0.3333 * (pos0[1] + pos1[1] + pos2[1]);
      incenter[2] = 0.3333 * (pos0[2] + pos1[2] + pos2[2]);

      return 0;
    }

    /**
     * @brief Get the center of four points
     * 
     * @param[in] pos0 Position 0
     * @param[in] pos1 Position 1
     * @param[in] pos2 Position 2
     * @param[in] pos3 Position 3
     * @param[out] incenter Resulting position
     * @return int 0 on success
     */
    inline int getCenter(
      float pos0[3], float pos1[3], float pos2[3], float pos3[3],
      float incenter[3]) const {
      incenter[0] = 0.25 * (pos0[0] + pos1[0] + pos2[0] + pos3[0]);
      incenter[1] = 0.25 * (pos0[1] + pos1[1] + pos2[1] + pos3[1]);
      incenter[2] = 0.25 * (pos0[2] + pos1[2] + pos2[2] + pos3[2]);

      return 0;
    }

    /**
     * @brief Interpolate between two points (lambda = 0 -> pos1 / 1 -> pos0)
     * 
     * @param[in] pos0 Position 0
     * @param[in] pos1 Position 1
     * @param[in] lambda Interpolation parameter
     * @param[out] result Resulting position
     * @return int 0 on success
     */
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
      std::vector<SimplexId> *outputSeparatrices2_cells_mscIds_{};
      std::vector<SimplexId> *outputSeparatrices2_cells_connectivity_{};
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
#endif // TTK_ENABLE_OPENMP

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
    getCriticalPoints(criticalPoints, orderArr, triangulation);
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
        computeAscDscSegmentations(
          critPointPtr, ascendingManifold, descendingManifold,
          ascendingNeighbor, descendingNeighbor, ascCritMap, dscCritMap,
          numberOfMinima, numberOfMaxima, orderArr, triangulation);

        computeFinalSegmentation(numberOfMSCRegions, morseSmaleManifold,
          numberOfMinima, ascendingManifold, descendingManifold, triangulation);

        sepManifold = static_cast<SimplexId *>(outputMorseSmaleManifold_);
        SepManifoldIsValid = true;
      }
      break;
    case SEPARATRICES_MANIFOLD::ASCENDING :
      if(ascendingManifold) {
        computeSingleSegmentation<triangulationType>(
          critPointPtr, ascendingManifold, ascendingNeighbor, ascCritMap,
          numberOfMinima, orderArr, triangulation, true);
        
        sepManifold = static_cast<SimplexId *>(outputAscendingManifold_);
        numberOfMSCRegions = numberOfMinima;
        SepManifoldIsValid = true;
      }
      break;
    case SEPARATRICES_MANIFOLD::DESCENDING :
      if(descendingManifold) {
        computeSingleSegmentation<triangulationType>(
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
    std::vector<mscq::IntegralLine> sadExtrConns;
    
    if(SepManifoldIsValid) {
      if(sep2Mode == SEPARATRICES2_MODE::WALLS) {
        computeSeparatrices2_3D_wall<triangulationType>(
          numberOfMSCRegions, sepManifold, triangulation);
      } else if(sep2Mode == SEPARATRICES2_MODE::SEPARATEBASINSFINE) {
        computeBasinSeparation_3D_fine<triangulationType>(
          sepManifold, triangulation);
      } else if(sep2Mode == SEPARATRICES2_MODE::SEPARATEBASINSFAST) {
        computeBasinSeparation_3D_fast<triangulationType>(
          sepManifold, triangulation);
      } else if(sep2Mode == SEPARATRICES2_MODE::EXPERIMENT) {
        computeSeparatrices2_3D_split<triangulationType>(
          numberOfMinima, numberOfMinima,
          ascendingManifold, descendingManifold, triangulation);
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
  } else if(dim == 2) {
    std::vector<std::array<float, 6>> edgePos;
    std::vector<long long> msLabels;
    std::vector<SimplexId> sep2VertSet;
    std::map<long long, SimplexId> msLabelMap;

    computeSeparatrices1Pieces_2D(
      sep2VertSet, numberOfMSCRegions, 
      sepManifold, computeSaddles, triangulation);

    if(computeSaddles) {
      findSaddlesFromCadidates(criticalPoints, sep2VertSet, triangulation);
    }
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
int ttk::MorseSmaleSegmentationPL::computeSingleSegmentation(
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

  numberOfExtrema = 0;

  const bool computeCriticalPoints = criticalPoints != nullptr;
  const bool computeNeighbor = (neighbor != nullptr);

  // -----------------------------------------------------------------------
  // Compute MSC variant
  // -----------------------------------------------------------------------
  {
    /* compute the Descending Maifold iterating over each vertex, searching
     * the biggest neighbor and compressing its path to its maximum */
    const SimplexId nVertices = triangulation.getNumberOfVertices();
    std::map<SimplexId, SimplexId> extremas;

#ifdef TTK_ENABLE_OPENMP
if(threadNumber_ == 1) {
#endif // TTK_ENABLE_OPENMP

    // vertices that may still be compressed
    std::vector<SimplexId>* activeVertices = new std::vector<SimplexId>();
    SimplexId nActiveVertices;

    // find maxima and intialize vector of not fully compressed vertices
    for(SimplexId i = 0; i < nVertices; i++) {
      SimplexId neighborId;
      const SimplexId numNeighbors =
        triangulation.getVertexNeighborNumber(i);
      
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
      } else {
        if(ascending) {
          if(computeCriticalPoints) {
            criticalPoints->push_back(
              std::make_pair(i, (char)(CriticalType::Local_minimum)));
          }
        } else {
          if(computeCriticalPoints) {
            criticalPoints->push_back(
              std::make_pair(i, (char)(CriticalType::Local_maximum)));
          }
        }

        critMap.insert({numberOfExtrema, i});
        extremas.insert({i, numberOfExtrema++});
      }

      if(computeNeighbor) {
        neighbor[i] = dmi;
      }
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

      delete activeVertices;
      activeVertices = newActiveVert;

      nActiveVertices = activeVertices->size();
    }

    // make indices dense
    for(SimplexId i = 0; i < nVertices; i++) {
      manifold[i] = extremas[manifold[i]];
    }

    delete activeVertices;

#ifdef TTK_ENABLE_OPENMP
} else {

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
            newActiveVert->push_back(v);
          }
        }

        delete lActiveVertices;
        lActiveVertices = newActiveVert;
        lnActiveVertices = lActiveVertices->size();

        newActiveVert = new std::vector<SimplexId>();
      }

      delete newActiveVert;

      #pragma omp barrier

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
}
#endif // TTK_ENABLE_OPENMP

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
int ttk::MorseSmaleSegmentationPL::computeAscDscSegmentations(
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

  std::map<SimplexId, int> minMap, maxMap;

#ifdef TTK_ENABLE_OPENMP
if(threadNumber_ == 1) {
#endif // TTK_ENABLE_OPENMP

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

#ifdef TTK_ENABLE_OPENMP
} else {

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

}
#endif // TTK_ENABLE_OPENMP

  this->printMsg("Computed Manifolds", 1.0,
                 localTimer.getElapsedTime(), this->threadNumber_);

  return 1; // return success
}

template <typename triangulationType>
int ttk::MorseSmaleSegmentationPL::computeFinalSegmentation(
  SimplexId &numManifolds,
  SimplexId *const morseSmaleManifold,
  const SimplexId numberOfMinima,
  const SimplexId *const ascendingManifold,
  const SimplexId *const descendingManifold,
  const triangulationType &triangulation) const {
  ttk::Timer localTimer;

  this->printMsg("Computing MSC Manifold",
                 0, localTimer.getElapsedTime(), this->threadNumber_,
                 ttk::debug::LineMode::REPLACE);
  
  const size_t nVerts = triangulation.getNumberOfVertices();
  size_t numSparseIDs = 0;
  std::vector<SimplexId> sparseIds;
  std::map<SimplexId, size_t> sparseToDenseRegionId{};
  

#ifdef TTK_ENABLE_OPENMP
if(threadNumber_ == 1) {
#endif // TTK_ENABLE_OPENMP

  std::unordered_set<SimplexId> idSet;
  for(size_t i = 0; i < nVerts; ++i) {
    morseSmaleManifold[i] = ascendingManifold[i] * numberOfMinima +
      descendingManifold[i];
    idSet.insert(morseSmaleManifold[i]);
  }

  numSparseIDs = idSet.size();
  sparseIds.reserve(numSparseIDs);
  std::copy(idSet.begin(), idSet.end(), std::back_inserter(sparseIds));
  TTK_PSORT(this->threadNumber_, sparseIds.begin(), sparseIds.end());

  // "sparse region id" -> "dense region id"
  for(size_t i = 0; i < numSparseIDs; ++i) {
    sparseToDenseRegionId[sparseIds[i]] = i;
  }

  for(size_t i = 0; i < nVerts; ++i) {
    morseSmaleManifold[i] = sparseToDenseRegionId[morseSmaleManifold[i]];
  }

#ifdef TTK_ENABLE_OPENMP
} else {
  std::unordered_set<SimplexId>* lIdSet[threadNumber_];

  // parallel conquer unordered_sets variables
  int numConquer = threadNumber_ / 2;
  int conquerStep = 1;
  int cS2pow = 1;
  int unmergedSet = -1;
  bool isConquerEven = !(bool)(threadNumber_ % 2);

  #pragma omp parallel num_threads(threadNumber_)
  {
    const int tid = omp_get_thread_num();
    lIdSet[tid] = new std::unordered_set<SimplexId>;
    std::unordered_set<SimplexId>* localSet = lIdSet[tid];

    // Create sparseIds
    #pragma omp for schedule(static)
    for(size_t i = 0; i < nVerts; ++i) {
      morseSmaleManifold[i] = ascendingManifold[i] * numberOfMinima +
        descendingManifold[i];
      localSet->insert(morseSmaleManifold[i]);
    }

    // parallel conquer unordered_sets
    while(numConquer > 0) {
      #pragma omp for schedule(static) nowait
      for(int i = 0; i < numConquer * 2; i += 2) {
        lIdSet[i * cS2pow]->insert(lIdSet[i * cS2pow + cS2pow]->begin(), lIdSet[i * cS2pow + cS2pow]->end());
      }

      #pragma omp single
      { // Fix border problems
        if(!isConquerEven) {
          if(unmergedSet == -1) {
            unmergedSet = 2 * numConquer * cS2pow;
          } else {
            int unmergedSet2 = 2 * numConquer * cS2pow;
            lIdSet[unmergedSet2]->insert(lIdSet[unmergedSet]->begin(), lIdSet[unmergedSet]->end());
            unmergedSet = unmergedSet2;
          }
        }
      }

      #pragma omp single
      {
        conquerStep += 1;
        cS2pow = std::pow(2, conquerStep - 1);
        isConquerEven = !(bool)(numConquer % 2);
        numConquer = numConquer / 2;
      }
    }

    // Merge unmerged Set with set 0, sort ids and create a map
    #pragma omp single
    {
      if(unmergedSet != -1) {
        lIdSet[0]->insert(lIdSet[unmergedSet]->begin(), lIdSet[unmergedSet]->end());
      }

      numSparseIDs = lIdSet[0]->size();
      sparseIds.reserve(numSparseIDs);
      std::copy(lIdSet[0]->begin(), lIdSet[0]->end(), std::back_inserter(sparseIds));
      TTK_PSORT(this->threadNumber_, sparseIds.begin(), sparseIds.end());

      // "sparse region id" -> "dense region id"
      for(size_t i = 0; i < numSparseIDs; ++i) {
        sparseToDenseRegionId[sparseIds[i]] = i;
      }
    }

    delete(localSet);

    #pragma omp for schedule(static)
    for(size_t i = 0; i < nVerts; ++i) {
      morseSmaleManifold[i] = sparseToDenseRegionId[morseSmaleManifold[i]];
    }
  }
}
#endif // TTK_ENABLE_OPENMP

  numManifolds = numSparseIDs;

  this->printMsg("Computed MSC Manifold",
                 1, localTimer.getElapsedTime(), this->threadNumber_);

  return 0;
}

template <typename triangulationType>
int ttk::MorseSmaleSegmentationPL::computeSeparatrices1Pieces_2D(
  std::vector<SimplexId> &sep2VertSet,
  const SimplexId &numMSCRegions,
  const SimplexId *const morseSmaleManifold,
  const bool computeSaddles,
  const triangulationType &triangulation) const {
  // start a local timer for this subprocedure
  ttk::Timer localTimer;

  this->printMsg("Computing 2-separatrices",
                 0, // progress form 0-1
                 0, // elapsed time so far
                 this->threadNumber_, ttk::debug::LineMode::REPLACE);

  const SimplexId numTriangles = triangulation.getNumberOfTriangles();
  size_t numSparseIDs = 0;
  std::vector<long long> sparseIdVect;
  std::map<long long, SimplexId> msLabelMap;

#ifdef TTK_ENABLE_OPENMP
if(threadNumber_ == 1) {
#endif // TTK_ENABLE_OPENMP

  std::vector<std::pair<SimplexId, unsigned char>> validCases;
  std::unordered_set<long long> msLabelSet;
  std::unordered_set<SimplexId> vertexSet;
  SimplexId numEdges = 0;

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
      validCases.push_back(std::make_pair(tri, lookupIndex));
      numEdges += triangleLookupEdgeCount[lookupIndex];

      if(triangleLookupIs2Label[lookupIndex]) {

        msLabelSet.insert(getSparseId(
          msm[triangleLookupEdgeVerts[lookupIndex][0]],
          msm[triangleLookupEdgeVerts[lookupIndex][1]], numMSCRegions));

        if(computeSaddles) {
          if(triangulation.isVertexOnBoundary(vertices[0]))
            vertexSet.insert(vertices[0]);
          if(triangulation.isVertexOnBoundary(vertices[1]))
            vertexSet.insert(vertices[1]);
          if(triangulation.isVertexOnBoundary(vertices[2]))
            vertexSet.insert(vertices[2]);
        }
      } else {
        long long sparseMSIds[3] = {
          getSparseId(msm[0], msm[1], numMSCRegions),
          getSparseId(msm[0], msm[2], numMSCRegions),
          getSparseId(msm[1], msm[2], numMSCRegions)
        };

        msLabelSet.insert({sparseMSIds[0], sparseMSIds[1], sparseMSIds[2]});

        if(computeSaddles) {
          vertexSet.insert({vertices[0], vertices[1], vertices[2]});
        }
      }
    }
  }

  numSparseIDs = msLabelSet.size();
  sparseIdVect.reserve(numSparseIDs);
  std::copy(msLabelSet.begin(), msLabelSet.end(),
    std::back_inserter(sparseIdVect));
  std::copy(vertexSet.begin(), vertexSet.end(),
    std::back_inserter(sep2VertSet));
  TTK_PSORT(this->threadNumber_, sparseIdVect.begin(), sparseIdVect.end());

  // "sparse id" -> "dense id"
  for(size_t i = 0; i < numSparseIDs; ++i) {
    msLabelMap[sparseIdVect[i]] = i;
  }

  outputSeparatrices2_points_->resize(6 * numEdges);
  outputSeparatrices2_cells_connectivity_->resize(2 * numEdges);
  outputSeparatrices2_cells_mscIds_->resize(numEdges);
  *outputSeparatrices2_numberOfPoints_ = 2 * numEdges;
  *outputSeparatrices2_numberOfCells_ = numEdges;

  auto &points = *outputSeparatrices2_points_;
  auto &cellsConn = *outputSeparatrices2_cells_connectivity_;
  auto &cellsMSCIds = *outputSeparatrices2_cells_mscIds_;

  float* p = points.data();
  SimplexId* c = cellsConn.data();
  SimplexId* m = cellsMSCIds.data();
  SimplexId cellIndex = 0;

  for (const auto& vCase : validCases)
  {
    const SimplexId& tri = vCase.first;
    const unsigned char& lookupIndex = vCase.second;

    const int * edgeVerts = triangleLookupEdgeVerts[lookupIndex];

    SimplexId vertices[3];
    triangulation.getCellVertex(tri, 0, vertices[0]);
    triangulation.getCellVertex(tri, 1, vertices[1]);
    triangulation.getCellVertex(tri, 2, vertices[2]);

    float vertPos[3][3];
    triangulation.getVertexPoint(
      vertices[0], vertPos[0][0], vertPos[0][1], vertPos[0][2]);
    triangulation.getVertexPoint(
      vertices[1], vertPos[1][0], vertPos[1][1], vertPos[1][2]);
    triangulation.getVertexPoint(
      vertices[2], vertPos[2][0], vertPos[2][1], vertPos[2][2]);

    const SimplexId msm[3] = {
      morseSmaleManifold[vertices[0]], morseSmaleManifold[vertices[1]],
      morseSmaleManifold[vertices[2]]};

    if(triangleLookupIs2Label[lookupIndex]) {
      float eC[2][3];
      getCenter(vertPos[edgeVerts[0]], vertPos[edgeVerts[1]], eC[0]);
      getCenter(vertPos[edgeVerts[2]], vertPos[edgeVerts[3]], eC[1]);
      
      const long long sparseID = msLabelMap[getSparseId(
        msm[edgeVerts[0]], msm[edgeVerts[1]], numMSCRegions)];

      p[0] = eC[0][0]; p[1] = eC[0][1]; p[2] = eC[0][2]; 
      p[3] = eC[1][0]; p[4] = eC[1][1]; p[5] = eC[1][2];
      p += 6;

      c[0] = cellIndex; c[1] = cellIndex + 1;
      c += 2;
      cellIndex += 2;

      m[0] = sparseID;
      m += 1;

      if(computeSaddles) {
        if(triangulation.isVertexOnBoundary(vertices[0]))
          sep2VertSet.push_back(vertices[0]);
        if(triangulation.isVertexOnBoundary(vertices[1]))
          sep2VertSet.push_back(vertices[1]);
        if(triangulation.isVertexOnBoundary(vertices[2]))
          sep2VertSet.push_back(vertices[2]);
      }
    } else {
      float eC[4][3];
      getCenter(vertPos[0], vertPos[1], eC[0]);
      getCenter(vertPos[0], vertPos[2], eC[1]);
      getCenter(vertPos[1], vertPos[2], eC[2]);
      getCenter(vertPos[0], vertPos[1], vertPos[2], eC[3]);

      const long long sparseID[3] = {
          msLabelMap[getSparseId(msm[0], msm[1], numMSCRegions)],
          msLabelMap[getSparseId(msm[0], msm[2], numMSCRegions)],
          msLabelMap[getSparseId(msm[1], msm[2], numMSCRegions)]};

      p[0] = eC[0][0]; p[1] = eC[0][1]; p[2] = eC[0][2]; 
      p[3] = eC[3][0]; p[4] = eC[3][1]; p[5] = eC[3][2];
      p[6] = eC[1][0]; p[7] = eC[1][1]; p[8] = eC[1][2]; 
      p[9] = eC[3][0]; p[10] = eC[3][1]; p[11] = eC[3][2];
      p[12] = eC[2][0]; p[13] = eC[2][1]; p[14] = eC[2][2]; 
      p[15] = eC[3][0]; p[16] = eC[3][1]; p[17] = eC[3][2];
      p += 18;

      c[0] = cellIndex + 0; c[1] = cellIndex + 1;
      c[2] = cellIndex + 2; c[3] = cellIndex + 3;
      c[4] = cellIndex + 4; c[5] = cellIndex + 5;
      c += 6;
      cellIndex += 6;

      m[0] = sparseID[0];
      m[1] = sparseID[1];
      m[2] = sparseID[2];
      m += 3;
    }
  }

#ifdef TTK_ENABLE_OPENMP
} else {

  size_t edgeStartIndex[threadNumber_ + 1];

  // parallel conquer unordered_sets variables
  std::unordered_set<long long>* msSets[threadNumber_];
  std::unordered_set<SimplexId>* vertSets[threadNumber_];
  int numConquer = threadNumber_ / 2;
  int conquerStep = 1;
  int cS2pow = 1;
  int unmergedSet = -1;
  bool isConquerEven = !(bool)(threadNumber_ % 2);

  #pragma omp parallel num_threads(threadNumber_)
  {
    const int tid = omp_get_thread_num();
    std::vector<std::pair<SimplexId, unsigned char>> validCases;
    size_t numThreadEdges = 0;
    size_t numThreadIndexEdges = 0;

    msSets[tid] = new std::unordered_set<long long>;
    std::unordered_set<long long>* localMSSet = msSets[tid];

    vertSets[tid] = new std::unordered_set<SimplexId>;
    std::unordered_set<SimplexId>* localVertSet = vertSets[tid];

    #pragma omp single nowait
    this->printMsg("Computing 2-separatrices",
                 0.1, // progress form 0-1
                 0, // elapsed time so far
                 this->threadNumber_, ttk::debug::LineMode::REPLACE);

    #pragma omp for schedule(static) nowait
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

        validCases.push_back(std::make_pair(tri, lookupIndex));
        numThreadEdges += triangleLookupEdgeCount[lookupIndex];

        if(triangleLookupIs2Label[lookupIndex]) {
          localMSSet->insert(getSparseId(
          msm[triangleLookupEdgeVerts[lookupIndex][0]],
          msm[triangleLookupEdgeVerts[lookupIndex][1]], numMSCRegions));

          if(computeSaddles) {
            if(triangulation.isVertexOnBoundary(vertices[0]))
              localVertSet->insert(vertices[0]);
            if(triangulation.isVertexOnBoundary(vertices[1]))
              localVertSet->insert(vertices[1]);
            if(triangulation.isVertexOnBoundary(vertices[2]))
              localVertSet->insert(vertices[2]);
          }
        } else {
          long long sparseMSIds[3] = {
            getSparseId(msm[0], msm[1], numMSCRegions),
            getSparseId(msm[0], msm[2], numMSCRegions),
            getSparseId(msm[1], msm[2], numMSCRegions)
          };

          localMSSet->insert({sparseMSIds[0], sparseMSIds[1], sparseMSIds[2]});

          if(computeSaddles) {
            localVertSet->insert({vertices[0], vertices[1], vertices[2]});
          }
        }
      }
    }

    #pragma omp single nowait
    this->printMsg("Computing 2-separatrices",
                 0.2, // progress form 0-1
                 0, // elapsed time so far
                 this->threadNumber_, ttk::debug::LineMode::REPLACE);

    /************ CONQUER SETS ************/

    while(numConquer > 0) {
      #pragma omp for schedule(static) nowait
      for(int i = 0; i < numConquer * 2; i += 2) {
        msSets[i * cS2pow]->insert(
          msSets[i * cS2pow + cS2pow]->begin(),
          msSets[i * cS2pow + cS2pow]->end());
        vertSets[i * cS2pow]->insert(
          vertSets[i * cS2pow + cS2pow]->begin(),
          vertSets[i * cS2pow + cS2pow]->end());
      }

      #pragma omp single
      { // Fix border problems
        if(!isConquerEven) {
          if(unmergedSet == -1) {
            unmergedSet = 2 * numConquer * cS2pow;
          } else {
            int unmergedSet2 = 2 * numConquer * cS2pow;
            msSets[unmergedSet2]->insert(
              msSets[unmergedSet]->begin(), msSets[unmergedSet]->end());
            vertSets[unmergedSet2]->insert(
              vertSets[unmergedSet]->begin(), vertSets[unmergedSet]->end());
            unmergedSet = unmergedSet2;
          }
        }
      }

      #pragma omp single
      {
        conquerStep += 1;
        cS2pow = std::pow(2, conquerStep - 1);
        isConquerEven = !(bool)(numConquer % 2);
        numConquer = numConquer / 2;
      }
    }

    // Merge unmerged set with set 0, sort ids and create map
    #pragma omp single
    {
      if(unmergedSet != -1) {
        msSets[0]->insert(
          msSets[unmergedSet]->begin(), msSets[unmergedSet]->end());
        vertSets[0]->insert(
          vertSets[unmergedSet]->begin(), vertSets[unmergedSet]->end());
      }

      numSparseIDs = msSets[0]->size();
      sparseIdVect.reserve(numSparseIDs);
      std::copy(msSets[0]->begin(), msSets[0]->end(),
        std::back_inserter(sparseIdVect));
      std::copy(vertSets[0]->begin(), vertSets[0]->end(),
        std::back_inserter(sep2VertSet));
      TTK_PSORT(this->threadNumber_, sparseIdVect.begin(), sparseIdVect.end());

      // "sparse id" -> "dense id"
      for(size_t i = 0; i < numSparseIDs; ++i) {
        msLabelMap[sparseIdVect[i]] = i;
      }
    }

    delete(localMSSet);

    /************ END CONQUER SETS ************/

    #pragma omp single nowait
    this->printMsg("Computing 2-separatrices",
                 0.4, // progress form 0-1
                 0, // elapsed time so far
                 this->threadNumber_, ttk::debug::LineMode::REPLACE);

    edgeStartIndex[tid + 1] = numThreadEdges;

    #pragma omp barrier

    #pragma omp single
    {
      edgeStartIndex[0] = 0;

      // Count triangle number and create iterator start indices
      for(int t = 1; t < threadNumber_ + 1; ++t) {
        edgeStartIndex[t] += edgeStartIndex[t-1];
      }
    }

    const size_t numEdges = edgeStartIndex[threadNumber_];
    numThreadIndexEdges = edgeStartIndex[tid];

    #pragma omp single nowait
    outputSeparatrices2_points_->resize(6 * numEdges);

    #pragma omp single nowait
    outputSeparatrices2_cells_connectivity_->resize(2 * numEdges);

    #pragma omp single nowait
    outputSeparatrices2_cells_mscIds_->resize(numEdges);

    #pragma omp single
    {
      *outputSeparatrices2_numberOfPoints_ = 2 * numEdges;
      *outputSeparatrices2_numberOfCells_ = numEdges;
    }
  
    auto &points = *outputSeparatrices2_points_;
    auto &cellsConn = *outputSeparatrices2_cells_connectivity_;
    auto &cellsMSCIds = *outputSeparatrices2_cells_mscIds_;

    float* p = points.data();
    SimplexId* c = cellsConn.data();
    SimplexId* m = cellsMSCIds.data();

    p += (numThreadIndexEdges * 6);
    c += (numThreadIndexEdges * 2);
    m += numThreadIndexEdges;

    numThreadIndexEdges = 2 * numThreadIndexEdges;

    #pragma omp single nowait
    this->printMsg("Computing 2-separatrices",
                 0.5, // progress form 0-1
                 0, // elapsed time so far
                 this->threadNumber_, ttk::debug::LineMode::REPLACE);

    
    for (const auto& vCase : validCases)
    {
      const SimplexId& tri = vCase.first;
      const unsigned char& lookupIndex = vCase.second;

      const int * edgeVerts = triangleLookupEdgeVerts[lookupIndex];

      SimplexId vertices[3];
      triangulation.getCellVertex(tri, 0, vertices[0]);
      triangulation.getCellVertex(tri, 1, vertices[1]);
      triangulation.getCellVertex(tri, 2, vertices[2]);

      float vertPos[3][3];
      triangulation.getVertexPoint(
        vertices[0], vertPos[0][0], vertPos[0][1], vertPos[0][2]);
      triangulation.getVertexPoint(
        vertices[1], vertPos[1][0], vertPos[1][1], vertPos[1][2]);
      triangulation.getVertexPoint(
        vertices[2], vertPos[2][0], vertPos[2][1], vertPos[2][2]);

      const SimplexId msm[3] = {
        morseSmaleManifold[vertices[0]], morseSmaleManifold[vertices[1]],
        morseSmaleManifold[vertices[2]]};

      if(triangleLookupIs2Label[lookupIndex]) {
        float eC[2][3];
        getCenter(vertPos[edgeVerts[0]], vertPos[edgeVerts[1]], eC[0]);
        getCenter(vertPos[edgeVerts[2]], vertPos[edgeVerts[3]], eC[1]);
        
        const long long sparseID = msLabelMap[getSparseId(
          msm[edgeVerts[0]], msm[edgeVerts[1]], numMSCRegions)];

        p[0] = eC[0][0]; p[1] = eC[0][1]; p[2] = eC[0][2]; 
        p[3] = eC[1][0]; p[4] = eC[1][1]; p[5] = eC[1][2];
        p += 6;

        c[0] = numThreadIndexEdges; c[1] = numThreadIndexEdges + 1;
        c += 2;
        numThreadIndexEdges += 2;

        m[0] = sparseID;
        m += 1;

        if(computeSaddles) {
          if(triangulation.isVertexOnBoundary(vertices[0]))
            sep2VertSet.push_back(vertices[0]);
          if(triangulation.isVertexOnBoundary(vertices[1]))
            sep2VertSet.push_back(vertices[1]);
          if(triangulation.isVertexOnBoundary(vertices[2]))
            sep2VertSet.push_back(vertices[2]);
        }
      } else {
        float eC[4][3];
        getCenter(vertPos[0], vertPos[1], eC[0]);
        getCenter(vertPos[0], vertPos[2], eC[1]);
        getCenter(vertPos[1], vertPos[2], eC[2]);
        getCenter(vertPos[0], vertPos[1], vertPos[2], eC[3]);

        const long long sparseID[3] = {
          msLabelMap[getSparseId(msm[0], msm[1], numMSCRegions)],
          msLabelMap[getSparseId(msm[0], msm[2], numMSCRegions)],
          msLabelMap[getSparseId(msm[1], msm[2], numMSCRegions)]};

        p[0] = eC[0][0]; p[1] = eC[0][1]; p[2] = eC[0][2]; 
        p[3] = eC[3][0]; p[4] = eC[3][1]; p[5] = eC[3][2];
        p[6] = eC[1][0]; p[7] = eC[1][1]; p[8] = eC[1][2]; 
        p[9] = eC[3][0]; p[10] = eC[3][1]; p[11] = eC[3][2];
        p[12] = eC[2][0]; p[13] = eC[2][1]; p[14] = eC[2][2]; 
        p[15] = eC[3][0]; p[16] = eC[3][1]; p[17] = eC[3][2];
        p += 18;

        c[0] = numThreadIndexEdges + 0; c[1] = numThreadIndexEdges + 1;
        c[2] = numThreadIndexEdges + 2; c[3] = numThreadIndexEdges + 3;
        c[4] = numThreadIndexEdges + 4; c[5] = numThreadIndexEdges + 5;
        c += 6;
        numThreadIndexEdges += 6;

        m[0] = sparseID[0];
        m[1] = sparseID[1];
        m[2] = sparseID[2];
        m += 3;
      }
    }
  }
}
#endif // TTK_ENABLE_OPENMP

  this->printMsg("Computed 1-separatrix pieces", 1, localTimer.getElapsedTime(),
                 this->threadNumber_);

  return 1;
}

template <typename triangulationType>
int ttk::MorseSmaleSegmentationPL::findSaddlesFromCadidates(
  std::vector<std::pair<SimplexId, char>> &criticalPoints,
  const std::vector<SimplexId> &sep2VertSet,
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

#ifdef TTK_ENABLE_OPENMP  
if(threadNumber_ == 1) {
#endif // TTK_ENABLE_OPENMP

  for(const SimplexId& candidate: sep2VertSet) {
    const CriticalType critC = (CriticalType)sfcp.getCriticalType(
      candidate, orderArr, (AbstractTriangulation*)&triangulation);

    if(critC == CriticalType::Saddle1) {
      criticalPoints.push_back(std::make_pair(candidate, 1));
    } else if(critC == CriticalType::Saddle2) {
      criticalPoints.push_back(std::make_pair(candidate, 2));
    }
  }

#ifdef TTK_ENABLE_OPENMP
} else {
  #pragma omp parallel num_threads(threadNumber_)
  {
    #pragma omp for schedule(static)
    for(const SimplexId& candidate: sep2VertSet)
    {
      const CriticalType critC = (CriticalType)sfcp.getCriticalType(
        candidate, orderArr, (AbstractTriangulation*)&triangulation);
      
      if(critC == CriticalType::Saddle1 || critC == CriticalType::Saddle2) {
        #pragma omp critical
        {
          if(critC == CriticalType::Saddle1) {
            criticalPoints.push_back(std::make_pair(candidate, 1));
         } else {
           criticalPoints.push_back(std::make_pair(candidate, 2));
          }
        }
      }
    }
  }
}
#endif // TTK_ENABLE_OPENMP

  this->printMsg("Computed saddles", 1, localTimer.getElapsedTime(),
    this->threadNumber_);

  return 1;
}

int ttk::MorseSmaleSegmentationPL::getCriticalPoints(
  std::vector<std::pair<SimplexId, char>> &criticalPoints,
  const SimplexId * const orderArr,
  const AbstractTriangulation &triangulation) const {

  ttk::Timer localTimer;

  this->printMsg("Computing Critical Points", 0, localTimer.getElapsedTime(),
                 this->threadNumber_, ttk::debug::LineMode::REPLACE);

  const int dim = triangulation.getDimensionality();

  ScalarFieldCriticalPoints sfcp;
  sfcp.setDomainDimension(dim);
  sfcp.setOutput(&criticalPoints);
  sfcp.setThreadNumber(this->threadNumber_);
  sfcp.executeProgressive(orderArr, &triangulation);

  this->printMsg("Computed Critical Points", 1, localTimer.getElapsedTime(),
                 this->threadNumber_);

  return 1;
  
}

template <typename triangulationType>
int ttk::MorseSmaleSegmentationPL::computeSeparatrices2_3D_wall(
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
  size_t numSparseIDs = 0;
  std::vector<long long> sparseIdVect;
  std::map<long long, SimplexId> msLabelMap;

#ifdef TTK_ENABLE_OPENMP
if(threadNumber_ == 1) {
#endif // TTK_ENABLE_OPENMP

  std::vector<std::pair<SimplexId, unsigned char>> validCases;
  std::unordered_set<long long> msLabelSet;
  SimplexId numTriangles = 0;

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
    if(tetLookupIsMultiLabel[lookupIndex]) {
      validCases.push_back(std::make_pair(tet, lookupIndex));
      numTriangles += tetLookupNumWallTriangles[lookupIndex];

      const int *tetEdgeIndices = tetLookupWall[lookupIndex];
      if(tetEdgeIndices[0] == 10) { // 4 labels on tetraeder
        long long sparseMSIds[6] = {
          getSparseId(msm[0], msm[1], numMSCRegions),
          getSparseId(msm[0], msm[2], numMSCRegions),
          getSparseId(msm[0], msm[3], numMSCRegions),
          getSparseId(msm[1], msm[2], numMSCRegions),
          getSparseId(msm[1], msm[3], numMSCRegions),
          getSparseId(msm[2], msm[3], numMSCRegions)
        };

        msLabelSet.insert({
          sparseMSIds[0], sparseMSIds[1], sparseMSIds[2], 
          sparseMSIds[3], sparseMSIds[4], sparseMSIds[5]});
      } else {
        const int * const unqLabels = tetraederUniqueLabels[lookupIndex];

        if(tetLookupIs2Label[lookupIndex]) {          
          msLabelSet.insert(getSparseId(
            msm[0], msm[unqLabels[0]], numMSCRegions));

        } else {
          long long sparseMSIds[3] = {
            getSparseId(msm[0], msm[unqLabels[0]], numMSCRegions),
            getSparseId(msm[0], msm[unqLabels[1]], numMSCRegions),
            getSparseId(msm[unqLabels[0]], msm[unqLabels[1]], numMSCRegions)};

          msLabelSet.insert({
            sparseMSIds[0], sparseMSIds[1], sparseMSIds[2]});
        }
      }
    }
  }

  numSparseIDs = msLabelSet.size();
  sparseIdVect.reserve(numSparseIDs);
  std::copy(msLabelSet.begin(), msLabelSet.end(),
    std::back_inserter(sparseIdVect));
  TTK_PSORT(this->threadNumber_, sparseIdVect.begin(), sparseIdVect.end());

  // "sparse id" -> "dense id"
  for(size_t i = 0; i < numSparseIDs; ++i) {
    msLabelMap[sparseIdVect[i]] = i;
  }

  outputSeparatrices2_points_->resize(9 * numTriangles);
  outputSeparatrices2_cells_connectivity_->resize(3 * numTriangles);
  outputSeparatrices2_cells_mscIds_->resize(numTriangles);
  *outputSeparatrices2_numberOfPoints_ = 3 * numTriangles;
  *outputSeparatrices2_numberOfCells_ = numTriangles;

  auto &points = *outputSeparatrices2_points_;
  auto &cellsConn = *outputSeparatrices2_cells_connectivity_;
  auto &cellsMSCIds = *outputSeparatrices2_cells_mscIds_;

  float* p = points.data();
  SimplexId* c = cellsConn.data();
  SimplexId* m = cellsMSCIds.data();
  SimplexId cellIndex = 0;

  for (const auto& vCase : validCases)
  {
    const SimplexId& tet = vCase.first;
    const unsigned char& lookupIndex = vCase.second;

    const int *tetEdgeIndices = tetLookupWall[lookupIndex];
    const int *tetVertLabel = tetLookupWallLabel[lookupIndex];

    SimplexId vertices[4];
    triangulation.getCellVertex(tet, 0, vertices[0]);
    triangulation.getCellVertex(tet, 1, vertices[1]);
    triangulation.getCellVertex(tet, 2, vertices[2]);
    triangulation.getCellVertex(tet, 3, vertices[3]);

    const SimplexId msm[4] = {
      morseSmaleManifold[vertices[0]], morseSmaleManifold[vertices[1]],
      morseSmaleManifold[vertices[2]], morseSmaleManifold[vertices[3]]};

    float vertPos[4][3];
    triangulation.getVertexPoint(
      vertices[0], vertPos[0][0], vertPos[0][1], vertPos[0][2]);
    triangulation.getVertexPoint(
      vertices[1], vertPos[1][0], vertPos[1][1], vertPos[1][2]);
    triangulation.getVertexPoint(
      vertices[2], vertPos[2][0], vertPos[2][1], vertPos[2][2]);
    triangulation.getVertexPoint(
      vertices[3], vertPos[3][0], vertPos[3][1], vertPos[3][2]);

    float eC[10][3];
    // 6 edge centers
    getCenter(vertPos[0], vertPos[1], eC[0]);
    getCenter(vertPos[0], vertPos[2], eC[1]);
    getCenter(vertPos[0], vertPos[3], eC[2]);
    getCenter(vertPos[1], vertPos[2], eC[3]);
    getCenter(vertPos[1], vertPos[3], eC[4]);
    getCenter(vertPos[2], vertPos[3], eC[5]);

    // 4 triangle centers
    getCenter(vertPos[0], vertPos[1], vertPos[2], eC[6]);
    getCenter(vertPos[0], vertPos[1], vertPos[3], eC[7]);
    getCenter(vertPos[0], vertPos[2], vertPos[3], eC[8]);
    getCenter(vertPos[1], vertPos[2], vertPos[3], eC[9]);

    if(tetEdgeIndices[0] == 10) { // 4 labels on tetraeder
      float tetCenter[3];
      getCenter(vertPos[0], vertPos[1], vertPos[2], vertPos[3], tetCenter);

      long long sparseMSIds[6] = {
        msLabelMap[getSparseId(msm[0], msm[1], numMSCRegions)],
        msLabelMap[getSparseId(msm[0], msm[2], numMSCRegions)],
        msLabelMap[getSparseId(msm[0], msm[3], numMSCRegions)],
        msLabelMap[getSparseId(msm[1], msm[2], numMSCRegions)],
        msLabelMap[getSparseId(msm[1], msm[3], numMSCRegions)],
        msLabelMap[getSparseId(msm[2], msm[3], numMSCRegions)]
      };

      p[0] = eC[7][0]; p[1] = eC[7][1]; p[2] = eC[7][2]; 
      p[3] = eC[0][0]; p[4] = eC[0][1]; p[5] = eC[0][2]; 
      p[6] = tetCenter[0]; p[7] = tetCenter[1]; p[8] = tetCenter[2];
      p[9] = eC[0][0]; p[10] = eC[0][1]; p[11] = eC[0][2]; 
      p[12] = eC[6][0]; p[13] = eC[6][1]; p[14] = eC[6][2]; 
      p[15] = tetCenter[0]; p[16] = tetCenter[1]; p[17] = tetCenter[2];
      p[18] = eC[8][0]; p[19] = eC[8][1]; p[20] = eC[8][2];
      p[21] = eC[1][0]; p[22] = eC[1][1]; p[23] = eC[1][2];
      p[24] = tetCenter[0]; p[25] = tetCenter[1]; p[26] = tetCenter[2];
      p[27] = eC[1][0]; p[28] = eC[1][1]; p[29] = eC[1][2];
      p[30] = eC[6][0]; p[31] = eC[6][1]; p[32] = eC[6][2];
      p[33] = tetCenter[0]; p[34] = tetCenter[1]; p[35] = tetCenter[2];
      p[36] = eC[8][0]; p[37] = eC[8][1]; p[38] = eC[8][2]; 
      p[39] = eC[2][0]; p[40] = eC[2][1]; p[41] = eC[2][2]; 
      p[42] = tetCenter[0]; p[43] = tetCenter[1]; p[44] = tetCenter[2];
      p[45] = eC[2][0]; p[46] = eC[2][1]; p[47] = eC[2][2]; 
      p[48] = eC[7][0]; p[49] = eC[7][1]; p[50] = eC[7][2]; 
      p[51] = tetCenter[0]; p[52] = tetCenter[1]; p[53] = tetCenter[2];
      p[54] = eC[6][0]; p[55] = eC[6][1]; p[56] = eC[6][2];
      p[57] = eC[3][0]; p[58] = eC[3][1]; p[59] = eC[3][2];
      p[60] = tetCenter[0]; p[61] = tetCenter[1]; p[62] = tetCenter[2];
      p[63] = eC[3][0]; p[64] = eC[3][1]; p[65] = eC[3][2];
      p[66] = eC[9][0]; p[67] = eC[9][1]; p[68] = eC[9][2];
      p[69] = tetCenter[0]; p[70] = tetCenter[1]; p[71] = tetCenter[2];
      p[72] = eC[7][0]; p[73] = eC[7][1]; p[74] = eC[7][2];
      p[75] = eC[4][0]; p[76] = eC[4][1]; p[77] = eC[4][2];
      p[78] = tetCenter[0]; p[79] = tetCenter[1]; p[80] = tetCenter[2];
      p[81] = eC[4][0]; p[82] = eC[4][1]; p[83] = eC[4][2];
      p[84] = eC[9][0]; p[85] = eC[9][1]; p[86] = eC[9][2];
      p[87] = tetCenter[0]; p[88] = tetCenter[1]; p[89] = tetCenter[2];
      p[90] = eC[9][0]; p[91] = eC[9][1]; p[92] = eC[9][2]; 
      p[93] = eC[5][0]; p[94] = eC[5][1]; p[95] = eC[5][2]; 
      p[96] = tetCenter[0]; p[97] = tetCenter[1]; p[98] = tetCenter[2];
      p[99] = eC[5][0]; p[100] = eC[5][1]; p[101] = eC[5][2];
      p[102] = eC[8][0]; p[103] = eC[8][1]; p[104] = eC[8][2];
      p[105] = tetCenter[0]; p[106] = tetCenter[1]; p[107] = tetCenter[2];
      p += 108;

      c[0] = cellIndex + 0; c[1] = cellIndex + 1;
      c[2] = cellIndex + 2; c[3] = cellIndex + 3;
      c[4] = cellIndex + 4; c[5] = cellIndex + 5;
      c[6] = cellIndex + 6; c[7] = cellIndex + 7;
      c[8] = cellIndex + 8; c[9] = cellIndex + 9;
      c[10] = cellIndex + 10; c[11] = cellIndex + 11;
      c[12] = cellIndex + 12; c[13] = cellIndex + 13;
      c[14] = cellIndex + 14; c[15] = cellIndex + 15;
      c[16] = cellIndex + 16; c[17] = cellIndex + 17;
      c[18] = cellIndex + 18; c[19] = cellIndex + 19;
      c[20] = cellIndex + 20; c[21] = cellIndex + 21;
      c[22] = cellIndex + 22; c[23] = cellIndex + 23;
      c[24] = cellIndex + 24; c[25] = cellIndex + 25;
      c[26] = cellIndex + 26; c[27] = cellIndex + 27;
      c[28] = cellIndex + 28; c[29] = cellIndex + 29;
      c[30] = cellIndex + 30; c[31] = cellIndex + 31;
      c[32] = cellIndex + 32; c[33] = cellIndex + 33;
      c[34] = cellIndex + 34; c[35] = cellIndex + 35;
      c += 36;
      cellIndex += 36;

      m[0] = sparseMSIds[0]; m[1] = sparseMSIds[0];        
      m[2] = sparseMSIds[1]; m[3] = sparseMSIds[1];
      m[4] = sparseMSIds[2]; m[5] = sparseMSIds[2];        
      m[6] = sparseMSIds[3]; m[7] = sparseMSIds[3];       
      m[8] = sparseMSIds[4]; m[9] = sparseMSIds[4];
      m[10] = sparseMSIds[5]; m[11] = sparseMSIds[5];
      m += 12;
    } else { // 2 or 3 labels on tetraeder
      const size_t numTris = tetLookupNumWallTriangles[lookupIndex];
      long long sparseIds[numTris];
      for(size_t t = 0; t < numTris; ++t) {
        sparseIds[t] = msLabelMap[getSparseId(msm[tetVertLabel[t * 2]],
                                     msm[tetVertLabel[(t * 2) + 1]],
                                     numMSCRegions)];

        p[0] = eC[tetEdgeIndices[(t * 3)    ]][0];
        p[1] = eC[tetEdgeIndices[(t * 3)    ]][1];
        p[2] = eC[tetEdgeIndices[(t * 3)    ]][2];
        p[3] = eC[tetEdgeIndices[(t * 3) + 1]][0];
        p[4] = eC[tetEdgeIndices[(t * 3) + 1]][1];
        p[5] = eC[tetEdgeIndices[(t * 3) + 1]][2];
        p[6] = eC[tetEdgeIndices[(t * 3) + 2]][0];
        p[7] = eC[tetEdgeIndices[(t * 3) + 2]][1];
        p[8] = eC[tetEdgeIndices[(t * 3) + 2]][2];
        p += 9;

        c[0] = cellIndex + 0;
        c[1] = cellIndex + 1;
        c[2] = cellIndex + 2;
        c += 3;
        cellIndex += 3;

        m[0] = sparseIds[t];
        m += 1;
      }
    }
  }

#ifdef TTK_ENABLE_OPENMP
} else {
  size_t triangleStartIndex[threadNumber_ + 1];

  // parallel conquer unordered_sets variables
  std::unordered_set<long long>* msSets[threadNumber_];
  int numConquer = threadNumber_ / 2;
  int conquerStep = 1;
  int cS2pow = 1;
  int unmergedSet = -1;
  bool isConquerEven = !(bool)(threadNumber_ % 2);

  #pragma omp parallel num_threads(threadNumber_)
  {
    const int tid = omp_get_thread_num();
    std::vector<std::pair<SimplexId, unsigned char>> validCases;
    size_t numThreadTriangles = 0;
    size_t numThreadIndex = 0;

    msSets[tid] = new std::unordered_set<long long>;
    std::unordered_set<long long>* localMSSet = msSets[tid];

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

      if(tetLookupIsMultiLabel[lookupIndex]) {
        validCases.push_back(std::make_pair(tet, lookupIndex));
        numThreadTriangles += tetLookupNumWallTriangles[lookupIndex];

        const int *tetEdgeIndices = tetLookupWall[lookupIndex];
        if(tetEdgeIndices[0] == 10) { // 4 labels on tetraeder
          long long sparseMSIds[6] = {
            getSparseId(msm[0], msm[1], numMSCRegions),
            getSparseId(msm[0], msm[2], numMSCRegions),
            getSparseId(msm[0], msm[3], numMSCRegions),
            getSparseId(msm[1], msm[2], numMSCRegions),
            getSparseId(msm[1], msm[3], numMSCRegions),
            getSparseId(msm[2], msm[3], numMSCRegions)
          };

          localMSSet->insert({
            sparseMSIds[0], sparseMSIds[1], sparseMSIds[2], 
            sparseMSIds[3], sparseMSIds[4], sparseMSIds[5]});
        } else {
          const int * const unqLabels = tetraederUniqueLabels[lookupIndex];

          if(tetLookupIs2Label[lookupIndex]) {            
            localMSSet->insert(getSparseId(
              msm[0], msm[unqLabels[0]], numMSCRegions));

          } else {
            long long sparseMSIds[3] = {
              getSparseId(msm[0], msm[unqLabels[0]], numMSCRegions),
              getSparseId(msm[0], msm[unqLabels[1]], numMSCRegions),
              getSparseId(msm[unqLabels[0]], msm[unqLabels[1]], numMSCRegions)};

            localMSSet->insert({
              sparseMSIds[0], sparseMSIds[1], sparseMSIds[2]});
          }
        }
      }  
    }

    /************ CONQUER SETS ************/

    while(numConquer > 0) {
      #pragma omp for schedule(static) nowait
      for(int i = 0; i < numConquer * 2; i += 2) {
        msSets[i * cS2pow]->insert(
          msSets[i * cS2pow + cS2pow]->begin(),
          msSets[i * cS2pow + cS2pow]->end());
      }

      #pragma omp single
      { // Fix border problems
        if(!isConquerEven) {
          if(unmergedSet == -1) {
            unmergedSet = 2 * numConquer * cS2pow;
          } else {
            int unmergedSet2 = 2 * numConquer * cS2pow;
            msSets[unmergedSet2]->insert(
              msSets[unmergedSet]->begin(), msSets[unmergedSet]->end());
            unmergedSet = unmergedSet2;
          }
        }
      }

      #pragma omp single
      {
        conquerStep += 1;
        cS2pow = std::pow(2, conquerStep - 1);
        isConquerEven = !(bool)(numConquer % 2);
        numConquer = numConquer / 2;
      }
    }

    // Merge unmerged set with set 0, sort ids and create map
    #pragma omp single
    {
      if(unmergedSet != -1) {
        msSets[0]->insert(
          msSets[unmergedSet]->begin(), msSets[unmergedSet]->end());
      }

      numSparseIDs = msSets[0]->size();
      sparseIdVect.reserve(numSparseIDs);
      std::copy(msSets[0]->begin(), msSets[0]->end(),
        std::back_inserter(sparseIdVect));
      TTK_PSORT(this->threadNumber_, sparseIdVect.begin(), sparseIdVect.end());

      // "sparse id" -> "dense id"
      for(size_t i = 0; i < numSparseIDs; ++i) {
        msLabelMap[sparseIdVect[i]] = i;
      }
    }

    delete(localMSSet);

    /************ END CONQUER SETS ************/

    triangleStartIndex[tid + 1] = numThreadTriangles;

    #pragma omp barrier

    #pragma omp single
    {
      triangleStartIndex[0] = 0;

      // Count triangle number and create iterator start indices
      for(int t = 1; t < threadNumber_ + 1; ++t) {
        triangleStartIndex[t] += triangleStartIndex[t-1];
      }
    }

    const size_t numTriangles = triangleStartIndex[threadNumber_];
    numThreadIndex = triangleStartIndex[tid];

    #pragma omp single nowait
    outputSeparatrices2_points_->resize(9 * numTriangles);

    #pragma omp single nowait
    outputSeparatrices2_cells_connectivity_->resize(3 * numTriangles);

    #pragma omp single nowait
    outputSeparatrices2_cells_mscIds_->resize(numTriangles);

    #pragma omp single
    {
      *outputSeparatrices2_numberOfPoints_ = 3 * numTriangles;
      *outputSeparatrices2_numberOfCells_ = numTriangles;
    }
  
    auto &points = *outputSeparatrices2_points_;
    auto &cellsConn = *outputSeparatrices2_cells_connectivity_;
    auto &cellsMSCIds = *outputSeparatrices2_cells_mscIds_;

    float* p = points.data();
    SimplexId* c = cellsConn.data();
    SimplexId* m = cellsMSCIds.data();

    p += (numThreadIndex * 9);
    c += (numThreadIndex * 3);
    m += numThreadIndex;

    numThreadIndex = 3 * numThreadIndex;

    for (const auto& vCase : validCases)
    {
      const SimplexId& tet = vCase.first;
      const unsigned char& lookupIndex = vCase.second;
      
      const int *tetEdgeIndices = tetLookupWall[lookupIndex];
      const int *tetVertLabel = tetLookupWallLabel[lookupIndex];

      SimplexId vertices[4];
      triangulation.getCellVertex(tet, 0, vertices[0]);
      triangulation.getCellVertex(tet, 1, vertices[1]);
      triangulation.getCellVertex(tet, 2, vertices[2]);
      triangulation.getCellVertex(tet, 3, vertices[3]);

      const SimplexId msm[4] = {
        morseSmaleManifold[vertices[0]], morseSmaleManifold[vertices[1]],
        morseSmaleManifold[vertices[2]], morseSmaleManifold[vertices[3]]};

      float vertPos[4][3];
      triangulation.getVertexPoint(
        vertices[0], vertPos[0][0], vertPos[0][1], vertPos[0][2]);
      triangulation.getVertexPoint(
        vertices[1], vertPos[1][0], vertPos[1][1], vertPos[1][2]);
      triangulation.getVertexPoint(
        vertices[2], vertPos[2][0], vertPos[2][1], vertPos[2][2]);
      triangulation.getVertexPoint(
        vertices[3], vertPos[3][0], vertPos[3][1], vertPos[3][2]);

      float eC[10][3];
      // 6 edge centers
      getCenter(vertPos[0], vertPos[1], eC[0]);
      getCenter(vertPos[0], vertPos[2], eC[1]);
      getCenter(vertPos[0], vertPos[3], eC[2]);
      getCenter(vertPos[1], vertPos[2], eC[3]);
      getCenter(vertPos[1], vertPos[3], eC[4]);
      getCenter(vertPos[2], vertPos[3], eC[5]);

      // 4 triangle centers
      getCenter(vertPos[0], vertPos[1], vertPos[2], eC[6]);
      getCenter(vertPos[0], vertPos[1], vertPos[3], eC[7]);
      getCenter(vertPos[0], vertPos[2], vertPos[3], eC[8]);
      getCenter(vertPos[1], vertPos[2], vertPos[3], eC[9]);

      if(tetEdgeIndices[0] == 10) { // 4 labels on tetraeder
        float tetCenter[3];
        getCenter(vertPos[0], vertPos[1], vertPos[2], vertPos[3], tetCenter);

        long long sparseMSIds[6] = {
          msLabelMap[getSparseId(msm[0], msm[1], numMSCRegions)],
          msLabelMap[getSparseId(msm[0], msm[2], numMSCRegions)],
          msLabelMap[getSparseId(msm[0], msm[3], numMSCRegions)],
          msLabelMap[getSparseId(msm[1], msm[2], numMSCRegions)],
          msLabelMap[getSparseId(msm[1], msm[3], numMSCRegions)],
          msLabelMap[getSparseId(msm[2], msm[3], numMSCRegions)]
        };

        p[0] = eC[7][0]; p[1] = eC[7][1]; p[2] = eC[7][2]; 
        p[3] = eC[0][0]; p[4] = eC[0][1]; p[5] = eC[0][2]; 
        p[6] = tetCenter[0]; p[7] = tetCenter[1]; p[8] = tetCenter[2];
        p[9] = eC[0][0]; p[10] = eC[0][1]; p[11] = eC[0][2]; 
        p[12] = eC[6][0]; p[13] = eC[6][1]; p[14] = eC[6][2]; 
        p[15] = tetCenter[0]; p[16] = tetCenter[1]; p[17] = tetCenter[2];
        p[18] = eC[8][0]; p[19] = eC[8][1]; p[20] = eC[8][2];
        p[21] = eC[1][0]; p[22] = eC[1][1]; p[23] = eC[1][2];
        p[24] = tetCenter[0]; p[25] = tetCenter[1]; p[26] = tetCenter[2];
        p[27] = eC[1][0]; p[28] = eC[1][1]; p[29] = eC[1][2];
        p[30] = eC[6][0]; p[31] = eC[6][1]; p[32] = eC[6][2];
        p[33] = tetCenter[0]; p[34] = tetCenter[1]; p[35] = tetCenter[2];
        p[36] = eC[8][0]; p[37] = eC[8][1]; p[38] = eC[8][2]; 
        p[39] = eC[2][0]; p[40] = eC[2][1]; p[41] = eC[2][2]; 
        p[42] = tetCenter[0]; p[43] = tetCenter[1]; p[44] = tetCenter[2];
        p[45] = eC[2][0]; p[46] = eC[2][1]; p[47] = eC[2][2]; 
        p[48] = eC[7][0]; p[49] = eC[7][1]; p[50] = eC[7][2]; 
        p[51] = tetCenter[0]; p[52] = tetCenter[1]; p[53] = tetCenter[2];
        p[54] = eC[6][0]; p[55] = eC[6][1]; p[56] = eC[6][2];
        p[57] = eC[3][0]; p[58] = eC[3][1]; p[59] = eC[3][2];
        p[60] = tetCenter[0]; p[61] = tetCenter[1]; p[62] = tetCenter[2];
        p[63] = eC[3][0]; p[64] = eC[3][1]; p[65] = eC[3][2];
        p[66] = eC[9][0]; p[67] = eC[9][1]; p[68] = eC[9][2];
        p[69] = tetCenter[0]; p[70] = tetCenter[1]; p[71] = tetCenter[2];
        p[72] = eC[7][0]; p[73] = eC[7][1]; p[74] = eC[7][2];
        p[75] = eC[4][0]; p[76] = eC[4][1]; p[77] = eC[4][2];
        p[78] = tetCenter[0]; p[79] = tetCenter[1]; p[80] = tetCenter[2];
        p[81] = eC[4][0]; p[82] = eC[4][1]; p[83] = eC[4][2];
        p[84] = eC[9][0]; p[85] = eC[9][1]; p[86] = eC[9][2];
        p[87] = tetCenter[0]; p[88] = tetCenter[1]; p[89] = tetCenter[2];
        p[90] = eC[9][0]; p[91] = eC[9][1]; p[92] = eC[9][2]; 
        p[93] = eC[5][0]; p[94] = eC[5][1]; p[95] = eC[5][2]; 
        p[96] = tetCenter[0]; p[97] = tetCenter[1]; p[98] = tetCenter[2];
        p[99] = eC[5][0]; p[100] = eC[5][1]; p[101] = eC[5][2];
        p[102] = eC[8][0]; p[103] = eC[8][1]; p[104] = eC[8][2];
        p[105] = tetCenter[0]; p[106] = tetCenter[1]; p[107] = tetCenter[2];
        p += 108;

        c[0] = numThreadIndex + 0; c[1] = numThreadIndex + 1;
        c[2] = numThreadIndex + 2; c[3] = numThreadIndex + 3;
        c[4] = numThreadIndex + 4; c[5] = numThreadIndex + 5;
        c[6] = numThreadIndex + 6; c[7] = numThreadIndex + 7;
        c[8] = numThreadIndex + 8; c[9] = numThreadIndex + 9;
        c[10] = numThreadIndex + 10; c[11] = numThreadIndex + 11;
        c[12] = numThreadIndex + 12; c[13] = numThreadIndex + 13;
        c[14] = numThreadIndex + 14; c[15] = numThreadIndex + 15;
        c[16] = numThreadIndex + 16; c[17] = numThreadIndex + 17;
        c[18] = numThreadIndex + 18; c[19] = numThreadIndex + 19;
        c[20] = numThreadIndex + 20; c[21] = numThreadIndex + 21;
        c[22] = numThreadIndex + 22; c[23] = numThreadIndex + 23;
        c[24] = numThreadIndex + 24; c[25] = numThreadIndex + 25;
        c[26] = numThreadIndex + 26; c[27] = numThreadIndex + 27;
        c[28] = numThreadIndex + 28; c[29] = numThreadIndex + 29;
        c[30] = numThreadIndex + 30; c[31] = numThreadIndex + 31;
        c[32] = numThreadIndex + 32; c[33] = numThreadIndex + 33;
        c[34] = numThreadIndex + 34; c[35] = numThreadIndex + 35;
        c += 36;
        numThreadIndex += 36;

        m[0] = sparseMSIds[0]; m[1] = sparseMSIds[0];
        m[2] = sparseMSIds[1]; m[3] = sparseMSIds[1];
        m[4] = sparseMSIds[2]; m[5] = sparseMSIds[2];
        m[6] = sparseMSIds[3]; m[7] = sparseMSIds[3];
        m[8] = sparseMSIds[4]; m[9] = sparseMSIds[4];
        m[10] = sparseMSIds[5]; m[11] = sparseMSIds[5];
        m += 12;
        
      } else { // 2 or 3 labels on tetraeder
        const size_t numTris = tetLookupNumWallTriangles[lookupIndex];
        long long sparseIds[numTris];
        for(size_t t = 0; t < numTris; ++t) {
          sparseIds[t] = msLabelMap[getSparseId(msm[tetVertLabel[t * 2]],
                                     msm[tetVertLabel[(t * 2) + 1]],
                                     numMSCRegions)];

          p[0] = eC[tetEdgeIndices[(t * 3)    ]][0];
          p[1] = eC[tetEdgeIndices[(t * 3)    ]][1];
          p[2] = eC[tetEdgeIndices[(t * 3)    ]][2];
          p[3] = eC[tetEdgeIndices[(t * 3) + 1]][0];
          p[4] = eC[tetEdgeIndices[(t * 3) + 1]][1];
          p[5] = eC[tetEdgeIndices[(t * 3) + 1]][2];
          p[6] = eC[tetEdgeIndices[(t * 3) + 2]][0];
          p[7] = eC[tetEdgeIndices[(t * 3) + 2]][1];
          p[8] = eC[tetEdgeIndices[(t * 3) + 2]][2];
          p += 9;

          c[0] = numThreadIndex + 0;
          c[1] = numThreadIndex + 1;
          c[2] = numThreadIndex + 2;
          c += 3;
          numThreadIndex += 3;

          m[0] = sparseIds[t];
          m += 1;
        }
      }
    }
  }
}
#endif // TTK_ENABLE_OPENMP

  this->printMsg("Computed 2-Separatrices 3D[Walls]",
                 1, localTimer.getElapsedTime(), this->threadNumber_);

  return 0;
}

template <typename triangulationType>
int ttk::MorseSmaleSegmentationPL::computeSeparatrices2_3D_split(
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
  size_t numSparseIDs = 0;
  std::vector<long long> sparseIdVect;
  std::map<long long, SimplexId> msLabelMap;

#ifdef TTK_ENABLE_OPENMP
if(threadNumber_ == 1) {
#endif // TTK_ENABLE_OPENMP

  std::vector<std::pair<SimplexId, unsigned char>> validCasesA;
  std::vector<std::pair<SimplexId, unsigned char>> validCasesD;
  std::unordered_set<long long> msLabelSet;
  SimplexId numTriangles = 0;

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
    const unsigned char lookupIndexD = index1D | index2D | index3D;

    if(tetLookupIsMultiLabel[lookupIndexA]) {
      validCasesA.push_back(std::make_pair(tet, lookupIndexA));
      numTriangles += tetLookupNumWallTriangles[lookupIndexA];

      const int *tetEdgeIndices = tetLookupWall[lookupIndexA];
      if(tetEdgeIndices[0] == 10) { // 4 labels on tetraeder
        long long sparseMSIds[6] = {
          getSparseId(asc[0], asc[1], numAscRegions),
          getSparseId(asc[0], asc[2], numAscRegions),
          getSparseId(asc[0], asc[3], numAscRegions),
          getSparseId(asc[1], asc[2], numAscRegions),
          getSparseId(asc[1], asc[3], numAscRegions),
          getSparseId(asc[2], asc[3], numAscRegions)
        };

        msLabelSet.insert({
          sparseMSIds[0], sparseMSIds[1], sparseMSIds[2], 
          sparseMSIds[3], sparseMSIds[4], sparseMSIds[5]});
      } else {
        const int * const unqLabels = tetraederUniqueLabels[lookupIndexA];

        if(tetLookupIs2Label[lookupIndexA]) {            
          msLabelSet.insert(getSparseId(
            asc[0], asc[unqLabels[0]], numAscRegions));

        } else {
          long long sparseMSIds[3] = {
            getSparseId(asc[0], asc[unqLabels[0]], numAscRegions),
            getSparseId(asc[0], asc[unqLabels[1]], numAscRegions),
            getSparseId(asc[unqLabels[0]], asc[unqLabels[1]], numAscRegions)};

          msLabelSet.insert({
            sparseMSIds[0], sparseMSIds[1], sparseMSIds[2]});
        }
      }
    }

    if(tetLookupIsMultiLabel[lookupIndexD]) {
      validCasesD.push_back(std::make_pair(tet, lookupIndexD));
      numTriangles += tetLookupNumWallTriangles[lookupIndexD];

      const int *tetEdgeIndices = tetLookupWall[lookupIndexD];
      if(tetEdgeIndices[0] == 10) { // 4 labels on tetraeder
        long long sparseMSIds[6] = {
          getSparseId(dsc[0], dsc[1], numDscRegions) + ascSparseIDOffset,
          getSparseId(dsc[0], dsc[2], numDscRegions) + ascSparseIDOffset,
          getSparseId(dsc[0], dsc[3], numDscRegions) + ascSparseIDOffset,
          getSparseId(dsc[1], dsc[2], numDscRegions) + ascSparseIDOffset,
          getSparseId(dsc[1], dsc[3], numDscRegions) + ascSparseIDOffset,
          getSparseId(dsc[2], dsc[3], numDscRegions) + ascSparseIDOffset
        };

        msLabelSet.insert({
          sparseMSIds[0], sparseMSIds[1], sparseMSIds[2], 
          sparseMSIds[3], sparseMSIds[4], sparseMSIds[5]});
      } else {
        const int * const unqLabels = tetraederUniqueLabels[lookupIndexD];

        if(tetLookupIs2Label[lookupIndexD]) {            
          msLabelSet.insert(getSparseId(
            dsc[0], dsc[unqLabels[0]], numDscRegions) + ascSparseIDOffset);

        } else {
          long long sparseMSIds[3] = {
            getSparseId(dsc[0], dsc[unqLabels[0]], numDscRegions)
            + ascSparseIDOffset,
            getSparseId(dsc[0], dsc[unqLabels[1]], numDscRegions)
            + ascSparseIDOffset,
            getSparseId(dsc[unqLabels[0]], dsc[unqLabels[1]], numDscRegions)
            + ascSparseIDOffset};

          msLabelSet.insert({
            sparseMSIds[0], sparseMSIds[1], sparseMSIds[2]});
        }
      }
    }
  }

  numSparseIDs = msLabelSet.size();
  sparseIdVect.reserve(numSparseIDs);
  std::copy(msLabelSet.begin(), msLabelSet.end(),
    std::back_inserter(sparseIdVect));
  TTK_PSORT(this->threadNumber_, sparseIdVect.begin(), sparseIdVect.end());

  // "sparse id" -> "dense id"
  for(size_t i = 0; i < numSparseIDs; ++i) {
    msLabelMap[sparseIdVect[i]] = i;
  }

  outputSeparatrices2_points_->resize(9 * numTriangles);
  outputSeparatrices2_cells_connectivity_->resize(3 * numTriangles);
  outputSeparatrices2_cells_mscIds_->resize(numTriangles);
  *outputSeparatrices2_numberOfPoints_ = 3 * numTriangles;
  *outputSeparatrices2_numberOfCells_ = numTriangles;

  auto &points = *outputSeparatrices2_points_;
  auto &cellsConn = *outputSeparatrices2_cells_connectivity_;
  auto &cellsMSCIds = *outputSeparatrices2_cells_mscIds_;

  float* p = points.data();
  SimplexId* c = cellsConn.data();
  SimplexId* m = cellsMSCIds.data();
  SimplexId cellIndex = 0;

  for (const auto& vCase : validCasesA)
  {
    const SimplexId& tet = vCase.first;
    const unsigned char& lookupIndex = vCase.second;

    const int *tetEdgeIndices = tetLookupWall[lookupIndex];
    const int *tetVertLabel = tetLookupWallLabel[lookupIndex];

    SimplexId vertices[4];
    triangulation.getCellVertex(tet, 0, vertices[0]);
    triangulation.getCellVertex(tet, 1, vertices[1]);
    triangulation.getCellVertex(tet, 2, vertices[2]);
    triangulation.getCellVertex(tet, 3, vertices[3]);

    const SimplexId msm[4] = {
      ascManifold[vertices[0]], ascManifold[vertices[1]],
      ascManifold[vertices[2]], ascManifold[vertices[3]]};

    float vertPos[4][3];
    triangulation.getVertexPoint(
      vertices[0], vertPos[0][0], vertPos[0][1], vertPos[0][2]);
    triangulation.getVertexPoint(
      vertices[1], vertPos[1][0], vertPos[1][1], vertPos[1][2]);
    triangulation.getVertexPoint(
      vertices[2], vertPos[2][0], vertPos[2][1], vertPos[2][2]);
    triangulation.getVertexPoint(
      vertices[3], vertPos[3][0], vertPos[3][1], vertPos[3][2]);

    float eC[10][3];
    // 6 edge centers
    getCenter(vertPos[0], vertPos[1], eC[0]);
    getCenter(vertPos[0], vertPos[2], eC[1]);
    getCenter(vertPos[0], vertPos[3], eC[2]);
    getCenter(vertPos[1], vertPos[2], eC[3]);
    getCenter(vertPos[1], vertPos[3], eC[4]);
    getCenter(vertPos[2], vertPos[3], eC[5]);

    // 4 triangle centers
    getCenter(vertPos[0], vertPos[1], vertPos[2], eC[6]);
    getCenter(vertPos[0], vertPos[1], vertPos[3], eC[7]);
    getCenter(vertPos[0], vertPos[2], vertPos[3], eC[8]);
    getCenter(vertPos[1], vertPos[2], vertPos[3], eC[9]);

    if(tetEdgeIndices[0] == 10) { // 4 labels on tetraeder
      float tetCenter[3];
      getCenter(vertPos[0], vertPos[1], vertPos[2], vertPos[3], tetCenter);

      long long sparseMSIds[6] = {
        msLabelMap[getSparseId(msm[0], msm[1], numAscRegions)],
        msLabelMap[getSparseId(msm[0], msm[2], numAscRegions)],
        msLabelMap[getSparseId(msm[0], msm[3], numAscRegions)],
        msLabelMap[getSparseId(msm[1], msm[2], numAscRegions)],
        msLabelMap[getSparseId(msm[1], msm[3], numAscRegions)],
        msLabelMap[getSparseId(msm[2], msm[3], numAscRegions)]
      };

      p[0] = eC[7][0]; p[1] = eC[7][1]; p[2] = eC[7][2]; 
      p[3] = eC[0][0]; p[4] = eC[0][1]; p[5] = eC[0][2]; 
      p[6] = tetCenter[0]; p[7] = tetCenter[1]; p[8] = tetCenter[2];
      p[9] = eC[0][0]; p[10] = eC[0][1]; p[11] = eC[0][2]; 
      p[12] = eC[6][0]; p[13] = eC[6][1]; p[14] = eC[6][2]; 
      p[15] = tetCenter[0]; p[16] = tetCenter[1]; p[17] = tetCenter[2];
      p[18] = eC[8][0]; p[19] = eC[8][1]; p[20] = eC[8][2];
      p[21] = eC[1][0]; p[22] = eC[1][1]; p[23] = eC[1][2];
      p[24] = tetCenter[0]; p[25] = tetCenter[1]; p[26] = tetCenter[2];
      p[27] = eC[1][0]; p[28] = eC[1][1]; p[29] = eC[1][2];
      p[30] = eC[6][0]; p[31] = eC[6][1]; p[32] = eC[6][2];
      p[33] = tetCenter[0]; p[34] = tetCenter[1]; p[35] = tetCenter[2];
      p[36] = eC[8][0]; p[37] = eC[8][1]; p[38] = eC[8][2]; 
      p[39] = eC[2][0]; p[40] = eC[2][1]; p[41] = eC[2][2]; 
      p[42] = tetCenter[0]; p[43] = tetCenter[1]; p[44] = tetCenter[2];
      p[45] = eC[2][0]; p[46] = eC[2][1]; p[47] = eC[2][2]; 
      p[48] = eC[7][0]; p[49] = eC[7][1]; p[50] = eC[7][2]; 
      p[51] = tetCenter[0]; p[52] = tetCenter[1]; p[53] = tetCenter[2];
      p[54] = eC[6][0]; p[55] = eC[6][1]; p[56] = eC[6][2];
      p[57] = eC[3][0]; p[58] = eC[3][1]; p[59] = eC[3][2];
      p[60] = tetCenter[0]; p[61] = tetCenter[1]; p[62] = tetCenter[2];
      p[63] = eC[3][0]; p[64] = eC[3][1]; p[65] = eC[3][2];
      p[66] = eC[9][0]; p[67] = eC[9][1]; p[68] = eC[9][2];
      p[69] = tetCenter[0]; p[70] = tetCenter[1]; p[71] = tetCenter[2];
      p[72] = eC[7][0]; p[73] = eC[7][1]; p[74] = eC[7][2];
      p[75] = eC[4][0]; p[76] = eC[4][1]; p[77] = eC[4][2];
      p[78] = tetCenter[0]; p[79] = tetCenter[1]; p[80] = tetCenter[2];
      p[81] = eC[4][0]; p[82] = eC[4][1]; p[83] = eC[4][2];
      p[84] = eC[9][0]; p[85] = eC[9][1]; p[86] = eC[9][2];
      p[87] = tetCenter[0]; p[88] = tetCenter[1]; p[89] = tetCenter[2];
      p[90] = eC[9][0]; p[91] = eC[9][1]; p[92] = eC[9][2]; 
      p[93] = eC[5][0]; p[94] = eC[5][1]; p[95] = eC[5][2]; 
      p[96] = tetCenter[0]; p[97] = tetCenter[1]; p[98] = tetCenter[2];
      p[99] = eC[5][0]; p[100] = eC[5][1]; p[101] = eC[5][2];
      p[102] = eC[8][0]; p[103] = eC[8][1]; p[104] = eC[8][2];
      p[105] = tetCenter[0]; p[106] = tetCenter[1]; p[107] = tetCenter[2];
      p += 108;

      c[0] = cellIndex + 0; c[1] = cellIndex + 1;
      c[2] = cellIndex + 2; c[3] = cellIndex + 3;
      c[4] = cellIndex + 4; c[5] = cellIndex + 5;
      c[6] = cellIndex + 6; c[7] = cellIndex + 7;
      c[8] = cellIndex + 8; c[9] = cellIndex + 9;
      c[10] = cellIndex + 10; c[11] = cellIndex + 11;
      c[12] = cellIndex + 12; c[13] = cellIndex + 13;
      c[14] = cellIndex + 14; c[15] = cellIndex + 15;
      c[16] = cellIndex + 16; c[17] = cellIndex + 17;
      c[18] = cellIndex + 18; c[19] = cellIndex + 19;
      c[20] = cellIndex + 20; c[21] = cellIndex + 21;
      c[22] = cellIndex + 22; c[23] = cellIndex + 23;
      c[24] = cellIndex + 24; c[25] = cellIndex + 25;
      c[26] = cellIndex + 26; c[27] = cellIndex + 27;
      c[28] = cellIndex + 28; c[29] = cellIndex + 29;
      c[30] = cellIndex + 30; c[31] = cellIndex + 31;
      c[32] = cellIndex + 32; c[33] = cellIndex + 33;
      c[34] = cellIndex + 34; c[35] = cellIndex + 35;
      c += 36;
      cellIndex += 36;

      m[0] = sparseMSIds[0]; m[1] = sparseMSIds[0];        
      m[2] = sparseMSIds[1]; m[3] = sparseMSIds[1];
      m[4] = sparseMSIds[2]; m[5] = sparseMSIds[2];        
      m[6] = sparseMSIds[3]; m[7] = sparseMSIds[3];       
      m[8] = sparseMSIds[4]; m[9] = sparseMSIds[4];
      m[10] = sparseMSIds[5]; m[11] = sparseMSIds[5];
      m += 12;
    } else { // 2 or 3 labels on tetraeder
      const size_t numTris = tetLookupNumWallTriangles[lookupIndex];
      long long sparseIds[numTris];
      for(size_t t = 0; t < numTris; ++t) {
        sparseIds[t] = msLabelMap[getSparseId(msm[tetVertLabel[t * 2]],
                                     msm[tetVertLabel[(t * 2) + 1]],
                                     numAscRegions)];

        p[0] = eC[tetEdgeIndices[(t * 3)    ]][0];
        p[1] = eC[tetEdgeIndices[(t * 3)    ]][1];
        p[2] = eC[tetEdgeIndices[(t * 3)    ]][2];
        p[3] = eC[tetEdgeIndices[(t * 3) + 1]][0];
        p[4] = eC[tetEdgeIndices[(t * 3) + 1]][1];
        p[5] = eC[tetEdgeIndices[(t * 3) + 1]][2];
        p[6] = eC[tetEdgeIndices[(t * 3) + 2]][0];
        p[7] = eC[tetEdgeIndices[(t * 3) + 2]][1];
        p[8] = eC[tetEdgeIndices[(t * 3) + 2]][2];
        p += 9;

        c[0] = cellIndex + 0;
        c[1] = cellIndex + 1;
        c[2] = cellIndex + 2;
        c += 3;
        cellIndex += 3;

        m[0] = sparseIds[t];
        m += 1;
      }
    }
  }

  for (const auto& vCase : validCasesD)
  {
    const SimplexId& tet = vCase.first;
    const unsigned char& lookupIndex = vCase.second;

    const int *tetEdgeIndices = tetLookupWall[lookupIndex];
    const int *tetVertLabel = tetLookupWallLabel[lookupIndex];

    SimplexId vertices[4];
    triangulation.getCellVertex(tet, 0, vertices[0]);
    triangulation.getCellVertex(tet, 1, vertices[1]);
    triangulation.getCellVertex(tet, 2, vertices[2]);
    triangulation.getCellVertex(tet, 3, vertices[3]);

    const SimplexId msm[4] = {
      dscManifold[vertices[0]], dscManifold[vertices[1]],
      dscManifold[vertices[2]], dscManifold[vertices[3]]};

    float vertPos[4][3];
    triangulation.getVertexPoint(
      vertices[0], vertPos[0][0], vertPos[0][1], vertPos[0][2]);
    triangulation.getVertexPoint(
      vertices[1], vertPos[1][0], vertPos[1][1], vertPos[1][2]);
    triangulation.getVertexPoint(
      vertices[2], vertPos[2][0], vertPos[2][1], vertPos[2][2]);
    triangulation.getVertexPoint(
      vertices[3], vertPos[3][0], vertPos[3][1], vertPos[3][2]);

    float eC[10][3];
    // 6 edge centers
    getCenter(vertPos[0], vertPos[1], eC[0]);
    getCenter(vertPos[0], vertPos[2], eC[1]);
    getCenter(vertPos[0], vertPos[3], eC[2]);
    getCenter(vertPos[1], vertPos[2], eC[3]);
    getCenter(vertPos[1], vertPos[3], eC[4]);
    getCenter(vertPos[2], vertPos[3], eC[5]);

    // 4 triangle centers
    getCenter(vertPos[0], vertPos[1], vertPos[2], eC[6]);
    getCenter(vertPos[0], vertPos[1], vertPos[3], eC[7]);
    getCenter(vertPos[0], vertPos[2], vertPos[3], eC[8]);
    getCenter(vertPos[1], vertPos[2], vertPos[3], eC[9]);

    if(tetEdgeIndices[0] == 10) { // 4 labels on tetraeder
      float tetCenter[3];
      getCenter(vertPos[0], vertPos[1], vertPos[2], vertPos[3], tetCenter);

      long long sparseMSIds[6] = {
        msLabelMap[getSparseId(msm[0], msm[1], numDscRegions)
        + ascSparseIDOffset],
        msLabelMap[getSparseId(msm[0], msm[2], numDscRegions)
        + ascSparseIDOffset],
        msLabelMap[getSparseId(msm[0], msm[3], numDscRegions)
        + ascSparseIDOffset],
        msLabelMap[getSparseId(msm[1], msm[2], numDscRegions)
        + ascSparseIDOffset],
        msLabelMap[getSparseId(msm[1], msm[3], numDscRegions)
        + ascSparseIDOffset],
        msLabelMap[getSparseId(msm[2], msm[3], numDscRegions)
        + ascSparseIDOffset]
      };

      p[0] = eC[7][0]; p[1] = eC[7][1]; p[2] = eC[7][2]; 
      p[3] = eC[0][0]; p[4] = eC[0][1]; p[5] = eC[0][2]; 
      p[6] = tetCenter[0]; p[7] = tetCenter[1]; p[8] = tetCenter[2];
      p[9] = eC[0][0]; p[10] = eC[0][1]; p[11] = eC[0][2]; 
      p[12] = eC[6][0]; p[13] = eC[6][1]; p[14] = eC[6][2]; 
      p[15] = tetCenter[0]; p[16] = tetCenter[1]; p[17] = tetCenter[2];
      p[18] = eC[8][0]; p[19] = eC[8][1]; p[20] = eC[8][2];
      p[21] = eC[1][0]; p[22] = eC[1][1]; p[23] = eC[1][2];
      p[24] = tetCenter[0]; p[25] = tetCenter[1]; p[26] = tetCenter[2];
      p[27] = eC[1][0]; p[28] = eC[1][1]; p[29] = eC[1][2];
      p[30] = eC[6][0]; p[31] = eC[6][1]; p[32] = eC[6][2];
      p[33] = tetCenter[0]; p[34] = tetCenter[1]; p[35] = tetCenter[2];
      p[36] = eC[8][0]; p[37] = eC[8][1]; p[38] = eC[8][2]; 
      p[39] = eC[2][0]; p[40] = eC[2][1]; p[41] = eC[2][2]; 
      p[42] = tetCenter[0]; p[43] = tetCenter[1]; p[44] = tetCenter[2];
      p[45] = eC[2][0]; p[46] = eC[2][1]; p[47] = eC[2][2]; 
      p[48] = eC[7][0]; p[49] = eC[7][1]; p[50] = eC[7][2]; 
      p[51] = tetCenter[0]; p[52] = tetCenter[1]; p[53] = tetCenter[2];
      p[54] = eC[6][0]; p[55] = eC[6][1]; p[56] = eC[6][2];
      p[57] = eC[3][0]; p[58] = eC[3][1]; p[59] = eC[3][2];
      p[60] = tetCenter[0]; p[61] = tetCenter[1]; p[62] = tetCenter[2];
      p[63] = eC[3][0]; p[64] = eC[3][1]; p[65] = eC[3][2];
      p[66] = eC[9][0]; p[67] = eC[9][1]; p[68] = eC[9][2];
      p[69] = tetCenter[0]; p[70] = tetCenter[1]; p[71] = tetCenter[2];
      p[72] = eC[7][0]; p[73] = eC[7][1]; p[74] = eC[7][2];
      p[75] = eC[4][0]; p[76] = eC[4][1]; p[77] = eC[4][2];
      p[78] = tetCenter[0]; p[79] = tetCenter[1]; p[80] = tetCenter[2];
      p[81] = eC[4][0]; p[82] = eC[4][1]; p[83] = eC[4][2];
      p[84] = eC[9][0]; p[85] = eC[9][1]; p[86] = eC[9][2];
      p[87] = tetCenter[0]; p[88] = tetCenter[1]; p[89] = tetCenter[2];
      p[90] = eC[9][0]; p[91] = eC[9][1]; p[92] = eC[9][2]; 
      p[93] = eC[5][0]; p[94] = eC[5][1]; p[95] = eC[5][2]; 
      p[96] = tetCenter[0]; p[97] = tetCenter[1]; p[98] = tetCenter[2];
      p[99] = eC[5][0]; p[100] = eC[5][1]; p[101] = eC[5][2];
      p[102] = eC[8][0]; p[103] = eC[8][1]; p[104] = eC[8][2];
      p[105] = tetCenter[0]; p[106] = tetCenter[1]; p[107] = tetCenter[2];
      p += 108;

      c[0] = cellIndex + 0; c[1] = cellIndex + 1;
      c[2] = cellIndex + 2; c[3] = cellIndex + 3;
      c[4] = cellIndex + 4; c[5] = cellIndex + 5;
      c[6] = cellIndex + 6; c[7] = cellIndex + 7;
      c[8] = cellIndex + 8; c[9] = cellIndex + 9;
      c[10] = cellIndex + 10; c[11] = cellIndex + 11;
      c[12] = cellIndex + 12; c[13] = cellIndex + 13;
      c[14] = cellIndex + 14; c[15] = cellIndex + 15;
      c[16] = cellIndex + 16; c[17] = cellIndex + 17;
      c[18] = cellIndex + 18; c[19] = cellIndex + 19;
      c[20] = cellIndex + 20; c[21] = cellIndex + 21;
      c[22] = cellIndex + 22; c[23] = cellIndex + 23;
      c[24] = cellIndex + 24; c[25] = cellIndex + 25;
      c[26] = cellIndex + 26; c[27] = cellIndex + 27;
      c[28] = cellIndex + 28; c[29] = cellIndex + 29;
      c[30] = cellIndex + 30; c[31] = cellIndex + 31;
      c[32] = cellIndex + 32; c[33] = cellIndex + 33;
      c[34] = cellIndex + 34; c[35] = cellIndex + 35;
      c += 36;
      cellIndex += 36;

      m[0] = sparseMSIds[0]; m[1] = sparseMSIds[0];        
      m[2] = sparseMSIds[1]; m[3] = sparseMSIds[1];
      m[4] = sparseMSIds[2]; m[5] = sparseMSIds[2];        
      m[6] = sparseMSIds[3]; m[7] = sparseMSIds[3];       
      m[8] = sparseMSIds[4]; m[9] = sparseMSIds[4];
      m[10] = sparseMSIds[5]; m[11] = sparseMSIds[5];
      m += 12;
    } else { // 2 or 3 labels on tetraeder
      const size_t numTris = tetLookupNumWallTriangles[lookupIndex];
      long long sparseIds[numTris];
      for(size_t t = 0; t < numTris; ++t) {
        sparseIds[t] = msLabelMap[getSparseId(msm[tetVertLabel[t * 2]],
                                  msm[tetVertLabel[(t * 2) + 1]],
                                  numDscRegions) + ascSparseIDOffset];

        p[0] = eC[tetEdgeIndices[(t * 3)    ]][0];
        p[1] = eC[tetEdgeIndices[(t * 3)    ]][1];
        p[2] = eC[tetEdgeIndices[(t * 3)    ]][2];
        p[3] = eC[tetEdgeIndices[(t * 3) + 1]][0];
        p[4] = eC[tetEdgeIndices[(t * 3) + 1]][1];
        p[5] = eC[tetEdgeIndices[(t * 3) + 1]][2];
        p[6] = eC[tetEdgeIndices[(t * 3) + 2]][0];
        p[7] = eC[tetEdgeIndices[(t * 3) + 2]][1];
        p[8] = eC[tetEdgeIndices[(t * 3) + 2]][2];
        p += 9;

        c[0] = cellIndex + 0;
        c[1] = cellIndex + 1;
        c[2] = cellIndex + 2;
        c += 3;
        cellIndex += 3;

        m[0] = sparseIds[t];
        m += 1;
      }
    }
  }

#ifdef TTK_ENABLE_OPENMP
} else {
  size_t triangleStartIndex[threadNumber_ + 1];

  // parallel conquer unordered_sets variables
  std::unordered_set<long long>* msSets[threadNumber_];
  int numConquer = threadNumber_ / 2;
  int conquerStep = 1;
  int cS2pow = 1;
  int unmergedSet = -1;
  bool isConquerEven = !(bool)(threadNumber_ % 2);

  #pragma omp parallel num_threads(threadNumber_)
  {
    const int tid = omp_get_thread_num();
    std::vector<std::pair<SimplexId, unsigned char>> validCasesA;
    std::vector<std::pair<SimplexId, unsigned char>> validCasesD;
    size_t numThreadTriangles = 0;
    size_t numThreadIndex = 0;

    msSets[tid] = new std::unordered_set<long long>;
    std::unordered_set<long long>* localMSSet = msSets[tid];

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
      const unsigned char lookupIndexD = index1D | index2D | index3D;

      if(tetLookupIsMultiLabel[lookupIndexA]) {
        validCasesA.push_back(std::make_pair(tet, lookupIndexA));
        numThreadTriangles += tetLookupNumFineTriangles[lookupIndexA];

        const int *tetEdgeIndices = tetLookupWall[lookupIndexA];
        if(tetEdgeIndices[0] == 10) { // 4 labels on tetraeder
          long long sparseMSIds[6] = {
            getSparseId(asc[0], asc[1], numAscRegions),
            getSparseId(asc[0], asc[2], numAscRegions),
            getSparseId(asc[0], asc[3], numAscRegions),
            getSparseId(asc[1], asc[2], numAscRegions),
            getSparseId(asc[1], asc[3], numAscRegions),
            getSparseId(asc[2], asc[3], numAscRegions)
          };

          localMSSet->insert({
            sparseMSIds[0], sparseMSIds[1], sparseMSIds[2], 
            sparseMSIds[3], sparseMSIds[4], sparseMSIds[5]});
        } else {
          const int * const unqLabels = tetraederUniqueLabels[lookupIndexA];

          if(tetLookupIs2Label[lookupIndexA]) {            
            localMSSet->insert(getSparseId(
              asc[0], asc[unqLabels[0]], numAscRegions));

          } else {
            long long sparseMSIds[3] = {
              getSparseId(asc[0], asc[unqLabels[0]], numAscRegions),
              getSparseId(asc[0], asc[unqLabels[1]], numAscRegions),
              getSparseId(asc[unqLabels[0]], asc[unqLabels[1]], numAscRegions)};

            localMSSet->insert({
              sparseMSIds[0], sparseMSIds[1], sparseMSIds[2]});
          }
        }
      }

      if(tetLookupIsMultiLabel[lookupIndexD]) {
        validCasesD.push_back(std::make_pair(tet, lookupIndexD));
        numThreadTriangles += tetLookupNumFineTriangles[lookupIndexD];

        const int *tetEdgeIndices = tetLookupWall[lookupIndexD];
        if(tetEdgeIndices[0] == 10) { // 4 labels on tetraeder
          long long sparseMSIds[6] = {
            getSparseId(dsc[0], dsc[1], numDscRegions)
            + ascSparseIDOffset,
            getSparseId(dsc[0], dsc[2], numDscRegions)
            + ascSparseIDOffset,
            getSparseId(dsc[0], dsc[3], numDscRegions)
            + ascSparseIDOffset,
            getSparseId(dsc[1], dsc[2], numDscRegions)
            + ascSparseIDOffset,
            getSparseId(dsc[1], dsc[3], numDscRegions)
            + ascSparseIDOffset,
            getSparseId(dsc[2], dsc[3], numDscRegions)
            + ascSparseIDOffset
          };

          localMSSet->insert({
            sparseMSIds[0], sparseMSIds[1], sparseMSIds[2], 
            sparseMSIds[3], sparseMSIds[4], sparseMSIds[5]});
        } else {
          const int * const unqLabels = tetraederUniqueLabels[lookupIndexD];

          if(tetLookupIs2Label[lookupIndexD]) {            
            localMSSet->insert(getSparseId(
              dsc[0], dsc[unqLabels[0]], numDscRegions) + ascSparseIDOffset);

          } else {
            long long sparseMSIds[3] = {
              getSparseId(dsc[0], dsc[unqLabels[0]], numDscRegions)
              + ascSparseIDOffset,
              getSparseId(dsc[0], dsc[unqLabels[1]], numDscRegions)
              + ascSparseIDOffset,
              getSparseId(dsc[unqLabels[0]], dsc[unqLabels[1]], numDscRegions)
              + ascSparseIDOffset};

            localMSSet->insert({
              sparseMSIds[0], sparseMSIds[1], sparseMSIds[2]});
          }
        }
      }
    }

    /************ CONQUER SETS ************/

    while(numConquer > 0) {
      #pragma omp for schedule(static) nowait
      for(int i = 0; i < numConquer * 2; i += 2) {
        msSets[i * cS2pow]->insert(
          msSets[i * cS2pow + cS2pow]->begin(),
          msSets[i * cS2pow + cS2pow]->end());
      }

      #pragma omp single
      { // Fix border problems
        if(!isConquerEven) {
          if(unmergedSet == -1) {
            unmergedSet = 2 * numConquer * cS2pow;
          } else {
            int unmergedSet2 = 2 * numConquer * cS2pow;
            msSets[unmergedSet2]->insert(
              msSets[unmergedSet]->begin(), msSets[unmergedSet]->end());
            unmergedSet = unmergedSet2;
          }
        }
      }

      #pragma omp single
      {
        conquerStep += 1;
        cS2pow = std::pow(2, conquerStep - 1);
        isConquerEven = !(bool)(numConquer % 2);
        numConquer = numConquer / 2;
      }
    }

    // Merge unmerged set with set 0, sort ids and create map
    #pragma omp single
    {
      if(unmergedSet != -1) {
        msSets[0]->insert(
          msSets[unmergedSet]->begin(), msSets[unmergedSet]->end());
      }

      numSparseIDs = msSets[0]->size();
      sparseIdVect.reserve(numSparseIDs);
      std::copy(msSets[0]->begin(), msSets[0]->end(),
        std::back_inserter(sparseIdVect));
      TTK_PSORT(this->threadNumber_, sparseIdVect.begin(), sparseIdVect.end());

      // "sparse id" -> "dense id"
      for(size_t i = 0; i < numSparseIDs; ++i) {
        msLabelMap[sparseIdVect[i]] = i;
      }
    }

    delete(localMSSet);

    /************ END CONQUER SETS ************/

    triangleStartIndex[tid + 1] = numThreadTriangles;

    #pragma omp barrier

    #pragma omp single
    {
      triangleStartIndex[0] = 0;

      // Count triangle number and create iterator start indices
      for(int t = 1; t < threadNumber_ + 1; ++t) {
        triangleStartIndex[t] += triangleStartIndex[t-1];
      }
    }

    const size_t numTriangles = triangleStartIndex[threadNumber_];
    numThreadIndex = triangleStartIndex[tid];

    #pragma omp single nowait
    outputSeparatrices2_points_->resize(9 * numTriangles);

    #pragma omp single nowait
    outputSeparatrices2_cells_connectivity_->resize(3 * numTriangles);

    #pragma omp single nowait
    outputSeparatrices2_cells_mscIds_->resize(numTriangles);

    #pragma omp single
    {
      *outputSeparatrices2_numberOfPoints_ = 3 * numTriangles;
      *outputSeparatrices2_numberOfCells_ = numTriangles;
    }
  
    auto &points = *outputSeparatrices2_points_;
    auto &cellsConn = *outputSeparatrices2_cells_connectivity_;
    auto &cellsMSCIds = *outputSeparatrices2_cells_mscIds_;

    float* p = points.data();
    SimplexId* c = cellsConn.data();
    SimplexId* m = cellsMSCIds.data();

    p += (numThreadIndex * 9);
    c += (numThreadIndex * 3);
    m += numThreadIndex;

    numThreadIndex = 3 * numThreadIndex;

    for (const auto& vCase : validCasesA)
    {
      const SimplexId& tet = vCase.first;
      const unsigned char& lookupIndex = vCase.second;
      
      const int *tetEdgeIndices = tetLookupWall[lookupIndex];
      const int *tetVertLabel = tetLookupWallLabel[lookupIndex];

      SimplexId vertices[4];
      triangulation.getCellVertex(tet, 0, vertices[0]);
      triangulation.getCellVertex(tet, 1, vertices[1]);
      triangulation.getCellVertex(tet, 2, vertices[2]);
      triangulation.getCellVertex(tet, 3, vertices[3]);

      const SimplexId asc[4] = {
        ascManifold[vertices[0]], ascManifold[vertices[1]],
        ascManifold[vertices[2]], ascManifold[vertices[3]]};

      float vertPos[4][3];
      triangulation.getVertexPoint(
        vertices[0], vertPos[0][0], vertPos[0][1], vertPos[0][2]);
      triangulation.getVertexPoint(
        vertices[1], vertPos[1][0], vertPos[1][1], vertPos[1][2]);
      triangulation.getVertexPoint(
        vertices[2], vertPos[2][0], vertPos[2][1], vertPos[2][2]);
      triangulation.getVertexPoint(
        vertices[3], vertPos[3][0], vertPos[3][1], vertPos[3][2]);

      float eC[10][3];
      // 6 edge centers
      getCenter(vertPos[0], vertPos[1], eC[0]);
      getCenter(vertPos[0], vertPos[2], eC[1]);
      getCenter(vertPos[0], vertPos[3], eC[2]);
      getCenter(vertPos[1], vertPos[2], eC[3]);
      getCenter(vertPos[1], vertPos[3], eC[4]);
      getCenter(vertPos[2], vertPos[3], eC[5]);

      // 4 triangle centers
      getCenter(vertPos[0], vertPos[1], vertPos[2], eC[6]);
      getCenter(vertPos[0], vertPos[1], vertPos[3], eC[7]);
      getCenter(vertPos[0], vertPos[2], vertPos[3], eC[8]);
      getCenter(vertPos[1], vertPos[2], vertPos[3], eC[9]);

      if(tetEdgeIndices[0] == 10) { // 4 labels on tetraeder
        float tetCenter[3];
        getCenter(vertPos[0], vertPos[1], vertPos[2], vertPos[3], tetCenter);

        long long sparseMSIds[6] = {
          msLabelMap[getSparseId(asc[0], asc[1], numAscRegions)],
          msLabelMap[getSparseId(asc[0], asc[2], numAscRegions)],
          msLabelMap[getSparseId(asc[0], asc[3], numAscRegions)],
          msLabelMap[getSparseId(asc[1], asc[2], numAscRegions)],
          msLabelMap[getSparseId(asc[1], asc[3], numAscRegions)],
          msLabelMap[getSparseId(asc[2], asc[3], numAscRegions)]
        };

        p[0] = eC[7][0]; p[1] = eC[7][1]; p[2] = eC[7][2]; 
        p[3] = eC[0][0]; p[4] = eC[0][1]; p[5] = eC[0][2]; 
        p[6] = tetCenter[0]; p[7] = tetCenter[1]; p[8] = tetCenter[2];
        p[9] = eC[0][0]; p[10] = eC[0][1]; p[11] = eC[0][2]; 
        p[12] = eC[6][0]; p[13] = eC[6][1]; p[14] = eC[6][2]; 
        p[15] = tetCenter[0]; p[16] = tetCenter[1]; p[17] = tetCenter[2];
        p[18] = eC[8][0]; p[19] = eC[8][1]; p[20] = eC[8][2];
        p[21] = eC[1][0]; p[22] = eC[1][1]; p[23] = eC[1][2];
        p[24] = tetCenter[0]; p[25] = tetCenter[1]; p[26] = tetCenter[2];
        p[27] = eC[1][0]; p[28] = eC[1][1]; p[29] = eC[1][2];
        p[30] = eC[6][0]; p[31] = eC[6][1]; p[32] = eC[6][2];
        p[33] = tetCenter[0]; p[34] = tetCenter[1]; p[35] = tetCenter[2];
        p[36] = eC[8][0]; p[37] = eC[8][1]; p[38] = eC[8][2]; 
        p[39] = eC[2][0]; p[40] = eC[2][1]; p[41] = eC[2][2]; 
        p[42] = tetCenter[0]; p[43] = tetCenter[1]; p[44] = tetCenter[2];
        p[45] = eC[2][0]; p[46] = eC[2][1]; p[47] = eC[2][2]; 
        p[48] = eC[7][0]; p[49] = eC[7][1]; p[50] = eC[7][2]; 
        p[51] = tetCenter[0]; p[52] = tetCenter[1]; p[53] = tetCenter[2];
        p[54] = eC[6][0]; p[55] = eC[6][1]; p[56] = eC[6][2];
        p[57] = eC[3][0]; p[58] = eC[3][1]; p[59] = eC[3][2];
        p[60] = tetCenter[0]; p[61] = tetCenter[1]; p[62] = tetCenter[2];
        p[63] = eC[3][0]; p[64] = eC[3][1]; p[65] = eC[3][2];
        p[66] = eC[9][0]; p[67] = eC[9][1]; p[68] = eC[9][2];
        p[69] = tetCenter[0]; p[70] = tetCenter[1]; p[71] = tetCenter[2];
        p[72] = eC[7][0]; p[73] = eC[7][1]; p[74] = eC[7][2];
        p[75] = eC[4][0]; p[76] = eC[4][1]; p[77] = eC[4][2];
        p[78] = tetCenter[0]; p[79] = tetCenter[1]; p[80] = tetCenter[2];
        p[81] = eC[4][0]; p[82] = eC[4][1]; p[83] = eC[4][2];
        p[84] = eC[9][0]; p[85] = eC[9][1]; p[86] = eC[9][2];
        p[87] = tetCenter[0]; p[88] = tetCenter[1]; p[89] = tetCenter[2];
        p[90] = eC[9][0]; p[91] = eC[9][1]; p[92] = eC[9][2]; 
        p[93] = eC[5][0]; p[94] = eC[5][1]; p[95] = eC[5][2]; 
        p[96] = tetCenter[0]; p[97] = tetCenter[1]; p[98] = tetCenter[2];
        p[99] = eC[5][0]; p[100] = eC[5][1]; p[101] = eC[5][2];
        p[102] = eC[8][0]; p[103] = eC[8][1]; p[104] = eC[8][2];
        p[105] = tetCenter[0]; p[106] = tetCenter[1]; p[107] = tetCenter[2];
        p += 108;

        c[0] = numThreadIndex + 0; c[1] = numThreadIndex + 1;
        c[2] = numThreadIndex + 2; c[3] = numThreadIndex + 3;
        c[4] = numThreadIndex + 4; c[5] = numThreadIndex + 5;
        c[6] = numThreadIndex + 6; c[7] = numThreadIndex + 7;
        c[8] = numThreadIndex + 8; c[9] = numThreadIndex + 9;
        c[10] = numThreadIndex + 10; c[11] = numThreadIndex + 11;
        c[12] = numThreadIndex + 12; c[13] = numThreadIndex + 13;
        c[14] = numThreadIndex + 14; c[15] = numThreadIndex + 15;
        c[16] = numThreadIndex + 16; c[17] = numThreadIndex + 17;
        c[18] = numThreadIndex + 18; c[19] = numThreadIndex + 19;
        c[20] = numThreadIndex + 20; c[21] = numThreadIndex + 21;
        c[22] = numThreadIndex + 22; c[23] = numThreadIndex + 23;
        c[24] = numThreadIndex + 24; c[25] = numThreadIndex + 25;
        c[26] = numThreadIndex + 26; c[27] = numThreadIndex + 27;
        c[28] = numThreadIndex + 28; c[29] = numThreadIndex + 29;
        c[30] = numThreadIndex + 30; c[31] = numThreadIndex + 31;
        c[32] = numThreadIndex + 32; c[33] = numThreadIndex + 33;
        c[34] = numThreadIndex + 34; c[35] = numThreadIndex + 35;
        c += 36;
        numThreadIndex += 36;

        m[0] = sparseMSIds[0]; m[1] = sparseMSIds[0];
        m[2] = sparseMSIds[1]; m[3] = sparseMSIds[1];
        m[4] = sparseMSIds[2]; m[5] = sparseMSIds[2];
        m[6] = sparseMSIds[3]; m[7] = sparseMSIds[3];
        m[8] = sparseMSIds[4]; m[9] = sparseMSIds[4];
        m[10] = sparseMSIds[5]; m[11] = sparseMSIds[5];
        m += 12;
      } else { // 2 or 3 labels on tetraeder
        const size_t numTris = tetLookupNumWallTriangles[lookupIndex];
        long long sparseIds[numTris];
        for(size_t t = 0; t < numTris; ++t) {
          sparseIds[t] = msLabelMap[getSparseId(asc[tetVertLabel[t * 2]],
                                     asc[tetVertLabel[(t * 2) + 1]],
                                     numAscRegions)];

          p[0] = eC[tetEdgeIndices[(t * 3)    ]][0];
          p[1] = eC[tetEdgeIndices[(t * 3)    ]][1];
          p[2] = eC[tetEdgeIndices[(t * 3)    ]][2];
          p[3] = eC[tetEdgeIndices[(t * 3) + 1]][0];
          p[4] = eC[tetEdgeIndices[(t * 3) + 1]][1];
          p[5] = eC[tetEdgeIndices[(t * 3) + 1]][2];
          p[6] = eC[tetEdgeIndices[(t * 3) + 2]][0];
          p[7] = eC[tetEdgeIndices[(t * 3) + 2]][1];
          p[8] = eC[tetEdgeIndices[(t * 3) + 2]][2];
          p += 9;

          c[0] = numThreadIndex + 0;
          c[1] = numThreadIndex + 1;
          c[2] = numThreadIndex + 2;
          c += 3;
          numThreadIndex += 3;

          m[0] = sparseIds[t];
          m += 1;
        }
      }
    }

    for (const auto& vCase : validCasesD)
    {
      const SimplexId& tet = vCase.first;
      const unsigned char& lookupIndex = vCase.second;
      
      const int *tetEdgeIndices = tetLookupWall[lookupIndex];
      const int *tetVertLabel = tetLookupWallLabel[lookupIndex];

      SimplexId vertices[4];
      triangulation.getCellVertex(tet, 0, vertices[0]);
      triangulation.getCellVertex(tet, 1, vertices[1]);
      triangulation.getCellVertex(tet, 2, vertices[2]);
      triangulation.getCellVertex(tet, 3, vertices[3]);

      const SimplexId dsc[4] = {
        dscManifold[vertices[0]], dscManifold[vertices[1]],
        dscManifold[vertices[2]], dscManifold[vertices[3]]};

      float vertPos[4][3];
      triangulation.getVertexPoint(
        vertices[0], vertPos[0][0], vertPos[0][1], vertPos[0][2]);
      triangulation.getVertexPoint(
        vertices[1], vertPos[1][0], vertPos[1][1], vertPos[1][2]);
      triangulation.getVertexPoint(
        vertices[2], vertPos[2][0], vertPos[2][1], vertPos[2][2]);
      triangulation.getVertexPoint(
        vertices[3], vertPos[3][0], vertPos[3][1], vertPos[3][2]);

      float eC[10][3];
      // 6 edge centers
      getCenter(vertPos[0], vertPos[1], eC[0]);
      getCenter(vertPos[0], vertPos[2], eC[1]);
      getCenter(vertPos[0], vertPos[3], eC[2]);
      getCenter(vertPos[1], vertPos[2], eC[3]);
      getCenter(vertPos[1], vertPos[3], eC[4]);
      getCenter(vertPos[2], vertPos[3], eC[5]);

      // 4 triangle centers
      getCenter(vertPos[0], vertPos[1], vertPos[2], eC[6]);
      getCenter(vertPos[0], vertPos[1], vertPos[3], eC[7]);
      getCenter(vertPos[0], vertPos[2], vertPos[3], eC[8]);
      getCenter(vertPos[1], vertPos[2], vertPos[3], eC[9]);

      if(tetEdgeIndices[0] == 10) { // 4 labels on tetraeder
        float tetCenter[3];
        getCenter(vertPos[0], vertPos[1], vertPos[2], vertPos[3], tetCenter);

        long long sparseMSIds[6] = {
          msLabelMap[getSparseId(dsc[0], dsc[1], numDscRegions) 
          + ascSparseIDOffset],
          msLabelMap[getSparseId(dsc[0], dsc[2], numDscRegions) 
          + ascSparseIDOffset],
          msLabelMap[getSparseId(dsc[0], dsc[3], numDscRegions) 
          + ascSparseIDOffset],
          msLabelMap[getSparseId(dsc[1], dsc[2], numDscRegions) 
          + ascSparseIDOffset],
          msLabelMap[getSparseId(dsc[1], dsc[3], numDscRegions) 
          + ascSparseIDOffset],
          msLabelMap[getSparseId(dsc[2], dsc[3], numDscRegions) 
          + ascSparseIDOffset]
        };

        p[0] = eC[7][0]; p[1] = eC[7][1]; p[2] = eC[7][2]; 
        p[3] = eC[0][0]; p[4] = eC[0][1]; p[5] = eC[0][2]; 
        p[6] = tetCenter[0]; p[7] = tetCenter[1]; p[8] = tetCenter[2];
        p[9] = eC[0][0]; p[10] = eC[0][1]; p[11] = eC[0][2]; 
        p[12] = eC[6][0]; p[13] = eC[6][1]; p[14] = eC[6][2]; 
        p[15] = tetCenter[0]; p[16] = tetCenter[1]; p[17] = tetCenter[2];
        p[18] = eC[8][0]; p[19] = eC[8][1]; p[20] = eC[8][2];
        p[21] = eC[1][0]; p[22] = eC[1][1]; p[23] = eC[1][2];
        p[24] = tetCenter[0]; p[25] = tetCenter[1]; p[26] = tetCenter[2];
        p[27] = eC[1][0]; p[28] = eC[1][1]; p[29] = eC[1][2];
        p[30] = eC[6][0]; p[31] = eC[6][1]; p[32] = eC[6][2];
        p[33] = tetCenter[0]; p[34] = tetCenter[1]; p[35] = tetCenter[2];
        p[36] = eC[8][0]; p[37] = eC[8][1]; p[38] = eC[8][2]; 
        p[39] = eC[2][0]; p[40] = eC[2][1]; p[41] = eC[2][2]; 
        p[42] = tetCenter[0]; p[43] = tetCenter[1]; p[44] = tetCenter[2];
        p[45] = eC[2][0]; p[46] = eC[2][1]; p[47] = eC[2][2]; 
        p[48] = eC[7][0]; p[49] = eC[7][1]; p[50] = eC[7][2]; 
        p[51] = tetCenter[0]; p[52] = tetCenter[1]; p[53] = tetCenter[2];
        p[54] = eC[6][0]; p[55] = eC[6][1]; p[56] = eC[6][2];
        p[57] = eC[3][0]; p[58] = eC[3][1]; p[59] = eC[3][2];
        p[60] = tetCenter[0]; p[61] = tetCenter[1]; p[62] = tetCenter[2];
        p[63] = eC[3][0]; p[64] = eC[3][1]; p[65] = eC[3][2];
        p[66] = eC[9][0]; p[67] = eC[9][1]; p[68] = eC[9][2];
        p[69] = tetCenter[0]; p[70] = tetCenter[1]; p[71] = tetCenter[2];
        p[72] = eC[7][0]; p[73] = eC[7][1]; p[74] = eC[7][2];
        p[75] = eC[4][0]; p[76] = eC[4][1]; p[77] = eC[4][2];
        p[78] = tetCenter[0]; p[79] = tetCenter[1]; p[80] = tetCenter[2];
        p[81] = eC[4][0]; p[82] = eC[4][1]; p[83] = eC[4][2];
        p[84] = eC[9][0]; p[85] = eC[9][1]; p[86] = eC[9][2];
        p[87] = tetCenter[0]; p[88] = tetCenter[1]; p[89] = tetCenter[2];
        p[90] = eC[9][0]; p[91] = eC[9][1]; p[92] = eC[9][2]; 
        p[93] = eC[5][0]; p[94] = eC[5][1]; p[95] = eC[5][2]; 
        p[96] = tetCenter[0]; p[97] = tetCenter[1]; p[98] = tetCenter[2];
        p[99] = eC[5][0]; p[100] = eC[5][1]; p[101] = eC[5][2];
        p[102] = eC[8][0]; p[103] = eC[8][1]; p[104] = eC[8][2];
        p[105] = tetCenter[0]; p[106] = tetCenter[1]; p[107] = tetCenter[2];
        p += 108;

        c[0] = numThreadIndex + 0; c[1] = numThreadIndex + 1;
        c[2] = numThreadIndex + 2; c[3] = numThreadIndex + 3;
        c[4] = numThreadIndex + 4; c[5] = numThreadIndex + 5;
        c[6] = numThreadIndex + 6; c[7] = numThreadIndex + 7;
        c[8] = numThreadIndex + 8; c[9] = numThreadIndex + 9;
        c[10] = numThreadIndex + 10; c[11] = numThreadIndex + 11;
        c[12] = numThreadIndex + 12; c[13] = numThreadIndex + 13;
        c[14] = numThreadIndex + 14; c[15] = numThreadIndex + 15;
        c[16] = numThreadIndex + 16; c[17] = numThreadIndex + 17;
        c[18] = numThreadIndex + 18; c[19] = numThreadIndex + 19;
        c[20] = numThreadIndex + 20; c[21] = numThreadIndex + 21;
        c[22] = numThreadIndex + 22; c[23] = numThreadIndex + 23;
        c[24] = numThreadIndex + 24; c[25] = numThreadIndex + 25;
        c[26] = numThreadIndex + 26; c[27] = numThreadIndex + 27;
        c[28] = numThreadIndex + 28; c[29] = numThreadIndex + 29;
        c[30] = numThreadIndex + 30; c[31] = numThreadIndex + 31;
        c[32] = numThreadIndex + 32; c[33] = numThreadIndex + 33;
        c[34] = numThreadIndex + 34; c[35] = numThreadIndex + 35;
        c += 36;
        numThreadIndex += 36;

        m[0] = sparseMSIds[0]; m[1] = sparseMSIds[0];
        m[2] = sparseMSIds[1]; m[3] = sparseMSIds[1];
        m[4] = sparseMSIds[2]; m[5] = sparseMSIds[2];
        m[6] = sparseMSIds[3]; m[7] = sparseMSIds[3];
        m[8] = sparseMSIds[4]; m[9] = sparseMSIds[4];
        m[10] = sparseMSIds[5]; m[11] = sparseMSIds[5];
        m += 12;
        
      } else { // 2 or 3 labels on tetraeder
        const size_t numTris = tetLookupNumWallTriangles[lookupIndex];
        long long sparseIds[numTris];
        for(size_t t = 0; t < numTris; ++t) {
          sparseIds[t] = msLabelMap[getSparseId(dsc[tetVertLabel[t * 2]],
                                     dsc[tetVertLabel[(t * 2) + 1]],
                                     numDscRegions) + ascSparseIDOffset];

          p[0] = eC[tetEdgeIndices[(t * 3)    ]][0];
          p[1] = eC[tetEdgeIndices[(t * 3)    ]][1];
          p[2] = eC[tetEdgeIndices[(t * 3)    ]][2];
          p[3] = eC[tetEdgeIndices[(t * 3) + 1]][0];
          p[4] = eC[tetEdgeIndices[(t * 3) + 1]][1];
          p[5] = eC[tetEdgeIndices[(t * 3) + 1]][2];
          p[6] = eC[tetEdgeIndices[(t * 3) + 2]][0];
          p[7] = eC[tetEdgeIndices[(t * 3) + 2]][1];
          p[8] = eC[tetEdgeIndices[(t * 3) + 2]][2];
          p += 9;

          c[0] = numThreadIndex + 0;
          c[1] = numThreadIndex + 1;
          c[2] = numThreadIndex + 2;
          c += 3;
          numThreadIndex += 3;

          m[0] = sparseIds[t];
          m += 1;
        }
      }
    }
  }
}
#endif // TTK_ENABLE_OPENMP

  this->printMsg("Computed 2-Separatrices 3D[Split]", 1,
    localTimer.getElapsedTime(), this->threadNumber_);

  return 0;
}

template <typename triangulationType>
int ttk::MorseSmaleSegmentationPL::computeBasinSeparation_3D_fine(
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

#ifdef TTK_ENABLE_OPENMP
if(threadNumber_ == 1) {
#endif // TTK_ENABLE_OPENMP

  std::vector<std::pair<SimplexId, unsigned char>> validCases;
  SimplexId numTriangles = 0;

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

    if(tetLookupIsMultiLabel[lookupIndex]) {
      validCases.push_back(std::make_pair(tet, lookupIndex));
      numTriangles += tetLookupNumFineTriangles[lookupIndex];
    }
  }

  outputSeparatrices2_points_->resize(9 * numTriangles);
  outputSeparatrices2_cells_connectivity_->resize(3 * numTriangles);
  outputSeparatrices2_cells_mscIds_->resize(numTriangles);
  *outputSeparatrices2_numberOfPoints_ = 3 * numTriangles;
  *outputSeparatrices2_numberOfCells_ = numTriangles;

  auto &points = *outputSeparatrices2_points_;
  auto &cellsConn = *outputSeparatrices2_cells_connectivity_;
  auto &cellsMSCIds = *outputSeparatrices2_cells_mscIds_;

  float* p = points.data();
  SimplexId* c = cellsConn.data();
  SimplexId* m = cellsMSCIds.data();
  SimplexId cellIndex = 0;

  for (const auto& vCase : validCases)
  {
    const SimplexId& tet = vCase.first;
    const unsigned char& lookupIndex = vCase.second;

    SimplexId vertices[4];
    triangulation.getCellVertex(tet, 0, vertices[0]);
    triangulation.getCellVertex(tet, 1, vertices[1]);
    triangulation.getCellVertex(tet, 2, vertices[2]);
    triangulation.getCellVertex(tet, 3, vertices[3]);

    const SimplexId msm[4] = {
      morseSmaleManifold[vertices[0]], morseSmaleManifold[vertices[1]],
      morseSmaleManifold[vertices[2]], morseSmaleManifold[vertices[3]]};
    
    float vPos[4][3];
    triangulation.getVertexPoint(
      vertices[0], vPos[0][0], vPos[0][1], vPos[0][2]);
    triangulation.getVertexPoint(
      vertices[1], vPos[1][0], vPos[1][1], vPos[1][2]);
    triangulation.getVertexPoint(
      vertices[2], vPos[2][0], vPos[2][1], vPos[2][2]);
    triangulation.getVertexPoint(
      vertices[3], vPos[3][0], vPos[3][1], vPos[3][2]);

    if(tetLookupIs2Label[lookupIndex]) { // 2 labels (eg. AAAB / AABB)
      const int *vIds = tetLookupSplitBasins2Label[lookupIndex];

      float vert00[3], vert01[3], vert02[3],
            vert10[3], vert11[3], vert12[3];

      interpolatePoints(vPos[vIds[0]], vPos[vIds[1]], d0, vert00);
      interpolatePoints(vPos[vIds[2]], vPos[vIds[3]], d0, vert01);
      interpolatePoints(vPos[vIds[4]], vPos[vIds[5]], d0, vert02);
      interpolatePoints(vPos[vIds[0]], vPos[vIds[1]], d1, vert10);
      interpolatePoints(vPos[vIds[2]], vPos[vIds[3]], d1, vert11);
      interpolatePoints(vPos[vIds[4]], vPos[vIds[5]], d1, vert12);

      p[0] = vert00[0]; p[1] = vert00[1]; p[2] = vert00[2]; 
      p[3] = vert01[0]; p[4] = vert01[1]; p[5] = vert01[2]; 
      p[6] = vert02[0]; p[7] = vert02[1]; p[8] = vert02[2];

      p[9] = vert10[0]; p[10] = vert10[1]; p[11] = vert10[2]; 
      p[12] = vert11[0]; p[13] = vert11[1]; p[14] = vert11[2]; 
      p[15] = vert12[0]; p[16] = vert12[1]; p[17] = vert12[2];
      p += 18;

      c[0] = cellIndex + 0; c[1] = cellIndex + 1;
      c[2] = cellIndex + 2; c[3] = cellIndex + 3;
      c[4] = cellIndex + 4; c[5] = cellIndex + 5;
      c += 6;
      cellIndex += 6;

      m[0] = msm[vIds[0]]; m[1] = msm[vIds[1]];
      m += 2;

      if(vIds[6] != -1) { // 2 vertices per label (e.g. AABB)
        float vert03[3], vert13[3];

        interpolatePoints(vPos[vIds[6]], vPos[vIds[7]], d0, vert03);
        interpolatePoints(vPos[vIds[6]], vPos[vIds[7]], d1, vert13);

        p[0] = vert00[0]; p[1] = vert00[1]; p[2] = vert00[2]; 
        p[3] = vert01[0]; p[4] = vert01[1]; p[5] = vert01[2]; 
        p[6] = vert03[0]; p[7] = vert03[1]; p[8] = vert03[2];
        p[9] = vert10[0]; p[10] = vert10[1]; p[11] = vert10[2]; 
        p[12] = vert11[0]; p[13] = vert11[1]; p[14] = vert11[2]; 
        p[15] = vert13[0]; p[16] = vert13[1]; p[17] = vert13[2];
        p += 18;

        c[0] = cellIndex + 0; c[1] = cellIndex + 1;
        c[2] = cellIndex + 2; c[3] = cellIndex + 3;
        c[4] = cellIndex + 4; c[5] = cellIndex + 5;
        c += 6;
        cellIndex += 6;

        m[0] = msm[vIds[0]]; m[1] = msm[vIds[1]];
        m += 2;
      }
    } else if(tetLookupIs3Label[lookupIndex]) {
      const int *vIds = tetLookupSplitBasisns3Label[lookupIndex];

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

      // Label 0
      p[0] = edge00[0]; p[1] = edge00[1]; p[2] = edge00[2];
      p[3] = edge02[0]; p[4] = edge02[1]; p[5] = edge02[2];
      p[6] = tri00[0]; p[7] = tri00[1]; p[8] = tri00[2];
      p[9] = edge02[0]; p[10] = edge02[1]; p[11] = edge02[2];
      p[12] = tri00[0]; p[13] = tri00[1]; p[14] = tri00[2];
      p[15] = tri01[0]; p[16] = tri01[1]; p[17] = tri01[2];        
      p[18] = edge01[0]; p[19] = edge01[1]; p[20] = edge01[2];
      p[21] = edge03[0]; p[22] = edge03[1]; p[23] = edge03[2];
      p[24] = tri00[0]; p[25] = tri00[1]; p[26] = tri00[2];
      p[27] = edge03[0]; p[28] = edge03[1]; p[29] = edge03[2];
      p[30] = tri00[0]; p[31] = tri00[1]; p[32] = tri00[2];
      p[33] = tri01[0]; p[34] = tri01[1]; p[35] = tri01[2];

      // Label 1           
      p[36] = edge10[0]; p[37] = edge10[1]; p[38] = edge10[2];
      p[39] = edge11[0]; p[40] = edge11[1]; p[41] = edge11[2];
      p[42] = tri10[0]; p[43] = tri10[1]; p[44] = tri10[2];
      p[45] = edge11[0]; p[46] = edge11[1]; p[47] = edge11[2];
      p[48] = tri10[0]; p[49] = tri10[1]; p[50] = tri10[2];
      p[51] = tri11[0]; p[52] = tri11[1]; p[53] = tri11[2];
      p[54] = edge12[0]; p[55] = edge12[1]; p[56] = edge12[2];
      p[57] = tri10[0]; p[58] = tri10[1]; p[59] = tri10[2];
      p[60] = tri11[0]; p[61] = tri11[1]; p[62] = tri11[2];

      // Label 2           
      p[63] = edge20[0]; p[64] = edge20[1]; p[65] = edge20[2];
      p[66] = edge21[0]; p[67] = edge21[1]; p[68] = edge21[2]; 
      p[69] = tri20[0]; p[70] = tri20[1]; p[71] = tri20[2];
      p[72] = edge21[0]; p[73] = edge21[1]; p[74] = edge21[2];
      p[75] = tri20[0]; p[76] = tri20[1]; p[77] = tri20[2];
      p[78] = tri21[0]; p[79] = tri21[1]; p[80] = tri21[2];
      p[81] = edge22[0]; p[82] = edge22[1]; p[83] = edge22[2]; 
      p[84] = tri20[0]; p[85] = tri20[1]; p[86] = tri20[2];
      p[87] = tri21[0]; p[88] = tri21[1]; p[89] = tri21[2];
      p += 90;

      c[0] = cellIndex + 0; c[1] = cellIndex + 1;
      c[2] = cellIndex + 2; c[3] = cellIndex + 3;
      c[4] = cellIndex + 4; c[5] = cellIndex + 5;
      c[6] = cellIndex + 6; c[7] = cellIndex + 7;
      c[8] = cellIndex + 8; c[9] = cellIndex + 9;
      c[10] = cellIndex + 10; c[11] = cellIndex + 11;
      c[12] = cellIndex + 12; c[13] = cellIndex + 13;
      c[14] = cellIndex + 14; c[15] = cellIndex + 15;
      c[16] = cellIndex + 16; c[17] = cellIndex + 17;
      c[18] = cellIndex + 18; c[19] = cellIndex + 19;
      c[20] = cellIndex + 20; c[21] = cellIndex + 21;
      c[22] = cellIndex + 22; c[23] = cellIndex + 23;
      c[24] = cellIndex + 24; c[25] = cellIndex + 25;
      c[26] = cellIndex + 26; c[27] = cellIndex + 27;
      c[28] = cellIndex + 28; c[29] = cellIndex + 29;
      c += 30;
      cellIndex += 30;

      m[0] = msm[vIds[0]]; m[1] = msm[vIds[0]];        
      m[2] = msm[vIds[0]]; m[3] = msm[vIds[0]];
      m[4] = msm[vIds[8]]; m[5] = msm[vIds[8]];        
      m[6] = msm[vIds[8]]; m[7] = msm[vIds[9]];       
      m[8] = msm[vIds[9]]; m[9] = msm[vIds[9]];
      m += 10;
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
      p[0] = vert00[0];   p[1] = vert00[1];   p[2] = vert00[2]; 
      p[3] = vert0t0[0];  p[4] = vert0t0[1];  p[5] = vert0t0[2]; 
      p[6] = vert0tet[0]; p[7] = vert0tet[1]; p[8] = vert0tet[2];
      p[9] = vert00[0];   p[10] = vert00[1];   p[11] = vert00[2]; 
      p[12] = vert0t1[0];  p[13] = vert0t1[1];  p[14] = vert0t1[2]; 
      p[15] = vert0tet[0]; p[16] = vert0tet[1]; p[17] = vert0tet[2];
      p[18] = vert01[0];   p[19] = vert01[1];   p[20] = vert01[2]; 
      p[21] = vert0t0[0];  p[22] = vert0t0[1];  p[23] = vert0t0[2]; 
      p[24] = vert0tet[0]; p[25] = vert0tet[1]; p[26] = vert0tet[2];
      p[27] = vert01[0];   p[28] = vert01[1];   p[29] = vert01[2]; 
      p[30] = vert0t2[0];  p[31] = vert0t2[1];  p[32] = vert0t2[2]; 
      p[33] = vert0tet[0]; p[34] = vert0tet[1]; p[35] = vert0tet[2];
      p[36] = vert02[0];   p[37] = vert02[1];   p[38] = vert02[2]; 
      p[39] = vert0t2[0];  p[40] = vert0t2[1];  p[41] = vert0t2[2]; 
      p[42] = vert0tet[0]; p[43] = vert0tet[1]; p[44] = vert0tet[2];
      p[45] = vert02[0];   p[46] = vert02[1];   p[47] = vert02[2]; 
      p[48] = vert0t1[0];  p[49] = vert0t1[1];  p[50] = vert0t1[2]; 
      p[51] = vert0tet[0]; p[52] = vert0tet[1]; p[53] = vert0tet[2];

      // Label Vert 1       
      p[54] = vert10[0];   p[55] = vert10[1];   p[56] = vert10[2]; 
      p[57] = vert1t0[0];  p[58] = vert1t0[1];  p[59] = vert1t0[2]; 
      p[60] = vert1tet[0]; p[61] = vert1tet[1]; p[62] = vert1tet[2];
      p[63] = vert10[0];   p[64] = vert10[1];   p[65] = vert10[2]; 
      p[66] = vert1t1[0];  p[67] = vert1t1[1];  p[68] = vert1t1[2]; 
      p[69] = vert1tet[0]; p[70] = vert1tet[1]; p[71] = vert1tet[2];
      p[72] = vert11[0];   p[73] = vert11[1];   p[74] = vert11[2]; 
      p[75] = vert1t0[0];  p[76] = vert1t0[1];  p[77] = vert1t0[2]; 
      p[78] = vert1tet[0]; p[79] = vert1tet[1]; p[80] = vert1tet[2];
      p[81] = vert11[0];   p[82] = vert11[1];   p[83] = vert11[2]; 
      p[84] = vert1t2[0];  p[85] = vert1t2[1];  p[86] = vert1t2[2]; 
      p[87] = vert1tet[0]; p[88] = vert1tet[1]; p[89] = vert1tet[2];
      p[90] = vert12[0];   p[91] = vert12[1];   p[92] = vert12[2]; 
      p[93] = vert1t2[0];  p[94] = vert1t2[1];  p[95] = vert1t2[2]; 
      p[96] = vert1tet[0]; p[97] = vert1tet[1]; p[98] = vert1tet[2];
      p[99] = vert12[0];   p[100] = vert12[1];   p[101] = vert12[2]; 
      p[102] = vert1t1[0];  p[103] = vert1t1[1];  p[104] = vert1t1[2]; 
      p[105] = vert1tet[0]; p[106] = vert1tet[1]; p[107] = vert1tet[2];


      // Label Vert 2     
      p[108] = vert20[0];   p[109] = vert20[1];   p[110] = vert20[2]; 
      p[111] = vert2t0[0];  p[112] = vert2t0[1];  p[113] = vert2t0[2]; 
      p[114] = vert2tet[0]; p[115] = vert2tet[1]; p[116] = vert2tet[2];
      p[117] = vert20[0];   p[118] = vert20[1];   p[119] = vert20[2]; 
      p[120] = vert2t1[0];  p[121] = vert2t1[1];  p[122] = vert2t1[2]; 
      p[123] = vert2tet[0]; p[124] = vert2tet[1]; p[125] = vert2tet[2];
      p[126] = vert21[0];   p[127] = vert21[1];   p[128] = vert21[2]; 
      p[129] = vert2t0[0];  p[130] = vert2t0[1];  p[131] = vert2t0[2]; 
      p[132] = vert2tet[0]; p[133] = vert2tet[1]; p[134] = vert2tet[2];
      p[135] = vert21[0];   p[136] = vert21[1];   p[137] = vert21[2]; 
      p[138] = vert2t2[0];  p[139] = vert2t2[1];  p[140] = vert2t2[2]; 
      p[141] = vert2tet[0]; p[142] = vert2tet[1]; p[143] = vert2tet[2];
      p[144] = vert22[0];   p[145] = vert22[1];   p[146] = vert22[2]; 
      p[147] = vert2t2[0];  p[148] = vert2t2[1];  p[149] = vert2t2[2]; 
      p[150] = vert2tet[0]; p[151] = vert2tet[1]; p[152] = vert2tet[2];
      p[153] = vert22[0];   p[154] = vert22[1];   p[155] = vert22[2]; 
      p[156] = vert2t1[0];  p[157] = vert2t1[1];  p[158] = vert2t1[2]; 
      p[159] = vert2tet[0]; p[160] = vert2tet[1]; p[161] = vert2tet[2];

      // Label Vert 3      
      p[162] = vert30[0];   p[163] = vert30[1];   p[164] = vert30[2]; 
      p[165] = vert3t0[0];  p[166] = vert3t0[1];  p[167] = vert3t0[2]; 
      p[168] = vert3tet[0]; p[169] = vert3tet[1]; p[170] = vert3tet[2];            
      p[171] = vert30[0];   p[172] = vert30[1];   p[173] = vert30[2]; 
      p[174] = vert3t1[0];  p[175] = vert3t1[1];  p[176] = vert3t1[2]; 
      p[177] = vert3tet[0]; p[178] = vert3tet[1]; p[179] = vert3tet[2];
      p[180] = vert31[0];   p[181] = vert31[1];   p[182] = vert31[2]; 
      p[183] = vert3t0[0];  p[184] = vert3t0[1];  p[185] = vert3t0[2]; 
      p[186] = vert3tet[0]; p[187] = vert3tet[1]; p[188] = vert3tet[2];
      p[189] = vert31[0];   p[190] = vert31[1];   p[191] = vert31[2]; 
      p[192] = vert3t2[0];  p[193] = vert3t2[1];  p[194] = vert3t2[2]; 
      p[195] = vert3tet[0]; p[196] = vert3tet[1]; p[197] = vert3tet[2];
      p[198] = vert32[0];   p[199] = vert32[1];   p[200] = vert32[2]; 
      p[201] = vert3t2[0];  p[202] = vert3t2[1];  p[203] = vert3t2[2]; 
      p[204] = vert3tet[0]; p[205] = vert3tet[1]; p[206] = vert3tet[2];
      p[207] = vert32[0];   p[208] = vert32[1];   p[209] = vert32[2]; 
      p[210] = vert3t1[0];  p[211] = vert3t1[1];  p[212] = vert3t1[2]; 
      p[213] = vert3tet[0]; p[214] = vert3tet[1]; p[215] = vert3tet[2];
      p += 216;

      c[0] = cellIndex + 0; c[1] = cellIndex + 1;
      c[2] = cellIndex + 2; c[3] = cellIndex + 3;
      c[4] = cellIndex + 4; c[5] = cellIndex + 5;
      c[6] = cellIndex + 6; c[7] = cellIndex + 7;
      c[8] = cellIndex + 8; c[9] = cellIndex + 9;
      c[10] = cellIndex + 10; c[11] = cellIndex + 11;
      c[12] = cellIndex + 12; c[13] = cellIndex + 13;
      c[14] = cellIndex + 14; c[15] = cellIndex + 15;
      c[16] = cellIndex + 16; c[17] = cellIndex + 17;
      c[18] = cellIndex + 18; c[19] = cellIndex + 19;
      c[20] = cellIndex + 20; c[21] = cellIndex + 21;
      c[22] = cellIndex + 22; c[23] = cellIndex + 23;
      c[24] = cellIndex + 24; c[25] = cellIndex + 25;
      c[26] = cellIndex + 26; c[27] = cellIndex + 27;
      c[28] = cellIndex + 28; c[29] = cellIndex + 29;
      c[30] = cellIndex + 30; c[31] = cellIndex + 31;
      c[32] = cellIndex + 32; c[33] = cellIndex + 33;
      c[34] = cellIndex + 34; c[35] = cellIndex + 35;
      c[36] = cellIndex + 36; c[37] = cellIndex + 37;
      c[38] = cellIndex + 38; c[39] = cellIndex + 39;
      c[40] = cellIndex + 40; c[41] = cellIndex + 41;
      c[42] = cellIndex + 42; c[43] = cellIndex + 43;
      c[44] = cellIndex + 44; c[45] = cellIndex + 45;
      c[46] = cellIndex + 46; c[47] = cellIndex + 47;
      c[48] = cellIndex + 48; c[49] = cellIndex + 49;
      c[50] = cellIndex + 50; c[51] = cellIndex + 51;
      c[52] = cellIndex + 52; c[53] = cellIndex + 53;
      c[54] = cellIndex + 54; c[55] = cellIndex + 55;
      c[56] = cellIndex + 56; c[57] = cellIndex + 57;
      c[58] = cellIndex + 58; c[59] = cellIndex + 59;
      c[60] = cellIndex + 60; c[61] = cellIndex + 61;
      c[62] = cellIndex + 62; c[63] = cellIndex + 63;
      c[64] = cellIndex + 64; c[65] = cellIndex + 65;
      c[66] = cellIndex + 66; c[67] = cellIndex + 67;
      c[68] = cellIndex + 68; c[69] = cellIndex + 69;
      c[70] = cellIndex + 70; c[71] = cellIndex + 71;
      c += 72;
      cellIndex += 72;

      m[0] = msm[0]; m[1] = msm[0];
      m[2] = msm[0]; m[3] = msm[0];
      m[4] = msm[0]; m[5] = msm[0];
      m[6] = msm[1]; m[7] = msm[1];
      m[8] = msm[1]; m[9] = msm[1];
      m[10] = msm[1]; m[11] = msm[1];
      m[12] = msm[2]; m[13] = msm[2];
      m[14] = msm[2]; m[15] = msm[2];
      m[16] = msm[2]; m[17] = msm[2];
      m[18] = msm[3]; m[19] = msm[3];
      m[20] = msm[3]; m[21] = msm[3];
      m[22] = msm[3]; m[23] = msm[3];
      m += 24;
    }
  }

#ifdef TTK_ENABLE_OPENMP
} else {
  size_t triangleStartIndex[threadNumber_ + 1];

  #pragma omp parallel num_threads(threadNumber_)
  {
    const int tid = omp_get_thread_num();
    std::vector<std::pair<SimplexId, unsigned char>> validCases;
    size_t numThreadTriangles = 0;
    size_t numThreadIndex = 0;

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

      if(tetLookupIsMultiLabel[lookupIndex]) {
        validCases.push_back(std::make_pair(tet, lookupIndex));
        numThreadTriangles += tetLookupNumFineTriangles[lookupIndex];
      }
    }

    triangleStartIndex[tid + 1] = numThreadTriangles;

    #pragma omp barrier

    #pragma omp single
    {
      triangleStartIndex[0] = 0;

      // Count triangle number and create iterator start indices
      for(int t = 1; t < threadNumber_ + 1; ++t) {
        triangleStartIndex[t] += triangleStartIndex[t-1];
      }
    }

    const size_t numTriangles = triangleStartIndex[threadNumber_];
    numThreadIndex = triangleStartIndex[tid];

    #pragma omp single nowait
    outputSeparatrices2_points_->resize(9 * numTriangles);

    #pragma omp single nowait
    outputSeparatrices2_cells_connectivity_->resize(3 * numTriangles);

    #pragma omp single nowait
    outputSeparatrices2_cells_mscIds_->resize(numTriangles);

    #pragma omp single
    {
      *outputSeparatrices2_numberOfPoints_ = 3 * numTriangles;
      *outputSeparatrices2_numberOfCells_ = numTriangles;
    }
  
    auto &points = *outputSeparatrices2_points_;
    auto &cellsConn = *outputSeparatrices2_cells_connectivity_;
    auto &cellsMSCIds = *outputSeparatrices2_cells_mscIds_;

    float* p = points.data();
    SimplexId* c = cellsConn.data();
    SimplexId* m = cellsMSCIds.data();

    p += (numThreadIndex * 9);
    c += (numThreadIndex * 3);
    m += numThreadIndex;

    numThreadIndex = 3 * numThreadIndex;

    for (const auto& vCase : validCases)
    {
      const SimplexId& tet = vCase.first;
      const unsigned char& lookupIndex = vCase.second;

      SimplexId vertices[4];
      triangulation.getCellVertex(tet, 0, vertices[0]);
      triangulation.getCellVertex(tet, 1, vertices[1]);
      triangulation.getCellVertex(tet, 2, vertices[2]);
      triangulation.getCellVertex(tet, 3, vertices[3]);

      const SimplexId msm[4] = {
        morseSmaleManifold[vertices[0]], morseSmaleManifold[vertices[1]],
        morseSmaleManifold[vertices[2]], morseSmaleManifold[vertices[3]]};

      float vPos[4][3];
      triangulation.getVertexPoint(
        vertices[0], vPos[0][0], vPos[0][1], vPos[0][2]);
      triangulation.getVertexPoint(
        vertices[1], vPos[1][0], vPos[1][1], vPos[1][2]);
      triangulation.getVertexPoint(
        vertices[2], vPos[2][0], vPos[2][1], vPos[2][2]);
      triangulation.getVertexPoint(
        vertices[3], vPos[3][0], vPos[3][1], vPos[3][2]);

      if(tetLookupIs2Label[lookupIndex]) { // 2 labels (eg. AAAB / AABB)
        const int *vIds = tetLookupSplitBasins2Label[lookupIndex];

        float vert00[3], vert01[3], vert02[3],
              vert10[3], vert11[3], vert12[3];

        interpolatePoints(vPos[vIds[0]], vPos[vIds[1]], d0, vert00);
        interpolatePoints(vPos[vIds[2]], vPos[vIds[3]], d0, vert01);
        interpolatePoints(vPos[vIds[4]], vPos[vIds[5]], d0, vert02);
        interpolatePoints(vPos[vIds[0]], vPos[vIds[1]], d1, vert10);
        interpolatePoints(vPos[vIds[2]], vPos[vIds[3]], d1, vert11);
        interpolatePoints(vPos[vIds[4]], vPos[vIds[5]], d1, vert12);

        p[0] = vert00[0]; p[1] = vert00[1]; p[2] = vert00[2]; 
        p[3] = vert01[0]; p[4] = vert01[1]; p[5] = vert01[2]; 
        p[6] = vert02[0]; p[7] = vert02[1]; p[8] = vert02[2];

        p[9] = vert10[0]; p[10] = vert10[1]; p[11] = vert10[2]; 
        p[12] = vert11[0]; p[13] = vert11[1]; p[14] = vert11[2]; 
        p[15] = vert12[0]; p[16] = vert12[1]; p[17] = vert12[2];
        p += 18;

        c[0] = numThreadIndex + 0; c[1] = numThreadIndex + 1;
        c[2] = numThreadIndex + 2; c[3] = numThreadIndex + 3;
        c[4] = numThreadIndex + 4; c[5] = numThreadIndex + 5;
        c += 6;
        numThreadIndex += 6;

        m[0] = msm[vIds[0]]; m[1] = msm[vIds[1]];
        m += 2;

        if(vIds[6] != -1) { // 2 vertices per label (e.g. AABB)
          float vert03[3], vert13[3];

          interpolatePoints(vPos[vIds[6]], vPos[vIds[7]], d0, vert03);
          interpolatePoints(vPos[vIds[6]], vPos[vIds[7]], d1, vert13);

          p[0] = vert00[0]; p[1] = vert00[1]; p[2] = vert00[2]; 
          p[3] = vert01[0]; p[4] = vert01[1]; p[5] = vert01[2]; 
          p[6] = vert03[0]; p[7] = vert03[1]; p[8] = vert03[2];
          p[9] = vert10[0]; p[10] = vert10[1]; p[11] = vert10[2]; 
          p[12] = vert11[0]; p[13] = vert11[1]; p[14] = vert11[2]; 
          p[15] = vert13[0]; p[16] = vert13[1]; p[17] = vert13[2];
          p += 18;

          c[0] = numThreadIndex + 0; c[1] = numThreadIndex + 1;
          c[2] = numThreadIndex + 2; c[3] = numThreadIndex + 3;
          c[4] = numThreadIndex + 4; c[5] = numThreadIndex + 5;
          c += 6;
          numThreadIndex += 6;

          m[0] = msm[vIds[0]]; m[1] = msm[vIds[1]];
          m += 2;
        }
      } else if(tetLookupIs3Label[lookupIndex]) {
        const int *vIds = tetLookupSplitBasisns3Label[lookupIndex];

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

        // Label 0
        p[0] = edge00[0]; p[1] = edge00[1]; p[2] = edge00[2];
        p[3] = edge02[0]; p[4] = edge02[1]; p[5] = edge02[2];
        p[6] = tri00[0]; p[7] = tri00[1]; p[8] = tri00[2];
        p[9] = edge02[0]; p[10] = edge02[1]; p[11] = edge02[2];
        p[12] = tri00[0]; p[13] = tri00[1]; p[14] = tri00[2];
        p[15] = tri01[0]; p[16] = tri01[1]; p[17] = tri01[2];        
        p[18] = edge01[0]; p[19] = edge01[1]; p[20] = edge01[2];
        p[21] = edge03[0]; p[22] = edge03[1]; p[23] = edge03[2];
        p[24] = tri00[0]; p[25] = tri00[1]; p[26] = tri00[2];
        p[27] = edge03[0]; p[28] = edge03[1]; p[29] = edge03[2];
        p[30] = tri00[0]; p[31] = tri00[1]; p[32] = tri00[2];
        p[33] = tri01[0]; p[34] = tri01[1]; p[35] = tri01[2];

        // Label 1           
        p[36] = edge10[0]; p[37] = edge10[1]; p[38] = edge10[2];
        p[39] = edge11[0]; p[40] = edge11[1]; p[41] = edge11[2];
        p[42] = tri10[0]; p[43] = tri10[1]; p[44] = tri10[2];
        p[45] = edge11[0]; p[46] = edge11[1]; p[47] = edge11[2];
        p[48] = tri10[0]; p[49] = tri10[1]; p[50] = tri10[2];
        p[51] = tri11[0]; p[52] = tri11[1]; p[53] = tri11[2];
        p[54] = edge12[0]; p[55] = edge12[1]; p[56] = edge12[2];
        p[57] = tri10[0]; p[58] = tri10[1]; p[59] = tri10[2];
        p[60] = tri11[0]; p[61] = tri11[1]; p[62] = tri11[2];

        // Label 2           
        p[63] = edge20[0]; p[64] = edge20[1]; p[65] = edge20[2];
        p[66] = edge21[0]; p[67] = edge21[1]; p[68] = edge21[2]; 
        p[69] = tri20[0]; p[70] = tri20[1]; p[71] = tri20[2];
        p[72] = edge21[0]; p[73] = edge21[1]; p[74] = edge21[2];
        p[75] = tri20[0]; p[76] = tri20[1]; p[77] = tri20[2];
        p[78] = tri21[0]; p[79] = tri21[1]; p[80] = tri21[2];
        p[81] = edge22[0]; p[82] = edge22[1]; p[83] = edge22[2]; 
        p[84] = tri20[0]; p[85] = tri20[1]; p[86] = tri20[2];
        p[87] = tri21[0]; p[88] = tri21[1]; p[89] = tri21[2];
        p += 90;

        c[0] = numThreadIndex + 0; c[1] = numThreadIndex + 1;
        c[2] = numThreadIndex + 2; c[3] = numThreadIndex + 3;
        c[4] = numThreadIndex + 4; c[5] = numThreadIndex + 5;
        c[6] = numThreadIndex + 6; c[7] = numThreadIndex + 7;
        c[8] = numThreadIndex + 8; c[9] = numThreadIndex + 9;
        c[10] = numThreadIndex + 10; c[11] = numThreadIndex + 11;
        c[12] = numThreadIndex + 12; c[13] = numThreadIndex + 13;
        c[14] = numThreadIndex + 14; c[15] = numThreadIndex + 15;
        c[16] = numThreadIndex + 16; c[17] = numThreadIndex + 17;
        c[18] = numThreadIndex + 18; c[19] = numThreadIndex + 19;
        c[20] = numThreadIndex + 20; c[21] = numThreadIndex + 21;
        c[22] = numThreadIndex + 22; c[23] = numThreadIndex + 23;
        c[24] = numThreadIndex + 24; c[25] = numThreadIndex + 25;
        c[26] = numThreadIndex + 26; c[27] = numThreadIndex + 27;
        c[28] = numThreadIndex + 28; c[29] = numThreadIndex + 29;
        c += 30;
        numThreadIndex += 30;

        m[0] = msm[vIds[0]]; m[1] = msm[vIds[0]];        
        m[2] = msm[vIds[0]]; m[3] = msm[vIds[0]];
        m[4] = msm[vIds[8]]; m[5] = msm[vIds[8]];        
        m[6] = msm[vIds[8]]; m[7] = msm[vIds[9]];       
        m[8] = msm[vIds[9]]; m[9] = msm[vIds[9]];
        m += 10;

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
        p[0] = vert00[0];   p[1] = vert00[1];   p[2] = vert00[2]; 
        p[3] = vert0t0[0];  p[4] = vert0t0[1];  p[5] = vert0t0[2]; 
        p[6] = vert0tet[0]; p[7] = vert0tet[1]; p[8] = vert0tet[2];
        p[9] = vert00[0];   p[10] = vert00[1];   p[11] = vert00[2]; 
        p[12] = vert0t1[0];  p[13] = vert0t1[1];  p[14] = vert0t1[2]; 
        p[15] = vert0tet[0]; p[16] = vert0tet[1]; p[17] = vert0tet[2];
        p[18] = vert01[0];   p[19] = vert01[1];   p[20] = vert01[2]; 
        p[21] = vert0t0[0];  p[22] = vert0t0[1];  p[23] = vert0t0[2]; 
        p[24] = vert0tet[0]; p[25] = vert0tet[1]; p[26] = vert0tet[2];
        p[27] = vert01[0];   p[28] = vert01[1];   p[29] = vert01[2]; 
        p[30] = vert0t2[0];  p[31] = vert0t2[1];  p[32] = vert0t2[2]; 
        p[33] = vert0tet[0]; p[34] = vert0tet[1]; p[35] = vert0tet[2];
        p[36] = vert02[0];   p[37] = vert02[1];   p[38] = vert02[2]; 
        p[39] = vert0t2[0];  p[40] = vert0t2[1];  p[41] = vert0t2[2]; 
        p[42] = vert0tet[0]; p[43] = vert0tet[1]; p[44] = vert0tet[2];
        p[45] = vert02[0];   p[46] = vert02[1];   p[47] = vert02[2]; 
        p[48] = vert0t1[0];  p[49] = vert0t1[1];  p[50] = vert0t1[2]; 
        p[51] = vert0tet[0]; p[52] = vert0tet[1]; p[53] = vert0tet[2];

        // Label Vert 1       
        p[54] = vert10[0];   p[55] = vert10[1];   p[56] = vert10[2]; 
        p[57] = vert1t0[0];  p[58] = vert1t0[1];  p[59] = vert1t0[2]; 
        p[60] = vert1tet[0]; p[61] = vert1tet[1]; p[62] = vert1tet[2];
        p[63] = vert10[0];   p[64] = vert10[1];   p[65] = vert10[2]; 
        p[66] = vert1t1[0];  p[67] = vert1t1[1];  p[68] = vert1t1[2]; 
        p[69] = vert1tet[0]; p[70] = vert1tet[1]; p[71] = vert1tet[2];
        p[72] = vert11[0];   p[73] = vert11[1];   p[74] = vert11[2]; 
        p[75] = vert1t0[0];  p[76] = vert1t0[1];  p[77] = vert1t0[2]; 
        p[78] = vert1tet[0]; p[79] = vert1tet[1]; p[80] = vert1tet[2];
        p[81] = vert11[0];   p[82] = vert11[1];   p[83] = vert11[2]; 
        p[84] = vert1t2[0];  p[85] = vert1t2[1];  p[86] = vert1t2[2]; 
        p[87] = vert1tet[0]; p[88] = vert1tet[1]; p[89] = vert1tet[2];
        p[90] = vert12[0];   p[91] = vert12[1];   p[92] = vert12[2]; 
        p[93] = vert1t2[0];  p[94] = vert1t2[1];  p[95] = vert1t2[2]; 
        p[96] = vert1tet[0]; p[97] = vert1tet[1]; p[98] = vert1tet[2];
        p[99] = vert12[0];   p[100] = vert12[1];   p[101] = vert12[2]; 
        p[102] = vert1t1[0];  p[103] = vert1t1[1];  p[104] = vert1t1[2]; 
        p[105] = vert1tet[0]; p[106] = vert1tet[1]; p[107] = vert1tet[2];


        // Label Vert 2     
        p[108] = vert20[0];   p[109] = vert20[1];   p[110] = vert20[2]; 
        p[111] = vert2t0[0];  p[112] = vert2t0[1];  p[113] = vert2t0[2]; 
        p[114] = vert2tet[0]; p[115] = vert2tet[1]; p[116] = vert2tet[2];
        p[117] = vert20[0];   p[118] = vert20[1];   p[119] = vert20[2]; 
        p[120] = vert2t1[0];  p[121] = vert2t1[1];  p[122] = vert2t1[2]; 
        p[123] = vert2tet[0]; p[124] = vert2tet[1]; p[125] = vert2tet[2];
        p[126] = vert21[0];   p[127] = vert21[1];   p[128] = vert21[2]; 
        p[129] = vert2t0[0];  p[130] = vert2t0[1];  p[131] = vert2t0[2]; 
        p[132] = vert2tet[0]; p[133] = vert2tet[1]; p[134] = vert2tet[2];
        p[135] = vert21[0];   p[136] = vert21[1];   p[137] = vert21[2]; 
        p[138] = vert2t2[0];  p[139] = vert2t2[1];  p[140] = vert2t2[2]; 
        p[141] = vert2tet[0]; p[142] = vert2tet[1]; p[143] = vert2tet[2];
        p[144] = vert22[0];   p[145] = vert22[1];   p[146] = vert22[2]; 
        p[147] = vert2t2[0];  p[148] = vert2t2[1];  p[149] = vert2t2[2]; 
        p[150] = vert2tet[0]; p[151] = vert2tet[1]; p[152] = vert2tet[2];
        p[153] = vert22[0];   p[154] = vert22[1];   p[155] = vert22[2]; 
        p[156] = vert2t1[0];  p[157] = vert2t1[1];  p[158] = vert2t1[2]; 
        p[159] = vert2tet[0]; p[160] = vert2tet[1]; p[161] = vert2tet[2];

        // Label Vert 3      
        p[162] = vert30[0];   p[163] = vert30[1];   p[164] = vert30[2]; 
        p[165] = vert3t0[0];  p[166] = vert3t0[1];  p[167] = vert3t0[2]; 
        p[168] = vert3tet[0]; p[169] = vert3tet[1]; p[170] = vert3tet[2];            
        p[171] = vert30[0];   p[172] = vert30[1];   p[173] = vert30[2]; 
        p[174] = vert3t1[0];  p[175] = vert3t1[1];  p[176] = vert3t1[2]; 
        p[177] = vert3tet[0]; p[178] = vert3tet[1]; p[179] = vert3tet[2];
        p[180] = vert31[0];   p[181] = vert31[1];   p[182] = vert31[2]; 
        p[183] = vert3t0[0];  p[184] = vert3t0[1];  p[185] = vert3t0[2]; 
        p[186] = vert3tet[0]; p[187] = vert3tet[1]; p[188] = vert3tet[2];
        p[189] = vert31[0];   p[190] = vert31[1];   p[191] = vert31[2]; 
        p[192] = vert3t2[0];  p[193] = vert3t2[1];  p[194] = vert3t2[2]; 
        p[195] = vert3tet[0]; p[196] = vert3tet[1]; p[197] = vert3tet[2];
        p[198] = vert32[0];   p[199] = vert32[1];   p[200] = vert32[2]; 
        p[201] = vert3t2[0];  p[202] = vert3t2[1];  p[203] = vert3t2[2]; 
        p[204] = vert3tet[0]; p[205] = vert3tet[1]; p[206] = vert3tet[2];
        p[207] = vert32[0];   p[208] = vert32[1];   p[209] = vert32[2]; 
        p[210] = vert3t1[0];  p[211] = vert3t1[1];  p[212] = vert3t1[2]; 
        p[213] = vert3tet[0]; p[214] = vert3tet[1]; p[215] = vert3tet[2];
        p += 216;

        c[0] = numThreadIndex + 0; c[1] = numThreadIndex + 1;
        c[2] = numThreadIndex + 2; c[3] = numThreadIndex + 3;
        c[4] = numThreadIndex + 4; c[5] = numThreadIndex + 5;
        c[6] = numThreadIndex + 6; c[7] = numThreadIndex + 7;
        c[8] = numThreadIndex + 8; c[9] = numThreadIndex + 9;
        c[10] = numThreadIndex + 10; c[11] = numThreadIndex + 11;
        c[12] = numThreadIndex + 12; c[13] = numThreadIndex + 13;
        c[14] = numThreadIndex + 14; c[15] = numThreadIndex + 15;
        c[16] = numThreadIndex + 16; c[17] = numThreadIndex + 17;
        c[18] = numThreadIndex + 18; c[19] = numThreadIndex + 19;
        c[20] = numThreadIndex + 20; c[21] = numThreadIndex + 21;
        c[22] = numThreadIndex + 22; c[23] = numThreadIndex + 23;
        c[24] = numThreadIndex + 24; c[25] = numThreadIndex + 25;
        c[26] = numThreadIndex + 26; c[27] = numThreadIndex + 27;
        c[28] = numThreadIndex + 28; c[29] = numThreadIndex + 29;
        c[30] = numThreadIndex + 30; c[31] = numThreadIndex + 31;
        c[32] = numThreadIndex + 32; c[33] = numThreadIndex + 33;
        c[34] = numThreadIndex + 34; c[35] = numThreadIndex + 35;
        c[36] = numThreadIndex + 36; c[37] = numThreadIndex + 37;
        c[38] = numThreadIndex + 38; c[39] = numThreadIndex + 39;
        c[40] = numThreadIndex + 40; c[41] = numThreadIndex + 41;
        c[42] = numThreadIndex + 42; c[43] = numThreadIndex + 43;
        c[44] = numThreadIndex + 44; c[45] = numThreadIndex + 45;
        c[46] = numThreadIndex + 46; c[47] = numThreadIndex + 47;
        c[48] = numThreadIndex + 48; c[49] = numThreadIndex + 49;
        c[50] = numThreadIndex + 50; c[51] = numThreadIndex + 51;
        c[52] = numThreadIndex + 52; c[53] = numThreadIndex + 53;
        c[54] = numThreadIndex + 54; c[55] = numThreadIndex + 55;
        c[56] = numThreadIndex + 56; c[57] = numThreadIndex + 57;
        c[58] = numThreadIndex + 58; c[59] = numThreadIndex + 59;
        c[60] = numThreadIndex + 60; c[61] = numThreadIndex + 61;
        c[62] = numThreadIndex + 62; c[63] = numThreadIndex + 63;
        c[64] = numThreadIndex + 64; c[65] = numThreadIndex + 65;
        c[66] = numThreadIndex + 66; c[67] = numThreadIndex + 67;
        c[68] = numThreadIndex + 68; c[69] = numThreadIndex + 69;
        c[70] = numThreadIndex + 70; c[71] = numThreadIndex + 71;
        c += 72;
        numThreadIndex += 72;

        m[0] = msm[0]; m[1] = msm[0];
        m[2] = msm[0]; m[3] = msm[0];
        m[4] = msm[0]; m[5] = msm[0];
        m[6] = msm[1]; m[7] = msm[1];
        m[8] = msm[1]; m[9] = msm[1];
        m[10] = msm[1]; m[11] = msm[1];
        m[12] = msm[2]; m[13] = msm[2];
        m[14] = msm[2]; m[15] = msm[2];
        m[16] = msm[2]; m[17] = msm[2];
        m[18] = msm[3]; m[19] = msm[3];
        m[20] = msm[3]; m[21] = msm[3];
        m[22] = msm[3]; m[23] = msm[3];
        m += 24;
      }
    }
  }
}
#endif // TTK_ENABLE_OPENMP

  this->printMsg("Computed 2-Separatrices 3D[Fine]",
                 1, localTimer.getElapsedTime(), this->threadNumber_);

  return 0;
}

template <typename triangulationType>
int ttk::MorseSmaleSegmentationPL::computeBasinSeparation_3D_fast(
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

#ifdef TTK_ENABLE_OPENMP
if(threadNumber_ == 1) {
#endif // TTK_ENABLE_OPENMP

  std::vector<std::pair<SimplexId, unsigned char>> validCases;
  SimplexId numTriangles = 0;

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

    if(tetLookupFast[lookupIndex]) {
      validCases.push_back(std::make_pair(tet, lookupIndex));
      numTriangles += tetLookupNumWallTriangles[lookupIndex];
    }
  }

  outputSeparatrices2_points_->resize(9 * numTriangles);
  outputSeparatrices2_cells_connectivity_->resize(3 * numTriangles);
  outputSeparatrices2_cells_mscIds_->resize(numTriangles);
  *outputSeparatrices2_numberOfPoints_ = 3 * numTriangles;
  *outputSeparatrices2_numberOfCells_ = numTriangles;

  auto &points = *outputSeparatrices2_points_;
  auto &cellsConn = *outputSeparatrices2_cells_connectivity_;
  auto &cellsMSCIds = *outputSeparatrices2_cells_mscIds_;

  float* p = points.data();
  SimplexId* c = cellsConn.data();
  SimplexId* m = cellsMSCIds.data();
  SimplexId cellIndex = 0;
  
  for (const auto& vCase : validCases) {
    const SimplexId& tet = vCase.first;
    const unsigned char& lookupIndex = vCase.second;

    SimplexId vertices[4];
    triangulation.getCellVertex(tet, 0, vertices[0]);
    triangulation.getCellVertex(tet, 1, vertices[1]);
    triangulation.getCellVertex(tet, 2, vertices[2]);
    triangulation.getCellVertex(tet, 3, vertices[3]);

    const SimplexId msm[4] = {
      morseSmaleManifold[vertices[0]], morseSmaleManifold[vertices[1]],
      morseSmaleManifold[vertices[2]], morseSmaleManifold[vertices[3]]};

    const int id0 =
      tetLookupFastTri[tetLookupFastCase[lookupIndex]][0];
    const int id1 =
      tetLookupFastTri[tetLookupFastCase[lookupIndex]][1];
    const int id2 =
      tetLookupFastTri[tetLookupFastCase[lookupIndex]][2];
    
    float vertPos[4][3];
    triangulation.getVertexPoint(
      vertices[0], vertPos[0][0], vertPos[0][1], vertPos[0][2]);
    triangulation.getVertexPoint(
      vertices[1], vertPos[1][0], vertPos[1][1], vertPos[1][2]);
    triangulation.getVertexPoint(
      vertices[2], vertPos[2][0], vertPos[2][1], vertPos[2][2]);
    triangulation.getVertexPoint(
      vertices[3], vertPos[3][0], vertPos[3][1], vertPos[3][2]);

    p[0] = vertPos[id0][0]; p[1] = vertPos[id0][1]; p[2] = vertPos[id0][2];
    p[3] = vertPos[id1][0]; p[4] = vertPos[id1][1]; p[5] = vertPos[id1][2];
    p[6] = vertPos[id2][0]; p[7] = vertPos[id2][1]; p[8] = vertPos[id2][2];
    p += 9;

    c[0] = cellIndex + 0;
    c[1] = cellIndex + 1;
    c[2] = cellIndex + 2;
    c += 3;
    cellIndex += 3;

    m[0] = msm[id0];
    m += 1;
  }

#ifdef TTK_ENABLE_OPENMP
} else {
  size_t triangleStartIndex[threadNumber_ + 1];
  #pragma omp parallel num_threads(threadNumber_)
  {
    const int tid = omp_get_thread_num();
    std::vector<std::pair<SimplexId, unsigned char>> validCases;
    size_t numThreadTriangles = 0;
    size_t numThreadIndex = 0;

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

      if(tetLookupFast[lookupIndex]) {
        validCases.push_back(std::make_pair(tet, lookupIndex));
        numThreadTriangles += 1;
      }  
    }

    triangleStartIndex[tid + 1] = numThreadTriangles;

    #pragma omp barrier

    #pragma omp single
    {
      triangleStartIndex[0] = 0;

      // Count triangle number and create iterator start indices
      for(int t = 1; t < threadNumber_ + 1; ++t) {
        triangleStartIndex[t] += triangleStartIndex[t-1];
      }
    }

    const size_t numTriangles = triangleStartIndex[threadNumber_];
    numThreadIndex = triangleStartIndex[tid];

    #pragma omp single nowait
    outputSeparatrices2_points_->resize(9 * numTriangles);

    #pragma omp single nowait
    outputSeparatrices2_cells_connectivity_->resize(3 * numTriangles);

    #pragma omp single nowait
    outputSeparatrices2_cells_mscIds_->resize(numTriangles);

    #pragma omp single
    {
      *outputSeparatrices2_numberOfPoints_ = 3 * numTriangles;
      *outputSeparatrices2_numberOfCells_ = numTriangles;
    }
  
    auto &points = *outputSeparatrices2_points_;
    auto &cellsConn = *outputSeparatrices2_cells_connectivity_;
    auto &cellsMSCIds = *outputSeparatrices2_cells_mscIds_;

    float* p = points.data();
    SimplexId* c = cellsConn.data();
    SimplexId* m = cellsMSCIds.data();

    p += (numThreadIndex * 9);
    c += (numThreadIndex * 3);
    m += numThreadIndex;

    numThreadIndex = 3 * numThreadIndex;

    for (const auto& vCase : validCases)
    {
      const SimplexId& tet = vCase.first;
      const unsigned char& lookupIndex = vCase.second;

      SimplexId vertices[4];
      triangulation.getCellVertex(tet, 0, vertices[0]);
      triangulation.getCellVertex(tet, 1, vertices[1]);
      triangulation.getCellVertex(tet, 2, vertices[2]);
      triangulation.getCellVertex(tet, 3, vertices[3]);

      const SimplexId msm[4] = {
        morseSmaleManifold[vertices[0]], morseSmaleManifold[vertices[1]],
        morseSmaleManifold[vertices[2]], morseSmaleManifold[vertices[3]]};

      const int id0 =
        tetLookupFastTri[tetLookupFastCase[lookupIndex]][0];
      const int id1 =
        tetLookupFastTri[tetLookupFastCase[lookupIndex]][1];
      const int id2 =
        tetLookupFastTri[tetLookupFastCase[lookupIndex]][2];
      
      float vertPos[4][3];
      triangulation.getVertexPoint(
        vertices[0], vertPos[0][0], vertPos[0][1], vertPos[0][2]);
      triangulation.getVertexPoint(
        vertices[1], vertPos[1][0], vertPos[1][1], vertPos[1][2]);
      triangulation.getVertexPoint(
        vertices[2], vertPos[2][0], vertPos[2][1], vertPos[2][2]);
      triangulation.getVertexPoint(
        vertices[3], vertPos[3][0], vertPos[3][1], vertPos[3][2]);

      p[0] = vertPos[id0][0]; p[1] = vertPos[id0][1]; p[2] = vertPos[id0][2];
      p[3] = vertPos[id1][0]; p[4] = vertPos[id1][1]; p[5] = vertPos[id1][2];
      p[6] = vertPos[id2][0]; p[7] = vertPos[id2][1]; p[8] = vertPos[id2][2];
      p += 9;

      c[0] = numThreadIndex + 0;
      c[1] = numThreadIndex + 1;
      c[2] = numThreadIndex + 2;
      c += 3;
      numThreadIndex += 3;

      m[0] = msm[id0];
      m += 1;
    }
  }
}
#endif // TTK_ENABLE_OPENMP

  this->printMsg("Computed 2-Separatrices 3D[Fast]",
                 1, localTimer.getElapsedTime(), this->threadNumber_);

  return 0;
}

template <typename triangulationType>
int ttk::MorseSmaleSegmentationPL::computeSaddleExtremaConnectors_3D(
  std::vector<mscq::IntegralLine> &sadExtrConns,
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

#ifdef TTK_ENABLE_OPENMP
if(threadNumber_ == 1) {
#endif // TTK_ENABLE_OPENMP

  std::vector<std::pair<SimplexId, SimplexId>> reachableExtrema[2];
  
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
    mscq::IntegralLine sadExtrConn(c.first);
    sadExtrConn.type_ = 0;

    SimplexId currentVert = c.second;

    while(ascNeighbor[currentVert] != currentVert) {
      sadExtrConn.geometry_.push_back(currentVert);
      currentVert = ascNeighbor[currentVert];
    }

    sadExtrConn.geometry_.push_back(currentVert);

    sadExtrConns.push_back(sadExtrConn);
  }

  std::map<long long, size_t> s2ToMax;
  for(const auto& c : reachableExtrema[1]) { // 2-saddle maximum
    mscq::IntegralLine sadExtrConn(c.first);
    sadExtrConn.type_ = 1;

    SimplexId currentVert = c.second;

    while(dscNeighbor[currentVert] != currentVert) {
      sadExtrConn.geometry_.push_back(currentVert);
      currentVert = dscNeighbor[currentVert];
    }

    sadExtrConn.geometry_.push_back(currentVert);

    sadExtrConns.push_back(sadExtrConn);
  }

#ifdef TTK_ENABLE_OPENMP
} else {
  size_t sadExtrStartIndex[threadNumber_ + 1];
  #pragma omp parallel num_threads(threadNumber_)
  {
    const int tid = omp_get_thread_num();
    std::vector<mscq::IntegralLine> sadExtrConnsLocal;

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
      mscq::IntegralLine sadExtrConn(c.first);
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
      mscq::IntegralLine sadExtrConn(c.first);
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
}
#endif // TTK_ENABLE_OPENMP

  this->printMsg("Computed crit point connectors",
                 1, localTimer.getElapsedTime(), this->threadNumber_);

  return 1;
}

template <typename triangulationType>
int ttk::MorseSmaleSegmentationPL::computeSaddleExtremaConnectors_3D_fast(
  std::vector<mscq::IntegralLine> &sadExtrConns,
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

#ifdef TTK_ENABLE_OPENMP
if(threadNumber_ == 1) {
#endif // TTK_ENABLE_OPENMP

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
          mscq::IntegralLine sadExtrConn(crit.first);
          if(isSaddle1) {
            sadExtrConn.type_ = 0;
          } else {
            sadExtrConn.type_ = 1;
          }

          sadExtrConn.geometry_.push_back(neighbor);
          sadExtrConns.push_back(sadExtrConn);

          neigborSet.insert(neighbor);
        }
      }
    }
  }

#ifdef TTK_ENABLE_OPENMP
} else {
  size_t sadExtrStartIndex[threadNumber_ + 1];
  #pragma omp parallel num_threads(threadNumber_)
  {
    const int tid = omp_get_thread_num();
    std::vector<mscq::IntegralLine> sadExtrConnsLocal;

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
            mscq::IntegralLine sadExtrConn(crit.first);
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
}
#endif // TTK_ENABLE_OPENMP

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

size_t numcritPoints = criticalPoints.size();

#ifdef TTK_ENABLE_OPENMP
if(threadNumber_ == 1) {
#endif // TTK_ENABLE_OPENMP

  outputCriticalPoints_points_->resize(numcritPoints);
  outputCriticalPoints_points_cellDimensions_->resize(numcritPoints);
  outputCriticalPoints_points_cellIds_->resize(numcritPoints);

  auto &points = *outputCriticalPoints_points_;
  auto &cellsDim = *outputCriticalPoints_points_cellDimensions_;
  auto &cellsIds = *outputCriticalPoints_points_cellIds_;

  float x, y, z;

  for(size_t i = 0; i < numcritPoints; ++i) {
    const auto& crit = criticalPoints[i];

    triangulation.getVertexPoint(crit.first, x, y, z);

    points[i] = {x, y, z};
    cellsDim[i] = crit.second;
    cellsIds[i] = crit.first;
  }

#ifdef TTK_ENABLE_OPENMP
} else {
  #pragma omp parallel num_threads(threadNumber_)
  {
    #pragma omp single nowait
    outputCriticalPoints_points_->resize(numcritPoints);
    #pragma omp single nowait
    outputCriticalPoints_points_cellDimensions_->resize(numcritPoints);
    #pragma omp single
    outputCriticalPoints_points_cellIds_->resize(numcritPoints);

    auto &points = *outputCriticalPoints_points_;
    auto &cellsDim = *outputCriticalPoints_points_cellDimensions_;
    auto &cellsIds = *outputCriticalPoints_points_cellIds_;

    float x, y, z;

    #pragma omp for schedule(static)
    for(size_t i = 0; i < numcritPoints; ++i) {
      const auto& crit = criticalPoints[i];

      triangulation.getVertexPoint(crit.first, x, y, z);

      points[i] = {x, y, z};
      cellsDim[i] = crit.second;
      cellsIds[i] = crit.first;
    }
  }
}
#endif // TTK_ENABLE_OPENMP

  this->printMsg("Wrote " + std::to_string(criticalPoints.size()) +
                 " critical points",
                 1, localTimer.getElapsedTime(), this->threadNumber_);

  return 1;
}

template <typename triangulationType>
int ttk::MorseSmaleSegmentationPL::setSeparatrices1_3D(
  const std::vector<mscq::IntegralLine> &separatrices,
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
#endif // TTK_ENABLE_KAMIKAZE

  this->printMsg("Writing 1-separatrices", 1, localTimer.getElapsedTime(),
                 this->threadNumber_, ttk::debug::LineMode::REPLACE);

  // number of separatrices
  const size_t numSep = separatrices.size();
  size_t npoints = 0;
  size_t numCells = 0;

#ifdef TTK_ENABLE_OPENMP
if(threadNumber_ == 1) {
#endif // TTK_ENABLE_OPENMP

  // count total number of points and cells
  for(size_t i = 0; i < numSep; ++i) {
    const auto &sep = separatrices[i];

    const int numSepPts = sep.geometry_.size();

    npoints += numSepPts;
    numCells += numSepPts - 1;
  }

  // resize arrays
  outputSeparatrices1_points_->resize(npoints * 3);
  auto &points = *outputSeparatrices1_points_;
  outputSeparatrices1_cells_connectivity_->resize(numCells * 2);
  auto &cellsConn = *outputSeparatrices1_cells_connectivity_;
  outputSeparatrices1_cells_separatrixTypes_->resize(numCells);
  auto &sepType = *outputSeparatrices1_cells_separatrixTypes_;
  outputseparatrices1_cells_extremaDistance_->resize(numCells);
  auto &extremaDistance = *outputseparatrices1_cells_extremaDistance_;
  outputseparatrices1_cells_extremaDistanceAbs_->resize(numCells);
  auto &extremaDistanceAbs = *outputseparatrices1_cells_extremaDistanceAbs_;

  SimplexId currPId = 0, currCId = 0;
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

      sepType[currCId] = sep.type_;

      const int jInv = geoSize - (j - 1);
      extremaDistance[currCId] = (float)(jInv) / (float)geoSize;
      extremaDistanceAbs[currCId] = jInv;

      currPId += 1;
      currCId += 1;
    }
  }

#ifdef TTK_ENABLE_OPENMP
} else {
  SimplexId StartPId[numSep];
  SimplexId StartCId[numSep];

  #pragma omp parallel num_threads(threadNumber_) shared(npoints, numCells, separatrices)
  {
     // count total number of points and cells
    #pragma omp for schedule(static) reduction(+: npoints) reduction(+: numCells) nowait
    for(size_t i = 0; i < numSep; ++i) {
      const auto &sep = separatrices[i];
      const int numSepPts = sep.geometry_.size();

      npoints += numSepPts;
      numCells += numSepPts - 1;
    }

    #pragma omp single
    {
      StartPId[0] = 0;
      StartCId[0] = 0;

      for(size_t i = 1; i < numSep; ++i) {
        const auto &sep = separatrices[i - 1];
        const int numSepPts = sep.geometry_.size();

        StartPId[i] = StartPId[i-1] + numSepPts;
        StartCId[i] = StartCId[i-1] + numSepPts - 1;
      }
    }

    // resize arrays
  #pragma omp single nowait
    outputSeparatrices1_points_->resize(npoints * 3);
  #pragma omp single nowait
    outputSeparatrices1_cells_connectivity_->resize(numCells * 2);
  #pragma omp single nowait
    outputSeparatrices1_cells_separatrixTypes_->resize(numCells);
  #pragma omp single nowait
    outputseparatrices1_cells_extremaDistance_->resize(numCells);
  #pragma omp single
    outputseparatrices1_cells_extremaDistanceAbs_->resize(numCells);

    auto &points = *outputSeparatrices1_points_;
    auto &cellsConn = *outputSeparatrices1_cells_connectivity_;
    auto &sepType = *outputSeparatrices1_cells_separatrixTypes_;
    auto &extremaDistance = *outputseparatrices1_cells_extremaDistance_;
    auto &extremaDistanceAbs = *outputseparatrices1_cells_extremaDistanceAbs_;

    #pragma omp for schedule(static)
    for(size_t i = 0; i < numSep; ++i) {
      const auto &sep = separatrices[i];

      SimplexId currPId = StartPId[i], currCId = StartCId[i];

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

        sepType[currCId] = sep.type_;

        const int jInv = geoSize - (j - 1);
        extremaDistance[currCId] = (float)(jInv) / (float)geoSize;
        extremaDistanceAbs[currCId] = jInv;

        currPId += 1;
        currCId += 1;
      }
    }
  }
}
#endif // TTK_ENABLE_OPENMP

  // update pointers
  *outputSeparatrices1_numberOfPoints_ = npoints;
  *outputSeparatrices1_numberOfCells_ = numCells;

  this->printMsg("Wrote " + std::to_string(numSep) + " 1-separatrices",
    1, localTimer.getElapsedTime(), this->threadNumber_);

  return 0;
}
