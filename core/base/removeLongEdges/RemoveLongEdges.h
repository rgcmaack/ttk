/// \ingroup base
/// \class ttk::RemoveLongEdges
/// \author Jonas Lukasczyk <jl@jluk.de>
/// \date 1.09.2019
///
/// TODO

#pragma once

// ttk common includes
#include <Debug.h>
#include <RoadRefinement.h>
#include <math.h>

namespace ttk {

  class RemoveLongEdges : virtual public Debug {

  public:
    RemoveLongEdges() {
      this->setDebugMsgPrefix(
        "RemoveLongEdges"); // inherited from Debug: prefix will be printed at
                            // the beginning of every msg
    };
    ~RemoveLongEdges(){};

    template <class idType>
    int getNumberOfnewEdges(const float *originalPointCoords,
                            const idType *originalConnectivityList,
                            const idType &nOriginalCells,
                            int &nNewEdges,
                            const float edgeThreshold) const;

    template <class idType>
    int removeEdges(const float *originalPointCoords,
                    const idType *originalConnectivityList,
                    const float edgeThreshold,
                    idType *cellConnectivityData,
                    const idType &nOriginalCells) const;
  };
} // namespace ttk

template <class idType>
int ttk::RemoveLongEdges::getNumberOfnewEdges(
  const float *originalPointCoords,
  const idType *originalConnectivityList,
  const idType &nOriginalCells,
  int &nNewEdges,
  const float edgeThreshold) const {
  // get the edge
  for(int index4edge = 0; index4edge < nOriginalCells; index4edge++) {
    auto p1 = originalConnectivityList[index4edge * 2];
    auto p2 = originalConnectivityList[index4edge * 2 + 1];
    float lat4p1 = originalPointCoords[p1 * 3 + 1];
    float long4p1 = originalPointCoords[p1 * 3];
    float lat4p2 = originalPointCoords[p2 * 3 + 1];
    float long4p2 = originalPointCoords[p2 * 3];
    //  std::cout << "lat4p1 " << lat4p1 << std::endl;
    //  std::cout << "long4p1 " << long4p1 << std::endl;
    auto distance = ttk::RoadRefinement::distance4TwoPoints(
      lat4p1, long4p1, lat4p2, long4p2);
    if((distance == 0) || (distance > edgeThreshold))
      continue;
    nNewEdges++;
  }

  return 0;
}

template <class idType>
int ttk::RemoveLongEdges::removeEdges(const float *originalPointCoords,
                                      const idType *originalConnectivityList,
                                      const float edgeThreshold,
                                      idType *cellConnectivityData,
                                      const idType &nOriginalCells) const {

  int nTotalEdges = 0;

  for(int cellIndex = 0; cellIndex < nOriginalCells; cellIndex++) {
    auto p1 = originalConnectivityList[cellIndex * 2];
    auto p2 = originalConnectivityList[cellIndex * 2 + 1];

    float p1Lat = originalPointCoords[p1 * 3 + 1];
    float p1Lon = originalPointCoords[p1 * 3];
    float p2Lat = originalPointCoords[p2 * 3 + 1];
    float p2Lon = originalPointCoords[p2 * 3];

    auto distance
      = ttk::RoadRefinement::distance4TwoPoints(p1Lat, p1Lon, p2Lat, p2Lon);
    if((distance == 0) || (distance > edgeThreshold))
      continue;

    cellConnectivityData[nTotalEdges * 2] = p1;
    cellConnectivityData[nTotalEdges * 2 + 1] = p2;
    nTotalEdges++;
  }

  std::cout << "nTotalEdges: " << nTotalEdges << std::endl;

  return 0;
}