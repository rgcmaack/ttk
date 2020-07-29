/// \ingroup base
/// \class ttk::MapEventToRefinedRoad
/// \author Jonas Lukasczyk <jl@jluk.de>
/// \date 1.09.2019
///
/// TODO

#pragma once

// ttk common includes
#include <Debug.h>
#include <RoadRefinement.h>
#include <limits>
#include <math.h>
#include <mutex>

std::mutex mtx;

namespace ttk {

  class MapEventToRefinedRoad : virtual public Debug {

  public:
    MapEventToRefinedRoad() {
      this->setDebugMsgPrefix(
        "MapEventToRefinedRoad"); // inherited from Debug: prefix will be
                                  // printed at the beginning of every msg
    };
    ~MapEventToRefinedRoad(){};

    template <class idType>
    int calculateClosetPointofRoadforEventData(
      const idType *originalEventPoints,
      const idType *refinedRoadPoints,
      idType *eventSample2rRoadPoint,
      int *gridCellIndex,
      int &nEventPoints,
      int &nRroadPoints,
      const idType &threshold4closeness) const;
  };
} // namespace ttk

template <class idType>
int ttk::MapEventToRefinedRoad::calculateClosetPointofRoadforEventData(
  const idType *originalEventPoints,
  const idType *refinedRoadPoints,
  idType *eventSample2rRoadPoint,
  int *gridCellIndex,
  int &nEventPoints,
  int &nRroadPoints,
  const idType &threshold4closeness) const {
  // initialize the eventSample property for the rRoadpoints and find the
  // boarding box of these points
  float maxLat = refinedRoadPoints[1], maxLong = refinedRoadPoints[0],
        minLat = refinedRoadPoints[1], minLong = refinedRoadPoints[0];
  for(int i = 0; i < nRroadPoints; i++) {
    eventSample2rRoadPoint[i] = 0.0;
    auto curLong = refinedRoadPoints[i * 3];
    auto curLat = refinedRoadPoints[i * 3 + 1];
    maxLat = std::max(maxLat, curLat);
    maxLong = std::max(maxLong, curLong);
    minLat = std::min(minLat, curLat);
    minLong = std::min(minLong, curLong);
  }
  for(int j = 0; j < nEventPoints; j++) {
    auto curLong = originalEventPoints[j * 3];
    auto curLat = originalEventPoints[j * 3 + 1];
    maxLat = std::max(maxLat, curLat);
    maxLong = std::max(maxLong, curLong);
    minLat = std::min(minLat, curLat);
    minLong = std::min(minLong, curLong);
  }
  // get the size of the bbx. m * n grid for the boarding box includding road
  // network and eventData
  float bbxWidth
    = ttk::RoadRefinement::distance4TwoPoints(minLat, maxLong, minLat, minLong);
  float bbxLong
    = ttk::RoadRefinement::distance4TwoPoints(maxLat, maxLong, minLat, maxLong);
  int n = ceil(bbxWidth / threshold4closeness),
      m = ceil(bbxLong / threshold4closeness);
  float colUnit = (maxLong - minLong) / n;
  float rowUnit = (maxLat - minLat) / m;

  std::cout << "M: " << m << std::endl;
  std::cout << "N: " << n << std::endl;

  // calculate the gridcell the rRoad points belong to and initialize the
  // gridVector
  std::vector<int> grids[m * n];
  for(int i = 0; i < nRroadPoints; i++) {
    auto curLong = refinedRoadPoints[i * 3];
    auto curLat = refinedRoadPoints[i * 3 + 1];
    int y = floor((curLong - minLong) / colUnit);
    int x = floor((curLat - minLat) / rowUnit);
    if(y == n)
      y = n - 1;
    if(x == m)
      x = m - 1;

    int index = n * x + y;
    grids[index].push_back(i);
    gridCellIndex[i] = index;
  }
  std::cout << "Finish generating grid for road network: " << m * n
            << std::endl;
  ///////////////////////////////////////print the nPoints in the gridCell
  // for(size_t iGird = 0; iGird < nRroadPoints; iGird++) {
  //   std::cout << "gridCellIndex: " << gridCellIndex[iGird] << std::endl;
  // }

// parallelly to
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif

  for(int eventPointIndex = 0; eventPointIndex < nEventPoints;
      eventPointIndex++) {

    int closetPoint = -1;
    double closetDist = threshold4closeness + 1.0;

    auto eventLon = originalEventPoints[eventPointIndex * 3];
    auto eventLat = originalEventPoints[eventPointIndex * 3 + 1];

    // find the current cell for the event data.
    int yIndex = floor((eventLon - minLong) / colUnit);
    int xIndex = floor((eventLat - minLat) / rowUnit);
    if(yIndex == n)
      yIndex = n - 1;
    if(xIndex == m)
      xIndex = m - 1;

    int delta[9][2] = {{0, 0}, {0, 1},  {0, -1}, {1, 0},  {-1, 0},
                       {1, 1}, {1, -1}, {-1, 1}, {-1, -1}};

    // calculate the 9 cells to find the close refined road Points
    for(size_t iCell = 0; iCell < 9; iCell++) {
      // get the cell x, y index
      int cellX = xIndex + delta[iCell][0], cellY = yIndex + delta[iCell][1];
      if(cellX < 0 || cellY < 0 || cellX >= m || cellY >= n)
        continue;

      // get the index of the cell in the grid vector
      auto indexCell4grids = n * cellX + cellY;
      auto targetCell = grids[indexCell4grids];
      // access the refined road points in the cell
      for(std::vector<int>::iterator it = targetCell.begin();
          it != targetCell.end(); ++it) {
        auto pIndex = *it;
        auto plong = refinedRoadPoints[pIndex * 3];
        auto plat = refinedRoadPoints[pIndex * 3 + 1];
        auto distance = ttk::RoadRefinement::distance4TwoPoints(
          eventLat, eventLon, plat, plong);

        if(distance < closetDist) {
          closetPoint = pIndex;
          closetDist = distance;
        }
      }
    } // end the search for the 9 cells
    // get the closet points for the event data

    if(closetPoint == -1)
      continue;

    if(closetDist > threshold4closeness)
      continue;

    mtx.lock();
    eventSample2rRoadPoint[closetPoint] += 1;
    mtx.unlock();
  }

  return 0;
}
