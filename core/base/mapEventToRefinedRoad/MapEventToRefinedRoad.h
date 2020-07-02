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
  int &nEventPoints,
  int &nRroadPoints,
  const idType &threshold4closeness) const {
  // initialize the eventSample property for the rRoadpoints
  for(int i = 0; i < nRroadPoints; i++) {
    eventSample2rRoadPoint[i] = 0.0;
  }

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif

  for(int eventPointIndex = 0; eventPointIndex < nEventPoints;
      eventPointIndex++) {
    int closetPoint = -1;
    int closetDist = std::numeric_limits<int>::max();

    auto eventLon = originalEventPoints[eventPointIndex * 3];
    auto eventLat = originalEventPoints[eventPointIndex * 3 + 1];

    // find out the closet rRoadPoints within threashould
    for(int rRoadPointInex = 0; rRoadPointInex < nRroadPoints;
        rRoadPointInex++) {
      auto rRoadPLon = refinedRoadPoints[rRoadPointInex * 3];
      auto rRoadPLat = refinedRoadPoints[rRoadPointInex * 3 + 1];
      auto distance = ttk::RoadRefinement::distance4TwoPoints(
        eventLat, eventLon, rRoadPLat, rRoadPLon);

      if(distance <= threshold4closeness and distance < closetDist) {
        closetPoint = rRoadPointInex;
        closetDist = distance;
      }
    }
    // calculate the event influence to each closet point and assign the value
    // to the eventSample property of each point.

    if(closetPoint == -1) {
      continue;
    } else {
      mtx.lock();
      eventSample2rRoadPoint[closetPoint] += 1;
      mtx.unlock();
    }
  }

  return 0;
}
