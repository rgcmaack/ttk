/// \ingroup base
/// \class ttk::RoadRefinement
/// \author Jonas Lukasczyk <jl@jluk.de>
/// \date 1.09.2019
///
/// TODO

#pragma once

// ttk common includes
#include <Debug.h>
#include <math.h>

namespace ttk {

  class RoadRefinement : virtual public Debug {

  public:
    RoadRefinement() {
      this->setDebugMsgPrefix(
        "RoadRefinement"); // inherited from Debug: prefix will be printed at
                           // the beginning of every msg
    };
    ~RoadRefinement(){};

    static inline double DegreeToRadian(const float &angle) {
      return M_PI * angle / 180.0;
    }

    static inline double RadianToDegree(const float &radian) {
      return 180.0 * radian / M_PI;
    }

    // common functions for geographic calculation
    static inline float distance4TwoPoints(float latStartP,
                                           float longStartP,
                                           float latEndP,
                                           float longEndP) {
      const static float EarthRadiusMeter = 6371e3;
      float diffLat = DegreeToRadian(latEndP - latStartP);
      float diffLong = DegreeToRadian(longEndP - longStartP);

      float a = sin(diffLat / 2) * sin(diffLat / 2)
                + cos(DegreeToRadian(latStartP)) * cos(DegreeToRadian(latEndP))
                    * sin(diffLong / 2) * sin(diffLong / 2);
      float computation = 2 * atan2(sqrt(a), sqrt(1 - a));
      float result = computation * EarthRadiusMeter;

      return result;
    }

    // the bearing should be calculated p1 -p0
    static inline float bearing4TwoPoints(float latStartP,
                                          float longStartP,
                                          float latEndP,
                                          float longEndP) {
      auto latSPRadian = DegreeToRadian(latStartP);
      auto latEPRadian = DegreeToRadian(latEndP);
      auto diffLong = DegreeToRadian(longEndP - longStartP);

      auto bearing
        = atan2(sin(diffLong) * cos(latEPRadian),
                cos(latSPRadian) * sin(latEPRadian)
                  - sin(latSPRadian) * cos(latEPRadian) * cos(diffLong));
      bearing = RadianToDegree(bearing);
      bearing = (((int)bearing + 360) % 360);
      return bearing;
    }

    static inline double *destinationPointGivenDistandBearingfromStartPoint(
      float latStartP, float longStartP, float bearing, float distance) {
      static double result[2];
      auto latSPRadian = DegreeToRadian(latStartP);
      auto longSPRadian = DegreeToRadian(longStartP);
      bearing = DegreeToRadian(bearing);

      const float EarthRadiusMeter = 6371e3;
      double ratio = distance / EarthRadiusMeter;

      double latNewPoint = asin(sin(latSPRadian) * cos(ratio)
                                + cos(latSPRadian) * sin(ratio) * cos(bearing));

      double longNewPoint
        = longSPRadian
          + atan2(sin(bearing) * sin(ratio) * cos(latSPRadian),
                  cos(ratio) - sin(latSPRadian) * sin(latNewPoint));

      latNewPoint = RadianToDegree(latNewPoint);
      longNewPoint = RadianToDegree(longNewPoint);

      result[0] = latNewPoint;
      result[1] = longNewPoint;
      return result;
    }

    template <class idType>
    int getNumberOfRefinedPointsAndEdges(const float *originalPointCoords,
                                         const idType *originalConnectivityList,
                                         const idType &nOriginalCells,
                                         int &nRefinedPoints,
                                         int &nRefinedEdges,
                                         const float segmentUnit) const;

    template <class idType>
    int refineRoads(const float *originalPointCoords,
                    const idType *originalConnectivityList,
                    const float refinementUnit,
                    float *refinedPointCoords,
                    idType *refinedConnectivityList,
                    const idType &nOriginalCells,
                    const idType &nOriginalPoints) const;
  };
} // namespace ttk

template <class idType>
int ttk::RoadRefinement::getNumberOfRefinedPointsAndEdges(
  const float *originalPointCoords,
  const idType *originalConnectivityList,
  const idType &nOriginalCells,
  int &nRefinedPoints,
  int &nRefinedEdges,
  const float segmentUnit) const {
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
    int nSegments = ceil(distance / segmentUnit);
    if(nSegments > 0) {
      // not dirty edge
      int nNewPoints = nSegments - 1;
      int nNewEdges = nSegments - 1;

      nRefinedPoints += nNewPoints;
      nRefinedEdges += nNewEdges;
    } else {
      // std::cout << "dirty edge!!!!!!!!!!!!!!!!!!!!" << std::endl;
      nRefinedEdges--;
    }
  }

  return 0;
}

template <class idType>
int ttk::RoadRefinement::refineRoads(const float *originalPointCoords,
                                     const idType *originalConnectivityList,
                                     const float refinementUnit,
                                     float *refinedPointCoords,
                                     idType *refinedConnectivityList,
                                     const idType &nOriginalCells,
                                     const idType &nOriginalPoints) const {
  int nTotalPoints = nOriginalPoints;
  int nTotalEdges = 0;

  for(int i = 0, j = nOriginalPoints * 3; i < j; i++) {
    refinedPointCoords[i] = originalPointCoords[i];
  }

  for(int cellIndex = 0; cellIndex < nOriginalCells; cellIndex++) {
    auto p1 = originalConnectivityList[cellIndex * 2];
    auto p2 = originalConnectivityList[cellIndex * 2 + 1];
    if(p1 == p2)
      continue;

    float p1Lat = originalPointCoords[p1 * 3 + 1];
    float p1Lon = originalPointCoords[p1 * 3];
    float p2Lat = originalPointCoords[p2 * 3 + 1];
    float p2Lon = originalPointCoords[p2 * 3];

    auto distance
      = ttk::RoadRefinement::distance4TwoPoints(p1Lat, p1Lon, p2Lat, p2Lon);
    int nSegments = ceil(distance / refinementUnit);

    if(nSegments == 1) {
      // std::cout << "I do not need to segment! " << std::endl;
      // if there is not new points and edges, insert the orignal edges to
      // the cell array
      refinedConnectivityList[nTotalEdges * 2] = p1;
      refinedConnectivityList[nTotalEdges * 2 + 1] = p2;
      nTotalEdges++;
      continue;
    }

    if(nSegments > 1) {
      // std::cout << "I need to segment! " << nSegments << std::endl;
      float segmentDist = distance / nSegments;
      auto bearing
        = ttk::RoadRefinement::bearing4TwoPoints(p1Lat, p1Lon, p2Lat, p2Lon);

      for(int index4seg = 0; index4seg < nSegments; index4seg++) {

        int c1 = index4seg == 0 ? p1 : nTotalPoints - 1;
        int c2 = index4seg == nSegments - 1 ? p2 : nTotalPoints;

        refinedConnectivityList[nTotalEdges * 2] = c1;
        refinedConnectivityList[nTotalEdges * 2 + 1] = c2;

        float distance2startpoint = (index4seg + 1) * segmentDist;
        double *newPoint = ttk::RoadRefinement::
          destinationPointGivenDistandBearingfromStartPoint(
            p1Lat, p1Lon, bearing, distance2startpoint);

        if(c2 != p2) {
          refinedPointCoords[c2 * 3] = newPoint[1];
          refinedPointCoords[c2 * 3 + 1] = newPoint[0];
          refinedPointCoords[c2 * 3 + 2] = 0;

          nTotalPoints += 1;
        }

        nTotalEdges += 1;
      }
    }
    }

  std::cout << "nTotalPoints: " << nTotalPoints << std::endl;
  std::cout << "nTotalEdges: " << nTotalEdges << std::endl;

  return 0;
}