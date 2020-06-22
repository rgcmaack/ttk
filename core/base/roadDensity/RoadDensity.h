/// \ingroup base
/// \class ttk::RoadDensity
/// \author Jonas Lukasczyk <jl@jluk.de>
/// \date 1.09.2019
///
/// TODO

#pragma once

// ttk common includes
#include <Debug.h>
#include <RoadRefinement.h>
#include <Triangulation.h>
#include <limits>
#include <math.h>

namespace ttk {

  class RoadDensity : virtual public Debug {

  public:
    RoadDensity() {
      this->setDebugMsgPrefix(
        "RoadDensity"); // inherited from Debug: prefix will be
                        // printed at the beginning of every msg
    };
    ~RoadDensity(){};

    // template <class idType>
    // int CalculateRRoadPointsEventSamples(
    //   std::vector<int> &rRoadPointtoEventSampleDictionary,
    //   int &nEventPoints,
    //   idType *closetPoints) const;

    template <class idType>
    int CalculateRoadWeight(idType *rRoadPointtoEventSampleDictionary,
                            const int &refinedRoadPointsNum,
                            ttk::Triangulation *triangulation,
                            idType *rRoadPointsCoors,
                            idType *rRoadPointsWeight,
                            idType &KernelBandwidth) const;
  };
} // namespace ttk
////eill remove later

int getNeighbors(std::vector<int> &visitedVertics,
                 std::vector<float> &visitedVShortestDistance,
                 int &rRoadPointIndex,
                 ttk::Triangulation *triangulation,
                 float *rRoadPointsCoors,
                 float KernelBandwidth) {
  std::vector<int> candidateVertics;
  std::vector<float> pilotDistance;

  candidateVertics.push_back(rRoadPointIndex);
  pilotDistance.push_back(0);

  while(!candidateVertics.empty()) {
    // get the smallest item in candidates
    auto minDist
      = *std::min_element(pilotDistance.begin(), pilotDistance.end());
    // std::cout << "minDist: " << minDist << std::endl;
    if(minDist > KernelBandwidth) {
      break;
    }

    auto minumDistIndex
      = std::min_element(pilotDistance.begin(), pilotDistance.end())
        - pilotDistance.begin();

    auto verticIndex2visited = candidateVertics[minumDistIndex] * 1;
    

    // push it the the visited array
    visitedVertics.push_back(verticIndex2visited);
    visitedVShortestDistance.push_back(minDist);

    // delete the item in the candidates
    candidateVertics.erase(candidateVertics.begin() + minumDistIndex);
    pilotDistance.erase(pilotDistance.begin() + minumDistIndex);

    // calculate the new neighbors based on the visited item
    std::vector<int> neighborVertics;
    std::vector<float> distancefromStartpointArray;
    auto nNeighbours
      = triangulation->getVertexNeighborNumber(verticIndex2visited);
    for(int j = 0; j < nNeighbours; j++) {
      ttk::SimplexId neighbourId = 0;
      triangulation->getVertexNeighbor(verticIndex2visited, j, neighbourId);

      // check the neighbor is not in the visited array
      if(std::find(visitedVertics.begin(), visitedVertics.end(), neighbourId)
         != visitedVertics.end()) {
        continue;
      }

      // get the distance betwen the visitedpoint to the neighbor
      auto distance = ttk::RoadRefinement::distance4TwoPoints(
        rRoadPointsCoors[verticIndex2visited * 3 + 1],
        rRoadPointsCoors[verticIndex2visited * 3],
        rRoadPointsCoors[neighbourId * 3 + 1],
        rRoadPointsCoors[neighbourId * 3]);
      float newDist = distance + minDist;

      // push or update the neighbor in the candidate pool
      std::vector<int>::iterator itr = std::find(
        candidateVertics.begin(), candidateVertics.end(), neighbourId);
      if(itr != candidateVertics.cend()) {
        // get the index of the neighbor in the candidate pool
        auto indexinCandidatePool
          = std::distance(candidateVertics.begin(), itr);

        // update the newDist fot the neighbor in the candidate pool
        if(pilotDistance[indexinCandidatePool] > newDist) {
          pilotDistance[indexinCandidatePool] = newDist;
        }
      } else {
        // push the neighbor to the candidate pool
        candidateVertics.push_back(neighbourId);
        pilotDistance.push_back(newDist);
      }
    }
  }

  return 0;
}

int kernelDensity(std::vector<int> &neighbors,
                  std::vector<float> &shortestDistances,
                  float &eventSamplesValue,
                  float *rRoadPointsWeight,
                  float KernelBandwidth) {
  int neighborsNumincludingitself = shortestDistances.size();
  for(int i = 0; i < neighborsNumincludingitself; i++) {
    auto vertexIndex = neighbors[i];
    auto vertexDist = shortestDistances[i];

    float vertexDensity
      = eventSamplesValue * (1 - vertexDist / KernelBandwidth);
    
    rRoadPointsWeight[vertexIndex] += vertexDensity;

  }

  return 0;
}

// template <class idType>
// int ttk::RoadDensity::CalculateRRoadPointsEventSamples(
//   std::vector<int> &rRoadPointtoEventSampleDictionary,
//   int &nEventPoints,
//   idType *closetPoints) const {

//   for(int eventPointIndex = 0; eventPointIndex < nEventPoints;
//       eventPointIndex++) {
//     int refinedRoadPointIndex = closetPoints->GetTuple1(eventPointIndex);
//     if(refinedRoadPointIndex != -1) {
//       rRoadPointtoEventSampleDictionary[refinedRoadPointIndex] += 1;
//     }
//   }
//   return 0;
// }

template <class idType>
int ttk::RoadDensity::CalculateRoadWeight(
  idType *rRoadPointtoEventSampleDictionary,
  const int &refinedRoadPointsNum,
  ttk::Triangulation *triangulation,
  idType *rRoadPointsCoors,
  idType *rRoadPointsWeight,
  idType &KernelBandwidth) const {

  for(int i = 0; i < refinedRoadPointsNum; i ++) {
    rRoadPointsWeight[i] = 0;
  }
    // int count4nonzero = 0;
    for(int rRoadPointIndex = 0; rRoadPointIndex < refinedRoadPointsNum;
        rRoadPointIndex++) {
      if(rRoadPointtoEventSampleDictionary[rRoadPointIndex] != 0) {
        // count4nonzero++;
        auto eventSamplesValue
          = rRoadPointtoEventSampleDictionary[rRoadPointIndex];
        
        std::vector<int> neighbors;
        std::vector<float> shortestDistances;
        getNeighbors(neighbors, shortestDistances, rRoadPointIndex,
                     triangulation, rRoadPointsCoors, KernelBandwidth);

        std::vector<float> neighborsDensity(shortestDistances.size());
        kernelDensity(neighbors, shortestDistances, eventSamplesValue,
                      rRoadPointsWeight, KernelBandwidth);
      }
    }
  // std::cout << "count4nonzero " << count4nonzero << std::endl;

  return 1;
}
