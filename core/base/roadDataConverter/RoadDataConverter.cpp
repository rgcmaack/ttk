#include <RoadDataConverter.h>
#include <algorithm>
#include <typeinfo>
#include <unordered_map>
#include <unordered_set>

std::vector<std::string> splitStr(std::string str2split, char letter4split) {
  std::stringstream splitFun(str2split);
  std::string segment;
  std::vector<std::string> seglist;
  while(std::getline(splitFun, segment, letter4split)) {
    seglist.push_back(segment);
  }
  return seglist;
}

int ttk::RoadDataConverter::getNumberOfPointsandEdges(const std::string &path,
                                                      int &npoints,
                                                      int &nedges) {
  this->printMsg("Code in base layer: getNumberOfPointsandEdges");
  std::unordered_set<std::string> pointSet;

  std::ifstream infile(path);
  for(std::string line; getline(infile, line);) {
    int curPoints = 0;
    int nonDuplicatePoints = 0;

    std::string::iterator end_pos = std::remove(line.begin(), line.end(), ' ');
    line.erase(end_pos, line.end());

    if(line.compare(0, 2, "{\"") != 0)
      continue;

    // get coordinates for each line
    std::size_t coordinatesInfor = line.find("\"coordinates\"");
    std::string extractedCoor = line.substr(coordinatesInfor);

    std::string parsed4coor;
    std::stringstream input_stringstream(extractedCoor);

    while(std::getline(input_stringstream, parsed4coor, '[')) {
      if(parsed4coor.compare(0, 1, "-") != 0)
        continue;
      curPoints++;
      // get each point
      std::vector<std::string> seglist4latlong = splitStr(parsed4coor, ']');

      std::string longlatstr = seglist4latlong[0];

      // check whether the point is already in the pointCoor
      /////if exist, do nothing
      if(pointSet.find(longlatstr) != pointSet.end()) {
        continue;
      }

      //// if not exist, add 1
      pointSet.insert(longlatstr);
      nonDuplicatePoints++;
    }
    npoints += nonDuplicatePoints;
    nedges += curPoints - 1;
  }
  std::cout << "npoints in getting: " << npoints << std::endl;
  std::cout << "nedges in getting: " << nedges << std::endl;
  return 1;
};

int ttk::RoadDataConverter::parsePointCoords(
  const std::string &path,
  float *pointCoords,
  long long int *cellConnectivityData,
  unsigned char *category4edgeArr,
  std::vector<std::string> &categoryDictionary) {
  this->printMsg("Code in base layer: parsePointCoords");
  std::unordered_map<std::string, int> pointMap;

  // record point index and edge index
  int npoints = 0, nedges = 0;

  std::ifstream infile(path);
  for(std::string line; getline(infile, line);) {
    // record the num of points in each line
    int nonDuplicatePoints = 0;

    std::string::iterator end_pos = std::remove(line.begin(), line.end(), ' ');
    line.erase(end_pos, line.end());

    if(line.compare(0, 2, "{\"") != 0)
      continue;

    // get property of road
    std::size_t osmIdInfor = line.find("\"highway\":");
    std::string extractedOsmId = line.substr(osmIdInfor);

    std::vector<std::string> seglist4osmid = splitStr(extractedOsmId, ',');
    std::string osmId = seglist4osmid[0];

    // get coordinates
    std::size_t coordinatesInfor = line.find("\"coordinates\"");
    std::string extractedCoor = line.substr(coordinatesInfor);

    std::string parsed4coor;
    std::stringstream input_stringstream(extractedCoor);

    std::vector<int> points4line;
    while(std::getline(input_stringstream, parsed4coor, '[')) {
      if(parsed4coor.compare(0, 1, "-") != 0)
        continue;

      // get each point
      std::vector<std::string> seglist4latlong = splitStr(parsed4coor, ']');

      std::string longlatstr = seglist4latlong[0];
      std::vector<std::string> longlatPair = splitStr(longlatstr, ',');

      // check whether the point is already existing.
      float longitude = std::stof(longlatPair[0]);
      float latitude = std::stof(longlatPair[1]);
      float hashCode = latitude * 31.0 - longitude;
      //// if exist, do nothing to the pointCoor array
      if(pointMap.find(longlatstr) != pointMap.end()) {
        points4line.push_back(pointMap[longlatstr]);
        continue;
      }

      //// if not exist, add 1 and add to the pointCoor array
      // insert point to pointCoords
      pointCoords[3 * (npoints + nonDuplicatePoints)] = longitude;
      pointCoords[3 * (npoints + nonDuplicatePoints) + 1] = latitude;
      pointCoords[3 * (npoints + nonDuplicatePoints) + 2] = 0;
      pointMap[longlatstr] = npoints + nonDuplicatePoints;
      points4line.push_back(npoints + nonDuplicatePoints);
      nonDuplicatePoints++;
    }

    // insert cell/line to cellConnectivityData and insert osmid to
    // category4edgeArr
    // TODO: Lookup
    int target_index = -1;
    auto target
      = std::find(categoryDictionary.begin(), categoryDictionary.end(), osmId);
    if(target != categoryDictionary.end()) {
      // if it's in the dictionary, get the index and store to the
      // categoryindex array
      target_index = std::distance(categoryDictionary.begin(), target);
    } else {
      // if it's not in the dictionary, add to the dictionary , add the last
      // index of the dictionary to the categoryindex array
      categoryDictionary.push_back(osmId);
      target_index = categoryDictionary.size() - 1;
    }

    ////fill the cell connectivity for each line
    for(int index = 0; index < points4line.size() - 1; index++) {
      // fill category
      category4edgeArr[nedges + index] = target_index;
      // fill cellconnectivity
      cellConnectivityData[(nedges + index) * 2] = points4line[index];
      cellConnectivityData[(nedges + index) * 2 + 1] = points4line[index + 1];
    }

    npoints += nonDuplicatePoints;
    nedges += points4line.size() - 1;
  }
  std::cout << "npoints in parseing: " << npoints << std::endl;
  std::cout << "nedges in parseing: " << nedges << std::endl;

  return 1;
};
