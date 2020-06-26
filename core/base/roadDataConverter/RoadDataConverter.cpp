#include <RoadDataConverter.h>
#include <algorithm>
#include <typeinfo>

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

  std::ifstream infile(path);
  for(std::string line; getline(infile, line);) {
    int curPoints = 0;

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
    }
    npoints += curPoints;
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
  std::string *category4edgeArr) {
  this->printMsg("Code in base layer: parsePointCoords");

  // record point index and edge index
  int npoints = 0, nedges = 0;

  std::ifstream infile(path);
  for(std::string line; getline(infile, line);) {
    // record the num of points in each line
    int curPoints = 0;

    std::string::iterator end_pos = std::remove(line.begin(), line.end(), ' ');
    line.erase(end_pos, line.end());

    if(line.compare(0, 2, "{\"") != 0)
      continue;

    // get osmId
    std::size_t osmIdInfor = line.find("\"osm_id\":");
    std::string extractedOsmId = line.substr(osmIdInfor);

    std::vector<std::string> seglist4osmid = splitStr(extractedOsmId, ',');
    std::string osmId = seglist4osmid[0];

    // get coordinates
    std::size_t coordinatesInfor = line.find("\"coordinates\"");
    std::string extractedCoor = line.substr(coordinatesInfor);

    std::string parsed4coor;
    std::stringstream input_stringstream(extractedCoor);

    while(std::getline(input_stringstream, parsed4coor, '[')) {
      if(parsed4coor.compare(0, 1, "-") != 0)
        continue;

      // get each point
      std::vector<std::string> seglist4latlong = splitStr(parsed4coor, ']');

      std::string longlatstr = seglist4latlong[0];
      std::vector<std::string> longlatPair = splitStr(longlatstr, ',');

      // insert point to pointCoords
      pointCoords[3 * (npoints + curPoints)] = std::stof(longlatPair[0]);
      pointCoords[3 * (npoints + curPoints) + 1] = std::stof(longlatPair[1]);
      pointCoords[3 * (npoints + curPoints) + 2] = 0;

      curPoints++;
    }

    // insert cell/line to cellConnectivityData and insert osmid to
    // category4edgeArr
    int startIndexofpoint4cell = npoints * 1;
    for(int index = 0; index < curPoints - 1; index++) {
      // fill category
      category4edgeArr[nedges + index] = osmId;
      // fill cellconnectivity
      cellConnectivityData[(nedges + index) * 2]
        = startIndexofpoint4cell + index;
      cellConnectivityData[(nedges + index) * 2 + 1]
        = startIndexofpoint4cell + index + 1;
    }

    npoints += curPoints;
    nedges += curPoints - 1;
  }
  std::cout << "npoints in parseing: " << npoints << std::endl;
  std::cout << "nedges in parseing: " << nedges << std::endl;

  return 1;
};
