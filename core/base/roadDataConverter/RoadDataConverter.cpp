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

int *ttk::RoadDataConverter::getNumberOfPointsandEdges(
  const std::string &path) {
  this->printMsg("Code in base layer: getNumberOfPointsandEdges");

  static int pointsandEdges[2] = {0, 0};

  std::ifstream infile(path);
  for(std::string line; getline(infile, line);) {
    std::string::iterator end_pos = std::remove(line.begin(), line.end(), ' ');
    line.erase(end_pos, line.end());

    if(line[0] == '{') {
      int lines4obj = 0;
      int points4obj = 0;
      std::size_t coordinatesInfor = line.find("\"coordinates\"");
      std::string extractedCoor = line.substr(coordinatesInfor);

      //// split for coorsArray
      std::vector<std::string> seglist = splitStr(extractedCoor, '[');

      for(std::vector<std::string>::iterator it = seglist.begin();
          it != seglist.end(); ++it) {
        // this->printMsg("tokens:" + *it);

        std::string segmentStr = *it;
        // std::cout << segmentStr[0] << std::endl;
        if(segmentStr[0] == '-') {
          // split lat and long
          std::vector<std::string> seglist4latlong = splitStr(segmentStr, ']');
          // this->printMsg("latlongs: " + seglist4latlong[0]);
          points4obj += 1;
        }
      }

      lines4obj = points4obj - 1;
      if(points4obj < 2) {
        lines4obj = 0;
        this->printMsg("no lines!!!!!");
      }

      pointsandEdges[0] += points4obj;
      pointsandEdges[1] += lines4obj;
    }
  }

  // std::cout << pointsandEdges << std::endl;
  // this->printMsg("points: " + std::to_string(pointsandEdges[0]));
  // this->printMsg("Edges: " + std::to_string(pointsandEdges[1]));
  return pointsandEdges;
};

int ttk::RoadDataConverter::parsePointCoords(const std::string &path,
                                             float *pointCoords,
                                             long long int *cellIds) {
  this->printMsg("Code in base layer: parsePointCoords");

  std::ifstream infile(path);
  int pointNum = 0;
  int lineNum = 0;

  for(std::string line; getline(infile, line);) {

    std::string::iterator end_pos = std::remove(line.begin(), line.end(), ' ');
    line.erase(end_pos, line.end());

    if(line[0] == '{') {
      int numberOfPointsPerLine = 0;
      int numberOfEdgesPerLine = 0;
      std::size_t coordinatesInfor = line.find("\"coordinates\"");
      std::string extractedCoor = line.substr(coordinatesInfor);

      //// split for coorsArray
      std::vector<std::string> seglist = splitStr(extractedCoor, '[');

      // std::cout << " -------------- "<< std::endl;

      for(std::vector<std::string>::iterator it = seglist.begin();
          it != seglist.end(); ++it) {
        std::string segmentStr = *it;
        if(segmentStr[0] == '-') {
          // split lat and long
          std::vector<std::string> seglist4latlong = splitStr(segmentStr, ']');

          std::string latlongstr = seglist4latlong[0];
          std::vector<std::string> latlongPair = splitStr(latlongstr, ',');

          //// insert point to pointCoords
          pointCoords[pointNum * 3 + 3 * numberOfPointsPerLine]
            = std::stof(latlongPair[0]);
          pointCoords[pointNum * 3 + 3 * numberOfPointsPerLine + 1]
            = std::stof(latlongPair[1]);
          pointCoords[pointNum * 3 + 3 * numberOfPointsPerLine + 2] = 0;

          numberOfPointsPerLine++;
        }
      }
      numberOfEdgesPerLine = numberOfPointsPerLine - 1;
      // //// insert cell/line to cellID
      // std::cout << "-----------------------------" << std::endl;
      int startIndexofpoint4cell = pointNum * 1;
      for(int index = 0; index < numberOfEdgesPerLine; index++) {
        cellIds[(lineNum + index) * 3] = 2;
        cellIds[(lineNum + index) * 3 + 1] = startIndexofpoint4cell + index;
        cellIds[(lineNum + index) * 3 + 2] = startIndexofpoint4cell + index + 1;
      }

      lineNum += numberOfEdgesPerLine;
      pointNum += numberOfPointsPerLine;

    } // end of the correct line to insert point and edge
  } // end of the line
  return 1;
};
