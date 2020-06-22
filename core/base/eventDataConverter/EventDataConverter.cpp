#include <EventDataConverter.h>
#include <algorithm>
#include <iterator>
#include <typeinfo>

int ttk::EventDataConverter::getNumberOfPoints(const std::string &path) {
  this->printMsg("Code in base layer: getNumberOfPoints");

  int eventNum = 0;

  std::ifstream infile(path);
  for(std::string line; getline(infile, line);) {
    // this->printMsg("'" + line + "'");
    std::string::iterator end_pos = std::remove(line.begin(), line.end(), ' ');
    line.erase(end_pos, line.end());
    // this->printMsg("after:' " + line + "'");
    if(line.compare(0, 1, "{") == 0) {
      eventNum++;
    }
  }

  //   this->printMsg("eventNum: " + std::to_string(eventNum));
  return eventNum - 1;
};

int ttk::EventDataConverter::parsePointCoords(
  const std::string &path,
  float *pointCoords,
  unsigned char *categoryIndex,
  std::vector<std::string> &categoryDictionary) {
  this->printMsg("Code in base layer: parsePointCoords");

  int pointIndex = 0;
  std::ifstream infile(path);
  for(std::string line; getline(infile, line);) {

    std::string::iterator end_pos = std::remove(line.begin(), line.end(), ' ');
    line.erase(end_pos, line.end());

    if(!(line.compare("[") == 0 || line.compare("},") == 0
         || line.compare("}") == 0 || line.compare("]") == 0
         || line.compare("{") == 0)) {
      //   this->printMsg("We have to parse: " + line);
      //// split line to two parts, key and value
      std::stringstream test(line);
      std::string segment;
      std::vector<std::string> seglist;

      while(std::getline(test, segment, ':')) {
        seglist.push_back(segment);
      }

      std::string key = seglist[0];
      std::string value = seglist[1];
      //   this->printMsg("seglist[0]: " + key);
      //   this->printMsg("seglist[1]: " + value);

      //// insert lat and long to buffer
      if(key.compare("\"latitude\"") == 0) {
        pointCoords[3 * pointIndex + 1] = std::stof(seglist[1]);
        // this->printMsg("key is latitude: " + key);
        // this->printMsg("seglist[1]: " + std::stof(seglist[1]));
      } else if(key.compare("\"longitude\"") == 0) {
        pointCoords[3 * pointIndex] = std::stof(seglist[1]);
        pointCoords[3 * pointIndex + 2] = 0;
      } else if(key.compare("\"category\"") == 0) {
        // std::cout << "value " << value << std::endl;
        // TODO: Lookup
        auto target = std::find(
          categoryDictionary.begin(), categoryDictionary.end(), value);
        if(target != categoryDictionary.end()) {
          // if it's in the dictionary, get the index and store to the categoryindex array
          auto target_index = std::distance(categoryDictionary.begin(), target);
          categoryIndex[pointIndex] = target_index;
        } else {
          // if it's not in the dictionary, add to the dictionary , add the last index of the dictionary to the categoryindex array
          categoryDictionary.push_back(value);
          categoryIndex[pointIndex] = categoryDictionary.size() - 1;
        }

        pointIndex++;
      }
    }
  }

  //   pointCoords[0] = p0x;
  //   pointCoords[1] = p0y;
  //   pointCoords[2] = p0z;
  //   pointCoords[3] = p1x;
  //   pointCoords[4] = p1y;
  //   pointCoords[5] = p1z;
  //   pointCoords[0] = p0x;

  return 1;
};