#include <EventDataConverter.h>
#include <algorithm>
#include <iterator>
#include <typeinfo>

int ttk::EventDataConverter::getNumberOfPoints(const std::string &path) {
  this->printMsg("Code in base layer: getNumberOfPoints");

  int eventNum = 0;

  std::ifstream infile(path);
  for(std::string line; getline(infile, line);) {
    std::string::iterator end_pos = std::remove(line.begin(), line.end(), ' ');
    line.erase(end_pos, line.end());
    if(line.compare(0, 2, "{\"") != 0)
      continue;

    // get the valid line and extact the lat and long info to remove the dirty
    // events.
    std::string parsed;
    std::stringstream input_stringstream(line);
    bool validLong = false, validLat = false;
    while(std::getline(input_stringstream, parsed, ',')) {
      std::string longit = "\"Longitude\":", lat = "\"Latitude\":";
      if(parsed.compare(0, longit.length(), longit) == 0) {
        std::string longitude = parsed.substr(longit.length());
        if(std::stof(longitude) == 0.0)
          continue;
        validLong = true;
      }
      if(parsed.compare(0, lat.length(), lat) == 0) {
        std::string latitude = parsed.substr(lat.length());
        latitude.pop_back();
        if(std::stof(latitude) == 0.0)
          continue;
        validLat = true;
      }
    }
    if(validLat and validLong) {
      eventNum++;
    }
  }

  //   this->printMsg("eventNum: " + std::to_string(eventNum));
  return eventNum;
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

    if(line.compare(0, 2, "{\"") != 0)
      continue;

    std::string parsed;
    std::stringstream input_stringstream(line);
    bool validLong = false, validLat = false;
    while(std::getline(input_stringstream, parsed, ',')) {
      std::string offense = "\"Offense\":", longit = "\"Longitude\":",
                  lat = "\"Latitude\":";

      if(parsed.compare(0, longit.length(), longit) == 0) {
        std::string longitude = parsed.substr(longit.length());
        // std::ostringstream out;
        // out << std::setprecision(8) << std::stof(longitude);
        // float longitF = std::stof(out.str());
        if(std::stof(longitude) == 0.0)
          continue;
        validLong = true;
        pointCoords[3 * pointIndex] = std::stof(longitude);

      } else if(parsed.compare(0, lat.length(), lat) == 0) {
        std::string latitude = parsed.substr(lat.length());
        latitude.pop_back();
        // std::ostringstream out;
        // out << std::setprecision(8) << std::stof(latitude);
        // float latF = std::stof(out.str());
        if(std::stof(latitude) == 0.0)
          continue;
        validLat = true;
        pointCoords[3 * pointIndex + 1] = std::stof(latitude);
        pointCoords[3 * pointIndex + 2] = 0;
      } else if(parsed.compare(0, offense.length(), offense) == 0) {
        std::string offenseCategory = parsed.substr(offense.length());
        // TODO: Lookup
        auto target = std::find(categoryDictionary.begin(),
                                categoryDictionary.end(), offenseCategory);
        if(target != categoryDictionary.end()) {
          // if it's in the dictionary, get the index and store to the
          // categoryindex array
          auto target_index = std::distance(categoryDictionary.begin(), target);
          categoryIndex[pointIndex] = target_index;
        } else {
          // if it's not in the dictionary, add to the dictionary , add the last
          // index of the dictionary to the categoryindex array
          categoryDictionary.push_back(offenseCategory);
          categoryIndex[pointIndex] = categoryDictionary.size() - 1;
        }
      }
    }

    if(validLat and validLong) {
      pointIndex++;
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