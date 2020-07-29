/// \ingroup base
/// \class ttk::EventDataConverter
/// \author Jonas Lukasczyk <jl@jluk.de>
/// \date 1.09.2019
///
/// TODO

#pragma once

// ttk common includes
#include <Debug.h>

namespace ttk {

    class EventDataConverter : virtual public Debug {

        public:
            EventDataConverter(){
                this->setDebugMsgPrefix("EventDataConverter"); // inherited from Debug: prefix will be printed at the beginning of every msg
            };
            ~EventDataConverter(){};

            int getNumberOfPoints( const std::string& path );
            int parsePointCoords(const std::string &path,
                                 float *pointCoords,
                                 unsigned char *categoryIndex,
                                 unsigned char *yearIndex,
                                 std::vector<std::string> &categoryDictionary,
                                 std::vector<std::string> &yearDictionary);
    };
}
