/// \ingroup base
/// \class ttk::RoadDataConverter
/// \author Jonas Lukasczyk <jl@jluk.de>
/// \date 1.09.2019
///
/// TODO

#pragma once

// ttk common includes
#include <Debug.h>

namespace ttk {

    class RoadDataConverter : virtual public Debug {

        public:
            RoadDataConverter(){
                this->setDebugMsgPrefix("RoadDataConverter"); // inherited from Debug: prefix will be printed at the beginning of every msg
            };
            ~RoadDataConverter(){};

            int* getNumberOfPointsandEdges(const std::string &path);
            int parsePointCoords(const std::string &path, float *pointCoords, long long int* cellIds);
    };
}
