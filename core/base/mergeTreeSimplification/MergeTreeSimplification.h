#pragma once

#include <Debug.h>

namespace ttk {

    class MergeTreeSimplification : virtual public Debug {

        public:

            MergeTreeSimplification(){
                this->setDebugMsgPrefix("MergeTreeSimplification");
            };
            ~MergeTreeSimplification(){};

    };
}