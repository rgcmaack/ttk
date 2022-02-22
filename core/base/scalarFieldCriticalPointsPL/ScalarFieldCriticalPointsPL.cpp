#include <ScalarFieldCriticalPointsPL.h>

ttk::ScalarFieldCriticalPointsPL::ScalarFieldCriticalPointsPL() {
  this->setDebugMsgPrefix("ScalarFieldCriticalPointsPL");
}

void ttk::ScalarFieldCriticalPointsPL::displayStats() {

  SimplexId minimumNumber = 0, maximumNumber = 0, saddleNumber = 0,
            oneSaddleNumber = 0, twoSaddleNumber = 0, monkeySaddleNumber = 0;

  if(debugLevel_ >= (int)debug::Priority::INFO) {
    if(dimension_ == 3) {
      for(size_t i = 0; i < criticalPoints_->size(); i++) {
        switch((*criticalPoints_)[i].second) {

          case(char)(CriticalType::Local_minimum):
            minimumNumber++;
            break;

          case(char)(CriticalType::Saddle1):
            oneSaddleNumber++;
            break;

          case(char)(CriticalType::Saddle2):
            twoSaddleNumber++;
            break;

          case(char)(CriticalType::Local_maximum):
            maximumNumber++;
            break;

          case(char)(CriticalType::Degenerate):
            monkeySaddleNumber++;
            break;
        }
      }
    } else if(dimension_ == 2) {
      for(size_t i = 0; i < criticalPoints_->size(); i++) {
        switch((*criticalPoints_)[i].second) {

          case(char)(CriticalType::Local_minimum):
            minimumNumber++;
            break;

          case(char)(CriticalType::Saddle1):
            saddleNumber++;
            break;

          case(char)(CriticalType::Local_maximum):
            maximumNumber++;
            break;

          case(char)(CriticalType::Degenerate):
            monkeySaddleNumber++;
            break;
        }
      }
    }

    {
      std::vector<std::vector<std::string>> stats;
      stats.push_back({"  #Minima", std::to_string(minimumNumber)});
      if(dimension_ == 3) {
        stats.push_back({"  #1-saddles", std::to_string(oneSaddleNumber)});
        stats.push_back({"  #2-saddles", std::to_string(twoSaddleNumber)});
      }
      if(dimension_ == 2) {
        stats.push_back({"  #Saddles", std::to_string(saddleNumber)});
      }
      stats.push_back({"  #Multi-saddles", std::to_string(monkeySaddleNumber)});
      stats.push_back({"  #Maxima", std::to_string(maximumNumber)});

      printMsg(stats);
    }
  }
}