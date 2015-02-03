#ifndef EventData_H
#define EventData_H


class EventData {
  public:
    int channel;
    int nJets;
    int nBJets;
    float dilMass;
    float met_Et;
    float htJets;
  
    EventData() {}

    ClassDef(EventData,1);
};




#endif
