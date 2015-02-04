#ifndef EventData_H
#define EventData_H

/**
 * Container class to keep useful data for each event and save to skimmed TTree
 * if required.
 */
class EventData {
  public:
    int run_number;
    int lumi_block;
    int event_number;
    int process_ID;  

    EventData () {}

    ClassDef(EventData, 1);
    
};

#endif
