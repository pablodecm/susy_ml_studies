#ifndef EventTopology_H
#define EventTopology_H

/**
 * Container class to keep the channel and other variables about the final
 * state topology of the event.
 */
class EventTopology {
  public:
    int channel;
    int n_lept;
    int n_elec;
    int n_muon; 
    int n_jet;
    int n_b_jet;
  
    EventTopology() {}

    ClassDef(EventTopology,1);
};




#endif
