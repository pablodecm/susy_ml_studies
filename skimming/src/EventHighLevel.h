#ifndef EventHighLevel_H
#define EventHighLevel_H

/**
 * Container class with derived high level features of the event
 */
class EventHighLevel {
  public:
    float dilept_inv_mass;
    float dilepton_MT2;
    float jets_ht;
    float max_CSV;

    EventHighLevel() {}

    ClassDef(EventHighLevel, 1);

};

#endif
