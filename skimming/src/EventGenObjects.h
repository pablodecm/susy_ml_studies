#ifndef EventGenObjects_H
#define EventGenObjects_H

/**
 * Container class with the generation level physical object information of
 * interest of the event
 */
class EventGenObjects {
  public:
    float gen_stop_mass;
    float gen_chi0_mass;
    float gen_chargino_mass;


    EventGenObjects() {}

    ClassDef(EventGenObjects, 1);

};

#endif
