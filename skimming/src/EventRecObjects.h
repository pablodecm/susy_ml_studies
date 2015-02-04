#ifndef EventRecObjects_H
#define EventRecObjects_H

/**
 * Container class with the physical object information of interest of the
 * event
 */
class EventRecObjects {
  public:
    float lepton_1_pt;
    float lepton_1_Eta;
    float lepton_1_Phi;
    float lepton_2_Pt;
    float lepton_2_Eta;
    float lepton_2_Phi;
    float pfmet_Et;
    float pfmet_Phi;
    float pfmet_Sig;
    float jet_1_Et;
    float jet_1_Eta;
    float jet_1_Phi;
    float jet_1_CSV;
    float jet_2_Et;
    float jet_2_Eta;
    float jet_2_Pt;
    float jet_2_CSV;


    EventRecObjects() {}

    ClassDef(EventRecObjects, 1);

};

#endif
