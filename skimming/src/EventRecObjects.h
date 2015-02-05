#ifndef EventRecObjects_H
#define EventRecObjects_H

#include <TLorentzVector.h>

/**
 * Container class with the physical object information of interest of the
 * event
 */
class EventRecObjects {
  public:
    float lepton_1_Pt;
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
    float jet_2_Phi;
    float jet_2_CSV;


    EventRecObjects() {}

    void SetLeadingLepton(const TLorentzVector& lepton_1);
    void SetTrailingLepton(const TLorentzVector& lepton_2);
    void SetLeadingJet(const TLorentzVector& jet_1, float CSV_value);
    void SetTrailingJet(const TLorentzVector& jet_2, float CSV_value);

    ClassDef(EventRecObjects, 1);

};

#endif
