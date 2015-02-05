
#include "EventRecObjects.h"

void EventRecObjects::SetLeadingLepton( const TLorentzVector& lepton_1) {
  // properly set members
  lepton_1_Pt = lepton_1.Pt();
  lepton_1_Eta = lepton_1.Eta();
  lepton_1_Phi = lepton_1.Phi();
}

void EventRecObjects::SetTrailingLepton( const TLorentzVector& lepton_2) {
  // properly set members
  lepton_2_Pt = lepton_2.Pt();
  lepton_2_Eta = lepton_2.Eta();
  lepton_2_Phi = lepton_2.Phi();
}

void EventRecObjects::SetLeadingJet( const TLorentzVector& jet_1, float jet_1_CSV ) {
  // properly set members
  jet_1_Et = jet_1.Et();
  jet_1_Eta = jet_1.Eta();
  jet_1_Phi = jet_1.Phi();
  jet_1_CSV = jet_1_CSV;
}

void EventRecObjects::SetTrailingJet( const TLorentzVector& jet_2, float jet_2_CSV ) {
  // properly set members
  jet_2_Et = jet_2.Et();
  jet_2_Eta = jet_2.Eta();
  jet_2_Phi = jet_2.Phi();
  jet_2_CSV = jet_2_CSV;
}


