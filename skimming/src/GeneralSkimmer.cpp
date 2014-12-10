#define GeneralSkimmer_cxx
// The class definition in GeneralSkimmer.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.

// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// Root > T->Process("GeneralSkimmer.C")
// Root > T->Process("GeneralSkimmer.C","some options")
// Root > T->Process("GeneralSkimmer.C+")
//

#include "GeneralSkimmer.h"
#include <TH2.h>
#include <TStyle.h>


void GeneralSkimmer::Begin(TTree * /*tree*/)
{
    // The Begin() function is called at the start of the query.
    // When running with PROOF Begin() is only called on the client.
    // The tree argument is deprecated (on PROOF 0 is passed).

    TString option = GetOption();

}

void GeneralSkimmer::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

}

Bool_t GeneralSkimmer::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // It can be passed to either GeneralSkimmer::GetEntry() or TBranch::GetEntry()
   // to read either all or the required parts of the data. When processing
   // keyed objects with PROOF, the object is already loaded and is available
   // via the fObject pointer.
   //
   // This function should contain the "body" of the analysis. It can contain
   // simple or elaborate selection criteria, run algorithms on the data
   // of the event and typically fill histograms.
   //
   // The processing can be stopped by calling Abort().
   //
   // Use fStatus to set the return value of TTree::Process().
   //
   // The return value is currently not used.

    int eSize = fChain->GetTree()->GetEntry(entry);
    if(entry % 50000 == 0){
    std::cout << "Events processed: " << entry << std::endl;
    std::cout << " - event size is: " << eSize << " bytes" << std::endl;
    }

    // Electron Selection
    int nElec = T_Elec_Pt->size();
    std::vector<TLorentzVector> vElec; // vector for valid electrons
    std::vector<int> cElec; // electron charge vector
    for (int i=0; i < nElec; i++) {

        // geometry dependent cuts
        bool gCuts = false;
        if( T_Elec_isEB->at(i) && // i.e. barrel electron
             ((fabs(T_Elec_deltaEtaIn->at(i)) < 0.004) &&
              (fabs(T_Elec_deltaPhiIn->at(i)) < 0.06 ) &&
              (T_Elec_sigmaIetaIeta->at(i) < 0.01 ) &&
              (T_Elec_HtoE->at(i) < 0.12 ))) gCuts = true;
        else if( T_Elec_isEE->at(i) && // i.e. endcap electron
             ((T_Elec_deltaEtaIn->at(i) < 0.007) &&
              (T_Elec_deltaPhiIn->at(i) < 0.03 ) &&
              (T_Elec_sigmaIetaIeta->at(i) < 0.03 ) &&
              (T_Elec_HtoE->at(i) < 0.10 ))) gCuts = true;

        // relative electron isolation
        double elecRelIso = (T_Elec_chargedHadronIso->at(i) + std::max(0.0,T_Elec_neutralHadronIso->at(i) + T_Elec_photonIso->at(i) - (T_Event_RhoIso * GetEffectiveArea( T_Elec_SC_Eta->at(i))) ) )/ T_Elec_Pt->at(i);

        if(gCuts &&
           T_Elec_simpleEleId80->at(i) &&
           T_Elec_passConversionVeto->at(i) &&
           T_Elec_isPF->at(i) &&
           (T_Elec_isEB->at(i) || T_Elec_isEE->at(i)) &&
           fabs( (1-T_Elec_eSuperClusterOverP->at(i))/T_Elec_ecalEnergy->at(i)) < 0.05   &&
           (fabs(T_Elec_SC_Eta->at(i) ) < 1.4442 || fabs(T_Elec_SC_Eta->at(i)) > 1.566 ) &&
           fabs( T_Elec_Pt->at(i) - T_Elec_PFElecPt->at(i) ) < 10. &&
           fabs( T_Elec_IPwrtPV->at(i) ) <  0.02 &&
           fabs( T_Elec_dzwrtPV->at(i) ) <  0.1  &&
           elecRelIso <= 0.15 &&
           T_Elec_nHits->at(i) >= 1 &&
           T_Elec_Pt->at(i) >= 20. &&
           fabs(T_Elec_SC_Eta->at(i)) <= 2.5
           ) {
            // add good electron to vector
            TLorentzVector gElec( T_Elec_Px->at(i), T_Elec_Py->at(i),
                                  T_Elec_Pz->at(i), T_Elec_Energy->at(i));
            vElec.push_back(gElec);
            cElec.push_back(T_Elec_Charge->at(i));
        } // end if good electron
    } // end electron loop


    // Muon Selection
    int nMuon = T_Muon_Pt->size();
    std::vector <TLorentzVector> vMuon; // vector for the good muons
    std::vector <int> cMuon; // muon charge vector
    for (int i = 0; i < nMuon; i++ ) {

        // relative muon isolation
        double muonRelIso = muonRelIso = (T_Muon_chargedHadronIsoR03->at(i) + std::max(0.0, T_Muon_neutralHadronIsoR03->at(i) + T_Muon_photonIsoR03->at(i) - 0.5*T_Muon_sumPUPtR03->at(i))) / T_Muon_Pt->at(i);

        if(T_Muon_IsGlobalMuon->at(i) &&
           T_Muon_IsGMPTMuons->at(i) &&
           T_Muon_isPFMuon->at(i) &&
           T_Muon_Chi2InTrk->at(i) < 10 &&  // not in FRs code
           T_Muon_NValidHitsInTrk->at(i) > 0 &&  // not in FRs code
           T_Muon_NumOfMatchedStations->at(i) > 1 &&
           T_Muon_IPwrtAveBSInTrack->at(i) < 0.2 &&
           fabs(T_Muon_vz->at(i) - T_Vertex_z->at(0)) < 0.5 &&
           T_Muon_NLayers->at(i) > 5 &&
           T_Muon_NValidPixelHitsInTrk->at(i) > 0 &&
           muonRelIso <=  0.15  &&
           T_Muon_Pt->at(i) >= 20. &&
           fabs(T_Muon_Eta->at(i)) <=  2.4   &&
           fabs( T_Muon_Pt->at(i) - T_Muon_PFMuonPt->at(i) ) < 5.
           ) {
            TLorentzVector gMuon(T_Muon_Px->at(i), T_Muon_Py->at(i),
                                 T_Muon_Pz->at(i), T_Muon_Energy->at(i));
            vMuon.push_back(gMuon);
            cMuon.push_back(T_Muon_Charge->at(i));
        } // end if good muon
    } // end muon loop

    return kTRUE;
}

void GeneralSkimmer::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void GeneralSkimmer::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

}

//------------------------------------------------------------------------------
// Effective areas for rho correction to electron isolation
//------------------------------------------------------------------------------
Double_t GetEffectiveArea(float eta)
{
  Double_t Aeff=0.;

  if( fabs(eta) < 1.0 )         Aeff = 0.13; // +/- 0.001
  else if( fabs(eta)<1.479 )    Aeff = 0.14; // +/- 0.002
  else if( fabs(eta)<2.0 )      Aeff = 0.07; // +/- 0.001
  else if( fabs(eta)<2.2 )      Aeff = 0.09; // +/- 0.001
  else if( fabs(eta)<2.3 )      Aeff = 0.11; // +/- 0.002
  else if( fabs(eta)<2.4 )      Aeff = 0.11; // +/- 0.003
  if( fabs(eta) > 2.4 )         Aeff = 0.14; // +/- 0.004

  return Aeff;
}
