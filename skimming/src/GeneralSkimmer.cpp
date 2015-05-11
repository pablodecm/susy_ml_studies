
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
#include <TFile.h>
#include <limits>





void GeneralSkimmer::Begin(TTree * /*tree*/)
{
    // The Begin() function is called at the start of the query.
    // When running with PROOF Begin() is only called on the client.
    // The tree argument is deprecated (on PROOF 0 is passed).

    std::string option = GetOption();

    std::size_t i_ofile = option.find("ofile="); 
    if (i_ofile != std::string::npos) {
      std::size_t length = (option.find(";", i_ofile) -  option.find("=", i_ofile) - 1);
      o_filename = option.substr(option.find("=", i_ofile)+1 , length );
    } else {
      o_filename = "output.root";
    }

    std::cout << "Output filename: " << o_filename << std::endl;

}

void GeneralSkimmer::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();


    std::cout << "Creating TTree: " << std::endl;
   _skimTree = new TTree("skimTree", "Skimmed TTree for output");
    std::cout << "tree created "  << std::endl;
   
    _ev_data = new EventData();
   _ev_topo = new EventTopology();
   _ev_reco = new EventRecObjects();
   _ev_geno = new EventGenObjects();
   _ev_high = new EventHighLevel();
   _skimTree->Branch("eventData",&_ev_data, 8000, 1);
   _skimTree->Branch("eventTopology",&_ev_topo, 8000, 1);
   _skimTree->Branch("eventRecObjects",&_ev_reco, 8000, 1);
   _skimTree->Branch("eventGenObjects",&_ev_geno, 8000, 1);
   _skimTree->Branch("eventHighLevel",&_ev_high, 8000, 1);
   fOutput->Add(_skimTree);

   h_gen_ch_all = new TH2I("h_gen_ch_all","Events per channel (generation level - no acceptance cuts)",
           3, 0.0, 3.0,3,0.0,3.0);
   h_gen_ch_all->GetXaxis()->SetBinLabel(1,"e");
   h_gen_ch_all->GetXaxis()->SetBinLabel(2,"#mu");
   h_gen_ch_all->GetXaxis()->SetBinLabel(3,"#tau");
   h_gen_ch_all->GetYaxis()->SetBinLabel(1,"e");
   h_gen_ch_all->GetYaxis()->SetBinLabel(2,"#mu");
   h_gen_ch_all->GetYaxis()->SetBinLabel(3,"#tau");
   fOutput->Add(h_gen_ch_all);

   h_gen_ch_acc = new TH2I("h_gen_ch_acc","Events per channel (generation level - with acceptance cuts)",
           3, 0.0, 3.0,3,0.0,3.0);
   h_gen_ch_acc->GetXaxis()->SetBinLabel(1,"e");
   h_gen_ch_acc->GetXaxis()->SetBinLabel(2,"#mu");
   h_gen_ch_acc->GetXaxis()->SetBinLabel(3,"#tau");
   h_gen_ch_acc->GetYaxis()->SetBinLabel(1,"e");
   h_gen_ch_acc->GetYaxis()->SetBinLabel(2,"#mu");
   h_gen_ch_acc->GetYaxis()->SetBinLabel(3,"#tau");
   fOutput->Add(h_gen_ch_acc);

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

    long eSize = fChain->GetTree()->GetEntry(entry);
    if(entry % 50000 == 0){
    std::cout << "Events processed: " << entry << std::endl;
    std::cout << " - event size is: " << eSize << " bytes" << std::endl;
    }

    // store event metadata
    _ev_data->run_number = T_Event_RunNumber;
    _ev_data->lumi_block = T_Event_LuminosityBlock;
    _ev_data->event_number = T_Event_EventNumber;
    _ev_data->process_ID = T_Event_processID;

    // store SUSY generation level  masses
    if (T_Gen_StopMass->size() > 0) {
      _ev_geno->gen_stop_mass = T_Gen_StopMass->at(0);
    }
    if (T_Gen_Chi0Mass->size() > 0) {
      _ev_geno->gen_chi0_mass = T_Gen_Chi0Mass->at(0);
    }
    if (T_Gen_CharginoMass->size() > 0) {
      _ev_geno->gen_chargino_mass = T_Gen_CharginoMass->at(0);
    }

    // Electron Selection
    int nElec = T_Elec_Pt->size();
    std::vector<TLorentzVector> vElec; // vector for valid electrons
    std::vector<int> cElec; // electron charge vector
    for (int i=0; i < nElec; i++) {

        if(isTightElec(i)) {
            // add good electron to vector
            TLorentzVector gElec( T_Elec_Px->at(i), T_Elec_Py->at(i),
                                  T_Elec_Pz->at(i), T_Elec_Energy->at(i));
            vElec.push_back(gElec);
            cElec.push_back(T_Elec_Charge->at(i));
        } // end if tight elec
    } // end electron loop


    // Muon Selection
    int nMuon = T_Muon_Pt->size();
    std::vector <TLorentzVector> vMuon; // vector for the good muons
    std::vector <int> cMuon; // muon charge vector
    for (int i = 0; i < nMuon; i++ ) {

        if(isTightMuon(i)) {
            TLorentzVector gMuon(T_Muon_Px->at(i), T_Muon_Py->at(i),
                                 T_Muon_Pz->at(i), T_Muon_Energy->at(i));
            vMuon.push_back(gMuon);
            cMuon.push_back(T_Muon_Charge->at(i));
        } // end if tight muon
    } // end muon loop

    // get number of loose leptons
    int nVetoLepton = 0;
    for (int i = 0; i < nElec; i++ ) {
        if(isVetoElec(i)) {
            nVetoLepton++;
        } // end if veto elec
    } // end muon loop
   for (int i = 0; i < nMuon; i++ ) {
        if(isVetoMuon(i)) {
            nVetoLepton++;
        } // end if veto muon
    } // end muon loop



    // get number of events at each channel at gen level (no acceptance cuts)
    int nSt3Elec = T_Gen_Elec_Energy->size();
    int nSt3Muon = T_Gen_Muon_Energy->size();
    int nSt3Tau = T_Gen_Tau_Energy->size();
    if(nSt3Elec == 2 && nSt3Muon == 0 && nSt3Tau == 0 ) {
        h_gen_ch_all->Fill("e","e",1);
    } else if (nSt3Elec == 0 && nSt3Muon == 2 && nSt3Tau == 0 ) {
        h_gen_ch_all->Fill("#mu","#mu",1);
    } else if (nSt3Elec == 0 && nSt3Muon == 0 && nSt3Tau == 2 ) {
        h_gen_ch_all->Fill("#tau","#tau",1);
    }


    // obtain channel and selected leptons  (veto any event with additional leptons)
    _ev_topo->n_elec = vElec.size();
    _ev_topo->n_muon = vMuon.size();
    _ev_topo->n_lept = _ev_topo->n_elec + _ev_topo->n_muon;
    bool isOppSign = false;
    _ev_topo->channel = -1;
    std::vector <TLorentzVector> vLept; // to save two selected leptons
    std::vector <TLorentzVector> aLept(vElec.size() + vMuon.size()); // to save all selected leptons
    aLept.insert( aLept.end(), vElec.begin(), vElec.end());
    aLept.insert( aLept.end(), vMuon.begin(), vMuon.end());
    // vElec and vMuon are ordered - get two highest pt leptons
    if (_ev_topo->n_elec > 0) {
      if (_ev_topo->n_muon > 0) { // at least 1 mu and 1 e
        if (vElec[0].Pt() > vMuon[0].Pt() ) { // e is leading lepton
          if(_ev_topo->n_elec > 1) { // at least 1 mu and 2 e
            if (vElec[1].Pt() > vMuon[0].Pt() ) { // e is trailing lepton
              _ev_topo->channel = 0; // ee channel
              vLept.push_back(vElec[0]);
              vLept.push_back(vElec[1]);
              if (cElec[0]*cElec[1] < 0) isOppSign =  true;
            } else { // mu is trailing lepton
              _ev_topo->channel = 2; // emu+mue channel
              vLept.push_back(vElec[0]);
              vLept.push_back(vMuon[0]);
              if (cElec[0]*cMuon[0] < 0) isOppSign =  true;
            }
          } else { // 1 e and 1 mu 
              _ev_topo->channel = 2; // emu+mue channel
              vLept.push_back(vElec[0]);
              vLept.push_back(vMuon[0]);
              if (cElec[0]*cMuon[0] < 0) isOppSign =  true;
          }
        }
      } else if ( _ev_topo->n_elec > 1) { // no muon and 2 or more elecs
        _ev_topo->channel = 0; // ee channel
        vLept.push_back(vElec[0]);
        vLept.push_back(vElec[1]);
        if (cElec[0]*cElec[1] < 0) isOppSign =  true;
      } else { // no muon and less than 2 elecs
        return kFALSE;
      } 
    }

    if (_ev_topo->n_muon > 0) {
      if (_ev_topo->n_elec > 0) { // at least 1 mu and 1 e
        if (vElec[0].Pt() < vMuon[0].Pt() ) { // mu is leading lepton
          if(_ev_topo->n_muon > 1) { // at least 2 mu and 1 e
            if (vElec[0].Pt() < vMuon[1].Pt() ) { // mu is trailing lepton
              _ev_topo->channel = 1; // mumu channel
              vLept.push_back(vMuon[0]);
              vLept.push_back(vMuon[1]);
              if (cMuon[0]*cMuon[1] < 0) isOppSign =  true;
            } else { // e is trailing lepton
              _ev_topo->channel = 2; // emu+mue channel
              vLept.push_back(vMuon[0]);
              vLept.push_back(vElec[0]);
              if (cElec[0]*cMuon[0] < 0) isOppSign =  true;
            }
          } else { // 1 e and 1 mu 
              _ev_topo->channel = 2; // emu+mue channel
              vLept.push_back(vMuon[0]);
              vLept.push_back(vElec[0]);
              if (cElec[0]*cMuon[0] < 0) isOppSign =  true;
          }
 
        }
      } else if ( _ev_topo->n_muon > 1) { // no elec and 2 or more muons
        _ev_topo->channel = 1; // mumu channel
        vLept.push_back(vMuon[0]);
        vLept.push_back(vMuon[1]);
        if (cMuon[0]*cMuon[1] < 0) isOppSign =  true;       
      } else { // no elecs and less than 2 muons
        return kFALSE;
      } 
    }

   
    // Jet Selection
    int nJet= T_JetAKCHS_Et->size();
    std::vector <TLorentzVector> vJet;
    std::vector <float> vJet_CSV;
    _ev_topo->n_jet = 0;
    _ev_topo->n_b_jet = 0;
    _ev_high->jets_ht = 0;
    _ev_high->max_CSV = std::numeric_limits<float>::lowest();
    for ( int i = 0 ; i < nJet; i++) {
        TLorentzVector jet( T_JetAKCHS_Px->at(i), T_JetAKCHS_Py->at(i),
                             T_JetAKCHS_Pz->at(i), T_JetAKCHS_Energy->at(i));
        // check if jets are cleaned of selected leptons
        bool isClean = true;
        for (auto lept : aLept ) {
          if (jet.DeltaR(lept) < 0.4) isClean = false;
        } // end cleaning check
        if (T_JetAKCHS_Et->at(i) > 30. &&
            fabs(T_JetAKCHS_Eta->at(i)) < 2.4 &&
            isClean && passJetID(i) ) {
            // add energy to ht
            _ev_high->jets_ht += T_JetAKCHS_Et->at(i);
            _ev_topo->n_jet++;
            // check if bjet ( CSV medium working point)
            bool isBJet = T_JetAKCHS_Tag_CombInclusiveSVtxV2->at(i) > 0.423;
            if (T_JetAKCHS_Tag_CombInclusiveSVtxV2->at(i) >_ev_high->max_CSV) 
               _ev_high->max_CSV = T_JetAKCHS_Tag_CombInclusiveSVtxV2->at(i);
            if (isBJet) _ev_topo->n_b_jet++;
            // keep only two higher pt jets
            if (vJet.size() < 3) {
              vJet.push_back(jet);
              vJet_CSV.push_back(T_JetAKCHS_Tag_CombInclusiveSVtxV2->at(i));
            } // end only saving two jets 
        } // end good jet loop
    } // end jet loop

    
    if (vLept.size() > 1) {
    // assign leptons
    _ev_reco->SetLeadingLepton(vLept[0]);
    _ev_reco->SetTrailingLepton(vLept[1]);
     // dilepton invariant mass
     _ev_high->dilept_inv_mass = (vLept[0]+vLept[1]).M();
     // dilepton MT2
     _ev_high->dilepton_MT2 = getMT2(vLept[0], vLept[1],
                                    _ev_reco->pfmet_Et, _ev_reco->pfmet_Phi);
    }

    // assign tranverse energy variables
    _ev_reco->pfmet_Et = T_METPF_ET; 
    _ev_reco->pfmet_Phi = T_METPF_Phi;
    // MET sig not avaliable for PHYS14
    _ev_reco->pfmet_Sig = -1; 


    if(vJet.size() > 1 ) {
    // assign jets
    _ev_reco->SetLeadingJet(vJet[0], vJet_CSV[0]);
    _ev_reco->SetTrailingJet(vJet[1], vJet_CSV[1]);
    
    }
         

    // basic selection for ee and mumu channels 
    if ( (_ev_topo->channel == 0 || _ev_topo->channel == 1 ) && 
         isOppSign  && (_ev_high->dilept_inv_mass > 20.)
         && (nVetoLepton < 3)
       ) {
      _skimTree->Fill();
    }

    // basic selection for emu and mue channel
    if ( (_ev_topo->channel == 2) &&
         isOppSign && (_ev_high->dilept_inv_mass > 20.)
         && (nVetoLepton < 3)
       ) {
      _skimTree->Fill();
    }
         
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

    TTree* skimTree = dynamic_cast<TTree *>(fOutput->FindObject(Form("skimTree")));
    std::cout << "Total number of selected entries: " << skimTree->GetEntries() << std::endl;
    std::cout << "  - ee      channel: " << skimTree->Draw("channel", "channel == 0", "goff") << std::endl;
    std::cout << "  - mumu    channel: " << skimTree->Draw("channel", "channel == 1", "goff") << std::endl;
    std::cout << "  - emu+mue channel: " << skimTree->Draw("channel", "channel == 2", "goff") << std::endl;

    TH2I* h_gen_ch_all = dynamic_cast<TH2I *>(fOutput->FindObject(Form("h_gen_ch_all")));
    std::cout << "Total number of signal entries: " << h_gen_ch_all->Integral() << std::endl;
    std::cout << "  - ee      channel: " << h_gen_ch_all->GetBinContent(1,1) << std::endl;
    std::cout << "  - mumu    channel: " << h_gen_ch_all->GetBinContent(2,2) << std::endl;
    std::cout << "  - tau-tau channel: " << h_gen_ch_all->GetBinContent(3,3) << std::endl;

    TFile *o_file = new TFile(o_filename.c_str(), "RECREATE");
    skimTree->Write();

}

bool GeneralSkimmer::isTightMuon( unsigned iMuon, float minPt) 
{
  TLorentzVector lep;
  lep.SetPxPyPzE(T_Muon_Px->at(iMuon), T_Muon_Py->at(iMuon),
                 T_Muon_Pz->at(iMuon), T_Muon_Energy->at(iMuon));
  if (lep.Pt() < minPt)            return false;
  if (TMath::Abs(lep.Eta()) > 2.4) return false;
  
  // POG Tight Muons definition              
  if (!T_Muon_IsGlobalMuon->at(iMuon)                            ) return false;
  if (!T_Muon_IsPFMuon->at(iMuon)                                ) return false;
  if (T_Muon_NormChi2GTrk->at(iMuon)                       >= 10.) return false;
  if (T_Muon_NValidHitsGTrk->at(iMuon)                     <= 0  ) return false;
  if (T_Muon_NumOfMatchedStations->at(iMuon)               <= 1  ) return false;                     
  //if (TMath::Abs(T_Muon_IPwrtAveBSInTrack->at(iMuon))      >= 0.2) return false; 
  //if (TMath::Abs(T_Muon_vz->at(iMuon) - T_Vertex_z->at(0)) >= 0.5) return false;
  if (TMath::Abs(T_Muon_dxyInTrack->at(iMuon))             >= 0.2) return false; 
  if (TMath::Abs(T_Muon_dzInTrack ->at(iMuon))             >= 0.5) return false;
  if (T_Muon_NLayers->at(iMuon)                            <= 5  ) return false;
  if (T_Muon_NValidPixelHitsInTrk->at(iMuon)               <= 0  ) return false;

  float relIso = muonIsolation(iMuon);
  
  if (relIso > 0.12) return false;
  
  return true;
}

bool GeneralSkimmer::isVetoMuon( unsigned iMuon, float minPt) 
{
  TLorentzVector lep;
  lep.SetPxPyPzE(T_Muon_Px->at(iMuon), T_Muon_Py->at(iMuon),
                 T_Muon_Pz->at(iMuon), T_Muon_Energy->at(iMuon));
  if (lep.Pt() < minPt)            return false;
  if (TMath::Abs(lep.Eta()) > 2.4) return false;
  float relIso = muonIsolation(iMuon);
  if (relIso > 0.20)               return false; 
  if (T_Muon_IsPFMuon->at(iMuon) == 0) return false;
  if (T_Muon_IsGlobalMuon->at(iMuon) == 0 && T_Muon_IsTrackerMuonArbitrated->at(iMuon) == 0) return false; 
  return true;
}

float GeneralSkimmer::muonIsolation(unsigned iMuon){
  if (iMuon < 0) return 9999.;
  if (iMuon >= (int)T_Muon_chargedHadronIsoR04->size()) return 9999.;

  TLorentzVector lep;
  lep.SetPxPyPzE(T_Muon_Px->at(iMuon), T_Muon_Py->at(iMuon),
                 T_Muon_Pz->at(iMuon), T_Muon_Energy->at(iMuon));
  
  return (T_Muon_chargedHadronIsoR04->at(iMuon) + std::max(0.0 , T_Muon_neutralHadronIsoR04->at(iMuon) + T_Muon_photonIsoR04->at(iMuon)- 0.5*T_Muon_sumPUPtR04->at(iMuon)))/lep.Pt();
}

bool GeneralSkimmer::isTightElec(unsigned int iElec, float minPt){
  
  TLorentzVector lep( T_Elec_Px->at(iElec), T_Elec_Py->at(iElec),
                      T_Elec_Pz->at(iElec), T_Elec_Energy->at(iElec));
  
  float pt  = lep.Pt();
  float sceta = TMath::Abs(T_Elec_SC_Eta->at(iElec));
  if (sceta > 1.4442 && sceta < 1.566) return false;
  if (lep.Pt() < minPt)                return false;
  if (TMath::Abs(lep.Eta()) > 2.5)     return false;

  float relIso =  elecIsolation(iElec);

  // Tight ID requirements
  bool passTightID = false;
  if(TMath::Abs(T_Elec_SC_Eta->at(iElec)) < 1.479) {  //Barrel electron
    if(T_Elec_sigmaIetaIetaFull5by5->at(iElec)              < 0.010181 && //
       TMath::Abs(T_Elec_deltaEtaIn->at(iElec))             < 0.006574 && //
       TMath::Abs(T_Elec_deltaPhiIn->at(iElec))             < 0.022868 && //       
       T_Elec_HtoE->at(iElec)                               < 0.037553 && //
       TMath::Abs( 1.0/T_Elec_ecalEnergy->at(iElec) -
                   T_Elec_eSuperClusterOverP->at(iElec) /
                   T_Elec_ecalEnergy->at(iElec))            < 0.131191 &&
       TMath::Abs(T_Elec_IPwrtPV->at(iElec))                < 0.009924 && //
       TMath::Abs(T_Elec_vz->at(iElec) - T_Vertex_z->at(0)) < 0.015310 &&
       T_Elec_nLost->at(iElec)                              <= 1       && //
       T_Elec_passConversionVeto->at(iElec)                 > 0        &&
       relIso                                               < 0.074355)
      passTightID = true;
  }
  else {
    if(T_Elec_sigmaIetaIetaFull5by5->at(iElec)              < 0.028766 && //
       TMath::Abs(T_Elec_deltaEtaIn->at(iElec))             < 0.005681 && //
       TMath::Abs(T_Elec_deltaPhiIn->at(iElec))             < 0.032046 && //       
       T_Elec_HtoE->at(iElec)                               < 0.081902 && //
       TMath::Abs( 1.0/T_Elec_ecalEnergy->at(iElec) -
                   T_Elec_eSuperClusterOverP->at(iElec) /
                   T_Elec_ecalEnergy->at(iElec))            < 0.106055 &&
       TMath::Abs(T_Elec_IPwrtPV->at(iElec))                < 0.027261 && //
       TMath::Abs(T_Elec_vz->at(iElec) - T_Vertex_z->at(0)) < 0.147154 &&
       T_Elec_nLost->at(iElec)                              <= 1       && //
       T_Elec_passConversionVeto->at(iElec)                 > 0        &&
       relIso                                               < 0.090185)
       passTightID = true;
  }  
  
  // Medium ID requirements
  bool passMediumID = false;
  if(TMath::Abs(T_Elec_SC_Eta->at(iElec)) < 1.479) {  //Barrel electron
    if(T_Elec_sigmaIetaIetaFull5by5->at(iElec)              < 0.010399 && 
       TMath::Abs(T_Elec_deltaEtaIn->at(iElec))             < 0.007641 && 
       TMath::Abs(T_Elec_deltaPhiIn->at(iElec))             < 0.032643 &&    
       T_Elec_HtoE->at(iElec)                               < 0.060662 && 
       relIso                                               < 0.097213 &&
       TMath::Abs( 1.0/T_Elec_ecalEnergy->at(iElec) - 
                  T_Elec_eSuperClusterOverP->at(iElec) / 
                  T_Elec_ecalEnergy->at(iElec))             < 0.153897 &&
       TMath::Abs(T_Elec_IPwrtPV->at(iElec))                < 0.011811 && 
       //TMath::Abs(T_Elec_vz->at(iElec) - T_Vertex_z->at(0)) < 0.070775 &&
       TMath::Abs(T_Elec_dzwrtPV->at(iElec))                < 0.070775 &&
       T_Elec_nLost->at(iElec)                              <= 1       && 
       T_Elec_passConversionVeto->at(iElec)                 > 0 )
      passMediumID = true;
  }
  else {
    if(T_Elec_sigmaIetaIetaFull5by5->at(iElec)              < 0.029524 && 
       TMath::Abs(T_Elec_deltaEtaIn->at(iElec))             < 0.009285 && 
       TMath::Abs(T_Elec_deltaPhiIn->at(iElec))             < 0.042447 &&    
       T_Elec_HtoE->at(iElec)                               < 0.104263 && 
       relIso                                               < 0.116708 &&
       TMath::Abs( 1.0/T_Elec_ecalEnergy->at(iElec) - 
       T_Elec_eSuperClusterOverP->at(iElec) /
       T_Elec_ecalEnergy->at(iElec))                        < 0.137468 &&
       TMath::Abs(T_Elec_IPwrtPV->at(iElec))                < 0.051682 && 
       //TMath::Abs(T_Elec_vz->at(iElec) - T_Vertex_z->at(0)) < 0.180720 &&
       TMath::Abs(T_Elec_dzwrtPV->at(iElec))                < 0.180720 &&
       T_Elec_nLost->at(iElec)                              <= 1       && 
       T_Elec_passConversionVeto->at(iElec)                 > 0 )
       passMediumID = true;
  }
  
  if (!passMediumID) return false;
   
  return true;
}

bool GeneralSkimmer::isVetoElec(unsigned int iElec, float minPt){
  
  TLorentzVector lep( T_Elec_Px->at(iElec), T_Elec_Py->at(iElec),
                      T_Elec_Pz->at(iElec), T_Elec_Energy->at(iElec));
  
  float pt  = lep.Pt();
  float sceta = TMath::Abs(T_Elec_SC_Eta->at(iElec));
  if (sceta > 1.4442 && sceta < 1.566) return false;
  if (lep.Pt() < minPt)                return false;
  if (TMath::Abs(lep.Eta()) > 2.5)     return false;

  float relIso =  elecIsolation(iElec);

  bool passVetoID = false;
  // veto ID requirements
  if(TMath::Abs(T_Elec_SC_Eta->at(iElec)) < 1.479) {  //Barrel electron
    if(T_Elec_sigmaIetaIetaFull5by5->at(iElec)              < 0.011100 &&
       TMath::Abs(T_Elec_deltaEtaIn->at(iElec))             < 0.016315 && 
       TMath::Abs(T_Elec_deltaPhiIn->at(iElec))             < 0.252044 &&       
       T_Elec_HtoE->at(iElec)                               < 0.345843 && 
       relIso                                               < 0.164369 &&
       ( 1.0/T_Elec_ecalEnergy->at(iElec) -
         T_Elec_eSuperClusterOverP->at(iElec) /
         T_Elec_ecalEnergy->at(iElec))                      < 0.248070 &&
       TMath::Abs(T_Elec_IPwrtPV->at(iElec))                < 0.060279 && 
       TMath::Abs(T_Elec_dzwrtPV->at(iElec))                < 0.800538 &&
       T_Elec_nLost->at(iElec)                              <= 2       && 
       T_Elec_passConversionVeto->at(iElec)                 > 0 )
      passVetoID = true;
  }
  else {
    if(T_Elec_sigmaIetaIetaFull5by5->at(iElec)              < 0.033987 && 
       TMath::Abs(T_Elec_deltaEtaIn->at(iElec))             < 0.010671 && 
       TMath::Abs(T_Elec_deltaPhiIn->at(iElec))             < 0.245263 &&    
       T_Elec_HtoE->at(iElec)                               < 0.134691 && 
       relIso                                               < 0.212604 &&
       ( 1.0/T_Elec_ecalEnergy->at(iElec) -
         T_Elec_eSuperClusterOverP->at(iElec) /
         T_Elec_ecalEnergy->at(iElec))                      < 0.157160 &&
       TMath::Abs(T_Elec_IPwrtPV->at(iElec))                < 0.273097 && 
       TMath::Abs(T_Elec_dzwrtPV->at(iElec))                < 0.885860 &&
       T_Elec_nLost->at(iElec)                              <= 3       && 
       T_Elec_passConversionVeto->at(iElec)                 > 0 )
       passVetoID = true;
  }  
  
  if (!passVetoID) return false;
  
  return true;
}


float GeneralSkimmer::elecIsolation(unsigned iElec) {

  if (iElec < 0) return 9999.;
  if (iElec >= (int)T_Elec_chargedHadronIso->size()) return 9999.;
  
  TLorentzVector lep( T_Elec_Px->at(iElec), T_Elec_Py->at(iElec),
                      T_Elec_Pz->at(iElec), T_Elec_Energy->at(iElec));

  float pt     = lep.Pt();
  float relIso = (T_Elec_sumChargedHadronPt->at(iElec) +
                  std::max(0.0 , T_Elec_sumNeutralHadronEt->at(iElec) +
                                 T_Elec_sumPhotonEt->at(iElec) -
                                 0.5*T_Elec_sumPUPt->at(iElec)))/pt;
  return relIso;

}

float GeneralSkimmer::elecEffArea(float eta){  // for a 0.3 CONE
  float abseta = TMath::Abs(eta);
  
  if      (abseta < 1.0)                      return 0.13; // +/- 0.001
  else if (abseta >= 1.0   && abseta < 1.479) return 0.14; // +/- 0.002
  else if (abseta >= 1.479 && abseta < 2.0)   return 0.07; // +/- 0.001
  else if (abseta >= 2.0   && abseta < 2.2)   return 0.09; // +/- 0.001
  else if (abseta >= 2.2   && abseta < 2.3)   return 0.11; // +/- 0.002
  else if (abseta >= 2.3   && abseta < 2.4)   return 0.11; // +/- 0.003
  else if (abseta >= 2.4)                     return 0.14; // +/- 0.004
  
  std::cout << "[ERROR] getEACorrection(): No correction factor found!!" << std::endl;
  return -999999.;
}

bool GeneralSkimmer::passJetID(unsigned iJet) {

  if ( !(T_JetAKCHS_nDaughters->at(iJet)        > 1   ) ) return false;
  if ( !(T_JetAKCHS_NeutHadEnergyFrac->at(iJet) < 0.99) ) return false;
  if ( !(T_JetAKCHS_NeutEmEnergyFrac ->at(iJet) < 0.99) ) return false;
  if (TMath::Abs(T_JetAKCHS_Eta->at(iJet)) < 2.5){
  if ( !(T_JetAKCHS_CharEmEnergyFrac->at(iJet)  < 0.99) ) return false;
  if ( !(T_JetAKCHS_CharHadEnergyFrac->at(iJet) > 0.00) ) return false;
  if ( !(T_JetAKCHS_ChargedMultiplicity->at(iJet) > 0 ) ) return false;
  return true;
  }
}



//-----------------------------------------------------------------------------
// MT2 calculation by mt2bisect.h
//-----------------------------------------------------------------------------
float getMT2(TLorentzVector lep0, TLorentzVector lep1, float met_Et, float met_Phi) {
    double pa[3];
    double pb[3];
    double pmiss[3];

    TLorentzVector pmet;
    pmet.SetPtEtaPhiM(met_Et, 0., met_Phi, 0.);
    pmiss[0] = 0.; // not required
    pmiss[1] = pmet.Px();
    pmiss[2] = pmet.Py();

    pa[0] = 0.;
    pa[1] = lep0.Px();
    pa[2] = lep0.Py();

    pb[0] = 0.;
    pb[1] = lep1.Px();
    pb[2] = lep1.Py();

    mt2bisect* MT2bisect = new mt2bisect();
    MT2bisect->set_verbose(0);
    MT2bisect->set_momenta(pa, pb, pmiss);
    MT2bisect->set_mn(0.); // test mass
    double mt2 = MT2bisect->get_mt2();
    delete MT2bisect;
    return mt2;
}
