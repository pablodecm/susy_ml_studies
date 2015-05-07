
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

        // geometry dependent cuts
        bool gCuts = false;
        if( T_Elec_isEB->at(i) && // i.e. barrel electron
             ((fabs(T_Elec_deltaEtaIn->at(i)) < 0.004) &&
              (fabs(T_Elec_deltaPhiIn->at(i)) < 0.06 ) &&
              (T_Elec_sigmaIetaIeta->at(i) < 0.01 ) &&
              (T_Elec_HtoE->at(i) < 0.12 ))) gCuts = true;
        else if( T_Elec_isEE->at(i) && // i.e. endcap electron
             ((fabs(T_Elec_deltaEtaIn->at(i)) < 0.007) &&
              (fabs(T_Elec_deltaPhiIn->at(i)) < 0.03 ) &&
              (T_Elec_sigmaIetaIeta->at(i) < 0.03 ) &&
              (T_Elec_HtoE->at(i) < 0.10 ))) gCuts = true;

        // relative electron isolation
        double elecRelIso = (T_Elec_chargedHadronIso->at(i) + std::max(0.0,
                    T_Elec_neutralHadronIso->at(i) + T_Elec_photonIso->at(i)
                    - (T_Event_Rho * GetEffectiveArea( T_Elec_SC_Eta->at(i)))))/ T_Elec_Pt->at(i);

        if(gCuts &&
           // T_Elec_MVAoutput->at(i) && // need to check working point
           T_Elec_passConversionVeto->at(i) &&
           // T_Elec_isPF->at(i) &&
           (T_Elec_isEB->at(i) || T_Elec_isEE->at(i)) &&
           fabs( (1-T_Elec_eSuperClusterOverP->at(i))/T_Elec_ecalEnergy->at(i)) < 0.05   &&
           (fabs(T_Elec_SC_Eta->at(i) ) < 1.4442 || fabs(T_Elec_SC_Eta->at(i)) > 1.566 ) &&
           // fabs( T_Elec_Pt->at(i) - T_Elec_PFElecPt->at(i) ) < 10. &&
           fabs( T_Elec_IPwrtPV->at(i) ) <  0.02 &&
           fabs( T_Elec_dzwrtPV->at(i) ) <  0.1  &&
           elecRelIso <= 0.15 &&
           // T_Elec_nHits->at(i) > 0 && // there seems to be a problem with this variable
           T_Elec_Pt->at(i) >= 10. &&
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
        double muonRelIso = (T_Muon_chargedHadronIsoR03->at(i) + std::max(0.0,
                    T_Muon_neutralHadronIsoR03->at(i) + T_Muon_photonIsoR03->at(i)
                    - 0.5*T_Muon_sumPUPtR03->at(i))) / T_Muon_Pt->at(i);

        if(T_Muon_IsGlobalMuon->at(i) &&
           T_Muon_IsGMPTMuons->at(i) &&
           T_Muon_IsPFMuon->at(i) && 
           T_Muon_NormChi2GTrk->at(i) < 10 &&  // not in FRs code ( normalized is important!)
           T_Muon_NValidHitsInTrk->at(i) > 0 &&  // not in FRs code
           T_Muon_NumOfMatchedStations->at(i) > 1 &&
           T_Muon_IPwrtAveBSInTrack->at(i) < 0.2 &&
           fabs(T_Muon_vz->at(i) - T_Vertex_z->at(0)) < 0.5 &&
           T_Muon_NLayers->at(i) > 5 &&
           T_Muon_NValidPixelHitsInTrk->at(i) > 0 &&
           muonRelIso <=  0.15  &&
           T_Muon_Pt->at(i) >= 10. &&
           fabs(T_Muon_Eta->at(i)) <=  2.4 //  &&
           // fabs( T_Muon_Pt->at(i) - T_Muon_PFMuonPt->at(i) ) < 5. // not avaliable PHYS14
           ) {
            TLorentzVector gMuon(T_Muon_Px->at(i), T_Muon_Py->at(i),
                                 T_Muon_Pz->at(i), T_Muon_Energy->at(i));
            vMuon.push_back(gMuon);
            cMuon.push_back(T_Muon_Charge->at(i));
        } // end if good muon
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

    // get number of events at each channel at gen level (with acceptance cuts)
    // int nAccElec = T_Gen_Elec_Energy->size();
    // int nAccMuon = T_Gen_Muon_Energy->size();
    // int nAccTau = T_Gen_Tau_Energy->size();

    // obtain channel and selected leptons  (veto any event with additional leptons)
    _ev_topo->n_elec = vElec.size();
    _ev_topo->n_muon = vMuon.size();
    _ev_topo->n_lept = _ev_topo->n_elec + _ev_topo->n_muon;
    bool isOppSign = false;
    _ev_topo->channel = -1;
    std::vector <TLorentzVector> vLept;
    if(_ev_topo->n_elec == 2 && _ev_topo->n_muon == 0) {
        _ev_topo->channel = 0; // ee channel
        vLept.push_back(vElec[0]);
        vLept.push_back(vElec[1]);
        if (cElec[0]*cElec[1] < 0) isOppSign =  true;
    } else if (_ev_topo->n_elec == 0 && _ev_topo->n_muon == 2 ) {
        _ev_topo->channel = 1; // mumu channel
        vLept.push_back(vMuon[0]);
        vLept.push_back(vMuon[1]);
        if (cMuon[0]*cMuon[1] < 0) isOppSign =  true;
    } else if (_ev_topo->n_elec == 1 && _ev_topo->n_muon == 1 ) {
        _ev_topo->channel = 2; // emu+mue channel
        if (vElec[0].Pt() > vMuon[0].Pt() ) {
            vLept.push_back(vElec[0]);
            vLept.push_back(vMuon[0]);
        }
        else if (vMuon[0].Pt() > vElec[0].Pt()) {
            vLept.push_back(vMuon[0]);
            vLept.push_back(vElec[0]);  
        }
        if (cElec[0]*cMuon[0] < 0) isOppSign =  true;
    } else {
      return kFALSE; // exit if not two leptons are present
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
        for (auto lept : vLept ) {
          if (jet.DeltaR(lept) < 0.5) isClean = false;
        } // end cleaning check
        if (T_JetAKCHS_Et->at(i) > 30 &&
            fabs(T_JetAKCHS_Eta->at(i)) < 2.4 &&
            isClean ) {
            // add energy to ht
            _ev_high->jets_ht += T_JetAKCHS_Et->at(i);
            _ev_topo->n_jet++;
            // check if bjet ( CSV medium working point)
            bool isBJet = T_JetAKCHS_Tag_pfCombinedSVtx->at(i) > 0.679;
            if (T_JetAKCHS_Tag_pfCombinedSVtx->at(i) >_ev_high->max_CSV) 
               _ev_high->max_CSV = T_JetAKCHS_Tag_pfCombinedSVtx->at(i);
            if (isBJet) _ev_topo->n_b_jet++;
            // keep only two higher pt jets
            if (vJet.size() < 3) {
              vJet.push_back(jet);
              vJet_CSV.push_back(T_JetAKCHS_Tag_pfCombinedSVtx->at(i));
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
         isOppSign 
       ) {
      _skimTree->Fill();
    }

    // basic selection for emu and mue channel
    if ( (_ev_topo->channel == 2) &&
         isOppSign 
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
