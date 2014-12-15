//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Dec  5 12:23:41 2014 by ROOT version 5.34/18
// from TTree Tree/
// found on file: /gpfs/csic_projects/tier3data/MC_Summer12_53X/SUSY/Tree_TTbar_Madgraph_0.root
//////////////////////////////////////////////////////////

#ifndef GeneralSkimmer_h
#define GeneralSkimmer_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TLorentzVector.h>

// Header file for the classes stored in the TTree if any.
#include <vector>
#include <iostream>

struct EventInfo{
    int runNumber;
    int lumiBlock;
    int eventNumber;
    int processID;
};

struct EventData{
    int channel;
};


// Fixed size dimensions of array or collections stored in the TTree if any.

class GeneralSkimmer : public TSelector {

protected:
   EventData      _eventData;
   EventData      _genData;
   TTree          *skimTree;
   TTree          *genTree;
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain

   // Declaration of leaf types
   Float_t         T_Event_Rho;
   Float_t         T_Event_RhoIso;
   Bool_t          T_EventF_EcalDeadCell;
   Bool_t          T_EventF_logErrorTooManyClusters;
   Bool_t          T_EventF_trackingFailure;
   Bool_t          T_EventF_hcalLaser;
   Bool_t          T_EventF_ecalLaserCorr;
   Bool_t          T_EventF_toomanystripclus;
   Bool_t          T_EventF_manystripclus;
   Bool_t          T_EventF_eeBadSc;
   Int_t           T_Event_RunNumber;
   Int_t           T_Event_EventNumber;
   Int_t           T_Event_LuminosityBlock;
   Int_t           T_Event_processID;
   Int_t           T_Event_nPU;
   Float_t         T_Event_nTruePU;
   Int_t           T_Event_nPUm;
   Int_t           T_Event_nPUp;
   Float_t         T_Event_AveNTruePU;
   Bool_t          T_HLT_Mu8_v16;
   Bool_t          T_HLT_Mu17_v3;
   Bool_t          T_HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v12;
   Bool_t          T_HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v13;
   Bool_t          T_HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v14;
   Bool_t          T_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v3;
   Bool_t          T_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v4;
   Bool_t          T_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v5;
   std::vector<float>   *T_Gen_StopMass;
   std::vector<float>   *T_Gen_Chi0Mass;
   std::vector<float>   *T_Gen_CharginoMass;
   std::vector<double>  *T_Gen_polWeights;
   std::vector<int>     *T_Gen_Muon_PID;
   std::vector<float>   *T_Gen_Muon_Px;
   std::vector<float>   *T_Gen_Muon_Py;
   std::vector<float>   *T_Gen_Muon_Pz;
   std::vector<float>   *T_Gen_Muon_Energy;
   std::vector<int>     *T_Gen_Muon_MPID;
   std::vector<float>   *T_Gen_Muon_MPx;
   std::vector<float>   *T_Gen_Muon_MPy;
   std::vector<float>   *T_Gen_Muon_MPz;
   std::vector<float>   *T_Gen_Muon_MEnergy;
   std::vector<int>     *T_Gen_Muon_MSt;
   std::vector<int>     *T_Gen_Elec_PID;
   std::vector<float>   *T_Gen_Elec_Px;
   std::vector<float>   *T_Gen_Elec_Py;
   std::vector<float>   *T_Gen_Elec_Pz;
   std::vector<float>   *T_Gen_Elec_Energy;
   std::vector<int>     *T_Gen_Elec_MPID;
   std::vector<float>   *T_Gen_Elec_MPx;
   std::vector<float>   *T_Gen_Elec_MPy;
   std::vector<float>   *T_Gen_Elec_MPz;
   std::vector<float>   *T_Gen_Elec_MEnergy;
   std::vector<int>     *T_Gen_Elec_MSt;
   std::vector<int>     *T_Gen_b_PID;
   std::vector<float>   *T_Gen_b_Px;
   std::vector<float>   *T_Gen_b_Py;
   std::vector<float>   *T_Gen_b_Pz;
   std::vector<float>   *T_Gen_b_Energy;
   std::vector<int>     *T_Gen_b_MPID;
   std::vector<float>   *T_Gen_b_MPx;
   std::vector<float>   *T_Gen_b_MPy;
   std::vector<float>   *T_Gen_b_MPz;
   std::vector<float>   *T_Gen_b_MEnergy;
   std::vector<int>     *T_Gen_b_MSt;
   std::vector<int>     *T_Gen_MuonSt3_pdgId;
   std::vector<int>     *T_Gen_MuonSt3_firstMother;
   std::vector<int>     *T_Gen_MuonSt3_i;
   std::vector<float>   *T_Gen_MuonSt3_energy;
   std::vector<float>   *T_Gen_MuonSt3_pt;
   std::vector<float>   *T_Gen_MuonSt3_eta;
   std::vector<float>   *T_Gen_MuonSt3_phi;
   std::vector<int>     *T_Gen_ElecSt3_pdgId;
   std::vector<int>     *T_Gen_ElecSt3_firstMother;
   std::vector<int>     *T_Gen_ElecSt3_i;
   std::vector<float>   *T_Gen_ElecSt3_energy;
   std::vector<float>   *T_Gen_ElecSt3_pt;
   std::vector<float>   *T_Gen_ElecSt3_eta;
   std::vector<float>   *T_Gen_ElecSt3_phi;
   std::vector<int>     *T_Gen_TauSt3_pdgId;
   std::vector<int>     *T_Gen_TauSt3_firstMother;
   std::vector<int>     *T_Gen_TauSt3_i;
   std::vector<float>   *T_Gen_TauSt3_energy;
   std::vector<float>   *T_Gen_TauSt3_pt;
   std::vector<float>   *T_Gen_TauSt3_eta;
   std::vector<float>   *T_Gen_TauSt3_phi;
   std::vector<int>     *T_Gen_StopSt3_pdgId;
   std::vector<int>     *T_Gen_StopSt3_firstMother;
   std::vector<int>     *T_Gen_StopSt3_i;
   std::vector<float>   *T_Gen_StopSt3_energy;
   std::vector<float>   *T_Gen_StopSt3_pt;
   std::vector<float>   *T_Gen_StopSt3_eta;
   std::vector<float>   *T_Gen_StopSt3_phi;
   std::vector<int>     *T_Gen_Chi0St3_pdgId;
   std::vector<int>     *T_Gen_Chi0St3_firstMother;
   std::vector<int>     *T_Gen_Chi0St3_i;
   std::vector<float>   *T_Gen_Chi0St3_energy;
   std::vector<float>   *T_Gen_Chi0St3_pt;
   std::vector<float>   *T_Gen_Chi0St3_eta;
   std::vector<float>   *T_Gen_Chi0St3_phi;
   std::vector<int>     *T_Gen_ChiPMSt3_pdgId;
   std::vector<int>     *T_Gen_ChiPMSt3_firstMother;
   std::vector<int>     *T_Gen_ChiPMSt3_i;
   std::vector<float>   *T_Gen_ChiPMSt3_energy;
   std::vector<float>   *T_Gen_ChiPMSt3_pt;
   std::vector<float>   *T_Gen_ChiPMSt3_eta;
   std::vector<float>   *T_Gen_ChiPMSt3_phi;
   std::vector<int>     *T_Gen_bSt3_pdgId;
   std::vector<int>     *T_Gen_bSt3_firstMother;
   std::vector<int>     *T_Gen_bSt3_i;
   std::vector<float>   *T_Gen_bSt3_energy;
   std::vector<float>   *T_Gen_bSt3_pt;
   std::vector<float>   *T_Gen_bSt3_eta;
   std::vector<float>   *T_Gen_bSt3_phi;
   std::vector<int>     *T_Gen_tSt3_pdgId;
   std::vector<int>     *T_Gen_tSt3_firstMother;
   std::vector<int>     *T_Gen_tSt3_i;
   std::vector<float>   *T_Gen_tSt3_energy;
   std::vector<float>   *T_Gen_tSt3_pt;
   std::vector<float>   *T_Gen_tSt3_eta;
   std::vector<float>   *T_Gen_tSt3_phi;
   std::vector<int>     *T_Gen_WSt3_pdgId;
   std::vector<int>     *T_Gen_WSt3_firstMother;
   std::vector<int>     *T_Gen_WSt3_i;
   std::vector<float>   *T_Gen_WSt3_energy;
   std::vector<float>   *T_Gen_WSt3_pt;
   std::vector<float>   *T_Gen_WSt3_eta;
   std::vector<float>   *T_Gen_WSt3_phi;
   std::vector<bool>    *T_Gen_TauSt3_IsLepDec;
   std::vector<int>     *T_Gen_TauSt3_LepDec_PID;
   std::vector<float>   *T_Gen_TauSt3_LepDec_Px;
   std::vector<float>   *T_Gen_TauSt3_LepDec_Py;
   std::vector<float>   *T_Gen_TauSt3_LepDec_Pz;
   std::vector<float>   *T_Gen_TauSt3_LepDec_Energy;
   std::vector<float>   *T_Muon_Eta;
   std::vector<bool>    *T_Muon_IsGlobalMuon;
   std::vector<bool>    *T_Muon_IsAllTrackerMuons;
   std::vector<bool>    *T_Muon_IsTrackerMuonArbitrated;
   std::vector<bool>    *T_Muon_IsGMPTMuons;
   std::vector<float>   *T_Muon_SegmentCompatibility;
   std::vector<float>   *T_Muon_trkKink;
   std::vector<float>   *T_Muon_Px;
   std::vector<float>   *T_Muon_Py;
   std::vector<float>   *T_Muon_Pz;
   std::vector<float>   *T_Muon_Pt;
   std::vector<float>   *T_Muon_deltaPt;
   std::vector<float>   *T_Muon_Energy;
   std::vector<int>     *T_Muon_Charge;
   std::vector<float>   *T_Muon_NormChi2GTrk;
   std::vector<int>     *T_Muon_NValidHitsInTrk;
   std::vector<int>     *T_Muon_NValidPixelHitsInTrk;
   std::vector<int>     *T_Muon_InnerTrackFound;
   std::vector<int>     *T_Muon_NValidHitsSATrk;
   std::vector<int>     *T_Muon_NValidHitsGTrk;
   std::vector<float>   *T_Muon_Chi2InTrk;
   std::vector<float>   *T_Muon_dofInTrk;
   std::vector<float>   *T_Muon_IPAbsGTrack;
   std::vector<float>   *T_Muon_IPAbsInTrack;
   std::vector<float>   *T_Muon_IPwrtAveBSInTrack;
   std::vector<float>   *T_Muon_chargedHadronIsoR04;
   std::vector<float>   *T_Muon_neutralHadronIsoR04;
   std::vector<float>   *T_Muon_photonIsoR04;
   std::vector<float>   *T_Muon_chargedParticleIsoR03;
   std::vector<float>   *T_Muon_chargedHadronIsoR03;
   std::vector<float>   *T_Muon_neutralHadronIsoR03;
   std::vector<float>   *T_Muon_photonIsoR03;
   std::vector<float>   *T_Muon_sumPUPtR04;
   std::vector<float>   *T_Muon_sumPUPtR03;
   std::vector<float>   *T_Muon_vz;
   std::vector<float>   *T_Muon_vy;
   std::vector<float>   *T_Muon_vx;
   std::vector<int>     *T_Muon_NumOfMatchedStations;
   std::vector<float>   *T_Muon_PFMuonPt;
   std::vector<float>   *T_Muon_PFMuonPx;
   std::vector<float>   *T_Muon_PFMuonPy;
   std::vector<float>   *T_Muon_PFMuonPz;
   std::vector<float>   *T_Muon_PFMuonE;
   std::vector<bool>    *T_Muon_isPFMuon;
   std::vector<int>     *T_Muon_NLayers;
   std::vector<float>   *T_Vertex_z;
   std::vector<float>   *T_Vertex_y;
   std::vector<float>   *T_Vertex_x;
   std::vector<float>   *T_Vertex_Chi2Prob;
   std::vector<float>   *T_Vertex_rho;
   std::vector<float>   *T_Vertex_ndof;
   std::vector<bool>    *T_Vertex_isFake;
   std::vector<int>     *T_Vertex_tracksSize;
   std::vector<float>   *T_Elec_Eta;
   std::vector<float>   *T_Elec_IPwrtAveBS;
   std::vector<float>   *T_Elec_IPwrtPV;
   std::vector<float>   *T_Elec_dzwrtPV;
   std::vector<float>   *T_Elec_Px;
   std::vector<float>   *T_Elec_Py;
   std::vector<float>   *T_Elec_Pz;
   std::vector<float>   *T_Elec_Pt;
   std::vector<float>   *T_Elec_Energy;
   std::vector<int>     *T_Elec_Charge;
   std::vector<int>     *T_Elec_nBrems;
   std::vector<float>   *T_Elec_fBrem;
   std::vector<float>   *T_Elec_eSuperClusterOverP;
   std::vector<float>   *T_Elec_ecalEnergy;
   std::vector<float>   *T_Elec_dr03TkSumPt;
   std::vector<float>   *T_Elec_dr03EcalSumEt;
   std::vector<float>   *T_Elec_dr03HcalSumEt;
   std::vector<bool>    *T_Elec_isEB;
   std::vector<bool>    *T_Elec_isEE;
   std::vector<float>   *T_Elec_MVA;
   std::vector<float>   *T_Elec_simpleEleId80;
   std::vector<float>   *T_Elec_chargedHadronIso;
   std::vector<float>   *T_Elec_neutralHadronIso;
   std::vector<float>   *T_Elec_photonIso;
   std::vector<float>   *T_Elec_puChargedHadronIso;
   std::vector<bool>    *T_Elec_passConversionVeto;
   std::vector<float>   *T_Elec_sigmaIetaIeta;
   std::vector<float>   *T_Elec_deltaPhiIn;
   std::vector<float>   *T_Elec_deltaEtaIn;
   std::vector<bool>    *T_Elec_isEcalDriven;
   std::vector<float>   *T_Elec_HtoE;
   std::vector<float>   *T_Elec_vz;
   std::vector<float>   *T_Elec_vy;
   std::vector<float>   *T_Elec_vx;
   std::vector<int>     *T_Elec_nLost;
   std::vector<int>     *T_Elec_nHits;
   std::vector<float>   *T_Elec_SC_Et;
   std::vector<float>   *T_Elec_SC_Eta;
   std::vector<float>   *T_Elec_PFElecPt;
   std::vector<float>   *T_Elec_PFElecPx;
   std::vector<float>   *T_Elec_PFElecPy;
   std::vector<float>   *T_Elec_PFElecPz;
   std::vector<float>   *T_Elec_PFElecE;
   std::vector<bool>    *T_Elec_isPF;
   std::vector<float>   *T_JetAKCHS_Px;
   std::vector<float>   *T_JetAKCHS_Py;
   std::vector<float>   *T_JetAKCHS_Pz;
   std::vector<float>   *T_JetAKCHS_Et;
   std::vector<float>   *T_JetAKCHS_Eta;
   std::vector<float>   *T_JetAKCHS_Energy;
   std::vector<float>   *T_JetAKCHS_Tag_HighEffTC;
   std::vector<float>   *T_JetAKCHS_Tag_CombSVtx;
   std::vector<float>   *T_JetAKCHS_Tag_CombSVtxMVA;
   std::vector<float>   *T_JetAKCHS_Tag_TauJet;
   std::vector<float>   *T_JetAKCHS_Tag_ImpParMVA;
   std::vector<float>   *T_JetAKCHS_Tag_JetBProb;
   std::vector<float>   *T_JetAKCHS_Tag_JetProb;
   std::vector<float>   *T_JetAKCHS_Tag_HighEffSimpSVtx;
   std::vector<float>   *T_JetAKCHS_Tag_HighPurSimpSVtx;
   std::vector<float>   *T_JetAKCHS_Tag_HighPurTC;
   std::vector<float>   *T_JetAKCHS_Parton_Px;
   std::vector<float>   *T_JetAKCHS_Parton_Py;
   std::vector<float>   *T_JetAKCHS_Parton_Pz;
   std::vector<float>   *T_JetAKCHS_Parton_Energy;
   std::vector<int>     *T_JetAKCHS_Parton_Flavour;
   std::vector<float>   *T_JetAKCHS_CharHadEnergyFrac;
   std::vector<float>   *T_JetAKCHS_NeutHadEnergyFrac;
   std::vector<float>   *T_JetAKCHS_CharEmEnergyFrac;
   std::vector<float>   *T_JetAKCHS_NeutEmEnergyFrac;
   std::vector<float>   *T_JetAKCHS_CharHadEnergy;
   std::vector<float>   *T_JetAKCHS_NeutHadEnergy;
   std::vector<float>   *T_JetAKCHS_CharEmEnergy;
   std::vector<float>   *T_JetAKCHS_NeutEmEnergy;
   std::vector<int>     *T_JetAKCHS_MuonMultiplicity;
   std::vector<int>     *T_JetAKCHS_NeutralMultiplicity;
   std::vector<int>     *T_JetAKCHS_ChargedMultiplicity;
   std::vector<bool>    *T_JetAKCHS_IDLoose;
   std::vector<int>     *T_JetAKCHS_nDaughters;
   std::vector<float>   *T_JetAKCHS_GenJet_InvisibleE;
   std::vector<float>   *T_JetAKCHS_GenJet_Px;
   std::vector<float>   *T_JetAKCHS_GenJet_Py;
   std::vector<float>   *T_JetAKCHS_GenJet_Pz;
   std::vector<float>   *T_JetAKCHS_GenJet_Et;
   std::vector<float>   *T_JetAKCHS_GenJet_Eta;
   std::vector<float>   *T_JetAKCHS_GenJet_Energy;
   std::vector<bool>    *T_JetAKCHS_IsGenJet;
   Float_t         T_METPF_ET;
   Float_t         T_METPF_Phi;
   Float_t         T_METPF_Sig;
   Float_t         T_METPFTypeI_ET;
   Float_t         T_METPFTypeI_Phi;
   Float_t         T_METPFTypeI_Sig;
   Float_t         T_METgen_ET;
   Float_t         T_METgen_Phi;
   Bool_t          T_passTriggerDoubleMu;
   Bool_t          T_passTriggerDoubleEl;
   Bool_t          T_passTriggerElMu;

   // List of branches
   TBranch        *b_T_Event_Rho;   //!
   TBranch        *b_T_Event_RhoIso;   //!
   TBranch        *b_T_EventF_EcalDeadCell;   //!
   TBranch        *b_T_EventF_logErrorTooManyClusters;   //!
   TBranch        *b_T_EventF_trackingFailure;   //!
   TBranch        *b_T_EventF_hcalLaser;   //!
   TBranch        *b_T_EventF_ecalLaserCorr;   //!
   TBranch        *b_T_EventF_toomanystripclus;   //!
   TBranch        *b_T_EventF_manystripclus;   //!
   TBranch        *b_T_EventF_eeBadSc;   //!
   TBranch        *b_T_Event_RunNumber;   //!
   TBranch        *b_T_Event_EventNumber;   //!
   TBranch        *b_T_Event_LuminosityBlock;   //!
   TBranch        *b_T_Event_processID;   //!
   TBranch        *b_T_Event_nPU;   //!
   TBranch        *b_T_Event_nTruePU;   //!
   TBranch        *b_T_Event_nPUm;   //!
   TBranch        *b_T_Event_nPUp;   //!
   TBranch        *b_T_Event_AveNTruePU;   //!
   TBranch        *b_T_HLT_Mu8_v16;   //!
   TBranch        *b_T_HLT_Mu17_v3;   //!
   TBranch        *b_T_HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v12;   //!
   TBranch        *b_T_HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v13;   //!
   TBranch        *b_T_HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v14;   //!
   TBranch        *b_T_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v3;   //!
   TBranch        *b_T_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v4;   //!
   TBranch        *b_T_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v5;   //!
   TBranch        *b_T_Gen_StopMass;   //!
   TBranch        *b_T_Gen_Chi0Mass;   //!
   TBranch        *b_T_Gen_CharginoMass;   //!
   TBranch        *b_T_Gen_polWeights;   //!
   TBranch        *b_T_Gen_Muon_PID;   //!
   TBranch        *b_T_Gen_Muon_Px;   //!
   TBranch        *b_T_Gen_Muon_Py;   //!
   TBranch        *b_T_Gen_Muon_Pz;   //!
   TBranch        *b_T_Gen_Muon_Energy;   //!
   TBranch        *b_T_Gen_Muon_MPID;   //!
   TBranch        *b_T_Gen_Muon_MPx;   //!
   TBranch        *b_T_Gen_Muon_MPy;   //!
   TBranch        *b_T_Gen_Muon_MPz;   //!
   TBranch        *b_T_Gen_Muon_MEnergy;   //!
   TBranch        *b_T_Gen_Muon_MSt;   //!
   TBranch        *b_T_Gen_Elec_PID;   //!
   TBranch        *b_T_Gen_Elec_Px;   //!
   TBranch        *b_T_Gen_Elec_Py;   //!
   TBranch        *b_T_Gen_Elec_Pz;   //!
   TBranch        *b_T_Gen_Elec_Energy;   //!
   TBranch        *b_T_Gen_Elec_MPID;   //!
   TBranch        *b_T_Gen_Elec_MPx;   //!
   TBranch        *b_T_Gen_Elec_MPy;   //!
   TBranch        *b_T_Gen_Elec_MPz;   //!
   TBranch        *b_T_Gen_Elec_MEnergy;   //!
   TBranch        *b_T_Gen_Elec_MSt;   //!
   TBranch        *b_T_Gen_b_PID;   //!
   TBranch        *b_T_Gen_b_Px;   //!
   TBranch        *b_T_Gen_b_Py;   //!
   TBranch        *b_T_Gen_b_Pz;   //!
   TBranch        *b_T_Gen_b_Energy;   //!
   TBranch        *b_T_Gen_b_MPID;   //!
   TBranch        *b_T_Gen_b_MPx;   //!
   TBranch        *b_T_Gen_b_MPy;   //!
   TBranch        *b_T_Gen_b_MPz;   //!
   TBranch        *b_T_Gen_b_MEnergy;   //!
   TBranch        *b_T_Gen_b_MSt;   //!
   TBranch        *b_T_Gen_MuonSt3_pdgId;   //!
   TBranch        *b_T_Gen_MuonSt3_firstMother;   //!
   TBranch        *b_T_Gen_MuonSt3_i;   //!
   TBranch        *b_T_Gen_MuonSt3_energy;   //!
   TBranch        *b_T_Gen_MuonSt3_pt;   //!
   TBranch        *b_T_Gen_MuonSt3_eta;   //!
   TBranch        *b_T_Gen_MuonSt3_phi;   //!
   TBranch        *b_T_Gen_ElecSt3_pdgId;   //!
   TBranch        *b_T_Gen_ElecSt3_firstMother;   //!
   TBranch        *b_T_Gen_ElecSt3_i;   //!
   TBranch        *b_T_Gen_ElecSt3_energy;   //!
   TBranch        *b_T_Gen_ElecSt3_pt;   //!
   TBranch        *b_T_Gen_ElecSt3_eta;   //!
   TBranch        *b_T_Gen_ElecSt3_phi;   //!
   TBranch        *b_T_Gen_TauSt3_pdgId;   //!
   TBranch        *b_T_Gen_TauSt3_firstMother;   //!
   TBranch        *b_T_Gen_TauSt3_i;   //!
   TBranch        *b_T_Gen_TauSt3_energy;   //!
   TBranch        *b_T_Gen_TauSt3_pt;   //!
   TBranch        *b_T_Gen_TauSt3_eta;   //!
   TBranch        *b_T_Gen_TauSt3_phi;   //!
   TBranch        *b_T_Gen_StopSt3_pdgId;   //!
   TBranch        *b_T_Gen_StopSt3_firstMother;   //!
   TBranch        *b_T_Gen_StopSt3_i;   //!
   TBranch        *b_T_Gen_StopSt3_energy;   //!
   TBranch        *b_T_Gen_StopSt3_pt;   //!
   TBranch        *b_T_Gen_StopSt3_eta;   //!
   TBranch        *b_T_Gen_StopSt3_phi;   //!
   TBranch        *b_T_Gen_Chi0St3_pdgId;   //!
   TBranch        *b_T_Gen_Chi0St3_firstMother;   //!
   TBranch        *b_T_Gen_Chi0St3_i;   //!
   TBranch        *b_T_Gen_Chi0St3_energy;   //!
   TBranch        *b_T_Gen_Chi0St3_pt;   //!
   TBranch        *b_T_Gen_Chi0St3_eta;   //!
   TBranch        *b_T_Gen_Chi0St3_phi;   //!
   TBranch        *b_T_Gen_ChiPMSt3_pdgId;   //!
   TBranch        *b_T_Gen_ChiPMSt3_firstMother;   //!
   TBranch        *b_T_Gen_ChiPMSt3_i;   //!
   TBranch        *b_T_Gen_ChiPMSt3_energy;   //!
   TBranch        *b_T_Gen_ChiPMSt3_pt;   //!
   TBranch        *b_T_Gen_ChiPMSt3_eta;   //!
   TBranch        *b_T_Gen_ChiPMSt3_phi;   //!
   TBranch        *b_T_Gen_bSt3_pdgId;   //!
   TBranch        *b_T_Gen_bSt3_firstMother;   //!
   TBranch        *b_T_Gen_bSt3_i;   //!
   TBranch        *b_T_Gen_bSt3_energy;   //!
   TBranch        *b_T_Gen_bSt3_pt;   //!
   TBranch        *b_T_Gen_bSt3_eta;   //!
   TBranch        *b_T_Gen_bSt3_phi;   //!
   TBranch        *b_T_Gen_tSt3_pdgId;   //!
   TBranch        *b_T_Gen_tSt3_firstMother;   //!
   TBranch        *b_T_Gen_tSt3_i;   //!
   TBranch        *b_T_Gen_tSt3_energy;   //!
   TBranch        *b_T_Gen_tSt3_pt;   //!
   TBranch        *b_T_Gen_tSt3_eta;   //!
   TBranch        *b_T_Gen_tSt3_phi;   //!
   TBranch        *b_T_Gen_WSt3_pdgId;   //!
   TBranch        *b_T_Gen_WSt3_firstMother;   //!
   TBranch        *b_T_Gen_WSt3_i;   //!
   TBranch        *b_T_Gen_WSt3_energy;   //!
   TBranch        *b_T_Gen_WSt3_pt;   //!
   TBranch        *b_T_Gen_WSt3_eta;   //!
   TBranch        *b_T_Gen_WSt3_phi;   //!
   TBranch        *b_T_Gen_TauSt3_IsLepDec;   //!
   TBranch        *b_T_Gen_TauSt3_LepDec_PID;   //!
   TBranch        *b_T_Gen_TauSt3_LepDec_Px;   //!
   TBranch        *b_T_Gen_TauSt3_LepDec_Py;   //!
   TBranch        *b_T_Gen_TauSt3_LepDec_Pz;   //!
   TBranch        *b_T_Gen_TauSt3_LepDec_Energy;   //!
   TBranch        *b_T_Muon_Eta;   //!
   TBranch        *b_T_Muon_IsGlobalMuon;   //!
   TBranch        *b_T_Muon_IsAllTrackerMuons;   //!
   TBranch        *b_T_Muon_IsTrackerMuonArbitrated;   //!
   TBranch        *b_T_Muon_IsGMPTMuons;   //!
   TBranch        *b_T_Muon_SegmentCompatibility;   //!
   TBranch        *b_T_Muon_trkKink;   //!
   TBranch        *b_T_Muon_Px;   //!
   TBranch        *b_T_Muon_Py;   //!
   TBranch        *b_T_Muon_Pz;   //!
   TBranch        *b_T_Muon_Pt;   //!
   TBranch        *b_T_Muon_deltaPt;   //!
   TBranch        *b_T_Muon_Energy;   //!
   TBranch        *b_T_Muon_Charge;   //!
   TBranch        *b_T_Muon_NormChi2GTrk;   //!
   TBranch        *b_T_Muon_NValidHitsInTrk;   //!
   TBranch        *b_T_Muon_NValidPixelHitsInTrk;   //!
   TBranch        *b_T_Muon_InnerTrackFound;   //!
   TBranch        *b_T_Muon_NValidHitsSATrk;   //!
   TBranch        *b_T_Muon_NValidHitsGTrk;   //!
   TBranch        *b_T_Muon_Chi2InTrk;   //!
   TBranch        *b_T_Muon_dofInTrk;   //!
   TBranch        *b_T_Muon_IPAbsGTrack;   //!
   TBranch        *b_T_Muon_IPAbsInTrack;   //!
   TBranch        *b_T_Muon_IPwrtAveBSInTrack;   //!
   TBranch        *b_T_Muon_chargedHadronIsoR04;   //!
   TBranch        *b_T_Muon_neutralHadronIsoR04;   //!
   TBranch        *b_T_Muon_photonIsoR04;   //!
   TBranch        *b_T_Muon_chargedParticleIsoR03;   //!
   TBranch        *b_T_Muon_chargedHadronIsoR03;   //!
   TBranch        *b_T_Muon_neutralHadronIsoR03;   //!
   TBranch        *b_T_Muon_photonIsoR03;   //!
   TBranch        *b_T_Muon_sumPUPtR04;   //!
   TBranch        *b_T_Muon_sumPUPtR03;   //!
   TBranch        *b_T_Muon_vz;   //!
   TBranch        *b_T_Muon_vy;   //!
   TBranch        *b_T_Muon_vx;   //!
   TBranch        *b_T_Muon_NumOfMatchedStations;   //!
   TBranch        *b_T_Muon_PFMuonPt;   //!
   TBranch        *b_T_Muon_PFMuonPx;   //!
   TBranch        *b_T_Muon_PFMuonPy;   //!
   TBranch        *b_T_Muon_PFMuonPz;   //!
   TBranch        *b_T_Muon_PFMuonE;   //!
   TBranch        *b_T_Muon_isPFMuon;   //!
   TBranch        *b_T_Muon_NLayers;   //!
   TBranch        *b_T_Vertex_z;   //!
   TBranch        *b_T_Vertex_y;   //!
   TBranch        *b_T_Vertex_x;   //!
   TBranch        *b_T_Vertex_Chi2Prob;   //!
   TBranch        *b_T_Vertex_rho;   //!
   TBranch        *b_T_Vertex_ndof;   //!
   TBranch        *b_T_Vertex_isFake;   //!
   TBranch        *b_T_Vertex_tracksSize;   //!
   TBranch        *b_T_Elec_Eta;   //!
   TBranch        *b_T_Elec_IPwrtAveBS;   //!
   TBranch        *b_T_Elec_IPwrtPV;   //!
   TBranch        *b_T_Elec_dzwrtPV;   //!
   TBranch        *b_T_Elec_Px;   //!
   TBranch        *b_T_Elec_Py;   //!
   TBranch        *b_T_Elec_Pz;   //!
   TBranch        *b_T_Elec_Pt;   //!
   TBranch        *b_T_Elec_Energy;   //!
   TBranch        *b_T_Elec_Charge;   //!
   TBranch        *b_T_Elec_nBrems;   //!
   TBranch        *b_T_Elec_fBrem;   //!
   TBranch        *b_T_Elec_eSuperClusterOverP;   //!
   TBranch        *b_T_Elec_ecalEnergy;   //!
   TBranch        *b_T_Elec_dr03TkSumPt;   //!
   TBranch        *b_T_Elec_dr03EcalSumEt;   //!
   TBranch        *b_T_Elec_dr03HcalSumEt;   //!
   TBranch        *b_T_Elec_isEB;   //!
   TBranch        *b_T_Elec_isEE;   //!
   TBranch        *b_T_Elec_MVA;   //!
   TBranch        *b_T_Elec_simpleEleId80;   //!
   TBranch        *b_T_Elec_chargedHadronIso;   //!
   TBranch        *b_T_Elec_neutralHadronIso;   //!
   TBranch        *b_T_Elec_photonIso;   //!
   TBranch        *b_T_Elec_puChargedHadronIso;   //!
   TBranch        *b_T_Elec_passConversionVeto;   //!
   TBranch        *b_T_Elec_sigmaIetaIeta;   //!
   TBranch        *b_T_Elec_deltaPhiIn;   //!
   TBranch        *b_T_Elec_deltaEtaIn;   //!
   TBranch        *b_T_Elec_isEcalDriven;   //!
   TBranch        *b_T_Elec_HtoE;   //!
   TBranch        *b_T_Elec_vz;   //!
   TBranch        *b_T_Elec_vy;   //!
   TBranch        *b_T_Elec_vx;   //!
   TBranch        *b_T_Elec_nLost;   //!
   TBranch        *b_T_Elec_nHits;   //!
   TBranch        *b_T_Elec_SC_Et;   //!
   TBranch        *b_T_Elec_SC_Eta;   //!
   TBranch        *b_T_Elec_PFElecPt;   //!
   TBranch        *b_T_Elec_PFElecPx;   //!
   TBranch        *b_T_Elec_PFElecPy;   //!
   TBranch        *b_T_Elec_PFElecPz;   //!
   TBranch        *b_T_Elec_PFElecE;   //!
   TBranch        *b_T_Elec_isPF;   //!
   TBranch        *b_T_JetAKCHS_Px;   //!
   TBranch        *b_T_JetAKCHS_Py;   //!
   TBranch        *b_T_JetAKCHS_Pz;   //!
   TBranch        *b_T_JetAKCHS_Et;   //!
   TBranch        *b_T_JetAKCHS_Eta;   //!
   TBranch        *b_T_JetAKCHS_Energy;   //!
   TBranch        *b_T_JetAKCHS_Tag_HighEffTC;   //!
   TBranch        *b_T_JetAKCHS_Tag_CombSVtx;   //!
   TBranch        *b_T_JetAKCHS_Tag_CombSVtxMVA;   //!
   TBranch        *b_T_JetAKCHS_Tag_TauJet;   //!
   TBranch        *b_T_JetAKCHS_Tag_ImpParMVA;   //!
   TBranch        *b_T_JetAKCHS_Tag_JetBProb;   //!
   TBranch        *b_T_JetAKCHS_Tag_JetProb;   //!
   TBranch        *b_T_JetAKCHS_Tag_HighEffSimpSVtx;   //!
   TBranch        *b_T_JetAKCHS_Tag_HighPurSimpSVtx;   //!
   TBranch        *b_T_JetAKCHS_Tag_HighPurTC;   //!
   TBranch        *b_T_JetAKCHS_Parton_Px;   //!
   TBranch        *b_T_JetAKCHS_Parton_Py;   //!
   TBranch        *b_T_JetAKCHS_Parton_Pz;   //!
   TBranch        *b_T_JetAKCHS_Parton_Energy;   //!
   TBranch        *b_T_JetAKCHS_Parton_Flavour;   //!
   TBranch        *b_T_JetAKCHS_CharHadEnergyFrac;   //!
   TBranch        *b_T_JetAKCHS_NeutHadEnergyFrac;   //!
   TBranch        *b_T_JetAKCHS_CharEmEnergyFrac;   //!
   TBranch        *b_T_JetAKCHS_NeutEmEnergyFrac;   //!
   TBranch        *b_T_JetAKCHS_CharHadEnergy;   //!
   TBranch        *b_T_JetAKCHS_NeutHadEnergy;   //!
   TBranch        *b_T_JetAKCHS_CharEmEnergy;   //!
   TBranch        *b_T_JetAKCHS_NeutEmEnergy;   //!
   TBranch        *b_T_JetAKCHS_MuonMultiplicity;   //!
   TBranch        *b_T_JetAKCHS_NeutralMultiplicity;   //!
   TBranch        *b_T_JetAKCHS_ChargedMultiplicity;   //!
   TBranch        *b_T_JetAKCHS_IDLoose;   //!
   TBranch        *b_T_JetAKCHS_nDaughters;   //!
   TBranch        *b_T_JetAKCHS_GenJet_InvisibleE;   //!
   TBranch        *b_T_JetAKCHS_GenJet_Px;   //!
   TBranch        *b_T_JetAKCHS_GenJet_Py;   //!
   TBranch        *b_T_JetAKCHS_GenJet_Pz;   //!
   TBranch        *b_T_JetAKCHS_GenJet_Et;   //!
   TBranch        *b_T_JetAKCHS_GenJet_Eta;   //!
   TBranch        *b_T_JetAKCHS_GenJet_Energy;   //!
   TBranch        *b_T_JetAKCHS_IsGenJet;   //!
   TBranch        *b_T_METPF_ET;   //!
   TBranch        *b_T_METPF_Phi;   //!
   TBranch        *b_T_METPF_Sig;   //!
   TBranch        *b_T_METPFTypeI_ET;   //!
   TBranch        *b_T_METPFTypeI_Phi;   //!
   TBranch        *b_T_METPFTypeI_Sig;   //!
   TBranch        *b_T_METgen_ET;   //!
   TBranch        *b_T_METgen_Phi;   //!
   TBranch        *b_T_passTriggerDoubleMu;   //!
   TBranch        *b_T_passTriggerDoubleEl;   //!
   TBranch        *b_T_passTriggerElMu;   //!

   GeneralSkimmer(TTree * /*tree*/ =0) : fChain(0), skimTree(0) { }
   virtual ~GeneralSkimmer() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

   ClassDef(GeneralSkimmer,0);
};

#endif

#ifdef GeneralSkimmer_cxx
void GeneralSkimmer::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   T_Gen_StopMass = 0;
   T_Gen_Chi0Mass = 0;
   T_Gen_CharginoMass = 0;
   T_Gen_polWeights = 0;
   T_Gen_Muon_PID = 0;
   T_Gen_Muon_Px = 0;
   T_Gen_Muon_Py = 0;
   T_Gen_Muon_Pz = 0;
   T_Gen_Muon_Energy = 0;
   T_Gen_Muon_MPID = 0;
   T_Gen_Muon_MPx = 0;
   T_Gen_Muon_MPy = 0;
   T_Gen_Muon_MPz = 0;
   T_Gen_Muon_MEnergy = 0;
   T_Gen_Muon_MSt = 0;
   T_Gen_Elec_PID = 0;
   T_Gen_Elec_Px = 0;
   T_Gen_Elec_Py = 0;
   T_Gen_Elec_Pz = 0;
   T_Gen_Elec_Energy = 0;
   T_Gen_Elec_MPID = 0;
   T_Gen_Elec_MPx = 0;
   T_Gen_Elec_MPy = 0;
   T_Gen_Elec_MPz = 0;
   T_Gen_Elec_MEnergy = 0;
   T_Gen_Elec_MSt = 0;
   T_Gen_b_PID = 0;
   T_Gen_b_Px = 0;
   T_Gen_b_Py = 0;
   T_Gen_b_Pz = 0;
   T_Gen_b_Energy = 0;
   T_Gen_b_MPID = 0;
   T_Gen_b_MPx = 0;
   T_Gen_b_MPy = 0;
   T_Gen_b_MPz = 0;
   T_Gen_b_MEnergy = 0;
   T_Gen_b_MSt = 0;
   T_Gen_MuonSt3_pdgId = 0;
   T_Gen_MuonSt3_firstMother = 0;
   T_Gen_MuonSt3_i = 0;
   T_Gen_MuonSt3_energy = 0;
   T_Gen_MuonSt3_pt = 0;
   T_Gen_MuonSt3_eta = 0;
   T_Gen_MuonSt3_phi = 0;
   T_Gen_ElecSt3_pdgId = 0;
   T_Gen_ElecSt3_firstMother = 0;
   T_Gen_ElecSt3_i = 0;
   T_Gen_ElecSt3_energy = 0;
   T_Gen_ElecSt3_pt = 0;
   T_Gen_ElecSt3_eta = 0;
   T_Gen_ElecSt3_phi = 0;
   T_Gen_TauSt3_pdgId = 0;
   T_Gen_TauSt3_firstMother = 0;
   T_Gen_TauSt3_i = 0;
   T_Gen_TauSt3_energy = 0;
   T_Gen_TauSt3_pt = 0;
   T_Gen_TauSt3_eta = 0;
   T_Gen_TauSt3_phi = 0;
   T_Gen_StopSt3_pdgId = 0;
   T_Gen_StopSt3_firstMother = 0;
   T_Gen_StopSt3_i = 0;
   T_Gen_StopSt3_energy = 0;
   T_Gen_StopSt3_pt = 0;
   T_Gen_StopSt3_eta = 0;
   T_Gen_StopSt3_phi = 0;
   T_Gen_Chi0St3_pdgId = 0;
   T_Gen_Chi0St3_firstMother = 0;
   T_Gen_Chi0St3_i = 0;
   T_Gen_Chi0St3_energy = 0;
   T_Gen_Chi0St3_pt = 0;
   T_Gen_Chi0St3_eta = 0;
   T_Gen_Chi0St3_phi = 0;
   T_Gen_ChiPMSt3_pdgId = 0;
   T_Gen_ChiPMSt3_firstMother = 0;
   T_Gen_ChiPMSt3_i = 0;
   T_Gen_ChiPMSt3_energy = 0;
   T_Gen_ChiPMSt3_pt = 0;
   T_Gen_ChiPMSt3_eta = 0;
   T_Gen_ChiPMSt3_phi = 0;
   T_Gen_bSt3_pdgId = 0;
   T_Gen_bSt3_firstMother = 0;
   T_Gen_bSt3_i = 0;
   T_Gen_bSt3_energy = 0;
   T_Gen_bSt3_pt = 0;
   T_Gen_bSt3_eta = 0;
   T_Gen_bSt3_phi = 0;
   T_Gen_tSt3_pdgId = 0;
   T_Gen_tSt3_firstMother = 0;
   T_Gen_tSt3_i = 0;
   T_Gen_tSt3_energy = 0;
   T_Gen_tSt3_pt = 0;
   T_Gen_tSt3_eta = 0;
   T_Gen_tSt3_phi = 0;
   T_Gen_WSt3_pdgId = 0;
   T_Gen_WSt3_firstMother = 0;
   T_Gen_WSt3_i = 0;
   T_Gen_WSt3_energy = 0;
   T_Gen_WSt3_pt = 0;
   T_Gen_WSt3_eta = 0;
   T_Gen_WSt3_phi = 0;
   T_Gen_TauSt3_IsLepDec = 0;
   T_Gen_TauSt3_LepDec_PID = 0;
   T_Gen_TauSt3_LepDec_Px = 0;
   T_Gen_TauSt3_LepDec_Py = 0;
   T_Gen_TauSt3_LepDec_Pz = 0;
   T_Gen_TauSt3_LepDec_Energy = 0;
   T_Muon_Eta = 0;
   T_Muon_IsGlobalMuon = 0;
   T_Muon_IsAllTrackerMuons = 0;
   T_Muon_IsTrackerMuonArbitrated = 0;
   T_Muon_IsGMPTMuons = 0;
   T_Muon_SegmentCompatibility = 0;
   T_Muon_trkKink = 0;
   T_Muon_Px = 0;
   T_Muon_Py = 0;
   T_Muon_Pz = 0;
   T_Muon_Pt = 0;
   T_Muon_deltaPt = 0;
   T_Muon_Energy = 0;
   T_Muon_Charge = 0;
   T_Muon_NormChi2GTrk = 0;
   T_Muon_NValidHitsInTrk = 0;
   T_Muon_NValidPixelHitsInTrk = 0;
   T_Muon_InnerTrackFound = 0;
   T_Muon_NValidHitsSATrk = 0;
   T_Muon_NValidHitsGTrk = 0;
   T_Muon_Chi2InTrk = 0;
   T_Muon_dofInTrk = 0;
   T_Muon_IPAbsGTrack = 0;
   T_Muon_IPAbsInTrack = 0;
   T_Muon_IPwrtAveBSInTrack = 0;
   T_Muon_chargedHadronIsoR04 = 0;
   T_Muon_neutralHadronIsoR04 = 0;
   T_Muon_photonIsoR04 = 0;
   T_Muon_chargedParticleIsoR03 = 0;
   T_Muon_chargedHadronIsoR03 = 0;
   T_Muon_neutralHadronIsoR03 = 0;
   T_Muon_photonIsoR03 = 0;
   T_Muon_sumPUPtR04 = 0;
   T_Muon_sumPUPtR03 = 0;
   T_Muon_vz = 0;
   T_Muon_vy = 0;
   T_Muon_vx = 0;
   T_Muon_NumOfMatchedStations = 0;
   T_Muon_PFMuonPt = 0;
   T_Muon_PFMuonPx = 0;
   T_Muon_PFMuonPy = 0;
   T_Muon_PFMuonPz = 0;
   T_Muon_PFMuonE = 0;
   T_Muon_isPFMuon = 0;
   T_Muon_NLayers = 0;
   T_Vertex_z = 0;
   T_Vertex_y = 0;
   T_Vertex_x = 0;
   T_Vertex_Chi2Prob = 0;
   T_Vertex_rho = 0;
   T_Vertex_ndof = 0;
   T_Vertex_isFake = 0;
   T_Vertex_tracksSize = 0;
   T_Elec_Eta = 0;
   T_Elec_IPwrtAveBS = 0;
   T_Elec_IPwrtPV = 0;
   T_Elec_dzwrtPV = 0;
   T_Elec_Px = 0;
   T_Elec_Py = 0;
   T_Elec_Pz = 0;
   T_Elec_Pt = 0;
   T_Elec_Energy = 0;
   T_Elec_Charge = 0;
   T_Elec_nBrems = 0;
   T_Elec_fBrem = 0;
   T_Elec_eSuperClusterOverP = 0;
   T_Elec_ecalEnergy = 0;
   T_Elec_dr03TkSumPt = 0;
   T_Elec_dr03EcalSumEt = 0;
   T_Elec_dr03HcalSumEt = 0;
   T_Elec_isEB = 0;
   T_Elec_isEE = 0;
   T_Elec_MVA = 0;
   T_Elec_simpleEleId80 = 0;
   T_Elec_chargedHadronIso = 0;
   T_Elec_neutralHadronIso = 0;
   T_Elec_photonIso = 0;
   T_Elec_puChargedHadronIso = 0;
   T_Elec_passConversionVeto = 0;
   T_Elec_sigmaIetaIeta = 0;
   T_Elec_deltaPhiIn = 0;
   T_Elec_deltaEtaIn = 0;
   T_Elec_isEcalDriven = 0;
   T_Elec_HtoE = 0;
   T_Elec_vz = 0;
   T_Elec_vy = 0;
   T_Elec_vx = 0;
   T_Elec_nLost = 0;
   T_Elec_nHits = 0;
   T_Elec_SC_Et = 0;
   T_Elec_SC_Eta = 0;
   T_Elec_PFElecPt = 0;
   T_Elec_PFElecPx = 0;
   T_Elec_PFElecPy = 0;
   T_Elec_PFElecPz = 0;
   T_Elec_PFElecE = 0;
   T_Elec_isPF = 0;
   T_JetAKCHS_Px = 0;
   T_JetAKCHS_Py = 0;
   T_JetAKCHS_Pz = 0;
   T_JetAKCHS_Et = 0;
   T_JetAKCHS_Eta = 0;
   T_JetAKCHS_Energy = 0;
   T_JetAKCHS_Tag_HighEffTC = 0;
   T_JetAKCHS_Tag_CombSVtx = 0;
   T_JetAKCHS_Tag_CombSVtxMVA = 0;
   T_JetAKCHS_Tag_TauJet = 0;
   T_JetAKCHS_Tag_ImpParMVA = 0;
   T_JetAKCHS_Tag_JetBProb = 0;
   T_JetAKCHS_Tag_JetProb = 0;
   T_JetAKCHS_Tag_HighEffSimpSVtx = 0;
   T_JetAKCHS_Tag_HighPurSimpSVtx = 0;
   T_JetAKCHS_Tag_HighPurTC = 0;
   T_JetAKCHS_Parton_Px = 0;
   T_JetAKCHS_Parton_Py = 0;
   T_JetAKCHS_Parton_Pz = 0;
   T_JetAKCHS_Parton_Energy = 0;
   T_JetAKCHS_Parton_Flavour = 0;
   T_JetAKCHS_CharHadEnergyFrac = 0;
   T_JetAKCHS_NeutHadEnergyFrac = 0;
   T_JetAKCHS_CharEmEnergyFrac = 0;
   T_JetAKCHS_NeutEmEnergyFrac = 0;
   T_JetAKCHS_CharHadEnergy = 0;
   T_JetAKCHS_NeutHadEnergy = 0;
   T_JetAKCHS_CharEmEnergy = 0;
   T_JetAKCHS_NeutEmEnergy = 0;
   T_JetAKCHS_MuonMultiplicity = 0;
   T_JetAKCHS_NeutralMultiplicity = 0;
   T_JetAKCHS_ChargedMultiplicity = 0;
   T_JetAKCHS_IDLoose = 0;
   T_JetAKCHS_nDaughters = 0;
   T_JetAKCHS_GenJet_InvisibleE = 0;
   T_JetAKCHS_GenJet_Px = 0;
   T_JetAKCHS_GenJet_Py = 0;
   T_JetAKCHS_GenJet_Pz = 0;
   T_JetAKCHS_GenJet_Et = 0;
   T_JetAKCHS_GenJet_Eta = 0;
   T_JetAKCHS_GenJet_Energy = 0;
   T_JetAKCHS_IsGenJet = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("T_Event_Rho", &T_Event_Rho, &b_T_Event_Rho);
   fChain->SetBranchAddress("T_Event_RhoIso", &T_Event_RhoIso, &b_T_Event_RhoIso);
   fChain->SetBranchAddress("T_EventF_EcalDeadCell", &T_EventF_EcalDeadCell, &b_T_EventF_EcalDeadCell);
   fChain->SetBranchAddress("T_EventF_logErrorTooManyClusters", &T_EventF_logErrorTooManyClusters, &b_T_EventF_logErrorTooManyClusters);
   fChain->SetBranchAddress("T_EventF_trackingFailure", &T_EventF_trackingFailure, &b_T_EventF_trackingFailure);
   fChain->SetBranchAddress("T_EventF_hcalLaser", &T_EventF_hcalLaser, &b_T_EventF_hcalLaser);
   fChain->SetBranchAddress("T_EventF_ecalLaserCorr", &T_EventF_ecalLaserCorr, &b_T_EventF_ecalLaserCorr);
   fChain->SetBranchAddress("T_EventF_toomanystripclus", &T_EventF_toomanystripclus, &b_T_EventF_toomanystripclus);
   fChain->SetBranchAddress("T_EventF_manystripclus", &T_EventF_manystripclus, &b_T_EventF_manystripclus);
   fChain->SetBranchAddress("T_EventF_eeBadSc", &T_EventF_eeBadSc, &b_T_EventF_eeBadSc);
   fChain->SetBranchAddress("T_Event_RunNumber", &T_Event_RunNumber, &b_T_Event_RunNumber);
   fChain->SetBranchAddress("T_Event_EventNumber", &T_Event_EventNumber, &b_T_Event_EventNumber);
   fChain->SetBranchAddress("T_Event_LuminosityBlock", &T_Event_LuminosityBlock, &b_T_Event_LuminosityBlock);
   fChain->SetBranchAddress("T_Event_processID", &T_Event_processID, &b_T_Event_processID);
   fChain->SetBranchAddress("T_Event_nPU", &T_Event_nPU, &b_T_Event_nPU);
   fChain->SetBranchAddress("T_Event_nTruePU", &T_Event_nTruePU, &b_T_Event_nTruePU);
   fChain->SetBranchAddress("T_Event_nPUm", &T_Event_nPUm, &b_T_Event_nPUm);
   fChain->SetBranchAddress("T_Event_nPUp", &T_Event_nPUp, &b_T_Event_nPUp);
   fChain->SetBranchAddress("T_Event_AveNTruePU", &T_Event_AveNTruePU, &b_T_Event_AveNTruePU);
   fChain->SetBranchAddress("T_HLT_Mu8_v16", &T_HLT_Mu8_v16, &b_T_HLT_Mu8_v16);
   fChain->SetBranchAddress("T_HLT_Mu17_v3", &T_HLT_Mu17_v3, &b_T_HLT_Mu17_v3);
   fChain->SetBranchAddress("T_HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v12", &T_HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v12, &b_T_HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v12);
   fChain->SetBranchAddress("T_HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v13", &T_HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v13, &b_T_HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v13);
   fChain->SetBranchAddress("T_HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v14", &T_HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v14, &b_T_HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v14);
   fChain->SetBranchAddress("T_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v3", &T_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v3, &b_T_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v3);
   fChain->SetBranchAddress("T_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v4", &T_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v4, &b_T_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v4);
   fChain->SetBranchAddress("T_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v5", &T_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v5, &b_T_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v5);
   fChain->SetBranchAddress("T_Gen_StopMass", &T_Gen_StopMass, &b_T_Gen_StopMass);
   fChain->SetBranchAddress("T_Gen_Chi0Mass", &T_Gen_Chi0Mass, &b_T_Gen_Chi0Mass);
   fChain->SetBranchAddress("T_Gen_CharginoMass", &T_Gen_CharginoMass, &b_T_Gen_CharginoMass);
   fChain->SetBranchAddress("T_Gen_polWeights", &T_Gen_polWeights, &b_T_Gen_polWeights);
   fChain->SetBranchAddress("T_Gen_Muon_PID", &T_Gen_Muon_PID, &b_T_Gen_Muon_PID);
   fChain->SetBranchAddress("T_Gen_Muon_Px", &T_Gen_Muon_Px, &b_T_Gen_Muon_Px);
   fChain->SetBranchAddress("T_Gen_Muon_Py", &T_Gen_Muon_Py, &b_T_Gen_Muon_Py);
   fChain->SetBranchAddress("T_Gen_Muon_Pz", &T_Gen_Muon_Pz, &b_T_Gen_Muon_Pz);
   fChain->SetBranchAddress("T_Gen_Muon_Energy", &T_Gen_Muon_Energy, &b_T_Gen_Muon_Energy);
   fChain->SetBranchAddress("T_Gen_Muon_MPID", &T_Gen_Muon_MPID, &b_T_Gen_Muon_MPID);
   fChain->SetBranchAddress("T_Gen_Muon_MPx", &T_Gen_Muon_MPx, &b_T_Gen_Muon_MPx);
   fChain->SetBranchAddress("T_Gen_Muon_MPy", &T_Gen_Muon_MPy, &b_T_Gen_Muon_MPy);
   fChain->SetBranchAddress("T_Gen_Muon_MPz", &T_Gen_Muon_MPz, &b_T_Gen_Muon_MPz);
   fChain->SetBranchAddress("T_Gen_Muon_MEnergy", &T_Gen_Muon_MEnergy, &b_T_Gen_Muon_MEnergy);
   fChain->SetBranchAddress("T_Gen_Muon_MSt", &T_Gen_Muon_MSt, &b_T_Gen_Muon_MSt);
   fChain->SetBranchAddress("T_Gen_Elec_PID", &T_Gen_Elec_PID, &b_T_Gen_Elec_PID);
   fChain->SetBranchAddress("T_Gen_Elec_Px", &T_Gen_Elec_Px, &b_T_Gen_Elec_Px);
   fChain->SetBranchAddress("T_Gen_Elec_Py", &T_Gen_Elec_Py, &b_T_Gen_Elec_Py);
   fChain->SetBranchAddress("T_Gen_Elec_Pz", &T_Gen_Elec_Pz, &b_T_Gen_Elec_Pz);
   fChain->SetBranchAddress("T_Gen_Elec_Energy", &T_Gen_Elec_Energy, &b_T_Gen_Elec_Energy);
   fChain->SetBranchAddress("T_Gen_Elec_MPID", &T_Gen_Elec_MPID, &b_T_Gen_Elec_MPID);
   fChain->SetBranchAddress("T_Gen_Elec_MPx", &T_Gen_Elec_MPx, &b_T_Gen_Elec_MPx);
   fChain->SetBranchAddress("T_Gen_Elec_MPy", &T_Gen_Elec_MPy, &b_T_Gen_Elec_MPy);
   fChain->SetBranchAddress("T_Gen_Elec_MPz", &T_Gen_Elec_MPz, &b_T_Gen_Elec_MPz);
   fChain->SetBranchAddress("T_Gen_Elec_MEnergy", &T_Gen_Elec_MEnergy, &b_T_Gen_Elec_MEnergy);
   fChain->SetBranchAddress("T_Gen_Elec_MSt", &T_Gen_Elec_MSt, &b_T_Gen_Elec_MSt);
   fChain->SetBranchAddress("T_Gen_b_PID", &T_Gen_b_PID, &b_T_Gen_b_PID);
   fChain->SetBranchAddress("T_Gen_b_Px", &T_Gen_b_Px, &b_T_Gen_b_Px);
   fChain->SetBranchAddress("T_Gen_b_Py", &T_Gen_b_Py, &b_T_Gen_b_Py);
   fChain->SetBranchAddress("T_Gen_b_Pz", &T_Gen_b_Pz, &b_T_Gen_b_Pz);
   fChain->SetBranchAddress("T_Gen_b_Energy", &T_Gen_b_Energy, &b_T_Gen_b_Energy);
   fChain->SetBranchAddress("T_Gen_b_MPID", &T_Gen_b_MPID, &b_T_Gen_b_MPID);
   fChain->SetBranchAddress("T_Gen_b_MPx", &T_Gen_b_MPx, &b_T_Gen_b_MPx);
   fChain->SetBranchAddress("T_Gen_b_MPy", &T_Gen_b_MPy, &b_T_Gen_b_MPy);
   fChain->SetBranchAddress("T_Gen_b_MPz", &T_Gen_b_MPz, &b_T_Gen_b_MPz);
   fChain->SetBranchAddress("T_Gen_b_MEnergy", &T_Gen_b_MEnergy, &b_T_Gen_b_MEnergy);
   fChain->SetBranchAddress("T_Gen_b_MSt", &T_Gen_b_MSt, &b_T_Gen_b_MSt);
   fChain->SetBranchAddress("T_Gen_MuonSt3_pdgId", &T_Gen_MuonSt3_pdgId, &b_T_Gen_MuonSt3_pdgId);
   fChain->SetBranchAddress("T_Gen_MuonSt3_firstMother", &T_Gen_MuonSt3_firstMother, &b_T_Gen_MuonSt3_firstMother);
   fChain->SetBranchAddress("T_Gen_MuonSt3_i", &T_Gen_MuonSt3_i, &b_T_Gen_MuonSt3_i);
   fChain->SetBranchAddress("T_Gen_MuonSt3_energy", &T_Gen_MuonSt3_energy, &b_T_Gen_MuonSt3_energy);
   fChain->SetBranchAddress("T_Gen_MuonSt3_pt", &T_Gen_MuonSt3_pt, &b_T_Gen_MuonSt3_pt);
   fChain->SetBranchAddress("T_Gen_MuonSt3_eta", &T_Gen_MuonSt3_eta, &b_T_Gen_MuonSt3_eta);
   fChain->SetBranchAddress("T_Gen_MuonSt3_phi", &T_Gen_MuonSt3_phi, &b_T_Gen_MuonSt3_phi);
   fChain->SetBranchAddress("T_Gen_ElecSt3_pdgId", &T_Gen_ElecSt3_pdgId, &b_T_Gen_ElecSt3_pdgId);
   fChain->SetBranchAddress("T_Gen_ElecSt3_firstMother", &T_Gen_ElecSt3_firstMother, &b_T_Gen_ElecSt3_firstMother);
   fChain->SetBranchAddress("T_Gen_ElecSt3_i", &T_Gen_ElecSt3_i, &b_T_Gen_ElecSt3_i);
   fChain->SetBranchAddress("T_Gen_ElecSt3_energy", &T_Gen_ElecSt3_energy, &b_T_Gen_ElecSt3_energy);
   fChain->SetBranchAddress("T_Gen_ElecSt3_pt", &T_Gen_ElecSt3_pt, &b_T_Gen_ElecSt3_pt);
   fChain->SetBranchAddress("T_Gen_ElecSt3_eta", &T_Gen_ElecSt3_eta, &b_T_Gen_ElecSt3_eta);
   fChain->SetBranchAddress("T_Gen_ElecSt3_phi", &T_Gen_ElecSt3_phi, &b_T_Gen_ElecSt3_phi);
   fChain->SetBranchAddress("T_Gen_TauSt3_pdgId", &T_Gen_TauSt3_pdgId, &b_T_Gen_TauSt3_pdgId);
   fChain->SetBranchAddress("T_Gen_TauSt3_firstMother", &T_Gen_TauSt3_firstMother, &b_T_Gen_TauSt3_firstMother);
   fChain->SetBranchAddress("T_Gen_TauSt3_i", &T_Gen_TauSt3_i, &b_T_Gen_TauSt3_i);
   fChain->SetBranchAddress("T_Gen_TauSt3_energy", &T_Gen_TauSt3_energy, &b_T_Gen_TauSt3_energy);
   fChain->SetBranchAddress("T_Gen_TauSt3_pt", &T_Gen_TauSt3_pt, &b_T_Gen_TauSt3_pt);
   fChain->SetBranchAddress("T_Gen_TauSt3_eta", &T_Gen_TauSt3_eta, &b_T_Gen_TauSt3_eta);
   fChain->SetBranchAddress("T_Gen_TauSt3_phi", &T_Gen_TauSt3_phi, &b_T_Gen_TauSt3_phi);
   fChain->SetBranchAddress("T_Gen_StopSt3_pdgId", &T_Gen_StopSt3_pdgId, &b_T_Gen_StopSt3_pdgId);
   fChain->SetBranchAddress("T_Gen_StopSt3_firstMother", &T_Gen_StopSt3_firstMother, &b_T_Gen_StopSt3_firstMother);
   fChain->SetBranchAddress("T_Gen_StopSt3_i", &T_Gen_StopSt3_i, &b_T_Gen_StopSt3_i);
   fChain->SetBranchAddress("T_Gen_StopSt3_energy", &T_Gen_StopSt3_energy, &b_T_Gen_StopSt3_energy);
   fChain->SetBranchAddress("T_Gen_StopSt3_pt", &T_Gen_StopSt3_pt, &b_T_Gen_StopSt3_pt);
   fChain->SetBranchAddress("T_Gen_StopSt3_eta", &T_Gen_StopSt3_eta, &b_T_Gen_StopSt3_eta);
   fChain->SetBranchAddress("T_Gen_StopSt3_phi", &T_Gen_StopSt3_phi, &b_T_Gen_StopSt3_phi);
   fChain->SetBranchAddress("T_Gen_Chi0St3_pdgId", &T_Gen_Chi0St3_pdgId, &b_T_Gen_Chi0St3_pdgId);
   fChain->SetBranchAddress("T_Gen_Chi0St3_firstMother", &T_Gen_Chi0St3_firstMother, &b_T_Gen_Chi0St3_firstMother);
   fChain->SetBranchAddress("T_Gen_Chi0St3_i", &T_Gen_Chi0St3_i, &b_T_Gen_Chi0St3_i);
   fChain->SetBranchAddress("T_Gen_Chi0St3_energy", &T_Gen_Chi0St3_energy, &b_T_Gen_Chi0St3_energy);
   fChain->SetBranchAddress("T_Gen_Chi0St3_pt", &T_Gen_Chi0St3_pt, &b_T_Gen_Chi0St3_pt);
   fChain->SetBranchAddress("T_Gen_Chi0St3_eta", &T_Gen_Chi0St3_eta, &b_T_Gen_Chi0St3_eta);
   fChain->SetBranchAddress("T_Gen_Chi0St3_phi", &T_Gen_Chi0St3_phi, &b_T_Gen_Chi0St3_phi);
   fChain->SetBranchAddress("T_Gen_ChiPMSt3_pdgId", &T_Gen_ChiPMSt3_pdgId, &b_T_Gen_ChiPMSt3_pdgId);
   fChain->SetBranchAddress("T_Gen_ChiPMSt3_firstMother", &T_Gen_ChiPMSt3_firstMother, &b_T_Gen_ChiPMSt3_firstMother);
   fChain->SetBranchAddress("T_Gen_ChiPMSt3_i", &T_Gen_ChiPMSt3_i, &b_T_Gen_ChiPMSt3_i);
   fChain->SetBranchAddress("T_Gen_ChiPMSt3_energy", &T_Gen_ChiPMSt3_energy, &b_T_Gen_ChiPMSt3_energy);
   fChain->SetBranchAddress("T_Gen_ChiPMSt3_pt", &T_Gen_ChiPMSt3_pt, &b_T_Gen_ChiPMSt3_pt);
   fChain->SetBranchAddress("T_Gen_ChiPMSt3_eta", &T_Gen_ChiPMSt3_eta, &b_T_Gen_ChiPMSt3_eta);
   fChain->SetBranchAddress("T_Gen_ChiPMSt3_phi", &T_Gen_ChiPMSt3_phi, &b_T_Gen_ChiPMSt3_phi);
   fChain->SetBranchAddress("T_Gen_bSt3_pdgId", &T_Gen_bSt3_pdgId, &b_T_Gen_bSt3_pdgId);
   fChain->SetBranchAddress("T_Gen_bSt3_firstMother", &T_Gen_bSt3_firstMother, &b_T_Gen_bSt3_firstMother);
   fChain->SetBranchAddress("T_Gen_bSt3_i", &T_Gen_bSt3_i, &b_T_Gen_bSt3_i);
   fChain->SetBranchAddress("T_Gen_bSt3_energy", &T_Gen_bSt3_energy, &b_T_Gen_bSt3_energy);
   fChain->SetBranchAddress("T_Gen_bSt3_pt", &T_Gen_bSt3_pt, &b_T_Gen_bSt3_pt);
   fChain->SetBranchAddress("T_Gen_bSt3_eta", &T_Gen_bSt3_eta, &b_T_Gen_bSt3_eta);
   fChain->SetBranchAddress("T_Gen_bSt3_phi", &T_Gen_bSt3_phi, &b_T_Gen_bSt3_phi);
   fChain->SetBranchAddress("T_Gen_tSt3_pdgId", &T_Gen_tSt3_pdgId, &b_T_Gen_tSt3_pdgId);
   fChain->SetBranchAddress("T_Gen_tSt3_firstMother", &T_Gen_tSt3_firstMother, &b_T_Gen_tSt3_firstMother);
   fChain->SetBranchAddress("T_Gen_tSt3_i", &T_Gen_tSt3_i, &b_T_Gen_tSt3_i);
   fChain->SetBranchAddress("T_Gen_tSt3_energy", &T_Gen_tSt3_energy, &b_T_Gen_tSt3_energy);
   fChain->SetBranchAddress("T_Gen_tSt3_pt", &T_Gen_tSt3_pt, &b_T_Gen_tSt3_pt);
   fChain->SetBranchAddress("T_Gen_tSt3_eta", &T_Gen_tSt3_eta, &b_T_Gen_tSt3_eta);
   fChain->SetBranchAddress("T_Gen_tSt3_phi", &T_Gen_tSt3_phi, &b_T_Gen_tSt3_phi);
   fChain->SetBranchAddress("T_Gen_WSt3_pdgId", &T_Gen_WSt3_pdgId, &b_T_Gen_WSt3_pdgId);
   fChain->SetBranchAddress("T_Gen_WSt3_firstMother", &T_Gen_WSt3_firstMother, &b_T_Gen_WSt3_firstMother);
   fChain->SetBranchAddress("T_Gen_WSt3_i", &T_Gen_WSt3_i, &b_T_Gen_WSt3_i);
   fChain->SetBranchAddress("T_Gen_WSt3_energy", &T_Gen_WSt3_energy, &b_T_Gen_WSt3_energy);
   fChain->SetBranchAddress("T_Gen_WSt3_pt", &T_Gen_WSt3_pt, &b_T_Gen_WSt3_pt);
   fChain->SetBranchAddress("T_Gen_WSt3_eta", &T_Gen_WSt3_eta, &b_T_Gen_WSt3_eta);
   fChain->SetBranchAddress("T_Gen_WSt3_phi", &T_Gen_WSt3_phi, &b_T_Gen_WSt3_phi);
   fChain->SetBranchAddress("T_Gen_TauSt3_IsLepDec", &T_Gen_TauSt3_IsLepDec, &b_T_Gen_TauSt3_IsLepDec);
   fChain->SetBranchAddress("T_Gen_TauSt3_LepDec_PID", &T_Gen_TauSt3_LepDec_PID, &b_T_Gen_TauSt3_LepDec_PID);
   fChain->SetBranchAddress("T_Gen_TauSt3_LepDec_Px", &T_Gen_TauSt3_LepDec_Px, &b_T_Gen_TauSt3_LepDec_Px);
   fChain->SetBranchAddress("T_Gen_TauSt3_LepDec_Py", &T_Gen_TauSt3_LepDec_Py, &b_T_Gen_TauSt3_LepDec_Py);
   fChain->SetBranchAddress("T_Gen_TauSt3_LepDec_Pz", &T_Gen_TauSt3_LepDec_Pz, &b_T_Gen_TauSt3_LepDec_Pz);
   fChain->SetBranchAddress("T_Gen_TauSt3_LepDec_Energy", &T_Gen_TauSt3_LepDec_Energy, &b_T_Gen_TauSt3_LepDec_Energy);
   fChain->SetBranchAddress("T_Muon_Eta", &T_Muon_Eta, &b_T_Muon_Eta);
   fChain->SetBranchAddress("T_Muon_IsGlobalMuon", &T_Muon_IsGlobalMuon, &b_T_Muon_IsGlobalMuon);
   fChain->SetBranchAddress("T_Muon_IsAllTrackerMuons", &T_Muon_IsAllTrackerMuons, &b_T_Muon_IsAllTrackerMuons);
   fChain->SetBranchAddress("T_Muon_IsTrackerMuonArbitrated", &T_Muon_IsTrackerMuonArbitrated, &b_T_Muon_IsTrackerMuonArbitrated);
   fChain->SetBranchAddress("T_Muon_IsGMPTMuons", &T_Muon_IsGMPTMuons, &b_T_Muon_IsGMPTMuons);
   fChain->SetBranchAddress("T_Muon_SegmentCompatibility", &T_Muon_SegmentCompatibility, &b_T_Muon_SegmentCompatibility);
   fChain->SetBranchAddress("T_Muon_trkKink", &T_Muon_trkKink, &b_T_Muon_trkKink);
   fChain->SetBranchAddress("T_Muon_Px", &T_Muon_Px, &b_T_Muon_Px);
   fChain->SetBranchAddress("T_Muon_Py", &T_Muon_Py, &b_T_Muon_Py);
   fChain->SetBranchAddress("T_Muon_Pz", &T_Muon_Pz, &b_T_Muon_Pz);
   fChain->SetBranchAddress("T_Muon_Pt", &T_Muon_Pt, &b_T_Muon_Pt);
   fChain->SetBranchAddress("T_Muon_deltaPt", &T_Muon_deltaPt, &b_T_Muon_deltaPt);
   fChain->SetBranchAddress("T_Muon_Energy", &T_Muon_Energy, &b_T_Muon_Energy);
   fChain->SetBranchAddress("T_Muon_Charge", &T_Muon_Charge, &b_T_Muon_Charge);
   fChain->SetBranchAddress("T_Muon_NormChi2GTrk", &T_Muon_NormChi2GTrk, &b_T_Muon_NormChi2GTrk);
   fChain->SetBranchAddress("T_Muon_NValidHitsInTrk", &T_Muon_NValidHitsInTrk, &b_T_Muon_NValidHitsInTrk);
   fChain->SetBranchAddress("T_Muon_NValidPixelHitsInTrk", &T_Muon_NValidPixelHitsInTrk, &b_T_Muon_NValidPixelHitsInTrk);
   fChain->SetBranchAddress("T_Muon_InnerTrackFound", &T_Muon_InnerTrackFound, &b_T_Muon_InnerTrackFound);
   fChain->SetBranchAddress("T_Muon_NValidHitsSATrk", &T_Muon_NValidHitsSATrk, &b_T_Muon_NValidHitsSATrk);
   fChain->SetBranchAddress("T_Muon_NValidHitsGTrk", &T_Muon_NValidHitsGTrk, &b_T_Muon_NValidHitsGTrk);
   fChain->SetBranchAddress("T_Muon_Chi2InTrk", &T_Muon_Chi2InTrk, &b_T_Muon_Chi2InTrk);
   fChain->SetBranchAddress("T_Muon_dofInTrk", &T_Muon_dofInTrk, &b_T_Muon_dofInTrk);
   fChain->SetBranchAddress("T_Muon_IPAbsGTrack", &T_Muon_IPAbsGTrack, &b_T_Muon_IPAbsGTrack);
   fChain->SetBranchAddress("T_Muon_IPAbsInTrack", &T_Muon_IPAbsInTrack, &b_T_Muon_IPAbsInTrack);
   fChain->SetBranchAddress("T_Muon_IPwrtAveBSInTrack", &T_Muon_IPwrtAveBSInTrack, &b_T_Muon_IPwrtAveBSInTrack);
   fChain->SetBranchAddress("T_Muon_chargedHadronIsoR04", &T_Muon_chargedHadronIsoR04, &b_T_Muon_chargedHadronIsoR04);
   fChain->SetBranchAddress("T_Muon_neutralHadronIsoR04", &T_Muon_neutralHadronIsoR04, &b_T_Muon_neutralHadronIsoR04);
   fChain->SetBranchAddress("T_Muon_photonIsoR04", &T_Muon_photonIsoR04, &b_T_Muon_photonIsoR04);
   fChain->SetBranchAddress("T_Muon_chargedParticleIsoR03", &T_Muon_chargedParticleIsoR03, &b_T_Muon_chargedParticleIsoR03);
   fChain->SetBranchAddress("T_Muon_chargedHadronIsoR03", &T_Muon_chargedHadronIsoR03, &b_T_Muon_chargedHadronIsoR03);
   fChain->SetBranchAddress("T_Muon_neutralHadronIsoR03", &T_Muon_neutralHadronIsoR03, &b_T_Muon_neutralHadronIsoR03);
   fChain->SetBranchAddress("T_Muon_photonIsoR03", &T_Muon_photonIsoR03, &b_T_Muon_photonIsoR03);
   fChain->SetBranchAddress("T_Muon_sumPUPtR04", &T_Muon_sumPUPtR04, &b_T_Muon_sumPUPtR04);
   fChain->SetBranchAddress("T_Muon_sumPUPtR03", &T_Muon_sumPUPtR03, &b_T_Muon_sumPUPtR03);
   fChain->SetBranchAddress("T_Muon_vz", &T_Muon_vz, &b_T_Muon_vz);
   fChain->SetBranchAddress("T_Muon_vy", &T_Muon_vy, &b_T_Muon_vy);
   fChain->SetBranchAddress("T_Muon_vx", &T_Muon_vx, &b_T_Muon_vx);
   fChain->SetBranchAddress("T_Muon_NumOfMatchedStations", &T_Muon_NumOfMatchedStations, &b_T_Muon_NumOfMatchedStations);
   fChain->SetBranchAddress("T_Muon_PFMuonPt", &T_Muon_PFMuonPt, &b_T_Muon_PFMuonPt);
   fChain->SetBranchAddress("T_Muon_PFMuonPx", &T_Muon_PFMuonPx, &b_T_Muon_PFMuonPx);
   fChain->SetBranchAddress("T_Muon_PFMuonPy", &T_Muon_PFMuonPy, &b_T_Muon_PFMuonPy);
   fChain->SetBranchAddress("T_Muon_PFMuonPz", &T_Muon_PFMuonPz, &b_T_Muon_PFMuonPz);
   fChain->SetBranchAddress("T_Muon_PFMuonE", &T_Muon_PFMuonE, &b_T_Muon_PFMuonE);
   fChain->SetBranchAddress("T_Muon_isPFMuon", &T_Muon_isPFMuon, &b_T_Muon_isPFMuon);
   fChain->SetBranchAddress("T_Muon_NLayers", &T_Muon_NLayers, &b_T_Muon_NLayers);
   fChain->SetBranchAddress("T_Vertex_z", &T_Vertex_z, &b_T_Vertex_z);
   fChain->SetBranchAddress("T_Vertex_y", &T_Vertex_y, &b_T_Vertex_y);
   fChain->SetBranchAddress("T_Vertex_x", &T_Vertex_x, &b_T_Vertex_x);
   fChain->SetBranchAddress("T_Vertex_Chi2Prob", &T_Vertex_Chi2Prob, &b_T_Vertex_Chi2Prob);
   fChain->SetBranchAddress("T_Vertex_rho", &T_Vertex_rho, &b_T_Vertex_rho);
   fChain->SetBranchAddress("T_Vertex_ndof", &T_Vertex_ndof, &b_T_Vertex_ndof);
   fChain->SetBranchAddress("T_Vertex_isFake", &T_Vertex_isFake, &b_T_Vertex_isFake);
   fChain->SetBranchAddress("T_Vertex_tracksSize", &T_Vertex_tracksSize, &b_T_Vertex_tracksSize);
   fChain->SetBranchAddress("T_Elec_Eta", &T_Elec_Eta, &b_T_Elec_Eta);
   fChain->SetBranchAddress("T_Elec_IPwrtAveBS", &T_Elec_IPwrtAveBS, &b_T_Elec_IPwrtAveBS);
   fChain->SetBranchAddress("T_Elec_IPwrtPV", &T_Elec_IPwrtPV, &b_T_Elec_IPwrtPV);
   fChain->SetBranchAddress("T_Elec_dzwrtPV", &T_Elec_dzwrtPV, &b_T_Elec_dzwrtPV);
   fChain->SetBranchAddress("T_Elec_Px", &T_Elec_Px, &b_T_Elec_Px);
   fChain->SetBranchAddress("T_Elec_Py", &T_Elec_Py, &b_T_Elec_Py);
   fChain->SetBranchAddress("T_Elec_Pz", &T_Elec_Pz, &b_T_Elec_Pz);
   fChain->SetBranchAddress("T_Elec_Pt", &T_Elec_Pt, &b_T_Elec_Pt);
   fChain->SetBranchAddress("T_Elec_Energy", &T_Elec_Energy, &b_T_Elec_Energy);
   fChain->SetBranchAddress("T_Elec_Charge", &T_Elec_Charge, &b_T_Elec_Charge);
   fChain->SetBranchAddress("T_Elec_nBrems", &T_Elec_nBrems, &b_T_Elec_nBrems);
   fChain->SetBranchAddress("T_Elec_fBrem", &T_Elec_fBrem, &b_T_Elec_fBrem);
   fChain->SetBranchAddress("T_Elec_eSuperClusterOverP", &T_Elec_eSuperClusterOverP, &b_T_Elec_eSuperClusterOverP);
   fChain->SetBranchAddress("T_Elec_ecalEnergy", &T_Elec_ecalEnergy, &b_T_Elec_ecalEnergy);
   fChain->SetBranchAddress("T_Elec_dr03TkSumPt", &T_Elec_dr03TkSumPt, &b_T_Elec_dr03TkSumPt);
   fChain->SetBranchAddress("T_Elec_dr03EcalSumEt", &T_Elec_dr03EcalSumEt, &b_T_Elec_dr03EcalSumEt);
   fChain->SetBranchAddress("T_Elec_dr03HcalSumEt", &T_Elec_dr03HcalSumEt, &b_T_Elec_dr03HcalSumEt);
   fChain->SetBranchAddress("T_Elec_isEB", &T_Elec_isEB, &b_T_Elec_isEB);
   fChain->SetBranchAddress("T_Elec_isEE", &T_Elec_isEE, &b_T_Elec_isEE);
   fChain->SetBranchAddress("T_Elec_MVA", &T_Elec_MVA, &b_T_Elec_MVA);
   fChain->SetBranchAddress("T_Elec_simpleEleId80", &T_Elec_simpleEleId80, &b_T_Elec_simpleEleId80);
   fChain->SetBranchAddress("T_Elec_chargedHadronIso", &T_Elec_chargedHadronIso, &b_T_Elec_chargedHadronIso);
   fChain->SetBranchAddress("T_Elec_neutralHadronIso", &T_Elec_neutralHadronIso, &b_T_Elec_neutralHadronIso);
   fChain->SetBranchAddress("T_Elec_photonIso", &T_Elec_photonIso, &b_T_Elec_photonIso);
   fChain->SetBranchAddress("T_Elec_puChargedHadronIso", &T_Elec_puChargedHadronIso, &b_T_Elec_puChargedHadronIso);
   fChain->SetBranchAddress("T_Elec_passConversionVeto", &T_Elec_passConversionVeto, &b_T_Elec_passConversionVeto);
   fChain->SetBranchAddress("T_Elec_sigmaIetaIeta", &T_Elec_sigmaIetaIeta, &b_T_Elec_sigmaIetaIeta);
   fChain->SetBranchAddress("T_Elec_deltaPhiIn", &T_Elec_deltaPhiIn, &b_T_Elec_deltaPhiIn);
   fChain->SetBranchAddress("T_Elec_deltaEtaIn", &T_Elec_deltaEtaIn, &b_T_Elec_deltaEtaIn);
   fChain->SetBranchAddress("T_Elec_isEcalDriven", &T_Elec_isEcalDriven, &b_T_Elec_isEcalDriven);
   fChain->SetBranchAddress("T_Elec_HtoE", &T_Elec_HtoE, &b_T_Elec_HtoE);
   fChain->SetBranchAddress("T_Elec_vz", &T_Elec_vz, &b_T_Elec_vz);
   fChain->SetBranchAddress("T_Elec_vy", &T_Elec_vy, &b_T_Elec_vy);
   fChain->SetBranchAddress("T_Elec_vx", &T_Elec_vx, &b_T_Elec_vx);
   fChain->SetBranchAddress("T_Elec_nLost", &T_Elec_nLost, &b_T_Elec_nLost);
   fChain->SetBranchAddress("T_Elec_nHits", &T_Elec_nHits, &b_T_Elec_nHits);
   fChain->SetBranchAddress("T_Elec_SC_Et", &T_Elec_SC_Et, &b_T_Elec_SC_Et);
   fChain->SetBranchAddress("T_Elec_SC_Eta", &T_Elec_SC_Eta, &b_T_Elec_SC_Eta);
   fChain->SetBranchAddress("T_Elec_PFElecPt", &T_Elec_PFElecPt, &b_T_Elec_PFElecPt);
   fChain->SetBranchAddress("T_Elec_PFElecPx", &T_Elec_PFElecPx, &b_T_Elec_PFElecPx);
   fChain->SetBranchAddress("T_Elec_PFElecPy", &T_Elec_PFElecPy, &b_T_Elec_PFElecPy);
   fChain->SetBranchAddress("T_Elec_PFElecPz", &T_Elec_PFElecPz, &b_T_Elec_PFElecPz);
   fChain->SetBranchAddress("T_Elec_PFElecE", &T_Elec_PFElecE, &b_T_Elec_PFElecE);
   fChain->SetBranchAddress("T_Elec_isPF", &T_Elec_isPF, &b_T_Elec_isPF);
   fChain->SetBranchAddress("T_JetAKCHS_Px", &T_JetAKCHS_Px, &b_T_JetAKCHS_Px);
   fChain->SetBranchAddress("T_JetAKCHS_Py", &T_JetAKCHS_Py, &b_T_JetAKCHS_Py);
   fChain->SetBranchAddress("T_JetAKCHS_Pz", &T_JetAKCHS_Pz, &b_T_JetAKCHS_Pz);
   fChain->SetBranchAddress("T_JetAKCHS_Et", &T_JetAKCHS_Et, &b_T_JetAKCHS_Et);
   fChain->SetBranchAddress("T_JetAKCHS_Eta", &T_JetAKCHS_Eta, &b_T_JetAKCHS_Eta);
   fChain->SetBranchAddress("T_JetAKCHS_Energy", &T_JetAKCHS_Energy, &b_T_JetAKCHS_Energy);
   fChain->SetBranchAddress("T_JetAKCHS_Tag_HighEffTC", &T_JetAKCHS_Tag_HighEffTC, &b_T_JetAKCHS_Tag_HighEffTC);
   fChain->SetBranchAddress("T_JetAKCHS_Tag_CombSVtx", &T_JetAKCHS_Tag_CombSVtx, &b_T_JetAKCHS_Tag_CombSVtx);
   fChain->SetBranchAddress("T_JetAKCHS_Tag_CombSVtxMVA", &T_JetAKCHS_Tag_CombSVtxMVA, &b_T_JetAKCHS_Tag_CombSVtxMVA);
   fChain->SetBranchAddress("T_JetAKCHS_Tag_TauJet", &T_JetAKCHS_Tag_TauJet, &b_T_JetAKCHS_Tag_TauJet);
   fChain->SetBranchAddress("T_JetAKCHS_Tag_ImpParMVA", &T_JetAKCHS_Tag_ImpParMVA, &b_T_JetAKCHS_Tag_ImpParMVA);
   fChain->SetBranchAddress("T_JetAKCHS_Tag_JetBProb", &T_JetAKCHS_Tag_JetBProb, &b_T_JetAKCHS_Tag_JetBProb);
   fChain->SetBranchAddress("T_JetAKCHS_Tag_JetProb", &T_JetAKCHS_Tag_JetProb, &b_T_JetAKCHS_Tag_JetProb);
   fChain->SetBranchAddress("T_JetAKCHS_Tag_HighEffSimpSVtx", &T_JetAKCHS_Tag_HighEffSimpSVtx, &b_T_JetAKCHS_Tag_HighEffSimpSVtx);
   fChain->SetBranchAddress("T_JetAKCHS_Tag_HighPurSimpSVtx", &T_JetAKCHS_Tag_HighPurSimpSVtx, &b_T_JetAKCHS_Tag_HighPurSimpSVtx);
   fChain->SetBranchAddress("T_JetAKCHS_Tag_HighPurTC", &T_JetAKCHS_Tag_HighPurTC, &b_T_JetAKCHS_Tag_HighPurTC);
   fChain->SetBranchAddress("T_JetAKCHS_Parton_Px", &T_JetAKCHS_Parton_Px, &b_T_JetAKCHS_Parton_Px);
   fChain->SetBranchAddress("T_JetAKCHS_Parton_Py", &T_JetAKCHS_Parton_Py, &b_T_JetAKCHS_Parton_Py);
   fChain->SetBranchAddress("T_JetAKCHS_Parton_Pz", &T_JetAKCHS_Parton_Pz, &b_T_JetAKCHS_Parton_Pz);
   fChain->SetBranchAddress("T_JetAKCHS_Parton_Energy", &T_JetAKCHS_Parton_Energy, &b_T_JetAKCHS_Parton_Energy);
   fChain->SetBranchAddress("T_JetAKCHS_Parton_Flavour", &T_JetAKCHS_Parton_Flavour, &b_T_JetAKCHS_Parton_Flavour);
   fChain->SetBranchAddress("T_JetAKCHS_CharHadEnergyFrac", &T_JetAKCHS_CharHadEnergyFrac, &b_T_JetAKCHS_CharHadEnergyFrac);
   fChain->SetBranchAddress("T_JetAKCHS_NeutHadEnergyFrac", &T_JetAKCHS_NeutHadEnergyFrac, &b_T_JetAKCHS_NeutHadEnergyFrac);
   fChain->SetBranchAddress("T_JetAKCHS_CharEmEnergyFrac", &T_JetAKCHS_CharEmEnergyFrac, &b_T_JetAKCHS_CharEmEnergyFrac);
   fChain->SetBranchAddress("T_JetAKCHS_NeutEmEnergyFrac", &T_JetAKCHS_NeutEmEnergyFrac, &b_T_JetAKCHS_NeutEmEnergyFrac);
   fChain->SetBranchAddress("T_JetAKCHS_CharHadEnergy", &T_JetAKCHS_CharHadEnergy, &b_T_JetAKCHS_CharHadEnergy);
   fChain->SetBranchAddress("T_JetAKCHS_NeutHadEnergy", &T_JetAKCHS_NeutHadEnergy, &b_T_JetAKCHS_NeutHadEnergy);
   fChain->SetBranchAddress("T_JetAKCHS_CharEmEnergy", &T_JetAKCHS_CharEmEnergy, &b_T_JetAKCHS_CharEmEnergy);
   fChain->SetBranchAddress("T_JetAKCHS_NeutEmEnergy", &T_JetAKCHS_NeutEmEnergy, &b_T_JetAKCHS_NeutEmEnergy);
   fChain->SetBranchAddress("T_JetAKCHS_MuonMultiplicity", &T_JetAKCHS_MuonMultiplicity, &b_T_JetAKCHS_MuonMultiplicity);
   fChain->SetBranchAddress("T_JetAKCHS_NeutralMultiplicity", &T_JetAKCHS_NeutralMultiplicity, &b_T_JetAKCHS_NeutralMultiplicity);
   fChain->SetBranchAddress("T_JetAKCHS_ChargedMultiplicity", &T_JetAKCHS_ChargedMultiplicity, &b_T_JetAKCHS_ChargedMultiplicity);
   fChain->SetBranchAddress("T_JetAKCHS_IDLoose", &T_JetAKCHS_IDLoose, &b_T_JetAKCHS_IDLoose);
   fChain->SetBranchAddress("T_JetAKCHS_nDaughters", &T_JetAKCHS_nDaughters, &b_T_JetAKCHS_nDaughters);
   fChain->SetBranchAddress("T_JetAKCHS_GenJet_InvisibleE", &T_JetAKCHS_GenJet_InvisibleE, &b_T_JetAKCHS_GenJet_InvisibleE);
   fChain->SetBranchAddress("T_JetAKCHS_GenJet_Px", &T_JetAKCHS_GenJet_Px, &b_T_JetAKCHS_GenJet_Px);
   fChain->SetBranchAddress("T_JetAKCHS_GenJet_Py", &T_JetAKCHS_GenJet_Py, &b_T_JetAKCHS_GenJet_Py);
   fChain->SetBranchAddress("T_JetAKCHS_GenJet_Pz", &T_JetAKCHS_GenJet_Pz, &b_T_JetAKCHS_GenJet_Pz);
   fChain->SetBranchAddress("T_JetAKCHS_GenJet_Et", &T_JetAKCHS_GenJet_Et, &b_T_JetAKCHS_GenJet_Et);
   fChain->SetBranchAddress("T_JetAKCHS_GenJet_Eta", &T_JetAKCHS_GenJet_Eta, &b_T_JetAKCHS_GenJet_Eta);
   fChain->SetBranchAddress("T_JetAKCHS_GenJet_Energy", &T_JetAKCHS_GenJet_Energy, &b_T_JetAKCHS_GenJet_Energy);
   fChain->SetBranchAddress("T_JetAKCHS_IsGenJet", &T_JetAKCHS_IsGenJet, &b_T_JetAKCHS_IsGenJet);
   fChain->SetBranchAddress("T_METPF_ET", &T_METPF_ET, &b_T_METPF_ET);
   fChain->SetBranchAddress("T_METPF_Phi", &T_METPF_Phi, &b_T_METPF_Phi);
   fChain->SetBranchAddress("T_METPF_Sig", &T_METPF_Sig, &b_T_METPF_Sig);
   fChain->SetBranchAddress("T_METPFTypeI_ET", &T_METPFTypeI_ET, &b_T_METPFTypeI_ET);
   fChain->SetBranchAddress("T_METPFTypeI_Phi", &T_METPFTypeI_Phi, &b_T_METPFTypeI_Phi);
   fChain->SetBranchAddress("T_METPFTypeI_Sig", &T_METPFTypeI_Sig, &b_T_METPFTypeI_Sig);
   fChain->SetBranchAddress("T_METgen_ET", &T_METgen_ET, &b_T_METgen_ET);
   fChain->SetBranchAddress("T_METgen_Phi", &T_METgen_Phi, &b_T_METgen_Phi);
   fChain->SetBranchAddress("T_passTriggerDoubleMu", &T_passTriggerDoubleMu, &b_T_passTriggerDoubleMu);
   fChain->SetBranchAddress("T_passTriggerDoubleEl", &T_passTriggerDoubleEl, &b_T_passTriggerDoubleEl);
   fChain->SetBranchAddress("T_passTriggerElMu", &T_passTriggerElMu, &b_T_passTriggerElMu);

}

Bool_t GeneralSkimmer::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

Double_t GetEffectiveArea(float eta);

#endif // #ifdef GeneralSkimmer_cxx
