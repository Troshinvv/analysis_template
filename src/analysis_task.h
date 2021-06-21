//
// Created by mikhail on 6/16/20.
//

#ifndef QUALITY_ASSURANCE_SRC_TREE_READER_H_
#define QUALITY_ASSURANCE_SRC_TREE_READER_H_

#include <TFile.h>
#include <TTree.h>
#include <TH2F.h>

#include <AnalysisTree/FillTask.hpp>
#include <AnalysisTree/Detector.hpp>
#include <AnalysisTree/EventHeader.hpp>
#include <AnalysisTree/DataHeader.hpp>
#include <AnalysisTree/Matching.hpp>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TProfile.h>
#include <memory>
#include <string>
#include <math.h>
#include <Math/SpecFuncMathMore.h>

class AnalysisTask : public AnalysisTree::FillTask{
public:
 AnalysisTask() = default;
  ~AnalysisTask() override = default;
  void Init( std::map<std::string, void*>& branch_map ) override;
  void Exec() override;
  void Finish() override;

private:
  /* pointers to link tree's branches with */
  AnalysisTree::EventHeader* event_header_{nullptr}; 		// event info
  AnalysisTree::TrackDetector* vtx_tracks_{nullptr}; 		        // reconstructed tracks
  AnalysisTree::Particles* mc_tracks_{nullptr};                             //Mc tracks
  AnalysisTree::ModuleDetector* fhcal_modules_{nullptr}; 		// modules of FhCal branch
  AnalysisTree::ModulePositions fhcal_modules_positions_;
  
  AnalysisTree::Matching* tpc_mc_matching{nullptr};
  TH1F* pT_distribution_;
  TH1F* fhcal_energy_distribution_;
  TH2F*fhcal_modules_xy_;
  TH2F* TpcpT_vs_eta_proton;
  TH2F* TpcpT_vs_eta_kaon;
TH2F* TpcpT_vs_eta_pion;
TH2F* TpcpT_vs_eta_all;

  TH2F* TpcpT_vs_phi;
  TH2F* Tpcphi_vs_eta;
  TH2F* McpT_vs_eta_proton;
  TH2F* McpT_vs_eta_kaon;
TH2F* McpT_vs_eta_pion;
TH2F* McpT_vs_eta_all;


  TH2F* McpT_vs_phi;
  TH2F* Mcphi_vs_eta;
  TH2F* Energy_vs_moduleId;
  TH1F *Qxall;
  TH1F *Qyall;
  TH1F *Qx44;
  TH1F *Qx90;
  TH1F *Qy44;
  TH1F *Qy90;
  TH2F *ntracksxb;
  TH1F *fnall;
  TH1F *fn44;
  TH1F *fn90;
  TH1F *Qx44all[17];
  TH1F *Qy44all[17];
  TH1F *Qx90all[17];
  TH1F *Qy90all[17];
  TH1F *fn44all[17];
  TH1F *fn90all[17];
  TH1F *Bt;
  TH1F *Btt;
  TH2F *energyvsb;
  TH2F *gr;
  TH2F *gr2;
  TH1F *PID_proton;
  TH2F *dEdX_vs_pz_proton;
  TH2F *mass2_vs_pz_proton;
  TH2F *dEdX_vs_pz_kaon;
  TH2F *mass2_vs_pz_kaon;
  TH2F *dEdX_vs_pz_pion;
  TH2F *mass2_vs_pz_pion;
  TH2F *dEdX_vs_pz_all;

    float ResPhin[8]={0,0,0,0,0,0,0,0};
   float ResN[8]={0,0,0,0,0,0,0,0};
   TProfile* TpcflowvsB[3];
   TProfile* TpcPhivsB[3];
   TProfile* McPhivsB[3];
   TProfile* McflowvsB[3];
 //  TProfile* simflowvsB;

   TH1F *PhiModule;
   TProfile* ResRPS;
   TProfile* ResRPN;
   TProfile* ResRP3;
   TProfile* ResHalf;
   TProfile* TpcflowvspT[3][3];
   TProfile* TpcflowvsEta[3][3];
   TProfile* McflowvspT[3][3];
   TProfile* McflowvsEta[3][3];
 //  TProfile* simflowvspT[4];
 //  TProfile* simflowvsEta[4];

   TProfile* TpcflowvspTEta[3][2];
   TProfile* TpcflowvsEtapT[3][3];
   TProfile* McflowvspTEta[3][2];
   TProfile* McflowvsEtapT[3][3];
 //  TProfile* simflowvspTEta[3];
 //  TProfile* simflowvsEtapT[4];



   TH1F *Phit;

};
#endif // QUALITY_ASSURANCE_SRC_TREE_READER_H_
