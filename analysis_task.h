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
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TProfile.h>
#include <memory>
#include <string>

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
  AnalysisTree::ModuleDetector* fhcal_modules_{nullptr}; 		// modules of FhCal branch
  AnalysisTree::ModulePositions fhcal_modules_positions_;
  TH1F* pT_distribution_;
  TH1F* fhcal_energy_distribution_;
  TH2F*fhcal_modules_xy_;
  TH2F* pT_vs_eta;
  TH2F* pT_vs_phi;
  TH2F* phi_vs_eta;
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
  TH1F *Qx44all[9];
  TH1F *Qy44all[9];
  TH1F *Qx90all[9];
  TH1F *Qy90all[9];
  TH1F *Bt;

};
#endif // QUALITY_ASSURANCE_SRC_TREE_READER_H_