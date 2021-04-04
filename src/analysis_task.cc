//
// Created by mikhail on 6/16/20.
//

#include "analysis_task.h"

void AnalysisTask::Init(std::map<std::string, void *> &branch_map) {
  // linking pointers with branch fields
  event_header_ = static_cast<AnalysisTree::EventHeader*>(branch_map.at("McEvent."));
  vtx_tracks_ = static_cast<AnalysisTree::TrackDetector *>(branch_map.at("TpcTracks."));
  fhcal_modules_ = static_cast<AnalysisTree::ModuleDetector *>(branch_map.at("FHCalModules."));
  fhcal_modules_positions_ = data_header_->GetModulePositions(0);
 



  pT_distribution_ = new TH1F( "pT_distribution", ";p_{T} [GeV/c];entries", 250, 0., 2.5 );
  fhcal_energy_distribution_ = new TH1F( "fhcal_energy_distribution", ";E [GeV];entries", 500, 0., 1.0 );
  fhcal_modules_xy_ = new TH2F( "fhcal_modules_xy", ";X;Y", 100, -50., 50.0, 100, -50., 50.0 );
  pT_vs_eta =new TH2F("pT_vs_eta",";eta;p_{T} [GeV/c]",100,-3.0,3.0,100,0,2.5);
  pT_vs_phi =new TH2F("pT_vs_phi",";phi [rad];p_{T} [GeV/c]",100,-7.0,7.0,100,0,2.5);
  phi_vs_eta =new TH2F("phi_vs_eta",";eta;phi [rad]",100,-2.5,2.5,100,-7.0,7.0);
  Energy_vs_moduleId =new TH2F("Energy_vs_moduleId",";moduleId;Energy [GeV]",100,-50,100,100,0,1.0);
  Qxall =new TH1F ("Qxall", "Distribution of Qx and Qy in FHCal ;value of Qx,Qy,GeV;number of pulses N", 90,-3,3);
  Qyall =new TH1F ("Qyall", "Distribution of Qy in FHCal;value of Qy,GeV;number of pulses N", 90,-3,3);
  Qx44 =new TH1F ("Qx44","Distribution of Qx and Qy in first part of FHCal;value of Qx,Qy,Gev;number of pulses N",90,-3,3);
  Qx90 =new TH1F ("Qx90","Distribution of Qx and Qy in second part of FHCal;value of Qx,Qy,GeV;number of pulses N",90,-3,3);
  Qy44 =new TH1F ("Qy44","Distribution of Qy in first part of FHCal;value of Qy,Gev;number of pulses N",90,-3,3);
  Qy90 =new TH1F ("Qy90", "Distribution of Qy in second part of FHCal;value of Qy,Gev;number of pulses N", 90,-3,3);
  ntracksxb =new TH2F ("ntracksxb","multiplicity versus impact parameter,impact parameter b,fermi;multiplicity,number",90,0,18,90,0,500);
  fnall =new TH1F ("fnall","Distribution of event plane angle in FHCal;event plane angle fn,rad;number of pulses N",90,-4,4);
  fn44 =new TH1F ("fn44","Distribution of event plane angle in first part of FHCal;event plane angle fn,rad;number of pulses N",90,-4,4);
  fn90 =new TH1F ("fn90","Distribution of event plane angle in second part of FHCal;event plane angle fn,rad;number of pulses N",90,-4,4);
    for(int j=1;j<9;j++)
    {
    Qx44all[j]=new TH1F(Form("Qx44all_%i",j),Form("Distribution of Qx in first part of FHCal with centrality %i0%;value of Qx,Gev;number of pulses N",j),90,-4,4);

    }
    for(int j=1;j<9;j++)
    {

    Qy44all[j]=new TH1F(Form("Qy44all_%i",j),Form("Distribution of Qy first part of FHCal with centrality %i0%;value of Qy,Gev;number of pulses N",j),90,-4,4);
    }
    for(int j=1;j<9;j++)
    {
    Qx90all[j]=new TH1F(Form("Qx90all_%i",j),Form("Distribution of Qx in second part of FHCal with centrality %i0%;value of Qx,GeV;number of pulses N",j),90,-4,4);
    }
    for(int j=1;j<9;j++)
    {
    Qy90all[j]=new TH1F(Form("Qy90all_%i",j),Form("Distribution of Qy in second part of FHCal with centrality %i0%;value of Qy,Gev;number of pulses N",j),90,-4,4);
    }


}

void AnalysisTask::Exec() {
	auto tracks1=0;
	auto Qx1=0;
	auto Qy1=0;
	auto Qx441=0;
	auto Qy441=0;
	auto Qx901=0;
	auto Qy901=0;
        float ResPhin[8]={0,0,0,0,0,0,0,0};
	float ResN[8]={0,0,0,0,0,0,0,0};
	auto mc_event=event_header_->GetId();
	auto stats=.GetChannel(mc_event);

	
  for( auto& track : *vtx_tracks_->GetChannels() ){
    auto mom3 = track.GetMomentum3();
    auto pT = mom3.Pt();
    auto Eta =track.GetEta();
    auto Phi =track.GetPhi();
    pT_distribution_->Fill(pT);
    pT_vs_eta->Fill(Eta,pT);
    pT_vs_phi->Fill(Phi,pT);
    phi_vs_eta->Fill(Eta,Phi);
    if(abs(track.GetEta())<0.5)
      {
      tracks1=tracks1+1;
      }
      
      ntracksxb->Fill(stats.B(),tracks1);

  }
  for( auto& module : *fhcal_modules_->GetChannels()){
    auto id = module.GetId();
    auto module_pos = fhcal_modules_positions_.GetChannel(id);
    auto signal = module.GetSignal();
    fhcal_energy_distribution_->Fill(signal);
    fhcal_modules_xy_->Fill( module_pos.GetX(), module_pos.GetY() );
    Energy_vs_moduleId->Fill(id, module.GetSignal());
    int w_sign = (id<45) ? 1: -1;
      Qx1=Qx1+ w_sign * signal * cos(module_pos.GetPhi());
      Qy1=Qy1+ w_sign * signal * sin(module_pos.GetPhi());
      if(id<45)
      {
      Qx441=Qx441 + w_sign * signal * cos(module_pos.GetPhi());
      Qy441=Qy441 + w_sign * signal * sin(module_pos.GetPhi());
      }
      else
      {
        Qx901=Qx901 + w_sign * signal * cos(module_pos.GetPhi());
        Qy901=Qy901 + w_sign * signal * sin(module_pos.GetPhi());
      }


  }
   for(int j=1;j<9;j++)
        {
        if((j-1)*10<= ((stats.B())*(stats.B())/4) && ((stats.B())*(stats.B())/4)<j*10)
	{

        Qx44all[j]->Fill(Qx441);
        Qy44all[j]->Fill(Qy441);
        Qx90all[j]->Fill(Qx901);
        Qy90all[j]->Fill(Qy901);
//	fn44all[j]->Fill(atan2(Qy441,Qx441));
//        fn90all[j]->Fill(atan2(Qy901,Qx901));
        ResPhin[j-1]=ResPhin[j-1]+cos(atan2(Qy901,Qx901)-atan2(Qy441,Qx441));
        ResN[j-1]=ResN[j-1]+1;
//      flowvsCen->Fill((j*10-5),(cos(abs(phi_mpd[i]) - atan2(Qy901 + Qy441 - RecQy44[j] - RecQy90[j],Qx901 + Qx441 - RecQx44[j] -RecQx90[j]))/Resolution[j]));
        }
	}
      Qxall->Fill(Qx1);
      Qyall->Fill(Qy1);
      Qx44->Fill(Qx441);
      Qx90->Fill(Qx901);
      Qy44->Fill(Qy441);
      Qy90->Fill(Qy901);
      fnall->Fill(atan2(Qy1,Qx1));
      fn44->Fill(atan2(Qy441,Qx441));
      fn90->Fill(atan2(Qy901,Qx901));



}

void AnalysisTask::Finish() {
  // Writing histograms to file
  out_file_->cd();
  pT_distribution_->Write();
  fhcal_energy_distribution_->Write();
  fhcal_modules_xy_->Write();
  pT_vs_eta->Write();
  pT_vs_phi->Write();
  phi_vs_eta->Write();
  Energy_vs_moduleId->Write();
  Qxall->Write();
  Qyall->Write();
  Qx44->Write();
  Qx90->Write();
  Qy44->Write();
  Qy90->Write();
  ntracksxb->Write();
  fnall->Write();
  fn44->Write();
  fn90->Write();
    for(int j=1;j<9;j++)
    {
    Qx44all[j]->Write();
    }
    for(int j=1;j<9;j++)
    {

    Qy44all[j]->Write();
    }
    for(int j=1;j<9;j++)
    {
    Qx90all[j]->Write();
    }
    for(int j=1;j<9;j++)
    {
    Qy90all[j]->Write();
    }


}
