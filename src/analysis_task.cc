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


  float b_edges[4]={0.,9.,13.,16.};
  energyvsb =new TH2F("Energy in FHCal vs impact parameter",";B [Fm];Energy [GeV]",100,0,16,100,0,30);
  Bt =new TH1F("Impact parameter",";B [Fm];entries",100,0,16);
  Btt =new TH1F("Impact parameter(3)",";B [Fm];entries",3,b_edges);
  pT_distribution_ = new TH1F( "pT_distribution", ";p_{T} [GeV/c];entries", 250, 0., 2.5 );
  fhcal_energy_distribution_ = new TH1F( "fhcal_energy_distribution", ";E [GeV];entries", 500, 0, 1.0 );
  fhcal_modules_xy_ = new TH2F( "fhcal_modules_xy", ";X;Y", 100, -50., 50.0, 100, -50., 50.0 );
  pT_vs_eta =new TH2F("pT_vs_eta","pT_vs_eta;eta;p_{T} [GeV/c]",100,-3.0,3.0,100,0,2.5);
  pT_vs_phi =new TH2F("pT_vs_phi","pT_vs_phi;phi [rad];p_{T} [GeV/c]",100,-4.0,4.0,100,0,2.5);
  phi_vs_eta =new TH2F("phi_vs_eta","phi_vs_eta;eta;phi [rad]",100,-2.5,2.5,100,-4.0,4.0);
  Energy_vs_moduleId =new TH2F("Energy_vs_moduleId",";moduleId;Energy [GeV]",100,0,90,100,0,2.5);
  Qxall =new TH1F ("Qxall", "Distribution of Qx in FHCal ;value of Qx,GeV;entries", 100,-3.0,3.0);
  Qyall =new TH1F ("Qyall", "Distribution of Qy in FHCal;value of Qy,GeV;entries", 100,-3.0,3.0);
  Qx44 =new TH1F ("Qx44","Distribution of Qx in first part of FHCal;value of Qx,Gev;entries",100,-3.0,3.0);
  Qx90 =new TH1F ("Qx90","Distribution of Qx in second part of FHCal;value of Qx,GeV;entries",100,-3.0,3.0);
  Qy44 =new TH1F ("Qy44","Distribution of Qy in first part of FHCal;value of Qy,Gev;entries",100,-3.0,3.0);
  Qy90 =new TH1F ("Qy90", "Distribution of Qy in second part of FHCal;value of Qy,Gev;entries", 100,-3.0,3.0);
  ntracksxb =new TH2F ("ntracksxb","multiplicity versus impact parameter;impact parameter b,fermi;multiplicity,number",100,0,16,90,0,500);
  fnall =new TH1F ("fnall","Distribution of event plane angle in FHCal;event plane angle fn,rad;entries",100,-4.0,4.0);
  fn44 =new TH1F ("fn44","Distribution of event plane angle in first part of FHCal;event plane angle fn,rad;entries",100,-4.0,4.0);
  fn90 =new TH1F ("fn90","Distribution of event plane angle in second part of FHCal;event plane angle fn,rad;entries",100,-4.0,4.0);
    for(int j=1;j<4;j++)
    {
    Qx44all[j]=new TH1F(Form("Qx44all_%i",j),Form("Distribution of Qx in first part of FHCal in part %i;value of Qx,Gev;entries",j),100,-4.0,4.0);

    }
    for(int j=1;j<4;j++)
    {

    Qy44all[j]=new TH1F(Form("Qy44all_%i",j),Form("Distribution of Qy first part of FHCal in part %i;value of Qy,Gev;entries",j),100,-4.0,4.0);
    }
    for(int j=1;j<4;j++)
    {
    Qx90all[j]=new TH1F(Form("Qx90all_%i",j),Form("Distribution of Qx in second part of FHCal in part %i;value of Qx,GeV;entries",j),100,-4.0,4.0);
    }
    for(int j=1;j<4;j++)
    {
    Qy90all[j]=new TH1F(Form("Qy90all_%i",j),Form("Distribution of Qy in second part of FHCal in part %i;value of Qy,Gev;entries",j),100,-4.0,4.0);
    }

    for(int j=1;j<4;j++)
    {
    fn44all[j]=new TH1F(Form("fn44all_%i",j),Form("Distribution of event plane angle in first part of FHCal in part %i;Event plane angle fn[rad];entries",j),90,-4,4);
    }
    for(int j=1;j<4;j++)
    {
    fn90all[j]=new TH1F(Form("fn90all_%i",j),Form("Distribution of event plane angle in second part of FHCal in part %i;Event plane angle fn[rad];entries",j),90,-4,4);
    }

    gr =new TH2F("gr","resolution for half detector vs impact parameter;B[Fm];Res",3,0,16,100,0,1);
    gr2 =new TH2F("gr2","resolution for all detector vs impact parameter;B[Fm];Res",3,0,16,100,0,1);
     flowvsB = new TProfile("flowvsB", "directed flow vs impact parameter;B[Fm];v1",100,0,16);
     flowvsB2 = new TProfile("flowvsB2", "directed flow vs impact parameter;B[Fm];v1",10,0,16);
     flowvspT = new TProfile("flowvspT","directed flow vs pT;pT[GeV/c^2];v1",10,0,2.5);
     flowvsPhi =new TProfile("flowvsPhi","directed flow vs Phi;Phi[rad];v1",10,-4,4);
     flowvsEta =new TProfile("flowvsEta","directed flow vs Eta;pseudorapidity;v1",10,-2.5,2.5);
     ResRP =new TProfile("ResRP","Resolution(RP) vs impact parameter;B[Fm];Res(RP)",10,0,16);
     PhiModule =new TH1F("PhiModule","Phi Distribution of Module;Phi[rad]",90,-4,4);

    


}

void AnalysisTask::Exec() {
	float en=0;
	auto tracks1=0;
	float Qx1=0;/*-0.054;*/
	float Qy1=0;/*-0.0082;*/
	float Qx441=0;/*-0.026;*/
	float Qy441=0;/*-0.005;*/
	float Qx901=0;/*-0.029;*/
	float Qy901=0;/*-0.003;*/
	float Recx44 [4] ={0,-0.012,-0.027,-0.035};
	float Recy44 [4] ={0,0.005,-0.008,-0.01};
	float Recx90 [4] ={0, -0.02,-0.02,-0.037};
	float Recy90 [4] ={0,-0.01,0.0015,-0.0036};
	auto hit = event_header_->GetChannel(0);
        auto B = hit.GetField<float>(0);
	auto PhiRp =hit.GetField<float>(1);
	Bt->Fill(B);
	Btt->Fill(B);
	float Res[3]={0.805,0.776,0.249};

	
  for( auto& track : *vtx_tracks_->GetChannels() ){
    auto mom3 = track.GetMomentum3();
    auto pT = mom3.Pt();
    auto Eta =track.GetEta();
    auto Phi =track.GetPhi();
    pT_distribution_->Fill(pT);
    pT_vs_eta->Fill(Eta,pT);
    pT_vs_phi->Fill(Phi,pT);
    phi_vs_eta->Fill(Eta,Phi);
    if(abs(track.GetEta())<0.5 && pT<2)
      {
      tracks1=tracks1+1;
      }
  }
  ntracksxb->Fill(B,tracks1);
  for( auto& module : *fhcal_modules_->GetChannels()){
    auto id = module.GetId();
    auto module_pos = fhcal_modules_positions_.GetChannel(id);
    auto signal = module.GetSignal();
    fhcal_energy_distribution_->Fill(signal);
    fhcal_modules_xy_->Fill( module_pos.GetX(), module_pos.GetY() );
    Energy_vs_moduleId->Fill(id, module.GetSignal());
    en=en+signal;
    PhiModule->Fill(module_pos.GetPhi());
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
    float QyR1=Qy1 - Recy90[1]-Recy44[1];
    float QxR1=Qx1 - Recx44[1] -Recx90[1];
    float QyR2=Qy1 - Recy90[2]-Recy44[2];
    float QxR2=Qx1 - Recx44[2] -Recx90[2];
    float QyR3=Qy1 - Recy90[3]-Recy44[3];
    float QxR3=Qx1 - Recx44[3] -Recx90[3];

  for( auto& track : *vtx_tracks_->GetChannels() ){
    auto mom3 = track.GetMomentum3();
    auto pT = mom3.Pt();
    auto Eta =track.GetEta();
    auto Phi =track.GetPhi();
    if(abs(Phi)<4 && Eta<=0)
    {
  if(B<9 && B>0)
{
      flowvsB2->Fill(B,(cos(Phi - atan2(QyR1,QxR1))));
      flowvsB->Fill(4.5,(cos(Phi - atan2(QyR1,QxR1))));
      flowvspT->Fill(pT,(cos(Phi - atan2(QyR1,QxR1))));
      flowvsPhi->Fill(Phi,(cos(Phi - atan2(QyR1,QxR1))));
      flowvsEta->Fill(Eta,(cos(Phi - atan2(QyR1,QxR1))));

}
if(B>=9 && B<13)
{
        flowvsB2->Fill(B,(cos(Phi - atan2(QyR2,QxR2))));
        flowvsB->Fill(11.0,(cos(Phi - atan2(QyR2,QxR2))));
	flowvspT->Fill(pT,(cos(Phi - atan2(QyR2,QxR2))));
      flowvsPhi->Fill(Phi,(cos(Phi - atan2(QyR2,QxR2))));
      flowvsEta->Fill(Eta,(cos(Phi - atan2(QyR2,QxR2))));

}
if(B>=13 && B<16)
{
        flowvsB2->Fill(B,(cos(Phi - atan2(QyR3,QxR3))));
        flowvsB->Fill(14.5,(cos(Phi - atan2(QyR3,QxR3))));
	flowvspT->Fill(pT,(cos(Phi - atan2(QyR3,QxR3))));
      flowvsPhi->Fill(Phi,(cos(Phi - atan2(QyR3,QxR3))));
      flowvsEta->Fill(Eta,(cos(Phi - atan2(QyR3,QxR3))));

}
}
}

        if( B<9)
	{

        Qx44all[1]->Fill(Qx441);
        Qy44all[1]->Fill(Qy441);
        Qx90all[1]->Fill(Qx901);
        Qy90all[1]->Fill(Qy901);
	fn44all[1]->Fill(atan2(Qy441,Qx441));
        fn90all[1]->Fill(atan2(Qy901,Qx901));
        ResPhin[0]=ResPhin[0]+cos(atan2(Qy901-Recy90[1],Qx901-Recx90[1])-atan2(Qy441-Recy44[1],Qx441-Recx44[1]));
        ResN[0]=ResN[0]+1;
	ResRP->Fill(4.5,cos(atan2(QyR1,QxR1)-PhiRp));
//      flowvsB->Fill(4,(cos(module_pos.GetPhi() - atan2(Qy901 + Qy441 - Recy44[1] - Recy90[1],Qx901 + Qx441 - Recx44[1] -Recx90[1]))));
        }
	if( B>=9 && B<13)
        {

        Qx44all[2]->Fill(Qx441);
        Qy44all[2]->Fill(Qy441);
        Qx90all[2]->Fill(Qx901);
        Qy90all[2]->Fill(Qy901);
      fn44all[2]->Fill(atan2(Qy441,Qx441));
        fn90all[2]->Fill(atan2(Qy901,Qx901));
        ResPhin[1]=ResPhin[1]+cos(atan2(Qy901-Recy90[2],Qx901-Recx90[2])-atan2(Qy441-Recy44[2],Qx441-Recx44[2]));
        ResN[1]=ResN[1]+1;
	ResRP->Fill(11,cos(atan2(QyR2,QxR2)-PhiRp));

  //    flowvsB->Fill(11.5,(cos(module_pos.GetPhi() - atan2(Qy901 + Qy441 - Recy44[2] - Recy90[2],Qx901 + Qx441 - Recx44[2] -Recx90[2]))));
        }
if(B>=13 && B<16)
        {

        Qx44all[3]->Fill(Qx441);
        Qy44all[3]->Fill(Qy441);
        Qx90all[3]->Fill(Qx901);
        Qy90all[3]->Fill(Qy901);
        fn44all[3]->Fill(atan2(Qy441,Qx441));
   /*     fn90all[3]->Fill(atan2(Qy901,Qx901));*/
        ResPhin[2]=ResPhin[2]+cos(atan2(Qy901-Recy90[3],Qx901-Recx90[3])-atan2(Qy441-Recy44[3],Qx441-Recx44[3]));
        ResN[2]=ResN[2]+1;
	ResRP->Fill(14.5,cos(atan2(QyR3,QxR3)-PhiRp));
    //  flowvsB->Fill(14.5,(cos(module_pos.GetPhi() - atan2(Qy901 + Qy441 - Recy44[3] - Recy90[3],Qx901 + Qx441 - Recx44[3] -Recx90[3]))));
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
      energyvsb->Fill(B,en);
      



}

void AnalysisTask::Finish() {
  // Writing histograms to file
  out_file_->cd();

  PhiModule->Write();
  ResRP->Write();
  flowvspT->Write();
  flowvsPhi->Write();
  flowvsEta->Write();

  flowvsB->Write();
  flowvsB2->Write();
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
    for(int j=1;j<4;j++)
    {
    Qx44all[j]->Write();
    }
    for(int j=1;j<4;j++)
    {

    Qy44all[j]->Write();
    }
    for(int j=1;j<4;j++)
    {
    Qx90all[j]->Write();
    }
    for(int j=1;j<4;j++)
    {
    Qy90all[j]->Write();
    }
    Bt->Write();
    Btt->Write();
    energyvsb->Write();
    for(int j=1;j<4;j++)
    {
    fn90all[j]->Write();
    }
for(int j=1;j<4;j++)
    {
    fn44all[j]->Write();
    }

     
for(int j=0;j<3;j++)
{
	std::cout<<ResPhin[j]<<" /// "<< ResN[j]<<"next";
}





}
