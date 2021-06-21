//
// Created by mikhail on 6/16/20.
//

#include "analysis_task.h"

void AnalysisTask::Init(std::map<std::string, void *> &branch_map) {
  // linking pointers with branch fields
  event_header_ = static_cast<AnalysisTree::EventHeader*>(branch_map.at("McEvent."));
  vtx_tracks_ = static_cast<AnalysisTree::TrackDetector *>(branch_map.at("TpcTracks."));
  fhcal_modules_ = static_cast<AnalysisTree::ModuleDetector *>(branch_map.at("FHCalModules."));
  mc_tracks_ = static_cast<AnalysisTree::Particles *>(branch_map.at("McTracks."));
  fhcal_modules_positions_ = data_header_->GetModulePositions(0);
  tpc_mc_matching = static_cast<AnalysisTree::Matching *>(branch_map.at("TpcTracks2McTracks"));


  float b_edges[17]={0.,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16};
  float b_edges3[4]={0.,4.,8.,11.};
  std::string part[3]={"proton","kaon","pion"};
  std::string partpT[3][3]={{"0.2<pT<0.5 GeV/c (proton)","0.5<pT<1 GeV/c (proton)","1<pT<2.5 GeV/c (proton)"},{"0.2<pT<0.5 GeV/c (positive kaon)","0.5<pT<1 GeV/c (positive kaon)","1<pT<2.5 GeV/c (positive kaon)"},{"0.2<pT<0.5 GeV/c (positive pion)","0.5<pT<1 GeV/c (positive pion)","1<pT<2.5 GeV/c (positive pion)"}};
  std::string partEta[3][2]={{"0<|Eta|<0.5 (proton)","0.5<|Eta|<1 (proton)"},{"0<|Eta|<0.5 (positive kaon)","0.5<|Eta|<1 (positive kaon)"},{"0<|Eta|<0.5 (positive pion)","0.5<|Eta|<1 (positive pion)"}};
  std::string b_edges3part[3][3]={{"central(proton)","semi-central(proton)","peripheral(proton)"},{"central(positive kaon)","semi-central(positive kaon)","peripheral(positive kaon)"},{"central(positive pion)","semi-central(positive pion)","peripheral(positive pion)"}};
  energyvsb =new TH2F("Energy in FHCal vs impact parameter",";B [Fm];Energy [GeV]",100,0,16,100,0,30);
  Bt =new TH1F("Impact parameter",";B [Fm];entries",100,0,16);
  Btt =new TH1F("Impact parameter(8)",";B [Fm];entries",8,b_edges);
  pT_distribution_ = new TH1F( "pT_distribution", ";p_{T} [GeV/c];entries", 250, 0., 2.5 );
  fhcal_energy_distribution_ = new TH1F( "fhcal_energy_distribution", ";E [GeV];entries", 500, 0, 1.0 );
  fhcal_modules_xy_ = new TH2F( "fhcal_modules_xy", ";X;Y", 100, -50., 50.0, 100, -50., 50.0 );
  McpT_vs_eta_proton =new TH2F("McpT_vs_eta_proton","pT_vs_eta_proton(Mc);eta;p_{T} [GeV/c]",100,-3.0,3.0,100,0,2.5);
  McpT_vs_eta_kaon =new TH2F("McpT_vs_eta_kaon","pT_vs_eta_kaon(Mc);eta;p_{T} [GeV/c]",100,-3.0,3.0,100,0,2.5);
  McpT_vs_eta_pion =new TH2F("McpT_vs_eta_pion","pT_vs_eta_pion(Mc);eta;p_{T} [GeV/c]",100,-3.0,3.0,100,0,2.5);
McpT_vs_eta_all =new TH2F("McpT_vs_eta_all","pT_vs_eta_all(Mc);eta;p_{T} [GeV/c]",100,-3.0,3.0,100,0,2.5);

  McpT_vs_phi =new TH2F("McpT_vs_phi","pT_vs_phi(Mc);phi [rad];p_{T} [GeV/c]",100,-4.0,4.0,100,0,2.5);

  Mcphi_vs_eta =new TH2F("Mcphi_vs_eta","phi_vs_eta(Mc);eta;phi [rad]",100,-2.5,2.5,100,-4.0,4.0);
  TpcpT_vs_eta_proton =new TH2F("TpcpT_vs_eta_proton","pT_vs_eta_proton(Tpc);eta;p_{T} [GeV/c]",100,-3,3.0,100,0,2.5);
  TpcpT_vs_eta_kaon =new TH2F("TpcpT_vs_eta_kaon","pT_vs_eta_kaon(Tpc);eta;p_{T} [GeV/c]",100,-3,3.0,100,0,2.5);
  TpcpT_vs_eta_pion =new TH2F("TpcpT_vs_eta_pion","pT_vs_eta_pion(Tpc);eta;p_{T} [GeV/c]",100,-3,3.0,100,0,2.5);
  TpcpT_vs_eta_all =new TH2F("TpcpT_vs_eta_all","pT_vs_eta_all(Tpc);eta;p_{T} [GeV/c]",100,-3,3.0,100,0,2.5);

  TpcpT_vs_phi =new TH2F("TpcpT_vs_phi","pT_vs_phi(Tpc);phi [rad];p_{T} [GeV/c]",100,-4.0,4.0,100,0,2.5);
  Tpcphi_vs_eta =new TH2F("Tpcphi_vs_eta","phi_vs_eta(Tpc);eta;phi [rad]",100,-2.5,2.5,100,-4.0,4.0);

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
    for(int j=1;j<17;j++)
    {
    Qx44all[j]=new TH1F(Form("Qx44all_%i",j),Form("Distribution of Qx in first part of FHCal in part %i;value of Qx,Gev;entries",j),100,-4.0,4.0);

    }
    for(int j=1;j<17;j++)
    {

    Qy44all[j]=new TH1F(Form("Qy44all_%i",j),Form("Distribution of Qy first part of FHCal in part %i;value of Qy,Gev;entries",j),100,-4.0,4.0);
    }
    for(int j=1;j<17;j++)
    {
    Qx90all[j]=new TH1F(Form("Qx90all_%i",j),Form("Distribution of Qx in second part of FHCal in part %i;value of Qx,GeV;entries",j),100,-4.0,4.0);
    }
    for(int j=1;j<17;j++)
    {
    Qy90all[j]=new TH1F(Form("Qy90all_%i",j),Form("Distribution of Qy in second part of FHCal in part %i;value of Qy,Gev;entries",j),100,-4.0,4.0);
    }

    for(int j=1;j<17;j++)
    {
    fn44all[j]=new TH1F(Form("fn44all_%i",j),Form("Distribution of event plane angle in first part of FHCal in part %i;Event plane angle fn[rad];entries",j),90,-4,4);
    }
    for(int j=1;j<17;j++)
    {
    fn90all[j]=new TH1F(Form("fn90all_%i",j),Form("Distribution of event plane angle in second part of FHCal in part %i;Event plane angle fn[rad];entries",j),90,-4,4);
    }
    for(int j=0;j<3;j++)
    {

    TpcflowvsB[j]= new TProfile(Form("TpcflowvsB_%p",part[j]),Form("Directed flow vs impact parameter %s;B[fm];v1",part[j]),16,b_edges);
    TpcPhivsB[j]= new TProfile(Form("TpcPhivsB_%s",part[j]),Form("TpcPhi vs impact parameter %s;B[fm];Phi",part[j]),32,-16,16);
    McPhivsB[j]= new TProfile(Form("McPhivsB_%s",part[j]),Form("McPhi vs impact parameter %s;B[fm];Phi",part[j]),32,-16,16);
    McflowvsB[j]= new TProfile(Form("McflowvsB_%s",part[j]),Form("Directed flow vs impact parameter %s;B[fm];v1",part[j]),16,b_edges);
    }
//    simflowvsB= new TProfile("simflowvsB","Directed flow vs impact parameter;B[fm];v1",16,b_edges);

     for(int i=0;i<3;i++)
     {
	     for(int j=0;j<3;j++)
	     {

	     TpcflowvspT[i][j]= new TProfile(Form("TpcflowvspT_%s",b_edges3part[i][j]),Form("Directed flow vs pT %s;pT[GeV/c];v1",b_edges3part[i][j]),10,0,2.5);
	     TpcflowvsEta[i][j]= new TProfile(Form("TpcflowvsEta_%s",b_edges3part[i][j]),Form("Directed flow vs pseudorapidity %s;pseudorapidity;v1",b_edges3part[i][j]),10,-1,1);
	     McflowvspT[i][j]= new TProfile(Form("McflowvspT_%s",b_edges3part[i][j]),Form("Directed flow vs pT %s;pT[GeV/c];v1",b_edges3part[i][j]),10,0,2.5);
             McflowvsEta[i][j]= new TProfile(Form("McflowvsEta_%s",b_edges3part[i][j]),Form("Directed flow vs pseudorapidity %s;pseudorapidity;v1",b_edges3part[i][j]),10,-1,1);
	     TpcflowvsEtapT[i][j]= new TProfile(Form("TpcflowvsEtapT_%s",partpT[i][j]),Form("Directed flow vs pseudorapidity %s;pseudorapidity;v1",partpT[i][j]),10,-1,1);
             McflowvsEtapT[i][j]= new TProfile(Form("McflowvsEtapT_%s",partpT[i][j]),Form("Directed flow vs pseudorapidity %s;pseudorapidity;v1",partpT[i][j]),10,-1,1);
  //           simflowvspT[j]=new TProfile(Form("simflowvspT_%i",j),Form("Directed flow vs pT %i;pT[GeV/c];v1",j),10,0,2.5);
//	     simflowvsEta[j]= new TProfile(Form("simflowvsEta_%i",j),Form("Directed flow vs pseudorapidity %i;pseudorapidity;v1",j),10,-1,1);
  //           simflowvsEtapT[j]= new TProfile(Form("simflowvsEtapT_%i",j),Form("Directed flow vs pseudorapidity %i;pseudorapidity;v1",j),1,1);

     }
     }
     for(int i=0;i<3;i++)
     {
	     for(int j=0;j<2;j++)
	     {
	     TpcflowvspTEta[i][j]= new TProfile(Form("TpcflowvspTEta_%s",partEta[i][j]),Form("Directed flow vs pT %s;pT[GeV/c];v1",partEta[i][j]),10,0,2.5);
             McflowvspTEta[i][j]= new TProfile(Form("McflowvspTEta_%s",partEta[i][j]),Form("Directed flow vs pT %s;pT[GeV/c];v1",partEta[i][j]),10,0,2.5);
     }
     }



     ResRPS =new TProfile("ResRPS","Resolution(RP) vs impact parameter;B[Fm];Res(RP)",16,b_edges);
     ResRPN =new TProfile("ResRPN","Resolution(RPN) vs impact parameter;B[Fm];Res(RP)",16,b_edges);

     ResRP3 =new TProfile("ResRP3","Resolution(RP3) vs impact parameter;B[Fm];Res(RP3)",3,b_edges3);
     ResHalf =new TProfile("ResHalf","Resolution(Half) vs impact parameter;B[Fm];ResHalf",16,b_edges);
     PhiModule =new TH1F("PhiModule","Phi Distribution of Module;Phi[rad]",90,-4,4);

     Phit =new TH1F("Phit","Phit",90,-4,4);

    PID_proton =new TH1F("PID_proton",";PI_proton;entries",100,0,1);
    dEdX_vs_pz_proton =new TH2F("dEdX_vs_p/z_proton","dEdX_vs_p/z proton;p/z [GeV/c];dEdX",1000,-5,5,1000,0,15000);
    mass2_vs_pz_proton =new TH2F("mass2_vs_p/z_proton","mass2_vs_p/z proton;p/z [GeV/c];mass2",1000,-5,5,1000,0,5);
    dEdX_vs_pz_pion =new TH2F("dEdX_vs_p/z_pion","dEdX_vs_p/z pion+;p/z [GeV/c];dEdX",1000,-5,5,1000,0,15000);
    mass2_vs_pz_pion =new TH2F("mass2_vs_p/z_pion","mass2_vs_p/z pion+;p/z [GeV/c];mass2",1000,-5,5,1000,0,1);
    dEdX_vs_pz_kaon =new TH2F("dEdX_vs_p/z_kaon","dEdX_vs_p/z kaon+;p/z [GeV/c];dEdX",1000,-5,5,1000,0,15000);
    mass2_vs_pz_kaon =new TH2F("mass2_vs_p/z_kaon","mass2_vs_p/z kaon+;p/z [GeV/c];mass2",1000,-5,5,1000,0,1); 
    dEdX_vs_pz_all =new TH2F("dEdX_vs_p/z_all","dEdX_vs_p/z all;p/z [GeV/c];dEdX",1000,-5,5,1000,0,15000);





}

void AnalysisTask::Exec() {
	float b_edges[17]={0.,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16};
	float b_edges3[4]={0.,4.,8.,11.};
	float en=0;
	auto tracks1=0;
	float Qx1=0;/*-0.054;*/
	float Qy1=0;/*-0.0082;*/
	float Qx441=0;/*-0.026;*/
	float Qy441=0;/*-0.005;*/
	float Qx901=0;/*-0.029;*/
	float Qy901=0;/*-0.003;*/
	float Recx44 [4] ={0.,0.0041,0.0143,0.0253};
	float Recy44 [4] ={0.,0.0016,0.0058,0.0062};
	float Recx90 [4] ={0., 0.0075,0.0133,0.0266};
	float Recy90 [4] ={0.,0.0024,0.0023,0.0011};

/*	float Recx44R [8] ={0.0075,0.0145,0.0206,0.0242,0.0331,0.0347,0.351,0.0364};
	float Recy44R [8] ={0.0025,0.0056,0.0078,0.0069,0.0063,0.0054,0.0033,0.0032};
	float Recx90R [8] ={0.0094,0.0127,0.0214,0.0265,0.0343,0.0353,0.0348,0.035};
	float Recy90R [8] ={0.0020,0.0015,0.0025,0.0019,0.0005,-0.0011,-0.0022,-0.0020}; */
	float Recx44R [16] ={0.0011,0.002,0.0059,0.0044,0.0122,0.0116,0.0145,0.0174,0.0214,0.0216,0.0313,0.0329,0.0340,0.0344,0.0361,0.0364};
	float Recy44R [16] ={-0.0042,0.0028,0.002,0.0019,0.0042,0.0034,0.0054,0.0085,0.0069,0.0076,0.0048,0.0078,0.0056,0.0034,0.0039,0.0032};
	float Recx90R [16] ={0.0004,0.0074,0.0065,0.0095,0.0123,0.0107,0.0126,0.0168,0.0224,0.0256,0.0302,0.0355,0.0346,0.0364,0.0339,0.035};
	float Recy90R [16] ={0.0017,0.0026,0.0040,0.0013,0.0021,0.0032,-0.0005,0.0043,0.0010,0.0032,-0.0003,0.0012,-0.001,-0.001,-0.0032,-0.002};

	auto hit = event_header_->GetChannel(0);
	
	
        auto B = hit.GetField<float>(0);
	auto PhiRp =hit.GetField<float>(1);
	Bt->Fill(B);
	Btt->Fill(B);
	float Res[3]={0.4299,0.8392,0.8629};
	float Res11[11]={0.1345,0.3509,0.55,0.6843,0.7723,0.832,0.8674,0.8851,0.8874,0.8745,0.8269};
	float Resrp[3]={0.8026,0.6808,0.0679};
	float pidcode[3]={2212,321,211};
	float catpT[4]={0.2,0.5,1,2.5};
	float catEta[3]={0,0.5,1};

	int w=0;
  for( auto& track : *vtx_tracks_->GetChannels() )
  {
	  
    auto mom3 = track.GetMomentum3();
    auto pT = mom3.Pt();
    auto Eta =track.GetEta();
    auto Phi =track.GetPhi();
    auto pid_pr = track.GetField<float>(9);
    auto pid_kaon = track.GetField<float>(8);
    auto pid_pion = track.GetField<float>(7);
    auto dedx =track.GetField<float>(6);
    auto p = track.GetField<float>(-7);
    auto charge = track.GetField<int>(2);
    auto mass2 = track.GetField<float>(4);
    PID_proton->Fill(pid_pr);
  /* if(pid_pion>0.95)
   { 
    dEdX_vs_pz->Fill(p/charge,dedx);
   mass2_vs_pz->Fill(p/charge,mass2);} */
    pT_distribution_->Fill(pT);
    if(pT>0.2 && Eta!=0 && pid_pr>0.95 && charge>0){
    TpcpT_vs_eta_proton->Fill(Eta,pT);}
    if(pT>0.2 && Eta!=0 && pid_kaon>0.95 && charge>0){
    TpcpT_vs_eta_kaon->Fill(Eta,pT);}
if(pT>0.2 && Eta!=0 && pid_pion>0.95 && charge>0){
    TpcpT_vs_eta_pion->Fill(Eta,pT);}
if(pT>0.2 && Eta!=0)
{	
TpcpT_vs_eta_all->Fill(Eta,pT);
}

    TpcpT_vs_phi->Fill(Phi,pT);
    Phit->Fill(Phi);
    Tpcphi_vs_eta->Fill(Eta,Phi);
    if(abs(track.GetEta())<0.5 && pT<2)
      {
      tracks1=tracks1+1;
      }

    //matching
   //из бранча с матчингом берётся прямое совпадение и записывается id совпавшего трека
      auto matched_track = tpc_mc_matching->GetMatchDirect(w); // matching simulated track's ID to current reconstructed track

      //из бранча с сим-треками(true-траекториями) берётся соответствующий трек (набор кин. характеристик)
      auto sim_track = mc_tracks_->GetChannel(matched_track);    // getting matched simulated track // changed i->matched_track
      auto sim_pid = sim_track.GetPid();                           // getting PID of matched simulated track

      auto sim_pT = sim_track.GetPt(); // getting transverse momentum of matched simulated track

      auto sim_Eta = sim_track.GetRapidity(); // getting rapidity of matched simulated track
      if(sim_pid==2212 && pT>0.2 && pT<2.5 && abs(Eta)<1)
      {
	      dEdX_vs_pz_proton->Fill(p/charge,dedx);
   mass2_vs_pz_proton->Fill(p/charge,mass2);
      }
      if(sim_pid==321 && pT>0.2 && pT<2.5 && abs(Eta)<1)
      {
	      dEdX_vs_pz_kaon->Fill(p/charge,dedx);
   mass2_vs_pz_kaon->Fill(p/charge,mass2);
      }
      if(sim_pid==211 && pT>0.2 && pT<2.5 && abs(Eta)<1)
      {
	      dEdX_vs_pz_pion->Fill(p/charge,dedx);
   mass2_vs_pz_pion->Fill(p/charge,mass2);
      }
      if(pT>0.2 && pT<2.5 && abs(Eta<1))
      {
      dEdX_vs_pz_all->Fill(p/charge,dedx);
      }


      w=w+1; 
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
  float QxR[4]={0.,Qx1 - Recx44[1] -Recx90[1],Qx1 - Recx44[2] -Recx90[2],Qx1 - Recx44[3] -Recx90[3]};
  float QyR[4]={0.,Qy1 - Recy90[1]-Recy44[1],Qy1 - Recy90[2]-Recy44[2],Qy1 - Recy90[3]-Recy44[3]};
float QxR11[11];
float QyR11[11];
for(int i=0;i<11;i++)
{
	QxR11[i]=Qx1-Recx44R[i] -Recx90R[i];
	QyR11[i]=Qy1-Recy44R[i] -Recy90R[i];
}
w=0;
  for( auto& track : *mc_tracks_->GetChannels() ){
	  auto McpT = track.GetField<float>(-2);
	  auto McEta = track.GetField<float>(-6);
	  auto McPhi = track.GetField<float>(-1);
	  auto Mcpid = track.GetField<int>(-4);
	  if(Mcpid==2212)
	  {
	  McpT_vs_eta_proton->Fill(McEta,McpT);}
	  if(Mcpid==321)
          {
          McpT_vs_eta_kaon->Fill(McEta,McpT);}
if(Mcpid==211)
          {
          McpT_vs_eta_pion->Fill(McEta,McpT);}
McpT_vs_eta_all->Fill(McEta,McpT);

    McpT_vs_phi->Fill(McPhi,McpT);
    Mcphi_vs_eta->Fill(McEta,McPhi);
for(int k=0;k<3;k++)
{
	if(Mcpid==pidcode[k])
	{

	  for(int i=0;i<3;i++)
	  {

		  if(abs(McPhi)<4 && B>b_edges3[i] && B<=b_edges3[i+1])
		  {

			  if(McEta>0 && McEta<1)
			  {
				  if(McpT>0.2 && McpT<2.5){
					  McflowvsB[k]->Fill(B,cos(McPhi-PhiRp));
				  McPhivsB[k]->Fill(B,McPhi);
				  }
				  McflowvspT[k][i]->Fill(McpT,cos(McPhi-PhiRp));
			  }
			  if(McEta<0 && McEta>-1)
			  {
				  if(McpT>0.2 && McpT<2.5){
					  McflowvsB[k]->Fill(B,-cos(McPhi-PhiRp));
				  McPhivsB[k]->Fill(-B,McPhi);
				  }
                                  McflowvspT[k][i]->Fill(McpT,-cos(McPhi-PhiRp));
			  }
			 if(McpT>0.2 && McpT<2.5)
                          {McflowvsEta[k][i]->Fill(McEta,cos(McPhi-PhiRp));}
		  }
		  if(B>0 && B<=11 && McpT>catpT[i] && McpT<=catpT[i+1])
		  {
	McflowvsEtapT[k][i]->Fill(McEta,cos(McPhi-PhiRp));
}
}
for(int i=0;i<2;i++)
{
	if(abs(McEta)>catEta[i] && abs(McEta)<=catEta[i+1] && B>0 && B<=11)
	{

        if(McEta>0)
        McflowvspTEta[k][i]->Fill(McpT,cos(McPhi-PhiRp));
        else
                McflowvspTEta[k][i]->Fill(McpT,-cos(McPhi-PhiRp));
}
}
}
}






}
w=0;
 for( auto& track : *vtx_tracks_->GetChannels() ){
    auto mom3 = track.GetMomentum3();
    auto TpcpT = mom3.Pt();
    auto TpcEta =track.GetEta();
    auto TpcPhi =track.GetPhi();
    float pid_pr = track.GetField<float>(9);
    float pid_kaon = track.GetField<float>(8);
    float pid_pion = track.GetField<float>(7);
    float probpid[3]={pid_pr,pid_kaon,pid_pion};
    auto charge = track.GetField<int>(2);
    auto matched_track = tpc_mc_matching->GetMatchDirect(w); // matching simulated track's ID to current reconstructed track

      //из бранча с сим-треками(true-траекториями) берётся соответствующий трек (набор кин. характеристик)
      auto sim_track = mc_tracks_->GetChannel(matched_track);    // getting matched simulated track // changed i->matched_track
      auto sim_pid = sim_track.GetPid();                           // getting PID of matched simulated track
      auto sim_Phi =sim_track.GetPhi();
     auto sim_pT = sim_track.GetPt(); // getting transverse momentum of matched simulated track

      auto sim_Eta = sim_track.GetRapidity(); // getting rapidity of matched simulated track

          for(int k=0;k<3;k++)
          {
		  if( probpid[k]>0.95 && charge>0)
		  {
			  for(int i=0;i<3;i++)
			  {
                  if(abs(TpcPhi)<4 && B>b_edges3[i] && B<=b_edges3[i+1])
                  {
                          if(TpcEta>0.2 && TpcEta<1)
                          {
                                  if(TpcpT>0.2 && TpcpT<2.5)
				  {
				  TpcPhivsB[k]->Fill(B,TpcPhi);
				  }
                                  TpcflowvspT[k][i]->Fill(TpcpT,(cos(TpcPhi-atan2(QyR[i],QxR[i])))/Res[i]);
                          }
                          if(TpcEta<-0.2 && TpcEta>-1)
                          {
                                  if(TpcpT>0.2 && TpcpT<2.5){
				  TpcPhivsB[k]->Fill(-B,TpcPhi);
				  }
                                  TpcflowvspT[k][i]->Fill(TpcpT,(-cos(TpcPhi-atan2(QyR[i],QxR[i])))/Res[i]);
                          }

                          if(TpcpT>0.2 && TpcpT<2.5)
                          {TpcflowvsEta[k][i]->Fill(TpcEta,(cos(TpcPhi-atan2(QyR[i],QxR[i])))/Res[i]);}

		  }
		  if(B>0 && B<=11 && TpcpT>catpT[i] && TpcpT<=catpT[i+1])
		 {
        TpcflowvsEtapT[k][i]->Fill(TpcEta,cos(TpcPhi-atan2(QyR[i],QxR[i]))/Res[i]);
}
}
for(int i=0;i<2;i++)
{
	if(abs(TpcEta)>catEta[i] && abs(TpcEta)<=catEta[i+1] && B>0 && B<=11)
	{
        if(TpcEta>0)
        TpcflowvspTEta[k][i]->Fill(TpcpT,cos(TpcPhi-atan2(QyR[i],QxR[i]))/Res[i]);
        else
                TpcflowvspTEta[k][i]->Fill(TpcpT,-cos(TpcPhi-atan2(QyR[i],QxR[i]))/Res[i]);
}


                  }

/*                  if( B>b_edges3[i-1] && B<=b_edges3[i] && sim_pid==211)
                  {
                          if(sim_Eta>0.2 && sim_Eta<1)
                          {
                                  if(sim_pT>0.2 && sim_pT<2.5){
                                  simflowvsB->Fill(B,(cos(sim_Phi-atan2(QyR[i],QxR[i])))/Res[i]);}
                                  simflowvspT[i]->Fill(sim_pT,(cos(sim_Phi-atan2(QyR[i],QxR[i])))/Res[i]);
                          }
                          if(sim_Eta<-0.2 && sim_Eta>-1)
                          {
                                  if(sim_pT>0.2 && sim_pT<2.5){
                                          simflowvsB->Fill(B,(-cos(sim_Phi-atan2(QyR[i],QxR[i])))/Res[i]);}
                                  simflowvspT[i]->Fill(sim_pT,(-cos(sim_Phi-atan2(QyR[i],QxR[i])))/Res[i]);
                          }
                          if(sim_pT>0.2 && sim_pT<2.5)
                          {simflowvsEta[i]->Fill(sim_Eta,(cos(sim_Phi-atan2(QyR[i],QxR[i])))/Res[i]);}


                        if(sim_pT>0.2 && sim_pT<0.5)
{
        simflowvsEtapT[1]->Fill(sim_Eta,cos(sim_Phi-atan2(QyR[i],QxR[i]))/Res[i]);
}
if(sim_pT>=0.5 && sim_pT<1)
{
        simflowvsEtapT[2]->Fill(sim_Eta,cos(sim_Phi-atan2(QyR[i],QxR[i]))/Res[i]);
}
if(sim_pT>=1 && sim_pT<2.5)
{
        simflowvsEtapT[3]->Fill(sim_Eta,cos(sim_Phi-atan2(QyR[i],QxR[i]))/Res[i]);
}

if(abs(sim_Eta>0) && abs(sim_Eta<0.5))
                {
        if(sim_Eta>0)
        simflowvspTEta[1]->Fill(sim_pT,cos(sim_Phi-atan2(QyR[i],QxR[i]))/Res[i]);
        else
                simflowvspTEta[1]->Fill(sim_pT,-cos(sim_Phi-atan2(QyR[i],QxR[i]))/Res[i]);
}
if(abs(sim_Eta>=0.5) && abs(sim_Eta<1))
                {
        if(sim_Eta>0)
        simflowvspTEta[2]->Fill(sim_pT,cos(sim_Phi-atan2(QyR[i],QxR[i]))/Res[i]);
        else
                simflowvspTEta[2]->Fill(sim_pT,-cos(sim_Phi-atan2(QyR[i],QxR[i]))/Res[i]);
} 


                  
} */

for(int i=0;i<11;i++)
          {
                  if(abs(TpcPhi)<4 && B>b_edges[i] && B<=b_edges[i+1])
                  {
                          if(TpcEta>0.2 && TpcEta<1)
                          {
                                  if(TpcpT>0.2 && TpcpT<2.5){
                                  TpcflowvsB[k]->Fill(B,(cos(TpcPhi-PhiRp)));}
                          }
                          if(TpcEta<-0.2 && TpcEta>-1)
                          {
                                  if(TpcpT>0.2 && TpcpT<2.5){
                                          TpcflowvsB[k]->Fill(B,(-cos(TpcPhi-PhiRp)));}
                                  }
                          }
		  }
}
}
	  
w=w+1;

}

for(int j=1;j<17;j++)
{
	if(B<b_edges[j] && B>b_edges[j-1])
	{
    /*    Qx44all[j]->Fill(Qx441);
        Qy44all[j]->Fill(Qy441);
        Qx90all[j]->Fill(Qx901);
        Qy90all[j]->Fill(Qy901);
	fn44all[j]->Fill(atan2(Qy441,Qx441));
        fn90all[j]->Fill(atan2(Qy901,Qx901)); */
       /* ResPhin[j-1]=ResPhin[j-1]+cos(atan2(Qy901-Recy90R[j-1],Qx901-Recx90R[j-1])-atan2(Qy441-Recy44R[j-1],Qx441-Recx44R[j-1]));
        ResN[j-1]=ResN[j-1]+1;  */
	ResHalf->Fill((b_edges[j-1]+b_edges[j])/2,cos(atan2(Qy901-Recy90R[j-1],Qx901-Recx90R[j-1])-atan2(Qy441-Recy44R[j-1],Qx441-Recx44R[j-1])));
	ResRPS->Fill((b_edges[j-1]+b_edges[j])/2,cos(atan2(Qy901-Recy90R[j-1]+Qy441-Recy44R[j-1],Qx901-Recx90R[j-1]+Qx441-Recx44R[j-1])-PhiRp));
	ResRPN->Fill((b_edges[j-1]+b_edges[j])/2,cos(atan2(/*Qy901-Recy90R[j-1]+*/Qy441-Recy44R[j-1],/*Qx901-Recx90R[j-1]+*/Qx441-Recx44R[j-1])-PhiRp));
//      flowvsB->Fill(4,(cos(module_pos.GetPhi() - atan2(Qy901 + Qy441 - Recy44[1] - Recy90[1],Qx901 + Qx441 - Recx44[1] -Recx90[1]))));
        } 


}
if(B>0 && B<=4)
{
	Qx44all[1]->Fill(Qx441);
        Qy44all[1]->Fill(Qy441);
        Qx90all[1]->Fill(Qx901);
        Qy90all[1]->Fill(Qy901);
}
if(B>4 && B<=8)
{
	Qx44all[2]->Fill(Qx441);
        Qy44all[2]->Fill(Qy441);
        Qx90all[2]->Fill(Qx901);
        Qy90all[2]->Fill(Qy901);
}
if(B>8 && B<=11)
{
	Qx44all[3]->Fill(Qx441);
        Qy44all[3]->Fill(Qy441);
        Qx90all[3]->Fill(Qx901);
        Qy90all[3]->Fill(Qy901);
}


/*
if(B<9 && B>0)
{
	ResRP3->Fill(4.5,cos(atan2(QyR1,QxR1)-PhiRp));
}
if(B>=9 && B<13)
{
	ResRP3->Fill(11,cos(atan2(QyR2,QxR2)-PhiRp));
}
if(B>=13 && B<16)
{
	ResRP3->Fill(14.5,cos(atan2(QyR3,QxR3)-PhiRp));
}

*/

	/*
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
        fn90all[3]->Fill(atan2(Qy901,Qx901));
        ResPhin[2]=ResPhin[2]+cos(atan2(Qy901-Recy90[3],Qx901-Recx90[3])-atan2(Qy441-Recy44[3],Qx441-Recx44[3]));
        ResN[2]=ResN[2]+1;
	ResRP->Fill(14.5,cos(atan2(QyR3,QxR3)-PhiRp));
    //  flowvsB->Fill(14.5,(cos(module_pos.GetPhi() - atan2(Qy901 + Qy441 - Recy44[3] - Recy90[3],Qx901 + Qx441 - Recx44[3] -Recx90[3]))));
        }*/
	
	
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

  Phit->Write();

  PID_proton->Write();
  dEdX_vs_pz_proton->Write();
  mass2_vs_pz_proton->Write();
  dEdX_vs_pz_kaon->Write();
  mass2_vs_pz_kaon->Write();
  dEdX_vs_pz_pion->Write();
  mass2_vs_pz_pion->Write();
  dEdX_vs_pz_all->Write();


  PhiModule->Write();
  ResRPS->Write();
  ResRPN->Write();
  ResRP3->Write();
  ResHalf->Write();
 for(int k=0;k<3;k++)
 {

TpcflowvsB[k]->Write();
TpcPhivsB[k]->Write();
McPhivsB[k]->Write();
McflowvsB[k]->Write();
 }
/*simflowvsB->Write();*/

for(int k=0;k<3;k++)
{
	for(int j=0;j<3;j++)
	{
	TpcflowvspT[k][j]->Write();
	TpcflowvsEta[k][j]->Write();
	McflowvspT[k][j]->Write();
	McflowvsEta[k][j]->Write();
//	simflowvsEta[j]->Write();
//	simflowvspT[j]->Write();
	TpcflowvsEtapT[k][j]->Write();
        McflowvsEtapT[k][j]->Write();
//	simflowvsEtapT[j]->Write();

}
for(int j=0;j<2;j++)
{
	TpcflowvspTEta[k][j]->Write();
        McflowvspTEta[k][j]->Write();
//	simflowvspTEta[j]->Write();

}
}

  pT_distribution_->Write();
  fhcal_energy_distribution_->Write();
  fhcal_modules_xy_->Write();
  TpcpT_vs_eta_proton->Write();
  TpcpT_vs_eta_kaon->Write();
TpcpT_vs_eta_pion->Write();
TpcpT_vs_eta_all->Write();

  TpcpT_vs_phi->Write();
  Tpcphi_vs_eta->Write();
McpT_vs_eta_proton->Write();
McpT_vs_eta_kaon->Write();
McpT_vs_eta_pion->Write();
McpT_vs_eta_all->Write();

  McpT_vs_phi->Write();
  Mcphi_vs_eta->Write();

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
    for(int j=1;j<17;j++)
    {
    Qx44all[j]->Write();
    }
    for(int j=1;j<17;j++)
    {

    Qy44all[j]->Write();
    }
    for(int j=1;j<17;j++)
    {
    Qx90all[j]->Write();
    }
    for(int j=1;j<17;j++)
    {
    Qy90all[j]->Write();
    }
    Bt->Write();
    Btt->Write();
    energyvsb->Write();
    for(int j=1;j<17;j++)
    {
    fn90all[j]->Write();
    }
for(int j=1;j<17;j++)
    {
    fn44all[j]->Write();
    }

     
for(int j=0;j<16;j++)
{
	std::cout<<ResPhin[j]<<" /// "<< ResN[j]<<"next";
}





}
//
// Created by mikhail on 6/16/20.
