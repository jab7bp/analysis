#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCut.h"
#include "TEventList.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TMath.h"
#include <chrono>

using namespace std::chrono;
#include "/w/halla-scshelf2102/sbs/jboyd/include/include_files.h"
#include "/w/halla-scshelf2102/sbs/jboyd/include/GEM_lookups.h"
#include "/w/halla-scshelf2102/sbs/jboyd/include/beam_variables.h"

//Run info and lookups
int runnum = 11581;
TString experiment = "gmn";
int pass = 0;

int kine = lookup_kinematic(runnum); //Kinematic setting: 4 for SBS-4, 9 for SBS-9, etc.
TString run_target = lookup_target(runnum); //Target --> "LH2" or "LD2"
double Ebeam = lookup_beam_energy(runnum); //Electron beam energy (electron energy) in GeV
double sbsfieldscale = lookup_run_info(runnum, "sbs_field"); //Strength (in percentage) of SBS magnet
double bbtheta=35.9;
double sbstheta=22.1;
double sbsdist=2.80;
double hcaldist=17.0;

// TString rootfile_dir = "/lustre19/expphy/volatile/halla/sbs/jboyd/Rootfiles/xtalk/11449";
TString rootfile_dir = Form("/lustre19/expphy/volatile/halla/sbs/sbs-gmn/GMN_REPLAYS/pass0/SBS%i/%s/rootfiles", kine, run_target.Data());
TString input_rootfile = Form("%s/*%i*", rootfile_dir.Data(), runnum);

const char *rootfilename = input_rootfile.Data();

void ElasticQuickAndDirty(){

  auto total_time_start = high_resolution_clock::now();

  sbstheta *= TMath::Pi()/180.0;
  bbtheta *= TMath::Pi()/180.0;

  TChain *C = new TChain("T");

  C->Add(rootfilename); 

  TCut globalcut = "sbs.hcal.nclus>0&&bb.ps.e>0.120&&bb.sh.nclus>0&&bb.ps.nclus>0&&sbs.hcal.e>0.05";
  
  TEventList *elist = new TEventList("elist","");
  
  C->Draw(">>elist",globalcut);

  int MAXNTRACKS=10;
  
  //variables we need are BigBite track px,py,pz and sbs hcal x, y, e

  double ntrack;

  double epx[MAXNTRACKS];
  double epy[MAXNTRACKS];
  double epz[MAXNTRACKS];
  double ep[MAXNTRACKS];

  double vx[MAXNTRACKS];
  double vy[MAXNTRACKS];
  double vz[MAXNTRACKS];

  double xhcal,yhcal,ehcal;

  C->SetBranchStatus("*",0);
  C->SetBranchStatus("bb.tr.n",1);
  C->SetBranchStatus("bb.tr.vz",1);
  C->SetBranchStatus("bb.tr.px",1);
  C->SetBranchStatus("bb.tr.py",1);
  C->SetBranchStatus("bb.tr.pz",1);
  C->SetBranchStatus("bb.tr.p",1);

  C->SetBranchStatus("sbs.hcal.x",1);
  C->SetBranchStatus("sbs.hcal.y",1);
  C->SetBranchStatus("sbs.hcal.e",1);

  C->SetBranchAddress("bb.tr.n",&ntrack);
  C->SetBranchAddress("bb.tr.vz",vz);
  C->SetBranchAddress("bb.tr.px",epx);
  C->SetBranchAddress("bb.tr.py",epy);
  C->SetBranchAddress("bb.tr.pz",epz);
  C->SetBranchAddress("bb.tr.p",ep);
  C->SetBranchAddress("sbs.hcal.x",&xhcal);
  C->SetBranchAddress("sbs.hcal.y",&yhcal);
  C->SetBranchAddress("sbs.hcal.e",&ehcal);

  TLorentzVector Pbeam(0,0,Ebeam,Ebeam);
  TLorentzVector Ptarg(0,0,0,0.5*(0.938272+0.939565));

  long nevent=0;
  double W2nom = 0.916;
  double W2var = 0.2;
  double W2min = W2nom-W2var;
  double W2max = W2nom+W2var;

  TFile *fout = new TFile("elastic_temp.root","RECREATE");

  TH2D *hdxdy_all = new TH2D("hdxdy_all","All events;#Deltay (m);#Deltax (m)",125,-2,2,125,-4,6);
  TH2D *hdxdy_Wcut = new TH2D("hdxdy_Wcut",Form("|W^{2}-%.3f|<%.2f;Deltay (m);#Deltax (m)", W2nom, W2var),125,-4,4,125,-4,6);
  TH2D *hdxdy_Wanticut = new TH2D("hdxdy_Wanticut",Form("|W^{2}-%.3f|<%.2f;Deltay (m);#Deltax (m)", W2nom*W2nom, W2var*W2var),125,-2,2,125,-4,6);

  TH1D *hW2_all = new TH1D("hW2_all","All events; W^{2} (GeV^{2});", 250, -1, 4 );

  TH1D *hEHCAL_all = new TH1D("hEHCAL_all","All events; HCAL energy sum (GeV);",250,0,0.5);
  TH1D *hEHCAL_Wcut = new TH1D("hEHCAL_Wcut",Form("|W^{2}-%.3f|<%.2f;HCAL energy sum (GeV);", W2nom, W2var),250,0,0.5);

  TVector3 hcal_origin( -hcaldist*sin(sbstheta), 0, hcaldist*cos(sbstheta) );

  TVector3 hcal_zaxis = hcal_origin.Unit();
  TVector3 hcal_xaxis(0,-1,0);
  TVector3 hcal_yaxis = hcal_zaxis.Cross( hcal_xaxis ).Unit();
  
  while( C->GetEntry(elist->GetEntry(nevent++) ) ){
    if( ntrack == 1.0 ){
      TLorentzVector kprime( epx[0], epy[0], epz[0], ep[0] );
      TLorentzVector q = Pbeam - kprime;

      TVector3 qdir = q.Vect().Unit();

      TVector3 vertex(0,0,vz[0]);

      double sintersect = (hcal_origin-vertex).Dot( hcal_zaxis )/qdir.Dot( hcal_zaxis );

      TVector3 hcal_intersect = vertex + sintersect * qdir; 

      double xhcal_expect = hcal_intersect.Dot( hcal_xaxis );
      double yhcal_expect = hcal_intersect.Dot( hcal_yaxis );

      hdxdy_all->Fill( yhcal - yhcal_expect, xhcal - xhcal_expect );
      double W2recon = (Ptarg + q).M2();

      hW2_all->Fill( W2recon );

      hEHCAL_all->Fill( ehcal );

      if( W2recon > W2min && W2recon < W2max ){
	hdxdy_Wcut->Fill( yhcal - yhcal_expect, xhcal - xhcal_expect );
	hEHCAL_Wcut->Fill( ehcal );
      } else {
	hdxdy_Wanticut->Fill( yhcal - yhcal_expect, xhcal - xhcal_expect );
      }
    }
  }

  elist->Delete(); 
  fout->Write();

  auto total_time_end = high_resolution_clock::now();
  auto total_time_duration = duration_cast<minutes>(total_time_end - total_time_start);
  cout << "Total time for analysis: " << total_time_duration.count() << " minutes. " << endl;

} 