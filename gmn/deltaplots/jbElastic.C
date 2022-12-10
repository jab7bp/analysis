#include <iostream>
#include <fstream>
#include <numeric>
#include <iomanip>
#include <math.h>
#include <algorithm>
#include <string>
#include <chrono>

using namespace std::chrono;
#include "/w/halla-scshelf2102/sbs/jboyd/include/include_files.h"
#include "/w/halla-scshelf2102/sbs/jboyd/include/GEM_lookups.h"
#include "/w/halla-scshelf2102/sbs/jboyd/include/beam_variables.h"

int runnum = 13581;
TString experiment = "gmn";
int pass = 1;

int kine = lookup_kinematic(runnum); //Kinematic setting: 4 for SBS-4, 9 for SBS-9, etc.
TString run_target = lookup_target(runnum); //Target --> "LH2" or "LD2"

TString rootfile_dir = Form("/work/halla/sbs/sbs-gmn/pass%i/SBS%i/%s/rootfiles", pass, kine, run_target.Data());
TString rootfile = Form("%s/*%i*", rootfile_dir.Data(), runnum);

void jbElastic(){

	const char *rootfilename = rootfile.Data();
	double Ebeam = lookup_beam_energy(runnum);
	double bbtheta = lookup_BB_angle(runnum);
	double sbstheta = lookup_SBS_angle(runnum);
	double sbsfieldscale = lookup_run_info(runnum, "sbs_field");
	double hcaldist = lookup_HCal_dist(runnum);

	cout << "----------------------------------" << endl;
	cout << "Run parameters: " << endl;
	cout << "Runnum = " << runnum << endl;
	cout << "Beam energy = " << Ebeam << endl;
	cout << "BB angle = " << bbtheta << endl;
	cout << "SBS angle = " << sbstheta << endl;
	cout << "SBS field scale = " << sbsfieldscale << endl;
	cout << "SBS distance = " << sbsdist << endl;
	cout << "HCal distance = " << hcaldist << endl;
	cout << "----------------------------------" << endl;

	sbstheta *= TMath::Pi()/180.0;
	bbtheta *= TMath::Pi()/180.0;

	TChain *C = new TChain("T");

	C->Add(rootfilename); 

	TCut globalcut = "bb.ps.e>0.15&&abs(bb.tr.vz)<0.075&&sbs.hcal.nclus>0&&bb.tr.n==1";

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

	double W2min = 0.88-0.4;
	double W2max = 0.88+0.4;

	TFile *fout = new TFile("jbElastic_out.root","RECREATE");

	TH2D *hdxdy_all = new TH2D("hdxdy_all","All events;#Deltay (m);#Deltax (m)",125,-2,2,125,-4,6);
	TH2D *hdxdy_Wcut = new TH2D("hdxdy_Wcut","|W^{2}-0.88|<0.4;Deltay (m);#Deltax (m)",125,-2,2,125,-4,6);
	TH2D *hdxdy_W_anticut = new TH2D("hdxdy_W_anticut","|W^{2}-0.88|<0.4;Deltay (m);#Deltax (m)",125,-2,2,125,-4,6);

	TH1D *hW2_all = new TH1D("hW2_all","All events; W^{2} (GeV^{2});", 250, -1, 4 );

	TH1D *hEHCAL_all = new TH1D("hEHCAL_all","All events; HCAL energy sum (GeV);",250,0,0.5);
	TH1D *hEHCAL_Wcut = new TH1D("hEHCAL_Wcut","|W^{2}-0.88|<0.4;HCAL energy sum (GeV);",250,0,0.5);

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
		  } 
		  else {
			hdxdy_W_anticut->Fill( yhcal - yhcal_expect, xhcal - xhcal_expect );
		  }
		}
	}

	TCanvas *c_dxdy_Wcut = new TCanvas("c_dxdy_Wcut", "dxdy with W cuts", 600, 500);
	hdxdy_Wcut->Draw("colz");

	elist->Delete(); 
	fout->Write();

} 