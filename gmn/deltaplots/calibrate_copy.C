#include <iostream>
#include <fstream>
#include <numeric>
#include <iomanip>
#include <math.h>
#include <algorithm>
#include <string>
#include <chrono>
#include <TF1.h>

using namespace std::chrono;
#include "/w/halla-scshelf2102/sbs/jboyd/include/include_files.h"
#include "/w/halla-scshelf2102/sbs/jboyd/include/GEM_lookups.h"
#include "/w/halla-scshelf2102/sbs/jboyd/include/beam_variables.h"

Double_t fit_gaus(Double_t * x, Double_t *par){

	Double_t g = 0.0;

	g = par[0]*exp((-0.5)*pow(((x[0] -  par[1])/par[2]),2));
	
	return g;
}

bool single_run = true;
bool multi_run = false;

bool calc_W = false;

//Run info and lookups
int runnum = 13585;
vector<int> runnum_vec = {13585, 13586, 13587};
TString experiment = "gmn";
int pass = 1;

int kine = lookup_kinematic(runnum); //Kinematic setting: 4 for SBS-4, 9 for SBS-9, etc.
TString run_target = lookup_target(runnum); //Target --> "LH2" or "LD2"
double E_beam = lookup_beam_energy(runnum); //Electron beam energy (electron energy) in GeV
double SBS_field = lookup_run_info(runnum, "sbs_field"); //Strength (in percentage) of SBS magnet
int sbsfieldscale = int(100*lookup_run_info(runnum, "sbs_field"));

// TString rootfile_dir = "/lustre19/expphy/volatile/halla/sbs/jboyd/Rootfiles/xtalk/11449";
TString rootfile_dir = Form("/lustre19/expphy/volatile/halla/sbs/sbs-gmn/GMN_REPLAYS/pass%i/SBS%i/%s/rootfiles", pass, kine, run_target.Data());
// TString input_rootfile;

TFile *outfile;
TChain *TC = new TChain("T");
vector<TString> master_cut_vec;
TString master_cut_string;
TCut master_cut = "";

//Experimental Constants, Thresholds, cuts, etc DEFINITIONS
const double pi = TMath::Pi();
const double Mp = 0.938272; //Mass of proton [GeV]
const double Mn = 0.939565; //Mass of neutron [GeV]
double W2_mean; //Invariant Mass-squared (mean val) {With perfect optics W2 = Mp. Can be calculated run-by-run}
double W2_sigma; //Invariant Mass-squared sigma {Reasonable default/guess. Can be calculated run-by-run from W plot}

//HCal constants and stuff
double tdiff = 510;		//Time difference between BBCal and HCal signals
double tdiff_max = 20;	//Maximum time difference from coincidences through tdctrig cut
double HCal_dist = lookup_HCal_dist( runnum ); 	//Distace from HCal face to target chamber
double Hcal_theta = lookup_HCal_angle( runnum );		//Theta angle for HCal from downstream beamline

const Double_t sampfrac = 0.077; 	//Estimate of the sampling fraction from MC
const Int_t kNcell = 288; // Total number of HCal modules
const Int_t kNrows = 24; // Total number of HCal rows
const Int_t kNcols = 12; // Total number of HCal columns
const Int_t kNtrack = 100; // Reasonable max number of tracks per event
const Int_t kNtdc = 1000; // Reasonable max number of tdc signals per event
const Double_t Xi = -2.20; // Distance from beam center to top of HCal in m
const Double_t Xf = 1.47; // Distance from beam center to bottom of HCal in m
const Double_t Yi = -0.853; // Distance from beam center to opposite-beam side of HCal in m
const Double_t Yf = 0.853; // Distance from beam center to beam side of HCal in m

//Declare vars
Double_t atime[kNcell], row[kNcell], col[kNcell], tdctime[kNcell], cblkid[kNcell], cblke[kNcell];
Double_t nblk, nclus, SH_nclus, PS_nclus, hcal_x, hcal_y, hcal_e;

Double_t par[3];

Double_t bb_tr_p[kNtrack], bb_tr_px[kNtrack], bb_tr_py[kNtrack], bb_tr_pz[kNtrack];
Double_t bb_tr_vz[kNtrack];
Double_t bb_tr_n, bb_ps_x, bb_ps_y, bb_ps_e, bb_sh_x, bb_sh_y, bb_sh_e;

Double_t TDCT_id[kNtdc], TDCT_tdc[kNtdc], hodo_tmean[kNtdc]; 
Int_t TDCTndata, elastic_yield;

Long64_t Nevents;

//INITIALIZE ALL HISTOGRAMS:
TH1D *h_Ep, *h_PS, *h_HCal_e, *h_SHPS, *h_W2, *h_W2recon, *h_W, *h_Wrecon;
double Ep_center, Ep_sigma, PS_center, PS_sigma, PS_min, SHPS_center, SHPS_sigma, HCal_e, HCal_e_min, W2_fit_center, W2_fit_sigma;
vector<double> cuts_from_fits;

// TH1D *h_atime, *h_E, *h_E_cut, *h_E_exp, *h_vert, *h_W2, *h_W2recon, *h_dx, *h_dy, *h_KE_p, *h_KElow, *h_diff, *h_X, *h_Y;
// TH2D *h2_CvCh, * h2_EvCh, *h2_dxdy, *h2_dxdy_all, *h2_dxVE, *h2_rowcol;

void calibrate(){
	auto total_time_start = high_resolution_clock::now();
	cout << "--------------------------------------" << endl;
	cout << "Analysis started. " << endl;
	cout << "--------------------------------------" << endl;

	outfile = new TFile(Form("rootfiles/%i_dxdy.root", runnum), "RECREATE");

	//Define Histograms:
	h_Ep = new TH1D("Ep", Form("E/p - SBS = %i%%, %s, Run %i", sbsfieldscale, run_target.Data(), runnum), 200, 0, 2);
	h_PS = new TH1D("PS", Form("Pre-Shower Clus. E - SBS = %i%%, %s, Run %i; GeV", sbsfieldscale, run_target.Data(), runnum), 300, 0, 3);
	h_HCal_e = new TH1D("HCal_e", Form("HCal Clus. E - SBS = %i%%, %s, Run %i; GeV", sbsfieldscale, run_target.Data(), runnum), 200, 0, 0.4);
	h_SHPS = new TH1D("SHPS", Form("SH + PS Clus. E - SBS = %i%%, %s, Run %i; GeV", sbsfieldscale, run_target.Data(), runnum), 500, 0, 5);
	h_W2 = new TH1D( "W2", Form("Invariant Mass, W^{2} (No Cuts) - SBS = %i%%, %s, Run %i; GeV", sbsfieldscale, run_target.Data(), runnum), 200, 0.0, 4.0 );
	h_W2recon = new TH1D( "W2recon", Form("Reconstructed W^{2} (No Cuts) - SBS = %i%%, %s, Run %i; GeV", sbsfieldscale, run_target.Data(), runnum), 200, 0.0, 4.0 );
	h_W = new TH1D( "W", Form("Invariant Mass, W (No Cuts) - SBS = %i%%, %s, Run %i; GeV", sbsfieldscale, run_target.Data(), runnum), 200, 0.0, 4.0 );
	h_Wrecon = new TH1D( "W2recon", Form("Reconstructed W (No Cuts) - SBS = %i%%, %s, Run %i; GeV", sbsfieldscale, run_target.Data(), runnum), 200, 0.0, 4.0 );


	if( single_run ){
		cout << "Running in single run mode: " << runnum << endl;
		cout << "--------------------------------------" << endl;
		cout << "Adding files to TChain from: " << rootfile_dir.Data() << endl;
		TC->Add(Form("%s/*%i*.root", rootfile_dir.Data(), runnum));
	}
cout << "--------------------------------------" << endl;
	if( multi_run ){
		runnum = runnum_vec[0];
		cout << "Running in multi-run mode for runs: "  << endl;
		for(size_t run = 0; run < runnum_vec.size(); run++){
			cout << runnum_vec[run] << " ";
		}
		cout << endl;
		cout << "--------------------------------------" << endl;
		cout << "Adding files to TChain from: " << rootfile_dir.Data() << endl;
		for(size_t run = 0; run < runnum_vec.size(); run++){
			TC->Add(Form("%s/*%i*.root", rootfile_dir.Data(), runnum_vec[run]));
		}
	}
	cout << "Finished adding files to TChain. " << endl;
	cout << "--------------------------------------" << endl;

	if( calc_W ){
		W2_mean = lookup_cut(runnum, "W2");
		W2_sigma = lookup_cut(runnum, "W2");
	}
	else{
		W2_mean = 0.92;
		W2_sigma = 0.3;
	}

	//DEFINE CUTS
	if( !calc_W ){
		master_cut_vec = {
			"sbs.hcal.nclus>0",
			"bb.ps.nclus>0",
			"bb.sh.nclus>0",
			"abs(bb.tr.vz[0])<=0.075",
			"bb.gem.track.nhits[0]>3",
			"bb.tr.n==1",
		};
	}
	if( calc_W ){
		master_cut_vec = {
			"sbs.hcal.nclus>0",
			"bb.ps.nclus>0",
			"bb.sh.nclus>0",
			"abs(bb.tr.vz[0])<=0.075",
			"bb.gem.track.nhits[0]>3",
			"bb.tr.n==1",
			Form("bb.ps.e>%f", lookup_cut(runnum, "PS_clus_e_cut")),
			Form("((bb.sh.e+bb.ps.e)/(bb.tr.p[0]))>%f", lookup_cut(runnum, "SH_PS_clus_e_cut") - lookup_cut(runnum,"SH_PS_sigma")),
			Form("sbs.hcal.e>%f", lookup_cut(runnum, "HCal_clus_e_cut"))
		};
	}
	for(size_t cut = 0; cut < master_cut_vec.size(); cut++){
		if(cut == master_cut_vec.size() - 1){
			master_cut_string.Append(Form("%s", master_cut_vec[cut].Data()));
		}
		else{
			master_cut_string.Append(Form("%s%s", master_cut_vec[cut].Data(), "&&"));
		}
	}

	master_cut = Form("%s", master_cut_string.Data());
	cout << "--------------------------------------" << endl;
	cout << "Setting up branches... ";
	//Setup Branches
	TC->SetMakeClass(1); //Allows for viewing of general detector params from Event_Branch
	TC->SetBranchStatus( "*", 0 );
	TC->SetBranchStatus( "sbs.hcal.clus_blk.atime", 1 );  
	TC->SetBranchStatus( "sbs.hcal.clus_blk.row", 1 );
	TC->SetBranchStatus( "sbs.hcal.clus_blk.col", 1 );
	TC->SetBranchStatus( "sbs.hcal.clus_blk.tdctime", 1 );
	TC->SetBranchStatus( "sbs.hcal.clus_blk.id", 1 );
	TC->SetBranchStatus( "sbs.hcal.clus_blk.e", 1 );
	TC->SetBranchStatus( "sbs.hcal.x", 1 );
	TC->SetBranchStatus( "sbs.hcal.y", 1 );
	TC->SetBranchStatus( "sbs.hcal.e", 1 );
	TC->SetBranchStatus( "sbs.hcal.nblk", 1 );
	TC->SetBranchStatus( "sbs.hcal.nclus", 1 );
	TC->SetBranchStatus( "bb.sh.nclus", 1 );
	TC->SetBranchStatus( "bb.ps.nclus", 1 );
	TC->SetBranchStatus( "bb.tr.n", 1 );
	TC->SetBranchStatus( "bb.tr.px", 1 );
	TC->SetBranchStatus( "bb.tr.py", 1 );
	TC->SetBranchStatus( "bb.tr.pz", 1 );
	TC->SetBranchStatus( "bb.tr.p", 1 );
	TC->SetBranchStatus( "bb.tr.vz", 1 );
	TC->SetBranchStatus( "bb.ps.e", 1 );
	TC->SetBranchStatus( "bb.ps.x", 1 );
	TC->SetBranchStatus( "bb.ps.y", 1 );
	TC->SetBranchStatus( "bb.sh.e", 1 );
	TC->SetBranchStatus( "bb.sh.x", 1 );
	TC->SetBranchStatus( "bb.sh.y", 1 );
	TC->SetBranchStatus( "bb.tdctrig.tdc", 1 );
	TC->SetBranchStatus( "bb.tdctrig.tdcelemID", 1 );
	TC->SetBranchStatus( "Ndata.bb.tdctrig.tdcelemID", 1 );
	TC->SetBranchStatus( "bb.gem.track.nhits", 1);

	TC->SetBranchAddress( "sbs.hcal.clus_blk.atime", atime );
	TC->SetBranchAddress( "sbs.hcal.clus_blk.row", row );
	TC->SetBranchAddress( "sbs.hcal.clus_blk.col", col );
	TC->SetBranchAddress( "sbs.hcal.clus_blk.tdctime", tdctime );
	TC->SetBranchAddress( "sbs.hcal.clus_blk.id", cblkid );
	TC->SetBranchAddress( "sbs.hcal.clus_blk.e", cblke );
	TC->SetBranchAddress( "sbs.hcal.x", &hcal_x );
	TC->SetBranchAddress( "sbs.hcal.y", &hcal_y );
	TC->SetBranchAddress( "sbs.hcal.e", &hcal_e );
	TC->SetBranchAddress( "sbs.hcal.nblk", &nblk );
	TC->SetBranchAddress( "sbs.hcal.nclus", &nclus );
	TC->SetBranchAddress( "bb.sh.nclus", &SH_nclus );
	TC->SetBranchAddress( "bb.ps.nclus", &PS_nclus );
	TC->SetBranchAddress( "bb.tr.n", &bb_tr_n );
	TC->SetBranchAddress( "bb.tr.px", bb_tr_px );
	TC->SetBranchAddress( "bb.tr.py", bb_tr_py );
	TC->SetBranchAddress( "bb.tr.pz", bb_tr_pz );
	TC->SetBranchAddress( "bb.tr.vz", bb_tr_vz );
	TC->SetBranchAddress( "bb.tr.p", bb_tr_p );
	TC->SetBranchAddress( "bb.ps.e", &bb_ps_e );
	TC->SetBranchAddress( "bb.ps.x", &bb_ps_x );
	TC->SetBranchAddress( "bb.ps.y", &bb_ps_y );
	TC->SetBranchAddress( "bb.sh.e", &bb_sh_e );
	TC->SetBranchAddress( "bb.sh.x", &bb_sh_x );
	TC->SetBranchAddress( "bb.sh.y", &bb_sh_y );
	TC->SetBranchAddress( "bb.tdctrig.tdcelemID", TDCT_id );
	TC->SetBranchAddress( "bb.tdctrig.tdc", TDCT_tdc );
	TC->SetBranchAddress( "Ndata.bb.tdctrig.tdcelemID", &TDCTndata );
	cout << " done. " << endl;
	cout << "--------------------------------------" << endl;
	cout << "Applying master cut: " << endl;
	cout << master_cut_string.Data() << endl;
	cout << "--------------------------------------" << endl;

	TEventList *ev_list = new TEventList("ev_list", "Elastic Event List");

//Event list with cuts:
	TC->Draw(">>ev_list", master_cut);
	cout << "--------------------------------------" << endl;
	cout << "--------------------------------------" << endl;
	Nevents = ev_list->GetN();

	cout << "Number of events to analyze: " << Nevents << endl;
	cout << "--------------------------------------" << endl;
	cout << "--------------------------------------" << endl;
	cout << "Starting analysis loop on events..... " << endl;

	for(Long64_t nevent = 0; nevent < Nevents; nevent++){
		TC->GetEntry( ev_list->GetEntry( nevent ));
		if( nevent%100000 == 0){
			cout << "Evt: " << nevent <<"/" << Nevents << Form("(%.1f/100%%)", 100.0*double(1.0*nevent/Nevents)) << endl;
		}
		cout.flush();

		Double_t Ep = (bb_ps_e + bb_sh_e)/(bb_tr_p[0]);
		Double_t PS = bb_ps_e;
		Double_t SHPS = bb_ps_e + bb_sh_e;
		Double_t HCal_e = hcal_e;

		h_Ep->Fill(Ep);
		h_PS->Fill(PS);
		h_SHPS->Fill(SHPS);
		h_HCal_e->Fill(hcal_e);

//----------------------------------------------------
//---------------INVARIANT MASS -----------------------
		if( calc_W ){
			TLorentzVector kprime( bb_tr_px[0], bb_tr_py[0], bb_tr_pz[0], bb_tr_p[0] );
			TLorentzVector q = Pbeam - kprime; //Standard q-vector
			TVector3 vertex( 0, 0, bb_tr_vz[0] ); //Location of scattering event in target
			TVector3 qunit = q.Vect().Unit(); //q-vector direction

			Double_t sintersect = (hcal_origin-vertex).Dot( hcal_zaxis )/qunit.Dot( hcal_zaxis );
			TVector3 hcal_intersect = vertex + sintersect * qunit;
			Double_t xexpect = hcal_intersect.Dot( hcal_xaxis );
			Double_t yexpect = hcal_intersect.Dot( hcal_yaxis );
			Double_t W2recon = (Ptarg + q).M2();
			Double_t Wrecon = (Ptarg + q).M();
			Double_t E_ep = bb_tr_p[0]; // Obtain the scattered electron energy, neglect mass e
			Double_t p_ep = bb_tr_p[0]; // Obtain the magnitude of scattered electron momentum
			Double_t Q2 = 2*E_e*E_ep*( 1-(bb_tr_pz[0]/p_ep) );
			Double_t nu = E_e-E_ep; // Obtain energy transfer
			Double_t W2 = pow( M_p,2 )+2*M_p*nu-Q2; // Obtain W2 from Q2 and nu
			Double_t W = sqrt(W2);
			Double_t eth = acos( bb_tr_pz[0]/bb_tr_p[0] );
			Double_t eph = atan2( bb_tr_py[0], bb_tr_px[0] );
			Double_t phinucleon = eph + TMath::Pi(); //assume coplanarity
			Double_t thetanucleon = acos( (E_e - bb_tr_pz[0])/p_ep ); //use elastic constraint on nucleon kinematics
			TVector3 pNhat( sin(thetanucleon)*cos(phinucleon),sin(thetanucleon)*sin(phinucleon),cos(thetanucleon));
			Double_t KE_p = nu; //For elastics
			hKE_p->Fill( KE_p );   

			//Cut on BBCal and HCal trigger coincidence
			Double_t bbcal_time=0., hcal_time=0., coin_time=0., rf_time=0.;

		    Double_t diff = hcal_time - bbcal_time; 
		    h_Diff->Fill( diff );

		    hdxdy_all->Fill( hcaly - yexpect, hcalx - xexpect );
		    h_vert->Fill(BBtr_vz[0]);

		    h_W2->Fill(W2);
		    h_W2recon->Fill(W2recon);

		    h_Wrecon->Fill(Wrecon);
		}

	}
	cout << "---------------------------------------" << endl;
	cout << "-----Finished going through events-----" << endl;
	cout << "---------------------------------------" << endl;

//---------------------------------------
//---------------- Ep -----------------------

	TCanvas *c_Ep = new TCanvas("c_Ep", "E/p", 600, 500);
	h_Ep->Draw();
	TF1 *fit_Ep = new TF1("fit_Ep", fit_gaus, 0.5, 1.5, 3);
	fit_Ep->SetParName(0, "Ep Norm");
	fit_Ep->SetParName(1, "Ep Center");
	fit_Ep->SetParName(2, "Ep Sigma");

	fit_Ep->SetParLimits(0, 0, h_Ep->GetMaximum());
	fit_Ep->SetParLimits(1, 0.8, 1.2);
	fit_Ep->SetParLimits(2, 0.01, 0.4);

	h_Ep->Fit("fit_Ep", "R+");
	fit_Ep->Draw("same");
	Ep_center = fit_Ep->GetParameter(1);
	Ep_sigma = fit_Ep->GetParameter(2);

	// par[0] = 0.0; par[1] = 0.0, par[2] = 0.0;


//---------------------------------------
//---------------- PS -----------------------
	TCanvas *c_PS = new TCanvas("c_PS", "PS", 600, 500);
	h_PS->Draw();

	TF1 *fit_PS = new TF1("fit_PS", fit_gaus, 0.15, 3, 3);
	fit_PS->SetParName(0, "PS Norm");
	fit_PS->SetParName(1, "PS Center");
	fit_PS->SetParName(2, "PS Sigma");

	fit_PS->SetParLimits(0, 0, h_PS->GetMaximum());
	fit_PS->SetParLimits(1, 0.8, 1.2);
	fit_PS->SetParLimits(2, 0.01, 0.8);

	h_PS->Fit("fit_PS", "R+");
	// fit_PS->Draw("same");
	PS_center = fit_PS->GetParameter(1);
	PS_sigma = fit_PS->GetParameter(2);

	// par[0] = 0.0; par[1] = 0.0, par[2] = 0.0;

	cout << "Setting range of PS to : " << (3.0/300.0)*h_PS->FindFirstBinAbove((0.75)*(h_PS->GetMaximum())) << " - " << PS_center << endl;
	h_PS->GetXaxis()->SetRangeUser((3.0/300.0)*h_PS->FindFirstBinAbove((0.75)*(h_PS->GetMaximum())), PS_center);
	PS_min = (3.0/300.0)*h_PS->GetMinimumBin();
	PS_min = ceil(100*PS_min)/100.0;

	h_PS->GetXaxis()->SetRangeUser(0, 3);
	TLine *tl_PS = new TLine(PS_min, 0, PS_min, h_PS->GetMaximum());
	tl_PS->SetLineColor(6);
	tl_PS->Draw("same");

	cout << "---------------------------------------" << endl;
	cout << "PS_min set to: " << PS_min << endl;
	cout << "---------------------------------------" << endl;

//---------------------------------------
//---------------- SHPS -----------------------
	TCanvas *c_SHPS = new TCanvas("c_SHPS", "SHPS", 600, 500);
	h_SHPS->Draw();
	TF1 *fit_SHPS = new TF1("fit_SHPS", fit_gaus, 0.5, 5, 3);
	fit_SHPS->SetParName(0, "SHPS Norm");
	fit_SHPS->SetParName(1, "SHPS Center");
	fit_SHPS->SetParName(2, "SHPS Sigma");

	fit_SHPS->SetParLimits(0, 0, h_SHPS->GetMaximum());
	fit_SHPS->SetParLimits(1, 2.6, 3.4);
	fit_SHPS->SetParLimits(2, 0.01, 0.5);

	h_SHPS->Fit("fit_SHPS", "R+");
	fit_SHPS->Draw("same");
	SHPS_center = fit_SHPS->GetParameter(1);
	SHPS_sigma = fit_SHPS->GetParameter(2);

	// par[0] = 0.0; par[1] = 0.0, par[2] = 0.0;

//---------------------------------------
//---------------- HCal_e -----------------------
	TCanvas *c_HCal_e = new TCanvas("c_HCal_e", "HCal_e", 600, 500);
	h_HCal_e->Draw();
	HCal_e_min = (0.4/200.0)*h_HCal_e->FindFirstBinAbove((0.25)*(h_HCal_e->GetMaximum()));
	HCal_e_min = ceil(100*HCal_e_min)/100.0;

	cout << "---------------------------------------" << endl;
	cout << "HCal_e_min set to: " << HCal_e_min << endl;
	cout << "---------------------------------------" << endl;

	TLine *tl_HCal_e = new TLine(HCal_e_min, 0, HCal_e_min, h_HCal_e->GetMaximum());
	tl_HCal_e->SetLineColor(6);
	tl_HCal_e->Draw("same");

	// par[0] = 0.0; par[1] = 0.0, par[2] = 0.0;
//---------------------------------------
//---------------- W2 Recon -----------------------	
	TCanvas *c_W2recon = new TCanvas("c_W2recon", "W2recon", 600, 500);
	h_W2recon->Draw();
	TF1 *fit_W2recon = new TF1("fit_W2recon", fit_gaus, 0.5, 1.25, 3);
	fit_W2recon->SetParName(0, "W2recon Norm");
	fit_W2recon->SetParName(1, "W2recon Center");
	fit_W2recon->SetParName(2, "W2recon Sigma");

	fit_W2recon->SetParLimits(0, 0, h_W2recon->GetMaximum());
	fit_W2recon->SetParLimits(1, 0.85, 1.15);
	fit_W2recon->SetParLimits(2, 0.01, 0.2);

	h_W2recon->Fit("fit_W2recon", "R+");
	fit_W2recon->Draw("same");
	W2_fit_center = fit_W2recon->GetParameter(1);
	W2_fit_sigma = fit_W2recon->GetParameter(2);

	// par[0] = 0.0; par[1] = 0.0, par[2] = 0.0;

//---------------------------------------
//---------------- CUT VECTOR -----------------------	
	cuts_from_fits.push_back(runnum);
	cuts_from_fits.push_back(PS_min);
	cuts_from_fits.push_back(SHPS_center);
	cuts_from_fits.push_back(SHPS_sigma);
	cuts_from_fits.push_back(HCal_e_min);
	cuts_from_fits.push_back(Ep_center);
	cuts_from_fits.push_back(Ep_sigma);
	cuts_from_fits.push_back(W2_fit_center);
	cuts_from_fits.push_back(W2_fit_sigma);

	cout << "---------------------------------------" << endl;
	cout << "---------------------------------------" << endl << endl;	
	cout << "Vector with cuts: " << endl;
	cout << cuts_from_fits << endl << endl;
	cout << "---------------------------------------" << endl;
	cout << "---------------------------------------" << endl;

	cout << "--- Writing to output file and saving ---" << endl;

	// outfile->Write();
	// outfile->Close();

	// TFile *inputfile = new TFile(Form("rootfiles/%i_dxdy", runnum), "READ");



	// TH1D *h_atime = new TH1D( "atime", Form("HCal ADC Time, All Channels - SBS = %i%%, %s, Run %i; ns", sbsfieldscale, run_target.Data(), runnum), 160, 0, 160 );
	// TH1D *h_E = new TH1D( "E", Form("HCal Cluster E, All Channels - SBS = %i%%, %s, Run %i; GeV", sbsfieldscale, run_target.Data(), runnum), 500, 0, 0.5 );
	// TH1D *h_E_cut = new TH1D( "E_cut", Form("HCal Cluster E All Cuts, All Channels - SBS = %i%%, %s, Run %i; GeV", sbsfieldscale, run_target.Data(), runnum), 100, 0, 0.5 );
	// TH1D *h_E_exp = new TH1D( "E_exp", Form("Expected Energy Dep in HCal - SBS = %i%%, %s, Run %i; GeV", sbsfieldscale, run_target.Data(), runnum), 100, 0, 0.2 );
	// TH1D *h_vert = new TH1D( "vert", Form("Vertex Position - SBS = %i%%, %s, Run %i; m", sbsfieldscale, run_target.Data(), runnum), 200, -1.0, 1.0 );
	// TH1D *hdx = new TH1D( "dx", Form("HCal dx - SBS = %i%%, %s, Run %i; m", sbsfieldscale, run_target.Data(), runnum), 200, -4.0, 2.0 );
	// TH1D *hdy = new TH1D( "dy", Form("HCal dy - SBS = %i%%, %s, Run %i; m", sbsfieldscale, run_target.Data(), runnum), 100, -1.2, 1.2 );
	// TH1D *hKE_p = new TH1D( "KE_p", Form("Scattered Proton Kinetic Energy - SBS = %i%%, %s, Run %i", sbsfieldscale, run_target.Data(), runnum), 500, 0.0, 5.0 );
	// TH1D *hKElow = new TH1D( "KElow", Form("Lowest Elastic E Sampled in HCal - SBS = %i%%, %s, Run %i; GeV", sbsfieldscale, run_target.Data(), runnum), 500, 0.0, 0.2 );
	// TH1D *hDiff = new TH1D( "hDiff",Form("HCal time - BBCal time - SBS = %i%%, %s, Run %i; ns", sbsfieldscale, run_target.Data(), runnum), 1300, -500, 800 );
	// TH1D *hX = new TH1D( "X", Form("HCal X - SBS = %i%%, %s, Run %i; m", sbsfieldscale, run_target.Data(), runnum), 100, -4.0, 2.0 );
	// TH1D *hY = new TH1D( "Y", Form("HCal Y - SBS = %i%%, %s, Run %i; m", sbsfieldscale, run_target.Data(), runnum), 100, -1.2, 1.2 );
	// TH2D *h2_CvCh = new TH2D( "CvCh", Form("HCal Coeff Single Block Clusters - SBS = %i%%, %s, Run %i; channel, GeV", sbsfieldscale, run_target.Data(), runnum), kNcell, 0, kNcell, 200, 0, 1.0 );
	// TH2D *h2_EvCh = new TH2D( "EvCh", Form("HCal Cluster E Single Block Clusters - SBS = %i%%, %s, Run %i; channel, GeV", sbsfieldscale, run_target.Data(), runnum), kNcell, 0, kNcell, 50, 0, 0.5 );
	// TH2D *h2_dxdy = new TH2D("dxdy",Form("Delta Plot (q-proj. dxdy) - SBS = %i%%, %s, Run %i;y_{HCAL}-y_{expect} (m); x_{HCAL}-x_{expect} (m)", sbsfieldscale, run_target.Data(), runnum), 50, -2, 2, 100, -4.0, 2.0 );
	// TH2D *h2_dxdy_all = new TH2D("dxdy_all",Form("Delta Plot (NO CUT) - SBS = %i%%, %s, Run %i;y_{HCAL}-y_{expect} (m); x_{HCAL}-x_{expect} (m)", sbsfieldscale, run_target.Data(), runnum),125,-2,2,125,-4,6);
	// TH2D *h2_dxVE = new TH2D("dxVE",Form("dxVE - SBS = %i%%, %s, Run %i;x_{HCAL}-x_{expect} (m); GeV", sbsfieldscale, run_target.Data(), runnum), 100, -4.0, 2.0, 100, 0, 0.5 );
	// TH2D *h2_rowcol = new TH2D( "hrowcol", Form("HCal Block Position Elastics, HCal - SBS = %i%%, %s, Run %i; Col; Row", sbsfieldscale, run_target.Data(), runnum), kNcols, 0, kNcols, kNrows, -kNrows, 0 );


	auto total_time_end = high_resolution_clock::now();
	auto total_time_duration = duration_cast<minutes>(total_time_end - total_time_start);
	cout << "Total time for analysis: " << total_time_duration.count() << " minutes. " << endl;

}

