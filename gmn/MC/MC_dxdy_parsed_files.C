#include <iostream>
#include <fstream>
#include <numeric>
#include <iomanip>
#include <math.h>
#include <algorithm>
#include <string>
#include <chrono>
#include <TF1.h>
#include "TStopwatch.h"
using namespace std;
using namespace std::chrono;

#include "/w/halla-scshelf2102/sbs/jboyd/include/include_files.h"
#include "/w/halla-scshelf2102/sbs/jboyd/include/GEM_lookups.h"
#include "/w/halla-scshelf2102/sbs/jboyd/include/beam_variables.h"

Double_t fit_gaus(Double_t * x, Double_t *par){

	Double_t g = 0.0;

	g = par[0]*exp((-0.5)*pow(((x[0] -  par[1])/par[2]),2));
	
	return g;
}

template<typename T>
double VectorMean(std::vector<T> const& v){
	if(v.empty()){
		return 0;
	}
	return std::accumulate(v.begin(), v.end(), 0.0)/v.size();
}

bool single_run = true;
bool multi_run = false;

bool calc_W = true;
bool use_heavy_cut = false;

bool correct_beam_energy = false;
bool fiducial_cut = false;

bool crosstalk = false;
TString XTALK_ONOFF = "XTALK_OFF";
int ratio_threshold = 4;

//Run info and lookups
TString run_target = "LD2";
int kine = 8;
int sbsfieldscale = 70;

int runnum = lookup_parsed_runnums(run_target.Data(), kine, sbsfieldscale, 0);

vector<int> runnum_vec;

TString experiment = "gmn";
int pass = 1;

//Experimental Lookup Parameters
double E_beam = lookup_beam_energy(runnum); //Electron beam energy (electron energy) in GeV
double SBS_field = lookup_run_info(runnum, "sbs_field"); //Strength (in percentage) of SBS magnet

double BB_dist, BB_theta, W_mean, W_sigma;
double dx_p, dx_p_sigma, dy_p, dy_p_sigma, dx_n, dx_n_sigma, dy_n, dy_n_sigma, dx_pn_max;

// TString rootfile_dir = "/lustre19/expphy/volatile/halla/sbs/jboyd/Rootfiles/xtalk/11449";
TString rootfile_dir;
// TString input_rootfile;

TFile *outfile;
TChain *TC = new TChain("T");
vector<TString> master_cut_vec;
TString master_cut_string;

TString elastic_yield_str = "";
TCut master_cut = "";

bool theta_pq_cut = true;
double theta_pq_p_thresh = 0.04;
double theta_pq_n_thresh = 0.04;

//Experimental Constants, Thresholds, cuts, etc DEFINITIONS
const double pi = TMath::Pi();
const double Mp = 0.938272; //Mass of proton [GeV]
const double Mn = 0.939565; //Mass of neutron [GeV]
const double Me = 0.00051; //Mass of electron [GeV]

//SBS Magnet
const Double_t Dgap = 48.0*2.54/100.0; //about 1.22 m
const Double_t maxSBSfield = 1.26; //Tesla
const Double_t SBSdist = 2.25; //m
const Double_t dipGap = 1.22; //m
const Double_t sbsmaxfield = 3.1 * atan( 0.85/(11.0 - 2.25 - 1.22/2.0 ))/0.3/1.22/0.7;

double W2_mean; //Invariant Mass-squared (mean val) {With perfect optics W2 = Mp. Can be calculated run-by-run}
double W2_sigma; //Invariant Mass-squared sigma {Reasonable default/guess. Can be calculated run-by-run from W plot}

//HCal constants and stuff
double tdiff = 510;		//Time difference between BBCal and HCal signals
double tdiff_max = 10;	//Maximum time difference from coincidences through tdctrig cut
double HCal_dist; 	//Distace from HCal face to target chamber
double HCal_theta;		//Theta angle for HCal from downstream beamline
double scint_intersect, x_expected_HCal, y_expected_HCal;
double hcal_y_fmin, hcal_y_fmax, hcal_x_fmin, hcal_x_fmax;

//Scattered kinematics
double e_prime_theta; //Scattered electron theta angle
double e_prime_phi; //Scattered electron phi angle
double p_el, nu, pp, nucleon_theta, nucleon_phi, E_ep, p_ep, Q2, W, E_pp, E_nucleon, KE_p, dx, dy;

//Static Detector Parameters
const int maxTracks = 1000; // Reasonable limit on tracks to be stored per event
const int maxTdcChan = 10; // Set to accomodate original 5 TDCTrig channels with buffer
const double hcalheight = -0.2897; // Height of HCal above beamline

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

//Static Target Parameters
const double l_tgt = 0.15; // Length of the target (m)
const double rho_tgt = 0.0723; // Density of target (g/cc)
const double rho_Al = 2.7; // Density of aluminum windows (g/cc)
const double cell_diameter = 1.6*2.54; //cm, right now this is a guess
const double Ztgt = 1.0;
const double Atgt = 1.0;
const double Mmol_tgt = 1.008; //g/mol

//For energy-loss correction to beam energy:
const double dEdx_tgt=0.00574; //According to NIST ESTAR, the collisional stopping power of hydrogen is about 5.74 MeV*cm2/g at 2 GeV energy
const double dEdx_Al = 0.0021; //According to NIST ESTAR, the collisional stopping power of Aluminum is about 2.1 MeV*cm2/g between 1-4 GeV
const double uwallthick_LH2 = 0.0145; //cm
const double dwallthick_LH2 = 0.015; //cm
const double cellthick_LH2 = 0.02; //cm, this is a guess;
const double Alshieldthick = 2.54/8.0; //= 1/8 inch * 2.54 cm/inch

double p_recon, nu_recon, E_loss, E_corr, theta_pq_n, theta_pq_p;

int useAlshield = 0;

//Declare vars
Double_t atime[kNcell], row[kNcell], col[kNcell], tdctime[kNcell], cblkid[kNcell], cblke[kNcell];
Double_t nblk, nclus, SH_nclus, PS_nclus, hcal_x, hcal_y, hcal_e;

Double_t par[3];

double bb_tr_p[maxTracks], bb_tr_px[maxTracks], bb_tr_py[maxTracks], bb_tr_pz[maxTracks];
double bb_tr_vx[maxTracks], bb_tr_vy[maxTracks], bb_tr_vz[maxTracks], bb_tr_chi2[maxTracks];
double bb_fp_x[maxTracks], bb_fp_y[maxTracks], bb_fp_th[maxTracks], bb_fp_ph[maxTracks];
double bb_tgt_x[maxTracks], bb_tgt_y[maxTracks], bb_tgt_th[maxTracks], bb_tgt_ph[maxTracks];
double bb_tr_n, bb_ps_x, bb_ps_y, bb_ps_e, bb_sh_x, bb_sh_y, bb_sh_e;

Double_t TDCT_id[kNtdc], TDCT_tdc[kNtdc], hodo_tmean[kNtdc]; 
Int_t TDCTndata;

Long64_t Nevents;

//INITIALIZE ALL HISTOGRAMS:
TH1D *h_atime, *h_W, *h_W2, *h_W2recon, *h_KE_p, *h_KE_low, *h_Diff, *h_X, *h_Y, *h_E_eloss;
TH1D *h_W_cut, *h_W_fcut, *h_vz_cut;
TH1D *h_Q2, *h_E_ep, *h_E_pp;
TH1D *h_dy, *h_dy_cut, *h_dy_wcut, *h_dx, *h_dx_cut, *h_dx_wcut, *h_dx_fcut;
 
TH2D *h_E_ecorr_vs_vert;
TH2D *h_dxdy, *h_dxdy_cut, *h_dxdy_wcut, *h_dxdy_ncut, *h_dxdy_pcut, *h_dxdy_fcut;
TH2D *h_xy, *h_xy_cut, *h_xy_fcut, *h_xy_cut_p, *h_xy_cut_n, *h_PAngleCorr_theta, *h_PAngleCorr_phi;

double n_integral, p_integral, n_center, n_sigma, p_center, p_sigma;
double p_Beam, E_loss_outgoing;
double Eloss, E_beam_final, p_corr;
int n_counts, p_counts, elastic_yield;

TString pq_cut_String = "";

void MC_dxdy_parsed_files(){
	auto total_time_start = high_resolution_clock::now();
	TStopwatch *StopWatch = new TStopwatch();

	gStyle->SetPalette(55);
	cout << "--------------------------------------" << endl;
	cout << "Analysis started. " << endl;
	cout << "--------------------------------------" << endl;
	for(int i = 0; i < lookup_parsed_runs_cnt(run_target.Data(), kine, sbsfieldscale); i++){
		runnum_vec.push_back(lookup_parsed_runnums(run_target.Data(), kine, sbsfieldscale, i));
	}

	if( theta_pq_cut ){
		pq_cut_String = "_pqCut";
	}

	if( !crosstalk ){
		outfile = new TFile(Form("rootfiles/%s_SBS%i_mag%i%s_dxdy_parsed.root", run_target.Data(), kine, sbsfieldscale, pq_cut_String.Data() ), "RECREATE");		
	}
	if( crosstalk ){
		outfile = new TFile(Form("rootfiles/%s_SBS%i_mag%i_%s%s_dxdy_parsed.root", run_target.Data(), kine, sbsfieldscale, XTALK_ONOFF.Data(), pq_cut_String.Data() ), "RECREATE");		
	}

	h_E_eloss = new TH1D("E_eloss", Form("Scattered Electron Energy Loss in Target - SBS = %i%%, %s, Run %i", sbsfieldscale, run_target.Data(), runnum), 500, 0.0, (0.1)*E_beam);
	h_E_ecorr_vs_vert = new TH2D("h_E_ecorr_vs_vert", Form("Corrected Beam Energy vs Vertex - SBS = %i%%, %s, Run %i; E_{e} (GeV); Z_{vertex} (m)", sbsfieldscale, run_target.Data(), runnum), 250, -0.125, 0.125, 500, 0, 0.001);
	h_Q2 = new TH1D("h_Q2", Form("Momentum Transfer Q^2 - SBS = %i%%, %s, Run %i; GeV", sbsfieldscale, run_target.Data(), runnum), 250, 0.5, 3.0);
	h_E_ep = new TH1D("h_E_ep", Form("Scattered Electron Energy - SBS = %i%%, %s, Run %i; GeV", sbsfieldscale, run_target.Data(), runnum), 500, 0.0, 1.5*E_beam);
	h_E_pp = new TH1D("h_E_pp", Form("Scattered Proton Energy - SBS = %i%%, %s, Run %i; GeV", sbsfieldscale, run_target.Data(), runnum), 500, 0.0, 1.5*E_beam);
	h_W = new TH1D("h_W", Form("Invariant Mass W - SBS = %i%%, %s, Run %i; GeV", sbsfieldscale, run_target.Data(), runnum), 300, 0.0, 3.0);
	h_W_cut = new TH1D("h_W_cut", Form("Invariant Mass W (Coin & Vert Cuts) - SBS = %i%%, %s, Run %i; GeV", sbsfieldscale, run_target.Data(), runnum), 300, 0.0, 3.0);
	h_W_fcut = new TH1D("h_W_fcut", Form("Invariant Mass W (Fiduc. Cuts) - SBS = %i%%, %s, Run %i; GeV", sbsfieldscale, run_target.Data(), runnum), 300, 0.0, 3.0);
	h_KE_p = new TH1D("h_KE_p", Form("Scattered Proton Kinetic Energy - SBS = %i%%, %s, Run %i; GeV", sbsfieldscale, run_target.Data(), runnum), 500, 0.0, 1.5*E_beam);

	h_dx = new TH1D("h_dx",Form("dx (NO CUTS) - SBS = %i%%, %s, Run %i; x_{HCal} - x_{exp} (m);", sbsfieldscale, run_target.Data(), runnum), 500, -2.5,2.5);
	h_dx_cut = new TH1D("h_dx_cut",Form("dx (Basic CUTS) - SBS = %i%%, %s, Run %i; x_{HCal} - x_{exp} (m);", sbsfieldscale, run_target.Data(), runnum), 500, -2.5,2.5);
	h_dx_wcut = new TH1D("h_dx_wcut",Form("dx (W cut) - SBS = %i%%, %s, Run %i; x_{HCal} - x_{exp} (m);", sbsfieldscale, run_target.Data(), runnum), 500, -2.5,2.5);
	h_dx_fcut = new TH1D("h_dx_fcut",Form("dx (f cut) - SBS = %i%%, %s, Run %i; x_{HCal} - x_{exp} (m);", sbsfieldscale, run_target.Data(), runnum), 500, -2.5,2.5);
	h_dy = new TH1D("h_dy",Form("dy (NO CUTS) - SBS = %i%%, %s, Run %i; y_{HCal} - y_{exp} (m);", sbsfieldscale, run_target.Data(), runnum), 500, -1.25,1.25);
	h_dy_cut = new TH1D("h_dy_cut",Form("dy (Basic Cuts) - SBS = %i%%, %s, Run %i; y_{HCal} - y_{exp} (m);", sbsfieldscale, run_target.Data(), runnum), 500, -1.25,1.25);  
	h_dy_wcut = new TH1D("h_dy_wcut",Form("dy (W Cuts) - SBS = %i%%, %s, Run %i; y_{HCal} - y_{exp} (m);", sbsfieldscale, run_target.Data(), runnum), 500, -1.25,1.25);  

	h_dxdy = new TH2D("h_dxdy", Form("Hadron Spot(s) on HCal (NO CUTS) - SBS = %i%%, %s, Run %i;y_{HCal}-y_{expect} (m); x_{HCal}-x_{expect} (m)", sbsfieldscale, run_target.Data(), runnum), 300, -1.5, 1.5, 500, -2.5, 2.5 );
	h_dxdy_wcut = new TH2D("h_dxdy_wcut", Form("Hadron Spot(s) on HCal (W cut) - SBS = %i%%, %s, Run %i;y_{HCal}-y_{expect} (m); x_{HCal}-x_{expect} (m)", sbsfieldscale, run_target.Data(), runnum), 300, -1.5, 1.5, 500, -2.5, 2.5 );
	h_dxdy_cut = new TH2D("h_dxdy_cut", Form("Hadron Spot(s) on HCal (Basic cuts) - SBS = %i%%, %s, Run %i;y_{HCal}-y_{expect} (m); x_{HCal}-x_{expect} (m)", sbsfieldscale, run_target.Data(), runnum), 300, -1.5, 1.5, 500, -2.5, 2.5 );
	h_dxdy_ncut = new TH2D("h_dxdy_ncut", Form("Hadron Spot(s) on HCal (n cut) - SBS = %i%%, %s, Run %i;y_{HCal}-y_{expect} (m); x_{HCal}-x_{expect} (m)", sbsfieldscale, run_target.Data(), runnum), 250, -1.25, 1.25, 250, -2.5, 2.5 );
	h_dxdy_pcut = new TH2D("h_dxdy_pcut", Form("Hadron Spot(s) on HCal (p cut) - SBS = %i%%, %s, Run %i;y_{HCal}-y_{expect} (m); x_{HCal}-x_{expect} (m)", sbsfieldscale, run_target.Data(), runnum), 250, -1.25, 1.25, 250, -2.5, 2.5 );
	h_dxdy_fcut = new TH2D("h_dxdy_fcut", Form("Hadron Spot(s) on HCal (f cut) - SBS = %i%%, %s, Run %i;y_{HCal}-y_{expect} (m); x_{HCal}-x_{expect} (m)", sbsfieldscale, run_target.Data(), runnum), 250, -1.25, 1.25, 250, -2.5, 2.5 );
	h_xy = new TH2D("h_xy",Form("HCal Hadron Spots (x, y) (NO CUTS) - SBS = %i%%, %s, Run %i;y_{HCal} (m); x_{HCal} (m)", sbsfieldscale, run_target.Data(), runnum),12,-0.9,0.9,24,-2.165,1.435);
	h_xy_cut = new TH2D("h_xy_cut", Form("HCal Hadron Spots (x, y) (BASIC CUTS) - SBS = %i%%, %s, Run %i;y_{HCal} (m); x_{HCal} (m)", sbsfieldscale, run_target.Data(), runnum),12,-0.9,0.9,24,-2.165,1.435);
	h_xy_fcut = new TH2D("h_xy_fcut", Form("HCal Hadron Spots (x, y) (Fiduc. CUTS) - SBS = %i%%, %s, Run %i;y_{HCal} (m); x_{HCal} (m)", sbsfieldscale, run_target.Data(), runnum),12,-0.9,0.9,24,-2.165,1.435);
	h_xy_cut_p = new TH2D("h_xy_cut_p", Form("HCal Hadron Spots (x, y) (p CUT) - SBS = %i%%, %s, Run %i;y_{HCal} (m); x_{HCal} (m)", sbsfieldscale, run_target.Data(), runnum),12,-0.9,0.9,24,-2.165,1.435);
	h_xy_cut_n = new TH2D("h_xy_cut_n", Form("HCal Hadron Spots (x, y) (n CUT) - SBS = %i%%, %s, Run %i;y_{HCal} (m); x_{HCal} (m)", sbsfieldscale, run_target.Data(), runnum),12,-0.9,0.9,24,-2.165,1.435);

	h_PAngleCorr_theta = new TH2D( "h_PAngCorr_theta",Form("BB theta vs HCal theta - SBS = %i%%, %s, Run %i", sbsfieldscale, run_target.Data(), runnum), 200, 0.55, 0.75, 300, 0.35, 0.65 );
	h_PAngleCorr_phi = new TH2D( "h_PAngCorr_phi",Form("BB phi vs HCal phi - SBS = %i%%, %s, Run %i", sbsfieldscale, run_target.Data(), runnum), 500, -0.4, 0.1, 500, 2.7, 3.2 );
	h_vz_cut = new TH1D("h_vz_cut",Form("BB phi vs HCal phi - SBS = %i%%, %s, Run %i; vertex z (m);", sbsfieldscale, run_target.Data(), runnum), 250,-0.125,0.125);

	BB_dist = lookup_BB_dist(runnum);
	BB_theta = (pi/180.0)*lookup_BB_angle(runnum);
	HCal_dist = lookup_HCal_dist(runnum);
	HCal_theta = (pi/180.0)*lookup_HCal_angle(runnum);
	W_mean = lookup_parsed_cut(run_target.Data(), kine, sbsfieldscale, "W_mean");
	W_sigma = lookup_parsed_cut(run_target.Data(), kine, sbsfieldscale, "W_sigma");

	dx_p = lookup_dxdy(runnum, "dx_p");
	dx_p_sigma = lookup_dxdy(runnum, "dx_p_sigma");
	dy_p = lookup_dxdy(runnum, "dy");
	dy_p_sigma = lookup_dxdy(runnum, "dy_sigma");
	dx_n = lookup_dxdy(runnum, "dx_n");
	dx_n_sigma = lookup_dxdy(runnum, "dx_n_sigma");
	dy_n = lookup_dxdy(runnum, "dy");
	dy_n_sigma = lookup_dxdy(runnum, "dy_sigma");
	dx_pn_max = lookup_dxdy(runnum, "dx_pn_max");

	cout << "Run parameters: " << endl;
	cout << "Run: " << runnum << endl;
	cout << "Kinematic: SBS" << kine << endl;
	cout << "Target: " << run_target.Data() << endl;
	cout << "Beam Energy: " << E_beam << endl;
	cout << "SBS Field: " << sbsfieldscale << "%" << endl;
	cout << "-----------------------------------" << endl;
	cout << "BB angle: " << (180/pi)*BB_theta << endl;
	cout << "SBS angle: " << lookup_SBS_angle(runnum) << endl;
	cout << "HCal angle: " << (180/pi)*HCal_theta << endl;
	cout << "HCal distance: " << HCal_dist << endl;
	cout << "-----------------------------------" << endl << endl;

	if( !crosstalk ){
		rootfile_dir = "/lustre/expphy/volatile/halla/sbs/jboyd/simulation/out_dir/MC_REPLAY_OUT_DIR";
	}

 	if( crosstalk ){
		rootfile_dir = "/lustre/expphy/volatile/halla/sbs/jboyd/simulation/out_dir/MC_REPLAY_OUT_DIR";
	}

	if( single_run ){
		cout << "Running in single run mode: " << runnum << endl;
		cout << "--------------------------------------" << endl;
		cout << "Adding files to TChain from: " << rootfile_dir.Data() << endl;
				TC->Add(Form("%s/replayed_gmn_sbs8_v2_job_1.root", rootfile_dir.Data()));
				// TC->Add(Form("%s/gmn_parsed_%s_SBS%i_mag%i.root", rootfile_dir.Data(), run_target.Data(), kine, sbsfieldscale));
		// TC->Add(Form("%s/gmn_parsed_fulltree_SBS%i_%s_mag%i.root", rootfile_dir.Data(), kine, run_target.Data(), sbsfieldscale));
		// TC->Add(Form("%s/*%i*.root", rootfile_dir.Data(), runnum));
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

	cout << "--------------------------------------" << endl;
	cout << "Setting up branches... ";
	// Declare root tree variables and set values to memory locations in root file
	// Switch them on
	TC->SetBranchStatus( "*", 0 );

	// HCal
	TC->SetBranchStatus( "sbs.hcal.x", 1 );
	TC->SetBranchStatus( "sbs.hcal.y", 1 );
	TC->SetBranchStatus( "sbs.hcal.e", 1 );
	TC->SetBranchStatus( "sbs.hcal.nclus", 1);

	// BB track
	TC->SetBranchStatus( "bb.tr.chi2", 1 );
	TC->SetBranchStatus( "bb.tr.n", 1 );
	TC->SetBranchStatus( "bb.tr.px", 1 );
	TC->SetBranchStatus( "bb.tr.py", 1 );
	TC->SetBranchStatus( "bb.tr.pz", 1 );    
	TC->SetBranchStatus( "bb.tr.p", 1 );
	TC->SetBranchStatus( "bb.tr.vx", 1 );
	TC->SetBranchStatus( "bb.tr.vy", 1 );
	TC->SetBranchStatus( "bb.tr.vz", 1 );
	TC->SetBranchStatus( "bb.tr.r_x", 1 );
	TC->SetBranchStatus( "bb.tr.r_y", 1 );
	TC->SetBranchStatus( "bb.tr.r_th", 1 );
	TC->SetBranchStatus( "bb.tr.r_ph", 1 );
	TC->SetBranchStatus( "bb.tr.tg_x", 1 );
	TC->SetBranchStatus( "bb.tr.tg_y", 1 );
	TC->SetBranchStatus( "bb.tr.tg_th", 1 );
	TC->SetBranchStatus( "bb.tr.tg_ph", 1 );
	TC->SetBranchStatus( "bb.gem.track.nhits", 1);

	// BBCal shower preshower
	TC->SetBranchStatus( "bb.ps.e", 1 );
	TC->SetBranchStatus( "bb.ps.x", 1 );
	TC->SetBranchStatus( "bb.ps.y", 1 );
	TC->SetBranchStatus( "bb.sh.e", 1 );
	TC->SetBranchStatus( "bb.sh.x", 1 );
	TC->SetBranchStatus( "bb.sh.y", 1 );
	TC->SetBranchStatus( "bb.sh.nclus", 1 );
	TC->SetBranchStatus( "bb.ps.nclus", 1 );

	// Trigger TDC
	TC->SetBranchStatus( "bb.tdctrig.tdc", 1 );
	TC->SetBranchStatus( "bb.tdctrig.tdcelemID", 1 );
	TC->SetBranchStatus( "Ndata.bb.tdctrig.tdcelemID", 1 );

// Set BRANCH ADDRESSES
	// HCal
	TC->SetBranchAddress( "sbs.hcal.x", &hcal_x );
	TC->SetBranchAddress( "sbs.hcal.y", &hcal_y );
	TC->SetBranchAddress( "sbs.hcal.e", &hcal_e );
	TC->SetBranchAddress( "sbs.hcal.nclus", &nclus );

	// BB track
	TC->SetBranchAddress( "bb.tr.chi2", bb_tr_chi2 );
	TC->SetBranchAddress( "bb.tr.n", &bb_tr_n );
	TC->SetBranchAddress( "bb.tr.px", bb_tr_px );
	TC->SetBranchAddress( "bb.tr.py", bb_tr_py );
	TC->SetBranchAddress( "bb.tr.pz", bb_tr_pz );
	TC->SetBranchAddress( "bb.tr.p", bb_tr_p );
	TC->SetBranchAddress( "bb.tr.vx", bb_tr_vx );
	TC->SetBranchAddress( "bb.tr.vy", bb_tr_vy );
	TC->SetBranchAddress( "bb.tr.vz", bb_tr_vz );
	TC->SetBranchAddress( "bb.tr.r_x", bb_fp_x );
	TC->SetBranchAddress( "bb.tr.r_y", bb_fp_y );
	TC->SetBranchAddress( "bb.tr.r_th", bb_fp_th );
	TC->SetBranchAddress( "bb.tr.r_ph", bb_fp_ph );
	TC->SetBranchAddress( "bb.tr.tg_x", bb_tgt_x );
	TC->SetBranchAddress( "bb.tr.tg_y", bb_tgt_y );
	TC->SetBranchAddress( "bb.tr.tg_th", bb_tgt_th );
	TC->SetBranchAddress( "bb.tr.tg_ph", bb_tgt_ph );

	// BBCal shower preshower
	TC->SetBranchAddress( "bb.ps.e", &bb_ps_e );
	TC->SetBranchAddress( "bb.ps.x", &bb_ps_x );
	TC->SetBranchAddress( "bb.ps.y", &bb_ps_y );
	TC->SetBranchAddress( "bb.sh.e", &bb_sh_e );
	TC->SetBranchAddress( "bb.sh.x", &bb_sh_x );
	TC->SetBranchAddress( "bb.sh.y", &bb_sh_y );
	TC->SetBranchAddress( "bb.sh.nclus", &SH_nclus );
	TC->SetBranchAddress( "bb.ps.nclus", &PS_nclus );

	// Trigger TDC
	TC->SetBranchAddress( "bb.tdctrig.tdcelemID", TDCT_id );
	TC->SetBranchAddress( "bb.tdctrig.tdc", TDCT_tdc );
	TC->SetBranchAddress( "Ndata.bb.tdctrig.tdcelemID", &TDCTndata );
	cout << " done. " << endl;
	cout << "--------------------------------------" << endl;

	p_Beam = E_beam/(1.0 + E_beam/Mp*(1.0 - cos(BB_theta)));

	E_loss_outgoing = cell_diameter/2.0/sin(BB_theta)*rho_tgt*dEdx_tgt; //Should be about 1 MeV
	if( useAlshield !=0 ) E_loss_outgoing += Alshieldthick*rho_Al*dEdx_Al;

//FIDUCIAL CUT:
	hcal_y_fmin = -0.75;
	hcal_y_fmax = 0.75;
	hcal_x_fmin = -2.015;
	hcal_x_fmax = 1.285;

	master_cut_vec = {
			"sbs.hcal.nclus>0",
			"bb.ps.nclus>0",
			"bb.sh.nclus>0",
			"abs(bb.tr.vz[0])<=0.075",
			"bb.gem.track.nhits[0]>3",
			"bb.tr.n==1",
			Form("bb.ps.e>%f", lookup_parsed_cut(run_target.Data(), kine, sbsfieldscale, "PS_min")),
			Form("((bb.sh.e+bb.ps.e)/(bb.tr.p[0]))>(%f)&&((bb.sh.e+bb.ps.e)/(bb.tr.p[0]))<(%f)", lookup_parsed_cut(run_target.Data(), kine, sbsfieldscale, "Ep_min"), lookup_parsed_cut(run_target.Data(), kine, sbsfieldscale, "Ep_max")),
			// Form("((abs(((bb.sh.e+bb.ps.e)/(bb.tr.p[0]))))-%f)<%f", lookup_parsed_cut(runnum, "Ep"), lookup_parsed_cut(runnum, "Ep_sigma")),
			Form("sbs.hcal.e>%f",lookup_parsed_cut(run_target.Data(), kine, sbsfieldscale, "HCal_min")),
			Form("(bb.sh.e+bb.ps.e)>%f", lookup_parsed_cut(run_target.Data(), kine, sbsfieldscale, "SH_PS_min"))
		};
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
	cout << "Applying Master Cut: " << endl;
	cout << master_cut << endl;

	TEventList *ev_list = new TEventList("ev_list", "Elastic Events List");
	TC->Draw(">>ev_list", master_cut);
	Long64_t Nevents = ev_list->GetN();

	cout << "--------------------------------------" << endl;
	cout << "Number of events to analyze: " << Nevents << endl;
	cout << "--------------------------------------" << endl;
	cout << "--------------------------------------" << endl;
	cout << "Starting analysis loop on events..... " << endl;

	elastic_yield = 0;

	int watch_cnt = 0;	
	int five_percent = int(0.05*Nevents);
	vector<double> time_for_five;
	double average_time = 0.0, time_remaining;
	StopWatch->Start();
	
	for(Long64_t nevent = 0; nevent < Nevents; nevent++){
		TC->GetEntry( ev_list->GetEntry( nevent ));

		if( nevent%five_percent == 0){		
			StopWatch->Stop();

			if( watch_cnt == 0){
				cout << "Evt: " << nevent <<"/" << Nevents << Form("(%.0f/100%%)", 100.0*double(1.0*nevent/Nevents)) << ". Elastic yield = " << elastic_yield << ". "<< endl;
			}

			if( watch_cnt > 0 ){
				time_for_five.push_back(StopWatch->RealTime());	
				average_time = VectorMean(time_for_five);
				cout << "average time for 5 = " << average_time << endl;
				time_remaining = average_time*( 1.0 - double(nevent)/double(Nevents));
				cout << "Evt: " << nevent <<"/" << Nevents << Form("(%.0f/100%%)", 100.0*double(1.0*nevent/Nevents)) << ". Elastic yield = " << elastic_yield << ". Time left: " << time_remaining << endl;
			}
			StopWatch->Reset();
			StopWatch->Continue();
		}


		double bbcal_time=0.0, hcal_time=0.0;

		for(int ihit=0; ihit<TDCTndata; ihit++){
			if(TDCT_id[ihit]==5){
				bbcal_time=TDCT_tdc[ihit];
			}
			if(TDCT_id[ihit]==0){
				hcal_time=TDCT_tdc[ihit];
			}
		}

		double diff = hcal_time - bbcal_time; 

		if( fabs(diff - tdiff)>tdiff_max ){
			continue;
		}

		//Sanity check
      	if( (int)bb_tr_n!=1 ){
      		cout << "**************************************************************" << endl;
      		cout << "--------------------------------------------------------------" << endl;
      		cout << endl << endl << "WARNING: Total tracks not as expected from global cut. Check globalcut for errors." << endl << endl;
      		cout << "--------------------------------------------------------------" << endl;
      		cout << "**************************************************************" << endl;
      	} 


      	if( correct_beam_energy ){
      		Eloss = (bb_tr_vz[0]+l_tgt/2.0) * rho_tgt * dEdx_tgt + uwallthick_LH2 * rho_Al * dEdx_Al; //approximately 3 MeV
      		h_E_eloss->Fill( Eloss );

      		E_beam_final = E_beam - Eloss;
      		h_E_ecorr_vs_vert->Fill( bb_tr_vz[0], E_beam_final);
      	}
      	if( !correct_beam_energy){
      		Eloss = (bb_tr_vz[0]+l_tgt/2.0) * rho_tgt * dEdx_tgt + uwallthick_LH2 * rho_Al * dEdx_Al; //approximately 3 MeV
      		h_E_eloss->Fill( Eloss );

      		E_beam_final = E_beam;
      		h_E_ecorr_vs_vert->Fill( bb_tr_vz[0], E_beam_final);
      	}

      	////Corrections
	    //Correct the beam energy with energy loss in target using vertex position
	    Double_t E_loss = (bb_tr_vz[0]+l_tgt/2.0) * rho_tgt * dEdx_tgt + uwallthick_LH2 * rho_Al * dEdx_Al; //approximately 3 MeV
	    Double_t E_corr = E_beam - Eloss;	

      	p_corr = bb_tr_p[0] - E_loss_outgoing; //Neglecting mass of e'


    //Proceed only if at least one track exists in BB arm - lowest chi2 track always first element
      	if( bb_tr_n > 1){
      		continue;
      	}
		
		p_Beam = E_beam/(1.0 + E_beam/Mp*(1.0 - cos(BB_theta)));
      	
      	e_prime_theta = acos( bb_tr_pz[0]/bb_tr_p[0] ); //Ucorrected track momenutm to reconstruct e' theta
      	e_prime_phi = atan2( bb_tr_py[0], bb_tr_px[0]);

      	TVector3 vertex( 0, 0, bb_tr_vz[0] ); // z location of vertex in hall coordinates
		TLorentzVector P_beam( 0, 0, E_beam_final, E_beam_final ); //Mass of e negligable
		TLorentzVector k_prime( bb_tr_px[0], bb_tr_py[0], bb_tr_pz[0], bb_tr_p[0] );
		TLorentzVector P_targ( 0, 0, 0, Mp );

		p_el = E_beam_final/( 1.0+E_beam_final/Mp*( 1.0-cos(e_prime_theta) ) );
		nu = E_beam_final - bb_tr_p[0];
		pp = sqrt( pow(nu,2)+2.*Mp*nu );
		nucleon_phi = e_prime_phi + pi; //assume coplanarity
		//double thetanucleon = acos( (E_corr - BBtr_p[0]*cos(etheta))/pp ); //use elastic constraint on nucleon kinematics
		nucleon_theta = acos( (E_beam_final - bb_tr_pz[0])/pp ); //use elastic constraint on nucleon kinematics

		TVector3 pNhat( sin(nucleon_theta)*cos(nucleon_phi), sin(nucleon_theta)*sin(nucleon_phi), cos(nucleon_theta) );

		//Define HCal coordinate system
		TVector3 HCAL_zaxis( sin(-HCal_theta ), 0, cos(-HCal_theta) );
		TVector3 HCAL_xaxis( 0, -1, 0 );
		TVector3 HCAL_yaxis = HCAL_zaxis.Cross(HCAL_xaxis).Unit();

		TVector3 HCAL_origin = HCal_dist*HCAL_zaxis + hcalheight*HCAL_xaxis;

		TVector3 HCAL_pos = HCAL_origin + (hcal_x*HCAL_xaxis) + (hcal_y*HCAL_yaxis);

		//Define intersection points for hadron vector
		scint_intersect = ( HCAL_origin - vertex ).Dot( HCAL_zaxis ) / (pNhat.Dot( HCAL_zaxis ) );
		TVector3 HCAL_intersect = vertex + scint_intersect * pNhat;

		//Define the expected position of hadron on HCal from BB track
		x_expected_HCal = (HCAL_intersect - HCAL_origin).Dot( HCAL_xaxis );
		y_expected_HCal = (HCAL_intersect - HCAL_origin).Dot( HCAL_yaxis );

//--------------------------------------------------------------
//Calculate theta pq variables
      	//Reconstructed momentum, corrected for mean loss exiting the target
		p_recon = bb_tr_p[0] + E_loss_outgoing; 

		TLorentzVector k_prime_recon(p_recon*bb_tr_px[0]/bb_tr_p[0], p_recon*bb_tr_py[0]/bb_tr_p[0], p_recon*bb_tr_pz[0]/bb_tr_p[0], p_recon);
		TLorentzVector q_recon = P_beam - k_prime_recon;
		TVector3 qvec_recon = q_recon.Vect();

	//Calculate q vector as beam momentum - scattered k
		TLorentzVector q = P_beam - k_prime;

	//Expected neutron direction
		TVector3 Neutron_Direction = (HCAL_pos - vertex).Unit();

	//Expected proton direction
		//Need to incorporate deflection due to SBS magnet
		double Bdl = sbsfieldscale*maxSBSfield*Dgap;
		double Proton_Deflection = tan( 0.3*Bdl/qvec_recon.Mag() )*(HCal_dist - (SBSdist + Dgap/2.0) ); 

		TVector3 Proton_Direction = (HCAL_pos + Proton_Deflection*HCAL_xaxis - vertex).Unit();

		theta_pq_n = acos( Neutron_Direction.Dot( qvec_recon.Unit() ) );
		theta_pq_p = acos( Proton_Direction.Dot( qvec_recon.Unit() ) );

		if( theta_pq_cut ){
			if( run_target == "LD2" ){
				if( theta_pq_p < theta_pq_p_thresh && theta_pq_n < theta_pq_n_thresh){
					continue;
				}
				// if( theta_pq_n > theta_pq_n_thresh ){
				// 	continue;
				// }
			}
		}

//--------------------------------------------------------------

		TLorentzVector P_gammaN = P_targ + q; //(-px, -py, ebeam - pz, Mp + ebeam - p)

		E_ep = sqrt( pow(Me,2) + pow(bb_tr_p[0],2) ); // Obtain the scattered electron energy
		h_E_ep->Fill( E_ep );

		p_ep = bb_tr_p[0];

		Q2 = 2*E_beam_final*E_ep*( 1-(bb_tr_pz[0]/p_ep) ); // Obtain Q2 from beam energy, outgoing electron energy, and momenta
		h_Q2->Fill( Q2 );

		//Get invariant mass transfer W from the four-momentum of the scattered nucleon
		W = P_gammaN.M(); 
		h_W->Fill( W );

		//Use the electron kinematics to predict the proton momedntum assuming elastic scattering on free proton at rest (will need to correct for fermi motion):
		E_pp = nu + Mp; // Get energy of the proton
		E_nucleon = sqrt(pow(pp,2)+pow(Mp,2)); // Check on E_pp, same
		h_E_pp->Fill( E_pp ); // Fill histogram

		KE_p = nu; // For elastics
		h_KE_p->Fill( KE_p );

		dx = hcal_x - x_expected_HCal;
		dy = hcal_y - y_expected_HCal;


		//Resolve the hadron spots without cuts
		h_dx->Fill( dx );
		h_dy->Fill( dy );
		h_dxdy->Fill( dy, dx );
		h_xy->Fill( hcal_y, hcal_x );

	// Coincidence timing cut and vertex cut to resolve W well
		if( fabs(diff - tdiff)<tdiff_max && fabs(bb_tr_vz[0])<=0.075 ){
			h_W_cut->Fill( W );
		} 

		// Preliminary HCal projections with single cut on W
		if( fabs(W - W_mean) < W_sigma ){
			h_dx_wcut->Fill( dx );
			h_dy_wcut->Fill ( dy );
			h_dxdy_wcut->Fill( dy, dx );
		}

		//Populate position histograms with cuts
		h_dxdy_cut->Fill( dy, dx );
		h_dx_cut->Fill( dx );
		h_dy_cut->Fill( dy );

		//Populate BB/HCal correlation histograms from elastics
		h_PAngleCorr_phi->Fill( e_prime_phi, nucleon_phi );
		h_PAngleCorr_theta->Fill( e_prime_theta, nucleon_theta );

		//Fill vertex position histogram for cut on tracks
    	h_vz_cut->Fill( bb_tr_vz[0] );

		//Check "elastic" events on center HCal for id with spot checks
		bool HCal_on = false;
		bool is_p = false;
		bool is_n = false;

	//FIDUCIAL Cut

		if( fiducial_cut ){
			if( hcal_y>hcal_y_fmin && hcal_y<hcal_y_fmax && hcal_x>hcal_x_fmin && hcal_x<hcal_x_fmax ){
				HCal_on = true;
			}
			if( pow( (hcal_x - x_expected_HCal - dx_p)/dx_p_sigma,2) + pow( (hcal_y - y_expected_HCal - dy_p)/dy_p_sigma,2) <= pow(2.5,2) ){
				is_p = true;
			}
			if( pow( (hcal_x - x_expected_HCal - dx_n)/dx_n_sigma,2) + pow( (hcal_y - y_expected_HCal - dy_n)/dy_n_sigma,2) <= pow(2.5,2) ){
				is_n = true;
			}

	//Fill respective histograms for these checks.
			if( HCal_on && is_n ) h_dxdy_ncut->Fill( dy, dx );
			if( HCal_on && is_p ) h_dxdy_pcut->Fill( dy, dx );

	//----------neutron
			if( HCal_on && is_n ){
				if( (hcal_x - dx_pn_max )>hcal_x_fmin ){
					h_dxdy_fcut->Fill( dy, dx );
					h_dx_fcut->Fill( dx );
					h_W_fcut->Fill( W );
					h_xy_fcut->Fill( hcal_y, hcal_x );
					h_xy_cut_n->Fill( hcal_y, hcal_x );
					elastic_yield++;
				}
			}
	//----------proton
			else if( HCal_on && is_p ){
				if( (hcal_x + dx_pn_max)<hcal_x_fmax ){
					h_dxdy_fcut->Fill( dy, dx );
					h_dx_fcut->Fill( dx );
					h_W_fcut->Fill( W );
					h_xy_fcut->Fill( hcal_y, hcal_x );
					h_xy_cut_p->Fill( hcal_y, hcal_x );
					elastic_yield++;
				}
			}
		}
	//END OF FIDUCIAL Cut
		//Still should count elastic yields if we got this far.....
		if( !fiducial_cut ){
			elastic_yield++;
		}
//end of events loop  
    }
	
	cout << "---------------------------------------" << endl;
	cout << "-----Finished going through events-----" << endl;
	cout << "---------------------------------------" << endl;

	outfile->Write();
	outfile->Close();

	cout << "------------------------------------------------------------------"<< endl;
	cout << "                       ANALYSIS FINISHED" << endl;
	cout << "------------------------------------------------------------------"<< endl;
	cout << "Run parameters: " << endl;
	cout << "Run: " << runnum << endl;
	cout << "Kinematic: SBS" << kine << endl;
	cout << "Target: " << run_target.Data() << endl;
	cout << "Beam Energy: " << E_beam << endl;
	cout << "SBS Field: " << sbsfieldscale << "%" << endl;
	cout << "-----------------------------------" << endl;
	cout << "BB angle: " << lookup_BB_angle(runnum) << endl;
	cout << "SBS angle: " << lookup_SBS_angle(runnum) << endl;
	cout << "HCal angle: " << HCal_theta << endl;
	cout << "HCal distance: " << HCal_dist << endl;
	cout << "-----------------------------------" << endl << endl;
	cout << "Elastic yield: " << elastic_yield << endl << endl;
	cout << "---------------------------------------" << endl << endl;	
	cout << "------------------------------------------------------------------"<< endl;
	cout << "Cut info:" << endl;
	cout << "Master cut: " << master_cut << endl << endl;
	cout << "Basic cuts applied: " << endl;
	cout << "Pre-shower & Shower clusters each > 0" << endl;
	cout << "Shower + Pre-Shower cut: "<< lookup_parsed_cut(run_target.Data(), kine, sbsfieldscale, "SH_PS_min") - lookup_parsed_cut(run_target.Data(), kine, sbsfieldscale, "SH_PS_max") << endl;
	cout << "Vertex <= 0.075" << endl;
	cout << "Number of GEM planes hit > 3" << endl;
	cout << "Number of tracks per event = 1" << endl;
	cout << "Pre-shower cut: " << lookup_parsed_cut(run_target.Data(), kine, sbsfieldscale, "PS_min") << endl;
	cout << "E/p: center = " << lookup_parsed_cut(run_target.Data(), sbsfieldscale, kine, "Ep_min") << ", sigma = " << lookup_parsed_cut(run_target.Data(), kine, sbsfieldscale, "Ep_max") << endl;
	cout << "HCal cluster energy: " << lookup_parsed_cut(run_target.Data(), kine, sbsfieldscale, "HCal_min");
	cout << "------------------------------------------------------------------"<< endl;
	cout << "Cut lookups: " << endl;
	cout << "PS_min: " << lookup_parsed_cut(run_target.Data(), kine, sbsfieldscale, "PS_min") << endl;
	cout << "SH_PS: " << lookup_parsed_cut(run_target.Data(), kine, sbsfieldscale, "SH_PS_min") << "; " << endl;
	cout << "HCal_clus_e: " << lookup_parsed_cut(run_target.Data(), kine, sbsfieldscale, "HCal_min") << endl;
	cout << "Ep: " << lookup_parsed_cut(run_target.Data(), kine, sbsfieldscale, "Ep_min") << "; " << lookup_parsed_cut(run_target.Data(), kine, sbsfieldscale, "Ep_max") << endl;
	cout << "W2: " << lookup_parsed_cut(run_target.Data(), kine, sbsfieldscale, "W2_min") << "; " << lookup_parsed_cut(run_target.Data(), kine, sbsfieldscale, "W2_max") << endl;
	cout << "------------------------------------------------------------------"<< endl << endl;
	cout << "Output file: " << outfile->GetName() << endl << endl;
	cout << "------------------------------------------------------------------"<< endl;
	
	auto total_time_end = high_resolution_clock::now();
	auto total_time_duration = duration_cast<minutes>(total_time_end - total_time_start);
	cout << "Total time for analysis: " << total_time_duration.count() << " minutes. " << endl;
}