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

//Run info and lookups
int runnum = 11493;
TString experiment = "gmn";
int pass = 0;

int kine = lookup_kinematic(runnum); //Kinematic setting: 4 for SBS-4, 9 for SBS-9, etc.
TString run_target = lookup_target(runnum); //Target --> "LH2" or "LD2"
double e_beam = lookup_beam_energy(runnum); //Electron beam energy (electron energy) in GeV
double SBS_field = lookup_run_info(runnum, "sbs_field"); //Strength (in percentage) of SBS magnet

//Experimental Constants, Thresholds, cuts, etc DEFINITIONS
const double pi = TMath::Pi();
const double Mp = 0.938; //Mass of proton [GeV]
const double Mn = 0.93956; //Mass of neutron [GeV]
const double W2_mean = 0.92; //Invariant Mass-squared (mean val) {With perfect optics W2 = Mp. Can be calculated run-by-run}
const double W2_sigma = 0.3; //Invariant Mass-squared sigma {Reasonable default/guess. Can be calculated run-by-run from W plot}

//HCal constants and stuff
double tdiff = 510;		//Time difference between BBCal and HCal signals
double tdiff_max = 20;	//Maximum time difference from coincidences through tdctrig cut
double HCal_dist = 17.0; 	//Distace from HCal face to target chamber
double Hcal_theta = 34.7*(TMath::Pi()/180);		//Theta angle for HCal from downstream beamline
const Double_t sampfrac = 0.077; 	//Estimate of the sampling fraction from MC

const Int_t Ncell = 288; 	//Total N HCal modules
const Int_t Nrows = 24;		//Total N HCal rows
const Int_t Ncols = 12;		//Total N HCal cols
const Int_t Ntrack = 100;	//Max number of tracks per event
const Int_t Ntdc = 1000;	//Max number of tdc signals per event
const Double_t X_i = -2.20;	//Distance from center of beam to top of HCal [m]
const Double_t X_f = 1.47;	//Distance from center of beam to bottom of HCal [m] 
const Double_t Y_i = -0.853;	//Distance from center of beam to "opposite-beam" side of HCal [m] (BB side??)
const Double_t Y_f = 0.853;	//Distance from center of beam to side "away" from BB??


//DEFINE CUTS
TString master_cut = "sbs.hcal.nclus>0 && sbs.hcal.e>0.025 && bb.ps.nclus>0 && ";

TString rootfile_dir = Form("/lustre19/expphy/volatile/halla/sbs/sbs-gmn/GMN_REPLAYS/pass%i/SBS%i/%s/rootfiles", pass, kine, run_target.Data());

TChain *TC = new TChain("T");

void plot_dxdy(){

	auto total_time_start = high_resolution_clock::now();

	cout << "Starting analysis..." << endl;

	//DECLARE VARIABLES
	Double_t atime[Ncell], row[Ncell], col[Ncell], tdctime[Ncell], cblkid[Ncell], cblke[Ncell];
	Double_t nblk, nclus, SHnclus, PSnclus, hcalx, hcaly, hcale;
	UInt_t TBits;

	Double_t BBtr_p[Ntrack], BBtr_px[Ntrack], BBtr_py[Ntrack], BBtr_pz[Ntrack], BBtr_vz[Ntrack];
	Double_t BBtr_n, BBps_x, BBps_y, BBps_e, BBsh_x, BBsh_y, BBsh_e;

	Double_t TDCT_id[Ntdc], TDCT_tdc[Ntdc], hodo_tmean[Ntdc];
	Int_t TDCTndata;

//SETUP Branches and leaves
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
	TC->SetBranchStatus( "fEvtHdr.fTrigBits", 1 );
	TC->SetBranchStatus( "bb.tr.vz", 1 );
	TC->SetBranchStatus( "bb.ps.e", 1 );
	TC->SetBranchStatus( "bb.ps.x", 1 );
	TC->SetBranchStatus( "bb.ps.y", 1 );
	TC->SetBranchStatus( "bb.sh.e", 1 );
	TC->SetBranchStatus( "bb.sh.x", 1 );
	TC->SetBranchStatus( "bb.sh.y", 1 );
	TC->SetBranchStatus( "g.trigbits", 1 );
	TC->SetBranchStatus( "bb.tdctrig.tdc", 1 );
	TC->SetBranchStatus( "bb.tdctrig.tdcelemID", 1 );
	TC->SetBranchStatus( "Ndata.bb.tdctrig.tdcelemID", 1 );

	TC->SetBranchAddress( "sbs.hcal.clus_blk.atime", atime );
	TC->SetBranchAddress( "sbs.hcal.clus_blk.row", row );
	TC->SetBranchAddress( "sbs.hcal.clus_blk.col", col );
	TC->SetBranchAddress( "sbs.hcal.clus_blk.tdctime", tdctime );
	TC->SetBranchAddress( "sbs.hcal.clus_blk.id", cblkid );
	TC->SetBranchAddress( "sbs.hcal.clus_blk.e", cblke );
	TC->SetBranchAddress( "sbs.hcal.x", &hcalx );
	TC->SetBranchAddress( "sbs.hcal.y", &hcaly );
	TC->SetBranchAddress( "sbs.hcal.e", &hcale );
	TC->SetBranchAddress( "sbs.hcal.nblk", &nblk );
	TC->SetBranchAddress( "sbs.hcal.nclus", &nclus );
	TC->SetBranchAddress( "bb.sh.nclus", &SHnclus );
	TC->SetBranchAddress( "bb.ps.nclus", &PSnclus );
	TC->SetBranchAddress( "bb.tr.n", &BBtr_n );
	TC->SetBranchAddress( "bb.tr.px", BBtr_px );
	TC->SetBranchAddress( "bb.tr.py", BBtr_py );
	TC->SetBranchAddress( "bb.tr.pz", BBtr_pz );
	TC->SetBranchAddress( "bb.tr.vz", BBtr_vz );
	TC->SetBranchAddress( "bb.tr.p", BBtr_p );
	TC->SetBranchAddress( "bb.ps.e", &BBps_e );
	TC->SetBranchAddress( "bb.ps.x", &BBps_x );
	TC->SetBranchAddress( "bb.ps.y", &BBps_y );
	TC->SetBranchAddress( "bb.sh.e", &BBsh_e );
	TC->SetBranchAddress( "bb.sh.x", &BBsh_x );
	TC->SetBranchAddress( "bb.sh.y", &BBsh_y );
	TC->SetBranchAddress( "fEvtHdr.fTrigBits", &TBits ); //For GEn, TBits==5 is coin
	TC->SetBranchAddress( "bb.tdctrig.tdcelemID", TDCT_id );
	TC->SetBranchAddress( "bb.tdctrig.tdc", TDCT_tdc );
	TC->SetBranchAddress( "Ndata.bb.tdctrig.tdcelemID", &TDCTndata );

	  // Declare outfile
	TFile *fout = new TFile( "deltaplot_.root", "RECREATE" );

	// Cut on global parameters defined in setup file
	TEventList *elist = new TEventList("elist","Elastic Event List");
	C->Draw(">>elist",globalcut);

	// Initialize histograms
	TH1D *h_atime = new TH1D( "atime", "HCal ADC Time, All Channels; ns", 160, 0, 160 );
	TH2D *h_CvCh = new TH2D( "CvCh", "HCal Coeff Single Block Clusters; channel, GeV", kNcell, 0, kNcell, 200, 0, 1.0 );
	TH1D *h_E = new TH1D( "E", "HCal Cluster E, All Channels; GeV", 500, 0, 0.5 );
	TH1D *h_E_cut = new TH1D( "E_cut", "HCal Cluster E All Cuts, All Channels; GeV", 100, 0, 0.5 );
	TH1D *h_E_exp = new TH1D( "E_exp", "Expected Energy Dep in HCal; GeV", 100, 0, 0.2 );
	TH1D *h_vert = new TH1D( "vert", "Vertex Position; m", 200, -1.0, 1.0 );
	TH2D *h_EvCh = new TH2D( "EvCh", "HCal Cluster E Single Block Clusters; channel, GeV", kNcell, 0, kNcell, 50, 0, 0.5 );
	TH1D *h_W2 = new TH1D( "W2", "W2 No Cuts; GeV", 200, 0.0, 4.0 );
	TH1D *h_W2recon = new TH1D( "W2recon", "W2 No Cuts Reconstructed; GeV", 200, 0.0, 4.0 );
	TH2D *hdxdy = new TH2D("dxdy",";y_{HCAL}-y_{expect} (m); x_{HCAL}-x_{expect} (m)", 50, -2, 2, 100, -4.0, 2.0 );
	TH2D *hdxdy_all = new TH2D("dxdy_all",";y_{HCAL}-y_{expect} (m); x_{HCAL}-x_{expect} (m)",125,-2,2,125,-4,6);
	TH1D *hdx = new TH1D( "dx", "HCal dx; m", 200, -4.0, 2.0 );
	TH1D *hdy = new TH1D( "dy", "HCal dy; m", 100, -1.2, 1.2 );
	TH1D *hKE_p = new TH1D( "KE_p", "Scattered Proton Kinetic Energy", 500, 0.0, 5.0 );
	TH2D *hdxVE = new TH2D("dxVE",";x_{HCAL}-x_{expect} (m); GeV", 100, -4.0, 2.0, 100, 0, 0.5 );
	TH1D *hKElow = new TH1D( "KElow", "Lowest Elastic E Sampled in HCal (GeV)", 500, 0.0, 0.2 );
	TH1D *hDiff = new TH1D( "hDiff","HCal time - BBCal time (ns)", 1300, -500, 800 );
	TH2D *hrowcol = new TH2D( "hrowcol", "HCal Block Position Elastics, HCal; Col; Row", kNcols, 0, kNcols, kNrows, -kNrows, 0 );
	TH1D *hX = new TH1D( "X", "HCal X; m", 100, -4.0, 2.0 );
	TH1D *hY = new TH1D( "Y", "HCal Y; m", 100, -1.2, 1.2 );

	auto total_time_end = high_resolution_clock::now();
  	auto total_time_duration = duration_cast<minutes>(total_time_end - total_time_start);
  	cout << "Total time for analysis: " << total_time_duration.count() << " minutes. " << endl;

}