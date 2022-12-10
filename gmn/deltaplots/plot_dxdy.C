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
int runnum = 11581;
TString experiment = "gmn";
int pass = 0;

int kine = lookup_kinematic(runnum); //Kinematic setting: 4 for SBS-4, 9 for SBS-9, etc.
TString run_target = lookup_target(runnum); //Target --> "LH2" or "LD2"
double E_beam = lookup_beam_energy(runnum); //Electron beam energy (electron energy) in GeV
double SBS_field = lookup_run_info(runnum, "sbs_field"); //Strength (in percentage) of SBS magnet

// TString rootfile_dir = "/lustre19/expphy/volatile/halla/sbs/jboyd/Rootfiles/xtalk/11449";
TString rootfile_dir = Form("/lustre19/expphy/volatile/halla/sbs/sbs-gmn/GMN_REPLAYS/pass0/SBS%i/%s/rootfiles", kine, run_target.Data());
TString input_rootfile = Form("%s/*%i*", rootfile_dir.Data(), runnum);

//Experimental Constants, Thresholds, cuts, etc DEFINITIONS
const double pi = TMath::Pi();
const double M_p = 0.938; //Mass of proton [GeV]
const double M_n = 0.93956; //Mass of neutron [GeV]
const double W2_mean = 0.916; //Invariant Mass-squared (mean val) {With perfect optics W2 = M_p. Can be calculated run-by-run}
const double W2_sigma = .2; //Invariant Mass-squared sigma {Reasonable default/guess. Can be calculated run-by-run from W plot}

//HCal constants and stuff
double tdiff = 510;		//Time difference between BBCal and HCal signals
double tdiff_max = 20;	//Maximum time difference from coincidences through tdctrig cut
double HCal_dist = 17.0; 	//Distace from HCal face to target chamber
double HCal_theta = 22.1*(TMath::Pi()/180);		//Theta angle for HCal from downstream beamline
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
TString master_cut = "sbs.hcal.nclus>0&&bb.ps.e>0.120&&bb.sh.nclus>0&&bb.ps.nclus>0&&sbs.hcal.e>0.05";
//&&sbs.hcal.e>0.025&&bb.tr.n==1
// TString master_cut = "bb.ps.e>0.15&&abs(bb.tr.vz)<0.075&&sbs.hcal.nclus>0&&bb.tr.n==1";

TH1D *h_atime;
TH2D *h_CvCh;
TH1D *h_E;
TH1D *h_E_cut;
TH1D *h_E_exp;
TH1D *h_vert;
TH2D *h_EvCh;
TH1D *h_W2;
TH1D *h_W2recon;
TH2D *h_dxdy;
TH2D *h_dxdy_all;
TH1D *h_dx;
TH1D *h_dy;
TH1D *h_KE_p;
TH2D *h_dxVE;
TH1D *h_KElow;
TH1D *h_Diff;
TH2D *h_rowcol;
TH1D *h_X;
TH1D *h_Y;

TChain *TC = new TChain("T");

void plot_dxdy(){

	auto total_time_start = high_resolution_clock::now();

	cout << "Starting analysis..." << endl;

	cout << "Adding files to TChain from: " << endl;
	cout << rootfile_dir.Data() << endl;

	TC->Add(input_rootfile.Data());

	//DECLARE VARIABLES
	Double_t atime[Ncell], row[Ncell], col[Ncell], tdctime[Ncell], cblkid[Ncell], cblke[Ncell];
	Double_t nblk, nclus, SHnclus, PSnclus, hcalx, hcaly, hcale;
	UInt_t TBits;

	Double_t BBtr_p[Ntrack], BBtr_px[Ntrack], BBtr_py[Ntrack], BBtr_pz[Ntrack], BBtr_vz[Ntrack];
	Double_t BBtr_n, BBps_x, BBps_y, BBps_e, BBsh_x, BBsh_y, BBsh_e;

	Double_t TDCT_id[Ntdc], TDCT_tdc[Ntdc], hodo_tmean[Ntdc];
	Int_t TDCTndata;

	cout << "Setting up branches" << endl;

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
	// TC->SetBranchStatus( "g.trigbits", 1 );
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

	cout << "Finished setting up branches. " << endl;

	// Declare outfile
	TFile *fout = new TFile( Form("%i_dxdy.root", runnum), "RECREATE" );

	// Master Cut on global parameters defined above
	TEventList *elas_evt_list = new TEventList("elas_evt_list","Elastic Event List");
	TC->Draw(">>elas_evt_list",master_cut.Data());

	cout << "Initializing histograms. " << endl;
	// Initialize histograms
	h_atime = new TH1D( "h_atime", "HCal ADC Time, All Channels; ns", 160, 0, 160 );
	h_CvCh = new TH2D( "h_CvCh", "HCal Coeff Single Block Clusters; channel, GeV", Ncell, 0, Ncell, 200, 0, 1.0 );
	h_E = new TH1D( "h_E", "HCal Cluster E, All Channels; GeV", 500, 0, 0.5 );
	h_E_cut = new TH1D( "h_E_cut", "HCal Cluster E All Cuts, All Channels; GeV", 100, 0, 0.5 );
	h_E_exp = new TH1D( "h_E_exp", "Expected Energy Dep in HCal; GeV", 100, 0, 0.2 );
	h_vert = new TH1D( "h_vert", "Vertex Position; m", 200, -1.0, 1.0 );
	h_EvCh = new TH2D( "h_EvCh", "HCal Cluster E Single Block Clusters; channel, GeV", Ncell, 0, Ncell, 50, 0, 0.5 );
	h_W2 = new TH1D( "h_W2", "W2 - No Cuts; GeV", 200, 0.0, 4.0 );
	h_W2recon = new TH1D( "h_W2recon", "W2 Recon - No Cuts ; GeV", 200, 0.0, 4.0 );
	h_dxdy = new TH2D("h_dxdy",";y_{HCAL}-y_{expect} (m); x_{HCAL}-x_{expect} (m)", 50, -8, 8, 100, -8.0, 8.0 );
	h_dxdy_all = new TH2D("h_dxdy_all",";y_{HCAL}-y_{expect} (m); x_{HCAL}-x_{expect} (m)",125,-2,2,125,-4,6);
	h_dx = new TH1D( "h_dx", "HCal dx; m", 200, -4.0, 2.0 );
	h_dy = new TH1D( "h_dy", "HCal dy; m", 100, -1.2, 1.2 );
	h_KE_p = new TH1D( "h_KE_p", "Scattered Proton Kinetic Energy", 500, 0.0, 5.0 );
	h_dxVE = new TH2D("h_dxVE",";x_{HCAL}-x_{expect} (m); GeV", 100, -4.0, 2.0, 100, 0, 0.5 );
	h_KElow = new TH1D( "h_KElow", "Lowest Elastic E Sampled in HCal (GeV)", 500, 0.0, 0.2 );
	h_Diff = new TH1D( "h_hDiff","HCal time - BBCal time (ns)", 1300, -500, 800 );
	h_rowcol = new TH2D( "h_hrowcol", "HCal Block Position Elastics, HCal; Col; Row", Ncols, 0, Ncols, Nrows, -Nrows, 0 );
	h_X = new TH1D( "h_X", "HCal X; m", 100, -4.0, 2.0 );
	h_Y = new TH1D( "h_Y", "HCal Y; m", 100, -1.2, 1.2 );

	cout << "Finished making histograms..." << endl;
//Variables for counting/managing loops, etc
	Long64_t Nevents = elas_evt_list->GetN();

	cout << endl << "Number of events from elastic event list: " << Nevents << endl << endl;
	Int_t elastic_yield = 0;

	for( Int_t row = 0; row < Nrows; row++){
		for( Int_t col = 0; col < Ncols; col++){
			h_rowcol->Fill( (col+1), -(row+1) );
		}
	}

//Variables for Beam, detectors, etc. with definition of HCal coordinate axis
	TLorentzVector P_beam( 0, 0, E_beam, E_beam );
	TLorentzVector P_targ( 0, 0, 0, 0.5*(M_p + M_n) ); //Average of proton and neutron rest mass
	TVector3 hcal_origin( -HCal_dist*sin(HCal_theta), 0, HCal_dist*cos(HCal_theta) );
	TVector3 hcal_zaxis = hcal_origin.Unit();
	TVector3 hcal_xaxis(0,-1,0);
	TVector3 hcal_yaxis = hcal_zaxis.Cross( hcal_xaxis ).Unit();

//Loop over events in elastic_events_list
	cout << "Starting loop over events. " << endl;
	for(Long64_t evt = 0; evt < Nevents; evt++){

		TC->GetEntry( elas_evt_list->GetEntry(evt) );
		if( evt%1000 ==0 ){
			cout << "Evt: " << evt << " of " << Nevents << " total events. Elastic yield: " << elastic_yield << ". " << endl;
			std::cout.flush();
		}

	//Define Lorentz Vectors for kinematics
		TLorentzVector k_prime( BBtr_px[0], BBtr_py[0], BBtr_pz[0], BBtr_p[0] );
		TLorentzVector q = P_beam - k_prime; //Interaction q-vector
		TVector3 vertex( 0, 0, BBtr_vz[0] ); //Vertex for scattering event at target
		TVector3 q_unit_vec = q.Vect().Unit(); //q unit vector 

		Double_t scint_intersect = (hcal_origin-vertex).Dot( hcal_zaxis )/q_unit_vec.Dot( hcal_zaxis ); //Location of scintillator coincidence
		TVector3 hcal_intersect = vertex + scint_intersect * q_unit_vec;
		Double_t x_expect = hcal_intersect.Dot( hcal_xaxis );
		Double_t y_expect = hcal_intersect.Dot( hcal_yaxis );

		Double_t W2_recon = (P_targ + q).M2();
		Double_t E_ep = BBtr_p[0]; // Obtain the scattered electron energy, neglect mass e
		Double_t p_ep = BBtr_p[0]; // Obtain the magnitude of scattered electron momentum
		Double_t Q2 = 2*E_beam*E_ep*( 1-(BBtr_pz[0]/p_ep) );
		Double_t nu = E_beam-E_ep; // Obtain energy transfer
		Double_t W2 = pow( M_p,2 )+2*M_p*nu-Q2; // Obtain W2 from Q2 and nu
		Double_t e_theta = acos( BBtr_pz[0]/BBtr_p[0] );
		Double_t e_phi = atan2( BBtr_py[0], BBtr_px[0] );
		Double_t nucleon_phi = e_phi + TMath::Pi(); //assume coplanarity
		Double_t nucleon_theta = acos( (E_beam - BBtr_pz[0])/p_ep ); //use elastic constraint on nucleon kinematics
		TVector3 pNhat( sin(nucleon_theta)*cos(nucleon_phi),sin(nucleon_theta)*sin(nucleon_phi),cos(nucleon_theta));
		Double_t KE_p = nu; //For elastics

		h_KE_p->Fill( KE_p ); 

//Cut on BBCal and HCal trigger coincidence
	    Double_t bbcal_time=0., hcal_time=0., coin_time=0., rf_time=0.;
	    bool cointrig = false;

	    for(Int_t ihit=0; ihit<TDCTndata; ihit++){
			if(TDCT_id[ihit]==5) bbcal_time=TDCT_tdc[ihit];
			if(TDCT_id[ihit]==0) hcal_time=TDCT_tdc[ihit];
			if(TDCT_id[ihit]==1) {
				coin_time=TDCT_tdc[ihit];
				cointrig=true;
	      	}
	      	if(TDCT_id[ihit]==4) rf_time=TDCT_tdc[ihit];
	    }

	    Double_t diff = hcal_time - bbcal_time; 
	    h_Diff->Fill( diff );

	    h_dxdy_all->Fill( hcaly - y_expect, hcalx - x_expect );
	    h_vert->Fill(BBtr_vz[0]);

	    h_W2->Fill(W2);
	    h_W2recon->Fill(W2_recon);

	    ///////////////
	    //PRIMARYCUTS//
	    //Cut on vertex inside target (15 cm --> +/- 7.5 cm)
	    if( abs(BBtr_vz[0])>0.075 ) continue;
	    //Cut on coin
	    // if( TBits==0 ) continue; //TBits 1 is BBCal trig, TBits 5 is coin
	    //Atime cut hcal
	    //if( atime[0]>60. && atime[0]<90. ) //Observed HCal ADC time peak

	    //Cut on  W2
	    if( abs(W2-W2_mean)>W2_sigma ) continue; //Observed W2 peak
	    
	    //Cut on BBCal HCal trig diff
	    if( abs(diff-tdiff)>tdiff_max ) continue; //BBCal/HCal trigger difference time
	    //Bigbite track / HCal angular correlation cut
	    //ENDCUTS//
	    ///////////
	    //elastic_yield++;

	    //Fill some histograms
	    h_KElow->Fill( KE_p*sampfrac );

	    h_E->Fill(hcale);
	    h_dxVE->Fill(hcalx - x_expect,hcale);

	    //Fill delta plots and others
	    h_dxdy->Fill( hcaly-y_expect, hcalx-x_expect );
	    h_dx->Fill( hcalx-x_expect );
	    h_dy->Fill( hcaly-y_expect );
	    h_atime->Fill( atime[0] );
	    h_Y->Fill( hcaly );

	    if( nblk==1 ){
	      h_EvCh->Fill( cblkid[0], hcale );
	    }
	    for( Int_t b=0; b<nblk; b++){      
	      //hrowcol->Fill( (Double_t)col[b]+1., -(Double_t)row[b]+1.);
	      h_rowcol->Fill( (Double_t)col[b], -(Double_t)row[b]-1); 
	    }


	    ///////////////
	    //SECONDARYCUTS

	    //Reject events where the primary block in the primary cluster is on the edge of the acceptance
	    if( row[0]==0 || row[0]==23 || col[0]==0 || col[0]==11 ){
	      continue;
	    }
	    elastic_yield++;

	    //ENDCUTS
	    /////////

	    Double_t E_exp = KE_p*sampfrac;

	    Double_t clusE = 0.0;
	    for( Int_t blk = 0; blk<(int)nblk; blk++ ){
	      Int_t blkid = int(cblkid[blk])-1; //-1 necessary since sbs.hcal.clus_blk.id ranges from 1 to 288 (different than just about everything else)
	      if( nblk==1 ) h_CvCh->Fill(blkid,cblke[blk]/E_exp);
	    }
	    
	    h_E_cut->Fill(clusE);	
	    h_E_exp->Fill(E_exp);

  	}
  	TCanvas *c_dxdy = new TCanvas("c_dxdy", "c_dxdy", 600, 500);
  	h_dxdy->Draw("colz");

  	TCanvas *c_W2recon = new TCanvas("c_W2recon", "c_W2recon", 600, 500);
  	h_W2recon->Draw();
    
    TCanvas *c_dy = new TCanvas("c_dy", "c_dy", 600, 500);
    h_dy->Draw();

  	cout << "Total calibration yield for run with current cuts: " << elastic_yield << "." << endl; 

  	fout->Write();
  	fout->Close();
  
  	cout << "Histograms populated and written to file." << endl;

	

	auto total_time_end = high_resolution_clock::now();
	auto total_time_duration = duration_cast<minutes>(total_time_end - total_time_start);
	cout << "Total time for analysis: " << total_time_duration.count() << " minutes. " << endl;

}