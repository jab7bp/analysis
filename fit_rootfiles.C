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

using namespace std::chrono;
#include "/w/halla-scshelf2102/sbs/jboyd/include/include_files.h"
#include "/w/halla-scshelf2102/sbs/jboyd/include/GEM_lookups.h"
#include "/w/halla-scshelf2102/sbs/jboyd/include/beam_variables.h"

//Run info and lookups
int runnum = 11449;
double SBS_field = lookup_run_info(runnum, "sbs_field"); //Strength (in percentage) of SBS magnet
int sbsfieldscale = int(100*lookup_run_info(runnum, "sbs_field"));

TString select_cut = "_wcut";

bool crosstalk = true;
TString XTALK_ONOFF = "XTALK_ON";
int ratio_threshold = 8;
bool overlay_crosstalk = true;

double pi = TMath::Pi();
Int_t num_par;

double dy_pn, dy_pn_sigma, dx_p, dx_p_sigma, dx_n, dx_n_sigma, dx_pn_max, BG_dx, BG_dy;
double dy_pn_xtalk_OFF, dy_pn_sigma_xtalk_OFF, dx_p_xtalk_OFF, dx_p_sigma_xtalk_OFF, dx_n_xtalk_OFF, dx_n_sigma_xtalk_OFF, dx_pn_max_xtalk_OFF, BG_dx_xtalk_OFF, BG_dy_xtalk_OFF;
double dy_pn_xtalk_ON, dy_pn_sigma_xtalk_ON, dx_p_xtalk_ON, dx_p_sigma_xtalk_ON, dx_n_xtalk_ON, dx_n_sigma_xtalk_ON, dx_pn_max_xtalk_ON, BG_dx_xtalk_ON, BG_dy_xtalk_ON;
int p_ellipse_cnt, n_ellipse_cnt;
double p_integral, n_integral, pn_integral, p_count_ellipse, n_count_ellipse, pn_count_ellipse;

double p_integral_xtalk_OFF, n_integral_xtalk_OFF, pn_integral_xtalk_OFF, p_count_ellipse_xtalk_OFF, n_count_ellipse_xtalk_OFF, pn_count_ellipse_xtalk_OFF;
double p_integral_xtalk_ON, n_integral_xtalk_ON, pn_integral_xtalk_ON, p_count_ellipse_xtalk_ON, n_count_ellipse_xtalk_ON, pn_count_ellipse_xtalk_ON;

int p_ellipse_cnt_xtalk_OFF, n_ellipse_cnt_xtalk_OFF;
int p_ellipse_cnt_xtalk_ON, n_ellipse_cnt_xtalk_ON;

vector<double> dxdy_vec;
TCutG *tcg_p, *tcg_n;
TCutG *tcg_p_xtalk_OFF, *tcg_n_xtalk_OFF;
TCutG *tcg_p_xtalk_ON, *tcg_n_xtalk_ON;

Double_t fit_gaus(Double_t * x, Double_t *par){

	Double_t g = 0.0;

	g = par[0]*exp((-0.5)*pow(((x[0] -  par[1])/par[2]),2));
	
	return g;
}

Double_t fit_dy(Double_t *x, Double_t *par){
	double total_fit = 0.0, dy_gaus = 0.0, BG_dy_gaus = 0.0;

	dy_gaus = par[0]*exp((-0.5)*pow(((x[0] -  par[1])/par[2]),2));
	BG_dy_gaus = par[3]*exp((-0.5)*pow(((x[0] -  par[4])/par[5]),2));

	total_fit = dy_gaus + BG_dy_gaus;

	return total_fit;
}

Double_t fit_dxdy(Double_t * x, Double_t *par){
	double total_fit = 0.0, p_gaus = 0.0, n_gaus = 0.0, BG_gaus = 0.0, pn_gaus = 0.0;
	
	if( sbsfieldscale == 0 ){
		pn_gaus = par[0]*exp((-0.5)*pow(((x[0] -  par[1])/par[2]),2));
		BG_gaus = par[3]*exp((-0.5)*pow(((x[0] -  par[4])/par[5]),2));

		total_fit = pn_gaus + BG_gaus;
	}

	else{
		p_gaus = par[0]*exp((-0.5)*pow(((x[0] -  par[1])/par[2]),2));
		n_gaus = par[3]*exp((-0.5)*pow(((x[0] -  par[4])/par[5]),2));
		BG_gaus = par[6]*exp((-0.5)*pow(((x[0] -  par[7])/par[8]),2));

		total_fit = p_gaus + n_gaus + BG_gaus;
	}

	return total_fit;
}

Double_t plot_dxdy_fit(Double_t *par){
	Double_t fit_x[1800], gaus1_y[1800], gaus2_y[1800], gaus3_y[1800];
	Int_t fit_n = 500;
	int fit_cnt = 0;

	if( sbsfieldscale == 0 ){
		for(double i = -3; i < 3; i+=.01){
			fit_x[fit_cnt] = i;
			gaus1_y[fit_cnt] = par[0]*exp((-0.5)*pow(((i -  par[1])/par[2]),2));
			gaus2_y[fit_cnt] = par[3]*exp((-0.5)*pow(((i -  par[4])/par[5]),2));
			fit_cnt++;
		}

		TGraph *gr_gaus1 = new TGraph(fit_n, fit_x, gaus1_y);
		gr_gaus1->SetLineColor(3);
		TGraph *gr_gaus2 = new TGraph(fit_n, fit_x, gaus2_y);
		gr_gaus2->SetLineColor(7);

		gr_gaus1->Draw("same");
		gr_gaus2->Draw("same");
	}
	else{
		for(double i = -3; i < 3; i+=.01){
			fit_x[fit_cnt] = i;
			gaus1_y[fit_cnt] = par[0]*exp((-0.5)*pow(((i -  par[1])/par[2]),2));
			gaus2_y[fit_cnt] = par[3]*exp((-0.5)*pow(((i -  par[4])/par[5]),2));
			gaus3_y[fit_cnt] = par[6]*exp((-0.5)*pow(((i -  par[7])/par[8]),2));
			fit_cnt++;
		}

		TGraph *gr_gaus1 = new TGraph(fit_n, fit_x, gaus1_y);
		gr_gaus1->SetLineColor(3);
		TGraph *gr_gaus2 = new TGraph(fit_n, fit_x, gaus2_y);
		gr_gaus2->SetLineColor(6);
		TGraph *gr_gaus3 = new TGraph(fit_n, fit_x, gaus3_y);
		gr_gaus3->SetLineColor(7);

		gr_gaus1->Draw("same");
		gr_gaus2->Draw("same");
		gr_gaus3->Draw("same");
	}
	return 1;

}
Double_t plot_dy_fit(Double_t *par){
	Double_t fit_x[1800], gaus1_y[1800], gaus2_y[1800];
	Int_t fit_n = 500;
	int fit_cnt = 0;

	for(double i = -3; i < 3; i+=.01){
		fit_x[fit_cnt] = i;
		gaus1_y[fit_cnt] = par[0]*exp((-0.5)*pow(((i -  par[1])/par[2]),2));
		gaus2_y[fit_cnt] = par[3]*exp((-0.5)*pow(((i -  par[4])/par[5]),2));
		fit_cnt++;
	}

	TGraph *gr_gaus1 = new TGraph(fit_n, fit_x, gaus1_y);
	gr_gaus1->SetLineColor(3);
	TGraph *gr_gaus2 = new TGraph(fit_n, fit_x, gaus2_y);
	gr_gaus2->SetLineColor(7);

	gr_gaus1->Draw("same");
	gr_gaus2->Draw("same");

	
	return 1;
}

Int_t cnt_ellipse_steps(double h, double k, double a, double b){
	Int_t cnt = 0;
	double y = 0.0;
	for(double x = -abs(h); x <= abs(h); x+=0.0001){
		y = k + (b)*sqrt( 1 - ((pow( (x - h), 2))/(a*a)) );
		cnt++;
	}
	
	for(double x = -abs(h); x <= abs(h); x+=0.0001){
		y = k - (b)*sqrt( 1 - ( (pow( (x - h), 2))/(a*a)) );
		cnt++;
	}
	return cnt;
}

Int_t set_tcg_ellipse(double h, double k, double a, double b, TCutG *tcg){
	Int_t cnt = 0;
	double y = 0.0;
	for(double x = -abs(h); x <= abs(h); x+=0.0001){
		y = k + (b)*sqrt( 1 - ((pow( (x - h), 2))/(a*a)) );
		tcg->SetPoint(cnt, x, y);
		cnt++;
	}
	
	for(double x = -abs(h); x <= abs(h); x+=0.0001){
		y = k - (b)*sqrt( 1 - ( (pow( (x - h), 2))/(a*a)) );
		tcg->SetPoint(cnt, x, y);
		cnt++;
	}
	return cnt;
}

TString select_file = "dxdy";

TString rootfile_dir = Form("/w/halla-scshelf2102/sbs/jboyd/analysis/gmn/deltaplots/rootfiles");

TFile *infile, *outfile, *crosstalk_infile;
TFile *infile_xtalk_OFF, *infile_xtalk_ON;
TString infile_dxdy, infile_dxdy_xtalk_OFF, infile_dxdy_xtalk_ON;
TString infile_calib;
TString infile_calib_W;
TEllipse *p_ellipse, *n_ellipse, *pn_ellipse;
TEllipse *p_ellipse_xtalk_OFF, *n_ellipse_xtalk_OFF, *pn_ellipse_xtalk_OFF;
TEllipse *p_ellipse_xtalk_ON, *n_ellipse_xtalk_ON, *pn_ellipse_xtalk_ON;

TH1D *hin_dx, *hin_dx_cut, *hin_dx_wcut, *hin_dy, *hin_dy_cut, *hin_dy_wcut;
TH1D *hin_dx_xtalk_OFF, *hin_dx_xtalk_ON, *hin_dy_xtalk_OFF, *hin_dy_xtalk_ON;
TH1D *hin_xtalk_OFF_dx, *hin_xtalk_ON_dx, h_xtalk_dx_diff;
TH2D *hin_dxdy, *hin_dxdy_cut, *hin_dxdy_wcut, *hin_dxdy_xtalk_OFF, *hin_dxdy_xtalk_ON;

Double_t par[9];

void fit_rootfiles(){
	auto total_time_start = high_resolution_clock::now();

	cout << "----------------------------" << endl;
	cout << "     Analysis Started" << endl;
	cout << "----------------------------" << endl;


	if( !crosstalk ){
		// infile_dxdy = Form("%i_dxdy.root", runnum);
		infile_dxdy = "LD2_SBS8_mag70_dxdy_parsed.root";
		infile_calib = Form("%i_prime_calibration.root", runnum);
		infile_calib_W = Form("%i_full_calibration.root", runnum);
	}
	if( crosstalk ){
		infile_dxdy = Form("%i_%s_dxdy.root", runnum, XTALK_ONOFF.Data());
		infile_calib = Form("%i_%s_prime_calibration.root", runnum, XTALK_ONOFF.Data());
		infile_calib_W = Form("%i_%s_full_calibration.root", runnum, XTALK_ONOFF.Data());
	}	
	if( crosstalk && overlay_crosstalk ){
		infile_dxdy_xtalk_OFF = Form("%i_XTALK_OFF_dxdy.root", runnum);
		infile_dxdy_xtalk_ON = Form("%i_XTALK_ON_dxdy.root", runnum);
	}

	if( select_file=="dxdy" ){
		infile = new TFile(Form("%s/%s", rootfile_dir.Data(), infile_dxdy.Data()), "READ");
	}
	if( select_file=="dxdy" && overlay_crosstalk){
		infile_xtalk_OFF = new TFile(Form("%s/%s", rootfile_dir.Data(), infile_dxdy_xtalk_OFF.Data()), "READ");
		infile_xtalk_ON = new TFile(Form("%s/%s", rootfile_dir.Data(), infile_dxdy_xtalk_ON.Data()), "READ");
	}
	if( select_file=="calib" ){
		infile = new TFile(Form("%s/%s", rootfile_dir.Data(), infile_calib.Data()), "READ");
	}
	if( select_file=="calib_W" ){
		infile = new TFile(Form("%s/%s", rootfile_dir.Data(), infile_calib_W.Data()), "READ");
	}

	cout << "Run Parameters: " << endl;
	cout << "Run: " << runnum << endl;
	cout << "SBS Field: " << sbsfieldscale << "%" << endl;
	cout << "----------------------------" << endl;
	cout << "Crosstalk ON?    " << crosstalk << endl;
	cout << "XTALK On/Off? ---- " << XTALK_ONOFF.Data() << endl;
	cout << "Ratio threshold: " << ratio_threshold << endl;
	cout << "----------------------------" << endl;
	cout << "Selected file: " << select_file << endl;
	cout << "Input file: " << infile->GetName() << endl;
	cout << "----------------------------" << endl;

	hin_dxdy = static_cast<TH2D*>(infile->Get(Form("h_dxdy%s", select_cut.Data())));

	if( crosstalk && overlay_crosstalk ){
		hin_dxdy_xtalk_OFF = static_cast<TH2D*>(infile_xtalk_OFF->Get(Form("h_dxdy%s", select_cut.Data())));
		hin_dx_xtalk_OFF = static_cast<TH1D*>(infile_xtalk_OFF->Get(Form("h_dx%s", select_cut.Data())));
		hin_dy_xtalk_OFF = static_cast<TH1D*>(infile_xtalk_OFF->Get(Form("h_dy%s", select_cut.Data())));

		hin_dxdy_xtalk_ON = static_cast<TH2D*>(infile_xtalk_ON->Get(Form("h_dxdy%s", select_cut.Data())));
		hin_dx_xtalk_ON = static_cast<TH1D*>(infile_xtalk_ON->Get(Form("h_dx%s", select_cut.Data())));
		hin_dy_xtalk_ON = static_cast<TH1D*>(infile_xtalk_ON->Get(Form("h_dy%s", select_cut.Data())));

	}
	
	// hin_dxdy_cut = static_cast<TH2D*>(infile->Get("h_dxdy_cut"));
	// hin_dxdy_wcut = static_cast<TH2D*>(infile->Get("h_dxdy_wcut"));

	hin_dx = static_cast<TH1D*>(infile->Get(Form("h_dx%s", select_cut.Data())));
	// hin_dx_cut = static_cast<TH1D*>(infile->Get("h_dx_cut"));
	// hin_dx_wcut = static_cast<TH1D*>(infile->Get("h_dx_wcut"));

	hin_dy = static_cast<TH1D*>(infile->Get(Form("h_dy%s", select_cut.Data())));
	// hin_dy_cut = static_cast<TH1D*>(infile->Get("h_dy_cut"));
	// hin_dy_wcut = static_cast<TH1D*>(infile->Get("h_dy_wcut"));



//--------------------------------------------------------------
//-----------------dx-------------------------------
	cout << "Fitting dx" << endl;
	TCanvas *c_dx = new TCanvas("c_dx", "c_dx", 600, 500);
	hin_dx->Draw();

	if( SBS_field != 0 ){num_par = 9;}
	if( SBS_field == 0 ){num_par = 6;}
	TF1 *fit_dx = new TF1("fit_dx", fit_dxdy, -3, 3, 9);

	if( sbsfieldscale == 0 ){
		cout << "sbsfield scale = " << sbsfieldscale << ". Num params: " << num_par << endl;
		fit_dx->SetParName(0, "pn - Norm");
		fit_dx->SetParName(1, "pn - Center");
		fit_dx->SetParName(2, "pn - Sigma");
		fit_dx->SetParName(3, "BG_pn - Norm");
		fit_dx->SetParName(4, "BG_pn - Center");
		fit_dx->SetParName(5, "BG_pn - Sigma");
		
		fit_dx->SetParLimits(0, 0, hin_dx->GetMaximum());
		fit_dx->SetParLimits(1, -1, 0);
		fit_dx->SetParLimits(2, 0.05, 0.3);
		fit_dx->SetParLimits(3, 0, (0.25)*hin_dx->GetMaximum());
		fit_dx->SetParLimits(4, -0.4, 0.2);
		fit_dx->SetParLimits(5, 1, 5);
	}
	else{
		cout << "sbsfield scale = " << sbsfieldscale << ". Num params: " << num_par << endl;
		fit_dx->SetParName(0, "p - Norm");
		fit_dx->SetParName(1, "p - Center");
		fit_dx->SetParName(2, "p - Sigma");
		fit_dx->SetParName(3, "n - Norm");
		fit_dx->SetParName(4, "n - Center");
		fit_dx->SetParName(5, "n - Sigma");
		fit_dx->SetParName(6, "BG_pn - Norm");
		fit_dx->SetParName(7, "BG_pn - Center");
		fit_dx->SetParName(8, "BG_pn - Sigma");
		
		fit_dx->SetParLimits(0, 0, hin_dx->GetMaximum());
		fit_dx->SetParLimits(1, -0.8, -0.2);
		fit_dx->SetParLimits(2, 0.05, 0.3);
		fit_dx->SetParLimits(3, (0.10)*hin_dx->GetMaximum(), (0.5)*hin_dx->GetMaximum());
		fit_dx->SetParLimits(4, 0.0, 0.6);
		fit_dx->SetParLimits(5, 0.1, 0.3);
		fit_dx->SetParLimits(6, 1, (0.05)*hin_dx->GetMaximum());
		fit_dx->SetParLimits(7, -0.4, 0.2);
		fit_dx->SetParLimits(8, 0.5, 5);
	}
	for(int i=0; i < 10; i++){
		hin_dx->Fit("fit_dx", "R+");
	
		fit_dx->GetParameters(&par[0]);
		fit_dx->SetParameters(par);
	}

	if( sbsfieldscale == 0 ){
		dx_p = fit_dx->GetParameter(1);
		dx_p_sigma = fit_dx->GetParameter(2);
		dx_n = dx_p;
		dx_n_sigma = dx_p_sigma;
		dx_pn_max = 0.05;
		BG_dx = fit_dx->GetParameter(3);
	}
	else{
		dx_p = fit_dx->GetParameter(1);
		cout << "-------------------------- dx_p: " << fit_dx->GetParameter(1) << endl;
		dx_p_sigma = fit_dx->GetParameter(2);
		dx_n = fit_dx->GetParameter(4);
		dx_n_sigma = fit_dx->GetParameter(5);
		dx_pn_max = (dx_n + dx_n_sigma) - (dx_p - dx_p_sigma);
		BG_dx = fit_dx->GetParameter(6);
	//-------------------------------------------------------------------------
		p_integral = (fit_dx->GetParameter(0))*(fit_dx->GetParameter(2))*sqrt(2*pi)/hin_dx->GetBinWidth(1);
		n_integral = (fit_dx->GetParameter(3))*(fit_dx->GetParameter(5))*sqrt(2*pi)/hin_dx->GetBinWidth(1);
	}

	hin_dx->Fit("fit_dx", "R+");
	fit_dx->Draw("same");
	plot_dxdy_fit(par);

//--------------------------------------------------------------
//-----------------dy-------------------------------
	cout << "Fitting dy" << endl;
	TCanvas *c_dy = new TCanvas("c_dy", "c_dy", 600, 500);
	hin_dy->Draw();

	num_par = 6;
cout << "--------------------------------- " << endl;
		cout << "373" << endl;
	TF1 *tf_fit_dy = new TF1("tf_fit_dy", fit_dy, -3, 3, 6);
	tf_fit_dy->SetParName(0, "dy - Norm");
	tf_fit_dy->SetParName(1, "dy - Center");
	tf_fit_dy->SetParName(2, "dy - Sigma");
	tf_fit_dy->SetParName(3, "BG_dy - Norm");
	tf_fit_dy->SetParName(4, "BG_dy - Center");
	tf_fit_dy->SetParName(5, "BG_dy - Sigma");

	tf_fit_dy->SetParLimits(0, 0, hin_dy->GetMaximum());
	tf_fit_dy->SetParLimits(1, -0.4, 0.2);
	tf_fit_dy->SetParLimits(2, 0.05, 0.3);
	tf_fit_dy->SetParLimits(3, 1, BG_dx);
	tf_fit_dy->SetParLimits(4, -0.2, 0.2);
	tf_fit_dy->SetParLimits(5, 0.5, 5);

	hin_dy->Fit("tf_fit_dy", "R+");
	tf_fit_dy->Draw("same");

	tf_fit_dy->GetParameters(&par[0]);
	tf_fit_dy->SetParameters(par);

	dy_pn = tf_fit_dy->GetParameter(1);
	dy_pn_sigma = tf_fit_dy->GetParameter(2);
	BG_dy = tf_fit_dy->GetParameter(3);
	plot_dy_fit(par);

//--------------------------------------------------------------
//-----------------dxdy---------------------------------
	cout << "Fitting dxdy" << endl;
	TCanvas *c_dxdy = new TCanvas("c_dxdy", "c_dxdy", 600, 500);
	hin_dxdy->Draw("colz");

//--------------------------------------------------------------
//-----------------ELLIPSE---------------------------------
	
	//ellipse constants...
	double p_ellipse_a = 1.5*dy_pn_sigma;
	double p_ellipse_b = 1.5*dx_p_sigma;
	double n_ellipse_a = 1.25*dy_pn_sigma;
	double n_ellipse_b = 1.25*dx_n_sigma;

	p_ellipse = new TEllipse(dy_pn, dx_p, p_ellipse_a, p_ellipse_b);
	p_ellipse->SetLineColor(2);
	p_ellipse->SetFillStyle(4005);
	p_ellipse->SetLineWidth(3);
	p_ellipse->Draw("same");

	n_ellipse = new TEllipse(dy_pn, dx_n, n_ellipse_a, n_ellipse_b);
	n_ellipse->SetLineColor(6);
	n_ellipse->SetFillStyle(4005);
	n_ellipse->SetLineWidth(3);
	n_ellipse->Draw("same");

	//count steps...
	p_ellipse_cnt = cnt_ellipse_steps(dy_pn, dx_p, p_ellipse_a, p_ellipse_b);
	n_ellipse_cnt = cnt_ellipse_steps(dy_pn, dx_n, n_ellipse_a, n_ellipse_b);

	tcg_p = new TCutG("tcg_p", p_ellipse_cnt);
	tcg_p->SetVarX("x");
	tcg_p->SetVarY("y");
	set_tcg_ellipse(dy_pn, dx_p, p_ellipse_a, p_ellipse_b, tcg_p);
	p_ellipse_cnt = tcg_p->IntegralHist(hin_dxdy)/hin_dx->GetBinWidth(1);

	tcg_n = new TCutG("tcg_n", n_ellipse_cnt);
	tcg_n->SetVarX("x");
	tcg_n->SetVarY("y");
	set_tcg_ellipse(dy_pn, dx_n, n_ellipse_a, n_ellipse_b, tcg_n);
	n_ellipse_cnt = tcg_n->IntegralHist(hin_dxdy)/hin_dx->GetBinWidth(1);

	//Overlay plots	
	if( crosstalk && overlay_crosstalk){

				//--------------------------------------------------------------
		//-----------------XTALK OFF ELLIPSE---------------------------------
			
		//ellipse constants...
		double p_ellipse_a_xtalk_OFF = 1.5*dy_pn_sigma_xtalk_OFF;
		double p_ellipse_b_xtalk_OFF = 1.5*dx_p_sigma_xtalk_OFF;
		double n_ellipse_a_xtalk_OFF = 1.25*dy_pn_sigma_xtalk_OFF;
		double n_ellipse_b_xtalk_OFF = 1.25*dx_n_sigma_xtalk_OFF;

		p_ellipse_xtalk_OFF = new TEllipse(dy_pn_xtalk_OFF, dx_p_xtalk_OFF, p_ellipse_a_xtalk_OFF, p_ellipse_b_xtalk_OFF);
		p_ellipse_xtalk_OFF->SetLineColor(2);
		p_ellipse_xtalk_OFF->SetFillStyle(4005);
		p_ellipse_xtalk_OFF->SetLineWidth(3);
		p_ellipse_xtalk_OFF->Draw("same");

		n_ellipse_xtalk_OFF = new TEllipse(dy_pn_xtalk_OFF, dx_n_xtalk_OFF, n_ellipse_a_xtalk_OFF, n_ellipse_b_xtalk_OFF);
		n_ellipse_xtalk_OFF->SetLineColor(6);
		n_ellipse_xtalk_OFF->SetFillStyle(4005);
		n_ellipse_xtalk_OFF->SetLineWidth(3);
		n_ellipse_xtalk_OFF->Draw("same");

		//count steps...
		p_ellipse_cnt_xtalk_OFF = cnt_ellipse_steps(dy_pn_xtalk_OFF, dx_p_xtalk_OFF, p_ellipse_a_xtalk_OFF, p_ellipse_b_xtalk_OFF);
		n_ellipse_cnt_xtalk_OFF = cnt_ellipse_steps(dy_pn_xtalk_OFF, dx_n_xtalk_OFF, n_ellipse_a_xtalk_OFF, n_ellipse_b_xtalk_OFF);

		tcg_p_xtalk_OFF = new TCutG("tcg_p", p_ellipse_cnt_xtalk_OFF);
		tcg_p_xtalk_OFF->SetVarX("x");
		tcg_p_xtalk_OFF->SetVarY("y");
		set_tcg_ellipse(dy_pn_xtalk_OFF, dx_p_xtalk_OFF, p_ellipse_a_xtalk_OFF, p_ellipse_b_xtalk_OFF, tcg_p_xtalk_OFF);
		p_ellipse_cnt_xtalk_OFF = tcg_p_xtalk_OFF->IntegralHist(hin_dxdy_xtalk_OFF)/hin_dx_xtalk_OFF->GetBinWidth(1);

		tcg_n_xtalk_OFF = new TCutG("tcg_n_xtalk_OFF", n_ellipse_cnt_xtalk_OFF);
		tcg_n_xtalk_OFF->SetVarX("x");
		tcg_n_xtalk_OFF->SetVarY("y");
		set_tcg_ellipse(dy_pn_xtalk_OFF, dx_n_xtalk_OFF, n_ellipse_a_xtalk_OFF, n_ellipse_b_xtalk_OFF, tcg_n_xtalk_OFF);
		n_ellipse_cnt_xtalk_OFF = tcg_n_xtalk_OFF->IntegralHist(hin_dxdy_xtalk_OFF)/hin_dx_xtalk_OFF->GetBinWidth(1);

		//--------------------------------------------------------------
	//----------------- XTALK OFF dx-------------------------------
		cout << "Fitting XTALK OFF dx" << endl;
		TCanvas *c_dx_xtalk_OFF = new TCanvas("c_dx_xtalk_OFF", "c_dx_xtalk_OFF", 600, 500);
		hin_dx_xtalk_OFF->Draw();

		if( SBS_field != 0 ){num_par = 9;}
		if( SBS_field == 0 ){num_par = 6;}
		cout << "--------------------------------- " << endl;
		cout << "492" << endl;
		TF1 *fit_dx_xtalk_OFF = new TF1("fit_dx_xtalk_OFF", fit_dxdy, -3, 3, 9);

		if( sbsfieldscale == 0 ){
			cout << "sbsfield scale = " << sbsfieldscale << ". Num params: " << num_par << endl;
			fit_dx_xtalk_OFF->SetParName(0, "pn - Norm (xtalk_OFF)");
			fit_dx_xtalk_OFF->SetParName(1, "pn - Center (xtalk_OFF)");
			fit_dx_xtalk_OFF->SetParName(2, "pn - Sigma (xtalk_OFF)");
			fit_dx_xtalk_OFF->SetParName(3, "BG_pn - Norm (xtalk_OFF)");
			fit_dx_xtalk_OFF->SetParName(4, "BG_pn - Center (xtalk_OFF)");
			fit_dx_xtalk_OFF->SetParName(5, "BG_pn - Sigma (xtalk_OFF)");
			
			fit_dx_xtalk_OFF->SetParLimits(0, 0, hin_dx_xtalk_OFF->GetMaximum());
			fit_dx_xtalk_OFF->SetParLimits(1, -1, 0);
			fit_dx_xtalk_OFF->SetParLimits(2, 0.05, 0.3);
			fit_dx_xtalk_OFF->SetParLimits(3, 0, (0.25)*hin_dx_xtalk_OFF->GetMaximum());
			fit_dx_xtalk_OFF->SetParLimits(4, -0.4, 0.2);
			fit_dx_xtalk_OFF->SetParLimits(5, 1, 5);
		}
		else{
			cout << "sbsfield scale = " << sbsfieldscale << ". Num params: " << num_par << endl;
			fit_dx_xtalk_OFF->SetParName(0, "p - Norm (xtalk_OFF)");
			fit_dx_xtalk_OFF->SetParName(1, "p - Center (xtalk_OFF)");
			fit_dx_xtalk_OFF->SetParName(2, "p - Sigma (xtalk_OFF)");
			fit_dx_xtalk_OFF->SetParName(3, "n - Norm (xtalk_OFF)");
			fit_dx_xtalk_OFF->SetParName(4, "n - Center (xtalk_OFF)");
			fit_dx_xtalk_OFF->SetParName(5, "n - Sigma (xtalk_OFF)");
			fit_dx_xtalk_OFF->SetParName(6, "BG_pn - Norm (xtalk_OFF)");
			fit_dx_xtalk_OFF->SetParName(7, "BG_pn - Center (xtalk_OFF)");
			fit_dx_xtalk_OFF->SetParName(8, "BG_pn - Sigma (xtalk_OFF)");
			
			fit_dx_xtalk_OFF->SetParLimits(0, 0, hin_dx_xtalk_OFF->GetMaximum());
			fit_dx_xtalk_OFF->SetParLimits(1, -0.8, -0.2);
			fit_dx_xtalk_OFF->SetParLimits(2, 0.05, 0.3);
			fit_dx_xtalk_OFF->SetParLimits(3, (0.10)*hin_dx_xtalk_OFF->GetMaximum(), (0.5)*hin_dx->GetMaximum());
			fit_dx_xtalk_OFF->SetParLimits(4, 0.0, 0.6);
			fit_dx_xtalk_OFF->SetParLimits(5, 0.1, 0.3);
			fit_dx_xtalk_OFF->SetParLimits(6, 1, (0.05)*hin_dx_xtalk_OFF->GetMaximum());
			fit_dx_xtalk_OFF->SetParLimits(7, -0.4, 0.2);
			fit_dx_xtalk_OFF->SetParLimits(8, 0.5, 5);
		}
		for(int i=0; i < 10; i++){
			hin_dx_xtalk_OFF->Fit("fit_dx_xtalk_OFF", "R+");
		
			fit_dx_xtalk_OFF->GetParameters(&par[0]);
			fit_dx_xtalk_OFF->SetParameters(par);
		}

		if( sbsfieldscale == 0 ){
			dx_p_xtalk_OFF = fit_dx_xtalk_OFF->GetParameter(1);
			dx_p_sigma_xtalk_OFF = fit_dx_xtalk_OFF->GetParameter(2);
			dx_n_xtalk_OFF = dx_p_xtalk_OFF;
			dx_n_sigma_xtalk_OFF = dx_p_sigma_xtalk_OFF;
			dx_pn_max_xtalk_OFF = 0.05;
			BG_dx_xtalk_OFF = fit_dx_xtalk_OFF->GetParameter(3);
		}
		else{
			dx_p_xtalk_OFF = fit_dx_xtalk_OFF->GetParameter(1);
			cout << "-------------------------- dx_p: " << fit_dx->GetParameter(1) << endl;
			dx_p_sigma_xtalk_OFF = fit_dx_xtalk_OFF->GetParameter(2);
			dx_n_xtalk_OFF = fit_dx_xtalk_OFF->GetParameter(4);
			dx_n_sigma_xtalk_OFF = fit_dx_xtalk_OFF->GetParameter(5);
			dx_pn_max_xtalk_OFF = (dx_n_xtalk_OFF + dx_n_sigma_xtalk_OFF) - (dx_p_xtalk_OFF - dx_p_sigma_xtalk_OFF);
			BG_dx_xtalk_OFF = fit_dx_xtalk_OFF->GetParameter(6);
		//-------------------------------------------------------------------------
			p_integral_xtalk_OFF = (fit_dx_xtalk_OFF->GetParameter(0))*(fit_dx_xtalk_OFF->GetParameter(2))*sqrt(2*pi)/hin_dx_xtalk_OFF->GetBinWidth(1);
			n_integral_xtalk_OFF = (fit_dx_xtalk_OFF->GetParameter(3))*(fit_dx_xtalk_OFF->GetParameter(5))*sqrt(2*pi)/hin_dx_xtalk_OFF->GetBinWidth(1);
		}

		hin_dx_xtalk_OFF->Fit("fit_dx_xtalk_OFF", "R+");
		fit_dx_xtalk_OFF->Draw("same");
		plot_dxdy_fit(par);

	//--------------------------------------------------------------
	//-----------------dy _xtalk_OFF-------------------------------
		cout << "Fitting _xtalk_OFF dy" << endl;
		TCanvas *c_dy_xtalk_OFF = new TCanvas("c_dy_xtalk_OFF", "c_dy_xtalk_OFF", 600, 500);
		hin_dy_xtalk_OFF->Draw();

		num_par = 6;
cout << "--------------------------------- " << endl;
		cout << "571" << endl;
		TF1 *tf_fit_dy_xtalk_OFF = new TF1("tf_fit_dy_xtalk_OFF", fit_dy, -3, 3, 6);
		tf_fit_dy_xtalk_OFF->SetParName(0, "dy - Norm (xtalk_OFF)");
		tf_fit_dy_xtalk_OFF->SetParName(1, "dy - Center (xtalk_OFF)");
		tf_fit_dy_xtalk_OFF->SetParName(2, "dy - Sigma (xtalk_OFF)");
		tf_fit_dy_xtalk_OFF->SetParName(3, "BG_dy - Norm (xtalk_OFF)");
		tf_fit_dy_xtalk_OFF->SetParName(4, "BG_dy - Center (xtalk_OFF)");
		tf_fit_dy_xtalk_OFF->SetParName(5, "BG_dy - Sigma (xtalk_OFF)");

		tf_fit_dy_xtalk_OFF->SetParLimits(0, 0, hin_dy_xtalk_OFF->GetMaximum());
		tf_fit_dy_xtalk_OFF->SetParLimits(1, -0.4, 0.2);
		tf_fit_dy_xtalk_OFF->SetParLimits(2, 0.05, 0.3);
		tf_fit_dy_xtalk_OFF->SetParLimits(3, 1, BG_dx_xtalk_OFF);
		tf_fit_dy_xtalk_OFF->SetParLimits(4, -0.2, 0.2);
		tf_fit_dy_xtalk_OFF->SetParLimits(5, 0.5, 5);

		hin_dy_xtalk_OFF->Fit("tf_fit_dy_xtalk_OFF", "R+");
		tf_fit_dy_xtalk_OFF->Draw("same");

		tf_fit_dy_xtalk_OFF->GetParameters(&par[0]);
		tf_fit_dy_xtalk_OFF->SetParameters(par);

		dy_pn_xtalk_OFF = tf_fit_dy_xtalk_OFF->GetParameter(1);
		dy_pn_sigma_xtalk_OFF = tf_fit_dy_xtalk_OFF->GetParameter(2);
		BG_dy_xtalk_OFF = tf_fit_dy_xtalk_OFF->GetParameter(3);
		plot_dy_fit(par);

	//--------------------------------------------------------------
	//-----------------dxdy---------------------------------
		cout << "Fitting _xtalk_OFF dxdy" << endl;
		TCanvas *c_dxdy_xtalk_OFF = new TCanvas("c_dxdy_xtalk_OFF", "c_dxdy_xtalk_OFF", 600, 500);
		hin_dxdy_xtalk_OFF->Draw("colz");

	
	//--------------------------------------------------------------
	//-----------------XTALK ON ELLIPSE---------------------------------
		
		//ellipse constants...
		double p_ellipse_a_xtalk_ON = 1.5*dy_pn_sigma_xtalk_ON;
		double p_ellipse_b_xtalk_ON = 1.5*dx_p_sigma_xtalk_ON;
		double n_ellipse_a_xtalk_ON = 1.25*dy_pn_sigma_xtalk_ON;
		double n_ellipse_b_xtalk_ON = 1.25*dx_n_sigma_xtalk_ON;

		p_ellipse_xtalk_ON = new TEllipse(dy_pn_xtalk_ON, dx_p_xtalk_ON, p_ellipse_a_xtalk_ON, p_ellipse_b_xtalk_ON);
		p_ellipse_xtalk_ON->SetLineColor(2);
		p_ellipse_xtalk_ON->SetFillStyle(4005);
		p_ellipse_xtalk_ON->SetLineWidth(3);
		p_ellipse_xtalk_ON->Draw("same");

		n_ellipse_xtalk_ON = new TEllipse(dy_pn_xtalk_ON, dx_n_xtalk_ON, n_ellipse_a_xtalk_ON, n_ellipse_b_xtalk_ON);
		n_ellipse_xtalk_ON->SetLineColor(6);
		n_ellipse_xtalk_ON->SetFillStyle(4005);
		n_ellipse_xtalk_ON->SetLineWidth(3);
		n_ellipse_xtalk_ON->Draw("same");

		//count steps...
		p_ellipse_cnt_xtalk_ON = cnt_ellipse_steps(dy_pn_xtalk_ON, dx_p_xtalk_ON, p_ellipse_a_xtalk_ON, p_ellipse_b_xtalk_ON);
		n_ellipse_cnt_xtalk_ON = cnt_ellipse_steps(dy_pn_xtalk_ON, dx_n_xtalk_ON, n_ellipse_a_xtalk_ON, n_ellipse_b_xtalk_ON);

		tcg_p_xtalk_ON = new TCutG("tcg_p_xtalk_ON", p_ellipse_cnt_xtalk_ON);
		tcg_p_xtalk_ON->SetVarX("x");
		tcg_p_xtalk_ON->SetVarY("y");
		set_tcg_ellipse(dy_pn_xtalk_ON, dx_p_xtalk_ON, p_ellipse_a_xtalk_ON, p_ellipse_b_xtalk_ON, tcg_p_xtalk_ON);
		p_ellipse_cnt_xtalk_ON = tcg_p_xtalk_ON->IntegralHist(hin_dxdy_xtalk_ON)/hin_dx_xtalk_ON->GetBinWidth(1);

		tcg_n_xtalk_ON = new TCutG("tcg_n_xtalk_ON", n_ellipse_cnt_xtalk_ON);
		tcg_n_xtalk_ON->SetVarX("x");
		tcg_n_xtalk_ON->SetVarY("y");
		set_tcg_ellipse(dy_pn_xtalk_ON, dx_n_xtalk_ON, n_ellipse_a_xtalk_ON, n_ellipse_b_xtalk_ON, tcg_n_xtalk_ON);
		n_ellipse_cnt_xtalk_ON = tcg_n_xtalk_ON->IntegralHist(hin_dxdy_xtalk_ON)/hin_dx_xtalk_ON->GetBinWidth(1);
	
		//--------------------------------------------------------------
	//----------------- XTALK ON dx-------------------------------
		cout << "Fitting XTALK ON dx" << endl;
		TCanvas *c_dx_xtalk_ON = new TCanvas("c_dx_xtalk_ON", "c_dx_xtalk_ON", 600, 500);
		hin_dx_xtalk_ON->Draw();

		if( SBS_field != 0 ){num_par = 9;}
		if( SBS_field == 0 ){num_par = 6;}
		cout << "--------------------------------- " << endl;
		cout << "649" << endl;
		TF1 *fit_dx_xtalk_ON = new TF1("fit_dx_xtalk_ON", fit_dxdy, -3, 3, 9);

		if( sbsfieldscale == 0 ){
			cout << "sbsfield scale = " << sbsfieldscale << ". Num params: " << num_par << endl;
			fit_dx_xtalk_ON->SetParName(0, "pn - Norm (xtalk_ON)");
			fit_dx_xtalk_ON->SetParName(1, "pn - Center (xtalk_ON)");
			fit_dx_xtalk_ON->SetParName(2, "pn - Sigma (xtalk_ON)");
			fit_dx_xtalk_ON->SetParName(3, "BG_pn - Norm (xtalk_ON)");
			fit_dx_xtalk_ON->SetParName(4, "BG_pn - Center (xtalk_ON)");
			fit_dx_xtalk_ON->SetParName(5, "BG_pn - Sigma (xtalk_ON)");
			
			fit_dx_xtalk_ON->SetParLimits(0, 0, hin_dx_xtalk_ON->GetMaximum());
			fit_dx_xtalk_ON->SetParLimits(1, -1, 0);
			fit_dx_xtalk_ON->SetParLimits(2, 0.05, 0.3);
			fit_dx_xtalk_ON->SetParLimits(3, 0, (0.25)*hin_dx_xtalk_ON->GetMaximum());
			fit_dx_xtalk_ON->SetParLimits(4, -0.4, 0.2);
			fit_dx_xtalk_ON->SetParLimits(5, 1, 5);
		}
		else{
			cout << "sbsfield scale = " << sbsfieldscale << ". Num params: " << num_par << endl;
			fit_dx_xtalk_ON->SetParName(0, "p - Norm (xtalk_ON)");
			fit_dx_xtalk_ON->SetParName(1, "p - Center (xtalk_ON)");
			fit_dx_xtalk_ON->SetParName(2, "p - Sigma (xtalk_ON)");
			fit_dx_xtalk_ON->SetParName(3, "n - Norm (xtalk_ON)");
			fit_dx_xtalk_ON->SetParName(4, "n - Center (xtalk_ON)");
			fit_dx_xtalk_ON->SetParName(5, "n - Sigma (xtalk_ON)");
			fit_dx_xtalk_ON->SetParName(6, "BG_pn - Norm (xtalk_ON)");
			fit_dx_xtalk_ON->SetParName(7, "BG_pn - Center (xtalk_ON)");
			fit_dx_xtalk_ON->SetParName(8, "BG_pn - Sigma (xtalk_ON)");
			
			fit_dx_xtalk_ON->SetParLimits(0, 0, hin_dx_xtalk_ON->GetMaximum());
			fit_dx_xtalk_ON->SetParLimits(1, -0.8, -0.2);
			fit_dx_xtalk_ON->SetParLimits(2, 0.05, 0.3);
			fit_dx_xtalk_ON->SetParLimits(3, (0.10)*hin_dx_xtalk_ON->GetMaximum(), (0.5)*hin_dx->GetMaximum());
			fit_dx_xtalk_ON->SetParLimits(4, 0.0, 0.6);
			fit_dx_xtalk_ON->SetParLimits(5, 0.1, 0.3);
			fit_dx_xtalk_ON->SetParLimits(6, 1, (0.05)*hin_dx_xtalk_ON->GetMaximum());
			fit_dx_xtalk_ON->SetParLimits(7, -0.4, 0.2);
			fit_dx_xtalk_ON->SetParLimits(8, 0.5, 5);
		}
		for(int i=0; i < 10; i++){
			hin_dx_xtalk_ON->Fit("fit_dx_xtalk_ON", "R+");
		
			fit_dx_xtalk_ON->GetParameters(&par[0]);
			fit_dx_xtalk_ON->SetParameters(par);
		}

		if( sbsfieldscale == 0 ){
			dx_p_xtalk_ON = fit_dx_xtalk_ON->GetParameter(1);
			dx_p_sigma_xtalk_ON = fit_dx_xtalk_ON->GetParameter(2);
			dx_n_xtalk_ON = dx_p_xtalk_ON;
			dx_n_sigma_xtalk_ON = dx_p_sigma_xtalk_ON;
			dx_pn_max_xtalk_ON = 0.05;
			BG_dx_xtalk_ON = fit_dx_xtalk_ON->GetParameter(3);
		}
		else{
			dx_p_xtalk_ON = fit_dx_xtalk_ON->GetParameter(1);
			cout << "-------------------------- dx_p: " << fit_dx->GetParameter(1) << endl;
			dx_p_sigma_xtalk_ON = fit_dx_xtalk_ON->GetParameter(2);
			dx_n_xtalk_ON = fit_dx_xtalk_ON->GetParameter(4);
			dx_n_sigma_xtalk_ON = fit_dx_xtalk_ON->GetParameter(5);
			dx_pn_max_xtalk_ON = (dx_n_xtalk_ON + dx_n_sigma_xtalk_ON) - (dx_p_xtalk_ON - dx_p_sigma_xtalk_ON);
			BG_dx_xtalk_ON = fit_dx_xtalk_ON->GetParameter(6);
		//-------------------------------------------------------------------------
			p_integral_xtalk_ON = (fit_dx_xtalk_ON->GetParameter(0))*(fit_dx_xtalk_ON->GetParameter(2))*sqrt(2*pi)/hin_dx_xtalk_ON->GetBinWidth(1);
			n_integral_xtalk_ON = (fit_dx_xtalk_ON->GetParameter(3))*(fit_dx_xtalk_ON->GetParameter(5))*sqrt(2*pi)/hin_dx_xtalk_ON->GetBinWidth(1);
		}

		hin_dx_xtalk_ON->Fit("fit_dx_xtalk_ON", "R+");
		fit_dx_xtalk_ON->Draw("same");
		plot_dxdy_fit(par);

	//--------------------------------------------------------------
	//-----------------dy _xtalk_ON-------------------------------
		cout << "Fitting _xtalk_ON dy" << endl;
		TCanvas *c_dy_xtalk_ON = new TCanvas("c_dy_xtalk_ON", "c_dy_xtalk_ON", 600, 500);
		hin_dy_xtalk_ON->Draw();

		num_par = 6;

		TF1 *tf_fit_dy_xtalk_ON = new TF1("tf_fit_dy_xtalk_ON", fit_dy, -3, 3, 6);
		tf_fit_dy_xtalk_ON->SetParName(0, "dy - Norm (xtalk_ON)");
		tf_fit_dy_xtalk_ON->SetParName(1, "dy - Center (xtalk_ON)");
		tf_fit_dy_xtalk_ON->SetParName(2, "dy - Sigma (xtalk_ON)");
		tf_fit_dy_xtalk_ON->SetParName(3, "BG_dy - Norm (xtalk_ON)");
		tf_fit_dy_xtalk_ON->SetParName(4, "BG_dy - Center (xtalk_ON)");
		tf_fit_dy_xtalk_ON->SetParName(5, "BG_dy - Sigma (xtalk_ON)");

		tf_fit_dy_xtalk_ON->SetParLimits(0, 0, hin_dy_xtalk_ON->GetMaximum());
		tf_fit_dy_xtalk_ON->SetParLimits(1, -0.4, 0.2);
		tf_fit_dy_xtalk_ON->SetParLimits(2, 0.05, 0.3);
		tf_fit_dy_xtalk_ON->SetParLimits(3, 1, BG_dx_xtalk_ON);
		tf_fit_dy_xtalk_ON->SetParLimits(4, -0.2, 0.2);
		tf_fit_dy_xtalk_ON->SetParLimits(5, 0.5, 5);

		hin_dy_xtalk_ON->Fit("tf_fit_dy_xtalk_ON", "R+");
		tf_fit_dy_xtalk_ON->Draw("same");

		tf_fit_dy_xtalk_ON->GetParameters(&par[0]);
		tf_fit_dy_xtalk_ON->SetParameters(par);

		dy_pn_xtalk_ON = tf_fit_dy_xtalk_ON->GetParameter(1);
		dy_pn_sigma_xtalk_ON = tf_fit_dy_xtalk_ON->GetParameter(2);
		BG_dy_xtalk_ON = tf_fit_dy_xtalk_ON->GetParameter(3);
		plot_dy_fit(par);

	//--------------------------------------------------------------
	//-----------------dxdy---------------------------------
		cout << "Fitting _xtalk_ON dxdy" << endl;
		TCanvas *c_dxdy_xtalk_ON = new TCanvas("c_dxdy_xtalk_ON", "c_dxdy_xtalk_ON", 600, 500);
		hin_dxdy_xtalk_ON->Draw("colz");

	
	//PLOT XTALK STUFF

		TCanvas *c_xtalk_OFF_dx = new TCanvas("c_xtalk_OFF_dx", "c_xtalk_OFF_dx", 600, 500);
		hin_dx_xtalk_OFF->Draw();

		TCanvas *c_xtalk_ON_dy = new TCanvas("c_xtalk_ON_dy", "c_xtalk_ON_dy", 600, 500);
		hin_dx_xtalk_ON->Draw();

		TCanvas *c_xtalk_overlay = new TCanvas("c_xtalk_overlay", "c_xtalk_overlay", 600, 500);
		hin_dx_xtalk_OFF->Draw();
		hin_dx_xtalk_ON->SetLineColor(6);
		hin_dx_xtalk_ON->Draw("same");

		TLegend *tl_overlay = new TLegend(0.15, 0.75, 0.35, 0.85, "", "NDC");
		tl_overlay->AddEntry(hin_dx_xtalk_OFF, "Xtalk OFF");
		tl_overlay->AddEntry(hin_dx_xtalk_ON, "Xtalk ON");
		tl_overlay->Draw("same");

//--------------- DIFF
		// h_xtalk_dx_diff = new TH1D("h_xtalk_dx_diff", "h_xtalK_dx_diff", 250, -2.5, 2.5);
		hin_dx_xtalk_OFF->Copy(h_xtalk_dx_diff);
		h_xtalk_dx_diff.SetNameTitle("h_xtalk_dx_diff", "Xtalk_OFF_dx - Xtalk_ON_dx (Difference)");
		h_xtalk_dx_diff.Add(hin_dx_xtalk_ON, -1);
		// for(int bin = 0; bin < hin_xtalk_OFF_dx->GetNbinsX(); bin++){
		// 	Double_t OFF_entry = hin_xtalk_OFF_dx->GetBinContent(bin);
		// 	Double_t ON_entry = hin_xtalk_ON_dx->GetBinContent(bin);

		// 	Double_t diff = OFF_entry - ON_entry;
		// 	h_xtalk_dx_diff.SetBinContent(bin, diff);
		// }


		TCanvas *c_xtalk_dx_diff = new TCanvas("c_xtalk_dx_diff", "c_xtalk_dx_diff", 600, 500);
		h_xtalk_dx_diff.Draw();

//--------------------------------------------------------------
//--------------------------------------------------------------

	}

//END OF XTALK STUFF



//Open up Crosstalk file and plot the histograms on top of each other
	if( false ){
		crosstalk_infile = new TFile(outfile->GetName(), "READ");
		hin_xtalk_OFF_dx = static_cast<TH1D*>(crosstalk_infile->Get("h_xtalk_OFF_dx"));
		hin_xtalk_ON_dx = static_cast<TH1D*>(crosstalk_infile->Get("h_xtalk_ON_dx"));

		TCanvas *c_xtalk_OFF_dx = new TCanvas("c_xtalk_OFF_dx", "c_xtalk_OFF_dx", 600, 500);
		hin_xtalk_OFF_dx->Draw();

		TCanvas *c_xtalk_ON_dy = new TCanvas("c_xtalk_ON_dy", "c_xtalk_ON_dy", 600, 500);
		hin_xtalk_ON_dx->Draw();

		TCanvas *c_xtalk_overlay = new TCanvas("c_xtalk_overlay", "c_xtalk_overlay", 600, 500);
		hin_xtalk_OFF_dx->Draw();
		hin_xtalk_ON_dx->SetLineColor(6);
		hin_xtalk_ON_dx->Draw("same");

		TLegend *tl_overlay = new TLegend(0.15, 0.75, 0.35, 0.85, "", "NDC");
		tl_overlay->AddEntry(hin_xtalk_OFF_dx, "Xtalk OFF");
		tl_overlay->AddEntry(hin_xtalk_ON_dx, "Xtalk ON");
		tl_overlay->Draw("same");

//--------------- DIFF
		// h_xtalk_dx_diff = new TH1D("h_xtalk_dx_diff", "h_xtalK_dx_diff", 250, -2.5, 2.5);
		hin_xtalk_OFF_dx->Copy(h_xtalk_dx_diff);
		h_xtalk_dx_diff.SetNameTitle("h_xtalk_dx_diff", "Xtalk_OFF_dx - Xtalk_ON_dx (Difference)");
		h_xtalk_dx_diff.Add(hin_xtalk_ON_dx, -1);
		// for(int bin = 0; bin < hin_xtalk_OFF_dx->GetNbinsX(); bin++){
		// 	Double_t OFF_entry = hin_xtalk_OFF_dx->GetBinContent(bin);
		// 	Double_t ON_entry = hin_xtalk_ON_dx->GetBinContent(bin);

		// 	Double_t diff = OFF_entry - ON_entry;
		// 	h_xtalk_dx_diff.SetBinContent(bin, diff);
		// }


		TCanvas *c_xtalk_dx_diff = new TCanvas("c_xtalk_dx_diff", "c_xtalk_dx_diff", 600, 500);
		h_xtalk_dx_diff.Draw();

	}
	

	if( sbsfieldscale == 0 ){
		dxdy_vec.push_back(runnum);
		dxdy_vec.push_back(dx_p);
		dxdy_vec.push_back(dx_p_sigma);
		dxdy_vec.push_back(-1);
		dxdy_vec.push_back(-1);
		dxdy_vec.push_back(dx_n);
		dxdy_vec.push_back(dx_n_sigma);
		dxdy_vec.push_back(-1);
		dxdy_vec.push_back(-1);
		dxdy_vec.push_back(dx_pn_max);	
	}

	// TCanvas *c_dx = new TCanvas("c_dx", "c_dx", 600, 500);

	cout << "-----------------------------" << endl;
	cout << "dxdy vector: " << endl;
	cout << "{ runnum, dx_p, dx_p_sigma, dx_n, dx_n_sigma, dy, dy_sigma, dx_pn_max }" << endl << endl;
	cout << "{";
	cout << runnum << ", " << dx_p << ", " << dx_p_sigma << ", " << dx_n << ", " << dx_n_sigma;
	cout << ", " << dy_pn <<", "<< dy_pn_sigma << ", " << dx_pn_max << "}" << endl;
	cout << "-----------------------------" << endl;
	cout << "p integral: " << int(p_integral) << endl;
	cout << "n_integral: " << int(n_integral) << endl;
	cout << "-----------------------------" << endl;
	cout << "p ellipse cnt: " << p_ellipse_cnt << endl;
	cout << "n ellipse cnt: " << n_ellipse_cnt << endl;

	cout << endl << "-----------------------------" << endl;
	cout << "Run Parameters: " << endl;
	cout << "Run: " << runnum << endl;
	cout << "SBS Field: " << sbsfieldscale << "%" << endl;
	cout << "----------------------------" << endl;
	cout << "Crosstalk ON?    " << crosstalk << endl;
	cout << "XTALK On/Off? ---- " << XTALK_ONOFF.Data() << endl;
	cout << "Ratio threshold: " << ratio_threshold << endl;
	cout << "----------------------------" << endl;
	cout << "Selected file: " << select_file << endl;
	cout << "Input file: " << infile->GetName() << endl;
	cout << "----------------------------" << endl;
	if( crosstalk ){
		if( overlay_crosstalk ){
			cout << "Crosstalk output file created: " << outfile->GetName() << endl;
		}
	}
	cout << "----------------------------" << endl << endl;
	auto total_time_end = high_resolution_clock::now();
	auto total_time_duration = duration_cast<minutes>(total_time_end - total_time_start);
	cout << "Total time for analysis: " << total_time_duration.count() << " minutes. " << endl;

}