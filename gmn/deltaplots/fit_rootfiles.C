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
int runnum = 13566;
double SBS_field = lookup_run_info(runnum, "sbs_field"); //Strength (in percentage) of SBS magnet
int sbsfieldscale = int(100*lookup_run_info(runnum, "sbs_field"));
double pi = TMath::Pi();
Int_t num_par;

double dy_pn, dy_pn_sigma, dx_p, dx_p_sigma, dx_n, dx_n_sigma, dx_pn_max, BG_dx, BG_dy;
int p_ellipse_cnt, n_ellipse_cnt;
double p_integral, n_integral, pn_integral, p_count_ellipse, n_count_ellipse, pn_count_ellipse;
vector<double> dxdy_vec;
TCutG *tcg_p, *tcg_n;

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

TString select_file = "calib_W";

TString rootfile_dir = Form("/w/halla-scshelf2102/sbs/jboyd/analysis/gmn/deltaplots/rootfiles");

TFile *infile;
TString infile_dxdy = Form("%i_dxdy.root", runnum);
TString infile_calib = Form("%i_calibrate.root", runnum);
TString infile_calib_W = Form("%i_calibrate_W.root", runnum);
TEllipse *p_ellipse, *n_ellipse, *pn_ellipse;

TH1D *hin_dx, *hin_dx_cut, *hin_dx_wcut, *hin_dy, *hin_dy_cut, *hin_dy_wcut;
TH2D *hin_dxdy, *hin_dxdy_cut, *hin_dxdy_wcut;

Double_t par[9];

void fit_rootfiles(){

	if( select_file=="dxdy" ){
		infile = new TFile(Form("%s/%s", rootfile_dir.Data(), infile_dxdy.Data()), "READ");
	}
	if( select_file=="calib" ){
		infile = new TFile(Form("%s/%s", rootfile_dir.Data(), infile_calib.Data()), "READ");
	}
	if( select_file=="calib_W" ){
		infile = new TFile(Form("%s/%s", rootfile_dir.Data(), infile_calib_W.Data()), "READ");
	}

	hin_dxdy = static_cast<TH2D*>(infile->Get("h_dxdy"));
	hin_dxdy_cut = static_cast<TH2D*>(infile->Get("h_dxdy_cut"));
	hin_dxdy_wcut = static_cast<TH2D*>(infile->Get("h_dxdy_wcut"));

	hin_dx = static_cast<TH1D*>(infile->Get("h_dx"));
	hin_dx_cut = static_cast<TH1D*>(infile->Get("h_dx_cut"));
	hin_dx_wcut = static_cast<TH1D*>(infile->Get("h_dx_wcut"));

	hin_dy = static_cast<TH1D*>(infile->Get("h_dy"));
	hin_dy_cut = static_cast<TH1D*>(infile->Get("h_dy_cut"));
	hin_dy_wcut = static_cast<TH1D*>(infile->Get("h_dy_wcut"));

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
		fit_dx->SetParLimits(4, 0.2, 0.6);
		fit_dx->SetParLimits(5, 0.1, 0.3);
		fit_dx->SetParLimits(6, 50, (0.25)*hin_dx->GetMaximum());
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
		p_integral = (fit_dx->GetParameter(0))*(fit_dx->GetParameter(2))*sqrt(2*pi);
		n_integral = (fit_dx->GetParameter(3))*(fit_dx->GetParameter(5))*sqrt(2*pi);
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
	tf_fit_dy->SetParLimits(3, (0.99)*BG_dx, (1.01)*BG_dx);
	tf_fit_dy->SetParLimits(4, -0.4, 0.4);
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
	TCanvas *c_dxdy = new TCanvas("c_dxdy", "c_dxdy", 600, 500);
	hin_dxdy_wcut->Draw("colz");

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
	p_ellipse_cnt = tcg_p->IntegralHist(hin_dxdy_wcut);

	tcg_n = new TCutG("tcg_n", n_ellipse_cnt);
	tcg_n->SetVarX("x");
	tcg_n->SetVarY("y");
	set_tcg_ellipse(dy_pn, dx_n, n_ellipse_a, n_ellipse_b, tcg_n);
	n_ellipse_cnt = tcg_n->IntegralHist(hin_dxdy_wcut);

//--------------------------------------------------------------
//--------------------------------------------------------------

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

	cout << "Fitting dxdy" << endl;
	// TCanvas *c_dx = new TCanvas("c_dx", "c_dx", 600, 500);

	cout << "-----------------------------" << endl;
	cout << "dxdy vector: " << endl;
	cout << "{ runnum, dx_p, dx_p_sigma, dx_n, dx_n_sigma, dy, dy_sigma, dx_pn_max }" << endl << endl;
	cout << "{";
	cout << runnum << ", " << dx_p << ", " << dx_p_sigma << dx_n << ", " << dx_n_sigma;
	cout << ", " << dy_pn <<", "<< dy_pn_sigma << ", " << dx_pn_max << "}" << endl;
	cout << "-----------------------------" << endl;
	cout << "p integral: " << int(p_integral) << endl;
	cout << "n_integral: " << int(n_integral) << endl;
	cout << "-----------------------------" << endl;
	cout << "p ellipse cnt: " << p_ellipse_cnt << endl;
	cout << "n ellipse cnt: " << n_ellipse_cnt << endl;

}