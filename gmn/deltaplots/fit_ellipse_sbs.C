#include <iostream>
#include <fstream>
#include <numeric>
#include <iomanip>
#include <math.h>
#include <algorithm>
#include <string>
#include <chrono>
#include "TCanvas.h"
#include "TFile.h"
#include "TEllipse.h"
#include "TH2F.h"
#include "TMath.h"
#include "TMinuit.h"
#include "TRandom.h"
#include "TStyle.h"
#include "iostream"

using namespace std::chrono;
#include "/w/halla-scshelf2102/sbs/jboyd/include/include_files.h"
#include "/w/halla-scshelf2102/sbs/jboyd/include/GEM_lookups.h"
#include "/w/halla-scshelf2102/sbs/jboyd/include/beam_variables.h"

//Run info and lookups
int runnum = 13479;
vector<double> integrals;

Double_t params[5];

TH2D *h2;

TFile *TF = new TFile(Form("%i_dxdy.root", runnum), "READ");
TEllipse * el;
TEllipse *el_cut;
TEllipse *el_p;
double pi = TMath::Pi();
Double_t dx0_p, dx_sig_p, dy0_p, dy_sig_p, first_area, last_area;
double total_integral = 0;
double first_cut_integral = 0;
double final_cut_integral = 0;

//First ellipse constants:
	
double x_0, y_0, r1_0, r2_0, phi_min, phi_max, theta, h_min, h_max, k_min, k_max, r1_0_min, r1_0_max, r2_0_min, r2_0_max;


TCutG *mycutg, *first_el_cut;

Float_t thresh;


double ellipse_integral(Double_t h, Double_t k, Double_t a, Double_t b, TH2D *histo)
{
	// cout << "Setting cut-ellipse with: par[0] = " << par[0] << ", par[1] = " << par[1] << ", par[2] = " << par[2] << ", par[3] = " << par[3] << endl;
	// el_cut = new TEllipse(par[0], par[1], par[2], par[3], 0, 360, 0);
	// el_cut->SetFillStyle(4000);

	// double x_calc = h + (a)*sqrt( 1 - ( pow( (y - k), 2)/(b*b) ) );
	// double y_calc = k + (b)*sqrt( 1 - ( pow( (x - h), 2)/(a*a) ) );
	int mycutg_cnt = 0;
	for(double x_i = (h - a); x_i < (h + a); x_i += 0.001){
	  mycutg_cnt++;
	}

	for(double x_i = (h - a); x_i < (h + a); x_i += 0.001){
	  mycutg_cnt++;
	}

	mycutg = new TCutG("mycutg", mycutg_cnt);
	mycutg->SetVarX("x");
	mycutg->SetVarY("y");

	mycutg_cnt = 0;

	for(double x_i = (h - a); x_i <= (h + a); x_i += 0.001){

		double y_calc = k + (b)*sqrt( 1 - ( pow( (x_i - h), 2)/(a*a) ) );
		mycutg->SetPoint(mycutg_cnt, x_i, y_calc);
		mycutg_cnt++;
	}

	for(double x_i = (h - a); x_i <= (h + a); x_i += 0.001){

		double y_calc = k - (b)*sqrt( 1 - ( pow( (x_i - h), 2)/(a*a) ) );
		mycutg->SetPoint(mycutg_cnt, x_i, y_calc);
		mycutg_cnt++;
	}

	Double_t cut_integral = mycutg->IntegralHist(histo);

	final_cut_integral = cut_integral;

	integrals.push_back(cut_integral);

	return cut_integral;

}


void ellipse_fcn(Int_t &/*npar*/, Double_t */*gin*/, Double_t &f, Double_t *par, Int_t /*iflag*/)
{

	double integ = ellipse_integral(par[0], par[1], par[2], par[3], h2); 
	double area = TMath::Pi()*(par[2])*(par[3]);

	f = -integ;

}

void fit_ellipse_sbs(){

//First ellipse constants:

	h2 = static_cast<TH2D*>(TF->Get("h_dxdy"));

	TH1D *h2_x = new TH1D("h2_x", "h2_x", 400, -4.0, 4.0);
	TH1D *h2_y = new TH1D("h2_y", "h2_y", 400, -4.0, 4.0);

	// TCanvas *c_x = new TCanvas("c_x", "c_x", 600, 500);
	h2_x = h2->ProjectionX();
	// h2_x->Draw();

	TF1 *fit_Px = new TF1("fit_Px", "gaus", -1.5, 1.5);
	fit_Px->SetParName(0, "Proj x norm");
	fit_Px->SetParName(1, "Proj x center");
	fit_Px->SetParName(2, "Proj x sigma");

	fit_Px->SetParLimits(0, .95*h2_x->GetMaximum(), h2_x->GetMaximum());
	fit_Px->SetParLimits(1, 0.05, 0.2);
	fit_Px->SetParLimits(2, 0.05, 0.15);

	h2_x->Fit("fit_Px", "R+");
	// fit_Px->Draw("same");

	// TCanvas *c_y = new TCanvas("c_y", "c_y", 600, 500);
	h2_y = h2->ProjectionY();
	// h2_y->Draw();

	TF1 *fit_Py = new TF1("fit_Py", "gaus", -1.5, 1.5);
	fit_Py->SetParName(0, "Proj y norm");
	fit_Py->SetParName(1, "Proj y center");
	fit_Py->SetParName(2, "Proj y sigma");

	fit_Py->SetParLimits(0, .95*h2_y->GetMaximum(), h2_y->GetMaximum());
	fit_Py->SetParLimits(1, 0.04, 0.5);
	fit_Py->SetParLimits(2, 0.05, 0.15);

	h2_y->Fit("fit_Py", "R+");
	// fit_Py->Draw("same");

	dx0_p = fit_Px->GetParameter(1);
	dy0_p = fit_Py->GetParameter(1);

	dx_sig_p = fit_Px->GetParameter(2);
	dy_sig_p = fit_Py->GetParameter(2);

	x_0 = dx0_p;
	y_0 = dy0_p;
	r1_0 = 1.5*dx_sig_p;
	r2_0 = 1.5*dy_sig_p;
	phi_min = 0;
	phi_max = 360;
	theta = 0;

	cout << "------------------------------------------------------" << endl;
	cout << "dx_sig_p = " << dx_sig_p << "; dy_sig_p = " << dy_sig_p << endl;
	cout << "------------------------------------------------------" << endl;

	TCanvas *c_h2 = new TCanvas("c_h2", "c_h2", 600, 500);
	h2->Draw("colz");
	el_p = new TEllipse(x_0, y_0, r1_0, r2_0, phi_min, phi_max, theta);
	el_p->SetFillStyle(4000);
	// el_p->Draw("same");

	first_area = pi*x_0*y_0;

	total_integral = h2->Integral();
	
	int first_cut_cnt = 0;
	// cout << "first point: " << (h - R1) << ", last point: " << (h + R1) << endl;
	
	for(double x_step = x_0 - r1_0; x_step < x_0 + r1_0; x_step += 0.001){

	  double y_step = y_0 + r2_0*sqrt( 1 - ( pow((x_step-x_0), 2)/(r1_0*r1_0) ) );

	  first_cut_cnt++;
	}

	for(double x_step = x_0 - r1_0; x_step < x_0 + r1_0; x_step += 0.001){

	  double y_step = y_0 - r2_0*sqrt( 1 - ( pow((x_step-x_0), 2)/(r1_0*r1_0) ) );

	  first_cut_cnt++;
	}

	// cout << "mycutg_cnt: " << mycutg_cnt << endl << endl;

	first_el_cut = new TCutG("first_el_cut", first_cut_cnt);
	first_el_cut->SetVarX("x");
	first_el_cut->SetVarY("y");

	first_cut_cnt = 0;

	for(double x_step = x_0 - r1_0; x_step < x_0 + r1_0; x_step += 0.001){

		double y_step = y_0 + r2_0*sqrt( 1 - ( pow((x_step-x_0), 2)/(r1_0*r1_0) ) );
		first_el_cut->SetPoint(first_cut_cnt, x_step, y_step);
		first_cut_cnt++;
	}

	for(double x_step = x_0 - r1_0; x_step < x_0 + r1_0; x_step += 0.001){

	  double y_step = y_0 - r2_0*sqrt( 1 - ( pow((x_step-x_0), 2)/(r1_0*r1_0) ) );
		first_el_cut->SetPoint(first_cut_cnt, x_step, y_step);
		first_cut_cnt++;
	}
	
	first_cut_integral = first_el_cut->IntegralHist(h2);

	// double el_int = ellipse_integral(el_p, h2);
	// cout << "Integral of cut area: " << el_int << endl;

	TMinuit *gMinuit = new TMinuit(4);  //ini TMinuit with a maximum of 5 params
	gMinuit->SetFCN(ellipse_fcn);

	Double_t arglist[10];
	Int_t ierflg = 0;

	arglist[0] = 1;
	gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);

	// Set starting values and step sizes for parameters
	// static Double_t vstart[5] = {-1, -.5 , 1, 0.5 , 0.5};
	h_min = 0;
	h_max = x_0 + dx_sig_p;
	k_min = -0.5;
	k_max = y_0 + dy_sig_p;
	r1_0_min = (2.0/3.0)*r1_0;
	r1_0_max = (4.0/3.0)*r1_0;
	r2_0_min = (2.0/3.0)*r2_0;
	r2_0_max = (4.0/3.0)*r2_0;

	static Double_t vstart[4] = { x_0, y_0, r1_0, r2_0};
	static Double_t step[5] = {0.001 , 0.001 , 0.001 , 0.001, 0.001};
	gMinuit->mnparm(0, "h", vstart[0], step[0], h_min, h_max, ierflg);
	gMinuit->mnparm(1, "k", vstart[1], step[1], k_min, k_max, ierflg);
	gMinuit->mnparm(2, "a", vstart[2], step[2], r1_0_min, r1_0_max,ierflg);
	gMinuit->mnparm(3, "b", vstart[3], step[3], r2_0_min, r2_0_max,ierflg);

// Set error Definition
   gMinuit->SetErrorDef(1);

// Now ready for minimization step
   arglist[0] = 500;
   arglist[1] = 1.;
   gMinuit->mnexcm("MIGRAD", arglist ,1,ierflg);

// Print results
   Double_t amin,edm,errdef;

   Int_t nvpar,nparx,icstat;
   gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
   gMinuit->mnprin(3,amin);


	Double_t x1, y1, r1, r2, ex1, ey1, er1, er2;;
	gMinuit->GetParameter(0, x1, ex1);
	gMinuit->GetParameter(1, y1, ey1);
	gMinuit->GetParameter(2, r1, er1);
	gMinuit->GetParameter(3, r2, er2);

	last_area = TMath::Pi()*(r1)*(r1);

	cout << "Minuit parameters: x1 = " << x1 << ", y1 = " << y1 << ", r1 = " << r1 << ", r2 = " << r2 << endl;


   cout << endl << "--------------------------------------------" << endl;
   cout << "Ellipse limits: " << endl;
   cout << "h: h_min = " << h_min << ", h_start = " << x_0 << ", h_max = " << h_max << endl;
   cout << "k: k_min = " << k_min << ", k_start = " << y_0 << ", k_max = " << k_max<< endl;
   cout << "a: a_min = " << r1_0_min << ", a_start = " << r1_0 << ", a_max = " << r1_0_max << endl;
   cout << "b: b_min = " << r2_0_min << ", b_start = " << r2_0 << ", b_max = " << r2_0_max << endl;
   cout << "--------------------------------------------" << endl << endl;
   cout << "Total Integral (NO CUTS) = " << total_integral << endl;
   cout << "First cut integral = " << first_cut_integral << endl;
   cout << "Final cut integral = " << final_cut_integral << endl;
   cout << "--------------------------------------------" << endl;
   cout << "First area = " << first_area << endl;
   cout << "Last area =" << last_area << endl;

	TEllipse *el_final = new TEllipse(x1, y1, r1, r2, phi_min, phi_max, theta);
	el_final->SetFillStyle(4000);
	el_final->Draw("same");
	el_p->SetFillStyle(4000);
	el_p->SetLineColor(6);
	el_p->Draw("same");



	// el_cut->Draw("same");
}