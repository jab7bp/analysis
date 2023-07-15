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

//RUN Info/Parameters
int kine = 8;
int sbsfieldscale = 70;
TString run_target = "LD2";
double I_beam = 5.00;
TString I_beam_str;
double E_beam = lookup_beam_energy_from_kine( kine );

double SBS_field = 1.0*sbsfieldscale/100.0;
double ngen_total = 600000.0;
// double ngen_total = 5000000.0;

TString MC_rootfile_dir, data_rootfile_dir, MC_filename, data_filename;

TFile *MC_inputfile, *data_inputfile;

TH1D *h_dx_MC, *h_dx_fcut_MC, *h_dx_data, *h_dx_fcut_data, *h_dx_wcut_MC, *h_dx_wcut_data;
TH1D *h_dx_wcut_fcut_MC, *h_dx_wcut_fcut_data;
TH1D *h_dx_overlay, *h_dx_fcut_overlay;

TH2D *h_dxdy_MC, *h_dxdy_fcut_MC, *h_dxdy_data, *h_dxdy_fcut_data;

double dx_scale, dx_fcut_scale;
double data_scale, mc_scale;

//histo_scaling: 
//	0 -> Integral to 1
//	1 -> Normlize MC to scale of data
//	2 -> Normalized per entry and per unit of X axis

int histo_scaling = 0;

void MC_plot_comparisons(){

	auto total_time_start = high_resolution_clock::now();
	TStopwatch *StopWatch = new TStopwatch();

	gStyle->SetPalette(55);
	cout << "--------------------------------------" << endl;
	cout << "Analysis started. " << endl;
	cout << "--------------------------------------" << endl;

	I_beam_str = Form("%0.2f", I_beam);
	I_beam_str.ReplaceAll(".", "");

	MC_rootfile_dir = "/w/halla-scshelf2102/sbs/jboyd/analysis/gmn/MC/rootfiles";
	// MC_filename = Form("MC_SBS%i_%s_mag%i_dxdy.root", kine, run_target.Data(), sbsfieldscale );
	MC_filename = Form("MC_SBS%i_%s_mag%i_%suA_dxdy.root", kine, run_target.Data(), sbsfieldscale, I_beam_str.Data() );
	MC_inputfile = new TFile(Form("%s/%s", MC_rootfile_dir.Data(), MC_filename.Data()), "READ");

	data_rootfile_dir = "/w/halla-scshelf2102/sbs/jboyd/analysis/gmn/deltaplots/rootfiles";
	data_filename = Form("%s_SBS%i_mag%i_pqCut_dxdy_parsed.root", run_target.Data(), kine, sbsfieldscale);
	data_inputfile = new TFile(Form("%s/%s", data_rootfile_dir.Data(), data_filename.Data()), "READ");

	//Grab histograms
	h_dx_MC = static_cast<TH1D*>(MC_inputfile->Get("h_dx"));
	h_dx_fcut_MC = static_cast<TH1D*>(MC_inputfile->Get("h_dx_fcut"));
	h_dx_wcut_MC = static_cast<TH1D*>(MC_inputfile->Get("h_dx_wcut"));
	h_dxdy_MC = static_cast<TH2D*>(MC_inputfile->Get("h_dxdy"));
	h_dxdy_fcut_MC = static_cast<TH2D*>(MC_inputfile->Get("h_dxdy_fcut"));

	h_dx_data = static_cast<TH1D*>(data_inputfile->Get("h_dx"));
	h_dx_fcut_data = static_cast<TH1D*>(data_inputfile->Get("h_dx_fcut"));
	h_dx_wcut_data = static_cast<TH1D*>(data_inputfile->Get("h_dx_wcut"));
	h_dx_wcut_fcut_data = static_cast<TH1D*>(data_inputfile->Get("h_dx_wcut_fcut"));
	h_dxdy_data = static_cast<TH2D*>(data_inputfile->Get("h_dxdy"));
	h_dxdy_fcut_data = static_cast<TH2D*>(data_inputfile->Get("h_dxdy_fcut"));


	// h_dx_overlay h_dx_fcut_overlay;

// //Normalization

	dx_scale = (h_dx_data->GetXaxis()->GetBinWidth(1))/(h_dx_data->Integral());

	// dx_scale = (h_dx_data->GetMaximum()/h_dx_MC->GetMaximum());
	dx_scale = (h_dx_wcut_fcut_data->GetMaximum()/h_dx_MC->GetMaximum());
	// dx_scale = (h_dx_fcut_data->GetMaximum()/h_dx_MC->GetMaximum());
	// dx_scale = (h_dx_wcut_data->GetMaximum()/h_dx_MC->GetMaximum());
	// h_dx_MC->Scale(1/(h_dx_MC->Integral()));
	// h_dx_data->Scale(1/(h_dx_data->Integral()));
	
	h_dx_MC->Scale(dx_scale);

	// // h_dx_fcut_MC->Scale(1.0/(h_dx_fcut_MC->GetEntries()));
	// h_dxdy_MC->Scale(1.0/(h_dxdy_MC->GetEntries()));
	// h_dxdy_fcut_MC->Scale(1.0/(h_dxdy_fcut_MC->GetEntries()));

	// // h_dx_MC->Scale(1.0/(h_dx_MC->Integral()));
	// // h_dx_fcut_MC->Scale(1.0/(h_dx_fcut_MC->Integral()));
	// // h_dxdy_MC->Scale(1.0/(h_dxdy_MC->Integral()));
	// // h_dxdy_fcut_MC->Scale(1.0/(h_dxdy_fcut_MC->Integral()));

	// // h_dx_data->Scale(h_dx_MC->Integral()/h_dx_data->Integral());
	// h_dx_fcut_data->Scale(1.0/(h_dx_fcut_data->Integral()));
	// h_dxdy_data->Scale(1.0/(h_dxdy_data->Integral()));
	// h_dxdy_fcut_data->Scale(1.0/(h_dxdy_fcut_data->Integral()));

// PLOT
	TCanvas *c_dx_overlay = new TCanvas("c_dx_overlay", "c_dx_overlay", 600, 500);
	h_dx_data->SetLineColor(4);
	h_dx_data->Draw();
	// h_dx_wcut_fcut_data->SetLineColor(4);
	// h_dx_wcut_fcut_data->Draw();
	// h_dx_fcut_data->SetLineColor(4);
	// h_dx_fcut_data->Draw();
	// h_dx_wcut_data->SetLineColor(4);
	// h_dx_wcut_data->Draw();	

	h_dx_MC->SetLineColor(6);
	h_dx_MC->Draw("hist+same");
	h_dx_data->SetLineColor(4);
	// h_dx_data->Draw("same");

	TLegend *tl_dx_overlay = new TLegend(0.15, 0.70, 0.30, 0.85);
	tl_dx_overlay->AddEntry(h_dx_fcut_data, "Data");
	tl_dx_overlay->AddEntry(h_dx_MC, "Simulation");

	tl_dx_overlay->Draw("same");

	auto total_time_end = high_resolution_clock::now();
	auto total_time_duration = duration_cast<minutes>(total_time_end - total_time_start);
	cout << "Total time for analysis: " << total_time_duration.count() << " minutes. " << endl;

}