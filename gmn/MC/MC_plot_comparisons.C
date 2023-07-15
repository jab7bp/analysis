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
int kine = 9;
int sbsfieldscale = 70;
TString run_target = "LD2";
double I_beam = 12.00;
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
TH1D *h_dx_data_select, *h_dx_MC_select;

TH2D *h_dxdy_MC, *h_dxdy_fcut_MC, *h_dxdy_data, *h_dxdy_fcut_data;

double dx_scale, dx_fcut_scale;
double data_scale, mc_scale;

//histo_scaling: 
//	0 -> Integral to 1
//	1 -> Normlize MC to scale of data
//	2 -> Normalized per entry and per unit of X axis

int histo_scaling;
TString plot_cut_type;

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
	if( kine == 9 ){
		MC_filename = Form("MC_SBS%i_%s_mag%i_%suA_dxdy.root", kine, run_target.Data(), sbsfieldscale, I_beam_str.Data() );		
	}
	else{
		MC_filename = Form("MC_SBS%i_%s_mag%imod3_%suA_dxdy.root", kine, run_target.Data(), sbsfieldscale, I_beam_str.Data() );		
	}


	cout << "MC infile: " << MC_filename.Data() << endl << endl;
	MC_inputfile = new TFile(Form("%s/%s", MC_rootfile_dir.Data(), MC_filename.Data()), "READ");

	data_rootfile_dir = "/w/halla-scshelf2102/sbs/jboyd/analysis/gmn/deltaplots/rootfiles";
	//SBS4
	if( kine == 4 ){
		data_filename = Form("%s_SBS%i_mag%i_pqCut_dxdy_parsed_new.root", run_target.Data(), kine, sbsfieldscale);		
	}	

	//SBS8
	if( kine == 8 ){
		data_filename = Form("%s_SBS%i_mag%i_pqCut_dxdy_parsed_new.root", run_target.Data(), kine, sbsfieldscale);		
	}
	//SBS9
	if( kine == 9 ){
		data_filename = Form("%s_SBS%i_mag%i_pqCut_dxdy_parsed_new.root", run_target.Data(), kine, sbsfieldscale);		
	}


	cout << "Data infile: " << data_filename.Data() << endl << endl;
	data_inputfile = new TFile(Form("%s/%s", data_rootfile_dir.Data(), data_filename.Data()), "READ");

	//Grab histograms
	h_dx_MC = static_cast<TH1D*>(MC_inputfile->Get("h_dx"));
	h_dx_fcut_MC = static_cast<TH1D*>(MC_inputfile->Get("h_dx_fcut"));
	h_dx_wcut_MC = static_cast<TH1D*>(MC_inputfile->Get("h_dx_wcut"));
	h_dx_wcut_fcut_MC = static_cast<TH1D*>(MC_inputfile->Get("h_dx_wcut_fcut"));
	h_dxdy_MC = static_cast<TH2D*>(MC_inputfile->Get("h_dxdy"));
	h_dxdy_fcut_MC = static_cast<TH2D*>(MC_inputfile->Get("h_dxdy_fcut"));

	h_dx_data = static_cast<TH1D*>(data_inputfile->Get("h_dx"));
	h_dx_fcut_data = static_cast<TH1D*>(data_inputfile->Get("h_dx_fcut"));
	h_dx_wcut_data = static_cast<TH1D*>(data_inputfile->Get("h_dx_wcut"));
	h_dx_wcut_fcut_data = static_cast<TH1D*>(data_inputfile->Get("h_dx_wcut_fcut"));
	h_dxdy_data = static_cast<TH2D*>(data_inputfile->Get("h_dxdy"));
	h_dxdy_fcut_data = static_cast<TH2D*>(data_inputfile->Get("h_dxdy_fcut"));

	// h_dx_overlay h_dx_fcut_overlay;

// PLOT
	h_dx_data->SetLineColor(4);
	h_dx_wcut_data->SetLineColor(4);
	h_dx_fcut_data->SetLineColor(4);
	h_dx_wcut_fcut_data->SetLineColor(4);

	h_dx_MC->SetLineColor(6);
	h_dx_wcut_MC->SetLineColor(6);
	h_dx_fcut_MC->SetLineColor(6);
	h_dx_wcut_fcut_MC->SetLineColor(6);

	TCanvas *c_dx_overlay = new TCanvas("c_dx_overlay", "c_dx_overlay", 600, 500);
// //Normalization
	histo_scaling = 0;
	vector<TString> histo_scaling_select_vec = {"Integral", "Scale-to-data", "Per entry and bin"};
	TString histo_Scaling_select = histo_scaling_select_vec[histo_scaling];

	plot_cut_type = ""; 

	if( plot_cut_type == "wcut" || plot_cut_type == "fcut" ){
		plot_cut_type.Prepend("_");
	}
	
	h_dx_MC_select = static_cast<TH1D*>(MC_inputfile->Get(Form("h_dx%s", plot_cut_type.Data())));
	Double_t nEntries_MC = h_dx_MC_select->GetEntries();

	h_dx_data_select = static_cast<TH1D*>(data_inputfile->Get(Form("h_dx%s", plot_cut_type.Data())));
	Double_t nEntries_data = h_dx_data_select->GetEntries();

	if( histo_scaling == 0 ){    //	0 -> Integral to 1
		h_dx_MC_select->Scale(1/(h_dx_MC_select->Integral()));
		h_dx_data_select->Scale(1/(h_dx_data_select->Integral()));
		dx_scale = (h_dx_data_select->GetMaximum()/h_dx_MC_select->GetMaximum());
		h_dx_MC_select->Scale(dx_scale);

	}

	if( histo_scaling == 1){ //	1 -> Normlize MC to scale of data
		cout << "Normalizing histograms to scale of Data Histogram" << endl;

		dx_scale = (h_dx_data_select->GetMaximum()/h_dx_MC_select->GetMaximum());
		h_dx_MC_select->Scale(dx_scale);
	}

	if( histo_scaling == 2){ //	2 -> Normalized per entry and per unit of X axis
		cout << "Normalizing histograms per entry and per unit of X axis" << endl;

		h_dx_data_select->Scale( (h_dx_data_select->GetXaxis()->GetBinWidth(1))/(h_dx_data_select->Integral()) );
		h_dx_MC_select->Scale( (h_dx_MC_select->GetXaxis()->GetBinWidth(1))/(h_dx_MC_select->Integral()) );
	
	}

	h_dx_data_select->SetMarkerStyle(8);
	h_dx_data_select->SetMarkerColor(4);
	h_dx_data_select->Draw("hist");	
	// h_dx_data_select->Draw("P + hist");	
	h_dx_MC_select->SetMarkerStyle(8);
	h_dx_MC_select->SetMarkerColor(6);
	h_dx_MC_select->Draw("hist + same");	
	// h_dx_MC_select->Draw("P + hist + same");	


	TLegend *tl_dx_overlay = new TLegend(0.15, 0.70, 0.30, 0.85);
	tl_dx_overlay->AddEntry(h_dx_data_select, "Data");
	tl_dx_overlay->AddEntry(h_dx_MC_select, "Simulation");

	tl_dx_overlay->Draw("same");

	cout << "---------------------------" << endl;
	cout << "Normalizing method: " << histo_Scaling_select.Data() << endl;
	cout << "---------------------------" << endl;
	cout << "MC histo entries: " << int(nEntries_MC) << endl;
	cout << "data histo entries: " << int(nEntries_data) << endl;
	cout << "---------------------------" << endl;

	auto total_time_end = high_resolution_clock::now();
	auto total_time_duration = duration_cast<minutes>(total_time_end - total_time_start);
	cout << "Total time for analysis: " << total_time_duration.count() << " minutes. " << endl;

}