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


int pass;
int kine = 9;
int sbsfieldscale = 70;
TString run_target = "LD2";
vector<int> runnum_vec;

TString rootfile_dir;
TString output_dir;
TString outfile_name;
TFile *infile, *parsed_file;

TTree *T, *P, *tree_holder;
vector<TString> input_filenames;

vector<TString> parser_cut_vec;
TString parser_cut_string;

TCut parser_cut = "";

TString experiment = "gmn";

TChain *TC = new TChain("T");
// TEventList *ev_list = new TEventList("ev_list", "Combined Event List");

void parse_gmn_data(){
	auto total_time_start = high_resolution_clock::now();

	cout << "--------------------------------------" << endl;
	cout << "Analysis started. " << endl;
	cout << "--------------------------------------" << endl;
	TTree::SetMaxTreeSize( 1000000000000LL );

	cout << "Run parameters: " << endl;
	cout << "Target: " << run_target << endl;
	cout << "Kinematic: SBS" << kine << endl;
	cout << "SBS Field: " << sbsfieldscale << "%" << endl;
	cout << "--------------------------------------" << endl;

	if( kine == 4 ){
		pass = 0;
	}

	if( kine == 8 ){
		pass = 1;
	}

	if( kine == 9 ){
		pass = 1;
	}

	// rootfile_dir = Form("/lustre19/expphy/volatile/halla/sbs/sbs-gmn/GMN_REPLAYS/pass%i/SBS%i/%s/rootfiles", pass, kine, run_target.Data());
	rootfile_dir = Form("/work/halla/sbs/sbs-gmn/pass%i/SBS%i/%s/rootfiles/", pass, kine, run_target.Data());
	output_dir = "/lustre19/expphy/volatile/halla/sbs/jboyd/analysis_rootfiles/jboyd_parsed";

	outfile_name = Form("gmn_parsed_%s_SBS%i_mag%i.root", run_target.Data(), kine, sbsfieldscale);
	parsed_file = new TFile(Form("%s/%s", output_dir.Data(), outfile_name.Data()), "RECREATE");

	cout << "Building runnum vector. " << endl;
	for(int i = 0; i < lookup_parsed_runs_cnt(run_target.Data(), kine, sbsfieldscale); i++){
		runnum_vec.push_back(lookup_parsed_runnums(run_target.Data(), kine, sbsfieldscale, i));
	}
	cout << "--------------------------------------" << endl;
	cout << "Number of runs: " << runnum_vec.size() << endl;
	cout << "Runs: " << runnum_vec[0] << " thru " << runnum_vec[runnum_vec.size() - 1] << endl;
	cout << "--------------------------------------" << endl;

	cout << "Adding rootfiles to TChain." << endl;
	for(size_t run = 0; run < runnum_vec.size(); run++){
		TC->Add(Form("%s/*%i*", rootfile_dir.Data(), runnum_vec[run]));
	}

	if( kine == 4 ){
		parser_cut_vec = {
			"sbs.hcal.nclus>0",
			"bb.ps.nclus>0",
			"bb.sh.nclus>0",
			"abs(bb.tr.vz)<0.08",
			"bb.gem.track.nhits[0]>3",
			"bb.tr.n==1",
			Form("bb.ps.e>%f", 0.135), //was 0.140
			Form("((bb.sh.e+bb.ps.e)/(bb.tr.p[0]))>(%f)&&((bb.sh.e+bb.ps.e)/(bb.tr.p[0]))<(%f)", 0.65, 1.30), //was 0.75, 12.0
			// Form("((abs(((bb.sh.e+bb.ps.e)/(bb.tr.p[0]))))-%f)<%f", lookup_cut(runnum, "Ep"), lookup_cut(runnum, "Ep_sigma")),
			Form("sbs.hcal.e>%f", 0.02), //was 0.025
			Form("(bb.sh.e+bb.ps.e)>%f", 1.55), //was 1.5
			Form("((e.kine.W2)>%f)&&((e.kine.W2)<%f)", 0.55, 1.25) //was 0.60, 1.20
		};
	}

	if( kine == 8 ){
		parser_cut_vec = {
			"sbs.hcal.nclus>0",
			"bb.ps.nclus>0",
			"bb.sh.nclus>0",
			"abs(bb.tr.vz)<0.08",
			"bb.gem.track.nhits[0]>3",
			"bb.tr.n==1",
			Form("bb.ps.e>%f", 0.175),
			Form("((bb.sh.e+bb.ps.e)/(bb.tr.p[0]))>(%f)&&((bb.sh.e+bb.ps.e)/(bb.tr.p[0]))<(%f)", 0.75, 1.25),
			// Form("((abs(((bb.sh.e+bb.ps.e)/(bb.tr.p[0]))))-%f)<%f", lookup_cut(runnum, "Ep"), lookup_cut(runnum, "Ep_sigma")),
			Form("sbs.hcal.e>%f", 0.025),
			Form("(bb.sh.e+bb.ps.e)>%f", 2.1),
			Form("((e.kine.W2)>%f)&&((e.kine.W2)<%f)", 0.6, 1.30)
		};
	}

	if( kine == 9 ){
		parser_cut_vec = {
			"sbs.hcal.nclus>0",
			"bb.ps.nclus>0",
			"bb.sh.nclus>0",
			"abs(bb.tr.vz)<0.08",
			"bb.gem.track.nhits[0]>3",
			"bb.tr.n==1",
			Form("bb.ps.e>%f", lookup_pre_parsed_cut( run_target.Data(), kine, "PS_min")),
			Form("((bb.sh.e+bb.ps.e)/(bb.tr.p[0]))>(%f)&&((bb.sh.e+bb.ps.e)/(bb.tr.p[0]))<(%f)", lookup_pre_parsed_cut( run_target.Data(), kine, "Ep_min"), lookup_pre_parsed_cut( run_target.Data(), kine, "Ep_max")),
			// Form("((abs(((bb.sh.e+bb.ps.e)/(bb.tr.p[0]))))-%f)<%f", lookup_cut(runnum, "Ep"), lookup_cut(runnum, "Ep_sigma")),
			Form("sbs.hcal.e>%f", lookup_pre_parsed_cut( run_target.Data(), kine, "HCal_min")),
			Form("(bb.sh.e+bb.ps.e)>%f", lookup_pre_parsed_cut( run_target.Data(), kine, "SH_PS_min")),
			Form("((e.kine.W2)>%f)&&((e.kine.W2)<%f)", lookup_pre_parsed_cut( run_target.Data(), kine, "W2_min"), lookup_pre_parsed_cut( run_target.Data(), kine, "W2_max"))


			// Form("bb.ps.e>%f", 0.21),
			// Form("((bb.sh.e+bb.ps.e)/(bb.tr.p[0]))>(%f)&&((bb.sh.e+bb.ps.e)/(bb.tr.p[0]))<(%f)", 0.58, 1.58),
			// // Form("((abs(((bb.sh.e+bb.ps.e)/(bb.tr.p[0]))))-%f)<%f", lookup_cut(runnum, "Ep"), lookup_cut(runnum, "Ep_sigma")),
			// Form("sbs.hcal.e>%f", .03),
			// Form("(bb.sh.e+bb.ps.e)>%f", 1.2),
			// Form("((e.kine.W2)>%f)&&((e.kine.W2)<%f)", 0.6, 1.30)
		};
	}
	

	for(size_t cut = 0; cut < parser_cut_vec.size(); cut++){
		if(cut == parser_cut_vec.size() - 1){
			parser_cut_string.Append(Form("%s", parser_cut_vec[cut].Data()));
		}
		else{
			parser_cut_string.Append(Form("%s%s", parser_cut_vec[cut].Data(), "&&"));
		}
	}

	parser_cut = Form("%s", parser_cut_string.Data());

	cout << "Applied cut: " << endl;
	cout << parser_cut << endl;
	cout << "--------------------" << endl << endl;

	cout << "Creating parsed tree, P, from TChain TC..." << endl;
	P = (TTree*)TC->CopyTree(parser_cut_string.Data());

	cout << "--------------------------------------" << endl;
	cout << "Finished parsing tree. " << endl;
	cout << "--------------------------------------" << endl;
	cout << "Target: " << run_target.Data() << endl;
	cout << "SBS" << kine << endl;
	cout << "Magnet: " << sbsfieldscale << "%" << endl;
	cout << "Entries in original tree: " << TC->GetEntries() << endl;
	cout << "Entries in parsed tree: " << P->GetEntries() << endl;
	cout << "--------------------------------------" << endl;
	cout << "--------------------------------------" << endl;
	cout << "Writing tree to file and saving to: " << endl;
	cout << parsed_file->GetName() << endl;
	parsed_file->Write();
	cout << "--------------------------------------" << endl;
	cout << "Finished writing to file. " << endl;
	cout << "--------------------------------------" << endl;
	cout  << "-------------------------" << endl;
	// cout << "Entries in trees: " << endl;
	// cout << "TC: " << TC->GetEntries() << endl;
	// cout << "tree_holder: " << tree_holder->GetEntries() << endl;


	auto total_time_end = high_resolution_clock::now();
	auto total_time_duration = duration_cast<minutes>(total_time_end - total_time_start);
	cout << "Total time for analysis: " << total_time_duration.count() << " minutes. " << endl;
}