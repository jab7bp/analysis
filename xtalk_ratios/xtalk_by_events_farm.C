#include <iostream>
#include <fstream>
#include <numeric>
#include <iomanip>
#include <math.h>
#include <algorithm>
#include <string>
#include <chrono>
#include <dirent.h>

using namespace std;
using namespace std::chrono;

#include "/w/halla-scshelf2102/sbs/jboyd/xtalk/include/include_files.h"
#include "/w/halla-scshelf2102/sbs/jboyd/xtalk/include/GEM_lookups.h"


//BOOLEANS
bool build_all = true;
bool build_last_mod_only = false;
bool build_custom = false;

bool build_U_strips = true;
bool build_V_strips = true;

int ADCcut = 0;
int ADC_min_thresh = 0;
int ADC_max_thresh = 7500;

int first_event = 0;
int last_event;

int nLayers = 5;
int first_layer = 0;
int last_layer = 5;

int first_module = 0;
int last_module;

//Custom run settings
int first_layer_custom = 0;
int last_layer_custom = 5;
int first_module_custom = 0;
int last_module_custom = -1;
int first_apv_custom = 0;
int last_apv_custom = 5;

int module_num_layer;
int num_modules_layer;

int module_num_global;
int num_modules_total;

int gemType;
TString gemType_name[3] = {"UVa UV", "UVa XY", "INFN XY"};

int first_apv;
int last_apv_U;
int last_apv_V;

int nAPVs_U;
int nAPVs_V;
int max_apvs = 30;

int max_strips = 16000;

int neighbor_ratio_cnt_U = 0;
int neighbor_ratio_cnt_V = 0;

//Root input/output file variables
TString IN_DIR = gSystem->Getenv("$SWIF_JOB_WORK_DIR");
TString OUT_DIR = gSystem->Getenv("$SWIF_JOB_WORK_DIR");

TChain *TC = new TChain("T");

TString inputfile;

TFile *RatioRootfile;
TString rootfilename;

TString protorootfile;

//VECTORS TO HOLD ADC VALUES
vector<vector<Double32_t>> APV_chans_ADCsum_U;
vector<vector<Double32_t>> APV_chans_isampmax_U;

vector<vector<Double32_t>> APV_chans_ADCsum_V;
vector<vector<Double32_t>> APV_chans_isampmax_V;

// //ARRAYS TO HOLD ADCSUM VALUES
// Double32_t ***APV_chans_ADCsum_U;
// Int_t ***APV_chans_isampmax_U;

// Double32_t ***APV_chans_ADCsum_V;
// Int_t ***APV_chans_isampmax_V;

int nEntries;
int ratio_bin_cnt = 300;
int ratio_bin_min = 0;
int ratio_bin_max = 30;

//@@@@@@@@@@@@@@@@
//BRANCHES AND VARIABLES
//@@@@@@@@@@@@@@@@

//------------------
//---- U-Strips ----
//------------------
TTree *tr_ratios_U;

Int_t layer_U[5];
Int_t module_U[5][4];

TBranch *br_layer_U;
TBranch *br_module_U;

//Branch Variables are 3D arrays structured as [Layers][Modules on Layer][APVs]
Double_t ratio_neighbor_channels_ADCsum_U[5][4][30];
Double_t larger_neighbor_ADCsum_U[5][4][30];
Double_t smaller_neighbor_ADCsum_U[5][4][30];
Int_t larger_neighbor_isampmax_U[5][4][30];
Int_t smaller_neighbor_isampmax_U[5][4][30];

TBranch *br_ratio_neighbor_channels_ADCsum_U[5][4][30];
TBranch *br_larger_neighbor_ADCsum_U[5][4][30];
TBranch *br_smaller_neighbor_ADCsum_U[5][4][30];
TBranch *br_larger_neighbor_isampmax_U[5][4][30];
TBranch *br_smaller_neighbor_isampmax_U[5][4][30];

//------------------
//---- V-Strips ----
//------------------

TTree *tr_ratios_V;

Int_t layer_V[5];
Int_t module_V[5][4];

TBranch *br_layer_V;
TBranch *br_module_V;

//Branches are 3D arrays structured as [Layers][Modules on Layer][APVs]
Double_t ratio_neighbor_channels_ADCsum_V[5][4][30];
Double_t larger_neighbor_ADCsum_V[5][4][30];
Double_t smaller_neighbor_ADCsum_V[5][4][30];
Int_t larger_neighbor_isampmax_V[5][4][30];
Int_t smaller_neighbor_isampmax_V[5][4][30];

TBranch *br_ratio_neighbor_channels_ADCsum_V[5][4][30];
TBranch *br_larger_neighbor_ADCsum_V[5][4][30];
TBranch *br_smaller_neighbor_ADCsum_V[5][4][30];
TBranch *br_larger_neighbor_isampmax_V[5][4][30];
TBranch *br_smaller_neighbor_isampmax_V[5][4][30];

void xtalk_by_events_farm(int runnum = 13770, int last_event_select = 5000, int segment = 0, int build_select = 0){

	auto total_time_start = high_resolution_clock::now();

	if( build_select == 0 ){
		build_all = true;
		build_last_mod_only = false;
		build_custom = false;
	}
	if( build_select == 1 ){
		build_all = false;
		build_last_mod_only = true;
		build_custom = false;
	}
	
	TTree::SetMaxTreeSize( 1000000000000LL );
	num_modules_total = lookup_nModules_total(runnum);

	cout << "Crosstalk analysis started. " << endl;

	cout << "--------------------------------------------------" << endl;
	cout << "Run number: " << runnum << ", ADCcut: " << ADCcut << endl;

	if( build_all ){
		cout << "****************************************************************" << endl;
		cout << "**************Building all layers and all modules.**************" << endl;		
		cout << "First Layer: " << first_layer << ", Last layer: " << last_layer << endl;
		cout << "****************************************************************" << endl;
	}

	if( build_last_mod_only && !build_all ){
		cout << "--------------------------------------------------" << endl;
		cout << "Building only the last module of the configuration." << endl;
		cout << "Layer: " << nLayers - 1 << ", Module: " << lookup_nModules_layer(runnum, (nLayers - 1)) - 1<< endl;
		cout << "--------------------------------------------------" << endl;
	}

	if( build_custom && !build_all && !build_last_mod_only ){
		cout << "--------------------------------------------------" << endl;
		cout << "----------RUNNING CUSTOM SELECTIONS -------------" << endl;
		cout << "First Layer: " << first_layer << ", Last layer: " << last_layer << endl;
		cout << "--------------------------------------------------" << endl;
	}

	cout << "Building filename and checking its existence..." << endl;
	
	// IN_DIR.Form("/lustre/expphy/volatile/halla/sbs/jboyd/Rootfiles/xtalk/histo_with_corr/%i", runnum);	
	// OUT_DIR.Form("/lustre19/expphy/volatile/halla/sbs/jboyd/Rootfiles/xtalk/xtalk_by_events/%i", runnum);

	inputfile.Form("%se1209019_fullreplay_%i_*", IN_DIR.Data(), runnum);

	if( build_all ){
		rootfilename.Form("%s%i_xtalk_ratios_%i_thru%i_events_U%i_V%i_seg_%i_ALL.root", OUT_DIR.Data(), runnum, first_event, last_event, build_U_strips, build_V_strips, segment);
	}
	if( build_last_mod_only && !build_all ){
		rootfilename.Form("%s%i_xtalk_ratios_%i_thru%i_events_U%i_V%i_seg_%i_last_module_only.root", OUT_DIR.Data(), runnum, first_event, last_event, build_U_strips, build_V_strips, segment);
	}
	if( build_custom && !build_all && !build_last_mod_only ){
		rootfilename.Form("%s%i_xtalk_ratios_%i_thru%i_events_U%i_V%i_seg_%i.root", OUT_DIR.Data(), runnum, first_event, last_event, build_U_strips, build_V_strips, segment);
	}

	cout << endl << "Input file: " << inputfile.Data() << endl << endl;
	cout << endl << "Output file: " << rootfilename.Data() << endl << endl;

	cout << "Adding following file(s) to TChain: " << inputfile.Data() << endl;

	TC->Add( inputfile.Data() );

	cout << endl << "Calculating total entries..." << endl << endl;
	nEntries = 0;
	nEntries = TC->GetEntries();

	std::cout.flush();

	cout << endl << "**************************************************" << endl;
	cout << "Number of entries: " << nEntries << endl;
	cout << "**************************************************" << endl;

	if( last_event_select == -1 ){
		cout << "Configured to analyze ALL events. " << endl;
		last_event = nEntries;
	}
	else{ 
		cout << "Only analyzing user-defined number of events." << endl;
		last_event = last_event_select;
	}
	cout << "Number of events to analyze: " << last_event << endl;
	cout << "**************************************************" << endl;


//CREATE ROOTFILES AND TREES
	RatioRootfile = new TFile(Form( "%s", rootfilename.Data() ), "RECREATE" );

	if( build_U_strips ){ tr_ratios_U = new TTree("tr_ratios_U", "tr_ratios_U"); }
	if( build_V_strips ){ tr_ratios_V = new TTree("tr_ratios_V", "tr_ratios_V"); }

	cout << "Starting main loop...." << endl;

	if( build_all ){
		first_layer = 0;
		last_layer = nLayers;
	}
	if( build_last_mod_only && !build_all ){
		first_layer = (nLayers - 1);
		last_layer = nLayers;
	}
	if( build_custom && !build_all & !build_last_mod_only ){
		first_layer = first_layer_custom;
		last_layer = last_layer_custom;
	}

	for(int layer = first_layer; layer < last_layer; layer++){

		num_modules_layer = lookup_nModules_layer(runnum, layer);

//Initilaize some things...
		layer_U[layer] = 0;
		layer_V[layer] = 0;

//START OF MODULES LOOP
		if( build_all ){
			first_module = 0;
			last_module = num_modules_layer;
		}
		if( build_last_mod_only && !build_all ){
			first_module = lookup_nModules_layer(13770, (nLayers - 1)) - 1;
			last_module = lookup_nModules_layer(13770, (nLayers - 1));
		}
		if( build_custom && !build_all && !build_last_mod_only ){
			first_module = first_module_custom;
			if( last_module_custom == -1 ){
				last_module = num_modules_layer;
			}
			else{
				last_module = last_module_custom;
			}
		}

		for( int module_on_layer = first_module; module_on_layer < last_module; module_on_layer++ ){
			cout << "*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*" << endl;
			cout << "Working on Layer " << layer << "(Last layer: " << last_layer << ")" << endl;
//Initilaize some things...
			module_U[layer][module_on_layer] = 0;
			module_V[layer][module_on_layer] = 0;

			for( int apv_init = 0; apv_init < max_apvs; apv_init++ ){
				ratio_neighbor_channels_ADCsum_U[layer][module_on_layer][apv_init] = 0.0;
				larger_neighbor_ADCsum_U[layer][module_on_layer][apv_init] = 0.0;

				ratio_neighbor_channels_ADCsum_V[layer][module_on_layer][apv_init] = 0.0;
				larger_neighbor_ADCsum_V[layer][module_on_layer][apv_init] = 0.0;
			}
			gemType = lookup_GEM_type(runnum, layer);
			module_num_global = lookup_global_mod_num(runnum, layer, module_on_layer);
			nAPVs_U = lookup_nAPVs(gemType, 0);
			nAPVs_V = lookup_nAPVs(gemType, 1);
			cout <<  "Module: " << module_on_layer << " (Global mod num: " << module_num_global << ") " << endl;
			cout << "Config: " << lookup_config(runnum) << ", Gem type: " << gemType << ", nAPVs_U: " << nAPVs_U << ", nAPVs_V: " << nAPVs_V << endl;
			cout << "*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*" << endl;

//RESIZE VECTORS FOR ADC VALUES:
			if( build_U_strips ){
				APV_chans_ADCsum_U.resize( nAPVs_U, vector<Double32_t>(128, 0) );
				APV_chans_isampmax_U.resize( nAPVs_U, vector<Double32_t>(128, 0) );
			}
			if( build_V_strips ){
				APV_chans_ADCsum_V.resize( nAPVs_V, vector<Double32_t>(128, 0) );
				APV_chans_isampmax_V.resize( nAPVs_V, vector<Double32_t>(128, 0) );
			}

			if( build_U_strips ){
					// cout << "Clearing U " << endl;
					for( int apv_clr = 0; apv_clr < nAPVs_U; apv_clr++ ){
						for( int chan_clr = 0; chan_clr < 128; chan_clr++ ){
							APV_chans_ADCsum_U[apv_clr][chan_clr] = 0.0;
							APV_chans_isampmax_U[apv_clr][chan_clr] = 0;
						}
					}
				}
			if( build_V_strips ){
				// cout << "Clearing V " << endl;
				for( int apv_clr = 0; apv_clr < nAPVs_V; apv_clr++ ){
					for( int chan_clr = 0; chan_clr < 128; chan_clr++ ){
						APV_chans_ADCsum_V[apv_clr][chan_clr] = 0.0;
						APV_chans_isampmax_V[apv_clr][chan_clr] = 0;
					}
				}
			}

			if( build_custom ){
				first_apv = first_apv_custom;
				last_apv_U = last_apv_custom;
				last_apv_V = last_apv_custom;
			}
			else{
				first_apv = 0;
				last_apv_U = nAPVs_U;
				last_apv_V = nAPVs_V;
			}

//DEFINE THE BRANCHES FOR EACH LAYER, MODULE, APV.....
			if( build_U_strips ){
				for(int apv_init = first_apv; apv_init < last_apv_U; apv_init++ ){
					br_layer_U = tr_ratios_U->Branch( Form( "Layer_U_%i", layer), &layer_U[layer], Form("Layer_U_%i/I", layer) );
					br_module_U = tr_ratios_U->Branch( Form( "Layer%i_Module%i_U", layer, module_on_layer), &module_U[layer][module_on_layer], Form("Layer%i_Module%i_U/I", layer, module_on_layer) );
					br_ratio_neighbor_channels_ADCsum_U[layer][module_on_layer][apv_init] = tr_ratios_U->Branch( Form( "Layer%i_Module%i_APV%i_ratio_neighbor_channels_ADCsum_U", layer, module_on_layer, apv_init), &ratio_neighbor_channels_ADCsum_U[layer][module_on_layer][apv_init], Form("Layer%i_Module%i_APV%i_ratio_neighbor_channels_ADCsum_U/D", layer, module_on_layer, apv_init) );
					br_larger_neighbor_ADCsum_U[layer][module_on_layer][apv_init] = tr_ratios_U->Branch( Form("Layer%i_Module%i_APV%i_larger_neighbor_ADCsum_U", layer, module_on_layer, apv_init), &larger_neighbor_ADCsum_U[layer][module_on_layer][apv_init], Form("Layer%i_Module%i_APV%i_larger_neigh_ADCsum_U/D", layer, module_on_layer, apv_init) );
					br_smaller_neighbor_ADCsum_U[layer][module_on_layer][apv_init] = tr_ratios_U->Branch( Form("Layer%i_Module%i_APV%i_smaller_neighbor_ADCsum_U", layer, module_on_layer, apv_init), &smaller_neighbor_ADCsum_U[layer][module_on_layer][apv_init], Form("Layer%i_Module%i_APV%i_smaller_neigh_ADCsum_U/D", layer, module_on_layer, apv_init) );
					br_larger_neighbor_isampmax_U[layer][module_on_layer][apv_init] = tr_ratios_U->Branch( Form("Layer%i_Module%i_APV%i_larger_neighbor_isampmax_U", layer, module_on_layer, apv_init), &larger_neighbor_isampmax_U[layer][module_on_layer][apv_init], Form("Layer%i_Module%i_APV%i_larger_neighbor_isampmax_U/D", layer, module_on_layer, apv_init) );
					br_smaller_neighbor_isampmax_U[layer][module_on_layer][apv_init] = tr_ratios_U->Branch( Form("Layer%i_Module%i_APV%i_smaller_neighbor_isampmax_U", layer, module_on_layer, apv_init), &smaller_neighbor_isampmax_U[layer][module_on_layer][apv_init], Form("Layer%i_Module%i_APV%i_smaller_neighbor_isampmax_U/D", layer, module_on_layer, apv_init) );
				}
			}
			if( build_V_strips ){
				for(int apv_init = first_apv; apv_init < last_apv_V; apv_init++ ){
					br_layer_V = tr_ratios_V->Branch( Form( "Layer_V_%i", layer), &layer_V[layer], Form("Layer_V_%i/I", layer) );
					br_module_V = tr_ratios_V->Branch( Form( "Layer%i_Module%i_V", layer, module_on_layer), &module_V[layer][module_on_layer], Form("Layer%i_Module%i_V/I", layer, module_on_layer) );
					br_ratio_neighbor_channels_ADCsum_V[layer][module_on_layer][apv_init] = tr_ratios_V->Branch( Form( "Layer%i_Module%i_APV%i_ratio_neighbor_channels_ADCsum_V", layer, module_on_layer, apv_init), &ratio_neighbor_channels_ADCsum_V[layer][module_on_layer][apv_init], Form("Layer%i_Module%i_APV%i_ratio_neighbor_channels_ADCsum_V/D", layer, module_on_layer, apv_init) );
					br_larger_neighbor_ADCsum_V[layer][module_on_layer][apv_init] = tr_ratios_V->Branch( Form("Layer%i_Module%i_APV%i_larger_neighbor_ADCsum_V", layer, module_on_layer, apv_init), &larger_neighbor_ADCsum_V[layer][module_on_layer][apv_init], Form("Layer%i_Module%i_APV%i_larger_neigh_ADCsum_V/D", layer, module_on_layer, apv_init) );
					br_smaller_neighbor_ADCsum_V[layer][module_on_layer][apv_init] = tr_ratios_V->Branch( Form("Layer%i_Module%i_APV%i_smaller_neighbor_ADCsum_V", layer, module_on_layer, apv_init), &smaller_neighbor_ADCsum_V[layer][module_on_layer][apv_init], Form("Layer%i_Module%i_APV%i_smaller_neigh_ADCsum_V/D", layer, module_on_layer, apv_init) );
					br_larger_neighbor_isampmax_V[layer][module_on_layer][apv_init] = tr_ratios_V->Branch( Form("Layer%i_Module%i_APV%i_larger_neighbor_isampmax_V", layer, module_on_layer, apv_init), &larger_neighbor_isampmax_V[layer][module_on_layer][apv_init], Form("Layer%i_Module%i_APV%i_larger_neighbor_isampmax_V/D", layer, module_on_layer, apv_init) );
					br_smaller_neighbor_isampmax_V[layer][module_on_layer][apv_init] = tr_ratios_V->Branch( Form("Layer%i_Module%i_APV%i_smaller_neighbor_isampmax_V", layer, module_on_layer, apv_init), &smaller_neighbor_isampmax_V[layer][module_on_layer][apv_init], Form("Layer%i_Module%i_APV%i_smaller_neighbor_isampmax_V/D", layer, module_on_layer, apv_init) );
				}
			}

			cout << "Setting up branches & variables." << endl;
//Turn OFF/ON all Branches
			TC->SetBranchStatus("*", false);

			//NData Branches
			cout << "Branches for Layer: " << layer <<"; Module on layer: " << module_on_layer << "; Global module: " << module_num_global << endl;
			TC->SetBranchStatus(Form("Ndata.bb.gem.m%i.strip.IsU", module_num_global), true);
			TC->SetBranchStatus(Form("Ndata.bb.gem.m%i.strip.IsV", module_num_global), true);

			//Variable Branches
			TC->SetBranchStatus(Form("bb.gem.m%i.strip.istrip", module_num_global), true);
			TC->SetBranchStatus(Form("bb.gem.m%i.strip.ADCsum", module_num_global), true);
			TC->SetBranchStatus(Form("bb.gem.m%i.strip.isampmax", module_num_global), true);
			TC->SetBranchStatus(Form("bb.gem.m%i.strip.IsU", module_num_global), true);
			TC->SetBranchStatus(Form("bb.gem.m%i.strip.IsV", module_num_global), true);

//Assign branches to variables
			Int_t NData_IsU;
			Int_t NData_IsV;
			Double_t IsU[max_strips];
			Double_t IsV[max_strips];
			Double_t isampmax[max_strips];
			Double_t istrip[max_strips];
			Double_t ADCsum[max_strips];

			//NData
			TC->SetBranchAddress(Form("Ndata.bb.gem.m%i.strip.IsU", module_num_global), &NData_IsU);
			TC->SetBranchAddress(Form("Ndata.bb.gem.m%i.strip.IsV", module_num_global), &NData_IsV);

			//Strip Variable Data
			TC->SetBranchAddress(Form("bb.gem.m%i.strip.istrip", module_num_global), &istrip);
			TC->SetBranchAddress(Form("bb.gem.m%i.strip.ADCsum", module_num_global), &ADCsum);
			TC->SetBranchAddress(Form("bb.gem.m%i.strip.isampmax", module_num_global), &isampmax);
			TC->SetBranchAddress(Form("bb.gem.m%i.strip.IsU", module_num_global), &(IsU) );
			TC->SetBranchAddress(Form("bb.gem.m%i.strip.IsV", module_num_global), &(IsV) );

			cout << "Finished setting up branches. " << endl << endl;
			cout << "Checking vector sizes for nAPVs_U = " << nAPVs_U << " and nAPVs_V = " << nAPVs_V << endl;
			if( build_U_strips ){
				cout << "ADCsum_U DIM 1: " << APV_chans_ADCsum_U.size() << ", ADCsum_U DIM 2: " << APV_chans_ADCsum_U[0].size() << endl;
				if( (APV_chans_ADCsum_U.size() == abs(nAPVs_U) ) && (APV_chans_ADCsum_U[0].size() == 128) ){ 
					cout << "ADCsum_U Dimension 1 & 2 check passed. " << endl;
				}
				else{
					cout << "ADCsum_U Dimension 1 & 2 check failed. Check the nAPVs_U and/or chan declaration. Exiting." << endl;
					exit(0);
				}
				
			}
			if( build_V_strips ){
				cout << "ADCsum_V DIM 1: " << APV_chans_ADCsum_V.size() << ", ADCsum_V DIM 2: " << APV_chans_ADCsum_V[0].size() << endl;
				if( (APV_chans_ADCsum_V.size() == abs(nAPVs_V) ) && (APV_chans_ADCsum_V[0].size() == 128) ){ 
					cout << "ADCsum_V Dimension 1 & 2 check passed. " << endl;
				}
				else{
					cout << "ADCsum_V Dimension 1 & 2 check failed. Check the nAPVs_V and/or chan declaration. Exiting." << endl;
					exit(0);
				}
			}
			cout << "All vectors passed checks." << endl;
			cout << "-------------------------------------------------------------------------------" << endl;
			cout << "Looping through events, pulling variables, and placing them into arrays...." << endl;

//*********************EVENTS LOOP*************************
//*********************EVENTS LOOP*************************
//*********************EVENTS LOOP*************************

			for( int evt = 0; evt < last_event; evt++ ){
				if( build_all ){
					if( evt%10000 == 0 ){
						cout << "Pulling layer " << layer << " - Module " << module_on_layer << " (Global: " << module_num_global << ") - Event: " << evt << " of " << last_event << " events (" << (int(100*(double(evt)/(double(last_event))))) << "/100%" << " for module)(" << int(100*(( double(((module_num_global)*nEntries) + evt) )/( double(nEntries*num_modules_total) ))) << "/100%" << " in total)" << endl;
					}
				}
				if( last_event <= 500000 ){
					if( evt%5000 == 0 ){
						cout << "Pulling layer " << layer << " - Module " << module_on_layer << " (Global: " << module_num_global << ") - Event: " << evt << " of " << last_event << " events (" << (int(100*(double(evt)/(double(last_event))))) << "/100%" << " for module)(" << (((layer)*last_event) + int(100*(double(evt)/(last_layer*(double(last_event)))))) << "/100%" << " in total)" << endl;
					}
				}
				if( build_last_mod_only ){
					if( evt%10000 == 0 ){
						cout << "Pulling layer " << layer << " - Module " << module_on_layer << " (Global: " << module_num_global << ") - Event: " << evt << " of " << last_event << " events (" << (int(100*(double(evt)/(double(last_event))))) << "/100%" << " in total)" << endl;
					}
				}
				
				TC->GetEntry(evt);

		//Initialize/Clear arrays here to prevent any issues with left-over values
				if( build_U_strips ){
					// cout << "Clearing U " << endl;
					for( int apv_clr = 0; apv_clr < nAPVs_U; apv_clr++ ){
						for( int chan_clr = 0; chan_clr < 128; chan_clr++ ){
							APV_chans_ADCsum_U[apv_clr][chan_clr] = 0.0;
							APV_chans_isampmax_U[apv_clr][chan_clr] = 0;
						}
					}
				}
				if( build_V_strips ){
					// cout << "Clearing V " << endl;
					for( int apv_clr = 0; apv_clr < nAPVs_V; apv_clr++ ){
						for( int chan_clr = 0; chan_clr < 128; chan_clr++ ){
							APV_chans_ADCsum_V[apv_clr][chan_clr] = 0.0;
							APV_chans_isampmax_V[apv_clr][chan_clr] = 0;
						}
					}
				}

				if( build_U_strips ){
					// cout << "Pulling U " << endl;
					for( int i = 0; i < NData_IsU; i++ ){

						if( IsU[i] == 1 ){
							int apv_U = strip_to_APV( int(istrip[i]) );
							int chan_U = GEM_strip_to_channel( gemType, int(istrip[i])%128 );

							APV_chans_ADCsum_U[apv_U][chan_U] = ADCsum[i];
							APV_chans_isampmax_U[apv_U][chan_U] = isampmax[i];
						}
					}
				}
				if( build_V_strips ){
					// cout << "Pulling V " << endl;
					for( int i = 0; i < NData_IsV; i++ ){
						if( IsV[i] == 1 ){
							int apv_V = strip_to_APV( int(istrip[i]) );
							int chan_V = GEM_strip_to_channel( gemType, int(istrip[i])%128 );

							APV_chans_ADCsum_V[apv_V][chan_V] = ADCsum[i];
							APV_chans_isampmax_V[apv_V][chan_V] = isampmax[i];
						}
					}
				}

//---------------------------------------------------------------------------------
//					USE THE VARIABLES AND CALCULATE RATIOS
//---------------------------------------------------------------------------------		
				std::cout.flush();

	//U-Strips
				if( build_U_strips ){
					// cout << " U ADC " << endl;
					double neighbor_ADCsum_ratio_U = -99999;
					double chan_i_ADCsum_U = 0.0;
					double chan_j_ADCsum_U = 0.0;
					int chan_i_isampmax_U = -1;
					int chan_j_isampmax_U = -1;

					for( int apv_U = first_apv; apv_U < last_apv_U; apv_U++ ){
						for( int chan_U = 0; chan_U < 127; chan_U++ ){
							// cout << "apv: " << apv_U << ", chan: " << chan_U << endl;
							//Initialize variables for U-Strip ADC values per channel
							neighbor_ADCsum_ratio_U = -99999;
							chan_i_ADCsum_U = 0.0;
							chan_j_ADCsum_U = 0.0;
							chan_i_isampmax_U = -1;
							chan_j_isampmax_U = -1;

							//STORE
							chan_i_ADCsum_U = (1.0)*APV_chans_ADCsum_U[apv_U][chan_U];
							chan_j_ADCsum_U = (1.0)*APV_chans_ADCsum_U[apv_U][chan_U + 1];
							chan_i_isampmax_U = APV_chans_isampmax_U[apv_U][chan_U];
							chan_j_isampmax_U = APV_chans_isampmax_U[apv_U][chan_U + 1];

				//Channel I is larger than Channel J:
				//Chan i & j DO NOT equal 0 AND they have matching time sample indices
							if( (chan_i_ADCsum_U > chan_j_ADCsum_U) && (chan_i_ADCsum_U > ADCcut) && (chan_i_ADCsum_U != 0) && (chan_j_ADCsum_U != 0) && (chan_i_isampmax_U == chan_j_isampmax_U) ){
								neighbor_ADCsum_ratio_U = ( (1.0)*chan_i_ADCsum_U )/( (1.0)*chan_j_ADCsum_U );

								ratio_neighbor_channels_ADCsum_U[layer][module_on_layer][apv_U] = neighbor_ADCsum_ratio_U;
								larger_neighbor_ADCsum_U[layer][module_on_layer][apv_U] = chan_i_ADCsum_U;
								smaller_neighbor_ADCsum_U[layer][module_on_layer][apv_U] = chan_j_ADCsum_U;
								larger_neighbor_isampmax_U[layer][module_on_layer][apv_U] = chan_i_isampmax_U;
								smaller_neighbor_isampmax_U[layer][module_on_layer][apv_U] = chan_j_isampmax_U;

								layer_U[layer] = layer;
								module_U[layer][module_on_layer] = module_on_layer;

								tr_ratios_U->Fill();
								neighbor_ratio_cnt_U++;
							}

				//Channel J is larger than Channel I:
				//Chan i & j DO NOT equal 0 AND they have matching time sample indices			
							if( (chan_j_ADCsum_U > chan_i_ADCsum_U) && (chan_j_ADCsum_U > ADCcut) && (chan_i_ADCsum_U != 0) && (chan_j_ADCsum_U != 0) && (chan_i_isampmax_U == chan_j_isampmax_U) ){
								neighbor_ADCsum_ratio_U = ( (1.0)*chan_j_ADCsum_U )/( (1.0)*chan_i_ADCsum_U );

								ratio_neighbor_channels_ADCsum_U[layer][module_on_layer][apv_U] = neighbor_ADCsum_ratio_U;
								larger_neighbor_ADCsum_U[layer][module_on_layer][apv_U] = chan_j_ADCsum_U;
								smaller_neighbor_ADCsum_U[layer][module_on_layer][apv_U] = chan_i_ADCsum_U;
								larger_neighbor_isampmax_U[layer][module_on_layer][apv_U] = chan_j_isampmax_U;
								smaller_neighbor_isampmax_U[layer][module_on_layer][apv_U] = chan_i_isampmax_U;

								layer_U[layer] = layer;
								module_U[layer][module_on_layer] = module_on_layer;

								tr_ratios_U->Fill();
								neighbor_ratio_cnt_U++;
							}			

						}
					}
				}

	//V-Strips
				if( build_V_strips ){
					// cout << " V ADC " << endl;
					double neighbor_ADCsum_ratio_V = -99999;
					double chan_i_ADCsum_V = 0.0;
					double chan_j_ADCsum_V = 0.0;
					int chan_i_isampmax_V = -1;
					int chan_j_isampmax_V = -1;

					for( int apv_V = first_apv; apv_V < last_apv_V; apv_V++ ){
						for( int chan_V = 0; chan_V < 127; chan_V++ ){
							// cout << "apv: " << apv_V << ", chan: " << chan_V << endl;
							//Initialize variables for U-Strip ADC values per channel
							neighbor_ADCsum_ratio_V = -99999;
							chan_i_ADCsum_V = 0.0;
							chan_j_ADCsum_V = 0.0;
							chan_i_isampmax_V = -1;
							chan_j_isampmax_V = -1;

							//STORE
							chan_i_ADCsum_V = (1.0)*APV_chans_ADCsum_V[apv_V][chan_V];
							chan_j_ADCsum_V = (1.0)*APV_chans_ADCsum_V[apv_V][chan_V + 1];
							chan_i_isampmax_V = APV_chans_isampmax_V[apv_V][chan_V];
							chan_j_isampmax_V = APV_chans_isampmax_V[apv_V][chan_V + 1];

				//Channel I is larger than Channel J:
				//Chan i & j DO NOT equal 0 AND they have matching time sample indices
							if( (chan_i_ADCsum_V > chan_j_ADCsum_V) && (chan_i_ADCsum_V > ADCcut) && (chan_i_ADCsum_V != 0) && (chan_j_ADCsum_V != 0) && (chan_i_isampmax_V == chan_j_isampmax_V) ){
								neighbor_ADCsum_ratio_V = ( (1.0)*chan_i_ADCsum_V )/( (1.0)*chan_j_ADCsum_V );

								ratio_neighbor_channels_ADCsum_V[layer][module_on_layer][apv_V] = neighbor_ADCsum_ratio_V;
								larger_neighbor_ADCsum_V[layer][module_on_layer][apv_V] = chan_i_ADCsum_V;
								smaller_neighbor_ADCsum_V[layer][module_on_layer][apv_V] = chan_j_ADCsum_V;
								larger_neighbor_isampmax_V[layer][module_on_layer][apv_V] = chan_i_isampmax_V;
								smaller_neighbor_isampmax_V[layer][module_on_layer][apv_V] = chan_j_isampmax_V;

								layer_V[layer] = layer;
								module_V[layer][module_on_layer] = module_on_layer;

								tr_ratios_V->Fill();
								neighbor_ratio_cnt_V++;
							}

				//Channel J is larger than Channel I:
				//Chan i & j DO NOT equal 0 AND they have matching time sample indices			
							if( (chan_j_ADCsum_V > chan_i_ADCsum_V) && (chan_j_ADCsum_V > ADCcut) && (chan_i_ADCsum_V != 0) && (chan_j_ADCsum_V != 0) && (chan_i_isampmax_V == chan_j_isampmax_V) ){
								neighbor_ADCsum_ratio_V = ( (1.0)*chan_j_ADCsum_V )/( (1.0)*chan_i_ADCsum_V );

								ratio_neighbor_channels_ADCsum_V[layer][module_on_layer][apv_V] = neighbor_ADCsum_ratio_V;
								larger_neighbor_ADCsum_V[layer][module_on_layer][apv_V] = chan_j_ADCsum_V;
								smaller_neighbor_ADCsum_V[layer][module_on_layer][apv_V] = chan_i_ADCsum_V;
								larger_neighbor_isampmax_V[layer][module_on_layer][apv_V] = chan_j_isampmax_V;
								smaller_neighbor_isampmax_V[layer][module_on_layer][apv_V] = chan_i_isampmax_V;

								layer_V[layer] = layer;
								module_V[layer][module_on_layer] = module_on_layer;

								tr_ratios_V->Fill();
								neighbor_ratio_cnt_V++;
							}			

						}
					}
				}		

		


//---------------------------------------------------------------------------------
//---------------------------END OF MAIN EVENTS LOOP--------------------------------
//---------------------------------------------------------------------------------

			}

			cout << endl << "#####   Clearing ADC values in preparation for next module...   #####" << endl << endl;			

			if( build_U_strips ){
				APV_chans_ADCsum_U.clear();
				APV_chans_isampmax_U.clear();
			}
			if( build_V_strips ){
				APV_chans_ADCsum_V.clear();
				APV_chans_isampmax_V.clear();
			}

			cout << "Writing data to TTrees." << endl;

			if( build_U_strips ){ tr_ratios_U->Write("", TObject::kOverwrite); }
			if( build_V_strips ){ tr_ratios_V->Write("", TObject::kOverwrite); }

	//Reset 3D arrays:

			for(int l = 0; l < nLayers; l++){
				layer_U[l] = 0;
				layer_V[l] = 0;
				for(int m = 0; m < 4; m++){
					module_U[l][m] = 0;
					module_V[l][m] = 0;
					for(int a = 0; a < 30; a++){
						ratio_neighbor_channels_ADCsum_U[l][m][a] = 0.0;
						larger_neighbor_ADCsum_U[l][m][a] = 0.0;
						smaller_neighbor_ADCsum_U[l][m][a] = 0.0;
						larger_neighbor_isampmax_U[l][m][a] = 0;
						smaller_neighbor_isampmax_U[l][m][a] = 0;

						ratio_neighbor_channels_ADCsum_V[l][m][a] = 0.0;
						larger_neighbor_ADCsum_V[l][m][a] = 0.0;
						smaller_neighbor_ADCsum_V[l][m][a] = 0.0;
						larger_neighbor_isampmax_V[l][m][a] = 0;
						smaller_neighbor_isampmax_V[l][m][a] = 0;
					}
				}
			}
			cout << "*******************************************************************************" << endl;
			cout << "-------------------------------------------------------------------------------" << endl;
			cout << "Finished Layer " << layer << " - Module " << module_on_layer << " (Global: " << module_num_global << ") loop. " << endl;
			cout << "-------------------------------------------------------------------------------" << endl;
			cout << "*******************************************************************************" << endl << endl;
			auto module_time_stop = high_resolution_clock::now();
			auto module_duration = duration_cast<minutes>(module_time_stop - total_time_start);
			cout << "-------------------------------------------------------------------------------" << endl;
			cout << "Time for analysis so far: " << module_duration.count() << " minutes. " << endl;
			cout << "-------------------------------------------------------------------------------" << endl << endl << endl;


		//END OF MODULES LOOP
		}

		if( build_U_strips ){
			APV_chans_ADCsum_U.clear();
			APV_chans_isampmax_U.clear();
		}
		if( build_V_strips ){
			APV_chans_ADCsum_V.clear();
			APV_chans_isampmax_V.clear();
		}
		cout << "*******************************************************************************" << endl;
		cout << "-------------------------------------------------------------------------------" << endl;
		cout << "Finished Layer " << layer << " loop. " << endl;
		cout << "-------------------------------------------------------------------------------" << endl << endl;
		auto layer_time_stop = high_resolution_clock::now();
		auto layer_duration = duration_cast<minutes>(layer_time_stop - total_time_start);
		cout << "-------------------------------------------------------------------------------" << endl;
		cout << "Time for analysis so far: " << layer_duration.count() << " minutes. " << endl;
		cout << "-------------------------------------------------------------------------------" << endl << endl;
	}
	cout << endl << "Finished with analysis on all layers. " << endl;
	cout << "Closing Rootfile." << endl;
	RatioRootfile->Close();
	cout << "Output Rootfile: " << endl;
	cout << rootfilename.Data() << endl;

	cout << "********************************************************" << endl;
	cout << "Run: " << runnum << ", ADCcut: " << ADCcut << ", First layer: " << first_layer << ", Last layer: " << last_layer << endl;
	cout << "Number of events: " << last_event << endl;
	cout << "Build U: " << build_U_strips << " - Build V: " << build_V_strips << endl;
	cout << "********************************************************" << endl;
}