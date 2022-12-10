#include <iostream>

int macro()
{
	gROOT->Reset();
	TH1F *hist = new TH1F();
	TFile *tf = new TFile("/Users/john/UVa/SBS/inv_mass/gmn_replayed_11207_stream0_seg0_0.root", "READ");

	TTree *tree = (TTree*)tf->Get("T");

	Double_t bb_tr_p;
	TBranch *p;

	tree->SetBranchAddress("bb.tr.p", &bb_tr_p);

	p = tree->GetBranch("bb.tr.p");

	int nEntries = tree->GetEntries();

	for (int iEnt = 0; iEnt < nEntries; iEnt++) {
		p->GetEvent(iEnt);
		
		for(int i = 0; i < 500; i++) {
			tree->GetEntry(i);

		}
		
	}

	

	cout << "entries: " << nEntries <<  endl;

	int a = 5;

	return 0;
}