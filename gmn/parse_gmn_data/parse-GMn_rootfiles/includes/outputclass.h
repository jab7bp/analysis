#include "/work/halla/sbs/adr/GMn_analysis/physics_analysis/ElasticEventsStudy/includes/eventclass.h"

#ifndef RESULTSCLASS_H
#define RESULTSCLASS_H

class Output
{
private:

	Event& m_event;
	const char* m_outputfilename;
	TFile* m_fout;
	TTree* m_resultstree;

	// Define the member variables to hold data for output TTree variables.
	double m_bbtrvz{0.};
	double m_bbtrth{0.};
	double m_bbtrx{0.};
	double m_bbtrp{0.};
	double m_bbtrpx{0.};
	double m_bbtrpy{0.};
	double m_bbshe{0.};
	double m_bbpse{0.};
	double m_bbshpse{0.};
	double m_bbeoverp{0.};
	double m_bbphitgt{0.};
	double m_sbshcalx{0.};
	double m_sbshcaly{0.};
	double m_sbshcale{0.};
	double m_timediffhcalbbcal{0.};

	// Output histograms
	TH1D* m_h1_bb_tr_th;
	TH2D* m_h2_bb_tr_th_vs_x;
	TH1D* m_h1_phi_tgt;
	TH2D* m_h2_hcalclusX_vs_phi;

public:

	Output(Event& event, const char* outputfilename) : m_event{event}, m_outputfilename{outputfilename}
	{
		m_fout = new TFile(Form("%s.root",outputfilename),"RECREATE");
		m_resultstree = new TTree("T","Down beding track analysis");

		m_resultstree->Branch("bb.tr.vz", &m_bbtrvz);
		m_resultstree->Branch("bb.tr.th", &m_bbtrth);
		m_resultstree->Branch("bb.tr.x", &m_bbtrx);
		m_resultstree->Branch("bb.tr.p", &m_bbtrp);
		m_resultstree->Branch("bb.tr.px", &m_bbtrpx);
		m_resultstree->Branch("bb.tr.py", &m_bbtrpy);
		m_resultstree->Branch("bb.sh.e", &m_bbshe);
		m_resultstree->Branch("bb.ps.e", &m_bbpse);
		m_resultstree->Branch("bb.shps.e", &m_bbshpse);
		m_resultstree->Branch("bb.eoverp", &m_bbeoverp);
		m_resultstree->Branch("bb.phitgt", &m_bbphitgt);
		m_resultstree->Branch("sbs.hcal.x", &m_sbshcalx);
		m_resultstree->Branch("sbs.hcal.y", &m_sbshcaly);
		m_resultstree->Branch("sbs.hcal.e", &m_sbshcale);
		m_resultstree->Branch("timediff.hcalbbcal", &m_timediffhcalbbcal);

		//Define the output histograms.
		m_h1_bb_tr_th = new TH1D("h1_bb_tr_th0", "BB track theta (dx/dz) distribution; dx/dz", 100, -0.2, 0.8);
		m_h2_bb_tr_th_vs_x = new TH2D("h2_bb_tr_th0_vs_x", "BB track dx/dz vs x distribution; x (m); dx/dz", 160, -0.8, 0.8, 120, -0.6, 0.6);
		m_h1_phi_tgt = new TH1D("h1_phi_tgt", "Azimuthal angle (#phi) distribution; #phi degrees", 700, -35, 35);
		m_h2_hcalclusX_vs_phi = new TH2D("h2_hcalclusX_vs_phi", "HCal cluster vertical pos vs #phi distribution; #phi degrees; sbs.hcal.x (m)", 700, -30, 35, 500, -3.0, 2.0);
	}


	void copyFromEvent()
	{
		m_bbtrvz = m_event.return_BBTrVz();
		m_bbtrth = m_event.return_BBTrth();
		m_bbtrx = m_event.return_BBTrx();
		m_bbtrp = m_event.return_BBTrP();
		m_bbtrpx = m_event.return_BBTrPx();
		m_bbtrpy = m_event.return_BBtrPy();
		m_bbshe = m_event.return_BBSHe();
		m_bbpse = m_event.return_BBPSe();
		m_bbshpse = m_event.return_BBSHPSe();
		m_bbeoverp = m_event.return_EoverP();
		m_bbphitgt = m_event.return_BBphitgt();
		m_sbshcalx = m_event.return_SBSHCalx();
		m_sbshcaly = m_event.return_SBSHCaly();
		m_sbshcale = m_event.return_SBSHCale();
		m_timediffhcalbbcal = m_event.return_timediffHCalBBCal();
	}

	
	void fillOutTree()
	{
		m_resultstree->Fill();
	}

	void fillHistos()
	{
		m_h1_bb_tr_th->Fill(m_bbtrth);
		m_h2_bb_tr_th_vs_x->Fill(m_bbtrx, m_bbtrth);
		m_h1_phi_tgt->Fill(m_bbphitgt);
		m_h2_hcalclusX_vs_phi->Fill(m_bbphitgt, m_sbshcalx);
	}

	void closeOutFile()
	{
		m_resultstree->Write(0, TObject::kWriteDelete, 0);
		m_h1_bb_tr_th->Write();
		m_h2_bb_tr_th_vs_x->Write();
		m_h1_phi_tgt->Write();
		m_h2_hcalclusX_vs_phi->Write();
	}
};

#endif