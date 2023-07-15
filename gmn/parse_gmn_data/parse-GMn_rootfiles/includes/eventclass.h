#ifndef EVENTCLASS_H
#define EVENTCLASS_H

class Event
{
private:
	
	// Read-in from the main script.
	TChain* m_C; // TChain that holds the "T" TTree(s) of input replayed ROOT file(s).
	const double m_cut_HCalE;
	const double m_cut_CoinCutLow;
	const double m_cut_CoinCutHigh;
	const double m_cut_BBTrPCut;

	// Variables needed to copy information from "T" tree.
	const static int m_MAXNTRACKS {100};
	double m_bbtrvz[m_MAXNTRACKS];
	double m_bbtrth[m_MAXNTRACKS];
	double m_bbtrx[m_MAXNTRACKS];
	double m_bbtrp[m_MAXNTRACKS];
	double m_bbtrpx[m_MAXNTRACKS];
	double m_bbtrpy[m_MAXNTRACKS];
	double m_sbshcale{0.};
	double m_sbshcalx{0.};
	double m_sbshcaly{0.};
	double m_bbshe{0.};
	double m_bbpse{0.};
	const static int m_MAXNTDC{1000};
	double m_tdc_time[m_MAXNTDC];
	double m_tdc_elemID[m_MAXNTDC];
	int m_ndata_tdc{0};
	double m_timediff_hcalbbcal{0.};
	double m_hcal_clusblk_ADC_time[15]; // Maximum number of blocks in a cluster is 15 as per S.Seeds.	

public:

	Event(TChain* C, double cut_HCalE, double cut_CoinCutLow, double cut_CoinCutHigh, double cut_BBTrPCut) 
	: m_C{C}, m_cut_HCalE{cut_HCalE}, m_cut_CoinCutLow{cut_CoinCutLow}, m_cut_CoinCutHigh{cut_CoinCutHigh}, m_cut_BBTrPCut{cut_BBTrPCut}
	{
		m_C->SetBranchStatus("*", 0);
		m_C->SetBranchStatus("bb.tr.vz",1);
		m_C->SetBranchStatus("bb.tr.th", 1);
		m_C->SetBranchStatus("bb.tr.x", 1);
		m_C->SetBranchStatus("bb.tr.p",1);
		m_C->SetBranchStatus("bb.tr.px", 1);
		m_C->SetBranchStatus("bb.tr.py", 1);
		m_C->SetBranchStatus("sbs.hcal.x", 1);
		m_C->SetBranchStatus("sbs.hcal.y", 1);
		m_C->SetBranchStatus("sbs.hcal.e",1);
		m_C->SetBranchStatus("bb.sh.e",1);
		m_C->SetBranchStatus("bb.ps.e", 1);
		m_C->SetBranchStatus("bb.tdctrig.tdc",1);
		m_C->SetBranchStatus("bb.tdctrig.tdcelemID",1);
		m_C->SetBranchStatus("Ndata.bb.tdctrig.tdcelemID",1);
		m_C->SetBranchStatus("sbs.hcal.clus_blk.atime",1);

		m_C->SetBranchAddress("bb.tr.vz", m_bbtrvz);
		m_C->SetBranchAddress("bb.tr.th", m_bbtrth);
		m_C->SetBranchAddress("bb.tr.x", m_bbtrx);
		m_C->SetBranchAddress("bb.tr.p", m_bbtrp);
		m_C->SetBranchAddress("bb.tr.px", m_bbtrpx); 
		m_C->SetBranchAddress("bb.tr.py", m_bbtrpy); 
		m_C->SetBranchAddress("sbs.hcal.x", &m_sbshcalx);
		m_C->SetBranchAddress("sbs.hcal.y", &m_sbshcaly);
		m_C->SetBranchAddress("sbs.hcal.e", &m_sbshcale);
		m_C->SetBranchAddress("bb.sh.e", &m_bbshe);
		m_C->SetBranchAddress("bb.ps.e", &m_bbpse);
		m_C->SetBranchAddress("bb.tdctrig.tdc", m_tdc_time);
		m_C->SetBranchAddress("bb.tdctrig.tdcelemID", m_tdc_elemID);
		m_C->SetBranchAddress("Ndata.bb.tdctrig.tdcelemID", &m_ndata_tdc);
		m_C->SetBranchAddress("sbs.hcal.clus_blk.atime", m_hcal_clusblk_ADC_time);
	}

	int getEntry(int n) //Copies the enries of the "T" to the above defined member variables.
	{
		return m_C->GetEntry(n);
	}

private:	

	void get_TDC_times(int ndata_tdc, double tdc_elemID[1000], double tdc_time[1000], double& bbcal_time, double& hcal_time, double& coin_time, double&rf_time)
	{
		for(int ihit=0; ihit<ndata_tdc; ihit++)
		{
	      if(tdc_elemID[ihit]==5) bbcal_time=tdc_time[ihit];
	      if(tdc_elemID[ihit]==0) hcal_time=tdc_time[ihit];
	      if(tdc_elemID[ihit]==1) coin_time=tdc_time[ihit];
	      if(tdc_elemID[ihit]==4) rf_time=tdc_time[ihit];
	    }
	}

public:

	bool passCuts()
	{	
		if ( m_sbshcale < m_cut_HCalE ) return false;

		// Calculate BBCal and HCal coincidence times
		double bbcal_time{0.};
		double hcal_time{0.};
		double coin_time{0.};
		double rf_time{0.};
		get_TDC_times(m_ndata_tdc, m_tdc_elemID, m_tdc_time, bbcal_time, hcal_time, coin_time, rf_time);
		m_timediff_hcalbbcal = hcal_time - bbcal_time;
		if ( m_timediff_hcalbbcal < m_cut_CoinCutLow || m_timediff_hcalbbcal > m_cut_CoinCutHigh ) return false;

		return true;
	}

	double return_BBTrVz() 
	{
		return m_bbtrvz[0];
	}
	
	double return_BBTrth()
	{
		return m_bbtrth[0];
	}

	double return_BBTrx() 
	{
		return m_bbtrx[0];
	}
	
	double return_BBTrP()
	{
		return m_bbtrp[0];
	}
	
	double return_BBTrPx()
	{
		return m_bbtrpx[0];
	}

	double return_BBtrPy()
	{
		return m_bbtrpy[0];
	}
	
	double return_BBSHe()
	{
		return m_bbshe;
	}
	double return_BBPSe()
	{
		return m_bbpse;
	}

	double return_SBSHCalx()
	{
		return m_sbshcalx;
	}

	double return_SBSHCaly()
	{
		return m_sbshcaly;
	}

	double return_SBSHCale()
	{
		return m_sbshcale;
	}
	double return_timediffHCalBBCal()
	{
		return m_timediff_hcalbbcal;
	}

private: 	

	// Variables to hold the secondary calculation results done using the above variables copied from the "T" tree.
	double m_bbshpse{0.}; // bb.sh.e + bb.ps.e = Total energy deposited in the sh and ps.
	double m_eoverp{0.};
	double m_bbphitgt{0.}; // Azimuthal angle of pi+ detected by the BigBite spectrometer.


public: 

	double return_BBSHPSe()
	{
		m_bbshpse = m_bbshe + m_bbpse;
		return m_bbshpse;
	}

	double return_EoverP() //Make sure to call the above "return_BBSHPSe()" function before calling this function.
	{
		m_eoverp = m_bbshpse / m_bbtrp[0];
		return m_eoverp;
	}

	double return_BBphitgt()
	{
		m_bbphitgt = atan2(m_bbtrpy[0], m_bbtrpx[0])*TMath::RadToDeg();
		return m_bbphitgt;
	}	 	

};

#endif