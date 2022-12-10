// This function makes the plots
Double_t calcW(TLorentzVector target, TLorentzVector beam, TLorentzVector scatter, double theta){
  Double_t W;
  TLorentzVector Q, p, q;
  scatter.SetTheta(theta);
  q = beam - scatter;
  p = target;
  Q = p + q;
  W = Q.Mag();
  
  return W;

}




////// This is the main function ///////
void calc_inv_mass_multi(){

// gSystem.Load("libPhysics.so");

// List of runs and currents
const int np = 5;
//TString runs[np] = {"11420","11421","11422","11423","11425","11432"};
TString runs[np] = {"11207", "11214", "11218", "11226", "11228"};
//TString currents[np] = {"1","2","3","4","5","6"};
TString currents[np] = {"1", "2", "3", "4", "5"};
TString rootfiles[np] = {};
TH1F* histW[np] = {};
int color = 1;
THStack *hWStack = new THStack("hWStack", "Stack of Histograms");


for(int i = 0; i < np; i++){
  rootfiles[i] = "/Users/john/UVa/SBS/inv_mass/gmn_replayed_" + runs[i] + "_stream0_seg0_0.root";
  histW[i] = new TH1F("histW_" + runs[i], "inv_mass", 200, 0.4, 1.8);
  }

for(int m = 0; m < (np); m++){
  cout << "currently on run: " << runs[m] << endl;

  TFile *tf = new TFile(rootfiles[m], "READ");
  TTree *t = (TTree*)tf->Get("T");


  // TH1F *hW = new TH1F("hW", "Invariant Mass", 200, 0.4, 1.8);
  // TH1F *hp = new TH1F("hp", "Invariant Mass", 200, 0, 3); 

  Double_t bb_sh_e, bb_ps_e, nhits, bb_tr_pz, bb_tr_p, theta, W, cval;

  double E = 1.92;
  double Mp = 0.938;

  int p_index[50000] = {};


  int pz_index[50000] = {};
  int cnt_p = 0;
  int cnt_pz = 0;
  int p_next = 0;




  TLorentzVector target, beam, scatter;
  target.SetPxPyPzE(0.0, 0.0, 0.0, Mp);
  beam.SetPxPyPzE(0.0, 0.0, E, E);

  // t->SetBranchAddress("bb.tr.pz", &bb_tr_pz);
  // t->SetBranchAddress("bb.tr.p", &bb_tr_p);
  // t->SetBranchAddress("bb.sh.e", &bb_sh_e);
  // t->SetBranchAddress("bb.ps.e", &bb_ps_e);
  // t->SetBranchAddress("bb.gem.track.nhits", &nhits);

    
  Int_t entries = 50000;
  Double_t p_arr [50000] = {};
  Double_t new_p[50000] = {};
  Double_t pz_arr [50000] = {};
  Double_t new_pz[50000] = {};
  Double_t E_sh_arr[50000] = {};
  Double_t new_E_sh[50000] = {};
  Double_t E_ps_arr[50000] = {};
  Double_t new_E_ps[50000] = {};
  Int_t nhits_arr[50000] = {};
  Int_t new_nhits[50000] = {};

  cout << "Gathering data and building arrays...." << endl;
  for (Int_t b = 0; b<entries; b++) {
    t->GetEntry(b);
    //cout << "i: " << i << "  getval 0 : " << t->GetLeaf("bb.tr.p")->GetValueLongDouble(0) << " pz: " << t->GetLeaf("bb.tr.pz")->GetValueLongDouble(0) <<endl;
    p_arr[b] = t->GetLeaf("bb.tr.p")->GetValueLongDouble(0);
    pz_arr[b] = t->GetLeaf("bb.tr.pz")->GetValueLongDouble(0);
    E_sh_arr[b] = t->GetLeaf("bb.sh.e")->GetValueLongDouble(0);
    E_ps_arr[b] = t->GetLeaf("bb.ps.e")->GetValueLongDouble(0);
    nhits_arr[b] = t->GetLeaf("bb.gem.track.nhits")->GetValueLongDouble(0);

    if(b == entries/4){cout << "1/4th of the data sifted through...." << endl;}
    if(b == entries/2){cout << "Halfway there...." << endl;}
    if(b == entries*(.75)){cout << "75 percent done." << endl;}
    if(b == entries - 10){cout << "Ten more to go....." << endl;}
    //cout << "i: " << i << "   p: " << p_arr[i] << "   pz: " << pz_arr[i] << endl;

    //cout << "p[i] " << p_arr[i] << endl;
  }

  for (Int_t i = 1; i<(entries+1); i++) {
    if(p_arr[i] != p_arr[i-1]){
      new_p[cnt_p] = p_arr[i-1];
      new_E_ps[cnt_p] = E_ps_arr[i-1];
      new_E_sh[cnt_p] = E_sh_arr[i-1];
      new_nhits[cnt_p] = nhits_arr[i-1];
      // cout << "i: " << i << "   cnt: " << cnt_p << "  p new " << new_p[cnt_p] << "   p arr - 1 "  << p_arr[i-1] << endl;
      p_index[cnt_p] = i;
      cnt_p ++;
    }
    else{
      p_arr[i-1] = 0;
      }
    
    if(pz_arr[i] != pz_arr[i-1]){
      new_pz[cnt_pz] = pz_arr[i-1];
      pz_index[cnt_p] = i;
      cnt_pz ++;
    }
    else{
      pz_arr[i-1] = 0;
      }
    }

  cout << "************" << endl << "Calculating W" << endl << "************" << endl;
  for(int i = 0; i < 50000; i++){
    if((new_p[i] > 0) && (new_pz[i] > 0)) {
      
      
      //theta = acos(bb_tr_pz/bb_tr_p);
      theta = acos(new_pz[i]/new_p[i]);
    
      //scatter.SetPxPyPzE(0.0, 0.0, bb_tr_p, bb_tr_p);
      scatter.SetPxPyPzE(0.0, 0.0, new_p[i], new_p[i]);

      //if ((nhits > 3) && (0.6 < ((bb_sh_e + bb_ps_e)/bb_tr_p)) && (0.8 > ((bb_sh_e + bb_ps_e)/bb_tr_p))){
      if ((new_nhits[i] > 3) && (0.6 < ((new_E_sh[i] + new_E_ps[i])/new_p[i])) && (0.8 > ((new_E_sh[i] + new_E_ps[i])/new_p[i]))) {
        W = calcW(target, beam, scatter, theta);
       //cout << "W: " << W << endl;
        histW[m]->Fill(W);
      //cout << "i: " << i << "  p : " << bb_tr_p << "  pz: " << bb_tr_pz << endl;
      }
    } 
  }

  hWStack->Add(histW[m]);
  hWStack->Draw("nostack");
  //histW[m]->Draw("SAME");
  if(m == 4){
    color = 6;
  }
  histW[m]->SetLineColor(color);

  color++;

}

// for(int i = 0; i<50000; i++){
//   cout << "i: " << i <<  "    E sh: " << E_sh_arr[i] << "    E_ps: " << E_ps_arr[i] << endl;



}


