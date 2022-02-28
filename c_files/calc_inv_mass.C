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
void calc_inv_mass(){

// gSystem.Load("libPhysics.so");

TString Rootfile = "/Users/john/UVa/SBS/inv_mass/gmn_replayed_11207_stream0_seg0_0.root";

TFile *tf = new TFile(Rootfile);
TTree *t = (TTree*)tf->Get("T");


TH1F *hW = new TH1F("hW", "Invariant Mass", 200, 0.4, 1.8);
TH1F *hp = new TH1F("hp", "Invariant Mass", 200, 0, 3); 

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

t->SetBranchAddress("bb.tr.pz", &bb_tr_pz);
t->SetBranchAddress("bb.tr.p", &bb_tr_p);
t->SetBranchAddress("bb.sh.e", &bb_sh_e);
t->SetBranchAddress("bb.ps.e", &bb_ps_e);
t->SetBranchAddress("bb.gem.track.nhits", &nhits);

  
Int_t entries = t->GetEntries();
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

for (Int_t i = 0; i<entries; i++) {
  t->GetEntry(i);
  //cout << "i: " << i << "  getval 0 : " << t->GetLeaf("bb.tr.p")->GetValueLongDouble(0) << " pz: " << t->GetLeaf("bb.tr.pz")->GetValueLongDouble(0) <<endl;
  p_arr[i] = t->GetLeaf("bb.tr.p")->GetValueLongDouble(0);
  pz_arr[i] = t->GetLeaf("bb.tr.pz")->GetValueLongDouble(0);
  E_sh_arr[i] = t->GetLeaf("bb.sh.e")->GetValueLongDouble(0);
  E_ps_arr[i] = t->GetLeaf("bb.ps.e")->GetValueLongDouble(0);
  nhits_arr[i] = t->GetLeaf("bb.gem.track.nhits")->GetValueLongDouble(0);

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




for(int i = 0; i < 50000; i++){
  if((new_p[i] > 0) && (new_pz[i] > 0)) {
    //cout << "i: " << i << "    new_p: " << new_p[i] << "  new_pz: " << new_pz[i] << " new_E_sh: " << new_E_sh[i] << "  new E_ps: " << new_E_ps[i] << endl;
  
  

    //theta = acos(bb_tr_pz/bb_tr_p);
    theta = acos(new_pz[i]/new_p[i]);
  
    //scatter.SetPxPyPzE(0.0, 0.0, bb_tr_p, bb_tr_p);
    scatter.SetPxPyPzE(0.0, 0.0, new_p[i], new_p[i]);

    //if ((nhits > 3) && (0.6 < ((bb_sh_e + bb_ps_e)/bb_tr_p)) && (0.8 > ((bb_sh_e + bb_ps_e)/bb_tr_p))){
    if ((new_nhits[i] > 3) && (0.6 < ((new_E_sh[i] + new_E_ps[i])/new_p[i])) && (0.8 > ((new_E_sh[i] + new_E_ps[i])/new_p[i]))) {
      W = calcW(target, beam, scatter, theta);
     //cout << "W: " << W << endl;
      hW->Fill(W);
    //cout << "i: " << i << "  p : " << bb_tr_p << "  pz: " << bb_tr_pz << endl;
    }
  } 
}

hW->Draw();
// for(int i = 0; i<50000; i++){
//   cout << "i: " << i <<  "    E sh: " << E_sh_arr[i] << "    E_ps: " << E_ps_arr[i] << endl;

}


