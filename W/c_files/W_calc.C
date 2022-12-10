// This function makes the plots
TH1F *make_plot(TString run, TString current){
  
  
  //TString Rootfile = "/chafs2/work1/sbs/Rootfiles/gmn_replayed_" + run + "_stream0_seg0_0.root";
  TString Rootfile = "/Users/john/UVa/SBS/inv_mass/gmn_replayed_" + run + "_stream0_seg0_0.root";

  TChain *t = new TChain("T");

  t->Add(Rootfile);
  
  TH1F *h = new TH1F("hW_I" + current,"BB Elastic Peak;W (GeV);Counts",200,0.9,3);

  double M_p = 0.938;  //proton mass
  double E = 1.92;     //Electron Beam energy

  TCut cut1 = "bb.gem.track.nhits[0] > 3"; //number of hits on tracks > 3
  //TCut cut2 = "(bb.sh.e + bb.ps.e)/bb.tr.p > 0.78 && (bb.ps.e + bb.ps.e)/bb.tr.p < 1.15";  //For this cut we plot the variable and cut on the peak (see below)
  TCut cut2 = "(bb.sh.e + bb.ps.e)/bb.tr.p > 0.6 && (bb.ps.e + bb.ps.e)/bb.tr.p < 0.8";

  TCut cut_tot = cut1 + cut2;
  
  //This calcualtion is W = sqrt(M^2 + 2M(E - E') - 2E*E'(1 - cos(theta))

  t->Draw(Form("sqrt(%f^2 + 2*%f*(%f - bb.tr.p[0]) - 2*%f*bb.tr.p[0]*(1 - cos(acos(bb.tr.pz[0]/bb.tr.p[0])))) >> hW_I" + current,M_p, M_p, E, E), cut_tot);

  
  //t->Draw("(bb.sh.e+bb.ps.e)/bb.tr.p", "bb.ps.e > 0.1");

  TLorentzVector *TLV = new TLorentzVector(0.0, 0.0, 0.0, 0.938);

  return h;

}




////// This is the main function ///////
void W_calc(){

  gStyle->SetOptStat(0);

  // Set lines where we are cutting on the elastic peak
  // double low_x = 1.55;
  // double high_x = 2.;

  double low_x = 0.96;
  double high_x = 1.5;

  // List of runs and currents
  const int np = 6;
  //TString runs[np] = {"11420","11421","11422","11423","11425","11432"};
  TString runs[np] = {"11207", "11214", "11218", "11226", "11228"};
  //TString currents[np] = {"1","2","3","4","5","6"};
  TString currents[np] = {"1", "2", "3", "4", "5"};

  //TString runs[np] = {"11207","11214","11218","11224","11228"};
  //TString currents[np] = {"1","2","3","4","5"};

  TCanvas *c = new TCanvas("c","",1000,800);

  TH1F *hW[np];

  //Fill our plots
  for(int i = 0; i < np; i++){
    hW[i] = make_plot(runs[i], currents[i]);
  }


  TLegend *leg = new TLegend(0.1, 0.6, 0.38, 0.9);
  leg->SetHeader("Beam Currents with # of Events in Cuts", "C");

  //Loop through plots and add them to the canvas
  double max = 0;
  int color = 0;
  for(int i = 0; i < np; i++){
    
    color++;
    if(color == 3 || color == 5) color++; //Skip annoying colors
    hW[i]->SetLineColor(color);
   
    hW[i]->Draw("same");   //Draw plots on same canvas
    
    //Get number of events in the elastic peak
    double nevents = hW[i]->Integral(hW[i]->GetXaxis()->FindBin(low_x),hW[i]->GetXaxis()->FindBin(high_x));

    //Add histogram to the levend
    leg->AddEntry("hW_I" + currents[i], Form(currents[i] + " #muA (%i events) Run " + runs[i],(int)nevents), "l");

    //Record the tallest histogram for scaling at the end
    if(hW[i]->GetMaximum() > max) max = hW[i]->GetMaximum();
  }


  leg->Draw("same");

  hW[np - 1]->GetYaxis()->SetRangeUser(0,max*1.1); //Scale histogram

  //Draw lines to mark elastic peak cut
  TLine *tl = new TLine(low_x,0,low_x,max*1.1);
  TLine *th = new TLine(high_x,0,high_x,max*1.1);
  tl->SetLineColor(kRed);
  th->SetLineColor(kRed);
  tl->Draw("same");
  th->Draw("same");

  //Write other cuts on canvas
  TPaveText *pt1 = new TPaveText(0.7,0.78,0.9,0.89,"nbNDC");
  pt1->AddText("LH2 Target");
  pt1->AddText("hits on track > 3");
  pt1->AddText("0.78 < E/p < 1.15");
  pt1->SetFillColor(0);

  pt1->Draw("same");

}
