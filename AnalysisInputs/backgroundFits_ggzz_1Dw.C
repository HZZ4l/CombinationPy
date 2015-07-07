/* 
 * Fit ggZZ background shapes and write parameters in a card fragment.
 * usage: 
 * -set all input variables in Config.h
 * -run with:
 * root -q -b backgroundFits_ggzz_1Dw.C
 * This runs on the 3 final states for 7 and 8 TeV and writes the output in a file (see stdout).
 *
 */

/*
 #ifndef __CINT__
 #include "RooGlobalFunc.h"
 #endif
 #include "RooRealVar.h"
 #include "RooDataSet.h"
 #include "RooGaussian.h"
 #include "RooConstVar.h"
 #include "RooChebychev.h"
 #include "RooAddPdf.h"
 #include "RooWorkspace.h"
 #include "RooPlot.h"
 #include "TCanvas.h"
 #include "TAxis.h"
 #include "TFile.h"
 #include "TH1.h"
*/

#include <iostream>
#include <iomanip>

using namespace RooFit;
using namespace std;



//----------> SET INPUT VARIABLES in Config.h
#include "Config.h"
//<----------

int category(
	     int nExtraLeptons,
	     float ZZPt,
	     float ZZMass,
	     int nJets, 
	     int nBTaggedJets,
	     //float* jetpt,
	     //float* jeteta,
	     //float* jetphi,
	     //float* jetmass,
	     float Fisher
	     )
{
  
  int category = -1;
  // 0 = Untagged  
  // 1 = 1-jet tagged  
  // 2 = VBF tagged  
  // 3 = VH-leptonic tagged  
  // 4 = VH-hadronic tagged  
  // 5 = ttH tagged  

  if( nExtraLeptons==0 && nJets>=2 && nBTaggedJets<=1 && Fisher>0.5 ){
    
    category = 2; // VBF tagged
    
  }else if( ( nExtraLeptons==0 && nJets>=2 && ZZPt>ZZMass && false /*flagDijetVH(nJets,jetpt,jeteta,jetphi,jetmass)*/ )
            || ( nExtraLeptons==0 && nJets==2 && nBTaggedJets==2 ) ){

    category = 4; // VH-hadronic tagged

  }else if( nExtraLeptons>=1 && nJets<=2 && nBTaggedJets==0 ){

    category = 3; // VH-leptonic tagged

  }else if( nExtraLeptons>=1 || (nJets>=3 && nBTaggedJets>=1) ){

    category = 5; // ttH tagged

  }else if(nJets>=1){

    category = 1; // 1-jet tagged

  }else{

    category = 0; // Untagged

  }

  return category;

}

void backgroundFits_ggzz_1Dw(int channel, int sqrts, int VBFtag);

// Run all final states and sqrts in one go
void backgroundFits_ggzz_1Dw() {
  // gSystem->Exec("mkdir -p bkgFigs7TeV");
  gSystem->Exec("mkdir -p bkgFigs13TeV");

  for(int icat=0;icat<4;icat++){
    backgroundFits_ggzz_1Dw(1,13,icat);
    backgroundFits_ggzz_1Dw(2,13,icat);
    backgroundFits_ggzz_1Dw(3,13,icat);
  }
}

// The actual job
void backgroundFits_ggzz_1Dw(int channel, int sqrts, int VBFtag)
{
  TString schannel;
  if      (channel == 1) schannel = "4mu";
  else if (channel == 2) schannel = "4e";
  else if (channel == 3) schannel = "2e2mu";
  else cout << "Not a valid channel: " << schannel << endl;

  TString ssqrts = (long) sqrts + TString("TeV");

  cout << "schannel = " << schannel << "  sqrts = " << sqrts << " VBFtag = "<< VBFtag << endl;

  TString outfile;
  outfile = "CardFragments/ggzzBackgroundFit_" + ssqrts + "_" + schannel + "_" + Form("%d",int(VBFtag)) + ".txt";
  ofstream of(outfile,ios_base::out);

  gSystem->AddIncludePath("-I$ROOFITSYS/include");
  gROOT->ProcessLine(".L ../CreateDatacards/include/tdrstyle.cc");
  setTDRStyle(false);
  gStyle->SetPadLeftMargin(0.16);
	
  TString filepath;filepath.Form("AAAOK/ZZ%s/ZZ4lAnalysis.root",schannel.Data());
  TFile *f = TFile::Open(filepath);
  TTree *tree = f->Get("ZZTree/candTree");

  RooRealVar* MC_weight = new RooRealVar("MC_weight","MC_weight",0.,2.) ; 
  RooRealVar* ZZMass = new RooRealVar("ZZMass","ZZMass",100,100.,1000.);
  RooRealVar* NJets30 = new RooRealVar("NJets30","NJets30",0.,5.);
  RooArgSet ntupleVarSet(*ZZMass,*NJets30,*MC_weight);
  RooDataSet *set = new RooDataSet("set","set",ntupleVarSet,WeightVar("MC_weight"));
  //RooArgSet ntupleVarSet(*ZZMass,*NJets30);  
  //RooDataSet *set = new RooDataSet("set","set",ntupleVarSet);

  Float_t myMC,myMass;
  Int_t myNJets;
  int nentries = tree->GetEntries();

  Float_t myPt,myJetPt,myJetEta,myJetPhi,myJetMass,myFisher;
  Int_t myExtralep,myBJets;
  tree->SetBranchAddress("ZZMass",&myMass);
  tree->SetBranchAddress("genHEPMCweight",&myMC);
  tree->SetBranchAddress("nCleanedJetsPt30",&myNJets);
  tree->SetBranchAddress("ZZPt",&myPt);
  tree->SetBranchAddress("nExtraLep",&myExtralep);
  tree->SetBranchAddress("nCleanedJetsPt30BTagged",&myBJets);
  tree->SetBranchAddress("DiJetDEta",&myFisher);

  for(int i =0;i<nentries;i++) {
    tree->GetEntry(i);
    if(myMass<100.)continue;
    int cat = category(myExtralep,myPt, myMass,myNJets, myBJets,/* jetpt, jeteta, jetphi, jetmass,*/myFisher);
    if(VBFtag != cat )continue;

    ntupleVarSet.setRealValue("ZZMass",myMass);
    ntupleVarSet.setRealValue("MC_weight",myMC);
    ntupleVarSet.setRealValue("NJets30",(double)cat);

    set->add(ntupleVarSet, myMC);
  }

  //RooRealVar* ZZLD = new RooRealVar("ZZLD","ZZLD",0.,1.);
  //char cut[10];
  //sprintf(cut,"ZZLD>0.5");
  //RooDataSet* set = new RooDataSet("set","set",tree,RooArgSet(*ZZMass,*MC_weight,*ZZLD),cut,"MC_weight");

  double totalweight = 0.;
  for (int i=0 ; i<set->numEntries() ; i++) { 
    set->get(i) ; 
    totalweight += set->weight();
    //cout << CMS_zz4l_mass->getVal() << " = " << set->weight() << endl ; 
  } 
  cout << "nEntries: " << set->numEntries() << ", totalweight: " << totalweight << endl;
		
  gSystem->Load("libHiggsAnalysisCombinedLimit.so");
	
  //// ---------------------------------------
  //Background
  RooRealVar CMS_qqzzbkg_a0("CMS_qqzzbkg_a0","CMS_qqzzbkg_a0",115.3,0.,200.);
  RooRealVar CMS_qqzzbkg_a1("CMS_qqzzbkg_a1","CMS_qqzzbkg_a1",21.96,0.,200.);
  RooRealVar CMS_qqzzbkg_a2("CMS_qqzzbkg_a2","CMS_qqzzbkg_a2",122.8,0.,200.);
  RooRealVar CMS_qqzzbkg_a3("CMS_qqzzbkg_a3","CMS_qqzzbkg_a3",0.03479,0.,1.);
  RooRealVar CMS_qqzzbkg_a4("CMS_qqzzbkg_a4","CMS_qqzzbkg_a4",185.5,0.,200.);
  RooRealVar CMS_qqzzbkg_a5("CMS_qqzzbkg_a5","CMS_qqzzbkg_a5",12.67,0.,200.);
  RooRealVar CMS_qqzzbkg_a6("CMS_qqzzbkg_a6","CMS_qqzzbkg_a6",34.81,0.,100.);
  RooRealVar CMS_qqzzbkg_a7("CMS_qqzzbkg_a7","CMS_qqzzbkg_a7",0.1393,0.,1.);
  RooRealVar CMS_qqzzbkg_a8("CMS_qqzzbkg_a8","CMS_qqzzbkg_a8",66.,0.,200.);
  RooRealVar CMS_qqzzbkg_a9("CMS_qqzzbkg_a9","CMS_qqzzbkg_a9",0.07191,0.,1.);
	
  RooggZZPdf_v2* bkg_ggzz = new RooggZZPdf_v2("bkg_ggzz","bkg_ggzz",*ZZMass,
					      CMS_qqzzbkg_a0,CMS_qqzzbkg_a1,CMS_qqzzbkg_a2,CMS_qqzzbkg_a3,CMS_qqzzbkg_a4,
					      CMS_qqzzbkg_a5,CMS_qqzzbkg_a6,CMS_qqzzbkg_a7,CMS_qqzzbkg_a8,CMS_qqzzbkg_a9);
	
  //// ---------------------------------------
	
  RooFitResult *r1 = bkg_ggzz->fitTo( *set, Save(kTRUE), SumW2Error(kTRUE) );//, Save(kTRUE), SumW2Error(kTRUE)) ;

  cout << endl;
  cout << "------- Parameters for " << schannel << " sqrts=" << sqrts << endl;
  cout << "  a0_bkgd = " << CMS_qqzzbkg_a0.getVal() << endl;
  cout << "  a1_bkgd = " << CMS_qqzzbkg_a1.getVal() << endl;
  cout << "  a2_bkgd = " << CMS_qqzzbkg_a2.getVal() << endl;
  cout << "  a3_bkgd = " << CMS_qqzzbkg_a3.getVal() << endl;
  cout << "  a4_bkgd = " << CMS_qqzzbkg_a4.getVal() << endl;
  cout << "  a5_bkgd = " << CMS_qqzzbkg_a5.getVal() << endl;
  cout << "  a6_bkgd = " << CMS_qqzzbkg_a6.getVal() << endl;
  cout << "  a7_bkgd = " << CMS_qqzzbkg_a7.getVal() << endl;
  cout << "  a8_bkgd = " << CMS_qqzzbkg_a8.getVal() << endl;
  cout << "  a9_bkgd = " << CMS_qqzzbkg_a9.getVal() << endl;
  cout << "---------------------------" << endl << endl;  

  of << "ggZZshape a0_bkgd  " << CMS_qqzzbkg_a0.getVal() << endl;
  of << "ggZZshape a1_bkgd  " << CMS_qqzzbkg_a1.getVal() << endl;
  of << "ggZZshape a2_bkgd  " << CMS_qqzzbkg_a2.getVal() << endl;
  of << "ggZZshape a3_bkgd  " << CMS_qqzzbkg_a3.getVal() << endl;
  of << "ggZZshape a4_bkgd  " << CMS_qqzzbkg_a4.getVal() << endl;
  of << "ggZZshape a5_bkgd  " << CMS_qqzzbkg_a5.getVal() << endl;
  of << "ggZZshape a6_bkgd  " << CMS_qqzzbkg_a6.getVal() << endl;
  of << "ggZZshape a7_bkgd  " << CMS_qqzzbkg_a7.getVal() << endl;
  of << "ggZZshape a8_bkgd  " << CMS_qqzzbkg_a8.getVal() << endl;
  of << "ggZZshape a9_bkgd  " << CMS_qqzzbkg_a9.getVal() << endl;
  of << endl;
  of.close();

  cout << endl << "Output written to: " << outfile << endl;

  int iLineColor = 1;
  string lab = "blah";
  if (channel == 1) { iLineColor = 2; lab = "4#mu"; }
  if (channel == 3) { iLineColor = 4; lab = "2e2#mu"; }
  if (channel == 2) { iLineColor = 6; lab = "4e"; }
  char lname[192];
  sprintf(lname,"gg #rightarrow ZZ #rightarrow %s", lab.c_str() );
  char lname2[192];
  sprintf(lname2,"Shape Model, %s", lab.c_str() );
  // dummy!                                                                                                                                               
  TF1* dummyF = new TF1("dummyF","1",0.,1.);
  TH1F* dummyH = new TH1F("dummyH","",1, 0.,1.);
  dummyF->SetLineColor( iLineColor );
  dummyF->SetLineWidth( 2 );
  
  TLegend * box2 = new TLegend(0.5,0.70,0.90,0.90);
  box2->SetFillColor(0);
  box2->SetBorderSize(0);
  box2->AddEntry(dummyH,"Simulation (GG2ZZ)  ","pe");
  box2->AddEntry(dummyH,lname,"");
  box2->AddEntry(dummyH,"","");
  box2->AddEntry(dummyF,lname2,"l");

  TPaveText *pt = new TPaveText(0.15,0.955,0.4,0.99,"NDC");
  pt->SetFillColor(0);
  pt->SetBorderSize(0);
  pt->AddText("CMS Preliminary 2012");
  TPaveText *pt2 = new TPaveText(0.84,0.955,0.99,0.99,"NDC");
  pt2->SetFillColor(0);
  pt2->SetBorderSize(0);

  // Plot m4l and 
  RooPlot* frameM4l = ZZMass->frame(Title("M4L"),Bins(200)) ;
  set->plotOn(frameM4l, MarkerStyle(24)) ;
  bkg_ggzz->plotOn(frameM4l,LineColor(iLineColor)) ;
  set->plotOn(frameM4l) ;

  //comaprison with different shape, if needed (uncommenting also the code above)
  //bkg_ggzz_bkgd->plotOn(frameM4l,LineColor(1),NormRange("largerange")) ;

  frameM4l->GetXaxis()->SetTitle("m_{4l} [GeV]");
  frameM4l->GetYaxis()->SetTitle("a.u.");
  //frameM4l->GetYaxis()->SetRangeUser(0,0.03);
  //if(channel == 3)frameM4l->GetYaxis()->SetRangeUser(0,0.05);
  //if(VBFtag<2){
  //  if(channel == 3)frameM4l->GetYaxis()->SetRangeUser(0,0.01);
  //  else frameM4l->GetYaxis()->SetRangeUser(0,0.005);
  //}
  frameM4l->GetXaxis()->SetRangeUser(100,1000);
  TCanvas *c = new TCanvas("c","c",800,600);
  c->cd();
  frameM4l->Draw();
  box2->Draw();
  pt->Draw();
  pt2->Draw();

  TString outputPath = "bkgFigs";
  outputPath = outputPath+ (long) sqrts + "TeV/";
  TString outputName;
  outputName =  outputPath + "bkgggzz_" + schannel + "_" + Form("%d",int(VBFtag));
  c->SaveAs(outputName + ".eps");
  c->SaveAs(outputName + ".png");
  c->SaveAs(outputName + ".root");
  delete c;

  frameM4l->GetXaxis()->SetRangeUser(100,200);
  TCanvas *c = new TCanvas("c","c",800,600);
  c->cd();
  frameM4l->Draw();
  box2->Draw();
  pt->Draw();
  pt2->Draw();

  TString outputPath = "bkgFigs";
  outputPath = outputPath+ (long) sqrts + "TeV/";
  TString outputName;
  outputName =  outputPath + "bkgggzz_lowZoom_" + schannel + "_" + Form("%d",int(VBFtag));
  c->SaveAs(outputName + ".eps");
  c->SaveAs(outputName + ".png");
  c->SaveAs(outputName + ".root");
  delete c;
} 

