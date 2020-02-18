#include "TTree.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TMath.h"
#include "TStopwatch.h"
#include "TStyle.h"
#include "TComplex.h"
using namespace std;
const Double_t pi = TMath::Pi();
Int_t n[15] = {1,1,2,1,2,2,1,2,2,2,1,2,2,2,2};
Double_t ne[15];

vector<double>* vvy[15];

Double_t FT(const Double_t* tlm){
  Double_t sum=0;
  for (UInt_t i=0;i<15;i++) sum+=tlm[i]*ne[i];
  
  for (UInt_t ev=0;ev<vvy[0]->size();ev++){
    Double_t i0 = 0;
    for (UInt_t i=0;i<15;i++) i0+=n[i]*tlm[i]*(*vvy[i])[ev];
    sum-=log(i0);
  }
  return sum;
}


void fit_moments(){
  ROOT::Math::Minimizer* minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit","");
  minimizer->SetMaxFunctionCalls(100000);
  minimizer->SetMaxIterations(1000);
  minimizer->SetTolerance(0.01);
  minimizer->SetPrintLevel(2);
  ROOT::Math::Functor f(&FT,15);
  minimizer->SetFunction(f);

  TFile* fEff = new TFile("fpwamc.root");
  TH1D* hEff[15];
  TH1D* hMcM = (TH1D*) fEff->Get("hMcM");
  for (Int_t i=0;i<15;i++) hEff[i] = (TH1D*) fEff->Get(Form("hRcM%i",i));
  for (Int_t i=0;i<15;i++) hEff[i]->Divide(hMcM);
  
  Int_t nm=hEff[0]->GetNbinsX();
  Double_t mmin = hEff[0]->GetXaxis()->GetXmin();
  Double_t mmax = hEff[0]->GetXaxis()->GetXmax();

  TH1D* hTLM[15];
  for (Int_t i=0;i<15;i++) {
    hTLM[i] = (TH1D*) hEff[0]->Clone(Form("hTLM%i",i));
    hTLM[i]->Reset();
  }
  
  TH1D* hM = (TH1D*) hEff[0]->Clone("hM");
  
  TFile* fpwa = new TFile("fpwadata.root");
  // TFile* fpwa = new TFile("fpwamc.root");
  TTree* tpwa = (TTree*) fpwa->Get("tpwa");
  for (Int_t i=0;i<15;i++) tpwa->SetBranchAddress(Form("y%02i",i),&(vvy[i]));
  // tpwa->SetBranchAddress("y00",&(vvy[0]));
  
  TCanvas* c1 = new TCanvas("c1","c1",2400,1000);
  c1->Divide(5,3,0.001,0.001);

  Double_t init[15]={0};
  for (Int_t i=0;i<15;i++)
    minimizer->SetVariable(i,Form("t%i",i),0,100);
  for (Int_t im=nm-1;im>=1;im--){
    if (im<nm-40) continue;
    printf("Reading data\n");
    tpwa->GetEntry(im);
    hM->SetBinContent(im+1,vvy[0]->size());
    if (hEff[0]->GetXaxis()->GetBinLowEdge(im+1)<0.8)
      continue;
    printf("Reading efficiency\n");
    for (Int_t i=0;i<15;i++)
      ne[i]=n[i]*hEff[i]->GetBinContent(im+1);
    printf("Set initial values\n");
    for (Int_t i=0;i<15;i++)
      minimizer->SetVariable(i,Form("t%i",i),0,100);
    minimizer->SetVariableValue(0,vvy[0]->size()/hEff[0]->GetBinContent(im+1));
    for (Int_t i=1;i<15;i++)
      minimizer->SetVariableValue(i,0);
    printf("Minimize\n");
    minimizer->Minimize();
    Int_t status = minimizer->Status();
    printf("Status: %i\n",status);
    const Double_t* xval = minimizer->X();
    const Double_t* xerr = minimizer->Errors();
    for (Int_t i=0;i<15;i++) {
  //    if (status==0) init[i]=xval[i];
      c1->cd(i+1);
      hTLM[i]->SetBinContent(im+1,xval[i]);
      hTLM[i]->SetBinError(im+1,xerr[i]);
      hTLM[i]->Draw();
    }
    c1->Update();
  }
  
  new TCanvas;
  hM->Draw();
  
  new TCanvas;
  hEff[0]->Draw();
  
  TFile* fmoments = new TFile("moments.root","recreate");
  for (Int_t i=0;i<15;i++) hTLM[i]->Write();
  fmoments->Close();
}

