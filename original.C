#include "TH1D.h"
#include "TFile.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/MultiRootFinder.h"
#include <math.h>
#include "TMath.h"
#include "TCanvas.h"
#include "TComplex.h"
#include "TF2.h"
#include "TStyle.h"
double s3 = sqrt(3.);
double s5 = sqrt(5.);
double s6 = sqrt(6.);
double s15 = sqrt(15.);
double s10 = sqrt(10.);
double s30 = sqrt(30.);
double s12 = sqrt(1./2.);
double s32 = sqrt(3./2.);
double s53 = sqrt(5./3.);
double s35 = sqrt(3./5.);
Double_t pi = TMath::Pi();


//void moments_from_waves(double* xw, double* xt){
//  double s0  = xw[ 0];
//  double p0  = xw[ 1];
//  double pm  = xw[ 2];
//  double d0  = xw[ 3];
//  double dm  = xw[ 4];
//  double pp  = xw[ 5];
//  double dp  = xw[ 6];
//  double ap0 = xw[ 7];
//  double apm = xw[ 8];
//  double ad0 = xw[ 9];
//  double adm = xw[10];
//  double adp = xw[11];
//  double s0s0 = s0*s0;
//  double p0p0 = p0*p0;
//  double pmpm = pm*pm;
//  double d0d0 = d0*d0;
//  double dmdm = dm*dm;
//  double pppp = pp*pp;
//  double dpdp = dp*dp;
//  double s0p0 = 2*s0*p0*cos(ap0);
//  double p0d0 = 2*p0*d0*cos(ap0-ad0);
//  double pmdm = 2*pm*dm*cos(apm-adm);
//  double s0pm = 2*s0*pm*cos(apm);
//  double p0dm = 2*p0*dm*cos(ap0-adm);
//  double pmd0 = 2*pm*d0*cos(apm-ad0);
//  double s0d0 = 2*s0*d0*cos(ad0);
//  double s0dm = 2*s0*dm*cos(adm);
//  double p0pm = 2*p0*pm*cos(ap0-apm);
//  double d0dm = 2*d0*dm*cos(ad0-adm);
//  double ppdp = 2*pp*dp*cos(adp);
//  double t00 = s0s0 + p0p0 + pmpm + d0d0 + dmdm + pppp + dpdp;
//  double t10 = 1/s3*s0p0 + 2/s15*p0d0 + 1/s5*(pmdm+ppdp);
//  double t11 = 1/s6*s0pm + 1/s10*p0dm - 1/s30*pmd0;
//  double t20 = 1/s5*s0d0 + 2/5.*p0p0 - 1/5.*(pmpm+pppp) + 2/7.*d0d0 + 1/7.*(dmdm+dpdp);
//  double t21 = 1/s10*s0dm+ 1/5.*s32*p0pm + 1/7.*s12*d0dm;
//  double t22 = 1/5.*s32*(pmpm-pppp) + 1/7.*s32*(dmdm-dpdp);
//  double t30 = 3/7./s5*(s3*p0d0-(pmdm+ppdp));
//  double t31 = 1/7.*s35*(2*p0dm+s3*pmd0);
//  double t32 = 1/7.*s32*(pmdm-ppdp);
//  double t40 = 2/7.*d0d0-4/21.*(dmdm+dpdp);
//  double t41 = 1/7.*s53*d0dm;
//  double t42 = s10/21.*(dmdm-dpdp);
//  xt[ 0] = t00;
//  xt[ 1] = t10;
//  xt[ 2] = t11;
//  xt[ 3] = t20;
//  xt[ 4] = t21;
//  xt[ 5] = t22;
//  xt[ 6] = t30;
//  xt[ 7] = t31;
//  xt[ 8] = t32;
//  xt[ 9] = t40;
//  xt[10] = t41;
//  xt[11] = t42;
//}

double t00;
double t10;
double t11;
double t20;
double t21;
double t22;
double t30;
double t31;
double t32;
double t33;
double t40;
double t41;
double t42;
double t43;
double t44;

double fw(const double* xw){
  double s0  = xw[0];
  double p0  = xw[1];
  double d0  = xw[2];
  double ap0 = xw[3];
  double apm = xw[4];
  double ad0 = xw[5];
  double adm = xw[6];
  double s0s0 = s0*s0;
  double p0p0 = p0*p0;
  double d0d0 = d0*d0;
  double p0d0 = 2*p0*d0*cos(ap0-ad0);
  double dmdm = (21/s10*t42 + 21/4.*(-t40 + 2/7.*d0d0))/2;
  //printf("%f\n",dmdm);
  double pmdm = (7/s32*t32 + s3*p0d0 - t30*7/3.*s5)/2.;
  if (dmdm<0) dmdm=0;
  double dm   = sqrt(dmdm);
  double pm   = pmdm/dm/2./cos(apm-adm);
  double pmpm = pm*pm;
  double s0p0 = 2*s0*p0*cos(ap0);
  double s0pm = 2*s0*pm*cos(apm);
  double s0d0 = 2*s0*d0*cos(ad0);
  double s0dm = 2*s0*dm*cos(adm);
  double p0dm = 2*p0*dm*cos(ap0-adm);
  double pmd0 = 2*pm*d0*cos(apm-ad0);
  double p0pm = 2*p0*pm*cos(ap0-apm);
  double d0dm = 2*d0*dm*cos(ad0-adm);
  
  double f00 = -t00 + s0s0 + p0p0 + 5/2.*d0d0 + 2*pmpm + 5./s32*(3*s32/s10*t42 - t22) - 21/4.*t40;
  double f10 = -t10 + 1/s3*s0p0 + s53*p0d0 - 7/3.*t30;
  double f11 = -t11 + 1/s6*s0pm + 1/s10*p0dm - 1/s30*pmd0;
  double f20 = -t20 + 1/s5*s0d0 + 2/5.*p0p0 - 1/5.*(2*pmpm + 5./s32*(3*s32/s10*t42 - t22)) + 1/2.*d0d0 - 3/4.*t40;
  double f21 = -t21 + 1/s10*s0dm+ 1/5.*s32*p0pm + 1/7.*s12*d0dm;
  double f31 = -t31 + 1/7.*s35*(2*p0dm+s3*pmd0);
  double f41 = -t41 + 1/7.*s53*d0dm;
  double sum2 = f00*f00+f10*f10+f11*f11+f20*f20+f21*f21+f31*f31+f41*f41;
  //printf("sum2=%f\n",sum2);
  return sum2;
}

double xt[15];

TComplex a0,a1,a2,a3,a4;

double fImG(const double* x){
  TComplex u(x[0],x[1],1);
  TComplex G = a4*u*u*u*u-a3*u*u*u+a2*u*u-a1*u+a0;
  return G.Im();
}

double fReG(const double* x){
  TComplex u(x[0],x[1],1);
  TComplex G = a4*u*u*u*u-a3*u*u*u+a2*u*u-a1*u+a0;
  return G.Re();
}

double f2RhoG(double* x, double* p=0){
  TComplex u(x[0],x[1],1);
  TComplex G = a4*u*u*u*u-a3*u*u*u+a2*u*u-a1*u+a0;
  double rho = G.Rho();
  if (rho>3000) return 3000;
  return rho;
}

void calc(TComplex u1, TComplex u2, TComplex u3, TComplex u4, TComplex &S0, TComplex &P0, TComplex &Pm, TComplex &Pp, TComplex &D0, TComplex &Dm, TComplex &Dp){
  S0 = 1/6.*a4*(2.*u1*u2*u3*u4+u1*u2+u1*u3+u1*u4+u2*u3+u2*u4+u3*u4+TComplex(2));
  P0 = 1/2./s3*a4*(u1*u2*u3*u4-TComplex(1));
  D0 = 1/6./s5*a4*(u1*u2*u3*u4-u1*u2-u1*u3-u1*u4-u2*u3-u2*u4-u3*u4+TComplex(1));
  Pm = 1/4./s3*a4*(u1*u2*u3+u2*u3*u4+u3*u4*u1+u4*u1*u2+u1+u2+u3+u4);
  Dm = 1/4./s15*a4*(u1*u2*u3+u2*u3*u4+u3*u4*u1+u4*u1*u2-u1-u2-u3-u4);
  // Rotate all
  double s0t = S0.Theta();
  double p0t = P0.Theta();
  double d0t = D0.Theta();
  double pmt = Pm.Theta();
  double dmt = Dm.Theta();
  S0 = TComplex(S0.Rho(),s0t-s0t,1);
  P0 = TComplex(P0.Rho(),p0t-s0t,1);
  D0 = TComplex(D0.Rho(),d0t-s0t,1);
  Pm = TComplex(Pm.Rho(),pmt-s0t,1);
  Dm = TComplex(Dm.Rho(),dmt-s0t,1);
  double dmdm = Dm.Rho2();
  double dmdm_m_dpdp = 21./s10*t42;
  double dpdp = dmdm-dmdm_m_dpdp;
  if (dpdp<0) { Dp = TComplex(0,0); Pp = TComplex(0,0); return; }
  double dp = sqrt(dpdp);
  Dp = TComplex(dp,0,1);
  double pmpm_m_pppp = 5./s32*(t22 - 1/7.*s32*dmdm_m_dpdp);
  double pmpm = Pm.Rho2();
  double pppp = pmpm-pmpm_m_pppp;
  if (pppp<0) { Pp = TComplex(0,0); return; }
  double pmdm_m_ppdp = 7./s32*t32;
  double pmdm = Pm*TComplex::Conjugate(Dm)+TComplex::Conjugate(Pm)*Dm;
  double ppdp = pmdm-pmdm_m_ppdp;
  if (ppdp*ppdp>pppp*dpdp) { Pp = TComplex(0,0); return; }
  double pp = sqrt(pppp);
  Dp = TComplex(dp,acos(ppdp/pp/dp),1);
  Pp = TComplex(pp,0,1);
}


void original(){
  gStyle->SetOptStat(0);
  ROOT::Math::Minimizer* minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit","Minuit2");
  minimizer->SetMaxFunctionCalls(1000000);
  minimizer->SetMaxIterations(10000);
  minimizer->SetTolerance(0.001);
  minimizer->SetPrintLevel(0);
  ROOT::Math::Functor f(&fw,7);
  minimizer->SetFunction(f);
  
  TFile* fWaves = new TFile("waves.root");
  TH1D* hS0m2 = (TH1D*) fWaves->Get("hS0m2");
  TH1D* hD0m2 = (TH1D*) fWaves->Get("hD0m2");
  TH1D* hD1m2 = (TH1D*) fWaves->Get("hD1m2");
  TH1D* hD1p2 = (TH1D*) fWaves->Get("hD1p2");
  TH1D* hPhi0 = (TH1D*) fWaves->Get("hPhi0");
  TH1D* hPhi1 = (TH1D*) fWaves->Get("hPhi1");
  
  TFile* fTLM = new TFile("moments.root");
  TH1D* hTLM[15];
  for (Int_t i=0;i<15;i++){
    hTLM[i] = (TH1D*) fTLM->Get(Form("hTLM%i",i));
  }

  TCanvas* c1 = new TCanvas("c1","",1500,800);
  c1->Divide(5,3);
  for (Int_t i=0;i<15;i++){
    c1->cd(i+1);
    hTLM[i]->Draw();
  }
  
  ROOT::Math::MultiRootFinder rf(ROOT::Math::MultiRootFinder::kHybridS);
  ROOT::Math::Functor f1(&fReG, 2);
  ROOT::Math::Functor f2(&fImG, 2);
  std::vector<ROOT::Math::IMultiGenFunction*> funlist;
  funlist.push_back(&f1);
  funlist.push_back(&f2);
  rf.SetFunctionList(funlist.begin(),funlist.end());
  rf.SetPrintLevel(0);
  double x1[] = {0.616,3.59};
  double x2[] = {1.170,3.59};
  double x3[] = {0.840,5.90};
  double x4[] = {1.640,5.78};

  TString wname[12] = {"|S_{0}|^{2}","|P_{0}|^{2}","|P_{-}|^{2}","|P_{+}|^{2}","|D_{0}|^{2}","|D_{-}|^{2}","|D_{+}|^{2}","arg P_{0}","arg P_{-}","arg D_{0}","arg D_{-}","arg D_{+}"}; 
  TH1D* hX[8][12];
  hX[0][0] = (TH1D*) hS0m2->Clone(Form("hX%i%02i",0,0));
  hX[0][0]->Reset();
  for (Int_t s=0;s<8;s++){
    for (Int_t i=0;i<12;i++) {
      hX[s][i] = (TH1D*) hX[0][0]->Clone(Form("hX%i%02i",s,i));
      hX[s][i]->SetTitle(wname[i].Data());
    }
  }
  Int_t nm = hTLM[0]->GetNbinsX();
  Int_t l[15]={0,1,1,2,2,2,3,3,3,3,4,4,4,4,4};
  
  for (Int_t im=nm-1;im>=0;im--){
    // if (im<nm-9) continue;
    // printf("Mass bin: %i\n",im);
    
    for (Int_t i=0;i<15;i++){
      xt[i] = sqrt(4.*TMath::Pi())/sqrt(2*l[i]+1)*hTLM[i]->GetBinContent(im+1);
    }
    
    t00 = xt[ 0];
    t10 = xt[ 1];
    t11 = xt[ 2];
    t20 = xt[ 3];
    t21 = xt[ 4];
    t22 = xt[ 5];
    t30 = xt[ 6];
    t31 = xt[ 7];
    t32 = xt[ 8];
    t33 = xt[ 9];
    t40 = xt[10];
    t41 = xt[11];
    t42 = xt[12];
    t43 = xt[13];
    t44 = xt[14];
    Double_t s0m2 = hS0m2->GetBinContent(im+1);
    Double_t d0m2 = hD0m2->GetBinContent(im+1);
    Double_t d1m2 = hD1m2->GetBinContent(im+1);
    Double_t d1p2 = hD1p2->GetBinContent(im+1);
    Double_t phi0 = hPhi0->GetBinContent(im+1);
    Double_t phi1 = hPhi1->GetBinContent(im+1);
    Double_t p0m2 = t00-d0m2-d0m2-d1m2-d1p2;
    if (p0m2<0) p0m2=10;
    Double_t s0m = sqrt(s0m2);
    Double_t d0m = sqrt(d0m2);
    Double_t d1m = sqrt(d1m2);
    Double_t d1p = sqrt(d1p2);
    Double_t p0m = sqrt(p0m2);
    Double_t xw[7];
    xw[0] = s0m;
    xw[1] = p0m;
    xw[2] = d0m;
    xw[3] = 0;
    xw[4] = 0;
    xw[5] = phi0;
    xw[6] = phi1;
    //printf("%f\n",fw(xw));
    //return;

    minimizer->SetLimitedVariable(0,"s0",s0m,100,0,1e7);
    minimizer->SetLimitedVariable(1,"p0",p0m,100,0,1e7);
    minimizer->SetLimitedVariable(2,"d0",d0m,100,0,1e7);
    minimizer->SetLimitedVariable(3,"ap0",0,0.01,-2*pi,2*pi);
    minimizer->SetLimitedVariable(4,"apm",0,0.01,-2*pi,2*pi);
    minimizer->SetLimitedVariable(5,"ad0",phi0,0.01,-2*pi,2*pi);
    minimizer->SetLimitedVariable(6,"adm",phi1,0.01,-2*pi,2*pi);
    minimizer->Minimize();
    Int_t status = minimizer->Status();
    //printf("Status: %i\n",status);
    const Double_t* xval = minimizer->X();
    const Double_t* xerr = minimizer->Errors();
    // for (Int_t i=0;i<7;i++) hX[i]->SetBinContent(im+1,xval[i]);
    //printf("%.0f %.0f %.0f\n",xval[0],xval[1],xval[2]);
    double s0  = xval[0];
    double p0  = xval[1];
    double d0  = xval[2];
    double ap0 = xval[3];
    double apm = xval[4];
    double ad0 = xval[5];
    double adm = xval[6];
    if (ap0<0) ap0+=2*pi;
    if (apm<0) apm+=2*pi;
    if (ad0<0) ad0+=2*pi;
    if (adm<0) adm+=2*pi;
    double s0s0 = s0*s0;
    double p0p0 = p0*p0;
    double d0d0 = d0*d0;
    double p0d0 = 2*p0*d0*cos(ap0-ad0);
    double dmdm = (21/s10*t42 + 21/4.*(-t40 + 2/7.*d0d0))/2;
    double pmdm = (7/s32*t32 + s3*p0d0 - t30*7/3.*s5)/2.;
    if (dmdm<0) dmdm=0;
    double dm   = sqrt(dmdm);
    double pm   = pmdm/dm/2./cos(apm-adm);
    if (pm<0) pm=-pm;
    
    TComplex S0(s0,0.,1);
    TComplex P0(p0,ap0,1);
    TComplex PM(pm,apm,1);
    TComplex D0(d0,ad0,1);
    TComplex DM(dm,adm,1);
    
    a4 = S0 - s3*P0 + s5*D0;
    a3 = 2.*s3*(PM-s5*DM);
    a2 = 2.*S0-4.*s5*D0;
    a1 = 2.*s3*(PM+s5*DM);
    a0 = S0 + s3*P0 + s5*D0;
    bool ret;
    
    printf("%.04f %.04f\n",x1[0],x1[1]);
    rf.Solve(x1);
    printf("%.04f %.04f\n",rf.X()[0],rf.X()[1]);
    TComplex u1(rf.X()[0],rf.X()[1],1);
    
    printf("%.04f %.04f\n",x2[0],x2[1]);
    rf.Solve(x2);
    printf("%.04f %.04f\n",rf.X()[0],rf.X()[1]);
    TComplex u2(rf.X()[0],rf.X()[1],1);
    
    printf("%.04f %.04f\n",x3[0],x3[1]);
    rf.Solve(x3);
    printf("%.04f %.04f\n",rf.X()[0],rf.X()[1]);
    TComplex u3(rf.X()[0],rf.X()[1],1);
    
    printf("%.04f %.04f\n",x4[0],x4[1]);
    rf.Solve(x4);
    printf("%.04f %.04f\n",rf.X()[0],rf.X()[1]);
    TComplex u4(rf.X()[0],rf.X()[1],1);
    
    x1[0]=u1.Rho(); x1[1]=u1.Theta();
    x2[0]=u2.Rho(); x2[1]=u2.Theta();
    x3[0]=u3.Rho(); x3[1]=u3.Theta();
    x4[0]=u4.Rho(); x4[1]=u4.Theta();
    
    TComplex u2c = TComplex::Conjugate(u2);
    TComplex u3c = TComplex::Conjugate(u3);
    TComplex u4c = TComplex::Conjugate(u4);
    
    TComplex vS0[8];
    TComplex vP0[8];
    TComplex vPm[8];
    TComplex vPp[8];
    TComplex vD0[8];
    TComplex vDm[8];
    TComplex vDp[8];
    calc(u1,u2 ,u3 ,u4 ,vS0[0],vP0[0],vPm[0],vPp[0],vD0[0],vDm[0],vDp[0]);
    calc(u1,u2c,u3 ,u4 ,vS0[1],vP0[1],vPm[1],vPp[1],vD0[1],vDm[1],vDp[1]);
    calc(u1,u2 ,u3c,u4 ,vS0[2],vP0[2],vPm[2],vPp[2],vD0[2],vDm[2],vDp[2]);
    calc(u1,u2 ,u3 ,u4c,vS0[3],vP0[3],vPm[3],vPp[3],vD0[3],vDm[3],vDp[3]);
    calc(u1,u2c,u3c,u4 ,vS0[4],vP0[4],vPm[4],vPp[4],vD0[4],vDm[4],vDp[4]);
    calc(u1,u2 ,u3c,u4c,vS0[5],vP0[5],vPm[5],vPp[5],vD0[5],vDm[5],vDp[5]);
    calc(u1,u2c,u3 ,u4c,vS0[6],vP0[6],vPm[6],vPp[6],vD0[6],vDm[6],vDp[6]);
    calc(u1,u2c,u3c,u4c,vS0[7],vP0[7],vPm[7],vPp[7],vD0[7],vDm[7],vDp[7]);
    for (Int_t s=0;s<8;s++){
      hX[s][ 0]->SetBinContent(im+1,vS0[s].Rho2());
      hX[s][ 1]->SetBinContent(im+1,vP0[s].Rho2());
      hX[s][ 2]->SetBinContent(im+1,vPm[s].Rho2());
      hX[s][ 3]->SetBinContent(im+1,vPp[s].Rho2());
      hX[s][ 4]->SetBinContent(im+1,vD0[s].Rho2());
      hX[s][ 5]->SetBinContent(im+1,vDm[s].Rho2());
      hX[s][ 6]->SetBinContent(im+1,vDp[s].Rho2());
      hX[s][ 7]->SetBinContent(im+1,vP0[s].Theta()>0 ? vP0[s].Theta() : vP0[s].Theta()+2*pi);
      hX[s][ 8]->SetBinContent(im+1,vPm[s].Theta()>0 ? vPm[s].Theta() : vPm[s].Theta()+2*pi);
      hX[s][ 9]->SetBinContent(im+1,vD0[s].Theta()>0 ? vD0[s].Theta() : vD0[s].Theta()+2*pi);
      hX[s][10]->SetBinContent(im+1,vDm[s].Theta()>0 ? vDm[s].Theta() : vDm[s].Theta()+2*pi);
      hX[s][11]->SetBinContent(im+1,vDp[s].Theta()>0 ? vDp[s].Theta() : vDp[s].Theta()+2*pi);
    }
  }
  
  TCanvas* c2 = new TCanvas("c2","c2",1800,800);
  c2->Divide(4,3,0.001,0.001);
  Float_t color[8] = {kBlack,kRed,kMagenta,kGray,kGreen+1,kYellow+1,kCyan,kBlue+1};
  for (Int_t i=0;i<12;i++){
    c2->cd(i+1);
    gPad->SetRightMargin(0.01);
    for (Int_t s=0;s<8;s++){
      hX[s][i]->SetMaximum(i<7 ? 3e7 : 2*pi);
      hX[s][i]->SetLineColor(color[s]);
      hX[s][i]->Draw(s==0 ? "e": "e same");
    }
    if (i== 0) hS0m2->Draw("same");
    if (i== 4) hD0m2->Draw("same");
    if (i== 5) hD1m2->Draw("same");
    if (i== 6) hD1p2->Draw("same");
    if (i== 9) hPhi0->Draw("same");
    if (i==10) hPhi1->Draw("same");
  }

  new TCanvas;
  TF2* ff2RhoG = new TF2("ff2RhoG",&f2RhoG,0,3,0,2*pi);
  ff2RhoG->SetNpx(1000);
  ff2RhoG->SetNpy(1000);
  ff2RhoG->Draw("colz");
//  
//
}
