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
#include "TGraph.h"
#include "fstream"

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
TComplex a0pol, a1pol, a2pol, a3pol, a4pol;

double fImG(const double* x){
  TComplex u(x[0],x[1],0);
  TComplex G = a4*u*u*u*u-a3*u*u*u+a2*u*u-a1*u+a0;
  return G.Im();
}

//--------------------------------------
// additional functions for root finding
//
double r1[2] = {0.0,0.0};
double r2[2] = {0.0,0.0};
double r3[2] = {0.0,0.0};
double r4[2] = {0.0,0.0};
//r1 -- 1st root, r2 -- 2nd root, etc.
double fImG1(const double* x){
  TComplex u(x[0],x[1],0);
  TComplex u1(r1[0], r1[1], 0);
  TComplex G = a4*u*u*u*u-a3*u*u*u+a2*u*u-a1*u+a0;
  return ( G / (u - u1) ).Im();
}

double fImG2(const double* x){
  TComplex u(x[0],x[1],0);
  TComplex u1(r1[0], r1[1], 0);
  TComplex u2(r2[0], r2[1], 0);
  TComplex G = a4*u*u*u*u-a3*u*u*u+a2*u*u-a1*u+a0;
  return ( G / ((u - u1)*(u - u2)) ).Im();
}

double fImG3(const double* x){
  TComplex u(x[0],x[1],0);
  TComplex u1(r1[0], r1[1], 0);
  TComplex u2(r2[0], r2[1], 0);
  TComplex u3(r3[0], r3[1], 0);
  TComplex G = a4*u*u*u*u-a3*u*u*u+a2*u*u-a1*u+a0;
  return ( G / ((u - u1)*(u - u2)*(u - u3)) ).Im();
}

double fReG(const double* x){
  TComplex u(x[0],x[1],0);
  TComplex G = a4*u*u*u*u-a3*u*u*u+a2*u*u-a1*u+a0;
  return G.Re();
}

double fReG1(const double* x){
  TComplex u(x[0],x[1],0);
  TComplex u1(r1[0], r1[1], 0);
  TComplex G = a4*u*u*u*u-a3*u*u*u+a2*u*u-a1*u+a0;
  return ( G / (u - u1) ).Re();
}

double fReG2(const double* x){
  TComplex u(x[0],x[1],0);
  TComplex u1(r1[0], r1[1], 0);
  TComplex u2(r2[0], r2[1], 0);
  TComplex G = a4*u*u*u*u-a3*u*u*u+a2*u*u-a1*u+a0;
  return ( G / ((u - u1)*(u - u2)) ).Re();
}

double fReG3(const double* x){
  TComplex u(x[0],x[1],0);
  TComplex u1(r1[0], r1[1], 0);
  TComplex u2(r2[0], r2[1], 0);
  TComplex u3(r3[0], r3[1], 0);
  TComplex G = a4*u*u*u*u-a3*u*u*u+a2*u*u-a1*u+a0;
  return ( G / ((u - u1)*(u - u2)*(u - u3)) ).Re();
}
//-------------------------------------

double f2RhoG(double* x, double* p=0){
  TComplex u(x[0],x[1],0);
  TComplex G = a4*u*u*u*u-a3*u*u*u+a2*u*u-a1*u+a0;
  double rho = G.Rho();
  if (rho>0.1) 
    return 0.1;
  return rho;
}

void calc(TComplex u1, TComplex u2, TComplex u3, TComplex u4, TComplex &S0, TComplex &P0, TComplex &Pm, TComplex &Pp, TComplex &D0, TComplex &Dm, TComplex &Dp){
  S0 = 1/6.*a4pol*(2.*u1*u2*u3*u4+u1*u2+u1*u3+u1*u4+u2*u3+u2*u4+u3*u4+TComplex(2));
  P0 = 1/2./s3*a4pol*(u1*u2*u3*u4-TComplex(1));
  D0 = 1/6./s5*a4pol*(u1*u2*u3*u4-u1*u2-u1*u3-u1*u4-u2*u3-u2*u4-u3*u4+TComplex(1));
  Pm = 1/4./s3*a4pol*(u1*u2*u3+u2*u3*u4+u3*u4*u1+u4*u1*u2+u1+u2+u3+u4);
  Dm = 1/4./s15*a4pol*(u1*u2*u3+u2*u3*u4+u3*u4*u1+u4*u1*u2-u1-u2-u3-u4);
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

void checkRe(TComplex vS0, TComplex vP0, TComplex vPm, TComplex vPp, TComplex vD0, TComplex vDm, TComplex vDp, double *R) {

  double s0  = vS0.Rho();
  double p0  = vP0.Rho();
  double pm  = vPm.Rho();
  double d0  = vD0.Rho();
  double dm  = vDm.Rho();
  double pp  = vPp.Rho();
  double dp  = vDp.Rho();
  double ap0 = vP0.Theta();
  double apm = vPm.Theta();
  double ad0 = vD0.Theta();
  double adm = vDm.Theta();
  double adp = vDp.Theta();
  double s0s0 = s0*s0;
  double p0p0 = p0*p0;
  double pmpm = pm*pm;
  double d0d0 = d0*d0;
  double dmdm = dm*dm;
  double pppp = pp*pp;
  double dpdp = dp*dp;
  double s0p0 = 2*s0*p0*cos(ap0);
  double p0d0 = 2*p0*d0*cos(ap0-ad0);
  double pmdm = 2*pm*dm*cos(apm-adm);
  double s0pm = 2*s0*pm*cos(apm);
  double p0dm = 2*p0*dm*cos(ap0-adm);
  double pmd0 = 2*pm*d0*cos(apm-ad0);
  double s0d0 = 2*s0*d0*cos(ad0);
  double s0dm = 2*s0*dm*cos(adm);
  double p0pm = 2*p0*pm*cos(ap0-apm);
  double d0dm = 2*d0*dm*cos(ad0-adm);
  double ppdp = 2*pp*dp*cos(adp);

  R[0] =  -t00 + s0s0 + p0p0 + pmpm + d0d0 + dmdm + pppp + dpdp;
  R[1] =  -t10 + 1/s3*s0p0 + 2/s15*p0d0 + 1/s5*(pmdm+ppdp);
  R[2] =  -t11 + 1/s6*s0pm + 1/s10*p0dm - 1/s30*pmd0;
  R[3] =  -t20 + 1/s5*s0d0 + 2/5.*p0p0 - 1/5.*(pmpm+pppp) + 2/7.*d0d0 + 1/7.*(dmdm+dpdp);
  R[4] =  -t21 + 1/s10*s0dm+ 1/5.*s32*p0pm + 1/7.*s12*d0dm;
  R[5] =  -t22 + 1/5.*s32*(pmpm-pppp) + 1/7.*s32*(dmdm-dpdp);
  R[6] =  -t30 + 3/7./s5*(s3*p0d0-(pmdm+ppdp));
  R[7] =  -t31 + 1/7.*s35*(2*p0dm+s3*pmd0);
  R[8] =  -t32 + 1/7.*s32*(pmdm-ppdp);
  R[9] =  -t40 + 2/7.*d0d0-4/21.*(dmdm+dpdp);
  R[10] = -t41 + 1/7.*s53*d0dm;
  R[11] = -t42 + s10/21.*(dmdm-dpdp);

  // R[0] =  -t00 + vS0*vS0 + vP0*vP0 + vPm*vPm + vD0*vD0 + vDm*vDm + vPp*vPp + vDp*vDp;
  // R[1] =  -t10 + TComplex(1)/s3*vS0*vS0 + TComplex(2)/s15*vP0*vD0 + TComplex(1)/s5*(vPm*vDm+vPp*vDp);
  // R[2] =  -t11 + TComplex(1)/s6*vS0*vPm + TComplex(1)/s10*vP0*vDm - TComplex(1)/s30*vPm*vD0;
  // R[3] =  -t20 + TComplex(1)/s5*vS0*vD0 + TComplex(2/5.)*vP0*vP0 - TComplex(1/5.)*(vPm*vPm+vPp*vPp) + TComplex(2/7.)*vD0*vD0 + TComplex(1/7.)*(vDm*vDm+vDp*vDp);
  // R[4] =  -t21 + TComplex(1)/s10*vS0*vDm+ TComplex(1/5.)*s32*vP0*vPm + TComplex(1/7.)*s12*vD0*vDm;
  // R[5] =  -t22 + TComplex(1/5.)*s32*(vPm*vPm-vPp*vPp) + TComplex(1/7.)*s32*(vDm*vDm-vDp*vDp);
  // R[6] =  -t30 + TComplex(3/7.)/s5*(s3*vP0*vD0-(vPm*vDm+vPp*vDp));
  // R[7] =  -t31 + TComplex(1/7.)*s35*(TComplex(2)*vP0*vDm+s3*vPm*vD0);
  // R[8] =  -t32 + TComplex(1/7.)*s32*(vPm*vDm-vPp*vDp);
  // R[9] =  -t40 + TComplex(2/7.)*vD0*vD0-TComplex(4/21.)*(vDm*vDm+vDp*vDp);
  // R[10] = -t41 + TComplex(1/7.)*s53*vD0*vDm;
  // R[11] = -t42 + s10/TComplex(21.)*(vDm*vDm-vDp*vDp);
}

void waves_from_moments(){
  gStyle->SetOptStat(0);
  ROOT::Math::Minimizer* minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit","Minuit2");
  minimizer->SetMaxFunctionCalls(1000000);
  minimizer->SetMaxIterations(10000);
  minimizer->SetTolerance(0.00001);
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

  TString wnameMoments[12] = {"t00","t10","t11","t20","t21","t22","t30","t31","t32","t40","t41","t42"}; 
  TH1D* hMoments[8][12];
  hMoments[0][0] = (TH1D*) hS0m2->Clone(Form("hX%i%02i",0,0));
  hMoments[0][0]->Reset();
  for (Int_t s=0;s<8;s++){
    for (Int_t i=0;i<12;i++) {
      hMoments[s][i] = (TH1D*) hMoments[0][0]->Clone(Form("hX%i%02i",s,i));
      hMoments[s][i]->SetTitle(wnameMoments[i].Data());
    }
  }
  
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

//additional root finders
  //finding 2nd root, excluding 1st root
  ROOT::Math::MultiRootFinder rf1(ROOT::Math::MultiRootFinder::kHybridS);
  ROOT::Math::Functor f11(&fReG1, 2);
  ROOT::Math::Functor f21(&fImG1, 2);
  std::vector<ROOT::Math::IMultiGenFunction*> funlist1;
  funlist1.push_back(&f11);
  funlist1.push_back(&f21);
  rf1.SetFunctionList(funlist1.begin(),funlist1.end());
  rf1.SetPrintLevel(0);

  //finding 3rd root, excluding 1st and 2nd
  ROOT::Math::MultiRootFinder rf2(ROOT::Math::MultiRootFinder::kHybridS);
  ROOT::Math::Functor f12(&fReG2, 2);
  ROOT::Math::Functor f22(&fImG2, 2);
  std::vector<ROOT::Math::IMultiGenFunction*> funlist2;
  funlist2.push_back(&f12);
  funlist2.push_back(&f22);
  rf2.SetFunctionList(funlist2.begin(),funlist2.end());
  rf2.SetPrintLevel(0);

  //finding 4th root, excluding 1st, 2nd and 3rd
  ROOT::Math::MultiRootFinder rf3(ROOT::Math::MultiRootFinder::kHybridS);
  ROOT::Math::Functor f13(&fReG3, 2);
  ROOT::Math::Functor f23(&fImG3, 2);
  std::vector<ROOT::Math::IMultiGenFunction*> funlist3;
  funlist3.push_back(&f13);
  funlist3.push_back(&f23);
  rf3.SetFunctionList(funlist3.begin(),funlist3.end());
  rf3.SetPrintLevel(0);
//--------------------

  double x1[] = {-0.9, 0.7}; //{0.616,3.59};
  double x2[] = { 0.8, 0.7}; //{1.170,3.59};
  double x3[] = { 1.5,-1.4}; //{0.840,5.90};
  double x4[] = {-0.3,-0.5}; //{1.640,5.78};

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

  std::vector<double> vx1, vx2, vx3, vx4;
  std::vector<double> vy1, vy2, vy3, vy4;

  for (Int_t im=nm-1;im>=0;im--){
    // if (im<nm-30)
    //   continue;
    printf("Mass bin: %i\n",im);
    
    //-----------------
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
    //-----------------

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
    //if (dmdm<0) dmdm=0;
    double dm   = sqrt(dmdm);
    double pm   = pmdm/dm/2./cos(apm-adm);
    //if (pm<0) pm=-pm;
    
    TComplex S0(s0,0., 1);
    TComplex P0(p0,ap0,1);
    TComplex PM(pm,apm,1);
    TComplex D0(d0,ad0,1);
    TComplex DM(dm,adm,1);
    
    a4pol = S0 - s3*P0 + s5*D0;
    a3pol = 2.*s3*(PM-s5*DM);
    a2pol = 2.*S0-4.*s5*D0;
    a1pol = 2.*s3*(PM+s5*DM);
    a0pol = S0 + s3*P0 + s5*D0;
    
    a4 = TComplex(a4pol.Re(), a4pol.Im(), 0);
    a3 = TComplex(a3pol.Re(), a3pol.Im(), 0);
    a2 = TComplex(a2pol.Re(), a2pol.Im(), 0);
    a1 = TComplex(a1pol.Re(), a1pol.Im(), 0);
    a0 = TComplex(a0pol.Re(), a0pol.Im(), 0);
    
    std::cout<<"-----------\n";
    printf("1stB %.04f %.04f\n",x1[0],x1[1]);
    rf.Solve(x1, 100000, 1e-6);
    r1[0] = rf.X()[0];
    r1[1] = rf.X()[1];
    printf("1stA %.04f %.04f\n",rf.X()[0],rf.X()[1]);
    TComplex u1cor(rf.X()[0],rf.X()[1],0);
    
    printf("2ndB %.04f %.04f\n",x2[0],x2[1]);
    rf1.Solve(x2, 100000, 1e-10); // x1 excluded, getting approximate root
    r2[0] = rf1.X()[0];
    r2[1] = rf1.X()[1];
    rf.Solve(r2, 100000, 1e-10); // starting from approx. root
    printf("2ndA %.04f %.04f\n",rf.X()[0],rf.X()[1]);
    TComplex u2cor(rf.X()[0],rf.X()[1],0);
    
    printf("3rdB %.04f %.04f\n",x3[0],x3[1]);
    rf2.Solve(x3, 100000, 1e-10);
    r3[0] = rf2.X()[0];
    r3[1] = rf2.X()[1];
    rf.Solve(r3, 100000, 1e-10);
    printf("3rdA %.04f %.04f\n",rf.X()[0],rf.X()[1]);
    TComplex u3cor(rf.X()[0],rf.X()[1],0);
    
    printf("4thB %.04f %.04f\n",x4[0],x4[1]);
    rf3.Solve(x4, 100000, 1e-10);
    r4[0] = rf3.X()[0];
    r4[1] = rf3.X()[1];
    rf.Solve(r4, 100000, 1e-10);
    printf("4thA %.04f %.04f\n",rf.X()[0],rf.X()[1]);
    TComplex u4cor(rf.X()[0],rf.X()[1],0);
    std::cout<<"-----------\n";

//-----------------------------------
    x1[0]=u1cor.Re(); x1[1]=u1cor.Im();
    x2[0]=u2cor.Re(); x2[1]=u2cor.Im();
    x3[0]=u3cor.Re(); x3[1]=u3cor.Im();
    x4[0]=u4cor.Re(); x4[1]=u4cor.Im();

    vx1.push_back(x1[0]);vy1.push_back(x1[1]);
    vx2.push_back(x2[0]);vy2.push_back(x2[1]);
    vx3.push_back(x3[0]);vy3.push_back(x3[1]);
    vx4.push_back(x4[0]);vy4.push_back(x4[1]);
//-----------------------------------

    TComplex u1(u1cor.Rho(), u1cor.Theta(),1);
    TComplex u2(u2cor.Rho(), u2cor.Theta(),1);
    TComplex u3(u3cor.Rho(), u3cor.Theta(),1);
    TComplex u4(u4cor.Rho(), u4cor.Theta(),1);

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

      double hp0 = vP0[s].Theta()>0 ? vP0[s].Theta() : vP0[s].Theta()+2*pi;
      double hpm = vPm[s].Theta()>0 ? vPm[s].Theta() : vPm[s].Theta()+2*pi;
      double hdm = vDm[s].Theta()>0 ? vDm[s].Theta() : vDm[s].Theta()+2*pi;
      double hdp = vDp[s].Theta()>0 ? vDp[s].Theta() : vDp[s].Theta()+2*pi;
      double hd0 = vD0[s].Theta()>0 ? vD0[s].Theta() : vD0[s].Theta()+2*pi;

      if (hdm < pi) {
        hdm = 2*pi - hdm;
      }

      if (hd0 < pi) {
        hd0 = 2*pi - hd0;
      }

      hX[s][ 7]->SetBinContent(im+1,hp0);
      hX[s][ 8]->SetBinContent(im+1,hpm);
      hX[s][ 9]->SetBinContent(im+1,hd0);
      hX[s][10]->SetBinContent(im+1,hdm);
      hX[s][11]->SetBinContent(im+1,hdp);
    }

    double tt[12];

    for (int i = 0; i < 8; i++){
      checkRe(vS0[i],vP0[i],vPm[i],vPp[i],vD0[i],vDm[i],vDp[i], tt);
      hMoments[i][ 0]->SetBinContent(im+1,tt[ 0]);
      hMoments[i][ 1]->SetBinContent(im+1,tt[ 1]);
      hMoments[i][ 2]->SetBinContent(im+1,tt[ 2]);
      hMoments[i][ 3]->SetBinContent(im+1,tt[ 3]);
      hMoments[i][ 4]->SetBinContent(im+1,tt[ 4]);
      hMoments[i][ 5]->SetBinContent(im+1,tt[ 5]);
      hMoments[i][ 6]->SetBinContent(im+1,tt[ 6]);
      hMoments[i][ 7]->SetBinContent(im+1,tt[ 7]);
      hMoments[i][ 8]->SetBinContent(im+1,tt[ 8]);
      hMoments[i][ 9]->SetBinContent(im+1,tt[ 9]);
      hMoments[i][10]->SetBinContent(im+1,tt[10]);
      hMoments[i][11]->SetBinContent(im+1,tt[11]);
    }

  }

  

  TCanvas* c2 = new TCanvas("c2","c2",1900,800);
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

  TCanvas* c3 = new TCanvas("c3","c3",1900,800);
  c3->Divide(4,3,0.001,0.001);
  for (Int_t i=0;i<12;i++){
    c3->cd(i+1);
    gPad->SetRightMargin(0.01);
    for (Int_t s=0;s<8;s++){
      hMoments[s][i]->SetLineColor(color[s]);
      hMoments[s][i]->SetMaximum(1e7);
      hMoments[s][i]->SetMinimum(-1e7);
      hMoments[s][i]->Draw(s==0 ? "e": "e same");
    }
  }

  new TCanvas;
  TF2* ff2RhoG = new TF2("ff2RhoG",&f2RhoG,-3,3,-2*pi,2*pi);
  ff2RhoG->SetNpx(1000);
  ff2RhoG->SetNpy(1000);
  ff2RhoG->Draw("colz");

  std::ofstream ofile;
  ofile.open ("roots.txt");
  for (int i=0; i<vx1.size(); i++)
    ofile << vx1[i] << " " << vy1[i] << " " << vx2[i] << " " << vy2[i] << " " << vx3[i] << " " << vy3[i] << " " << vx4[i] << " " << vy4[i] << std::endl;
  ofile.close();
}
