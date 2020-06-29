#ifndef cutsAndBinZb_h
#define cutsAndBinZb_h

#include <TF1.h>
#include <TCut.h>
#include <TChain.h>
#include <TGraphAsymmErrors.h>
#include <TCanvas.h>
#include <TNtuple.h>
#include <iostream>
#include <TLine.h>
#include <TMath.h>
#include <TTree.h>
#include <TMath.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TROOT.h>
#include <TClonesArray.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>

int kMuGlb = 1;
int kMuGlbTrk = 2;
  
const int Ntrig = 4;
int kL1DoubleMuOpen = 0;
int kTrigUps = 13;
int kTrigL1DBOS40100 = 1;
int kTrigL1DB50100 = 2;
TString fTrigName[Ntrig] = {"L1DoubleMuOpen", "Ups", "L1DoubleMuOpenOS40100", "L1DoubleMuOpen50100"};

bool IsAcceptanceQQ(double pt, double eta)
{
  return ( (fabs(eta)<1.2 && pt>=3.5) ||
      (1.2<=fabs(eta) && fabs(eta)<2.1 && pt>=5.47-1.89*fabs(eta)) ||
      (2.1<=fabs(eta) && fabs(eta)<2.4 && pt>=1.5) );
}

bool IsAcceptanceNoTrig(double pt, double eta)
{
  return ( (fabs(eta)<0.3 && pt>=3.4) ||
      (0.3<=fabs(eta) && fabs(eta)<1.1 && pt>=3.3) ||
      (1.1<=fabs(eta) && fabs(eta)<1.4 && pt>=7.7-4.0*fabs(eta)) ||
      (1.4<=fabs(eta) && fabs(eta)<1.55 && pt>=2.1) ||
      (1.55<=fabs(eta) && fabs(eta)<2.2 && pt>=4.25-1.39*fabs(eta)) ||
      (2.2<=fabs(eta) && fabs(eta)<2.4 && pt>=1.2) );
}


// lumi Unc 
double lumi_unc_pp = 0.023;
double nMB_unc = TMath::Sqrt(0.02*0.02+0.01*0.01);

struct ParticleMass { double JPsi, Psi2S, Y1S, Y2S, Y3S, Z, PiPlMi, KaPlus; };
ParticleMass pdgMass = {3.096, 3.686, 9.460, 10.023, 10.355, 91.188, 0.13957061, 0.49367 };

int fpol1 = 1;
int fpol2 = 2;
int fpol3 = 3;
int fpol4 = 4;
int fExp = 10;


struct valErr { float val, err; } ; 

int kPPDATA = 0 ;
int kPADATA = 1 ;
int kAADATA = 2 ; // L1 doubleMu 0
int kPPMC = 3 ;
int kPAMC = 4 ;
int kAAMC = 5 ;
int kAADATAPeri = 6 ;
int kAADATACentL3 = 7 ;
int kPPMCUps1S = 8 ;
int kPPMCUps2S = 9 ;
int kPPMCUps3S = 10 ;
int kAAMCUps1S = 11 ;
int kAAMCUps2S = 12 ;
int kAAMCUps3S = 13 ;
int kPPAADATASIMUL = 20 ; // 2 and 0 simultaneous fit
int kPPAADATAPeriSIMUL = 60 ; // 6 and 0 simultaneous fit

TString getCollID( int collid ) {
  if ( collid == kPPDATA ) return "PP_DATA";
  else if ( collid == kPADATA ) return "PA_DATA";
  else if ( collid == kAADATA ) return "AA_DATA";
  else if ( collid == kPPMC ) return "PP_MC";
  else if ( collid == kPAMC ) return "PA_MC";
  else if ( collid == kAAMC ) return "AA_MC";
  else if ( collid == kAADATAPeri ) return "AA_DATA";
  else if ( collid == kAADATACentL3 ) return "AA_DATA_CentL3";
  else if ( collid == kPPMCUps1S ) return "PP_MC_Ups1S";
  else if ( collid == kPPMCUps2S ) return "PP_MC_Ups2S";
  else if ( collid == kPPMCUps3S ) return "PP_MC_Ups3S";
  else if ( collid == kAAMCUps1S ) return "AA_MC_Ups1S";
  else if ( collid == kAAMCUps2S ) return "AA_MC_Ups2S";
  else if ( collid == kAAMCUps3S ) return "AA_MC_Ups3S";
  else if ( collid == kPPAADATASIMUL ) return "PP_AA_DATA_SIMUL";
  else if ( collid == kPPAADATAPeriSIMUL ) return "PP_AA_DATA_PeriL1_SIMUL";

  else return "none";
}

int kEPl2HF = 0;
int kEPOppositeHF = 1;
int kEPSameSideHF = 2;


TString getEPSel( int eventPln) {
  if ( eventPln == kEPl2HF)  return "BothHFs";
  else if ( eventPln == kEPOppositeHF ) return "OppositeHF" ;
  else if ( eventPln == kEPSameSideHF ) return "SameSideHF" ;
  else return "none";
}


int kSoftMuCut = 0;
int kHighPtMuCut = 0;

class DiMuon {
 public:
 DiMuon() :
    run(0), lumi(0), event(0), cBin(0), ep2(0), dphiEp2(0),
    vz(-99),  mass(-1), pt(-1), y(999), phi(999), eta(999),
    pt1(-1), eta1(-1), phi1(-1),        
    pt2(-1), eta2(-1), phi2(-1), weight0(0), weight(0),       
    oniaIndex(-1), softFlag(0), highPtFlag(0),
    qxa(0), qya(0),
    qxb(0), qyb(0),
    qxc(0), qyc(0),
    qxdimu(0), qydimu(0)

    {}
  
  int run;
  int lumi;
  int event;
  int cBin;
  float ep2;
  float dphiEp2;
  float vz;
  float mass;
  float pt;
  float y;
  float phi;    
  float eta;
  float pt1; 
  float eta1;
  float phi1;
  float pt2;
  float eta2;
  float phi2;    
  float weight0;
  float weight;
  int oniaIndex;
  int softFlag;
  int highPtFlag;
  float qxa;
  float qya;
  float qxb;
  float qyb;
  float qxc;
  float qyc;
  float qxdimu;
  float qydimu;

  void clear() {
    run = -99;  lumi=-99; event=-99; cBin=-99; ep2=-99, dphiEp2=-99; 
    vz=-99;     mass = -99; pt=-99; y=-99; phi=-99; eta=-99;      
    pt1=-99; eta1=-99; phi1=-99; pt2=-99; eta2=-99; phi2=-99; weight0=-99, weight=-99;
    oniaIndex=-1; softFlag=-1; highPtFlag=-1; 
    qxa=0;  qya=0; 
    qxb=0;  qyb=0; 
    qxc=0;  qyc=0; 
    qxdimu=0; qydimu=0;
  }

};
TString branchString = "run/I:lumi:event:cBin:ep2/F:dphiEp2:vz:mass:pt:y:phi:eta:pt1:eta1:phi1:pt2:eta2:phi2:weight0:weight:oniaIndex/I:softFlag:highPtFlag:qxa/F:qya:qxb:qyb:qxc:qyc:qxdimu:qydimu";

TString getKineLabel(float ptLow, float ptHigh, float yLow, float yHigh, float muPtCut_, int cLow, int cHigh) {
  TString kineLabel = Form("pt%.1f-%.1f_y%.1f-%.1f_muPt%.1f",ptLow,ptHigh, yLow, yHigh, (float)muPtCut_) ;
    kineLabel = kineLabel+ Form("_centrality%d-%d",(int)cLow, (int)cHigh) ;
  return kineLabel;
}

#endif
