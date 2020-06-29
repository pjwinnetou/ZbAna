#include <ctime>

#include <TLorentzVector.h>
#include "commonUtility.h"
#include "HiEvtPlaneList.h"
#include "cutsAndBinChiB.h"
#include "TreeSetting.h"

static const long MAXTREESIZE = 1000000000000;

void SkimTree(int nevt=-1, bool isMC = false, int kTrigSel = kL1DoubleMuOpen) 
{

  using namespace std;

  //Read Input File
  TString fnameData = "/home/samba.old/UpsilonAnalysis/DataFiles/Onia2018/Conversion/oniaTreeConv_pPb_GTV19_merged.root";
  TString fnameMC = "/eos/cms/store/group/phys_heavyions/dileptons/MC2018/PbPb502TeV/TTrees/Upsi1S_TuneCP5_HydjetDrumMB_officialPythia8MC*_v20190801.root";

  TChain *mytree = new TChain("hionia/myTree");
  if(!isMC){
    mytree->Add(fnameData.Data());
  }
  else if(isMC){
    mytree->Add(fnameMC.Data());
  }

  //SetBranchAddress
  SetTree settree_;
  settree_.TreeSetting(mytree,isMC,1);
//  SetTree::TreeSetting(mytree,isMC,1);
  

  //TriggerIndex
  int trigIndx=0;
  if(kTrigSel == kL1DoubleMuOpen) trigIndx=0;
  else if(kTrigSel == kTrigUps) trigIndx=1;
  else if(kTrigSel == kTrigL1DBOS40100) trigIndx=2;
  else if(kTrigSel == kTrigL1DB50100) trigIndx=3;
  
  //For New Writing File
  TFile* newfile;
  if(isMC) newfile = new TFile(Form("ZBSkim_%sTrig_isMC%d.root",fTrigName[trigIndx].Data(),isMC),"recreate");
  else if(!isMC) newfile = new TFile(Form("ZBSkim_%sTrig_isMC%d.root",fTrigName[trigIndx].Data(),isMC),"recreate");

  //For New Tree
  const static int nMaxDimu = 1000;
  int evt;
  int runN;
  int lumi;
  int cBin;
  int nCand;
  float vz;
  float Dimu_mass[nMaxDimu];
  float Dimu_pt[nMaxDimu];
  float Dimu_y[nMaxDimu];
  float Dimu_phi[nMaxDimu];
  float Dimu_eta[nMaxDimu];
  float Simu_eta1[nMaxDimu];
  float Simu_eta2[nMaxDimu];
  float Simu_phi1[nMaxDimu];
  float Simu_phi2[nMaxDimu];
  float Simu_pt1[nMaxDimu];
  float Simu_pt2[nMaxDimu];
  int recoQQsign[nMaxDimu];
  
  float trk_mass[nMaxDimu];
  float trk_eta[nMaxDimu];
  float trk_y[nMaxDimu];
  float trk_pt[nMaxDimu];

  float Zb_mass[nMaxDimu];
  float Zb_pt[nMaxDimu];
  float Zb_y[nMaxDimu];
  
  float QQtrk_mass[nMaxDimu];
  float QQtrk_pt[nMaxDimu];
  float QQtrk_y[nMaxDimu];
  float Qval[nMaxDimu];

  TTree* mmevttree = new TTree("mmepevt","dimuonAndEventPlanes in event based");
  mmevttree->SetMaxTreeSize(MAXTREESIZE);
  mmevttree->Branch("event",&evt,"event/I");
  mmevttree->Branch("runN",&runN,"runN/I");
  mmevttree->Branch("lumi",&lumi,"lumi/I");
  mmevttree->Branch("vz",&vz,"vz/F");
  mmevttree->Branch("nCand",&nCand,"nCand/I");
  mmevttree->Branch("Dimu_mass",Dimu_mass,"Dimu_mass[nCand]/F");
  mmevttree->Branch("Dimu_y",Dimu_y,"Dimu_y[nCand]/F");
  mmevttree->Branch("Dimu_pt",Dimu_pt,"Dimu_pt[nCand]/F");
  mmevttree->Branch("Dimu_eta",Dimu_eta,"Dimu_eta[nCand]/F");
  mmevttree->Branch("Simu_pt1",Simu_pt1,"Simu_pt1[nCand]/F");
  mmevttree->Branch("Simu_pt2",Simu_pt2,"Simu_pt2[nCand]/F");
  mmevttree->Branch("Simu_eta1",Simu_eta1,"Simu_eta1[nCand]/F");
  mmevttree->Branch("Simu_eta2",Simu_eta2,"Simu_eta2[nCand]/F");
  mmevttree->Branch("recoQQsign",recoQQsign,"recoQQsign[nCand]/I");

  mmevttree->Branch("trk_mass",trk_mass,"trk_mass[nCand]/F");
  mmevttree->Branch("trk_eta",trk_eta,"trk_eta[nCand]/F");
  mmevttree->Branch("trk_pt",trk_pt,"trk_pt[nCand]/F");

  mmevttree->Branch("Zb_mass",Zb_mass,"Zb_mass[nCand]/F");
  mmevttree->Branch("Zb_pt",Zb_pt,"Zb_pt[nCand]/F");
  mmevttree->Branch("Zb_y",Zb_y,"Zb_y[nCand]/F");

  mmevttree->Branch("QQtrk_mass",QQtrk_mass,"QQtrk_mass[nCand]/F");
  mmevttree->Branch("QQtrk_pt",QQtrk_pt,"QQtrk_pt[nCand]/F");
  mmevttree->Branch("QQtrk_y",QQtrk_y,"QQtrk_y[nCand]/F");


      


  ////////////////////////////////////////////////////////////////////////
  ////////////////// TLorentzVector dummies & conversion photon
  ////////////////////////////////////////////////////////////////////////
  TLorentzVector* JP_Reco = new TLorentzVector;
  TLorentzVector* mupl_Reco = new TLorentzVector;
  TLorentzVector* mumi_Reco = new TLorentzVector;
  
  TLorentzVector* trk_4mom = new TLorentzVector;
  TLorentzVector* trkvtx_4mom = new TLorentzVector;
  TLorentzVector* trk_pi_4mom = new TLorentzVector;
  TLorentzVector* Zb_4mom = new TLorentzVector;
  TLorentzVector* QQtrk_4mom = new TLorentzVector;


  // event loop start
  if(nevt == -1) nevt = mytree->GetEntries();
  for(int iev=0; iev<nevt ; ++iev)
  {
    if(iev%100000==0) cout << ">>>>> EVENT " << iev << " / " << mytree->GetEntries() <<  " ("<<(int)(100.*iev/mytree->GetEntries()) << "%)" << endl;

    mytree->GetEntry(iev);

    //Trigger Matching
    if(!( (HLTriggers&((ULong64_t)pow(2, kTrigSel))) == ((ULong64_t)pow(2, kTrigSel)) ) ) continue;

    //Dimuon Loop
    nCand = 0;
    for (Int_t irqq=0; irqq<Reco_QQ_size; ++irqq) 
    {
      //Dimuon Filter Matching
      if(!( (Reco_QQ_trig[irqq]&((ULong64_t)pow(2, kTrigSel))) == ((ULong64_t)pow(2, kTrigSel)) ) ) continue;


      runN = runNb;
      evt = eventNb;
      lumi = LS;
      vz = zVtx;

      JP_Reco = (TLorentzVector*) Reco_QQ_4mom->At(irqq);
      mupl_Reco = (TLorentzVector*) Reco_mu_4mom->At(Reco_QQ_mupl_idx[irqq]);
      mumi_Reco = (TLorentzVector*) Reco_mu_4mom->At(Reco_QQ_mumi_idx[irqq]);
      
   
      //Soft Muon Id cut
      if(!settree_.SoftMuIdCut(irqq)) continue;
     
      //Dimuon Vertex probability cut 
      if ( Reco_QQ_VtxProb[irqq]  < 0.01 ) continue;

      //Acceptance cut --> Based on Jpsi trigger but just for skimming at this stage
      if(! (IsAcceptanceQQ(mupl_Reco->Pt(), mupl_Reco->Eta()) && IsAcceptanceQQ(mumi_Reco->Pt(), mumi_Reco->Eta()) ) ) continue;



      //Loop for Converted photon 

      for(int itrk =0; itrk<Reco_trk_size; itrk++)
      {
        trk_4mom = (TLorentzVector*) Reco_trk_4mom->At(itrk);
        trkvtx_4mom = (TLorentzVector*) Reco_trk_vtx->At(itrk); 
     
        //trk selection cut
        if(trk_4mom->Pt()<=0.5) continue;

        //Set pion 4vector
        trk_pi_4mom ->XYZM(trk_4mom->X(), trk_4mom->Y(), trk_4mom->Z(), pdgMass.PiPlMi);

        //make Z_b candidate and QQtrk Candidate with Q-value cut
        *Zb_4mom = *JP_Reco+*trk_pi_4mom;
        *QQtrk_4mom = *JP_Reco+*trk_4mom;
        Qval = QQtrk_4mom->M()-JP_Reco->M();
        

        recoQQsign[nCand] = Reco_QQ_sign[irqq];     
        Dimu_mass[nCand] = JP_Reco->M();
        Dimu_phi[nCand] = JP_Reco->Phi();
        Simu_phi1[nCand] = mupl_Reco->Phi();
        Simu_phi2[nCand] = mumi_Reco->Phi();
        Dimu_eta[nCand] = JP_Reco->Eta();
        Dimu_y[nCand] = JP_Reco->Rapidity();
        Dimu_pt[nCand] = JP_Reco->Pt();
        Simu_pt1[nCand] = mupl_Reco->Pt();
        Simu_pt2[nCand] = mumi_Reco->Pt();
        Simu_eta1[nCand] = mupl_Reco->Eta();
        Simu_eta2[nCand] = mumi_Reco->Eta();
        
        trk_mass[nCand] = trk_4mom->M(); 
        trk_pt[nCand] = trk_4mom->Pt(); 
        trk_eta[nCand] = trk_4mom->Eta(); 

        Zb_mass[nCand] = Zb_4mom->M();
        Zb_y[nCand] = Zb_4mom->Rapidity();
        Zb_pt[nCand] = Zb_4mom->Pt();
        
        QQtrk_mass[nCand] = QQtrk_4mom->M();
        QQtrk_y[nCand] = QQtrk_4mom->Rapidity();
        QQtrk_pt[nCand] = QQtrk_4mom->Pt();
        nCand++;
        
      } // end of photon loop


    } // end of dimuon loop

    if(nCand>0) mmevttree->Fill();
    
  } //end of event loop

  newfile->cd();
  mmevttree->Write();
  newfile->Close();
  
} 
