// This class has been automatically generated on
// Thu Apr 29 13:15:02 2021 by ROOT version 6.14/04
// from TTree myTree/My TTree of dimuons
// found on file: /data_CMS/cms/diab/JpsiJet/MC/pp/prompt/v3/HiForestAOD_merged.root
//////////////////////////////////////////////////////////

#ifndef myTree_h
#define myTree_h
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include <TLorentzVector.h>
#include <TH1.h>
#include <TH2.h>
#include <cmath>

using namespace std;

// Header file for the classes stored in the TTree if any.
#include "TClonesArray.h"

class myTree {
 public :
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain

  // Fixed size dimensions of array or collections stored in the TTree if any.

  // Declaration of leaf types
  UInt_t          eventNb;
  UInt_t          runNb;
  UInt_t          LS;
  Float_t         zVtx;
  Float_t         nPV;
  Int_t           nTrig;
  Int_t           trigPrescale[99];   //[nTrig]
  ULong64_t       HLTriggers;
  Int_t           Reco_QQ_size;
  Int_t           Reco_QQ_type[99];   //[Reco_QQ_size]
  Int_t           Reco_QQ_sign[99];   //[Reco_QQ_size]
  TClonesArray    *Reco_QQ_4mom;
  Int_t           Reco_QQ_mupl_idx[99];   //[Reco_QQ_size]
  Int_t           Reco_QQ_mumi_idx[99];   //[Reco_QQ_size]
  ULong64_t       Reco_QQ_trig[99];   //[Reco_QQ_size]
  Bool_t          Reco_QQ_isCowboy[99];   //[Reco_QQ_size]
  Float_t         Reco_QQ_ctau[99];   //[Reco_QQ_size]
  Float_t         Reco_QQ_ctauErr[99];   //[Reco_QQ_size]
  Float_t         Reco_QQ_cosAlpha[99];   //[Reco_QQ_size]
  Float_t         Reco_QQ_ctau3D[99];   //[Reco_QQ_size]
  Float_t         Reco_QQ_ctauErr3D[99];   //[Reco_QQ_size]
  Float_t         Reco_QQ_cosAlpha3D[99];   //[Reco_QQ_size]
  Int_t           Reco_QQ_whichGen[99];   //[Reco_QQ_size]
  Float_t         Reco_QQ_VtxProb[99];   //[Reco_QQ_size]
  Float_t         Reco_QQ_dca[99];   //[Reco_QQ_size]
  Float_t         Reco_QQ_MassErr[99];   //[Reco_QQ_size]
  TClonesArray    *Reco_QQ_vtx;
  Int_t           Reco_QQ_Ntrk[99];   //[Reco_QQ_size]
  Float_t         Reco_QQ_mupl_dxy_muonlessVtx[99];   //[Reco_QQ_size]
  Float_t         Reco_QQ_mumi_dxy_muonlessVtx[99];   //[Reco_QQ_size]
  Float_t         Reco_QQ_mupl_dxyErr_muonlessVtx[99];   //[Reco_QQ_size]
  Float_t         Reco_QQ_mumi_dxyErr_muonlessVtx[99];   //[Reco_QQ_size]
  Float_t         Reco_QQ_mupl_dz_muonlessVtx[99];   //[Reco_QQ_size]
  Float_t         Reco_QQ_mumi_dz_muonlessVtx[99];   //[Reco_QQ_size]
  Float_t         Reco_QQ_mupl_dzErr_muonlessVtx[99];   //[Reco_QQ_size]
  Float_t         Reco_QQ_mumi_dzErr_muonlessVtx[99];   //[Reco_QQ_size]
  Int_t           Reco_mu_size;
  Int_t           Reco_mu_type[99];   //[Reco_mu_size]
  Int_t           Reco_mu_whichGen[99];   //[Reco_mu_size]
  Int_t           Reco_mu_SelectionType[99];   //[Reco_mu_size]
  Int_t           Reco_mu_charge[99];   //[Reco_mu_size]
  TClonesArray    *Reco_mu_4mom;
  ULong64_t       Reco_mu_trig[99];   //[Reco_mu_size]
  Bool_t          Reco_mu_highPurity[99];   //[Reco_mu_size]
  Bool_t          Reco_mu_TrkMuArb[99];   //[Reco_mu_size]
  Bool_t          Reco_mu_TMOneStaTight[99];   //[Reco_mu_size]
  Int_t           Reco_mu_nPixValHits[99];   //[Reco_mu_size]
  Int_t           Reco_mu_nMuValHits[99];   //[Reco_mu_size]
  Int_t           Reco_mu_nTrkHits[99];   //[Reco_mu_size]
  Float_t         Reco_mu_normChi2_inner[99];   //[Reco_mu_size]
  Float_t         Reco_mu_normChi2_global[99];   //[Reco_mu_size]
  Int_t           Reco_mu_nPixWMea[99];   //[Reco_mu_size]
  Int_t           Reco_mu_nTrkWMea[99];   //[Reco_mu_size]
  Int_t           Reco_mu_StationsMatched[99];   //[Reco_mu_size]
  Float_t         Reco_mu_dxy[99];   //[Reco_mu_size]
  Float_t         Reco_mu_dxyErr[99];   //[Reco_mu_size]
  Float_t         Reco_mu_dz[99];   //[Reco_mu_size]
  Float_t         Reco_mu_dzErr[99];   //[Reco_mu_size]
  Float_t         Reco_mu_pt_inner[99];   //[Reco_mu_size]
  Float_t         Reco_mu_pt_global[99];   //[Reco_mu_size]
  Float_t         Reco_mu_ptErr_inner[99];   //[Reco_mu_size]
  Float_t         Reco_mu_ptErr_global[99];   //[Reco_mu_size]
  Int_t           Gen_QQ_size;
  Int_t           Gen_QQ_type[99];   //[Gen_QQ_size]
  TClonesArray    *Gen_QQ_4mom;
  Int_t           Gen_QQ_momId[99];   //[Gen_QQ_size]
  Float_t         Gen_QQ_ctau[99];   //[Gen_QQ_size]
  Float_t         Gen_QQ_ctau3D[99];   //[Gen_QQ_size]
  Int_t           Gen_QQ_mupl_idx[99];   //[Gen_QQ_size]
  Int_t           Gen_QQ_mumi_idx[99];   //[Gen_QQ_size]
  Int_t           Gen_QQ_whichRec[99];   //[Gen_QQ_size]
  Int_t           Gen_mu_size;
  Int_t           Gen_mu_type[99];   //[Gen_mu_size]
  Int_t           Gen_mu_charge[99];   //[Gen_mu_size]
  TClonesArray    *Gen_mu_4mom;
  Int_t           Gen_mu_whichRec[99];   //[Gen_mu_size]
  Float_t         pthat;
  Bool_t          isPbPb;
  Bool_t          isAcc;
  //centrality
  Int_t           hiBin;
  Float_t         hiHF;
  Int_t           Centrality;
  

  // List of branches
  TBranch        *b_eventNb;   //!
  TBranch        *b_runNb;   //!
  TBranch        *b_LS;   //!
  TBranch        *b_zVtx;   //!
  TBranch        *b_nPV;   //!
  TBranch        *b_nTrig;   //!
  TBranch        *b_trigPrescale;   //!
  TBranch        *b_HLTriggers;   //!
  TBranch        *b_Reco_QQ_size;   //!
  TBranch        *b_Reco_QQ_type;   //!
  TBranch        *b_Reco_QQ_sign;   //!
  TBranch        *b_Reco_QQ_4mom;   //!
  TBranch        *b_Reco_QQ_mupl_idx;   //!
  TBranch        *b_Reco_QQ_mumi_idx;   //!
  TBranch        *b_Reco_QQ_trig;   //!
  TBranch        *b_Reco_QQ_isCowboy;   //!
  TBranch        *b_Reco_QQ_ctau;   //!
  TBranch        *b_Reco_QQ_ctauErr;   //!
  TBranch        *b_Reco_QQ_cosAlpha;   //!
  TBranch        *b_Reco_QQ_ctau3D;   //!
  TBranch        *b_Reco_QQ_ctauErr3D;   //!
  TBranch        *b_Reco_QQ_cosAlpha3D;   //!
  TBranch        *b_Reco_QQ_whichGen;   //!
  TBranch        *b_Reco_QQ_VtxProb;   //!
  TBranch        *b_Reco_QQ_dca;   //!
  TBranch        *b_Reco_QQ_MassErr;   //!
  TBranch        *b_Reco_QQ_vtx;   //!
  TBranch        *b_Reco_QQ_Ntrk;   //!
  TBranch        *b_Reco_QQ_mupl_dxy_muonlessVtx;   //!
  TBranch        *b_Reco_QQ_mumi_dxy_muonlessVtx;   //!
  TBranch        *b_Reco_QQ_mupl_dxyErr_muonlessVtx;   //!
  TBranch        *b_Reco_QQ_mumi_dxyErr_muonlessVtx;   //!
  TBranch        *b_Reco_QQ_mupl_dz_muonlessVtx;   //!
  TBranch        *b_Reco_QQ_mumi_dz_muonlessVtx;   //!
  TBranch        *b_Reco_QQ_mupl_dzErr_muonlessVtx;   //!
  TBranch        *b_Reco_QQ_mumi_dzErr_muonlessVtx;   //!
  TBranch        *b_Reco_mu_size;   //!
  TBranch        *b_Reco_mu_type;   //!
  TBranch        *b_Reco_mu_whichGen;   //!
  TBranch        *b_Reco_mu_SelectionType;   //!
  TBranch        *b_Reco_mu_charge;   //!
  TBranch        *b_Reco_mu_4mom;   //!
  TBranch        *b_Reco_mu_trig;   //!
  TBranch        *b_Reco_mu_highPurity;   //!
  TBranch        *b_Reco_mu_TrkMuArb;   //!
  TBranch        *b_Reco_mu_TMOneStaTight;   //!
  TBranch        *b_Reco_mu_nPixValHits;   //!
  TBranch        *b_Reco_mu_nMuValHits;   //!
  TBranch        *b_Reco_mu_nTrkHits;   //!
  TBranch        *b_Reco_mu_normChi2_inner;   //!
  TBranch        *b_Reco_mu_normChi2_global;   //!
  TBranch        *b_Reco_mu_nPixWMea;   //!
  TBranch        *b_Reco_mu_nTrkWMea;   //!
  TBranch        *b_Reco_mu_StationsMatched;   //!
  TBranch        *b_Reco_mu_dxy;   //!
  TBranch        *b_Reco_mu_dxyErr;   //!
  TBranch        *b_Reco_mu_dz;   //!
  TBranch        *b_Reco_mu_dzErr;   //!
  TBranch        *b_Reco_mu_pt_inner;   //!
  TBranch        *b_Reco_mu_pt_global;   //!
  TBranch        *b_Reco_mu_ptErr_inner;   //!
  TBranch        *b_Reco_mu_ptErr_global;   //!
  TBranch        *b_Gen_QQ_size;   //!
  TBranch        *b_Gen_QQ_type;   //!
  TBranch        *b_Gen_QQ_4mom;   //!
  TBranch        *b_Gen_QQ_momId;   //!
  TBranch        *b_Gen_QQ_ctau;   //!
  TBranch        *b_Gen_QQ_ctau3D;   //!
  TBranch        *b_Gen_QQ_mupl_idx;   //!
  TBranch        *b_Gen_QQ_mumi_idx;   //!
  TBranch        *b_Gen_QQ_whichRec;   //!
  TBranch        *b_Gen_mu_size;   //!
  TBranch        *b_Gen_mu_type;   //!
  TBranch        *b_Gen_mu_charge;   //!
  TBranch        *b_Gen_mu_4mom;   //!
  TBranch        *b_Gen_mu_whichRec;   //!
  TBranch        *b_hiBin;
  TBranch        *b_hiHF;
  TBranch        *b_pthat;
  TBranch        *b_Centrality;
  myTree(Bool_t pbpb, Bool_t acc);
  virtual ~myTree();
  virtual Int_t      Cut(Long64_t entry);
  virtual Int_t      GetEntry(Long64_t entry);
  virtual Long64_t   LoadTree(Long64_t entry);
  virtual void       Init(TTree *tree);
  virtual Bool_t     Notify();
  virtual void       Show(Long64_t entry = -1);
  virtual Bool_t     isMuonInAccept2019 (TLorentzVector* Muon);
  virtual Bool_t     passQualityCuts2019 (Int_t iRecoMu);
  virtual Bool_t     isTriggerMatch (Int_t iRecoMu, Int_t TriggerBit);
  virtual TObjArray* eff (string effType,string Pbp,TH2D *num, TH2D *denom, int j);
  virtual void       AccCalc();
  virtual void       EffCalc();
  virtual Int_t      getMCHiBinFromhiHF(const Double_t hiHF);
  virtual Double_t   findNcoll(int hiBin);
};
#endif

#ifdef myTree_cxx
myTree::myTree(Bool_t pbpb, Bool_t acc) : fChain(0)
{
  TFile* f(0x0);
  isPbPb = pbpb;
  isAcc = acc;
  
  TString inputFiles [4] = {
    //Acc files
    //pp prompt 
    "/data_CMS/cms/mnguyen/jPsiJet/mc/prompt/acc/merged_acc.root",
    //PbPb prompt
    "/data_CMS/cms/mnguyen/jPsiJet/mc/prompt/acc/merged_acc.root",
    //Eff files
    //pp prompt
    "/data_CMS/cms/diab/JpsiJet/MC/pp/prompt/v3/HiForestAOD_merged.root", 
    //PbPb prompt
    //"/data_CMS/cms/diab/JpsiJet/MC/PbPb/prompt/v7/HiForestAOD_merged.root",
    "/data_CMS/cms/diab/tnpOniaTree/OniaTree_HINPbPbAutumn18DR-mva98_JPsi_MC_0.root", 
  };
  
  int inin = 0; // inin = input index
  
  if (isAcc && !isPbPb) inin = 0;
  else if (isAcc && isPbPb) inin = 1;
  else if (!isAcc && !isPbPb) inin = 2;
  else if (!isAcc && isPbPb) inin = 3;
  
  f = TFile::Open(inputFiles[inin]);
  
  cout<<Form("[INFO] %s tree in file ",(isPbPb?"PbPb":"pp"))<<f->GetName()<<endl;
  
  TTree * tree (0x0); TTree * skimtree (0x0); TTree * centtree (0x0);
  tree = (TTree*) f->Get("hionia/myTree");
  skimtree = (TTree*) f->Get("skimanalysis/HltTree");
  centtree = (TTree*) f->Get("hiEvtAnalyzer/HiTree");
  if(!tree){
    cout <<"[ERROR] Cannot find the onia tree"<<endl;
    return;
  }
  if(!skimtree)
    cout <<"[WARNING] Cannot find the skimanalysis tree"<<endl;
  else
    tree->AddFriend(skimtree);
  if (!centtree)
    cout<<"[WARNING] Cannot find the hiEvtAnalyzer tree"<<endl;
  else
    tree->AddFriend(centtree);
  Init(tree);
}

myTree::~myTree()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

Int_t myTree::GetEntry(Long64_t entry)
{
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}
Long64_t myTree::LoadTree(Long64_t entry)
{
  // Set the environment to read one entry
  if (!fChain) return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) return centry;
  if (fChain->GetTreeNumber() != fCurrent) {
    fCurrent = fChain->GetTreeNumber();
    Notify();
  }
  return centry;
}

void myTree::Init(TTree *tree)
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses and branch
  // pointers of the tree will be set.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).
  
  // Set object pointer  Reco_QQ_4mom = 0;
  Reco_QQ_vtx = 0;
  Reco_mu_4mom = 0;
  Gen_QQ_4mom = 0;
  Gen_mu_4mom = 0;
  // Set branch addresses and branch pointers
  if (!tree) return;
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);
  
  
  if (fChain->GetBranch("Centrality")) fChain->SetBranchAddress("Centrality", &Centrality, &b_Centrality);
  if (fChain->GetBranch("eventNb")) fChain->SetBranchAddress("eventNb", &eventNb, &b_eventNb);
  if (fChain->GetBranch("runNb")) fChain->SetBranchAddress("runNb", &runNb, &b_runNb);
  if (fChain->GetBranch("LS")) fChain->SetBranchAddress("LS", &LS, &b_LS);
  if (fChain->GetBranch("zVtx")) fChain->SetBranchAddress("zVtx", &zVtx, &b_zVtx);
  if (fChain->GetBranch("nPV")) fChain->SetBranchAddress("nPV", &nPV, &b_nPV);
  if (fChain->GetBranch("nTrig")) fChain->SetBranchAddress("nTrig", &nTrig, &b_nTrig);
  if (fChain->GetBranch("trigPrescale")) fChain->SetBranchAddress("trigPrescale", trigPrescale, &b_trigPrescale);
  if (fChain->GetBranch("HLTriggers")) fChain->SetBranchAddress("HLTriggers", &HLTriggers, &b_HLTriggers);
  if (fChain->GetBranch("Reco_QQ_size")) fChain->SetBranchAddress("Reco_QQ_size", &Reco_QQ_size, &b_Reco_QQ_size);
  if (fChain->GetBranch("Reco_QQ_type")) fChain->SetBranchAddress("Reco_QQ_type", Reco_QQ_type, &b_Reco_QQ_type);
  if (fChain->GetBranch("Reco_QQ_sign")) fChain->SetBranchAddress("Reco_QQ_sign", Reco_QQ_sign, &b_Reco_QQ_sign);
  if (fChain->GetBranch("Reco_QQ_4mom")) fChain->SetBranchAddress("Reco_QQ_4mom", &Reco_QQ_4mom, &b_Reco_QQ_4mom);
  if (fChain->GetBranch("Reco_QQ_mupl_idx")) fChain->SetBranchAddress("Reco_QQ_mupl_idx", Reco_QQ_mupl_idx, &b_Reco_QQ_mupl_idx);
  if (fChain->GetBranch("Reco_QQ_mumi_idx")) fChain->SetBranchAddress("Reco_QQ_mumi_idx", Reco_QQ_mumi_idx, &b_Reco_QQ_mumi_idx);
  if (fChain->GetBranch("Reco_QQ_trig")) fChain->SetBranchAddress("Reco_QQ_trig", Reco_QQ_trig, &b_Reco_QQ_trig);
  if (fChain->GetBranch("Reco_QQ_isCowboy")) fChain->SetBranchAddress("Reco_QQ_isCowboy", Reco_QQ_isCowboy, &b_Reco_QQ_isCowboy);
  if (fChain->GetBranch("Reco_QQ_ctau")) fChain->SetBranchAddress("Reco_QQ_ctau", Reco_QQ_ctau, &b_Reco_QQ_ctau);
  if (fChain->GetBranch("Reco_QQ_ctauErr")) fChain->SetBranchAddress("Reco_QQ_ctauErr", Reco_QQ_ctauErr, &b_Reco_QQ_ctauErr);
  if (fChain->GetBranch("Reco_QQ_cosAlpha")) fChain->SetBranchAddress("Reco_QQ_cosAlpha", Reco_QQ_cosAlpha, &b_Reco_QQ_cosAlpha);
  if (fChain->GetBranch("Reco_QQ_ctau3D")) fChain->SetBranchAddress("Reco_QQ_ctau3D", Reco_QQ_ctau3D, &b_Reco_QQ_ctau3D);
  if (fChain->GetBranch("Reco_QQ_ctauErr3D")) fChain->SetBranchAddress("Reco_QQ_ctauErr3D", Reco_QQ_ctauErr3D, &b_Reco_QQ_ctauErr3D);
  if (fChain->GetBranch("Reco_QQ_cosAlpha3D")) fChain->SetBranchAddress("Reco_QQ_cosAlpha3D", Reco_QQ_cosAlpha3D, &b_Reco_QQ_cosAlpha3D);
  if (fChain->GetBranch("Reco_QQ_whichGen")) fChain->SetBranchAddress("Reco_QQ_whichGen", Reco_QQ_whichGen, &b_Reco_QQ_whichGen);
  if (fChain->GetBranch("Reco_QQ_VtxProb")) fChain->SetBranchAddress("Reco_QQ_VtxProb", Reco_QQ_VtxProb, &b_Reco_QQ_VtxProb);
  if (fChain->GetBranch("Reco_QQ_dca")) fChain->SetBranchAddress("Reco_QQ_dca", Reco_QQ_dca, &b_Reco_QQ_dca);
  if (fChain->GetBranch("Reco_QQ_MassErr")) fChain->SetBranchAddress("Reco_QQ_MassErr", Reco_QQ_MassErr, &b_Reco_QQ_MassErr);
  if (fChain->GetBranch("Reco_QQ_vtx")) fChain->SetBranchAddress("Reco_QQ_vtx", &Reco_QQ_vtx, &b_Reco_QQ_vtx);
  if (fChain->GetBranch("Reco_QQ_Ntrk")) fChain->SetBranchAddress("Reco_QQ_Ntrk", Reco_QQ_Ntrk, &b_Reco_QQ_Ntrk);
  if (fChain->GetBranch("Reco_QQ_mupl_dxy_muonlessVtx")) fChain->SetBranchAddress("Reco_QQ_mupl_dxy_muonlessVtx", Reco_QQ_mupl_dxy_muonlessVtx
, &b_Reco_QQ_mupl_dxy_muonlessVtx);
  if (fChain->GetBranch("Reco_QQ_mumi_dxy_muonlessVtx")) fChain->SetBranchAddress("Reco_QQ_mumi_dxy_muonlessVtx", Reco_QQ_mumi_dxy_muonlessVtx
, &b_Reco_QQ_mumi_dxy_muonlessVtx);
  if (fChain->GetBranch("Reco_QQ_mupl_dxyErr_muonlessVtx")) fChain->SetBranchAddress("Reco_QQ_mupl_dxyErr_muonlessVtx", Reco_QQ_mupl_dxyErr_mu
onlessVtx, &b_Reco_QQ_mupl_dxyErr_muonlessVtx);
  if (fChain->GetBranch("Reco_QQ_mumi_dxyErr_muonlessVtx")) fChain->SetBranchAddress("Reco_QQ_mumi_dxyErr_muonlessVtx", Reco_QQ_mumi_dxyErr_mu
onlessVtx, &b_Reco_QQ_mumi_dxyErr_muonlessVtx);
  if (fChain->GetBranch("Reco_QQ_mupl_dz_muonlessVtx")) fChain->SetBranchAddress("Reco_QQ_mupl_dz_muonlessVtx", Reco_QQ_mupl_dz_muonlessVtx, &
b_Reco_QQ_mupl_dz_muonlessVtx);
  if (fChain->GetBranch("Reco_QQ_mumi_dz_muonlessVtx")) fChain->SetBranchAddress("Reco_QQ_mumi_dz_muonlessVtx", Reco_QQ_mumi_dz_muonlessVtx, &
b_Reco_QQ_mumi_dz_muonlessVtx);
  if (fChain->GetBranch("Reco_QQ_mupl_dzErr_muonlessVtx")) fChain->SetBranchAddress("Reco_QQ_mupl_dzErr_muonlessVtx", Reco_QQ_mupl_dzErr_muonl
essVtx, &b_Reco_QQ_mupl_dzErr_muonlessVtx);
  if (fChain->GetBranch("Reco_QQ_mumi_dzErr_muonlessVtx")) fChain->SetBranchAddress("Reco_QQ_mumi_dzErr_muonlessVtx", Reco_QQ_mumi_dzErr_muonl
essVtx, &b_Reco_QQ_mumi_dzErr_muonlessVtx);
  if (fChain->GetBranch("Reco_mu_size")) fChain->SetBranchAddress("Reco_mu_size", &Reco_mu_size, &b_Reco_mu_size);
  if (fChain->GetBranch("Reco_mu_type")) fChain->SetBranchAddress("Reco_mu_type", Reco_mu_type, &b_Reco_mu_type);
  if (fChain->GetBranch("Reco_mu_whichGen")) fChain->SetBranchAddress("Reco_mu_whichGen", Reco_mu_whichGen, &b_Reco_mu_whichGen);
  if (fChain->GetBranch("Reco_mu_SelectionType")) fChain->SetBranchAddress("Reco_mu_SelectionType", Reco_mu_SelectionType, &b_Reco_mu_Selectio
nType);
  if (fChain->GetBranch("Reco_mu_charge")) fChain->SetBranchAddress("Reco_mu_charge", Reco_mu_charge, &b_Reco_mu_charge);
  if (fChain->GetBranch("Reco_mu_4mom")) fChain->SetBranchAddress("Reco_mu_4mom", &Reco_mu_4mom, &b_Reco_mu_4mom);
  if (fChain->GetBranch("Reco_mu_trig")) fChain->SetBranchAddress("Reco_mu_trig", Reco_mu_trig, &b_Reco_mu_trig);
  if (fChain->GetBranch("Reco_mu_highPurity")) fChain->SetBranchAddress("Reco_mu_highPurity", Reco_mu_highPurity, &b_Reco_mu_highPurity);
  if (fChain->GetBranch("Reco_mu_TrkMuArb")) fChain->SetBranchAddress("Reco_mu_TrkMuArb", Reco_mu_TrkMuArb, &b_Reco_mu_TrkMuArb);
  if (fChain->GetBranch("Reco_mu_TMOneStaTight")) fChain->SetBranchAddress("Reco_mu_TMOneStaTight", Reco_mu_TMOneStaTight, &b_Reco_mu_TMOneSta
Tight);
  if (fChain->GetBranch("Reco_mu_nPixValHits")) fChain->SetBranchAddress("Reco_mu_nPixValHits", Reco_mu_nPixValHits, &b_Reco_mu_nPixValHits);
  if (fChain->GetBranch("Reco_mu_nMuValHits")) fChain->SetBranchAddress("Reco_mu_nMuValHits", Reco_mu_nMuValHits, &b_Reco_mu_nMuValHits);
  if (fChain->GetBranch("Reco_mu_nTrkHits")) fChain->SetBranchAddress("Reco_mu_nTrkHits", Reco_mu_nTrkHits, &b_Reco_mu_nTrkHits);
  if (fChain->GetBranch("Reco_mu_normChi2_inner")) fChain->SetBranchAddress("Reco_mu_normChi2_inner", Reco_mu_normChi2_inner, &b_Reco_mu_normC
hi2_inner);
  if (fChain->GetBranch("Reco_mu_normChi2_global")) fChain->SetBranchAddress("Reco_mu_normChi2_global", Reco_mu_normChi2_global, &b_Reco_mu_no
rmChi2_global);
  if (fChain->GetBranch("Reco_mu_nPixWMea")) fChain->SetBranchAddress("Reco_mu_nPixWMea", Reco_mu_nPixWMea, &b_Reco_mu_nPixWMea);
  if (fChain->GetBranch("Reco_mu_nTrkWMea")) fChain->SetBranchAddress("Reco_mu_nTrkWMea", Reco_mu_nTrkWMea, &b_Reco_mu_nTrkWMea);
  if (fChain->GetBranch("Reco_mu_StationsMatched")) fChain->SetBranchAddress("Reco_mu_StationsMatched", Reco_mu_StationsMatched, &b_Reco_mu_St
ationsMatched);
  if (fChain->GetBranch("Reco_mu_dxy")) fChain->SetBranchAddress("Reco_mu_dxy", Reco_mu_dxy, &b_Reco_mu_dxy);
  if (fChain->GetBranch("Reco_mu_dxyErr")) fChain->SetBranchAddress("Reco_mu_dxyErr", Reco_mu_dxyErr, &b_Reco_mu_dxyErr);
  if (fChain->GetBranch("Reco_,u_dz")) fChain->SetBranchAddress("Reco_mu_dz", Reco_mu_dz, &b_Reco_mu_dz);
  if (fChain->GetBranch("Reco_mu_dzErr")) fChain->SetBranchAddress("Reco_mu_dzErr", Reco_mu_dzErr, &b_Reco_mu_dzErr);
  if (fChain->GetBranch("Reco_mu_pt_inner")) fChain->SetBranchAddress("Reco_mu_pt_inner", Reco_mu_pt_inner, &b_Reco_mu_pt_inner);
  if (fChain->GetBranch("Reco_mu_pt_global")) fChain->SetBranchAddress("Reco_mu_pt_global", Reco_mu_pt_global, &b_Reco_mu_pt_global);
  if (fChain->GetBranch("Reco_mu_ptErr_inner")) fChain->SetBranchAddress("Reco_mu_ptErr_inner", Reco_mu_ptErr_inner, &b_Reco_mu_ptErr_inner);
  if (fChain->GetBranch("Reco_mu_ptErr_global")) fChain->SetBranchAddress("Reco_mu_ptErr_global", Reco_mu_ptErr_global, &b_Reco_mu_ptErr_globa
l);
  if (fChain->GetBranch("Gen_QQ_size")) fChain->SetBranchAddress("Gen_QQ_size", &Gen_QQ_size, &b_Gen_QQ_size);
  if (fChain->GetBranch("Gen_QQ_type")) fChain->SetBranchAddress("Gen_QQ_type", Gen_QQ_type, &b_Gen_QQ_type);
  if (fChain->GetBranch("Gen_QQ_4mom")) fChain->SetBranchAddress("Gen_QQ_4mom", &Gen_QQ_4mom, &b_Gen_QQ_4mom);
  if (fChain->GetBranch("Gen_QQ_momId")) fChain->SetBranchAddress("Gen_QQ_momId", Gen_QQ_momId, &b_Gen_QQ_momId);
  if (fChain->GetBranch("Gen_QQ_ctau")) fChain->SetBranchAddress("Gen_QQ_ctau", Gen_QQ_ctau, &b_Gen_QQ_ctau);
  if (fChain->GetBranch("Gen_QQ_ctau3D")) fChain->SetBranchAddress("Gen_QQ_ctau3D", Gen_QQ_ctau3D, &b_Gen_QQ_ctau3D);
  if (fChain->GetBranch("Gen_QQ_mupl_idx")) fChain->SetBranchAddress("Gen_QQ_mupl_idx", Gen_QQ_mupl_idx, &b_Gen_QQ_mupl_idx);
  if (fChain->GetBranch("Gen_QQ_mumi_idx")) fChain->SetBranchAddress("Gen_QQ_mumi_idx", Gen_QQ_mumi_idx, &b_Gen_QQ_mumi_idx);
  if (fChain->GetBranch("Gen_QQ_whichRec")) fChain->SetBranchAddress("Gen_QQ_whichRec", Gen_QQ_whichRec, &b_Gen_QQ_whichRec);
  if (fChain->GetBranch("Gen_mu_size")) fChain->SetBranchAddress("Gen_mu_size", &Gen_mu_size, &b_Gen_mu_size);
  if (fChain->GetBranch("Gen_mu_type")) fChain->SetBranchAddress("Gen_mu_type", Gen_mu_type, &b_Gen_mu_type);
  if (fChain->GetBranch("Gen_mu_charge")) fChain->SetBranchAddress("Gen_mu_charge", Gen_mu_charge, &b_Gen_mu_charge);
  if (fChain->GetBranch("Gen_mu_4mom")) fChain->SetBranchAddress("Gen_mu_4mom", &Gen_mu_4mom, &b_Gen_mu_4mom);
  if (fChain->GetBranch("Gen_mu_whichRec")) fChain->SetBranchAddress("Gen_mu_whichRec", Gen_mu_whichRec, &b_Gen_mu_whichRec);
  if (fChain->GetBranch("hiBin")) fChain->SetBranchAddress("hiBin", &hiBin, &b_hiBin);
  if (fChain->GetBranch("hiHF")) fChain->SetBranchAddress("hiHF", &hiHF, &b_hiHF);
  if (fChain->GetBranch("pthat")) fChain->SetBranchAddress("pthat", &pthat, &b_pthat);
  Notify();
}

Bool_t myTree::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.

  return kTRUE;
}

void myTree::Show(Long64_t entry)
{
  // Print contents of entry.
  // If entry is not specified, print current entry
  if (!fChain) return;
  fChain->Show(entry);
}
Int_t myTree::Cut(Long64_t entry)
{
  // This function may be called from Loop.
  // returns  1 if entry is accepted.
  // returns -1 otherwise.
  return 1;
}
Bool_t myTree::isMuonInAccept2019 (TLorentzVector* Muon)
{
  return ( fabs(Muon->Eta()) < 2.4 &&
	   ((fabs(Muon->Eta()) < 1.2 && Muon->Pt() >= 3.5) ||
	    (1.2 <= fabs(Muon->Eta()) && fabs(Muon->Eta()) < 2.1 && Muon->Pt() >= 5.47-1.89*fabs(Muon->Eta())) ||
	    (2.1 <= fabs(Muon->Eta()) && Muon->Pt() >= 1.5)));
};
Bool_t myTree::passQualityCuts2019 (Int_t iRecoMu)
{
  
  if ( ! (Reco_mu_SelectionType[iRecoMu]&((ULong64_t)pow(2, 1))) ) return false; // require the muons to be global muons
  // cout<<"1"<<endl;
  if ( ! (Reco_mu_SelectionType[iRecoMu]&((ULong64_t)pow(2, 3))) ) return false; // require the muons to be tracker muons
  //cout<<"2"<<endl;
  if ( ! (Reco_mu_nTrkWMea[iRecoMu] > 5) ) return false;
  //cout<<"3"<<endl;
  if ( ! (Reco_mu_nPixWMea[iRecoMu] > 0) ) return false;
  //cout<<"4"<<endl;
  if ( ! (fabs(Reco_mu_dxy[iRecoMu]) < 0.3) ) return false;
  if ( ! (fabs(Reco_mu_dz[iRecoMu]) < 20.) ) return false;
  
  return true;  
};

Bool_t myTree::isTriggerMatch (Int_t iRecoMu, Int_t TriggerBit)
{
  Bool_t cond = true;
  cond = cond && ( (HLTriggers&((ULong64_t)pow(2, TriggerBit))) == ((ULong64_t)pow(2, TriggerBit)) );
  cond = cond && ( (Reco_mu_trig[iRecoMu]&((ULong64_t)pow(2, TriggerBit))) == ((ULong64_t)pow(2, TriggerBit)) );
  return cond;
};

Int_t myTree::getMCHiBinFromhiHF(const Double_t hiHF) {
  const Int_t nBins = 200; // table of bin edges
  const Double_t binTable[nBins+1] = {0, 12.2187, 13.0371, 13.7674, 14.5129, 15.2603, 16.0086, 16.7623, 17.5335, 18.3283, 19.1596, 19.9989, 20
.8532, 21.7297, 22.6773, 23.6313, 24.6208, 25.6155, 26.6585, 27.7223, 28.8632, 30.041, 31.2865, 32.5431, 33.8655, 35.2539, 36.6912, 38.2064, 3
9.7876, 41.4818, 43.2416, 45.0605, 46.9652, 48.9918, 51.1, 53.2417, 55.5094, 57.9209, 60.3817, 62.9778, 65.6099, 68.4352, 71.3543, 74.4154, 77
.6252, 80.8425, 84.1611, 87.7395, 91.3973, 95.1286, 99.0571, 103.185, 107.482, 111.929, 116.45, 121.178, 126.081, 130.995, 136.171, 141.612, 1
47.298, 153.139, 159.419, 165.633, 172.114, 178.881, 185.844, 192.845, 200.244, 207.83, 215.529, 223.489, 231.878, 240.254, 249.319, 258.303, 
267.508, 277.037, 286.729, 296.845, 307.458, 317.882, 328.787, 340.074, 351.295, 362.979, 375.125, 387.197, 399.604, 412.516, 425.683, 439.001
, 452.667, 466.816, 481.007, 495.679, 510.588, 526.138, 541.782, 557.641, 574.141, 591.071, 608.379, 626.068, 643.616, 661.885, 680.288, 699.4
49, 718.925, 738.968, 758.983, 779.459, 800.376, 821.638, 843.555, 865.771, 888.339, 911.031, 934.979, 958.56, 982.582, 1007.02, 1031.9, 1057.
81, 1084.01, 1111.71, 1138.21, 1165.72, 1193.73, 1221.65, 1251.51, 1281.23, 1311.01, 1341.1, 1372.4, 1404.29, 1436.52, 1468.65, 1501.91, 1535.
56, 1569.69, 1604.69, 1640.65, 1676.05, 1712.62, 1749.28, 1787.43, 1825.89, 1866.07, 1906.58, 1947.84, 1989.66, 2031.4, 2072.8, 2115.32, 2159.
5, 2205.23, 2252.68, 2298.58, 2345.65, 2393.36, 2442.87, 2491.45, 2541.04, 2592.81, 2645.52, 2699.1, 2753.29, 2807.93, 2864.37, 2922.6, 2979.4
2, 3038.68, 3098.72, 3159.29, 3221.66, 3285.9, 3350.95, 3415.81, 3482.69, 3552.62, 3623.61, 3694.63, 3767.25, 3840.28, 3917.04, 3993.66, 4073.
36, 4154.33, 4238.13, 4322.21, 4409.83, 4498.89, 4589.72, 4681.56, 4777.09, 4877.95, 4987.05, 5113.04, 5279.58, 6242.82};
  
  Int_t binPos = -1;
  for(int i = 0; i < nBins; ++i){
    if(hiHF >= binTable[i] && hiHF < binTable[i+1]){
      binPos = i;
      break;
    }
  }
  binPos = nBins - 1 - binPos;
  return (Int_t)(200*((Double_t)binPos)/((Double_t)nBins));
}

Double_t myTree::findNcoll(int hiBin) {
  const int nbins = 200;
  const Double_t Ncoll[nbins] = {1976.95, 1944.02, 1927.29, 1891.9, 1845.3, 1807.2, 1760.45, 1729.18, 1674.8, 1630.3, 1590.52, 1561.72, 1516.1
, 1486.5, 1444.68, 1410.88, 1376.4, 1347.32, 1309.71, 1279.98, 1255.31, 1219.89, 1195.13, 1165.96, 1138.92, 1113.37, 1082.26, 1062.42, 1030.6,
 1009.96, 980.229, 955.443, 936.501, 915.97, 892.063, 871.289, 847.364, 825.127, 806.584, 789.163, 765.42, 751.187, 733.001, 708.31, 690.972, 
677.711, 660.682, 640.431, 623.839, 607.456, 593.307, 576.364, 560.967, 548.909, 530.475, 519.575, 505.105, 490.027, 478.133, 462.372, 451.115
, 442.642, 425.76, 416.364, 405.154, 392.688, 380.565, 371.167, 360.28, 348.239, 340.587, 328.746, 320.268, 311.752, 300.742, 292.172, 281.361
, 274.249, 267.025, 258.625, 249.931, 240.497, 235.423, 228.63, 219.854, 214.004, 205.425, 199.114, 193.618, 185.644, 180.923, 174.289, 169.64
1, 161.016, 157.398, 152.151, 147.425, 140.933, 135.924, 132.365, 127.017, 122.127, 117.817, 113.076, 109.055, 105.16, 101.323, 98.098, 95.054
8, 90.729, 87.6495, 84.0899, 80.2237, 77.2201, 74.8848, 71.3554, 68.7745, 65.9911, 63.4136, 61.3859, 58.1903, 56.4155, 53.8486, 52.0196, 49.29
21, 47.0735, 45.4345, 43.8434, 41.7181, 39.8988, 38.2262, 36.4435, 34.8984, 33.4664, 31.8056, 30.351, 29.2074, 27.6924, 26.7754, 25.4965, 24.2
802, 22.9651, 22.0059, 21.0915, 19.9129, 19.1041, 18.1487, 17.3218, 16.5957, 15.5323, 14.8035, 14.2514, 13.3782, 12.8667, 12.2891, 11.61, 11.0
026, 10.3747, 9.90294, 9.42648, 8.85324, 8.50121, 7.89834, 7.65197, 7.22768, 6.7755, 6.34855, 5.98336, 5.76555, 5.38056, 5.11024, 4.7748, 4.59
117, 4.23247, 4.00814, 3.79607, 3.68702, 3.3767, 3.16309, 2.98282, 2.8095, 2.65875, 2.50561, 2.32516, 2.16357, 2.03235, 1.84061, 1.72628, 1.62
305, 1.48916, 1.38784, 1.28366, 1.24693, 1.18552, 1.16085, 1.12596, 1.09298, 1.07402, 1.06105, 1.02954};
  return Ncoll[hiBin];
};

#endif // #ifdef myTree_cxx
-bash-4.2$   

