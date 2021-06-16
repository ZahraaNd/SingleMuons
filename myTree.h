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
  TBranch        *b_pthat;

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
  virtual TObjArray* eff (string effType,string Pbp,TH2D *num, TH2D *denom);
  virtual void       AccCalc();
  virtual void       EffCalc();
};
#endif

#ifdef myTree_cxx
/*myTree::myTree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
// TTree* jTree= NULL;
  
if (tree == 0) {
TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/data_CMS/cms/diab/JpsiJet/MC/pp/prompt/v3/HiForestAOD_ext_merged.root");
if (!f || !f->IsOpen()) {
f = new TFile("/data_CMS/cms/diab/JpsiJet/MC/pp/prompt/v3/HiForestAOD_ext_merged.root");
}
      
TDirectory * dir = (TDirectory*)f->Get("/data_CMS/cms/diab/JpsiJet/MC/pp/prompt/v3/HiForestAOD_ext_merged.root:/hionia");
TDirectory * dir1 = (TDirectory*)f->Get("/data_CMS/cms/diab/JpsiJet/MC/pp/prompt/v3/HiForestAOD_ext_merged.root:/ak3PFJetAnalyzer");
dir->GetObject("myTree",tree);      
//dir1->GetObject("t",jTree);
}
// tree->AddFriend(jTree);   
Init(tree);
}*/

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
    "/data_CMS/cms/diab/JpsiJet/MC/PbPb/prompt/v7/HiForestAOD_merged.root", 
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
  if (fChain->GetBranch("Reco_QQ_mupl_dxy_muonlessVtx")) fChain->SetBranchAddress("Reco_QQ_mupl_dxy_muonlessVtx", Reco_QQ_mupl_dxy_muonlessVtx, &b_Reco_QQ_mupl_dxy_muonlessVtx);
  if (fChain->GetBranch("Reco_QQ_mumi_dxy_muonlessVtx")) fChain->SetBranchAddress("Reco_QQ_mumi_dxy_muonlessVtx", Reco_QQ_mumi_dxy_muonlessVtx, &b_Reco_QQ_mumi_dxy_muonlessVtx);
  if (fChain->GetBranch("Reco_QQ_mupl_dxyErr_muonlessVtx")) fChain->SetBranchAddress("Reco_QQ_mupl_dxyErr_muonlessVtx", Reco_QQ_mupl_dxyErr_muonlessVtx, &b_Reco_QQ_mupl_dxyErr_muonlessVtx);
  if (fChain->GetBranch("Reco_QQ_mumi_dxyErr_muonlessVtx")) fChain->SetBranchAddress("Reco_QQ_mumi_dxyErr_muonlessVtx", Reco_QQ_mumi_dxyErr_muonlessVtx, &b_Reco_QQ_mumi_dxyErr_muonlessVtx);
  if (fChain->GetBranch("Reco_QQ_mupl_dz_muonlessVtx")) fChain->SetBranchAddress("Reco_QQ_mupl_dz_muonlessVtx", Reco_QQ_mupl_dz_muonlessVtx, &b_Reco_QQ_mupl_dz_muonlessVtx);
  if (fChain->GetBranch("Reco_QQ_mumi_dz_muonlessVtx")) fChain->SetBranchAddress("Reco_QQ_mumi_dz_muonlessVtx", Reco_QQ_mumi_dz_muonlessVtx, &b_Reco_QQ_mumi_dz_muonlessVtx);
  if (fChain->GetBranch("Reco_QQ_mupl_dzErr_muonlessVtx")) fChain->SetBranchAddress("Reco_QQ_mupl_dzErr_muonlessVtx", Reco_QQ_mupl_dzErr_muonlessVtx, &b_Reco_QQ_mupl_dzErr_muonlessVtx);
  if (fChain->GetBranch("Reco_QQ_mumi_dzErr_muonlessVtx")) fChain->SetBranchAddress("Reco_QQ_mumi_dzErr_muonlessVtx", Reco_QQ_mumi_dzErr_muonlessVtx, &b_Reco_QQ_mumi_dzErr_muonlessVtx);
  if (fChain->GetBranch("Reco_mu_size")) fChain->SetBranchAddress("Reco_mu_size", &Reco_mu_size, &b_Reco_mu_size);
  if (fChain->GetBranch("Reco_mu_type")) fChain->SetBranchAddress("Reco_mu_type", Reco_mu_type, &b_Reco_mu_type);
  if (fChain->GetBranch("Reco_mu_whichGen")) fChain->SetBranchAddress("Reco_mu_whichGen", Reco_mu_whichGen, &b_Reco_mu_whichGen);
  if (fChain->GetBranch("Reco_mu_charge")) fChain->SetBranchAddress("Reco_mu_charge", Reco_mu_charge, &b_Reco_mu_charge);
  if (fChain->GetBranch("Reco_mu_4mom")) fChain->SetBranchAddress("Reco_mu_4mom", &Reco_mu_4mom, &b_Reco_mu_4mom);
  if (fChain->GetBranch("Reco_mu_trig")) fChain->SetBranchAddress("Reco_mu_trig", Reco_mu_trig, &b_Reco_mu_trig);
  if (fChain->GetBranch("Reco_mu_highPurity")) fChain->SetBranchAddress("Reco_mu_highPurity", Reco_mu_highPurity, &b_Reco_mu_highPurity);
  if (fChain->GetBranch("Reco_mu_TrkMuArb")) fChain->SetBranchAddress("Reco_mu_TrkMuArb", Reco_mu_TrkMuArb, &b_Reco_mu_TrkMuArb);
  if (fChain->GetBranch("Reco_mu_TMOneStaTight")) fChain->SetBranchAddress("Reco_mu_TMOneStaTight", Reco_mu_TMOneStaTight, &b_Reco_mu_TMOneStaTight);
  if (fChain->GetBranch("Reco_mu_nPixValHits")) fChain->SetBranchAddress("Reco_mu_nPixValHits", Reco_mu_nPixValHits, &b_Reco_mu_nPixValHits);
  if (fChain->GetBranch("Reco_mu_nMuValHits")) fChain->SetBranchAddress("Reco_mu_nMuValHits", Reco_mu_nMuValHits, &b_Reco_mu_nMuValHits);
  if (fChain->GetBranch("Reco_mu_nTrkHits")) fChain->SetBranchAddress("Reco_mu_nTrkHits", Reco_mu_nTrkHits, &b_Reco_mu_nTrkHits);
  if (fChain->GetBranch("Reco_mu_normChi2_inner")) fChain->SetBranchAddress("Reco_mu_normChi2_inner", Reco_mu_normChi2_inner, &b_Reco_mu_normChi2_inner);
  if (fChain->GetBranch("Reco_mu_normChi2_global")) fChain->SetBranchAddress("Reco_mu_normChi2_global", Reco_mu_normChi2_global, &b_Reco_mu_normChi2_global);
  if (fChain->GetBranch("Reco_mu_nPixWMea")) fChain->SetBranchAddress("Reco_mu_nPixWMea", Reco_mu_nPixWMea, &b_Reco_mu_nPixWMea);
  if (fChain->GetBranch("Reco_mu_nTrkWMea")) fChain->SetBranchAddress("Reco_mu_nTrkWMea", Reco_mu_nTrkWMea, &b_Reco_mu_nTrkWMea);
  if (fChain->GetBranch("Reco_mu_StationsMatched")) fChain->SetBranchAddress("Reco_mu_StationsMatched", Reco_mu_StationsMatched, &b_Reco_mu_StationsMatched);
  if (fChain->GetBranch("Reco_mu_dxy")) fChain->SetBranchAddress("Reco_mu_dxy", Reco_mu_dxy, &b_Reco_mu_dxy);
  if (fChain->GetBranch("Reco_mu_dxyErr")) fChain->SetBranchAddress("Reco_mu_dxyErr", Reco_mu_dxyErr, &b_Reco_mu_dxyErr);
  if (fChain->GetBranch("Reco_,u_dz")) fChain->SetBranchAddress("Reco_mu_dz", Reco_mu_dz, &b_Reco_mu_dz);
  if (fChain->GetBranch("Reco_mu_dzErr")) fChain->SetBranchAddress("Reco_mu_dzErr", Reco_mu_dzErr, &b_Reco_mu_dzErr);
  if (fChain->GetBranch("Reco_mu_pt_inner")) fChain->SetBranchAddress("Reco_mu_pt_inner", Reco_mu_pt_inner, &b_Reco_mu_pt_inner);
  if (fChain->GetBranch("Reco_mu_pt_global")) fChain->SetBranchAddress("Reco_mu_pt_global", Reco_mu_pt_global, &b_Reco_mu_pt_global);
  if (fChain->GetBranch("Reco_mu_ptErr_inner")) fChain->SetBranchAddress("Reco_mu_ptErr_inner", Reco_mu_ptErr_inner, &b_Reco_mu_ptErr_inner);
  if (fChain->GetBranch("Reco_mu_ptErr_global")) fChain->SetBranchAddress("Reco_mu_ptErr_global", Reco_mu_ptErr_global, &b_Reco_mu_ptErr_global);
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
  //fChain->SetBranchAddress("pthat", &pthat, &b_pthat);
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
  if ( ! (Reco_mu_SelectionType[iRecoMu]&((ULong64_t)pow(2, 3))) ) return false; // require the muons to be tracker muons
  if ( ! (Reco_mu_nTrkWMea[iRecoMu] > 5) ) return false;
  if ( ! (Reco_mu_nPixWMea[iRecoMu] > 0) ) return false;
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



#endif // #ifdef myTree_cxx
