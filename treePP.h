//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Jul 13 19:13:46 2021 by ROOT version 6.14/04
// from TTree tree/jet tree with JP tagger
// found on file: /data_CMS/cms/mnguyen/bJetRaa/makeFlatTrees/fillPP/data/jetTree_JP_PP_DATA_SingleMuonMerged.root
//////////////////////////////////////////////////////////

#ifndef treePP_h
#define treePP_h

#include <TH1.h>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class treePP {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           nref;
   Int_t           run;
   Int_t           evt;
   Int_t           lumi;
   Float_t         vz;
   Float_t         mupt;
   Float_t         mueta;
   Float_t         rawpt;
   Float_t         jtpt;
   Float_t         jteta;
   Float_t         jtphi;
   Float_t         discr_csvV2;
   Float_t         discr_deepCSV;
   Float_t         discr_prob;
   Int_t           HLT_Mu3;
   Int_t           HLT_Mu5;
   Int_t           HLT_Mu7;
   Int_t           HLT_Mu5xJet30;
   Float_t         singleMuWeight;
   Float_t         onlineMuPt;
   
   // List of branches
   TBranch        *b_nref;   //!
   TBranch        *b_run;   //!
   TBranch        *b_evt;   //!
   TBranch        *b_lumi;   //!
   TBranch        *b_vz;   //!
   TBranch        *b_mupt;   //!
   TBranch        *b_mueta;   //!
   TBranch        *b_rawpt;   //!
   TBranch        *b_jtpt;   //!
   TBranch        *b_jteta;   //!
   TBranch        *b_jtphi;   //!
   TBranch        *b_discr_csvV2;   //!
   TBranch        *b_discr_deepCSV;   //!
   TBranch        *b_discr_prob;   //!
   TBranch        *b_HLT_Mu3;   //!
   TBranch        *b_HLT_Mu5;   //!
   TBranch        *b_HLT_Mu7;   //!
   TBranch        *b_HLT_Mu5xJet30;   //!
   TBranch        *b_singleMuWeight;   //!
   TBranch        *b_onlineMuPt;   //!

   treePP(TTree *treePP=0);
   virtual ~treePP();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     jetSp();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef treePP_cxx
treePP::treePP(TTree *treePP) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (treePP == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/data_CMS/cms/mnguyen/bJetRaa/makeFlatTrees/fillPP/data/jetTree_JP_PP_DATA_Singl
eMuonMerged.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/data_CMS/cms/mnguyen/bJetRaa/makeFlatTrees/fillPP/data/jetTree_JP_PP_DATA_SingleMuonMerged.root");
      }
      f->GetObject("tree",treePP);

   }
   Init(treePP);
}

treePP::~treePP()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t treePP::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t treePP::LoadTree(Long64_t entry)
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

void treePP::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("nref", &nref, &b_nref);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("evt", &evt, &b_evt);
   fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
   fChain->SetBranchAddress("vz", &vz, &b_vz);
   fChain->SetBranchAddress("mupt", &mupt, &b_mupt);
   fChain->SetBranchAddress("mueta", &mueta, &b_mueta);
   fChain->SetBranchAddress("rawpt", &rawpt, &b_rawpt);
   fChain->SetBranchAddress("jtpt", &jtpt, &b_jtpt);
   fChain->SetBranchAddress("jteta", &jteta, &b_jteta);
   fChain->SetBranchAddress("jtphi", &jtphi, &b_jtphi);
   fChain->SetBranchAddress("discr_csvV2", &discr_csvV2, &b_discr_csvV2);
   fChain->SetBranchAddress("discr_deepCSV", &discr_deepCSV, &b_discr_deepCSV);
   fChain->SetBranchAddress("discr_prob", &discr_prob, &b_discr_prob);
   fChain->SetBranchAddress("HLT_Mu3", &HLT_Mu3, &b_HLT_Mu3);
   fChain->SetBranchAddress("HLT_Mu5", &HLT_Mu5, &b_HLT_Mu5);
   fChain->SetBranchAddress("HLT_Mu7", &HLT_Mu7, &b_HLT_Mu7);
   fChain->SetBranchAddress("HLT_Mu5xJet30", &HLT_Mu5xJet30, &b_HLT_Mu5xJet30);
   fChain->SetBranchAddress("singleMuWeight", &singleMuWeight, &b_singleMuWeight);
   fChain->SetBranchAddress("onlineMuPt", &onlineMuPt, &b_onlineMuPt);
   Notify();
}

Bool_t treePP::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void treePP::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t treePP::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef tree_cxx
