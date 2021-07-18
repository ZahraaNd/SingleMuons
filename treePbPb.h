//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Jul 13 19:07:21 2021 by ROOT version 6.14/04
// from TTree tree/jet tree for LT method
// found on file: /data_CMS/cms/mnguyen/bJetRaa/makeFlatTrees/fillPbPb/data/jetTree_JP_PbPb_DATA_HISingleMuonPD.root
//////////////////////////////////////////////////////////

#ifndef treePbPb_h
#define treePbPb_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class treePbPb {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           nref;
   Int_t           run;
   Int_t           evt;
   Int_t           lumi;
   Int_t           hiBin;
   Float_t         vz;
   Float_t         hiHF;
   Float_t         rawpt;
   Float_t         jtpt;
   Float_t         jteta;
   Float_t         jtphi;
   Float_t         discr_csvV2;
   Float_t         discr_deepCSV;
   Float_t         pdiscr_csvV2;
   Float_t         ndiscr_csvV2;
   Float_t         discr_prob;
   Float_t         mupt;
   Float_t         mueta;
   Float_t         singleMuWeight;
   Float_t         jetWeight40;
   Float_t         jetWeight60;
   Int_t           HLT_Mu3;
   Int_t           HLT_Mu5;
   Int_t           HLT_Mu7;
   Int_t           HLT_Mu7_PS;
   Int_t           HLT_MuxJet40;
   Int_t           HLT_MuxJet60;
   Int_t           HLT_MuxJet80;
   Int_t           HLT_MuxJet100;
   Int_t           HLT_Jet40;
   Int_t           HLT_Jet60;
   Int_t           HLT_Jet80;
   Int_t           HLT_Jet100;
   Float_t         onlineMuPt;
   Float_t         onlineJetPt;
   Float_t         onlineJetEta;
   Float_t         onlineJetPhi;
   Float_t         onlineJetMuPt;
   Float_t         onlineJetMuEta;
   Float_t         onlineJetMuPhi;

   // List of branches
   TBranch        *b_nref;   //!
   TBranch        *b_run;   //!
   TBranch        *b_evt;   //!
   TBranch        *b_lumi;   //!
   TBranch        *b_hiBin;   //!
   TBranch        *b_vz;   //!
   TBranch        *b_hiHF;   //!
   TBranch        *b_rawpt;   //!
   TBranch        *b_jtpt;   //!
   TBranch        *b_jteta;   //!
   TBranch        *b_jtphi;   //!
   TBranch        *b_discr_csvV2;   //!
   TBranch        *b_discr_deepCSV;   //!
   TBranch        *b_pdiscr_csvV2;   //!
   TBranch        *b_ndiscr_csvV2;   //!
   TBranch        *b_discr_prob;   //!
   TBranch        *b_mupt;   //!
   TBranch        *b_mueta;   //!
   TBranch        *b_singleMuWeight;   //!
   TBranch        *b_jetWeight40;   //!
   TBranch        *b_jetWeight60;   //!
   TBranch        *b_HLT_Mu3;   //!
   TBranch        *b_HLT_Mu5;   //!
   TBranch        *b_HLT_Mu7;   //!
   TBranch        *b_HLT_Mu7_PS;   //!
   TBranch        *b_HLT_MuxJet40;   //!
   TBranch        *b_HLT_MuxJet60;   //!
   TBranch        *b_HLT_MuxJet80;   //!
   TBranch        *b_HLT_MuxJet100;   //!
   TBranch        *b_HLT_Jet40;   //!
   TBranch        *b_HLT_Jet60;   //!
   TBranch        *b_HLT_Jet80;   //!
   TBranch        *b_HLT_Jet100;   //!
   TBranch        *b_onlineMuPt;   //!
   TBranch        *b_onlineJetPt;   //!
   TBranch        *b_onlineJetEta;   //!
   TBranch        *b_onlineJetPhi;   //!
   TBranch        *b_onlineJetMuPt;   //!
   TBranch        *b_onlineJetMuEta;   //!
   TBranch        *b_onlineJetMuPhi;   //!

   treePbPb(TTree *tree=0);
   virtual ~treePbPb();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     jetSp();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   virtual Int_t    getMCHiBinFromhiHF(const Double_t hiHF); 
   virtual Double_t findNcoll(int hiBin);
};

#endif

#ifdef treePbPb_cxx
treePbPb::treePbPb(TTree *treePbPb) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (treePbPb == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/data_CMS/cms/mnguyen/bJetRaa/makeFlatTrees/fillPbPb/data/jetTree_JP_PbPb_DATA_H
ISingleMuonPD.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/data_CMS/cms/mnguyen/bJetRaa/makeFlatTrees/fillPbPb/data/jetTree_JP_PbPb_DATA_HISingleMuonPD.root");
      }
      f->GetObject("tree",treePbPb);

   }
   Init(treePbPb);
}

treePbPb::~treePbPb()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t treePbPb::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t treePbPb::LoadTree(Long64_t entry)
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

void treePbPb::Init(TTree *tree)
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
   fChain->SetBranchAddress("hiBin", &hiBin, &b_hiBin);
   fChain->SetBranchAddress("vz", &vz, &b_vz);
   fChain->SetBranchAddress("hiHF", &hiHF, &b_hiHF);
   fChain->SetBranchAddress("rawpt", &rawpt, &b_rawpt);
   fChain->SetBranchAddress("jtpt", &jtpt, &b_jtpt);
   fChain->SetBranchAddress("jteta", &jteta, &b_jteta);
   fChain->SetBranchAddress("jtphi", &jtphi, &b_jtphi);
   fChain->SetBranchAddress("discr_csvV2", &discr_csvV2, &b_discr_csvV2);
   fChain->SetBranchAddress("discr_deepCSV", &discr_deepCSV, &b_discr_deepCSV);
   fChain->SetBranchAddress("pdiscr_csvV2", &pdiscr_csvV2, &b_pdiscr_csvV2);
   fChain->SetBranchAddress("ndiscr_csvV2", &ndiscr_csvV2, &b_ndiscr_csvV2);
   fChain->SetBranchAddress("discr_prob", &discr_prob, &b_discr_prob);
   fChain->SetBranchAddress("mupt", &mupt, &b_mupt);
   fChain->SetBranchAddress("mueta", &mueta, &b_mueta);
   fChain->SetBranchAddress("singleMuWeight", &singleMuWeight, &b_singleMuWeight);
   fChain->SetBranchAddress("jetWeight40", &jetWeight40, &b_jetWeight40);
   fChain->SetBranchAddress("jetWeight60", &jetWeight60, &b_jetWeight60);
   fChain->SetBranchAddress("HLT_Mu3", &HLT_Mu3, &b_HLT_Mu3);
   fChain->SetBranchAddress("HLT_Mu5", &HLT_Mu5, &b_HLT_Mu5);
   fChain->SetBranchAddress("HLT_Mu7", &HLT_Mu7, &b_HLT_Mu7);
   fChain->SetBranchAddress("HLT_Mu7_PS", &HLT_Mu7_PS, &b_HLT_Mu7_PS);
   fChain->SetBranchAddress("HLT_MuxJet40", &HLT_MuxJet40, &b_HLT_MuxJet40);
   fChain->SetBranchAddress("HLT_MuxJet60", &HLT_MuxJet60, &b_HLT_MuxJet60);
   fChain->SetBranchAddress("HLT_MuxJet80", &HLT_MuxJet80, &b_HLT_MuxJet80);
   fChain->SetBranchAddress("HLT_MuxJet100", &HLT_MuxJet100, &b_HLT_MuxJet100);
   fChain->SetBranchAddress("HLT_Jet40", &HLT_Jet40, &b_HLT_Jet40);
   fChain->SetBranchAddress("HLT_Jet60", &HLT_Jet60, &b_HLT_Jet60);
   fChain->SetBranchAddress("HLT_Jet80", &HLT_Jet80, &b_HLT_Jet80);
   fChain->SetBranchAddress("HLT_Jet100", &HLT_Jet100, &b_HLT_Jet100);
   fChain->SetBranchAddress("onlineMuPt", &onlineMuPt, &b_onlineMuPt);
   fChain->SetBranchAddress("onlineJetPt", &onlineJetPt, &b_onlineJetPt);
   fChain->SetBranchAddress("onlineJetEta", &onlineJetEta, &b_onlineJetEta);
   fChain->SetBranchAddress("onlineJetPhi", &onlineJetPhi, &b_onlineJetPhi);
   fChain->SetBranchAddress("onlineJetMuPt", &onlineJetMuPt, &b_onlineJetMuPt);
   fChain->SetBranchAddress("onlineJetMuEta", &onlineJetMuEta, &b_onlineJetMuEta);
   fChain->SetBranchAddress("onlineJetMuPhi", &onlineJetMuPhi, &b_onlineJetMuPhi);
   Notify();
}

Bool_t treePbPb::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void treePbPb::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t treePbPb::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

Int_t treePbPb::getMCHiBinFromhiHF(const Double_t hiHF) {
  const Int_t nBins = 200; // table of bin edges     
  const Double_t binTable[nBins+1] = {0, 12.2187, 13.0371, 13.7674, 14.5129, 15.2603, 16.0086, 16.7623, 17.5335, 18.3283, 19.1596, 19.9989, 20
.8532, 21.7297, 22.6773, 23.6313, 24.6208, 25.6155, 26.6585, 27.7223, 28.8632, 30.041, 31.2865, 32.5431, 33.8655, 35.2539, 36.6912, 38.2064,39
.7876, 41.4818, 43.2416, 45.0605, 46.9652, 48.9918, 51.1, 53.2417, 55.5094, 57.9209, 60.3817, 62.9778, 65.6099, 68.4352, 71.3543, 74.4154,77.6
252, 80.8425, 84.1611, 87.7395, 91.3973, 95.1286, 99.0571, 103.185, 107.482, 111.929, 116.45, 121.178, 126.081, 130.995, 136.171, 141.612, 147
.298, 153.139, 159.419, 165.633, 172.114, 178.881, 185.844, 192.845, 200.244, 207.83, 215.529, 223.489, 231.878, 240.254, 249.319, 258.303, 26
7.508, 277.037, 286.729, 296.845, 307.458, 317.882, 328.787, 340.074, 351.295, 362.979, 375.125, 387.197, 399.604, 412.516, 425.683, 439.001, 
452.667, 466.816, 481.007, 495.679, 510.588, 526.138, 541.782, 557.641, 574.141, 591.071, 608.379, 626.068, 643.616, 661.885, 680.288, 699.449
, 718.925, 738.968, 758.983, 779.459, 800.376, 821.638, 843.555, 865.771, 888.339, 911.031, 934.979, 958.56, 982.582, 1007.02, 1031.9, 1057.81
, 1084.01, 1111.71, 1138.21, 1165.72, 1193.73, 1221.65, 1251.51, 1281.23, 1311.01, 1341.1, 1372.4, 1404.29, 1436.52, 1468.65, 1501.91, 1535.56
, 1569.69, 1604.69, 1640.65, 1676.05, 1712.62, 1749.28, 1787.43, 1825.89, 1866.07, 1906.58, 1947.84, 1989.66, 2031.4, 2072.8, 2115.32, 2159.5,
 2205.23, 2252.68, 2298.58, 2345.65, 2393.36, 2442.87, 2491.45, 2541.04, 2592.81, 2645.52, 2699.1, 2753.29, 2807.93, 2864.37, 2922.6, 2979.42,
 3038.68, 3098.72, 3159.29, 3221.66, 3285.9, 3350.95, 3415.81, 3482.69, 3552.62, 3623.61, 3694.63, 3767.25, 3840.28, 3917.04, 3993.66, 4073.36
, 4154.33, 4238.13, 4322.21, 4409.83, 4498.89, 4589.72, 4681.56, 4777.09, 4877.95, 4987.05, 5113.04, 5279.58, 6242.82};

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

Double_t treePbPb::findNcoll(int hiBin) {
  const int nbins = 200;
  const Double_t Ncoll[nbins] = {1976.95, 1944.02, 1927.29, 1891.9, 1845.3, 1807.2, 1760.45, 1729.18, 1674.8, 1630.3, 1590.52, 1561.72, 1516.1
, 1486.5, 1444.68, 1410.88, 1376.4, 1347.32, 1309.71, 1279.98, 1255.31, 1219.89, 1195.13, 1165.96, 1138.92, 1113.37, 1082.26, 1062.42, 1030.6,
 1009.96, 980.229, 955.443, 936.501, 915.97, 892.063, 871.289, 847.364, 825.127, 806.584, 789.163, 765.42, 751.187, 733.001, 708.31, 690.972, 
677.711, 660.682, 640.431, 623.839, 607.456, 593.307, 576.364, 560.967, 548.909, 530.475, 519.575, 505.105, 490.027, 478.133, 462.372, 451.115
, 442.642, 425.76, 416.364, 405.154, 392.688, 380.565, 371.167, 360.28, 348.239, 340.587, 328.746, 320.268, 311.752, 300.742, 292.172, 281.361
, 274.249, 267.025, 258.625, 249.931, 240.497, 235.423, 228.63, 219.854, 214.004, 205.425, 199.114, 193.618, 185.644, 180.923, 174.289, 169.64
1, 161.016, 157.398, 152.151, 147.425, 140.933, 135.924, 132.365, 127.017, 122.127, 117.817, 113.076, 109.055, 105.16, 101.323, 98.098,95.0548
, 90.729, 87.6495, 84.0899, 80.2237, 77.2201, 74.8848, 71.3554, 68.7745, 65.9911, 63.4136, 61.3859, 58.1903, 56.4155, 53.8486, 52.0196, 49.292
1, 47.0735, 45.4345, 43.8434, 41.7181, 39.8988, 38.2262, 36.4435, 34.8984, 33.4664, 31.8056, 30.351, 29.2074, 27.6924, 26.7754, 25.4965, 24.28
02, 22.9651, 22.0059, 21.0915, 19.9129, 19.1041, 18.1487, 17.3218, 16.5957, 15.5323, 14.8035, 14.2514, 13.3782, 12.8667, 12.2891, 11.61, 11.00
26, 10.3747, 9.90294, 9.42648, 8.85324, 8.50121, 7.89834, 7.65197, 7.22768, 6.7755, 6.34855, 5.98336, 5.76555, 5.38056, 5.11024, 4.7748, 4.591
17, 4.23247, 4.00814, 3.79607, 3.68702, 3.3767, 3.16309, 2.98282, 2.8095, 2.65875, 2.50561, 2.32516, 2.16357, 2.03235, 1.84061, 1.72628, 1.623
05, 1.48916, 1.38784, 1.28366, 1.24693, 1.18552, 1.16085, 1.12596, 1.09298, 1.07402, 1.06105, 1.02954};
  return Ncoll[hiBin];
};

#endif // #ifdef tree_cxx
