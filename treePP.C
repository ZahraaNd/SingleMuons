#define treePP_cxx
#include "treePP.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>

double etaBins[]={0,1.2,1.8,2.1,2.4};
int nEtaBins=sizeof(etaBins)/sizeof(double)-1;

double a,b;
int etaBin=-1;

TH1D* jetPt=new TH1D("jetPt","jet pt spectrum without correction",200,0,200);
  TH1D* jetPtCr=new TH1D("jetPtCr","jet pt spectrum with correction",200,0,200);
double eff=1;
TFile* f= TFile::Open("PPEff.root");
  
void treePP::jetSp()
{
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    for(int i=0;i<5;i++){
      a=etaBins[i];
      b=etaBins[(i+1)%nEtaBins];
      if(b<a) {
	int c=a;
	a=b;
	b=c;
      }
      if(a<mueta<b) etaBin=i;
    }
    if(etaBin!=-1) {
      TH1D* histoCr=(TH1D*) f->Get(Form("Res%d0",etaBin));
      eff=histoCr->GetBinContent(histoCr->FindBin(mupt));
      jetPt->Fill(jtpt,singleMuWeight);
      if(eff!=0) jetPtCr->Fill(jtpt,singleMuWeight*1/eff);
    }
  }
  TFile* PPjetPt=new TFile("PPjetPt.root","recreate");
  jetPt->Write();
  jetPtCr->Write();
  PPjetPt->Close();
}

