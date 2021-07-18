#define treePbPb_cxx
#include "treePbPb.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>

double etaBins[]={0,1.2,1.8,2.1,2.4};
int nEtaBins=sizeof(etaBins)/sizeof(double);
double centBins[]={0,10,20};//,30,40,60,80,100,140,200};    
int nCentBins=sizeof(centBins)/sizeof(double)-1;

int etaBin=0;
double a,b;

char histoName[50];
double cent;
double weight;
  
TH1D** jetPt= new TH1D*[nCentBins];
TH1D** jetPtCr= new TH1D*[nCentBins];
double eff;

TFile* f= TFile::Open("PbPbEff.root");
  
void treePbPb::jetSp()
{  
  TFile* PbPbjetPt=new TFile("PbPbjetPt.root","recreate");
  
  for(int k=0;k<nCentBins;k++){  
    cout<<"k="<<k<<endl;
    sprintf(histoName,"jetPt%d",k);
    jetPt[k]= new TH1D(histoName, "jet pt spectrum without correction PbPb",200,0,200);
    sprintf(histoName,"jetPtCr%d",k);
    jetPtCr[k]= new TH1D(histoName, "jet pt spectrum with correction PbPb",200,0,200);

    if (fChain == 0) return;
    Long64_t nentries = fChain->GetEntriesFast();
    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      
      cent = getMCHiBinFromhiHF(hiHF);
      weight = findNcoll(cent);
      if(!(cent>centBins[k] && cent<centBins[k+1])) continue;
      
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

      if(etaBin!=-1){
	TH1D* histoCr=(TH1D*) f->Get(Form("Qua%d%d",etaBin,k));
	eff=histoCr->GetBinContent(histoCr->FindBin(mupt));;
	jetPt[k]->Fill(jtpt,singleMuWeight*weight);
	if (eff!=0) jetPtCr[k]->Fill(jtpt,singleMuWeight*weight*1/eff);
      }
    }
    jetPt[k]->setXTitle("Pt[GeV]");
    jetPtCr[k]->setXTitle("Pt[GeV]");
    
    jetPt[k]->setYTitle("Number of Jet");
    jetPtCr[k]->setYTitle("Number of Jets");
    
    jetPt[k]->Write();
    jetPtCr[k]->Write();
  }
  PbPbjetPt->Close();
}
