#define myTree_cxx
#include "myTree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "tnp_weight_lowptPbPb.h"
#include "tnp_weight_lowptpp.h"

double ptBins []= {0,0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.25,2.5,2.75,3,3.25,3.5,3.75,4,4.5,5,7,10,15,30};
int nPtBins=sizeof(ptBins)/sizeof(double)-1;
double etaBins[]={0,1.2,1.8,2.1,2.4};
int nEtaBins=sizeof(etaBins)/sizeof(double);
double centBins[]={0,10,20};//,30,40,60,80,100,140,200};
int nCentBins=sizeof(centBins)/sizeof(double)-1;

int cent;
double weight;  

double tnp_weight_trk;
double tnp_weight;

double a,b; 
    
char histoName[50];
    
TObjArray* myTree::eff (string effType,string Pbp,TH2D* num,TH2D* denom, int cent){
  TObjArray* arrEff=new TObjArray();
  double BinNb;
  double sumNum=0;
  double sumDenom=0;
  double sumNumErr=0;
  double sumDenomErr=0;
  double numErr;
  double denomErr;
  for(int i=0;i<nEtaBins;i++){ 
    char histoName[50];
    char histoTitle[50];
    TH1D *histoEff=new TH1D("","",nPtBins,ptBins);
    
    if(etaBins[(i+1)%nEtaBins]>etaBins[i]){
      a=etaBins[i];
      b=etaBins[(i+1)%nEtaBins];
    }
    else{
      a=etaBins[(i+1)%nEtaBins];
      b=etaBins[i];
    }
    
    sprintf(histoName,"%s%d%d",effType.c_str(),i,cent);
    //sprintf(histoTitle,"%s Efficiency For %s Collisions Over Abs(Eta): %f,%f",effType.c_str(),Pbp.c_str(),a,b);
    histoEff->SetName(histoName);
    //histoEff->SetTitle(histoTitle);
    
    for(int j=0;j<nPtBins;j++){
      sumNum=0;
      sumDenom=0;
      sumNumErr=0;
      sumDenomErr=0;
    
      for(BinNb=num->FindBin(-b,ptBins[j]);BinNb<num->FindBin(-a,ptBins[j]);BinNb++){
	sumNum+=+num->GetBinContent(BinNb);
	sumDenom+=denom->GetBinContent(BinNb);
	sumNumErr+=pow(num->GetBinError(BinNb),2);
	sumDenomErr+=pow(denom->GetBinError(BinNb),2);
      }

      for(BinNb=num->FindBin(a,ptBins[j]);BinNb<num->FindBin(b,ptBins[j]);BinNb++){
	sumNum+=num->GetBinContent(BinNb);
	sumDenom+=denom->GetBinContent(BinNb);
	sumNumErr+=pow(num->GetBinError(BinNb),2);
	sumDenomErr+=pow(denom->GetBinError(BinNb),2);
      }
      numErr=sqrt(sumNumErr);
      denomErr=sqrt(sumDenomErr);
      
      if((sumNum==0) || (sumDenom==0) || (numErr==0) || (denomErr==0)) continue;
      histoEff->SetBinContent(histoEff->FindBin(ptBins[j]),sumNum/sumDenom);
      histoEff->SetBinError(histoEff->FindBin(ptBins[j]),sqrt((pow(numErr/sumNum,2)+pow(denomErr/sumDenom,2))*pow(sumNum/sumDenom,2)));
    }
    gStyle->SetOptStat(false);
    histoEff->SetXTitle("p_{T}^{#mu}[GeV/c]");
    histoEff->SetYTitle("Single #mu Efficiency");      
    histoEff->GetXaxis()->CenterTitle();
    arrEff->Add(histoEff);
  }
  return arrEff;
};

void myTree::EffCalc()
{
  if (fChain == 0) return;
  if (isAcc) {cout<<"[ERROR] you're trying to make Efficiency with Acceptance trees."<<endl; return;}
  
  Long64_t nentries = fChain->GetEntries();
  Long64_t nbytes = 0, nb = 0;
  //nentries=1000000;
  cout<<nentries<<endl;
  int HLT_HIL3Mu5_v1; 
  HLT_HIL3Mu5_v1=(isPbPb?24:17);
  cout<<"HLT_HIL3Mu5_v1="<<HLT_HIL3Mu5_v1<<endl;
    
  if(!isPbPb) nCentBins=1;
  cout<<"nCentBins="<<nCentBins<<endl;

  TH2D** histoAcc= new TH2D*[nCentBins];
  
  TH2D** histoReco= new TH2D*[nCentBins];
  TH2D** histoQua= new TH2D*[nCentBins];
  
  TH2D** histoRecoTnP= new TH2D*[nCentBins];
  TH2D** histoQuaTnP= new TH2D*[nCentBins];
  
  TH2D** histoTrigTnP= new TH2D*[nCentBins];  
  TH2D** histoResTnP= new TH2D*[nCentBins];
  
  TObjArray** histoRecoEff=new TObjArray*[nCentBins];
  TObjArray** histoQuaEff=new TObjArray*[nCentBins];
  TObjArray** histoTrigEff=new TObjArray*[nCentBins];
  TObjArray** histoEff=new TObjArray*[nCentBins];
  TObjArray** histoResEff=new TObjArray*[nCentBins];
  
  TFile* Efficiency= new TFile((isPbPb?"PbPbEff.root":"PPEff.root"),"recreate");

  for(int k=0;k<nCentBins;k++){

    sprintf(histoName,"histoAcc%d",k);
    histoAcc[k]= new TH2D(histoName, "Accepted muons",48,-2.4,2.4,nPtBins,ptBins);//accepted gen muons
    sprintf(histoName,"histoReco%d",k);
    histoReco[k]= new TH2D(histoName, "Reconstructed Muons",48,-2.4,2.4,nPtBins,ptBins);//accepted matched gen muons
    sprintf(histoName,"histoRecoTnP%d",k);
    histoRecoTnP[k]= new TH2D(histoName, "Reconstructed Muons",48,-2.4,2.4,nPtBins,ptBins);//accepted matched gen muons
    sprintf(histoName,"histoQua%d",k);
    histoQua[k]= new TH2D(histoName, "Muons Passing Quality Cuts",48,-2.4,2.4,nPtBins,ptBins);
    sprintf(histoName,"histoQuaTnP%d",k);
    histoQuaTnP[k]= new TH2D(histoName, "Muons Passing Quality Cuts",48,-2.4,2.4,nPtBins,ptBins);
    sprintf(histoName,"histoTrig%d",k);
    histoTrigTnP[k]= new TH2D(histoName, "Muons Passing Trigger",48,-2.4,2.4,nPtBins,ptBins);
    sprintf(histoName,"histoResTnP%d",k);
    histoResTnP[k]= new TH2D(histoName, "Total Muon Efficiency With Resolution",48,-2.4,2.4,nPtBins,ptBins);

    for (Long64_t jentry=0; jentry<nentries;jentry++) {
      if(jentry==747864) continue;
      if (LoadTree(jentry) < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      if(jentry%1000000==0) cout<<"jentry="<<jentry<<"/"<<nentries<<endl;
      
      if(isPbPb) {
	cent=Centrality;
	weight = findNcoll(cent); 
	if(!(cent>centBins[k] && cent<centBins[(k+1)])) continue;
      }
      else weight=1;
      
      for(int iMu=0; iMu<Gen_mu_size; iMu++){
	
	TLorentzVector* GenMu4mom = (TLorentzVector*) Gen_mu_4mom->At(iMu);
	
	if(isPbPb){
	  tnp_weight_trk=tnp_weight_trk_pbpb(GenMu4mom->Eta(),0);
	  tnp_weight = tnp_weight_trk_pbpb(GenMu4mom->Eta(),0)*tnp_weight_muid_pbpb(GenMu4mom->Pt(),GenMu4mom->Eta(),0);
	}
	else{
	  auto muGlb = tnp_weight_GlobalMuon_TightAcceptance_pp(GenMu4mom->Pt(),GenMu4mom->Eta());
	  auto muIdTrg = tnp_weight_HybridSoftID_LooseAcceptance_pp(GenMu4mom->Pt(),GenMu4mom->Eta());

	  tnp_weight_trk= std::get<0>(muGlb);
	  tnp_weight = std::get<0>(muGlb) * std::get<0>(muIdTrg);
	}
    
	if(!isMuonInAccept2019(GenMu4mom)) continue;
	histoAcc[k]->Fill(GenMu4mom->Eta(),GenMu4mom->Pt(),weight);
	if (Gen_mu_whichRec[iMu]<0) continue;
	TLorentzVector* RecoMu4mom = (TLorentzVector*) Reco_mu_4mom->At(Gen_mu_whichRec[iMu]);
	if(!(isMuonInAccept2019(RecoMu4mom))) continue;
	histoReco[k]->Fill(GenMu4mom->Eta(),GenMu4mom->Pt(),weight);
	histoRecoTnP[k]->Fill(GenMu4mom->Eta(),GenMu4mom->Pt(),weight*tnp_weight_trk);
	if(!(passQualityCuts2019(Gen_mu_whichRec[iMu]))) continue;
	histoQua[k]->Fill(GenMu4mom->Eta(),GenMu4mom->Pt(),weight);
	histoQuaTnP[k]->Fill(GenMu4mom->Eta(),GenMu4mom->Pt(),weight*tnp_weight);
	if(!(isTriggerMatch(Gen_mu_whichRec[iMu],HLT_HIL3Mu5_v1))) continue;
	histoTrigTnP[k]->Fill(GenMu4mom->Eta(),GenMu4mom->Pt(),weight*tnp_weight);
	histoResTnP[k]->Fill(RecoMu4mom->Eta(),RecoMu4mom->Pt(),weight*tnp_weight);
      
      }//GenMusizeloop
    }//jentryloop
   
    histoRecoEff[k]=eff("Reco",(isPbPb?"PbPb":"pp"),histoRecoTnP[k],histoAcc[k],k);
    histoQuaEff[k]=eff("Qua",(isPbPb?"PbPb":"pp"),histoQuaTnP[k],histoReco[k],k);
    histoTrigEff[k]=eff("Trig",(isPbPb?"PbPb":"pp"),histoTrigTnP[k],histoQua[k],k);
    histoEff[k]=eff("Tot",(isPbPb?"PbPb":"pp"),histoTrigTnP[k],histoAcc[k],k);
    histoResEff[k]=eff("Res",(isPbPb?"PbPb":"pp"),histoResTnP[k],histoAcc[k],k);
    
    histoRecoEff[k]->Write();
    histoQuaEff[k]->Write();
    histoTrigEff[k]->Write();
    histoEff[k]->Write();
    histoResEff[k]->Write();
  }
  Efficiency->Close();
} 

void myTree::AccCalc()
{  
  if (fChain == 0) return;
  if (!isAcc) {cout<<"[ERROR] you're trying to make Acceptance with Efficiency trees."<<endl; return;}

  Long64_t nentries = fChain->GetEntries();
  Long64_t nbytes = 0, nb = 0;
  //nentries=1000000;
  TH2D* histoTot= new TH2D(histoName, "All muons",48,-2.4,2.4,nPtBins,ptBins);
  TH2D* histoAcc= new TH2D(histoName, "Accepted muons",48,-2.4,2.4,nPtBins,ptBins);
  
  TFile* Acceptance= new TFile("Acc.root","recreate");
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    
    if (LoadTree(jentry) < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if(jentry%1000000==0) cout<<"jentry="<<jentry<<"/"<<nentries<<endl;
    
    for(int iMu=0; iMu<Gen_mu_size; iMu++){
      
      TLorentzVector* GenMu4mom = (TLorentzVector*) Gen_mu_4mom->At(iMu);
      histoTot->Fill(GenMu4mom->Eta(),GenMu4mom->Pt());
      if(!isMuonInAccept2019(GenMu4mom)) continue;
      histoAcc->Fill(GenMu4mom->Eta(),GenMu4mom->Pt());
    }
  }
  TObjArray* histoAccEff=eff("Acceptance","",histoAcc,histoTot,0);
  histoAccEff->Write();
  Acceptance->Close();
};
