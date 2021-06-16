#define myTree_cxx
#include "myTree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

double ptBins []= {0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 17.5, 18.0, 18.5, 19.0, 19.
5, 20.0, 20.25, 20.5, 20.75, 21.0, 21.25, 21.5, 21.75, 22.0, 22.25, 22.5, 22.75, 23.0, 23.25, 23.5, 23.75, 24.0, 24.5, 25.0};
int nPtBins=sizeof(ptBins)/sizeof(double)-1;
//double etaBins []= {-2.4,-2.1,-1.8,-1.2,0,1.2,1.8,2.1,2.4};
//int nEtaBins=sizeof(etaBins)/sizeof(double)-1;
double etaBins[5]={0,1.2,1.8,2.1,2.4};

TObjArray* myTree::eff (string effType,string Pbp,TH2D *num,TH2D *denom)
{
  TObjArray *arrEff=new TObjArray();
  double BinNb;
  double sumNum=0;
  double sumDenom=0;
  double sumNumErr=0;
  double sumDenomErr=0;
  double numErr;
  double denomErr;
  
  //cout<<"before eta loop"<<endl;
  for(int i=0;i<5;i++){ 
    char histoName[50];
    char histoTitle[50];
    TH1D *histoEff=new TH1D("","",nPtBins,ptBins);
    
    if(etaBins[(i+1)%5]<etaBins[i]){
      double a=etaBins[(i+1)%5];
      etaBins[(i+1)%5]=etaBins[i];
      etaBins[i]=a;
    }
    
    sprintf(histoName,"%sEff%s%d",effType.c_str(),Pbp.c_str(),i);
    sprintf(histoTitle,"%s Efficiency For %s Collisions Over Abs(Eta): %f,%f",effType.c_str(),Pbp.c_str(),etaBins[i],etaBins[(i+1)%5]);
    histoEff->SetName(histoName);
    histoEff->SetTitle(histoTitle);
    cout<<"eta range "<<i<<endl;
    cout<<"Name:"<<histoEff->GetName()<<endl;
    cout<<"Title:"<<histoEff->GetTitle()<<endl;
    for(int j=0;j<nPtBins;j++){
      sumNum=0;
      sumDenom=0;
      sumNumErr=0;
      sumDenomErr=0;
      cout<<"Pt"<<ptBins[j]<<endl;
      for(BinNb=num->FindBin(-etaBins[(i+1)%5],ptBins[j]);BinNb<num->FindBin(-etaBins[i],ptBins[j]);BinNb++){
	cout<<"summing over negative eta range"<<endl;
	sumNum+=+num->GetBinContent(BinNb);
	sumDenom+=denom->GetBinContent(BinNb);
	sumNumErr+=pow(num->GetBinError(BinNb),2);
	sumDenomErr+=pow(denom->GetBinError(BinNb),2);
      }
      for(BinNb=num->FindBin(etaBins[i],ptBins[j]);BinNb<num->FindBin(etaBins[(i+1)%5],ptBins[j]);BinNb++){
	cout<<"Summing over positive eta range"<<endl;
	sumNum+=num->GetBinContent(BinNb);
	sumDenom+=denom->GetBinContent(BinNb);
	sumNumErr+=pow(num->GetBinError(BinNb),2);
	sumDenomErr+=pow(denom->GetBinError(BinNb),2);
      }
      numErr=sqrt(sumNumErr);
      denomErr=sqrt(sumDenomErr);
      cout<<"filling"<<endl;
    
      if((sumNum==0) || (sumDenom==0) || (numErr==0) || (denomErr==0)) continue;
      histoEff->SetBinContent(histoEff->FindBin(ptBins[j]),sumNum/sumDenom);
      cout<<sumNum/sumDenom<<endl;
      histoEff->SetBinError(histoEff->FindBin(ptBins[j]),sqrt((pow(numErr/sumNum,2)+pow(denomErr/sumDenom,2))*pow(sumNum/sumDenom,2)));
      cout<<sqrt((pow(numErr/sumNum,2)+pow(denomErr/sumDenom,2))*pow(sumNum/sumDenom,2))<<endl;        
    }
    cout<<"addingto arrEff"<<endl;
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
  TH2D* histoAcc= new TH2D("histoAcc", "Accepted muons",48,-2.4,2.4,nPtBins,ptBins);//accepted gen muons
  TH2D* histoReco= new TH2D("histoReco", "Reconstructed Muons",48,-2.4,2.4,nPtBins,ptBins);//accepted matched gen muons
  TH2D* histoQua= new TH2D("histoQua", "Muons Passing Quality Cuts",48,-2.4,2.4,nPtBins,ptBins);
  TH2D* histoTrig= new TH2D("histoTrig", "Muons Passing Trigger",48,-2.4,2.4,nPtBins,ptBins);
  TH2D* histoRes= new TH2D("histoRes", "Total Muon Efficiency With Resolution",48,-2.4,2.4,nPtBins,ptBins);
  
  //double Gen_weight=1;
  
  int HLT_HIL3Mu5_v1;
  HLT_HIL3Mu5_v1=(isPbPb?24:17);
  cout<<"HLT_HIL3Mu5_v1="<<HLT_HIL3Mu5_v1<<endl;
  
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    if (LoadTree(jentry) < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if(jentry%1000000==0) cout<<"jentry="<<jentry<<"/"<<nentries<<endl;
    for(int iMu=0; iMu<Gen_mu_size; iMu++){
      TLorentzVector *GenMu4mom = (TLorentzVector*) Gen_mu_4mom->At(iMu);
      if(!isMuonInAccept2019(GenMu4mom)) continue;
      histoAcc->Fill(GenMu4mom->Eta(),GenMu4mom->Pt());
      if (Gen_mu_whichRec[iMu]<0) continue;
      TLorentzVector *RecoMu4mom = (TLorentzVector*) Reco_mu_4mom->At(Gen_mu_whichRec[iMu]);
      if(!(isMuonInAccept2019(RecoMu4mom))) continue;
      histoReco->Fill(GenMu4mom->Eta(),GenMu4mom->Pt());
      if(!(passQualityCuts2019(Gen_mu_whichRec[iMu]))) continue;
      histoQua->Fill(GenMu4mom->Eta(),GenMu4mom->Pt());  
      if(!(isTriggerMatch(Gen_mu_whichRec[iMu],HLT_HIL3Mu5_v1))) continue;
      histoTrig->Fill(GenMu4mom->Eta(),GenMu4mom->Pt());
      histoRes->Fill(RecoMu4mom->Eta(),RecoMu4mom->Pt());
    }//GenMusizeloop
  }//jentryloop
  
  //construct 1D efficiency histograms
  TObjArray *histoRecoEff=eff("Reco",(isPbPb?"PbPb":"pp"),histoReco,histoAcc);
  TObjArray *histoQuaEff=eff("Qua",(isPbPb?"PbPb":"pp"),histoQua,histoReco);
  TObjArray *histoTrigEff=eff("Trig",(isPbPb?"PbPb":"pp"),histoTrig,histoQua);
  TObjArray *histoEff=eff("Tot",(isPbPb?"PbPb":"pp"),histoTrig,histoAcc);
  TObjArray *histoResEff=eff("Res",(isPbPb?"PbPb":"pp"),histoRes,histoAcc);
  cout<<"7";
  
  //save histox 
  TFile* Efficiency= new TFile("Efficiency.root","recreate");
  histoAcc->Write();
  histoReco->Write();
  histoQua->Write();
  histoTrig->Write();
  histoRes->Write();
  histoRecoEff->Write();
  histoQuaEff->Write();
  histoTrigEff->Write();
  histoEff->Write();
  histoResEff->Write();
  Efficiency->Close();
};

void myTree::AccCalc(){
  
  if (fChain == 0) return;
  
  if (!isAcc) {cout<<"[ERROR] you're trying to make Acceptance with Efficiency trees."<<endl; return;}
   
  Long64_t nentries = fChain->GetEntries();
  Long64_t nbytes = 0, nb = 0;
  TH2D* histoTot= new TH2D("histoTot", "All muons",48,-2.4,2.4,nPtBins,ptBins);//accepted gen muons
  TH2D* histoAcc= new TH2D("histoAcc", "Accepted muons",48,-2.4,2.4,nPtBins,ptBins);//accepted gen muons

  nentries=1000000;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    if (LoadTree(jentry) < 0) break;
    nentries=1000000;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if(jentry%1000000==0) cout<<"jentry="<<jentry<<"/"<<nentries<<endl;
    for(int iMu=0; iMu<Gen_mu_size; iMu++){
      TLorentzVector *GenMu4mom = (TLorentzVector*) Gen_mu_4mom->At(iMu);
      histoTot->Fill(GenMu4mom->Eta(),GenMu4mom->Pt());
      if(!isMuonInAccept2019(GenMu4mom)) continue;
      histoAcc->Fill(GenMu4mom->Eta(),GenMu4mom->Pt());
    }
  }
  TObjArray *histoAccEff=eff("Acceptance",(isPbPb?"PbPb":"pp"),histoAcc,histoTot);
  TFile* Acceptance= new TFile("Acceptance.root","recreate");
  histoAccEff->Write();
  Acceptance->Close();
};
