#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>
#include <map>
#include <sstream> 
#include <cmath>
#include <exception>
#include <ctime>
#include <cmath>
#include "TNtuple.h"
#include "TROOT.h"
#include "TString.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TF1.h"
#include "TF2.h"
#include "TCanvas.h"
#include "TFitResult.h"
#include "TPaveText.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TObject.h"
#include "TObjArray.h"
#include "THStack.h"
#include "TStyle.h"
#include "TROOT.h"
#include "THashList.h"
#include "TApplication.h"
#include "TGraph.h"
#include "TMath.h"
#include "TLegend.h"
#include "TProfile.h"

using namespace std;

int main(int argc, char *argv[]){

  //gStyle->SetOptStat(0);
  //gROOT->ForceStyle();
  TH1::SetDefaultSumw2();
  if(argc != 5)
    throw std::runtime_error("Bad number of arguments!");
    
  TString file1 = argv[1];
  TString subDet = argv[2];
  TString slayer = argv[3];
  TString dir = argv[4];

  int layer = slayer.Atoi();

  uint16_t sdId = 0;

  if(subDet == "TIB")
    sdId = 3;
  else if(subDet == "TID")
    sdId = 4;
  else if(subDet == "TOB")
    sdId = 5;
  else if(subDet == "TEC")
    sdId = 6;
  else
    throw std::runtime_error("Wrong partition entered");
  TFile* f1 = NULL;
  TTree* t1 = NULL;
  f1 = TFile::Open(file1); 
  if(f1==NULL)
    throw std::runtime_error("File 1 address not set");
  t1 = dynamic_cast<TTree*>(f1->Get("testTree/tree"));
  if(t1==NULL)
    throw std::runtime_error("Tree 1 address not set");

  vector<float>* partition = 0;
  vector<float>* partition2 = 0;
  vector<float>* clustercharge = 0;
  vector<float>* clustercharge2 =0;
  vector<float>* clusterwidth = 0;
  vector<float>* clusterwidth2 = 0;
  vector<float>* clusterlayerwheel = 0;
  vector<float>* clusterlayerwheel2 = 0;
  vector<float>* clusterstripChargeSubdetid = 0;
  vector<float>* clusterstripChargeSubdetid2 = 0;
  vector<float>* clusterstripCharge = 0;
  vector<float>* clusterstripCharge2 = 0;
  vector<float>* clusterstripChargeLayerwheel = 0;
  vector<float>* clusterstripChargeLayerwheel2 = 0;
  vector<float>* clusterstripChargeStripNr = 0;
  vector<float>* clusterstripChargeStripNr2 = 0;
  vector<float>* clusterstripChargeTotWidth = 0;
  vector<float>* clusterstripChargeTotWidth2 = 0;
  vector<float>* clusterstripChargeTotCharge = 0;
  vector<float>* clusterstripChargeTotCharge2 = 0;
  vector<float>* clusterstripChargeLocalTrackPhi = 0;
  vector<float>* clusterstripChargeLocalTrackPhi2 = 0;
  vector<float>* clusterstripChargeGlobalTrackPhi = 0;
  vector<float>* clusterstripChargeGlobalTrackPhi2 = 0;
  vector<float>* clusterstripChargeLocalTrackTheta = 0;
  vector<float>* clusterstripChargeLocalTrackTheta2 = 0;
  vector<float>* clusterstripChargeGlobalTrackTheta = 0;
  vector<float>* clusterstripChargeGlobalTrackTheta2 = 0;
  vector<unsigned>* clusterstripChargeDetid = 0;
  vector<unsigned>* clusterstripChargeDetid2 = 0;
  vector<float>* clusterstripChargeLocalX = 0;
  vector<float>* clusterstripChargeLocalX2 = 0;
  vector<float>* clusterstripChargeLocalY = 0;
  vector<float>* clusterstripChargeLocalY2 = 0;
  vector<float>* tsostrackPt = 0;
  vector<float>* tsostrackPt2 = 0;
  vector<float>* clusterstripChargetrackPt = 0;
  vector<float>* clusterstripChargetrackPt2 = 0;
  vector<float>* tsoslocalphi = 0;
  vector<float>* tsoslocalphi2 = 0;
  vector<float>* tsoslocaltheta = 0;
  vector<float>* tsoslocaltheta2 = 0;
  vector<float>* tsoslocalpitch = 0;
  vector<float>* tsoslocalpitch2 = 0;
  vector<float>* clustersensorThickness = 0;
  vector<float>* clustersensorThickness2 = 0;
  vector<float>* tsosBdotY = 0;
  vector<float>* tsosBdotY2 = 0;
  vector<float>* clusterstripChargelocalpitch = 0;
  vector<float>* clusterstripChargelocalpitch2  = 0;
  vector<float>* clusterstripChargesensorThickness = 0;
  vector<float>* clusterstripChargesensorThickness2  = 0;
  vector<float>* clusterstripChargeBdotY = 0;
  vector<float>* clusterstripChargeBdotY2  = 0;


  vector<float> subpartition;
  vector<float> subpartition2;
  vector<float> subclustercharge;
  vector<float> subclustercharge2;
  vector<float> subclusterwidth;
  vector<float> subclusterwidth2;
  vector<float> subclusterlayerwheel;
  vector<float> subclusterlayerwheel2;
  vector<float> subclusterstripChargeSubdetid;
  vector<float> subclusterstripChargeSubdetid2;
  vector<float> subclusterstripCharge;
  vector<float> subclusterstripCharge2;
  vector<float> subclusterstripChargeLayerwheel;
  vector<float> subclusterstripChargeLayerwheel2;
  vector<float> subclusterstripChargeStripNr;
  vector<float> subclusterstripChargeStripNr2;
  vector<float> subclusterstripChargeTotWidth;
  vector<float> subclusterstripChargeTotWidth2;
  vector<float> subclusterstripChargeTotCharge;
  vector<float> subclusterstripChargeTotCharge2;
  vector<float> subclusterstripChargeLocalTrackPhi;
  vector<float> subclusterstripChargeLocalTrackPhi2;
  vector<float> subclusterstripChargeGlobalTrackPhi;
  vector<float> subclusterstripChargeGlobalTrackPhi2;
  vector<float> subclusterstripChargeLocalTrackTheta;
  vector<float> subclusterstripChargeLocalTrackTheta2;
  vector<float> subclusterstripChargeGlobalTrackTheta;
  vector<float> subclusterstripChargeGlobalTrackTheta2;
  vector<unsigned> subclusterstripChargeDetid;
  vector<unsigned> subclusterstripChargeDetid2;
  vector<float> subclusterstripChargeLocalX;
  vector<float> subclusterstripChargeLocalX2;
  vector<float> subclusterstripChargeLocalY;
  vector<float> subclusterstripChargeLocalY2;
  vector<float> subtsostrackPt;
  vector<float> subtsostrackPt2;
  vector<float> subclusterstripChargetrackPt;
  vector<float> subclusterstripChargetrackPt2;
  vector<float> subtsoslocalphi;
  vector<float> subtsoslocalphi2;
  vector<float> subtsoslocaltheta;
  vector<float> subtsoslocaltheta2;
  vector<float> subtsoslocalpitch;
  vector<float> subtsoslocalpitch2;
  vector<float> subclustersensorThickness;
  vector<float> subclustersensorThickness2;
  vector<float> subtsosBdotY;
  vector<float> subtsosBdotY2;
  vector<float> subclusterstripChargelocalpitch;
  vector<float> subclusterstripChargelocalpitch2;
  vector<float> subclusterstripChargesensorThickness;
  vector<float> subclusterstripChargesensorThickness2;
  vector<float> subclusterstripChargeBdotY;
  vector<float> subclusterstripChargeBdotY2;


  vector<float> subclusterstripChargeBdotYMC2;



  t1->SetBranchAddress("clustersubdetid",  &partition );
  t1->SetBranchAddress("clustercharge",  &clustercharge );
  t1->SetBranchAddress("clusterwidth",  &clusterwidth );
  t1->SetBranchAddress("clusterlayerwheel",  &clusterlayerwheel );
  t1->SetBranchAddress("clusterstripChargeSubdetid",  &clusterstripChargeSubdetid );
  t1->SetBranchAddress("clusterstripCharge",  &clusterstripCharge );
  t1->SetBranchAddress("clusterstripChargeLayerwheel",  &clusterstripChargeLayerwheel );
  t1->SetBranchAddress("tsostrackPt",  &tsostrackPt );
  t1->SetBranchAddress("tsoslocalphi",  &tsoslocalphi );
  t1->SetBranchAddress("tsoslocaltheta",  &tsoslocaltheta );
  t1->SetBranchAddress("tsosBdotY",  &tsosBdotY );
  t1->SetBranchAddress("clusterstripChargeStripNr",  &clusterstripChargeStripNr );
  t1->SetBranchAddress("clusterstripChargeTotWidth",  &clusterstripChargeTotWidth );
  t1->SetBranchAddress("clusterstripChargeTotCharge",  &clusterstripChargeTotCharge );
  t1->SetBranchAddress("clusterstripChargeLocalTrackPhi",  &clusterstripChargeLocalTrackPhi );
  t1->SetBranchAddress("clusterstripChargeGlobalTrackPhi",  &clusterstripChargeGlobalTrackPhi );
  t1->SetBranchAddress("clusterstripChargeLocalTrackTheta",  &clusterstripChargeLocalTrackTheta );
  t1->SetBranchAddress("clusterstripChargeGlobalTrackTheta",  &clusterstripChargeGlobalTrackTheta );
  t1->SetBranchAddress("clusterstripChargeDetid",  &clusterstripChargeDetid );
  t1->SetBranchAddress("clusterstripChargetrackPt",  &clusterstripChargetrackPt );


  uint32_t evCount=0;
   
  cout << "in here a" << endl;
  Int_t nentries = (Int_t)t1->GetEntries();
  cout << "entries " << nentries << endl;

  ///fill variables from tree 1
  for (Int_t e=0; e<nentries; e++) 
    {
      t1->GetEntry(e);
          
      //per cluster
      uint32_t up = partition->size();
      for(uint32_t k=0; k<up;k++)
	{
	  if( partition->at(k) == sdId )
	    {
	      if(clusterlayerwheel->at(k) == layer)
		{
		  subpartition.push_back(partition->at(k));
		  subclustercharge.push_back(clustercharge->at(k));
		  subclusterwidth.push_back(clusterwidth->at(k));
		  subclusterlayerwheel.push_back(clusterlayerwheel->at(k));
		  subtsoslocalphi.push_back(tsoslocalphi->at(k));
		  subtsoslocaltheta.push_back(tsoslocaltheta->at(k));
		  subtsosBdotY.push_back(tsosBdotY->at(k));
		}
	    }
	}
      //per strip
      uint32_t upStrip = clusterstripChargeSubdetid->size();
      std::cout<<"strip "<<upStrip<<" partition "<<up<<std::endl;
      for(uint32_t k=0; k<upStrip ; k++)
	{
	  if( clusterstripChargeSubdetid->at(k) == sdId)
	    {
	      if(clusterstripChargeLayerwheel->at(k)== layer)
		{
		  subclusterstripChargeSubdetid.push_back(clusterstripChargeSubdetid->at(k));
		  subclusterstripCharge.push_back(clusterstripCharge->at(k));
		  subclusterstripChargeLayerwheel.push_back(clusterstripChargeLayerwheel->at(k));
		  subclusterstripChargeLocalTrackPhi.push_back(clusterstripChargeLocalTrackPhi->at(k));
		  subclusterstripChargeLocalTrackTheta.push_back(clusterstripChargeLocalTrackTheta->at(k));
		  subclusterstripChargeGlobalTrackTheta.push_back(clusterstripChargeGlobalTrackTheta->at(k));
		  subclusterstripChargeTotWidth.push_back(clusterstripChargeTotWidth->at(k));
		  subclusterstripChargeStripNr.push_back(clusterstripChargeStripNr->at(k));
		  subclusterstripChargeTotCharge.push_back(clusterstripChargeTotCharge->at(k));

		}
	    }
	}
    }
  TH1F* chargeForAllWidthsData = new TH1F("chargeForAllWidths", "chargeForAllWidths" , 100, 0, 1000 );
  TH1F* clusterAllWidthsData = new TH1F("clusterAllWidths", "chargeForAllWidths" , 20, 0, 20 );
  TH1F* clusterAllWidthsfinal = new TH1F("clusterAllWidthsfinal", "chargeForAllWidthsfinal" , 20, 0, 20 );

  for(uint32_t m = 0; m<subclustercharge.size(); m++)
    {
      chargeForAllWidthsData->Fill(subclustercharge.at(m));
      clusterAllWidthsData->Fill(subclusterwidth.at(m));
      int factor = subtsosBdotY.at(m) > 0 ? 1 : -1;
      if(factor*tan(subtsoslocaltheta.at(m))> 0.78 && subclusterwidth.at(m)> 3){
	clusterAllWidthsfinal->Fill(subclusterwidth.at(m));

      }
    }
  clusterAllWidthsfinal->SetMarkerStyle(kFullCircle);
  TCanvas c0("chargeForAllWidths","chargeForAllWidths");
  clusterAllWidthsfinal->DrawNormalized("P");
  c0.SaveAs("plot2.eps");

  vector<float> clusterVector;
  TProfile* leadingStripProfileData = new TProfile("leadingcluschargeData", "leadingcluschargeData", 550, 0, 550, 0, 300);
  TH1D* leadingStripSumData = new TH1D("leadingclusSumData", "leadingclusSumData", 800, 0, 800);
  float lstripNr = 0; 
  float lstripCh = 0; 
  float lClusWidth = 0 ; 
  float lClusCharge = 0 ; 
  uint32_t stripCounter = 0; 
  bool clusterEnd = true;  
  vector<double> stripChargeSum;
  stripChargeSum.resize(512,0);

  TProfile* clusterStripCahrgeAsFceTanThetaData = new TProfile("clusterStripCahrgeAsFceTanThetaData", "clusterStripCahrgeAsFceTanThetaData" , 50, -4, 4, 0, 300 );

  TProfile* clusterShapeData = new TProfile("clusterShapeData", "clusterShapeData", 40, -20, 20, 0, 300);
  TProfile* clusterShape2Data = new TProfile("clusterShape2Data", "clusterShape2Data", 40, -20, 20, 0, 300);
  TProfile* clusterShape3Data = new TProfile("clusterShape3Data", "clusterShape3Data", 40, -20, 20, 0, 300);
  TProfile* clusterShape4Data = new TProfile("clusterShape4Data", "clusterShape4Data", 40, -20, 20, 0, 300);
  TProfile* clusterShape5Data = new TProfile("clusterShape5Data", "clusterShape5Data", 40, -20, 20, 0, 300);
  TProfile* clusterShape6Data = new TProfile("clusterShape6Data", "clusterShape6Data", 40, -20, 20, 0, 300);
  TProfile* clusterShape10Data = new TProfile("clusterShape10Data", "clusterShape10Data", 40, -20, 20, 0, 300);
  TProfile* clusterShape15Data = new TProfile("clusterShape15Data", "clusterShape15Data", 40, -20, 20, 0, 300);
  TProfile* clusterShape20Data = new TProfile("clusterShape20Data", "clusterShape20Data", 40, -20, 20, 0, 300);

  TProfile* clusterShapeDataPositive = new TProfile("clusterShapeDataPositive", "clusterShapeDataPositive", 40, 0, 20, 0, 300);
  TProfile* clusterShape2DataPositive = new TProfile("clusterShape2DataPositive", "clusterShape2DataPositive", 40, 0, 40, 0, 300);
  TProfile* clusterShape3DataPositive = new TProfile("clusterShape3DataPositive", "clusterShape3DataPositive", 40, 0, 40, 0, 300);
  TProfile* clusterShape4DataPositive = new TProfile("clusterShape4DataPositive", "clusterShape4DataPositive", 40, 0, 40, 0, 300);
  TProfile* clusterShape5DataPositive = new TProfile("clusterShape5DataPositive", "clusterShape5DataPositive", 40, 0, 40, 0, 300);
  TProfile* clusterShape6DataPositive = new TProfile("clusterShape6DataPositive", "clusterShape6DataPositive", 40, 0, 40, 0, 300);
  TProfile* clusterShape10DataPositive = new TProfile("clusterShape10DataPositive", "clusterShape10DataPositive", 40, 0, 40, 0, 300);
  TProfile* clusterShape15DataPositive = new TProfile("clusterShape15DataPositive", "clusterShape15DataPositive", 40, 0, 40, 0, 300);
  TProfile* clusterShape20DataPositive = new TProfile("clusterShape20DataPositive", "clusterShape20DataPositive", 40, 0, 40, 0, 300);


  TProfile* clusterShapeDataNegative = new TProfile("clusterShapeDataNegative", "clusterShapeDataNegative", 40, 0, 40, 0, 300);
  TProfile* clusterShape2DataNegative = new TProfile("clusterShape2DataNegative", "clusterShape2DataNegative", 40, 0, 40, 0, 300);
  TProfile* clusterShape3DataNegative = new TProfile("clusterShape3DataNegative", "clusterShape3DataNegative", 40, 0, 40, 0, 300);
  TProfile* clusterShape4DataNegative = new TProfile("clusterShape4DataNegative", "clusterShape4DataNegative", 40, 0, 40, 0, 300);
  TProfile* clusterShape5DataNegative = new TProfile("clusterShape5DataNegative", "clusterShape5DataNegative", 40, 0, 40, 0, 300);
  TProfile* clusterShape6DataNegative = new TProfile("clusterShape6DataNegative", "clusterShape6DataNegative", 40, 0, 40, 0, 300);
  TProfile* clusterShape10DataNegative = new TProfile("clusterShape10DataNegative", "clusterShape10DataNegative", 40, 0, 40, 0, 300);
  TProfile* clusterShape15DataNegative = new TProfile("clusterShape15DataNegative", "clusterShape15DataNegative", 40, 0, 40, 0, 300);
  TProfile* clusterShape20DataNegative = new TProfile("clusterShape20DataNegative", "clusterShape20DataNegative", 40, 0, 40, 0, 300);
  TProfile* clusterShape25DataNegative = new TProfile("clusterShape25DataNegative", "clusterShape25DataNegative", 40, 0, 40, 0, 300);


  TH1F* clusterChargeWOSaturationData = new TH1F("clusterChargeWOSaturationData", "clusterChargeWOSaturationData", 100, 0, 1000 );
  TH1F* clusterWidthWOSaturationData = new TH1F("clusterWidthWOSaturationData", "clusterWidthWOSaturationData", 15, 0, 15 );
  TH1F* clusterChargePerStripWOSaturationData = new TH1F("clusterChargePerStripWOSaturationData", "clusterChargePerStripWOSaturationData", 300, 0, 300 );
 
  cout<<"in line 254"<<endl;
  for(uint32_t m = 0; m<subclusterstripChargeLayerwheel.size(); m++)
    {
      if(subclusterstripChargeLayerwheel.at(m) == 3)
	{

	  if(lClusWidth==0 || clusterEnd==true )
	    {
	      lstripNr = subclusterstripChargeStripNr.at(m);
	      lstripCh = subclusterstripCharge.at(m);
	      lClusWidth = subclusterstripChargeTotWidth.at(m);
	      lClusCharge = subclusterstripChargeTotCharge.at(m);
	      stripCounter = 1; 
	      clusterVector.clear();
	    }
	  if(stripCounter <= lClusWidth)
	    {

	      clusterEnd = false;
	      if(lstripCh<subclusterstripCharge.at(m))
		{
		  lstripNr = subclusterstripChargeStripNr.at(m);
		  lstripCh = subclusterstripCharge.at(m);
		}
	      clusterVector.push_back(subclusterstripCharge.at(m)) ;
	    
	  if(stripCounter == lClusWidth)
	    {

	      leadingStripProfileData->Fill(lstripNr, lstripCh);
	      clusterStripCahrgeAsFceTanThetaData->Fill(tan(subclusterstripChargeLocalTrackTheta.at(m)), lstripCh);
	      int factor = subtsosBdotY.at(m) > 0 ? 1 : -1;
	      double prev = stripChargeSum.at(lstripNr);
	      stripChargeSum.at(lstripNr) = prev+lstripCh;
	      std::cout<<"line after factor "<<std::endl;
	      float maxChrg = 0;
	      int32_t maxChrgCtr = -1;
	      for(uint32_t ch=0; ch<clusterVector.size();ch++)
		{
		  if(maxChrg < clusterVector.at(ch))
		    {
		      maxChrg = clusterVector.at(ch);
		      maxChrgCtr = ch;
		    }
		}
	      if(maxChrg < 253)
		{
		  clusterWidthWOSaturationData->Fill(lClusWidth);
		  clusterChargeWOSaturationData->Fill(lClusCharge);
		}
	      for(uint32_t ch=0; ch<clusterVector.size();ch++)
		{
		  int32_t positionValue = ch-(maxChrgCtr);

		  clusterShapeData->Fill( positionValue ,clusterVector.at(ch) );
		  if( subclusterstripChargeLocalTrackTheta.at(m)>0.785 && subclusterstripChargeLocalTrackTheta.at(m)<1.57)
		    clusterShapeDataPositive->Fill( ch+1 ,clusterVector.at(ch) );
		  if( subclusterstripChargeLocalTrackTheta.at(m)>1.57 && subclusterstripChargeLocalTrackTheta.at(m)<2.355)
		    clusterShapeDataNegative->Fill( ch+1 ,clusterVector.at(ch) );

		  if(lClusWidth>3)
		    {
		      clusterShape2Data->Fill( positionValue ,clusterVector.at(ch) );
		      if( subclusterstripChargeLocalTrackTheta.at(m)>0.785 && subclusterstripChargeLocalTrackTheta.at(m)<1.57)
			clusterShape2DataPositive->Fill( ch+1 ,clusterVector.at(ch) );
		      if( subclusterstripChargeLocalTrackTheta.at(m)>1.57 && subclusterstripChargeLocalTrackTheta.at(m)<2.355)
			clusterShape2DataNegative->Fill( ch+1 ,clusterVector.at(ch) );
		    }
		  if(lClusWidth==3)
		    {
		      clusterShape3Data->Fill( positionValue ,clusterVector.at(ch) );
		      if( subclusterstripChargeLocalTrackTheta.at(m)>0.785 && subclusterstripChargeLocalTrackTheta.at(m)<1.57)
			clusterShape3DataPositive->Fill( ch+1 ,clusterVector.at(ch) );
		      if( subclusterstripChargeLocalTrackTheta.at(m)>1.57 && subclusterstripChargeLocalTrackTheta.at(m)<2.355)
			clusterShape3DataNegative->Fill( ch+1 ,clusterVector.at(ch) );
		    }
		  if(lClusWidth==4)
                    {
                      clusterShape4Data->Fill( positionValue ,clusterVector.at(ch) );
                      //if( subclusterstripChargeLocalTrackTheta.at(m)>0.785 && subclusterstripChargeLocalTrackTheta.at(m)<1.57)
		      if(factor*tan( subclusterstripChargeLocalTrackTheta.at(m))>1)
		      clusterShape4DataPositive->Fill( ch+1 ,clusterVector.at(ch) );
                      //if( subclusterstripChargeLocalTrackTheta.at(m)>1.57 && subclusterstripChargeLocalTrackTheta.at(m)<2.355)
		      if(factor*tan( subclusterstripChargeLocalTrackTheta.at(m))<-1)
                        clusterShape4DataNegative->Fill( ch+1 ,clusterVector.at(ch) );
                    }

		  if(lClusWidth==5)
                    {
                      clusterShape5Data->Fill( positionValue ,clusterVector.at(ch) );
                      if( subclusterstripChargeLocalTrackTheta.at(m)>0.785 && subclusterstripChargeLocalTrackTheta.at(m)<1.57)
                        clusterShape5DataPositive->Fill( ch+1 ,clusterVector.at(ch) );
                      if( subclusterstripChargeLocalTrackTheta.at(m)>1.57 && subclusterstripChargeLocalTrackTheta.at(m)<2.355)
                        clusterShape5DataNegative->Fill( ch+1 ,clusterVector.at(ch) );
                    }
		  
		  if(lClusWidth==6)
		    {
		  clusterShape6Data->Fill( positionValue ,clusterVector.at(ch) );
		  if( subclusterstripChargeLocalTrackTheta.at(m)>0.785 && subclusterstripChargeLocalTrackTheta.at(m)<1.57)
		    clusterShape6DataPositive->Fill( ch+1 ,clusterVector.at(ch) );
		  if( subclusterstripChargeLocalTrackTheta.at(m)>1.57 && subclusterstripChargeLocalTrackTheta.at(m)<2.355)
		    clusterShape6DataNegative->Fill( ch+1 ,clusterVector.at(ch) );
		    }

		  if(lClusWidth==10)
                    {
                      clusterShape10Data->Fill( positionValue ,clusterVector.at(ch) );
                      if( subclusterstripChargeLocalTrackTheta.at(m)>0.785 && subclusterstripChargeLocalTrackTheta.at(m)<1.57)
                        clusterShape10DataPositive->Fill( ch+1 ,clusterVector.at(ch) );
                      if( subclusterstripChargeLocalTrackTheta.at(m)>1.57 && subclusterstripChargeLocalTrackTheta.at(m)<2.355)
                        clusterShape10DataNegative->Fill( ch+1 ,clusterVector.at(ch) );
                    }
		  if(lClusWidth==15)
                    {
                      clusterShape15Data->Fill( positionValue ,clusterVector.at(ch) );
                      if( subclusterstripChargeLocalTrackTheta.at(m)>0.785 && subclusterstripChargeLocalTrackTheta.at(m)<1.57)
                        clusterShape15DataPositive->Fill( ch+1 ,clusterVector.at(ch) );
                      if( subclusterstripChargeLocalTrackTheta.at(m)>1.57 && subclusterstripChargeLocalTrackTheta.at(m)<2.355)
                        clusterShape15DataNegative->Fill( ch+1 ,clusterVector.at(ch) );
                    }

		  if(lClusWidth==20)
                    {
                      clusterShape20Data->Fill( positionValue ,clusterVector.at(ch) );
                      if( subclusterstripChargeLocalTrackTheta.at(m)>0.785 && subclusterstripChargeLocalTrackTheta.at(m)<1.57)
                        clusterShape20DataPositive->Fill( ch+1 ,clusterVector.at(ch) );
                      if( subclusterstripChargeLocalTrackTheta.at(m)>1.57 && subclusterstripChargeLocalTrackTheta.at(m)<2.355)
                        clusterShape20DataNegative->Fill( ch+1 ,clusterVector.at(ch) );
                    }

		}
	      clusterVector.clear();
	      clusterEnd =true;
	    }
	  stripCounter++;
	    }
	}
    
    }
  TFile* fout=TFile::Open("outputMC_327_6_new.root","RECREATE");
      clusterShape2Data->Write();
      clusterShape3Data->Write();
      clusterShape4Data->Write();
      clusterShape5Data->Write();
      clusterShape6Data->Write();
      clusterShape10Data->Write();
      clusterShape15Data->Write();
      clusterShape20Data->Write();
      clusterShape2DataPositive->Write();
      clusterShape3DataPositive->Write();
      clusterShape4DataPositive->Write();
      clusterShape5DataPositive->Write();
      clusterShape6DataPositive->Write();
      clusterShape10DataPositive->Write();
      clusterShape15DataPositive->Write();
      clusterShape20DataPositive->Write();
      clusterShape2DataNegative->Write();
      clusterShape3DataNegative->Write();
      clusterShape4DataNegative->Write();
      clusterShape5DataNegative->Write();
      clusterShape6DataNegative->Write();
      clusterShape10DataNegative->Write();
      clusterShape15DataNegative->Write();
      clusterShape20DataNegative->Write();
      return 0;
      
    }
