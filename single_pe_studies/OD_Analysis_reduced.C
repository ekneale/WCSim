// C++ Includes
#include <iostream>
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <utility>


// ROOT Includes
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TEllipse.h"
// Defines
#define PI 3.141592654


/*
 * To run this script just type:  root -l -x 'OD_Analysis_reduced.C("wcsim.root", false)' 
 * or replace wcsim.root for your input filename.
 * The final true can be changed to true to turn on verbosity.
 */



// A function to convert radians to degress
float RadToDeg(float x){
  return x*180/PI;
}
// A function to convert degress to radians
float DegToRad(float x){
  return x*PI/180;
}

void OD_Analysis_reduced( const char *inFileName = "wcsim.root", bool verbosity = 0){ 
	

  // Some nicely formatted text options
  std::cout << std::scientific; // This causes all numbers to be displayed in scientific notation.
  std::cout << std::setprecision(2); // Sets the decimal precision (no more than two decimal places)
  std::cout << std::left; // Sets the text justification to left
  const int txtW = 20; // Width of "box" holding text
  const int numW = 10; // Width of "box" holding numbers


  // Open the WCSim file
  TFile *inFile = new TFile(inFileName, "READ"); 
  if ( !inFile->IsOpen() ){
    std::cout << "Error: could not open input file \"" << inFileName << "\"." <<std::endl; 
	
  } else if (verbosity) {
    std::cout << "Input file: " << inFileName << std::endl;
  }


  //-------------------- Reading in ----------------------

  //inFile->ls();
  
  // TTree *inTree = Events;  
  TTree *inTree = (TTree*)inFile->Get("Events");

  Long64_t nentries = inTree->GetEntries();

  std::cerr << "Entries " << nentries << std::endl;

  Int_t eventId;
  Int_t tubeId;
  //Float_t trueTime;
  Float_t recTime;
  Float_t trueNumPhotons;
  Float_t recCharge;
  Int_t cylLoc;
  Float_t pmtX;
  Float_t pmtY;
  Float_t pmtZ;
  Float_t nearestInjDist;

  inTree->SetBranchAddress("eventId",&eventId);
  inTree->SetBranchAddress("tubeId",&tubeId);
  //inTree->SetBranchAddress("trueTime",&trueTime);
  inTree->SetBranchAddress("recTime",&recTime);
  inTree->SetBranchAddress("trueNumPhotons",&trueNumPhotons);
  inTree->SetBranchAddress("recCharge",&recCharge);
  inTree->SetBranchAddress("cylLoc",&cylLoc);
  inTree->SetBranchAddress("pmtX",&pmtX);
  inTree->SetBranchAddress("pmtY",&pmtY);
  inTree->SetBranchAddress("pmtZ",&pmtZ);
  inTree->SetBranchAddress("nearestInjDist",&nearestInjDist);

  TLeaf *eventid = inTree->GetLeaf("eventId");
  eventid->GetBranch()->GetEntry(nentries-1);
  int numevents = eventid->GetValue();

  std::cerr << "num eventids " << numevents << std::endl;

  int num_OD_PMT = 0;
  /*  TLeaf *tubeid = inTree->GetLeaf("tubeId");
      tubeid->GetBranch()->GetEntry(nentries-1);*/
  for (Long64_t i=0; i < nentries; i++) {
    TLeaf *tubeid = inTree->GetLeaf("tubeId");
    tubeid->GetBranch()->GetEntry(i);
    int temp_num_OD_PMT = tubeid->GetValue();
    if (temp_num_OD_PMT > num_OD_PMT) num_OD_PMT = temp_num_OD_PMT;
  }  

  std::cerr << "num PMTs " << num_OD_PMT << std::endl;

  int nbPEMaxByPMT = 5;
  //double WLSLength = 48;
  //  double WLSLength = 30; //correct setting but harder to see in ODQpoly plot
    double WLSLength = 60;

  //  TF1 *poiss = new TF1("poiss","[0]*TMath::PoissonI(x,[1])",0,nbPEMaxByPMT);
  TF1 *poissI = new TF1("poiss","[0]*TMath::PoissonI(x,[1])",0,nbPEMaxByPMT);
  TF1 *poiss = new TF1("poiss","[0]*TMath::Poisson(x,[1])",0,nbPEMaxByPMT);
  TF1 *poiss1 = new TF1("poiss1","[0]*TMath::Poisson(x,[1])",0,nbPEMaxByPMT/4);
  TF1 *poiss2 = new TF1("poiss2","[0]*TMath::Poisson(x,[1])",0,nbPEMaxByPMT/8);

  //--------------------- Tree ready to use -------------------

  // Basic histograms
  
  TH1F** PPerPMT = new TH1F*[(int)num_OD_PMT];
  TH1F** QPerPMT = new TH1F*[(int)num_OD_PMT];
  for (int i=0; i<num_OD_PMT; i++) {
    PPerPMT[i] = new TH1F(Form("h%d",i),"test",(nbPEMaxByPMT+10),-10,nbPEMaxByPMT);
    QPerPMT[i] = new TH1F(Form("tq%d",i),"test1",5*(nbPEMaxByPMT+20),-10,nbPEMaxByPMT+10);
  }

  TH2F* poiss_mean_P_per_pmt = new TH2F("poissmeanP","poissmeanP", num_OD_PMT, 0, num_OD_PMT,nbPEMaxByPMT,0,nbPEMaxByPMT);
  TProfile2D* poiss_mean_P_per_pmt_vs_events = new TProfile2D("poissmeanPvsevents","poissmeanPvsevents", num_OD_PMT, 0, num_OD_PMT,nbPEMaxByPMT,0,nbPEMaxByPMT,0,1000);

  TH2F* poiss_mean_Q_per_pmt = new TH2F("poissQmean","poissQmean", num_OD_PMT, 0, num_OD_PMT,100,0,10);
  TProfile2D* poiss_mean_Q_per_pmt_vs_events = new TProfile2D("poissmeanQvsevents","poissmeanQvsevents", num_OD_PMT, 0, num_OD_PMT,(10*nbPEMaxByPMT),0,nbPEMaxByPMT,0,1000);

  TH2F* poiss_mean_Q_vs_dist = new TH2F("poissQmeandist","poissQmeandist",100,0,1000,(2*nbPEMaxByPMT),0,nbPEMaxByPMT);
  TProfile2D* poiss_mean_Q_vs_dist_vs_events = new TProfile2D("poissmeanQvseventsdist","poissmeanQvseventsdist", 100,0,1000,(2*nbPEMaxByPMT),0,nbPEMaxByPMT,0,1000);

  TH2F* QFilledEventsVSDist = new TH2F("qevtsvsdist","qevtsvsdist",500,0,1000,100,0,100);
  TH2F* QFilledEventsVSPMT = new TH2F("qevtsvspmt","qevtsvspmt",num_OD_PMT,0,num_OD_PMT,100,0,100);

  //  inTree->SetEstimate(inTree->GetEntries()+1);

  int n = inTree->Draw("eventId:tubeId:trueTime:recTime:trueNumPhotons:recCharge:cylLoc:pmtX:pmtY:pmtZ:nearestInjDist","","goff");

  std::vector<double> tube_inj_dists(num_OD_PMT);
  std::vector<double> xpos(num_OD_PMT);
  std::vector<double> ypos(num_OD_PMT);
  std::vector<double> zpos(num_OD_PMT);
  std::vector<double> cyl(num_OD_PMT);

  double barrelposx = -99999;
  double barrelposy = -99999;
  double toppos = -99999;

  for (Long64_t i=0; i < n; i++) {
    //
    inTree->GetEntry(i);

    if (cylLoc == 4) {
      barrelposx = pmtX;
      barrelposy = pmtY;
    }
    if (cylLoc == 5) {
      toppos = pmtZ;
    }
    //    std::cerr << tubeId << "  " << trueNumPhotons << std::endl;
    if (recCharge > 0) QPerPMT[tubeId-1]->Fill(recCharge);
    tube_inj_dists.at(tubeId-1) = nearestInjDist;
    xpos.at(tubeId-1) = pmtX;
    ypos.at(tubeId-1) = pmtY;
    zpos.at(tubeId-1) = pmtZ;
    cyl.at(tubeId-1) = cylLoc;
    //    if (recCharge > 1500) std::cerr << pmtX << ", " << pmtY << ", " << pmtZ << std::endl;
    
  }


  double RadiusOD = sqrt( pow(barrelposx,2) + pow(barrelposy,2) );
  double HeightOD = 2*(abs(toppos));

  TH2Poly *ODQpoly = new TH2Poly("ODQ","Charge in OD PMTs",-12000,12000,-12000,12000);
  TH2Poly *ODPEpoly = new TH2Poly("ODPE","Photons incident on OD PMTs",-12000,12000,-12000,12000);

  for (int polybin = 0; polybin < num_OD_PMT; polybin++) {
    //
    double tube[3];
    int cylLoc = cyl[polybin];
    tube[0] = xpos[polybin];
    tube[1] = ypos[polybin];
    tube[2] = zpos[polybin];

    if ( cylLoc == 5){
      ODQpoly->AddBin(tube[0] - WLSLength/2, tube[1] + RadiusOD + HeightOD/2 + 100 - WLSLength/2, tube[0] + WLSLength/2, tube[1] + RadiusOD + HeightOD/2 + 100 + WLSLength/2);
      ODPEpoly->AddBin(tube[0] - WLSLength/2, tube[1] + RadiusOD + HeightOD/2 + 100 - WLSLength/2, tube[0] + WLSLength/2, tube[1] + RadiusOD + HeightOD/2 + 100 + WLSLength/2);
    }
    //Bot OD
    else if ( cylLoc == 3){
      ODQpoly->AddBin(tube[0] - WLSLength/2,-(HeightOD/2 +RadiusOD +tube[1] + 100) - WLSLength/2, tube[0] + WLSLength/2,-(HeightOD/2 +RadiusOD +tube[1] + 100) + WLSLength/2);
      ODPEpoly->AddBin(tube[0] - WLSLength/2,-(HeightOD/2 +RadiusOD +tube[1] + 100) - WLSLength/2, tube[0] + WLSLength/2,-(HeightOD/2 +RadiusOD +tube[1] + 100) + WLSLength/2);
    }
    //Barrel OD
    else {

      double pmtRad = sqrt((pow(tube[0], 2.0))+(pow(tube[1], 2.0)));
      double pmtTheta = std::atan(tube[1]/tube[0]);
      if (tube[0]>0) pmtTheta += (std::acos(-1)/2);
      else pmtTheta -= (std::acos(-1)/2);
      double pmtRadTheta = pmtRad*pmtTheta;

      ODQpoly->AddBin(pmtRadTheta - WLSLength/2, tube[2] - WLSLength/2, pmtRadTheta + WLSLength/2, tube[2] + WLSLength/2);
      ODPEpoly->AddBin(pmtRadTheta - WLSLength/2, tube[2] - WLSLength/2, pmtRadTheta + WLSLength/2, tube[2] + WLSLength/2);
    }
  }

  for (Long64_t i=0; i < n; i++) {
    //
    inTree->GetEntry(i);

    double pmtRad = sqrt((pow(pmtX, 2.0))+(pow(pmtY, 2.0)));
    double pmtTheta = std::atan(pmtY/pmtX);
    double l = sqrt((pow((0-pmtX),2)) + pow((-(HeightOD/2)-pmtY),2));
    double angle = 2*asin(l/2*RadiusOD);
    double length = angle*RadiusOD;
    if(recCharge > 0) {
      if (cylLoc == 3) {
	double ODpolyX = pmtX;
	double ODpolyY = -(pmtY + RadiusOD + HeightOD/2 + 100);
	//	std::cerr << "cylLoc " << cylLoc << ", ODpolyX " << ODpolyX << ", ODpolyY " << ODpolyY << ", recCharge " << recCharge << std::endl;
	//	ODPEpoly->Fill(ODpolyX,ODpolyY,trueNumPhotons);
	ODQpoly->Fill(ODpolyX,ODpolyY,recCharge);
      }
      else if (cylLoc == 5) {
	double ODpolyX = pmtX;
	double ODpolyY = (HeightOD/2 +RadiusOD + pmtY + 100);
	//ODPEpoly->Fill(ODpolyX,ODpolyY,trueNumPhotons);
	ODQpoly->Fill(ODpolyX,ODpolyY,recCharge);
      }
      else if (cylLoc == 4) {
	if (pmtX>0) pmtTheta += (std::acos(-1)/2);
	else pmtTheta -= (std::acos(-1)/2);
	double pmtRadTheta = pmtRad*pmtTheta;
	//	ODPEpoly->Fill(pmtRadTheta,pmtZ,trueNumPhotons);
	ODQpoly->Fill(pmtRadTheta,pmtZ,recCharge);
      }
    }
   
  }

  for (Int_t i=0; i<num_OD_PMT; i++) {
    
    if (QPerPMT[i]->GetEntries() > 0) {

      QFilledEventsVSDist->Fill(tube_inj_dists.at(i),(QPerPMT[i]->GetEntries())/10);
      QFilledEventsVSPMT->Fill(i,(QPerPMT[i]->GetEntries())/10);
      poiss->SetParameter(0,QPerPMT[i]->GetBinContent(QPerPMT[i]->GetMaximumBin()));
      poiss->SetParameter(1,QPerPMT[i]->GetMean());
      TFitResultPtr r = QPerPMT[i]->Fit(poiss,"SRQ");
      Double_t par1 = r->Parameter(1);
      Double_t par0 = r->Parameter(0);
      //std::cerr << "distance " << tube_inj_dists.at(i) << ", fit mean : " << par1 << ", hist mean : " << QPerPMT[i]->GetMean() << std::endl;      
      //std::cerr << i+1 << " QPerPMT[i]->GetEntries() " << QPerPMT[i]->GetEntries() << ", QPerPMT[i]->GetBinContent(QPerPMT[i]->GetMaximumBin()) " << QPerPMT[i]->GetBinContent(QPerPMT[i]->GetMaximumBin()) << " , QPerPMT[i]->GetMean() " << QPerPMT[i]->GetMean() << ", fit mean " << par1 << std::endl;
      //std::cerr << i+1 << " normalisation " << par0 << std::endl;
      poiss_mean_Q_per_pmt->Fill(i+1,par1);
      poiss_mean_Q_per_pmt_vs_events->Fill(i+1,par1,QPerPMT[i]->GetEntries(),1);
      poiss_mean_Q_vs_dist->Fill(tube_inj_dists.at(i),par1);
      poiss_mean_Q_vs_dist_vs_events->Fill(tube_inj_dists.at(i),par1,QPerPMT[i]->GetEntries(),1);
  
    }
  }


  double num_OD_PMTs_with_hits = 0;
  for (Int_t i = 0; i < num_OD_PMT; i++) {
    if (QFilledEventsVSPMT->GetBinContent(i)) num_OD_PMTs_with_hits += 1;
  }

  std::cerr << "num_OD_PMTs_with_hits " << (num_OD_PMTs_with_hits) << ", num_OD_PMT " << num_OD_PMT << std::endl;
  std::cerr << "PMT COVERAGE : " << (num_OD_PMTs_with_hits/num_OD_PMT)*100 << std::endl;

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(kFALSE);
  gStyle->SetPalette(55);

  TCanvas *qevtPMT = new TCanvas("qevtPMT","qevtPMT",800,600);
  QFilledEventsVSPMT->SetXTitle("PMT ID");
  QFilledEventsVSPMT->SetYTitle("% Events with digitised charge in PMT");
  QFilledEventsVSPMT->SetMarkerStyle(6);
  QFilledEventsVSPMT->Draw();
  //qevtPMT->SaveAs("QEventCovPMTID.png");

  TCanvas *qevtdist = new TCanvas("qevtdist","qevtdist ",800,600);
  QFilledEventsVSDist->SetXTitle("Distance between PMT and nearest injector (cm)");
  QFilledEventsVSDist->SetYTitle("% Events with digitised charge in PMT");
  QFilledEventsVSDist->SetMarkerStyle(6);
  QFilledEventsVSDist->Draw();
  ///  qevtdist->SaveAs("80refl_QEventCovDist.png");

  TCanvas *poissmeanQvev = new TCanvas("pmQv","pmQv",800,600);
  poiss_mean_Q_per_pmt_vs_events->SetXTitle("PMT ID");
  poiss_mean_Q_per_pmt_vs_events->SetYTitle("Mean Digitised Charge");
  poiss_mean_Q_per_pmt_vs_events->SetMarkerStyle(6);
  poiss_mean_Q_per_pmt_vs_events->Draw("COLZ");
  //poissmeanQvev->SaveAs("80refl_MeanQPerPMTPerEvt.png");


  TCanvas *poissmeanQdistvev = new TCanvas("pmQdv","pmQdv",800,600);
  poiss_mean_Q_vs_dist_vs_events->GetXaxis()->SetRangeUser(0,850);
  //  poiss_mean_Q_vs_dist_vs_events->GetYaxis()->SetRangeUser(0,1000);
  poiss_mean_Q_vs_dist_vs_events->SetXTitle("Distance between PMT and nearest injector (cm)");
  poiss_mean_Q_vs_dist_vs_events->SetYTitle("Mean Digitised Charge");
  poiss_mean_Q_vs_dist_vs_events->SetMarkerStyle(6);
  poiss_mean_Q_vs_dist_vs_events->Draw("COLZ");
  TLine *line = new TLine(0,200,850,200);
  line->SetLineColor(kRed);
  line->Draw();
  //poissmeanQdistvev->SaveAs("feb_sat_diffuse_facing_out_MeanQVsDistPerEvt.png");

  //  TCanvas *cpq = new TCanvas("cpq","cpq",1500,1000);
  TCanvas *cpq = new TCanvas("cpq","cpq",1100,1000);
    ODQpoly->GetZaxis()->SetTitle("Digitised Charge");
    // ODQpoly->GetZaxis()->SetRangeUser(300.0,4000);
  // ODQpoly->GetXaxis()->SetRangeUser(4800.0,6800);
  //ODQpoly->GetYaxis()->SetRangeUser(-1000,1000);
  gPad->SetLogz(1);
  //cpq->SetFrameFillColor(0);
  //cpq->SetFrameFillStyle(0);
  //cpq->SetFrameLineColor(0);
  //cpq->SetFrameBorderMode(0);
  cpq->SetRightMargin(0.15);
  ODQpoly->Draw("COLZ 0");
  //cpq->SaveAs("current_diffuse_charge_map.png");

}

