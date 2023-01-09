// C++ Includes
#include <iostream>
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <utility>
#include <string>

// ROOT Includes
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TEllipse.h"
// Defines
#define PI 3.141592654



// A function to convert radians to degress
float RadToDeg(float x){
  return x*180/PI;
}
// A function to convert degress to radians
float DegToRad(float x){
  return x*PI/180;
}

void OD_Saturation_Analysis_reduced( const char *inFileName = "wcsim.root", const char *barrel_or_top = "barrel", const char *in_or_out = "out", const char *bare_or_diffuse = "diffuse",/* const char *outFileName = "calib_analysis.root",*/ int MaxPE = 2000, bool verbosity = 0){ 
	

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
	
  // Create an output file
  //TFile *outFile = new TFile(outFileName, "RECREATE");


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
  double max_charge = 0;
  /*  TLeaf *tubeid = inTree->GetLeaf("tubeId");
      tubeid->GetBranch()->GetEntry(nentries-1);*/
  for (Long64_t i=0; i < nentries; i++) {
    TLeaf *tubeid = inTree->GetLeaf("tubeId");
    tubeid->GetBranch()->GetEntry(i);
    int temp_num_OD_PMT = tubeid->GetValue();
    if (temp_num_OD_PMT > num_OD_PMT) num_OD_PMT = temp_num_OD_PMT;
    TLeaf *recCharge = inTree->GetLeaf("recCharge");
    recCharge->GetBranch()->GetEntry(i);
    double charge = recCharge->GetValue();
    if (charge > max_charge) max_charge = charge;
  }  

  std::cerr << "num PMTs " << num_OD_PMT << std::endl;
  std::cerr << "max charge " << max_charge << std::endl;

  //    int nbPEMaxByPMT = MaxPE;
  int nbPEMaxByPMT = max_charge*0.75;
  //double WLSLength = 48;
  double WLSLength = 30;

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

  TH2F* poiss_mean_Q_per_pmt = new TH2F("poissQmean","poissQmean", num_OD_PMT, 0, num_OD_PMT,nbPEMaxByPMT,0,nbPEMaxByPMT);
  TProfile2D* poiss_mean_Q_per_pmt_vs_events = new TProfile2D("poissmeanQvsevents","poissmeanQvsevents", num_OD_PMT, 0, num_OD_PMT,(nbPEMaxByPMT),0,nbPEMaxByPMT,0,1000);

  TH2F* poiss_mean_Q_vs_dist = new TH2F("poissQmeandist","poissQmeandist",100,0,1000,(2*nbPEMaxByPMT),0,nbPEMaxByPMT);
  TProfile2D* poiss_mean_Q_vs_dist_vs_events = new TProfile2D("poissmeanQvseventsdist","poissmeanQvseventsdist", 100,0,1000,nbPEMaxByPMT,0,nbPEMaxByPMT,0,1000);

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

      double l = sqrt( pow((0 - tube[0]),2) + pow((-RadiusOD - tube[1]),2));
      double angle = 2*asin(l/(2*RadiusOD));
      double length = angle*RadiusOD ;
      if (tube[0]<0) length *= -1;
      ODQpoly->AddBin(length - WLSLength/2, tube[2] - WLSLength/2, length + WLSLength/2, tube[2] + WLSLength/2);
      ODPEpoly->AddBin(length - WLSLength/2, tube[2] - WLSLength/2, length + WLSLength/2, tube[2] + WLSLength/2);
    }
  }

  double num_pmts_nearby = 0;
  double num_pmts_nearby_200pe = 0;
  double sum_charge_in_radius = 0;
  double err_sum_charge_in_radius = 0;
  double source_x = -9999;
  double source_y = -9999;
  double source_z = -9999;

  if (strcmp(barrel_or_top,"barrel")==0) {
    source_z = 0.0;
    if (strcmp(in_or_out,"in")==0) {
      source_x = 3384.69; // facing in
      source_y = 673.257; // facing in
    } else if (strcmp(in_or_out,"out")==0) {
      source_x = 3296.42; // facing out
      source_y = 655.7; // facing out
    } else {
      std::cerr << "You need to specifiy if the sources are facing in or out" << std::endl;
      return 0;
    }

    for (Int_t i=0; i<num_OD_PMT; i++) {
    
      if (QPerPMT[i]->GetEntries() > 0) {
	QFilledEventsVSDist->Fill(tube_inj_dists.at(i),(QPerPMT[i]->GetEntries())/10);
	QFilledEventsVSPMT->Fill(i,(QPerPMT[i]->GetEntries())/10);
	double pmtRad = sqrt((pow(xpos.at(i), 2.0))+(pow(ypos.at(i), 2.0)));
	double pmtTheta = std::atan(ypos.at(i)/xpos.at(i));
	double l = sqrt((pow((0-xpos.at(i)),2)) + pow((-(HeightOD/2)-ypos.at(i)),2));
	double angle = 2*asin(l/2*RadiusOD);
	double length = angle*RadiusOD;
      
	double distance = sqrt((source_x-xpos.at(i))*(source_x-xpos.at(i)) + (source_y-ypos.at(i))*(source_y-ypos.at(i)) + (source_z-zpos.at(i))*(source_z-zpos.at(i)));
	if (distance <= 150) {
	  poiss->SetParameter(0,QPerPMT[i]->GetBinContent(QPerPMT[i]->GetMaximumBin()));
	  poiss->SetParameter(1,QPerPMT[i]->GetMean());
	  TFitResultPtr r = QPerPMT[i]->Fit(poiss,"SRQ");
	  Double_t par1 = r->Parameter(1);
	  Double_t par0 = r->Parameter(0);
	  Double_t charge_err = r->Error(1);
	  //	  std::cerr << "QPerPMT[i]->GetBinContent(QPerPMT[i]->GetMaximumBin()) " << QPerPMT[i]->GetBinContent(QPerPMT[i]->GetMaximumBin()) << ", QPerPMT[i]->GetMean() " << QPerPMT[i]->GetMean() << ", par1 " << par1 << std::endl;
	  //	  std::cerr << "distance " << distance << ", tube_inj_dists.at(i) " << tube_inj_dists.at(i) << std::endl;
	  //std::cerr << "distance " << distance << ", fit mean : " << par1 << ", hist mean : " << QPerPMT[i]->GetMean() << std::endl;
	  //std::cerr << i+1 << " QPerPMT[i]->GetEntries() " << QPerPMT[i]->GetEntries() << ", QPerPMT[i]->GetBinContent(QPerPMT[i]->GetMaximumBin()) " << QPerPMT[i]->GetBinContent(QPerPMT[i]->GetMaximumBin()) << " , QPerPMT[i]->GetMean() " << QPerPMT[i]->GetMean() << ", fit mean " << par1 << std::endl;
	  //std::cerr << i+1 << " normalisation " << par0 << std::endl;
	  poiss_mean_Q_per_pmt->Fill(i+1,par1);
	  poiss_mean_Q_per_pmt_vs_events->Fill(i+1,par1,QPerPMT[i]->GetEntries(),1);
	  //	  poiss_mean_Q_vs_dist->Fill(tube_inj_dists.at(i),par1);
	  poiss_mean_Q_vs_dist->Fill(distance,par1);
	  //poiss_mean_Q_vs_dist_vs_events->Fill(tube_inj_dists.at(i),par1,QPerPMT[i]->GetEntries(),1);
	  poiss_mean_Q_vs_dist_vs_events->Fill(distance,par1,QPerPMT[i]->GetEntries(),1);
	  if (xpos.at(i)>0) pmtTheta += (std::acos(-1)/2);
	  else pmtTheta -= (std::acos(-1)/2);
	  double pmtRadTheta = pmtRad*pmtTheta;
	  ODQpoly->Fill(pmtRadTheta,zpos.at(i),par1);
	  num_pmts_nearby += 1;
	  if (par1 > 199.9) num_pmts_nearby_200pe += 1;
	  sum_charge_in_radius += par1;
	  err_sum_charge_in_radius += charge_err*charge_err;
	}
	else if (distance > 150 && distance <= 450) {
	  poiss1->SetParameter(0,QPerPMT[i]->GetBinContent(QPerPMT[i]->GetMaximumBin()));
	  poiss1->SetParameter(1,QPerPMT[i]->GetMean());
	  TFitResultPtr r = QPerPMT[i]->Fit(poiss1,"SRQ");
	  Double_t par1 = r->Parameter(1);
	  Double_t par0 = r->Parameter(0);
	  Double_t charge_err = r->Error(1);      
	  //std::cerr << i+1 << " QPerPMT[i]->GetEntries() " << QPerPMT[i]->GetEntries() << ", QPerPMT[i]->GetBinContent(QPerPMT[i]->GetMaximumBin()) " << QPerPMT[i]->GetBinContent(QPerPMT[i]->GetMaximumBin()) << " , QPerPMT[i]->GetMean() " << QPerPMT[i]->GetMean() << ", fit mean " << par1 << std::endl;
	  //std::cerr << i+1 << " normalisation " << par0 << std::endl;
	  poiss_mean_Q_per_pmt->Fill(i+1,par1);
	  poiss_mean_Q_per_pmt_vs_events->Fill(i+1,par1,QPerPMT[i]->GetEntries(),1);
	  //poiss_mean_Q_vs_dist->Fill(tube_inj_dists.at(i),par1);
	  //poiss_mean_Q_vs_dist_vs_events->Fill(tube_inj_dists.at(i),par1,QPerPMT[i]->GetEntries(),1);
	  poiss_mean_Q_vs_dist->Fill(distance,par1);
	  poiss_mean_Q_vs_dist_vs_events->Fill(distance,par1,QPerPMT[i]->GetEntries(),1);
	  if (xpos.at(i)>0) pmtTheta += (std::acos(-1)/2);
	  else pmtTheta -= (std::acos(-1)/2);
	  double pmtRadTheta = pmtRad*pmtTheta;
	  ODQpoly->Fill(pmtRadTheta,zpos.at(i),par1);
	  num_pmts_nearby += 1;
	  if (par1 > 199.9) num_pmts_nearby_200pe += 1;
	  sum_charge_in_radius += par1;
	  err_sum_charge_in_radius += charge_err*charge_err;
	}
	else if (distance > 450 && distance < 850) { 
	  poiss2->SetParameter(0,QPerPMT[i]->GetBinContent(QPerPMT[i]->GetMaximumBin()));
	  poiss2->SetParameter(1,QPerPMT[i]->GetMean());
	  TFitResultPtr r = QPerPMT[i]->Fit(poiss2,"SRQ");
	  Double_t par1 = r->Parameter(1);
	  Double_t par0 = r->Parameter(0);
	  Double_t charge_err = r->Error(1);      
      
	  //std::cerr << i+1 << " QPerPMT[i]->GetEntries() " << QPerPMT[i]->GetEntries() << ", QPerPMT[i]->GetBinContent(QPerPMT[i]->GetMaximumBin()) " << QPerPMT[i]->GetBinContent(QPerPMT[i]->GetMaximumBin()) << " , QPerPMT[i]->GetMean() " << QPerPMT[i]->GetMean() << ", fit mean " << par1 << std::endl;
	  //std::cerr << i+1 << " normalisation " << par0 << std::endl;
	  poiss_mean_Q_per_pmt->Fill(i+1,par1);
	  poiss_mean_Q_per_pmt_vs_events->Fill(i+1,par1,QPerPMT[i]->GetEntries(),1);
	  //poiss_mean_Q_vs_dist->Fill(tube_inj_dists.at(i),par1);
	  //poiss_mean_Q_vs_dist_vs_events->Fill(tube_inj_dists.at(i),par1,QPerPMT[i]->GetEntries(),1);
	  poiss_mean_Q_vs_dist->Fill(distance,par1);
	  poiss_mean_Q_vs_dist_vs_events->Fill(distance,par1,QPerPMT[i]->GetEntries(),1);
	  if (xpos.at(i)>0) pmtTheta += (std::acos(-1)/2);
	  else pmtTheta -= (std::acos(-1)/2);
	  double pmtRadTheta = pmtRad*pmtTheta;
	  ODQpoly->Fill(pmtRadTheta,zpos.at(i),par1);
	  num_pmts_nearby += 1;
	  if (par1 > 199.9) num_pmts_nearby_200pe += 1;
	  sum_charge_in_radius += par1;
	  err_sum_charge_in_radius += charge_err*charge_err;
	}
      }
    }
  } else if (strcmp(barrel_or_top,"top")==0) {
    source_x = 800*(std::cos(3*(std::acos(-1)/12)));
    source_y = 800*(std::sin(3*(std::acos(-1)/12)));
    if (strcmp(in_or_out,"in")==0) {
      source_z = 3650;
    } else if (strcmp(in_or_out,"out")==0) {
      source_z = 3460;
    } else {
      std::cerr << "You need to specifiy if the sources are facing in or out" << std::endl;
      return 0;
    }
    for (Int_t i=0; i<num_OD_PMT; i++) {
    
      if (QPerPMT[i]->GetEntries() > 0) {
	double distance = sqrt((source_x-xpos.at(i))*(source_x-xpos.at(i)) + (source_y-ypos.at(i))*(source_y-ypos.at(i)) + (source_z-zpos.at(i))*(source_z-zpos.at(i)));
	QFilledEventsVSDist->Fill(distance,(QPerPMT[i]->GetEntries())/10);
	
	//	std::cerr << "check QFilledEventsVSPMT->Fill, i " << i << " (QPerPMT[i]->GetEntries())/10 " << (QPerPMT[i]->GetEntries())/10 << std::endl;
	QFilledEventsVSPMT->Fill(i,(QPerPMT[i]->GetEntries())/10);
      
	if (distance <= 150) {
	  poiss->SetParameter(0,QPerPMT[i]->GetBinContent(QPerPMT[i]->GetMaximumBin()));
	  poiss->SetParameter(1,QPerPMT[i]->GetMean());
	  TFitResultPtr r = QPerPMT[i]->Fit(poiss,"SRQ");
	  Double_t par1 = r->Parameter(1);
	  Double_t par0 = r->Parameter(0);
	  Double_t charge_err = r->Error(1);
      
	  poiss_mean_Q_per_pmt->Fill(i+1,par1);
	  poiss_mean_Q_per_pmt_vs_events->Fill(i+1,par1,QPerPMT[i]->GetEntries(),1);
	  //poiss_mean_Q_vs_dist->Fill(tube_inj_dists.at(i),par1);
	  //poiss_mean_Q_vs_dist_vs_events->Fill(tube_inj_dists.at(i),par1,QPerPMT[i]->GetEntries(),1);
	  poiss_mean_Q_vs_dist->Fill(distance,par1);
          poiss_mean_Q_vs_dist_vs_events->Fill(distance,par1,QPerPMT[i]->GetEntries(),1);
	  double ODpolyX = xpos.at(i);
	  double ODpolyY = (HeightOD/2 +RadiusOD + ypos.at(i) + 100);
	  ODQpoly->Fill(ODpolyX,ODpolyY,par1);
	  num_pmts_nearby += 1;
	  if (par1 > 199.9) num_pmts_nearby_200pe += 1;
	  sum_charge_in_radius += par1;
	  err_sum_charge_in_radius += charge_err*charge_err;
	}
	else if (distance > 150 && distance <= 450) {
	  poiss1->SetParameter(0,QPerPMT[i]->GetBinContent(QPerPMT[i]->GetMaximumBin()));
	  poiss1->SetParameter(1,QPerPMT[i]->GetMean());
	  TFitResultPtr r = QPerPMT[i]->Fit(poiss1,"SRQ");
	  Double_t par1 = r->Parameter(1);
	  Double_t par0 = r->Parameter(0);
	  Double_t charge_err = r->Error(1);
      
	  //std::cerr << i+1 << " QPerPMT[i]->GetEntries() " << QPerPMT[i]->GetEntries() << ", QPerPMT[i]->GetBinContent(QPerPMT[i]->GetMaximumBin()) " << QPerPMT[i]->GetBinContent(QPerPMT[i]->GetMaximumBin()) << " , QPerPMT[i]->GetMean() " << QPerPMT[i]->GetMean() << ", fit mean " << par1 << std::endl;
	  //std::cerr << i+1 << " normalisation " << par0 << std::endl;
	  poiss_mean_Q_per_pmt->Fill(i+1,par1);
	  poiss_mean_Q_per_pmt_vs_events->Fill(i+1,par1,QPerPMT[i]->GetEntries(),1);
	  //poiss_mean_Q_vs_dist->Fill(tube_inj_dists.at(i),par1);
	  //poiss_mean_Q_vs_dist_vs_events->Fill(tube_inj_dists.at(i),par1,QPerPMT[i]->GetEntries(),1);
	  poiss_mean_Q_vs_dist->Fill(distance,par1);
          poiss_mean_Q_vs_dist_vs_events->Fill(distance,par1,QPerPMT[i]->GetEntries(),1);
	  double ODpolyX = xpos.at(i);
	  double ODpolyY = (HeightOD/2 +RadiusOD + ypos.at(i) + 100);
	  ODQpoly->Fill(ODpolyX,ODpolyY,par1);
	  num_pmts_nearby += 1;
	  if (par1 > 199.9) num_pmts_nearby_200pe += 1;
	  sum_charge_in_radius += par1;
	  err_sum_charge_in_radius += charge_err*charge_err;
	}
	else if (distance > 450 && distance < 850) { 
	  poiss2->SetParameter(0,QPerPMT[i]->GetBinContent(QPerPMT[i]->GetMaximumBin()));
	  poiss2->SetParameter(1,QPerPMT[i]->GetMean());
	  TFitResultPtr r = QPerPMT[i]->Fit(poiss2,"SRQ");
	  Double_t par1 = r->Parameter(1);
	  Double_t par0 = r->Parameter(0);
	  Double_t charge_err = r->Error(1);
      
	  //std::cerr << i+1 << " QPerPMT[i]->GetEntries() " << QPerPMT[i]->GetEntries() << ", QPerPMT[i]->GetBinContent(QPerPMT[i]->GetMaximumBin()) " << QPerPMT[i]->GetBinContent(QPerPMT[i]->GetMaximumBin()) << " , QPerPMT[i]->GetMean() " << QPerPMT[i]->GetMean() << ", fit mean " << par1 << std::endl;
	  //std::cerr << i+1 << " normalisation " << par0 << std::endl;
	  poiss_mean_Q_per_pmt->Fill(i+1,par1);
	  poiss_mean_Q_per_pmt_vs_events->Fill(i+1,par1,QPerPMT[i]->GetEntries(),1);
	  //poiss_mean_Q_vs_dist->Fill(tube_inj_dists.at(i),par1);
	  //poiss_mean_Q_vs_dist_vs_events->Fill(tube_inj_dists.at(i),par1,QPerPMT[i]->GetEntries(),1);
	  poiss_mean_Q_vs_dist->Fill(distance,par1);
          poiss_mean_Q_vs_dist_vs_events->Fill(distance,par1,QPerPMT[i]->GetEntries(),1);
	  double ODpolyX = xpos.at(i);
	  double ODpolyY = (HeightOD/2 +RadiusOD + ypos.at(i) + 100);
	  ODQpoly->Fill(ODpolyX,ODpolyY,par1);
	  num_pmts_nearby += 1;
	  if (par1 > 199.9) num_pmts_nearby_200pe += 1;
	  sum_charge_in_radius += par1;
	  err_sum_charge_in_radius += charge_err*charge_err;
	}
      }
    }
  } else {
    std::cerr << "You need to specifiy if barrel or top cap source" << std::endl;
    return 0;
  }

  //  std::cerr << "QFilledEventsVSDist->GetEntries() " << QFilledEventsVSDist->GetEntries() << std::endl;

  double percentage_PMT_200pe_plus = (num_pmts_nearby_200pe/num_pmts_nearby)*100;
  std::cerr << "% of PMTs within 8.5m with 200+ pe " << percentage_PMT_200pe_plus << std::endl;

  double num_OD_PMTs_with_hits = 0;
  for (Int_t i = 0; i < num_OD_PMT; i++) {
    if (QFilledEventsVSPMT->GetBinContent(i)) num_OD_PMTs_with_hits += 1;
  }

  std::cerr << "Mean charge deposited within 8.5m of source " << sum_charge_in_radius << " error : " << sqrt(err_sum_charge_in_radius) << std::endl;

  std::cerr << "num_OD_PMTs_with_hits " << (num_OD_PMTs_with_hits) << ", num_OD_PMT " << num_OD_PMT << std::endl;
  std::cerr << "PMT COVERAGE : " << (num_OD_PMTs_with_hits/num_OD_PMT)*100 << std::endl;

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(kFALSE);
  gStyle->SetPalette(55);
  /*
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
  */

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
  
  //  poiss_mean_Q_vs_dist_vs_events->Write();
  //  poissmeanQdistvev->SaveAs("feb_sat_diffuse_facing_out_MeanQVsDistPerEvt.png");

  //  TCanvas *cpq = new TCanvas("cpq","cpq",1500,1000);
  TCanvas *cpq = new TCanvas("cpq","cpq",1100,1000);
  ODQpoly->GetZaxis()->SetRangeUser(10.0,nbPEMaxByPMT);
  if (strcmp(barrel_or_top,"barrel")==0) {
    ODQpoly->GetXaxis()->SetRangeUser(4800.0,6800);
    ODQpoly->GetYaxis()->SetRangeUser(-1000,1000);
  } else if (strcmp(barrel_or_top,"top")==0) {
    ODQpoly->GetXaxis()->SetRangeUser(-565,1565);
    ODQpoly->GetYaxis()->SetRangeUser((-565 + HeightOD/2 + RadiusOD + 100),(1565 + HeightOD/2 + RadiusOD + 100));
  }
  gPad->SetLogz(1);
  //cpq->SetFrameFillColor(0);
  //cpq->SetFrameFillStyle(0);
  //cpq->SetFrameLineColor(0);
  //cpq->SetFrameBorderMode(0);
  cpq->SetRightMargin(0.15);
  //  ODQpoly->Draw("COLZ TEXTx0");

  //  double sourceRad = sqrt((pow(source_x, 2.0))+(pow(source_y, 2.0)));
  double sourceRad = sqrt((pow(3296.42, 2.0))+(pow(655.7, 2.0)));
  double sourceTheta = std::atan(source_y/source_x) + (std::acos(-1)/2); 
  double sourceRadTheta = sourceRad*sourceTheta;
        
  double source_theta = -999;
  if (strcmp(bare_or_diffuse,"bare")==0) {
    source_theta = (12*std::acos(-1))/180;
  } else if (strcmp(bare_or_diffuse,"diffuse")==0) {
    source_theta = (40*std::acos(-1))/180;
  }

  double refl_1 = 9999;
  double refl_2 = 9999;
  double refl_3 = 9999;
  double refl_4 = 9999;
  
  double centre_x = -9999;
  double centre_y = -9999;
  if (strcmp(barrel_or_top,"top")==0) {
    centre_x = source_x;
    centre_y = HeightOD/2 +RadiusOD + source_y + 100;
    if (strcmp(in_or_out,"in")==0) {
      refl_1 = 2*(std::tan(source_theta));
    } else if (strcmp(in_or_out,"out")==0) {
      refl_1 = ((0.9 + 2)*std::tan(source_theta));
    }
    refl_2 = (refl_1 + 4*std::tan(source_theta));
    refl_3 = (refl_2 + 4*std::tan(source_theta));
    refl_4 = (refl_3 + 4*std::tan(source_theta));
  } else if (strcmp(barrel_or_top,"barrel")==0) {
    centre_x = sourceRadTheta;
    centre_y = 0;
    if (strcmp(in_or_out,"in")==0) {
      refl_1 = (std::tan(source_theta));
    } else if (strcmp(in_or_out,"out")==0) {
      refl_1 = ((0.9 + 1)*std::tan(source_theta));
    }
    refl_2 = (refl_1 + 2*std::tan(source_theta));
    refl_3 = (refl_2 + 2*std::tan(source_theta));
    refl_4 = (refl_3 + 2*std::tan(source_theta));
  }

  TEllipse *el1 = new TEllipse(centre_x, centre_y,refl_1*100,refl_1*100);
  TEllipse *el2 = new TEllipse(centre_x, centre_y,refl_2*100,refl_2*100);
  TEllipse *el3 = new TEllipse(centre_x, centre_y,refl_3*100,refl_3*100);
  TEllipse *el4 = new TEllipse(centre_x, centre_y,refl_4*100,refl_4*100);

  //std::cerr << "centrex " << centre_x << " rad1 " << refl_1*100 << " rad2 " << refl_2*100 << " rad3 " << refl_3*100 << " rad4 " << refl_4*100 <<std::endl;
  

  ODQpoly->Draw("COLZ0");
  el1->SetFillStyle(0);
  el1->SetLineColor(7);
  el1->SetLineWidth(3);
  el1->Draw();
  el2->SetFillStyle(0);
  el2->SetLineColor(6);
  el2->SetLineWidth(3);
  el2->Draw();
  el3->SetFillStyle(0);
  el3->SetLineColor(3);
  el3->SetLineWidth(3);
  el3->Draw();
  el4->SetFillStyle(0);
  el4->SetLineColor(2);
  el4->SetLineWidth(3);
  if (strcmp(barrel_or_top,"barrel")==0) el4->Draw();
  //  ODQpoly->Write();
  //  cpq->SaveAs("feb_sat_diffuse_facing_out_beamspot.pdf");
  poissmeanQdistvev->Update();
  
  //  outFile->Close();

}

