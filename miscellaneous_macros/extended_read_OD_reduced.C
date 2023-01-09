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

// Defines
#define PI 3.141592654


/*
 * To run this script just type:  root -l -x llib.C 'ODAnalysis.C("wcsim.root", "ODHits.root", false)' 
 * or replace wcsim.root and ODHits.root for your input and output filenames respectively
 * The final true can be changed to true to turn on verbosity.
 * The default: root -l -x llib.C ODAnalysis.C will run with the variables mentioned above
 */



// A function to convert radians to degress
float RadToDeg(float x){
  return x*180/PI;
}
// A function to convert degress to radians
float DegToRad(float x){
  return x*PI/180;
}

void extended_read_OD_reduced( const char *inFileName = "wcsim.root", bool verbosity = 0){ 
	

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
  
  inTree->SetBranchAddress("eventId",&eventId);
  inTree->SetBranchAddress("tubeId",&tubeId);
  inTree->SetBranchAddress("recTime",&recTime);
  inTree->SetBranchAddress("recCharge",&recCharge);
  inTree->SetBranchAddress("cylLoc",&cylLoc);
  inTree->SetBranchAddress("pmtX",&pmtX);
  inTree->SetBranchAddress("pmtY",&pmtY);
  inTree->SetBranchAddress("pmtZ",&pmtZ);
  
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

  int nbPEMaxByPMT = 25;
  //double WLSLength = 48;
  double WLSLength = 60;

  //  inTree->SetEstimate(inTree->GetEntries()+1);

  int n = inTree->Draw("eventId:tubeId:recTime:recCharge:cylLoc:pmtX:pmtY:pmtZ","","goff");

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
    if (cylLoc == 3) {
      toppos = pmtZ;
    }
    //    std::cerr << tubeId << "  " << trueNumPhotons << std::endl;
    xpos.at(tubeId-1) = pmtX;
    ypos.at(tubeId-1) = pmtY;
    zpos.at(tubeId-1) = pmtZ;
    cyl.at(tubeId-1) = cylLoc;

    //if (cylLoc == 3) std::cerr << "pmtX " << pmtX << ", pmtY " << pmtY << ", pmtZ " << pmtZ << std::endl;
  }

  // ---------------------------- Light coverage plots

  
  TH2Poly *ODQpoly = new TH2Poly("ODQ","Charge in OD PMTs",-12000,12000,-12000,12000);

  double RadiusOD = sqrt( pow(barrelposx,2) + pow(barrelposy,2) );
  double HeightOD = 2*(abs(toppos));

  std::cerr << "radius " << RadiusOD << " height " << HeightOD << std::endl;

  for (int polybin = 0; polybin < num_OD_PMT; polybin++) {
    //
    double tube[3];
    int cylLoc = cyl[polybin];
    tube[0] = xpos[polybin];
    tube[1] = ypos[polybin];
    tube[2] = zpos[polybin];

    if ( cylLoc == 5){
      ODQpoly->AddBin(tube[0] - WLSLength/2, tube[1] + RadiusOD + HeightOD/2 + 100 - WLSLength/2, tube[0] + WLSLength/2, tube[1] + RadiusOD + HeightOD/2 + 100 + WLSLength/2);
    }
    //Bot OD
    else if ( cylLoc == 3){
      //      std::cerr << "tube[0] - WLSLength/2 " << tube[0] - WLSLength/2 << ", -(HeightOD/2 +RadiusOD +tube[1] + 100) - WLSLength/2" << -(HeightOD/2 +RadiusOD +tube[1] + 100) - WLSLength/2 << ", tube[0] + WLSLength/2 " << tube[0] + WLSLength/2 << ", -(HeightOD/2 +RadiusOD +tube[1] + 100) + WLSLength/2" << -(HeightOD/2 +RadiusOD +tube[1] + 100) + WLSLength/2 << std::endl;
      ODQpoly->AddBin(tube[0] - WLSLength/2,-(HeightOD/2 +RadiusOD +tube[1] + 100) - WLSLength/2, tube[0] + WLSLength/2,-(HeightOD/2 +RadiusOD +tube[1] + 100) + WLSLength/2);
    }
    //Barrel OD
    else {
      double pmtRad = sqrt((pow(tube[0], 2.0))+(pow(tube[1], 2.0)));
      double pmtTheta = std::atan(tube[1]/tube[0]);
      if (tube[0]>0) pmtTheta += (std::acos(-1)/2);
      else pmtTheta -= (std::acos(-1)/2);
      double pmtRadTheta = pmtRad*pmtTheta;

      /*      double l = sqrt( pow((0 - tube[0]),2) + pow((-RadiusOD - tube[1]),2));
      double angle = 2*asin(l/(2*RadiusOD));
      double length = angle*RadiusOD ;
      if (tube[0]<0) length *= -1;*/
      ODQpoly->AddBin(pmtRadTheta - WLSLength/2, tube[2] - WLSLength/2, pmtRadTheta + WLSLength/2, tube[2] + WLSLength/2);
    }
  }

  // need another loop through branches to fill light cov plots

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
	//	continue;
	double ODpolyX = pmtX;
	double ODpolyY = (HeightOD/2 +RadiusOD + pmtY + 100);
	//ODPEpoly->Fill(ODpolyX,ODpolyY,trueNumPhotons);
	ODQpoly->Fill(ODpolyX,ODpolyY,recCharge);
      }
      else if (cylLoc == 4) {
	//continue;
	if (pmtX>0) pmtTheta += (std::acos(-1)/2);
	else pmtTheta -= (std::acos(-1)/2);
	double pmtRadTheta = pmtRad*pmtTheta;
	//	ODPEpoly->Fill(pmtRadTheta,pmtZ,trueNumPhotons);
	ODQpoly->Fill(pmtRadTheta,pmtZ,recCharge);
      }
    }
    /*    if(recCharge > 0) {
      if (cylLoc == 3) {
	double ODpolyX = pmtX;
	double ODpolyY = pmtY + RadiusOD + HeightOD/2 + 100;
	ODQpoly->Fill(ODpolyX,ODpolyY,recCharge);
      }
      else if (cylLoc == 5) {
	double ODpolyX = pmtX;
	double ODpolyY = -(HeightOD/2 +RadiusOD + pmtY + 100);
	ODQpoly->Fill(ODpolyX,ODpolyY,recCharge);
      }
      else if (cylLoc == 4) {
	if (pmtX>0) pmtTheta += (std::acos(-1)/2);
	else pmtTheta -= (std::acos(-1)/2);
	double pmtRadTheta = pmtRad*pmtTheta;
	ODQpoly->Fill(pmtRadTheta,pmtZ,recCharge);
      }
      }*/
    
  }



  gStyle->SetOptStat(0);
  gStyle -> SetOptTitle(kFALSE);
  gStyle->SetPalette(55);

  TCanvas *cpq = new TCanvas("cpq","cpq",650,600);
  ODQpoly->GetZaxis()->SetRangeUser(200,4000);
  ODQpoly->GetZaxis()->SetTitle("Digitised Charge");
  //  ODQpoly->GetXaxis()->SetRangeUser(-3500,3500);
  //   ODQpoly->GetYaxis()->SetRangeUser(-3500, 3500);
  //  ODQpoly->GetYaxis()->SetRangeUser(-3500, 3500);
  gPad->SetLogz(1);
    cpq->SetFrameFillColor(0);
  cpq->SetFrameFillStyle(0);
  cpq->SetFrameLineColor(0);
  cpq->SetFrameBorderMode(0);
  cpq->SetRightMargin(0.15);
  //ODQpoly->Draw("COLZ 0 AH");
  ODQpoly->Draw("COLZ 0");
  // cpq->SaveAs("diffuse_10k_QCoverage.png");


}

