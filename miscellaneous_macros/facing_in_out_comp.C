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

void facing_in_out_comp( const char *inFileName1 = "facing_out.root", const char *inFileName2 = "facing_in.root", bool verbosity = 0){ 
	

  // Some nicely formatted text options
  std::cout << std::scientific; // This causes all numbers to be displayed in scientific notation.
  std::cout << std::setprecision(2); // Sets the decimal precision (no more than two decimal places)
  std::cout << std::left; // Sets the text justification to left
  const int txtW = 20; // Width of "box" holding text
  const int numW = 10; // Width of "box" holding numbers


  // Open the WCSim file
  TFile *inFile1 = new TFile(inFileName1, "READ"); 
  if ( !inFile1->IsOpen() ){
    std::cout << "Error: could not open input file \"" << inFileName1 << "\"." <<std::endl; 
	
  } else if (verbosity) {
    std::cout << "Input file: " << inFileName1 << std::endl;
  }

  TFile *inFile2 = new TFile(inFileName2, "READ"); 
  if ( !inFile2->IsOpen() ){
    std::cout << "Error: could not open input file \"" << inFileName2 << "\"." <<std::endl; 
	
  } else if (verbosity) {
    std::cout << "Input file: " << inFileName2 << std::endl;
  }

  std::vector<Double_t> distance1;
  std::vector<Double_t> distance2;
  std::vector<Double_t> mean1;
  std::vector<Double_t> mean2;
  std::vector<Double_t> mean_err1;
  std::vector<Double_t> mean_err2;

	
  for (int i = 0; i < 2; i++) {
    if (i==0) {
      TTree *inTree = (TTree*)inFile1->Get("Events");

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

      int num_OD_PMT = 0;
      for (Long64_t i=0; i < nentries; i++) {
	TLeaf *tubeid = inTree->GetLeaf("tubeId");
	tubeid->GetBranch()->GetEntry(i);
	int temp_num_OD_PMT = tubeid->GetValue();
	if (temp_num_OD_PMT > num_OD_PMT) num_OD_PMT = temp_num_OD_PMT;
      }  

      std::cerr << "num PMTs " << num_OD_PMT << std::endl;

      int nbPEMaxByPMT = 2000;
      TF1 *poiss = new TF1("poiss","[0]*TMath::Poisson(x,[1])",0,nbPEMaxByPMT);
      TF1 *poiss1 = new TF1("poiss1","[0]*TMath::Poisson(x,[1])",0,nbPEMaxByPMT/4);
      TF1 *poiss2 = new TF1("poiss2","[0]*TMath::Poisson(x,[1])",0,nbPEMaxByPMT/8);

      //--------------------- Tree ready to use -------------------

      // Basic histograms
  
      TH1F** QPerPMT = new TH1F*[(int)num_OD_PMT];
      for (int i=0; i<num_OD_PMT; i++) {
	QPerPMT[i] = new TH1F(Form("tq%d",i),"test1",5*(nbPEMaxByPMT+20),-10,nbPEMaxByPMT+10);
      }


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


      double source_x = 3296.42;
      double source_y = 655.7;
      double source_z = 0.0;

      for (Int_t i=0; i<num_OD_PMT; i++) {
    
	if (QPerPMT[i]->GetEntries() > 0) {
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
	    Double_t err_par1 = r->ParError(1);
	    distance1.push_back(distance);
	    mean1.push_back(par1);
	    mean_err1.push_back(err_par1);
	    //	  std::cerr << "QPerPMT[i]->GetBinContent(QPerPMT[i]->GetMaximumBin()) " << QPerPMT[i]->GetBinContent(QPerPMT[i]->GetMaximumBin()) << ", QPerPMT[i]->GetMean() " << QPerPMT[i]->GetMean() << ", par1 " << par1 << std::endl;
	    //	  std::cerr << "distance " << distance << ", tube_inj_dists.at(i) " << tube_inj_dists.at(i) << std::endl;
	    //std::cerr << "distance " << distance << ", fit mean : " << par1 << ", hist mean : " << QPerPMT[i]->GetMean() << std::endl;
	    //std::cerr << i+1 << " QPerPMT[i]->GetEntries() " << QPerPMT[i]->GetEntries() << ", QPerPMT[i]->GetBinContent(QPerPMT[i]->GetMaximumBin()) " << QPerPMT[i]->GetBinContent(QPerPMT[i]->GetMaximumBin()) << " , QPerPMT[i]->GetMean() " << QPerPMT[i]->GetMean() << ", fit mean " << par1 << std::endl;
	    //std::cerr << i+1 << " normalisation " << par0 << std::endl;
	  
	  }
	  else if (distance > 150 && distance <= 450) {
	    poiss1->SetParameter(0,QPerPMT[i]->GetBinContent(QPerPMT[i]->GetMaximumBin()));
	    poiss1->SetParameter(1,QPerPMT[i]->GetMean());
	    TFitResultPtr r = QPerPMT[i]->Fit(poiss1,"SRQ");
	    Double_t par1 = r->Parameter(1);
	    Double_t par0 = r->Parameter(0);
	    Double_t err_par1 = r->ParError(1);
	    distance1.push_back(distance);
	    mean1.push_back(par1);
	    mean_err1.push_back(err_par1);
	    //std::cerr << i+1 << " QPerPMT[i]->GetEntries() " << QPerPMT[i]->GetEntries() << ", QPerPMT[i]->GetBinContent(QPerPMT[i]->GetMaximumBin()) " << QPerPMT[i]->GetBinContent(QPerPMT[i]->GetMaximumBin()) << " , QPerPMT[i]->GetMean() " << QPerPMT[i]->GetMean() << ", fit mean " << par1 << std::endl;
	    //std::cerr << i+1 << " normalisation " << par0 << std::endl;

	  }
	  else if (distance > 450 && distance < 850) { 
	    poiss2->SetParameter(0,QPerPMT[i]->GetBinContent(QPerPMT[i]->GetMaximumBin()));
	    poiss2->SetParameter(1,QPerPMT[i]->GetMean());
	    TFitResultPtr r = QPerPMT[i]->Fit(poiss2,"SRQ");
	    Double_t par1 = r->Parameter(1);
	    Double_t par0 = r->Parameter(0);
	    Double_t err_par1 = r->ParError(1);
	    distance1.push_back(distance);
	    mean1.push_back(par1);
	    mean_err1.push_back(err_par1);
	  }
	}
      }
    } else if (i==1) {
      TTree *inTree = (TTree*)inFile2->Get("Events");

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

      int num_OD_PMT = 0;
      for (Long64_t i=0; i < nentries; i++) {
	TLeaf *tubeid = inTree->GetLeaf("tubeId");
	tubeid->GetBranch()->GetEntry(i);
	int temp_num_OD_PMT = tubeid->GetValue();
	if (temp_num_OD_PMT > num_OD_PMT) num_OD_PMT = temp_num_OD_PMT;
      }  

      std::cerr << "num PMTs " << num_OD_PMT << std::endl;

      int nbPEMaxByPMT = 2000;
      TF1 *poiss = new TF1("poiss","[0]*TMath::Poisson(x,[1])",0,nbPEMaxByPMT);
      TF1 *poiss1 = new TF1("poiss1","[0]*TMath::Poisson(x,[1])",0,nbPEMaxByPMT/4);
      TF1 *poiss2 = new TF1("poiss2","[0]*TMath::Poisson(x,[1])",0,nbPEMaxByPMT/8);

      //--------------------- Tree ready to use -------------------

      // Basic histograms
  
      TH1F** QPerPMT = new TH1F*[(int)num_OD_PMT];
      for (int i=0; i<num_OD_PMT; i++) {
	QPerPMT[i] = new TH1F(Form("tq%d",i),"test1",5*(nbPEMaxByPMT+20),-10,nbPEMaxByPMT+10);
      }


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


      double source_x = 3296.42;
      double source_y = 655.7;
      double source_z = 0.0;

      for (Int_t i=0; i<num_OD_PMT; i++) {
    
	if (QPerPMT[i]->GetEntries() > 0) {
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
	    Double_t err_par1 = r->ParError(1);
	    distance2.push_back(distance);
	    mean2.push_back(par1);
	    mean_err2.push_back(err_par1);
	    //	  std::cerr << "QPerPMT[i]->GetBinContent(QPerPMT[i]->GetMaximumBin()) " << QPerPMT[i]->GetBinContent(QPerPMT[i]->GetMaximumBin()) << ", QPerPMT[i]->GetMean() " << QPerPMT[i]->GetMean() << ", par1 " << par1 << std::endl;
	    //	  std::cerr << "distance " << distance << ", tube_inj_dists.at(i) " << tube_inj_dists.at(i) << std::endl;
	    //std::cerr << "distance " << distance << ", fit mean : " << par1 << ", hist mean : " << QPerPMT[i]->GetMean() << std::endl;
	    //std::cerr << i+1 << " QPerPMT[i]->GetEntries() " << QPerPMT[i]->GetEntries() << ", QPerPMT[i]->GetBinContent(QPerPMT[i]->GetMaximumBin()) " << QPerPMT[i]->GetBinContent(QPerPMT[i]->GetMaximumBin()) << " , QPerPMT[i]->GetMean() " << QPerPMT[i]->GetMean() << ", fit mean " << par1 << std::endl;
	    //std::cerr << i+1 << " normalisation " << par0 << std::endl;
	  
	  }
	  else if (distance > 150 && distance <= 450) {
	    poiss1->SetParameter(0,QPerPMT[i]->GetBinContent(QPerPMT[i]->GetMaximumBin()));
	    poiss1->SetParameter(1,QPerPMT[i]->GetMean());
	    TFitResultPtr r = QPerPMT[i]->Fit(poiss1,"SRQ");
	    Double_t par1 = r->Parameter(1);
	    Double_t par0 = r->Parameter(0);
	    Double_t err_par1 = r->ParError(1);
	    distance2.push_back(distance);
	    mean2.push_back(par1);
	    mean_err2.push_back(err_par1);
	    //std::cerr << i+1 << " QPerPMT[i]->GetEntries() " << QPerPMT[i]->GetEntries() << ", QPerPMT[i]->GetBinContent(QPerPMT[i]->GetMaximumBin()) " << QPerPMT[i]->GetBinContent(QPerPMT[i]->GetMaximumBin()) << " , QPerPMT[i]->GetMean() " << QPerPMT[i]->GetMean() << ", fit mean " << par1 << std::endl;
	    //std::cerr << i+1 << " normalisation " << par0 << std::endl;

	  }
	  else if (distance > 450 && distance < 850) { 
	    poiss2->SetParameter(0,QPerPMT[i]->GetBinContent(QPerPMT[i]->GetMaximumBin()));
	    poiss2->SetParameter(1,QPerPMT[i]->GetMean());
	    TFitResultPtr r = QPerPMT[i]->Fit(poiss2,"SRQ");
	    Double_t par1 = r->Parameter(1);
	    Double_t par0 = r->Parameter(0);
	    Double_t err_par1 = r->ParError(1);
	    distance2.push_back(distance);
	    mean2.push_back(par1);
	    mean_err2.push_back(err_par1);
	  }
	}
      }
    }
  }

  std::cerr << "distance1.size() " << distance1.size() << " distance2.size() " << distance2.size() << std::endl;

  /*  for (int i = 0; i < distance1.size(); i++) {
    std::cerr << "dist " << distance1.at(i) << ", mean " << mean1.at(i) << ", mean_err " << mean_err1.at(i) << std::endl;
  }
  */
  TGraphErrors *ge1 = new TGraphErrors(distance1.size(), &(distance1[0]), &(mean1[0]), 0, &(mean_err1[0]) );
  TGraphErrors *ge2 = new TGraphErrors(distance2.size(), &(distance2[0]), &(mean2[0]), 0, &(mean_err2[0]) );
  
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(kFALSE);
  gStyle->SetPalette(55);

  std::vector<Double_t> ratio_fout_fin;
  std::vector<Double_t> err_ratio_fout_fin;
  for (int i = 0; i < distance1.size(); i++) {
    Double_t facing_out_mean = ge1->GetPointY(i);
    Double_t facing_in_mean = ge2->GetPointY(i);
    Double_t err_facing_out_mean = ge1->GetErrorY(i);
    Double_t err_facing_in_mean = ge2->GetErrorY(i);
    std::cerr << "facing_in_mean err_facing_in_mean " << facing_in_mean << " " << err_facing_in_mean << std::endl;
    Double_t ratio = facing_out_mean/facing_in_mean;
    Double_t comb_error = ratio*(sqrt((err_facing_out_mean/facing_out_mean)*(err_facing_out_mean/facing_out_mean) + (err_facing_in_mean/facing_in_mean)*(err_facing_in_mean/facing_in_mean)));
    ratio_fout_fin.push_back(ratio);
    err_ratio_fout_fin.push_back(comb_error);
  }

  TGraphErrors *ratio_plot = new TGraphErrors(distance1.size(), &(distance1[0]), &(ratio_fout_fin[0]), 0, &(err_ratio_fout_fin[0]) );
  TCanvas *c3 = new TCanvas("c3","Mean charge facing out/mean charge facing in",800,600);
  ratio_plot->SetMarkerColor(4);
  ratio_plot->SetMarkerStyle(7);
  ratio_plot->GetXaxis()->SetRangeUser(0,850);
  ratio_plot->GetYaxis()->SetTitle("Q(source facing out) / Q(source facing in)");
  ratio_plot->GetXaxis()->SetTitle("Distance between PMT and injector (cm)");
  ratio_plot->Draw("AP");
  TLine *line = new TLine(0,1,850,1);
  line->SetLineColor(14);
  line->SetLineStyle(8);
  line->Draw();

}

