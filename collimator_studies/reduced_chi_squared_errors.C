#include <stdio.h>     
#include <stdlib.h>
#include "TMath.h"
#include <vector>
#include <utility>
#include <iostream>
#include <fstream>

// #include <libgen.h>

void reduced_chi_squared_errors(char *filename1=NULL,char *filename2=NULL) {
  /* A simple script to plot aspects of phototube hits 
   * 
   * I like to run this macro as 
   * $ root -l -x 'poisson_test.C("OD.root")'
   */

 // Open the WCSim file
  TFile *inFile1 = new TFile(filename1, "READ"); 
  if ( !inFile1->IsOpen() ){
    std::cout << "Error: could not open input file \"" << filename1 << "\"." <<std::endl; 
	
  }	
  TFile *inFile2 = new TFile(filename2, "READ"); 
  if ( !inFile2->IsOpen() ){
    std::cout << "Error: could not open input file \"" << filename2 << "\"." <<std::endl; 
	
  }	
  
  TTree *inTree1 = (TTree*)inFile1->Get("Events");

  Long64_t nentries1 = inTree1->GetEntries();

  std::cerr << "Entries " << nentries1 << std::endl;
  /*
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

  inTree1->SetBranchAddress("eventId",&eventId);
  inTree1->SetBranchAddress("tubeId",&tubeId);
  //inTree->SetBranchAddress("trueTime",&trueTime);
  inTree1->SetBranchAddress("recTime",&recTime);
  inTree1->SetBranchAddress("trueNumPhotons",&trueNumPhotons);
  inTree1->SetBranchAddress("recCharge",&recCharge);
  inTree1->SetBranchAddress("cylLoc",&cylLoc);
  inTree1->SetBranchAddress("pmtX",&pmtX);
  inTree1->SetBranchAddress("pmtY",&pmtY);
  inTree1->SetBranchAddress("pmtZ",&pmtZ);

*/
  TLeaf *eventid = inTree1->GetLeaf("eventId");
  eventid->GetBranch()->GetEntry(nentries1-1);
  int numevents = eventid->GetValue();

  //  std::cerr << "num eventids " << numevents << std::endl;

  int num_OD_PMT = 0;
  /*  TLeaf *tubeid = inTree->GetLeaf("tubeId");
      tubeid->GetBranch()->GetEntry(nentries-1);*/
  for (Long64_t i=0; i < nentries1; i++) {
    TLeaf *tubeid = inTree1->GetLeaf("tubeId");
    tubeid->GetBranch()->GetEntry(i);
    int temp_num_OD_PMT = tubeid->GetValue();
    if (temp_num_OD_PMT > num_OD_PMT) num_OD_PMT = temp_num_OD_PMT;
  }  

  TTree *inTree2 = (TTree*)inFile2->Get("Events");

  Long64_t nentries2 = inTree2->GetEntries();

  std::cerr << "Entries " << nentries2 << std::endl;
  /*
  Int_t eventId2;
  Int_t tubeId2;
  //Float_t trueTime;
  Float_t recTime2;
  Float_t trueNumPhotons2;
  Float_t recCharge2;
  Int_t cylLoc2;
  Float_t pmtX2;
  Float_t pmtY2;
  Float_t pmtZ2;

  inTree2->SetBranchAddress("eventId",&eventId2);
  inTree2->SetBranchAddress("tubeId",&tubeId2);
  //inTree->SetBranchAddress("trueTime",&trueTime);
  inTree2->SetBranchAddress("recTime",&recTime2);
  inTree2->SetBranchAddress("trueNumPhotons",&trueNumPhotons2);
  inTree2->SetBranchAddress("recCharge",&recCharge2);
  inTree2->SetBranchAddress("cylLoc",&cylLoc2);
  inTree2->SetBranchAddress("pmtX",&pmtX2);
  inTree2->SetBranchAddress("pmtY",&pmtY2);
  inTree2->SetBranchAddress("pmtZ",&pmtZ2);

*/
  //////////////////////////////////////////
  // HISTOGRAMS DEFINITION /////////////////
  //////////////////////////////////////////

  //  double num_OD_PMT = wcsimrootgeom1->GetODWCNumPMT();
  
  /*TH2D *QvsTubeID1 = new TH2D("QvsTubesID1","charge per PMT", num_OD_PMT,0,num_OD_PMT,);
  QvsTubeID1->SetXTitle("PMT ID");
  QvsTubeID1->SetYTitle("Q");

  TH1D *QvsTubeID2 = new TH1D("QvsTubesID2","charge per PMT", num_OD_PMT,0,num_OD_PMT,);
  QvsTubeID2->SetXTitle("PMT ID");
  QvsTubeID2->SetYTitle("Q");*/

  int num_elements = num_OD_PMT;
  int default_value = 0;

  std::vector<double> chi2values;
  double sum_chi2 = 0;
  double ndof = 0;
  std::vector<double> charge_vec_1(num_elements,default_value);
  std::vector<double> charge_vec_2(num_elements,default_value);
  
  // END HISTOGRAMS DEFINITION /////////////
  //////////////////////////////////////////

  int n1 = inTree1->Draw("eventId:tubeId:trueTime:recTime:trueNumPhotons:recCharge:cylLoc:pmtX:pmtY:pmtZ","","goff");
  int n2 = inTree2->Draw("eventId:tubeId:trueTime:recTime:trueNumPhotons:recCharge:cylLoc:pmtX:pmtY:pmtZ","","goff");

  int n_cycles = 10;
  int n1_per_cycle = n1/n_cycles;
  int n2_per_cycle = n2/n_cycles;

  std::vector<double> total_charge_vec_1(num_elements,default_value);
  std::vector<double> total_charge_vec_2(num_elements,default_value);

  for(long unsigned int cycle = 0; cycle < n_cycles; cycle++) {

    std::vector<double> charge_vec_1(num_elements,default_value);
    std::vector<double> charge_vec_2(num_elements,default_value);

    for (Long64_t i=(cycle*n1_per_cycle); i < (cycle+1)*n1_per_cycle; i++) {
      //
      inTree1->GetEntry(i);
      TLeaf *_cylLoc = inTree1->GetLeaf("cylLoc");
      _cylLoc->GetBranch()->GetEntry(i);
      int cylLoc = _cylLoc->GetValue();
      TLeaf *_tubeId = inTree1->GetLeaf("tubeId");
      _tubeId->GetBranch()->GetEntry(i);
      int tubeId = _tubeId->GetValue();
      TLeaf *_recCharge = inTree1->GetLeaf("recCharge");
      _recCharge->GetBranch()->GetEntry(i);
      double recCharge = _recCharge->GetValue();
      TLeaf *_pmtX = inTree1->GetLeaf("pmtX");
      _pmtX->GetBranch()->GetEntry(i);
      double pmtX = _pmtX->GetValue();
      TLeaf *_pmtY = inTree1->GetLeaf("pmtY");
      _pmtY->GetBranch()->GetEntry(i);
      double pmtY = _pmtX->GetValue();
      if (i%(n1/10) == 0) std::cerr << i << " out of " << n1 << " entries" << std::endl;

      if (cylLoc==4) {
	ndof = 6481;
	double pmtRad = sqrt((pow(pmtX, 2.0))+(pow(pmtY, 2.0)));
	double pmtTheta = std::atan(pmtY/pmtX);
	if (pmtX>0) pmtTheta += (std::acos(-1)/2);
	else pmtTheta -= (std::acos(-1)/2);
	double pmtRadTheta = pmtRad*pmtTheta;
	if (pmtRadTheta > 1000 && pmtRadTheta < 9000) {
	  //	  std::cerr << "z pos " << wcsimrootgeom1->GetPMT(tubeId-1).GetPosition(2) << std::endl;
	  //if (charge_vec_1[tubeId-1] == 0) charge_vec_1[tubeId-1] = testcharge;
	  //else {
	  double old_charge = charge_vec_1[tubeId-1];
	  double plus_new_charge = old_charge + recCharge;
	  charge_vec_1[tubeId-1] = plus_new_charge;
	  double old_total_charge = total_charge_vec_1[tubeId-1];
	  double total_plus_new_charge = old_total_charge + recCharge;
	  total_charge_vec_1[tubeId-1] = total_plus_new_charge;
	}
	//      if (plus_new_charge != old_charge) {
	//std::cerr << "CV1 charge " << plus_new_charge << std::endl;
	//}//}
      }
    }



    for (Long64_t i=(cycle*n2_per_cycle); i < (cycle+1)*n2_per_cycle; i++) {
      //
      inTree2->GetEntry(i);
      TLeaf *_cylLoc = inTree2->GetLeaf("cylLoc");
      _cylLoc->GetBranch()->GetEntry(i);
      int cylLoc = _cylLoc->GetValue();
      TLeaf *_tubeId = inTree2->GetLeaf("tubeId");
      _tubeId->GetBranch()->GetEntry(i);
      int tubeId = _tubeId->GetValue();
      TLeaf *_recCharge = inTree2->GetLeaf("recCharge");
      _recCharge->GetBranch()->GetEntry(i);
      double recCharge = _recCharge->GetValue();
      TLeaf *_pmtX = inTree2->GetLeaf("pmtX");
      _pmtX->GetBranch()->GetEntry(i);
      double pmtX = _pmtX->GetValue();
      TLeaf *_pmtY = inTree2->GetLeaf("pmtY");
      _pmtY->GetBranch()->GetEntry(i);
      double pmtY = _pmtX->GetValue();

      if (i%(n2/10) == 0) std::cerr << i << " out of " << n2 << " entries" << std::endl;

      if (cylLoc==4) {
	double pmtRad = sqrt((pow(pmtX, 2.0))+(pow(pmtY, 2.0)));
	double pmtTheta = std::atan(pmtY/pmtX);
	if (pmtX>0) pmtTheta += (std::acos(-1)/2);
	else pmtTheta -= (std::acos(-1)/2);
	double pmtRadTheta = pmtRad*pmtTheta;
	if (pmtRadTheta > 1000 && pmtRadTheta < 9000) {
	  //	  std::cerr << "z pos " << wcsimrootgeom1->GetPMT(tubeId-1).GetPosition(2) << std::endl;
	  //if (charge_vec_1[tubeId-1] == 0) charge_vec_1[tubeId-1] = testcharge;
	  //else {
	  double old_charge = charge_vec_2[tubeId-1];
	  double plus_new_charge = old_charge + recCharge;
	  //if (plus_new_charge != old_charge) {
	  //std::cerr << "CV2 charge " << plus_new_charge << std::endl;
	  //}
	  charge_vec_2[tubeId-1] = plus_new_charge;
	  double old_total_charge = total_charge_vec_2[tubeId-1];
	  double total_plus_new_charge = old_total_charge + recCharge;
	  total_charge_vec_2[tubeId-1] = total_plus_new_charge;
	}
	//}
      }
    }

    double PMT_count = 0;
    for (int i = 0; i < total_charge_vec_1.size(); i++) {
      if (total_charge_vec_1[i] > 0) PMT_count++;
    }

    ndof = PMT_count - 1;
    std::cerr << "ndof " << ndof << std::endl;

    //    std::cerr << "charge_vec_1.size() " << charge_vec_1.size() << ", charge_vec_2.size() " << charge_vec_2.size() << std::endl;

    double chi_squared = -99.0;

    for (int i = 0; i < charge_vec_1.size(); i++) {
      double observed = charge_vec_1[i];
      double expected = charge_vec_2[i];
      if (observed > 0 && expected > 0) {
	double chi_i = pow((observed-expected),2.0) / expected;
	//std::cerr << "chi_squred " << chi_squared << ", chi_i " << chi_i << ", expected " << expected << ", observed " << observed << std::endl;
	if (chi_squared == -99.0) chi_squared = chi_i;
	else chi_squared += chi_i;
      }
    }
    //    std::cerr << " chi_sq " << chi_squared << " ndof " << ndof << " chi/ndof " << chi_squared/ndof << std::endl;
    chi_squared /= ndof;

    sum_chi2 += chi_squared;
    chi2values.push_back(chi_squared);

  }

  double mean_chi2 = sum_chi2/chi2values.size();
  std::cerr << "mean " << mean_chi2 << std::endl;
  double var_chi2 = 0;
  double sigma_chi2 = -9999;

  for (int i = 0; i<chi2values.size(); i++) {
    std::cerr << "chi 2 value for " << i << ", " << chi2values.at(i) << ", var " << pow((chi2values.at(i)-mean_chi2),2) << std::endl;
    var_chi2 += pow((chi2values.at(i)-mean_chi2),2);
  }

  sigma_chi2 = sqrt(var_chi2/n_cycles);

  std::cerr << "Std dev: " << sigma_chi2 << std::endl;
  double relative_sigma = sigma_chi2/mean_chi2;

  double total_chi_squared = -99;

  for (int i = 0; i < total_charge_vec_1.size(); i++) {
    double observed = total_charge_vec_1[i];
    double expected = total_charge_vec_2[i];
    if (observed > 0 && expected > 0) {
      double chi_i = pow((observed-expected),2.0) / expected;
      //                  std::cerr << "chi_squred " << total_chi_squared << ", chi_i " << chi_i << ", expected " << expected << ", observed " << observed << std::endl;
      if (total_chi_squared == -99.0) total_chi_squared = chi_i;
      else total_chi_squared += chi_i;
    }
  }

  total_chi_squared /= ndof;

  double error_chi_squared = total_chi_squared*relative_sigma;

  std::cerr << "CHI2 " << total_chi_squared << ", error " << error_chi_squared << std::endl;

  ofstream myfile;
  myfile.open("cylLoc4_rayff_chi2_errors.txt", std::ios_base::app);
  //myfile << sigma_chi2 << ", ";
  myfile << error_chi_squared << ", ";
  myfile.close();


  ofstream myfile2;
  myfile2.open("rayff_cylLoc4_reduced_chi_sq.txt", std::ios_base::app);
  myfile2 << total_chi_squared << ", ";
  myfile2.close();


} // END MACRO
