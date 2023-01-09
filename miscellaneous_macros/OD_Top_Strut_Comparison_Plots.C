// C++ Includes
#include <iostream>
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <utility>
#include <string>
#include <map>
// ROOT Includes
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TEllipse.h"
// Defines
#define PI 3.141592654

std::vector<double> Reduce_Vector(std::vector<double> input_vector, int factor) {

  int count = 0;
  int new_index = 0;
  std::vector<double> reduced_vector;

  for (int i = 0; i < input_vector.size(); i++) {
    if (count == 0) {
      reduced_vector.push_back(input_vector.at(i));
      count++;
    } else if (count < (factor-1)) {
      double temp_vector = reduced_vector.at(new_index) + input_vector.at(i);
      reduced_vector.at(new_index) = temp_vector;
      count++;
    } else if (count == (factor-1)) {
      double temp_vector = reduced_vector.at(new_index) + input_vector.at(i);
      double average_vector = temp_vector/factor;
      reduced_vector.at(new_index) = average_vector;
      new_index++;
      count = 0;
    }
  }
 return reduced_vector;
}

std::vector<double> Reduce_Error_Vector(std::vector<double> input_vector, int factor) {

  int count = 0;
  int new_index = 0;
  std::vector<double> reduced_vector;
  
  for (int i = 0; i < input_vector.size(); i++) {
    if (count == 0) {
      reduced_vector.push_back((input_vector.at(i)*input_vector.at(i)));
      count++;
    } else if (count < (factor-1)) {
      double temp_vector = reduced_vector.at(new_index) + (input_vector.at(i)*input_vector.at(i));
      reduced_vector.at(new_index) = temp_vector;
      count++;
    } else if (count == (factor-1)) {
      double temp_vector = reduced_vector.at(new_index) + (input_vector.at(i)*input_vector.at(i));
      double average_vector = temp_vector;
      reduced_vector.at(new_index) = average_vector;
      new_index++;
      count = 0;
    }
  }
  return reduced_vector;
}

void Sort_And_Reduce_Vectors(std::vector<double> distance_vec, std::vector<double> mean_charge_vec, std::vector<double> eff_vec, std::vector<double> mean_charge_err_vec, std::vector<double> eff_err_vec, std::vector<double>* reduced_distance_vec, std::vector<double>* reduced_mean_charge_vec, std::vector<double>* reduced_eff_vec, std::vector<double>* reduced_mean_charge_err_vec, std::vector<double>* reduced_eff_err_vec) {

  std::vector< std::pair <double,double> > dist_charge_vec;
  std::vector< std::pair <double,double> > dist_charge_err_vec;
  std::vector< std::pair <double,double> > dist_eff_vec;
  std::vector< std::pair <double,double> > dist_eff_err_vec;

  for(int i=0; i < distance_vec.size(); ++i) {
    double d = distance_vec.at(i);
    double q = mean_charge_vec.at(i);
    double e = eff_vec.at(i);
    double q_err = mean_charge_err_vec.at(i);
    double e_err = eff_err_vec.at(i);
    dist_charge_vec.push_back(std::make_pair(d,q));
    dist_charge_err_vec.push_back(std::make_pair(d,q_err));
    dist_eff_vec.push_back(std::make_pair(d,e));
    dist_eff_err_vec.push_back(std::make_pair(d,e_err));
  }
 
  // using modified sort() function to sort
  std::sort(dist_charge_vec.begin(), dist_charge_vec.end());

  std::sort(dist_charge_err_vec.begin(), dist_charge_err_vec.end());

  std::sort(dist_eff_vec.begin(), dist_eff_vec.end());

  std::sort(dist_eff_err_vec.begin(), dist_eff_err_vec.end());

  // Now order separate vectors
  std::vector<double> distance;
  std::vector<double> charge;
  std::vector<double> charge_error;
  std::vector<double> efficiency;
  std::vector<double> efficiency_error;

  for (int i=0; i<dist_charge_vec.size(); i++) {
    distance.push_back(dist_charge_vec[i].first);
    charge.push_back(dist_charge_vec[i].second);
    charge_error.push_back(dist_charge_err_vec[i].second);
    efficiency.push_back(dist_eff_vec[i].second);
    efficiency_error.push_back(dist_eff_err_vec[i].second);
  }

  //Now reduce the vectors
  *reduced_distance_vec = Reduce_Vector(distance, 10);
  *reduced_mean_charge_vec = Reduce_Vector(charge, 10);
  *reduced_mean_charge_err_vec = Reduce_Error_Vector(charge_error, 10);
  *reduced_eff_vec = Reduce_Vector(efficiency, 10);
  *reduced_eff_err_vec = Reduce_Error_Vector(efficiency_error, 10);
}


Double_t PMT_Injector_Distance(double pmt[3]) {
    //------------------------------------------------------------- INJ POS

  std::vector<std::vector<double>> toppos;

  // Injector positions
  double maxr = 3317.01;
  double rad1 = 800;
  double rad2 = 1810;
  double rad3 = 2800;

  // top cap
  for (int i = 0; i < 3; i++) {
    double zpos = 3362.01;
    double xpos = 0.;
    double ypos = 0.;
    if (i==0) {
      for (int angfrac1 = 0; angfrac1 < 3; angfrac1++) {
	std::vector<double> single_top_pos (3);
	xpos = rad1*(std::cos((angfrac1*8 + 3)*(std::acos(-1)/12)));
	ypos = rad1*(std::sin((angfrac1*8 + 3)*(std::acos(-1)/12)));
	single_top_pos.at(0) = xpos;
	single_top_pos.at(1) = ypos;
	single_top_pos.at(2) = zpos;
	toppos.push_back(single_top_pos);

      }
    }
    if (i==1) {
      for (int angfrac2 = 0; angfrac2 < 6; angfrac2++) {
	std::vector<double> single_top_pos (3);
	xpos = rad2*(std::cos((angfrac2*4 + 1)*(std::acos(-1)/12)));
	ypos = rad2*(std::sin((angfrac2*4 + 1)*(std::acos(-1)/12)));
	single_top_pos.at(0) = xpos;
	single_top_pos.at(1) = ypos;
	single_top_pos.at(2) = zpos;
	toppos.push_back(single_top_pos);
      }
    }
    if (i==2) {
      for (int angfrac3 = 0; angfrac3 < 12; angfrac3++) {
	std::vector<double> single_top_pos (3);
	xpos = rad3*(std::cos((angfrac3)*(std::acos(-1)/6)));
	ypos = rad3*(std::sin((angfrac3)*(std::acos(-1)/6)));
	single_top_pos.at(0) = xpos;
	single_top_pos.at(1) = ypos;
	single_top_pos.at(2) = zpos;
	toppos.push_back(single_top_pos);
      }
    }
  }
  
  //------------------------------------------------------------- INJ POS

  Double_t distance = 999999;
  for (int i = 0; i < toppos.size(); i++) {
    std::vector<double> temp_inj = toppos.at(i);
    double injx = temp_inj[0];
    double injy = temp_inj[1];
    double injz = temp_inj[2];

    //std::cerr << " IX " << injx << " IY " << injy << " IZ " << injz << std::endl;
    //std::cerr << " PX " << pmt[0] << " PY " << pmt[1] << " PZ " << pmt[2] << std::endl;
    
    Double_t xsquared = pow((pmt[0] - injx),2.0);
    Double_t ysquared = pow((pmt[1] - injy),2.0);
    Double_t zsquared = pow((pmt[2] - injz),2.0);
    Double_t temp_dist = sqrt(xsquared + ysquared + zsquared);
    if (temp_dist < distance) distance = temp_dist;
  }
  return distance;   
}

void Fill_Q_Eff_Vectors(std::vector<double>* distance_vec, std::vector<double>* mean_charge_vec, std::vector<double>* eff_vec, std::vector<double>* mean_charge_err_vec, std::vector<double>* eff_err_vec, int num_PMTs, std::vector<int> cyl, std::vector<TH1F*> QPerPMT, int nbPEMaxByPMT, std::vector<double> tube_inj_dists) {

  TF1 *poissI = new TF1("poiss","[0]*TMath::PoissonI(x,[1])",0,nbPEMaxByPMT);
  TF1 *poiss = new TF1("poiss","[0]*TMath::Poisson(x,[1])",0,nbPEMaxByPMT);
  TF1 *poiss1 = new TF1("poiss1","[0]*TMath::Poisson(x,[1])",0,nbPEMaxByPMT/4);
  TF1 *poiss2 = new TF1("poiss2","[0]*TMath::Poisson(x,[1])",0,nbPEMaxByPMT/8);

  for (Int_t i=0; i<num_PMTs; i++) {
    
    if (cyl.at(i) == 5 && QPerPMT[i]->GetEntries() > 0) {
      double distance = tube_inj_dists.at(i);
      distance_vec->push_back(distance);
      // Efficiency assumes we had 1000 events , make that better later
      // Also don't make these percentages, it is screwing you over
      double efficiency = (QPerPMT[i]->GetEntries()/1000);
      // error on counts is sqrt of counts
      // efficiency is counts*(1/1000), error on efficiency is err_counts*(1/1000)
      double efficiency_error = sqrt(QPerPMT[i]->GetEntries())/1000;
      eff_vec->push_back(efficiency);
      eff_err_vec->push_back(efficiency_error);
      if (distance <= 150) {
	poiss->SetParameter(0,QPerPMT[i]->GetBinContent(QPerPMT[i]->GetMaximumBin()));
	poiss->SetParameter(1,QPerPMT[i]->GetMean());
	TFitResultPtr r = QPerPMT[i]->Fit(poiss,"SRQ");
	Double_t par1 = r->Parameter(1);
	Double_t par0 = r->Parameter(0);
	Double_t err_par1 = r->ParError(1);

	mean_charge_vec->push_back(par1);
	mean_charge_err_vec->push_back(err_par1);
      }
      else if (distance > 150 && distance <= 450) {
	poiss1->SetParameter(0,QPerPMT[i]->GetBinContent(QPerPMT[i]->GetMaximumBin()));
	poiss1->SetParameter(1,QPerPMT[i]->GetMean());
	TFitResultPtr r = QPerPMT[i]->Fit(poiss1,"SRQ");
	Double_t par1 = r->Parameter(1);
	Double_t par0 = r->Parameter(0);
	Double_t err_par1 = r->ParError(1);      
	mean_charge_vec->push_back(par1);
	mean_charge_err_vec->push_back(err_par1);

      }
      else if (distance > 450) { 
	poiss2->SetParameter(0,QPerPMT[i]->GetBinContent(QPerPMT[i]->GetMaximumBin()));
	poiss2->SetParameter(1,QPerPMT[i]->GetMean());
	TFitResultPtr r = QPerPMT[i]->Fit(poiss2,"SRQ");
	Double_t par1 = r->Parameter(1);
	Double_t par0 = r->Parameter(0);
	Double_t err_par1 = r->ParError(1);
	mean_charge_vec->push_back(par1);
	mean_charge_err_vec->push_back(err_par1);
      }
    }
  }
}

void Process_Input_File(TFile* inFile,std::vector<double>* distance_vec, std::vector<double>* mean_charge_vec, std::vector<double>* eff_vec, std::vector<double>* mean_charge_err_vec, std::vector<double>* eff_err_vec) {

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

  int nbPEMaxByPMT = max_charge*0.75;

  /*TH1F** QPerPMT = new TH1F*[(int)num_OD_PMT];
  for (int i=0; i<num_OD_PMT; i++) {
    QPerPMT[i] = new TH1F(Form("tq%d",i),"test1",5*(nbPEMaxByPMT+20),-10,nbPEMaxByPMT+10);
    }*/
  std::vector<TH1F*> QPerPMT;
  for (int i=0; i<num_OD_PMT; i++) {
    TH1F* Qhist = new TH1F(Form("tq%d",i),"test1",5*(nbPEMaxByPMT+20),-10,nbPEMaxByPMT+10);
    QPerPMT.push_back(Qhist);
  }

  int n = inTree->Draw("eventId:tubeId:trueTime:recTime:trueNumPhotons:recCharge:cylLoc:pmtX:pmtY:pmtZ:nearestInjDist","","goff");

  std::vector<double> tube_inj_dists(num_OD_PMT);
  std::vector<double> xpos(num_OD_PMT);
  std::vector<double> ypos(num_OD_PMT);
  std::vector<double> zpos(num_OD_PMT);
  std::vector<int> cyl(num_OD_PMT);

  double barrelposx = -99999;
  double barrelposy = -99999;
  double topcappos = -99999;

  for (Long64_t i=0; i < n; i++) {
    //
    inTree->GetEntry(i);
    if (cylLoc == 4) {
      barrelposx = pmtX;
      barrelposy = pmtY;
    }

    cyl.at(tubeId-1) = cylLoc;
    if (cylLoc == 5) {
      topcappos = pmtZ;
      if (recCharge > 0) QPerPMT.at(tubeId-1)->Fill(recCharge);
      xpos.at(tubeId-1) = pmtX;
      ypos.at(tubeId-1) = pmtY;
      zpos.at(tubeId-1) = pmtZ;
      double pmtXYZ[3] = {pmtX, pmtY, pmtZ};
      tube_inj_dists.at(tubeId-1) = PMT_Injector_Distance(pmtXYZ);
    } 
  }

  double RadiusOD = sqrt( pow(barrelposx,2) + pow(barrelposy,2) );
  double HeightOD = 2*(abs(topcappos));

  Fill_Q_Eff_Vectors(distance_vec, mean_charge_vec, eff_vec, mean_charge_err_vec, eff_err_vec, num_OD_PMT, cyl, QPerPMT, nbPEMaxByPMT, tube_inj_dists);
}

// A function to convert radians to degress
float RadToDeg(float x){
  return x*180/PI;
}
// A function to convert degress to radians
float DegToRad(float x){
  return x*PI/180;
}

void OD_Top_Strut_Comparison_Plots( const char *inFileName1 = "wcsim.root", const char *inFileName2 = "wcsim.root", const char *in_or_out = "out", const char *bare_or_diffuse = "diffuse", int MaxPE = 20, bool verbosity = 0){
  /*Produce mean charge and illumination efficiency ratio plots (as a function of distance),
    comparing two input files.
    Ratios will be file 1 / file 2
  */

  // Open the WCSim files
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


  //define the vectors we're going to use
  std::vector<double> distance1, mean_charge1, efficiency1, charge_error1, efficiency_error1;
  std::vector<double> reduced_distance1, reduced_mean_charge1, reduced_efficiency1, reduced_charge_error1, reduced_efficiency_error1;

  std::vector<double> distance2, mean_charge2, efficiency2, charge_error2, efficiency_error2;
  std::vector<double> reduced_distance2, reduced_mean_charge2, reduced_efficiency2, reduced_charge_error2, reduced_efficiency_error2;

  //Read in and process the mean charge and efficiency per PMT
  Process_Input_File(inFile1, &distance1, &mean_charge1, &efficiency1, &charge_error1, &efficiency_error1);
  Process_Input_File(inFile2, &distance2, &mean_charge2, &efficiency2, &charge_error2, &efficiency_error2);

  //Sort by distance from nearest injector and reduce data down by factor (10)
  Sort_And_Reduce_Vectors(distance1, mean_charge1, efficiency1, charge_error1, efficiency_error1, &reduced_distance1, &reduced_mean_charge1, &reduced_efficiency1, &reduced_charge_error1, &reduced_efficiency_error1);
  Sort_And_Reduce_Vectors(distance2, mean_charge2, efficiency2, charge_error2, efficiency_error2, &reduced_distance2, &reduced_mean_charge2, &reduced_efficiency2, &reduced_charge_error2, &reduced_efficiency_error2);
  //Now make comparison plots
   std::vector<Double_t> charge_ratio;
   std::vector<Double_t> charge_ratio_error;
   std::vector<Double_t> efficiency_ratio;
   std::vector<Double_t> efficiency_ratio_error;
   for (int i = 0; i < reduced_distance1.size(); i++) {
     Double_t charge_1 = reduced_mean_charge1.at(i);
     Double_t charge_2 = reduced_mean_charge2.at(i);
     Double_t charge_err_1 = reduced_charge_error1.at(i);
     Double_t charge_err_2 = reduced_charge_error2.at(i);
     Double_t charge12_ratio = charge_1/charge_2;
     Double_t charge12_comb_error = charge12_ratio*(sqrt((charge_err_1/charge_1)*(charge_err_1/charge_1) + (charge_err_2/charge_2)*(charge_err_2/charge_2)));
     charge_ratio.push_back(charge12_ratio);
     charge_ratio_error.push_back(charge12_comb_error);

     Double_t efficiency_1 = reduced_efficiency1.at(i);
     Double_t efficiency_2 = reduced_efficiency2.at(i);
     Double_t efficiency_err_1 = reduced_efficiency_error1.at(i);
     Double_t efficiency_err_2 = reduced_efficiency_error2.at(i);
     Double_t efficiency12_ratio = efficiency_1/efficiency_2;
     Double_t efficiency12_comb_error = efficiency12_ratio*(sqrt((efficiency_err_1/efficiency_1)*(efficiency_err_1/efficiency_1) + (efficiency_err_2/efficiency_2)*(efficiency_err_2/efficiency_2)));
     efficiency_ratio.push_back(efficiency12_ratio);
     efficiency_ratio_error.push_back(efficiency12_comb_error);
   }

   /*   std::cerr << "Q size " << mean_charge1.size() << " dist size " << distance1.size() << std::endl;
   TH1F* h=new TH1F("h", "You're a bastard, ROOT", 1000, 0, 10000);
   TH2* h2 = new TH2D("h2", "You're a rat bastard, ROOT", 1000,0,10000,10,0,5);

   for (int i = 0; i<distance1.size(); i++) {
     h->Fill(distance1.at(i));
     h2->Fill(distance1.at(i), mean_charge1.at(i));
     if (distance1.at(i) > 1000) std::cerr << distance1.at(i) << std::endl;
   }
   TCanvas *ch = new TCanvas("ch","die in a fire",800,600);
   h->Draw();

   TCanvas *ch2 = new TCanvas("ch2","go to hell",800,600);
   h2->Draw();

   
   
   TGraphErrors *mean_charge_dist_plot1 = new TGraphErrors(distance1.size(), &(reduced_distance1[0]), &(mean_charge1[0]), 0, &(charge_error1[0]) );
   TCanvas *c4 = new TCanvas("c4","Mean charge (file 1)  vs dist",800,600);
   mean_charge_dist_plot1->Draw("AP");
   std::cerr << mean_charge_dist_plot1->GetN() << std::endl;
   TGraphErrors *efficiency_dist_plot1 = new TGraphErrors(distance1.size(), &(reduced_distance1[0]), &(efficiency1[0]), 0, &(efficiency_error1[0]) );
   TCanvas *c5 = new TCanvas("c5","Efficiency (file 2) vs dist",800,600);
   efficiency_dist_plot1->Draw("AP");
   std::cerr << efficiency_dist_plot1->GetN() << std::endl;
   TGraphErrors *mean_charge_dist_plot2 = new TGraphErrors(distance2.size(), &(reduced_distance2[0]), &(mean_charge2[0]), 0, &(charge_error2[0]) );
   TCanvas *c42 = new TCanvas("c42","Mean charge (file 1)  vs dist",800,600);
   mean_charge_dist_plot2->Draw("AP");
   std::cerr << mean_charge_dist_plot2->GetN() << std::endl;
   TGraphErrors *efficiency_dist_plot2 = new TGraphErrors(distance2.size(), &(reduced_distance2[0]), &(efficiency2[0]), 0, &(efficiency_error2[0]) );
   TCanvas *c52 = new TCanvas("c52","Efficiency (file 2) vs dist",800,600);
   efficiency_dist_plot2->Draw("AP");
   std::cerr << efficiency_dist_plot2->GetN() << std::endl;
   */
   TGraphErrors *charge_ratio_plot = new TGraphErrors(reduced_distance1.size(), &(reduced_distance1[0]), &(charge_ratio[0]), 0, &(charge_ratio_error[0]) );
   TCanvas *c3 = new TCanvas("c3","Mean charge ratio",800,600);
   charge_ratio_plot->SetMarkerColor(4);
   charge_ratio_plot->SetMarkerStyle(8);
   charge_ratio_plot->GetXaxis()->SetRangeUser(0,1000);
   charge_ratio_plot->SetTitle("");
   charge_ratio_plot->GetYaxis()->SetTitle("Charge ratio");
   charge_ratio_plot->GetXaxis()->SetTitle("Distance between PMT and injector (cm)");
   charge_ratio_plot->Draw("AP");
   TLine *line1 = new TLine(0,1,1000,1);
   line1->SetLineColor(2);
   line1->SetLineStyle(9);
   line1->Draw();

   TGraphErrors *efficiency_ratio_plot = new TGraphErrors(reduced_distance1.size(), &(reduced_distance1[0]), &(efficiency_ratio[0]), 0, &(efficiency_ratio_error[0]) );
   TCanvas *c2 = new TCanvas("c2","Efficiency ratio",800,600);
   efficiency_ratio_plot->SetMarkerColor(4);
   efficiency_ratio_plot->SetMarkerStyle(8);
   efficiency_ratio_plot->GetXaxis()->SetRangeUser(0,1000);
   efficiency_ratio_plot->GetYaxis()->SetTitle("Efficiency ratio");
   efficiency_ratio_plot->GetXaxis()->SetTitle("Distance between PMT and injector (cm)");
   efficiency_ratio_plot->SetTitle("");
   efficiency_ratio_plot->Draw("AP");
   TLine *line2 = new TLine(0,1,1000,1);
   line2->SetLineColor(2);
   line2->SetLineStyle(9);
   line2->Draw();

}
