#include <iostream>
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <utility>


// ROOT Includes
#include "TGraph.h"
#include "TMultiGraph.h"


void refl_chi2_plot(){
gROOT->Reset();
//  gStyle->SetOptFit();
TCanvas *c1 = new TCanvas("c1","multigraph",200,10,700,500);
//c1->SetGrid();
// double bottom_err = sqrt(2*1975);
// draw a frame to define the range
TMultiGraph *mg = new TMultiGraph();

// create first graph REFL
Int_t n1 = 18;
 Double_t x1[]  = {-30,-20,-15,-10,-9,-8,-7,-6,-5,-4.5,-4,-3.5,-3,-2.5,-2,-1.5,-1,-0.5,};
 Double_t ex1[] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  //  Double_t y1[] = {2717.68, 1921.14, 1430.81, 884.079, 765.562, 655.318, 541.308, 435.18, 328.247, 278.434, 231.3, 187.235, 143.213, 105.88, 72.8843, 44.7674, 22.972, 8.74484};
 Double_t y1[] = {7723.75/2.5, 5535.84/2.5, 4153.88/2.5, 2589.37/2.5, 2253.26/2.5, 1927.32/2.5, 1598.44/2.5, 1271.81/2.5, 963.771/2.5, 820.681/2.5, 675.201/2.5, 543.658/2.5, 421.297/2.5, 308.574/2.5, 208.565/2.5, 125.275/2.5, 62.7479/2.5, 18.6919/2.5};
  //Double_t ey1[18] = {7.28444, 10.6847, 8.12795, 4.56112, 5.61642, 6.18975, 4.92185, 4.64027, 3.25093, 3.3053, 2.82838, 2.66562, 3.26694, 2.30859, 1.23134, 0.716781, 0.760094, 0.295704};
 Double_t ey1[] = {31.2844/2.5, 30.4626/2.5, 33.4265/2.5, 32.3336/2.5, 25.9555/2.5, 20.2298/2.5, 15.3798/2.5, 17.792/2.5, 15.9858/2.5, 13.9338/2.5, 12.3094/2.5, 11.82/2.5, 12.1411/2.5, 6.79811/2.5, 9.16343/2.5, 4.67516/2.5, 2.15503/2.5, 1.09664/2.5};

TGraphErrors *gr1 = new TGraphErrors(n1,x1,y1,ex1,ey1);
gr1->SetMarkerColor(kGreen+2);
gr1->SetMarkerStyle(52);
 gr1->SetMarkerSize(1.3);
 gr1->SetName("Reflectivity parameter");
//gr1->Fit("pol6","q");
mg->Add(gr1);
   
mg->GetHistogram()->GetXaxis()->SetTitle("Percentage change in reflectivity");
mg->GetHistogram()->GetYaxis()->SetTitle("#chi^{2}/ndof");

mg->Draw("ap");

/*
 TLegend *mglegend = new TLegend(0.6,0.75,0.9,0.9);
 mglegend->AddEntry(gr1, "Reflectivity parameter","p");
 mglegend->Draw();
   */
//force drawing of canvas to generate the fit TPaveStats
c1->Update();
 mg->GetYaxis()->SetRangeUser(6,3500);
 gPad->SetLogy(1);
 gPad->SetGridy(1);


  c1->SaveAs("SDO_refl_chi2_barrel.png");
}
