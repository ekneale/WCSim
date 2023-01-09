#include <iostream>
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <utility>


// ROOT Includes
#include "TGraph.h"
#include "TMultiGraph.h"


void Barrel_Water_Params_chi2_plot(){
gROOT->Reset();
//  gStyle->SetOptFit();
TCanvas *c1 = new TCanvas("c1","multigraph",200,10,700,500);
//c1->SetGrid();
// double barrel_err = sqrt(2*1975);
// draw a frame to define the range
TMultiGraph *mg = new TMultiGraph();

// create first graph ABWFF
Int_t n1 = 20;
//Double_t x1[]  = {-30,-25,-20,-15,-10,-5,-4,-3,-2,-1,0,1,2,3,4,5,10,15,20,25,30};
//Double_t ex1[] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
 Double_t x1[]  = {-30,-25,-20,-15,-10,-5,-4,-3,-2,-1,1,2,3,4,5,10,15,20,25,30};
 Double_t ex1[] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
 Double_t y1[] = {26.3313/2.5, 17.6586/2.5, 11.9361/2.5, 8.64058/2.5, 6.20661/2.5, 4.37403/2.5, 4.88391/2.5, 4.41495/2.5, 4.32038/2.5, 4.4696/2.5, 4.77374/2.5, 4.52944/2.5, 4.55365/2.5, 4.55956/2.5, 4.6255/2.5, 5.42336/2.5, 7.26046/2.5, 8.29323/2.5, 9.95565/2.5, 11.9897/2.5};
 Double_t ey1[] = {1.74051/2.5, 1.12192/2.5, 0.616571/2.5, 0.556766/2.5, 0.37754/2.5, 0.255563/2.5, 0.265368/2.5, 0.285464/2.5, 0.265428/2.5, 0.195626/2.5, 0.278404/2.5, 0.312292/2.5, 0.334248/2.5, 0.265734/2.5, 0.253861/2.5, 0.22094/2.5, 0.270918/2.5, 0.560449/2.5, 0.849519/2.5, 0.457802/2.5};

TGraphErrors *gr1 = new TGraphErrors(n1,x1,y1,ex1,ey1);
gr1->SetMarkerColor(kBlue);
 gr1->SetMarkerStyle(4);
 gr1->SetMarkerSize(0.9);
 gr1->SetName("Absorption parameter");
//gr1->Fit("pol6","q");
mg->Add(gr1);

// create second graph RAYFF
Int_t n2 = 20;
Double_t x2[]  = {-30,-25,-20,-15,-10,-5,-4,-3,-2,-1,1,2,3,4,5,10,15,20,25,30};
Double_t ex2[] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
 Double_t y2[]  = {7630.3/2.5, 895.349/2.5, 562.494/2.5, 326.582/2.5, 167.306/2.5, 70.0387/2.5, 20.0063/2.5, 13.9664/2.5, 9.28972/2.5, 6.83962/2.5, 5.11507/2.5, 4.84045/2.5, 6.37489/2.5, 9.89094/2.5, 13.4513/2.5, 16.0542/2.5, 52.1575/2.5, 105.129/2.5, 167.385/2.5, 240.999/2.5, 322.435/2.5};
 Double_t ey2[]  = {256.319/2.5, 10.2631/2.5, 14.246/2.5, 7.29874/2.5, 5.97562/2.5, 4.33149/2.5, 1.57242/2.5, 0.834105/2.5, 0.554619/2.5, 0.450913/2.5, 0.284505/2.5, 0.238873/2.5, 0.39012/2.5, 0.497306/2.5, 0.64318/2.5, 0.930081/2.5, 2.54884/2.5, 3.28662/2.5, 6.58651/2.5, 6.57027/2.5, 6.75077/2.5};

TGraphErrors *gr2 = new TGraphErrors(n2,x2,y2,ex2,ey2);
gr2->SetMarkerColor(kRed);
gr2->SetMarkerStyle(20);
 gr2->SetMarkerSize(0.9);
 gr2->SetName("Rayleigh parameter");
//gr2->Fit("pol5","q");
   
mg->Add(gr2);
   
mg->GetHistogram()->GetXaxis()->SetTitle("Percentage change from default value");
mg->GetHistogram()->GetYaxis()->SetTitle("#chi^{2}/ndof");

mg->Draw("ap");


 TLegend *mglegend = new TLegend(0.6,0.75,0.9,0.9);
 mglegend->AddEntry(gr1, "Absorption parameter","p");
 mglegend->AddEntry(gr2, "Rayleigh parameter","p");
 mglegend->Draw();
   
//force drawing of canvas to generate the fit TPaveStats
c1->Update();
  mg->GetYaxis()->SetRangeUser(1,300);
 gPad->SetLogy(1);
 gPad->SetGridy(1);


 c1->SaveAs("SDO_abs_ray_chi2_barrel.png");
}
