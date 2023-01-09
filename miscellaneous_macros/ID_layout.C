#include <stdio.h>     
#include <stdlib.h>    
// C++ Includes
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <utility>


// ROOT Includes
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TEllipse.h"
#include "TBox.h"

void ID_layout() {

  TH2F* layout = new TH2F("layout","layout", 2400,-12000,12000,2200,-11000,11000);
  double RadiusOD = 3290.0;
  double HeightOD = 6780.0;
  double rad1 = 400;
  double xposb = -999;
  double yposb = -999;
  double zposb = -999;
    // barrel
  for(int ip = 0; ip < 7; ip++) {
    if (ip==0) zposb = 1800;
    if (ip==1) zposb = 900;
    if (ip==2) zposb = 0;
    if (ip==3) zposb = -900;
    if (ip==4) zposb = -1800;
    if (ip==5) zposb = 2700;
    if (ip==6) zposb = -2700;
    for(int theta = 0; theta < 4; theta++) {
      std::vector<double> single_barrel_pos (3);
      xposb = RadiusOD*(std::cos((theta+0.5)*(std::acos(-1))/2));
      yposb = RadiusOD*(std::sin((theta+0.5)*(std::acos(-1))/2));
      double pmtRad = sqrt((pow(xposb, 2.0))+(pow(yposb, 2.0)));
      double pmtTheta = std::atan(yposb/xposb);
      double l = sqrt((pow((0-xposb),2)) + pow((-(HeightOD/2)-yposb),2));
      double angle = 2*asin(l/2*RadiusOD);
      double length = angle*RadiusOD;
      if (xposb>0) pmtTheta += (std::acos(-1)/2);
      else pmtTheta -= (std::acos(-1)/2);
      double pmtRadTheta = pmtRad*pmtTheta;
      layout->Fill(pmtRadTheta,zposb);
      layout->Fill(pmtRadTheta,zposb);
    }
  }

  // bottom
  double xpos = 0.;
  double ypos = 0.;
  //G4ThreeVector postemp(0, 0, zpos);
  for (int angfrac1 = 0; angfrac1 < 4; angfrac1++) {
    xpos = rad1*(std::cos(angfrac1*(std::acos(-1))/2));
    ypos = rad1*(std::sin(angfrac1*(std::acos(-1))/2));

    double ODpolyX = xpos;
    double ODpolyY = ypos + RadiusOD + HeightOD/2 + 100;
    layout->Fill(ODpolyX,-ODpolyY);
    layout->Fill(ODpolyX,-ODpolyY);

  }

  // top
  for (int angfrac2 = 0; angfrac2 < 4; angfrac2++) {
    xpos = rad1*(std::cos(angfrac2*(std::acos(-1))/2));
    ypos = rad1*(std::sin(angfrac2*(std::acos(-1))/2));

    double ODpolyX = xpos;
    double ODpolyY = ypos + RadiusOD + HeightOD/2 + 100;
    layout->Fill(ODpolyX,ODpolyY);
  }

gStyle->SetOptStat(0);

double pi = 2*acos(0.0);
double length = pi*2*RadiusOD;

  TCanvas *lay = new TCanvas("lay","lay",800,700);
  layout->SetXTitle("cm");
  layout->SetYTitle("cm");
  layout->SetTitle("");
  layout->SetMarkerStyle(43);
  layout->SetMarkerSize(1);
  layout->Draw("");
TEllipse *el1 = new TEllipse(0,(RadiusOD + HeightOD/2 + 100),RadiusOD,RadiusOD);
TEllipse *el2 = new TEllipse(0,-(RadiusOD + HeightOD/2 + 100),RadiusOD,RadiusOD);
TBox* b1 = new TBox(-(length/2),-(HeightOD/2),(length/2),(HeightOD/2));
//el1->SetFillColorAlpha(kBlue,1.0);
el1->SetFillStyle(0);
el1->Draw();
el2->SetFillStyle(0);
el2->Draw();
b1->SetFillStyle(0);
b1->Draw();




} // END MACRO
