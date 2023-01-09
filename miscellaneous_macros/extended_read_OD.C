#include <stdio.h>     
#include <stdlib.h>    
// C++ Includes
#include <iostream>


// ROOT Includes
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"

// WCSim Includes
#include "include/WCSimRootEvent.hh"
#include "include/WCSimRootGeom.hh"

// #include <libgen.h>

void extended_read_OD(char *filename=NULL) {
  /* A simple script to plot aspects of phototube hits 
   * 
   * I like to run this macro as 
   * $ root -l -x 'extended_read_OD.C("OD.root")'
   */

  gROOT->Reset();
  char* wcsimdirenv;
  wcsimdirenv = getenv ("WCSIMDIR");
  if(wcsimdirenv !=  NULL){
    gSystem->Load("${WCSIMDIR}/libWCSimRoot.so");
    // gSystem->Load("${WCSIMDIR}/libWCSimRoot.rootmap");
    gSystem->Load("${WCSIMDIR}/WCSimRootDict_rdict.pcm");
  }else{
    std::cout << "Can't load WCSim ROOT dictionaries" << std::endl;
  }
  gStyle->SetOptStat(1);

  TFile *f;
  char fTest[128];
  if (filename==NULL){
    std::cout << "Please provide filename in option" << std::endl;
    std::cout << "Will load auto wcsim.root in WCSIMDIR ..." << std::endl;
    char* name = "wcsim.root";
    strncpy(fTest, wcsimdirenv, sizeof(fTest));
    strncat(fTest, "/", (sizeof(fTest) - strlen(fTest)) );
    strncat(fTest, name, (sizeof(fTest) - strlen(fTest)) );
    f = new TFile(fTest);
  }else{
    f = new TFile(filename);
  }
  if (!f->IsOpen()){
    cout << "Error, could not open input file: " << filename << endl;
    return -1;
  }else{
    if (filename==NULL) cout << "File open bro: " << fTest << endl;
    else cout << "File open bro: " << filename << endl;
  }

  char *dirc, *basec, *bname, *dname;

  double WLSLength = 80;
  //double WLSLength = 120;
  //  double WLSLength = 40;

  // dirc = strdup(filename);
  // basec = strdup(filename);
  // dname = dirname(dirc);
  // bname = basename(basec);
  // printf("dirname=%s, basename=%s\n", dname, bname);

  // Get a pointer to the tree from the input file
  TTree *wcsimT    = (TTree*) f->Get("wcsimT");
	
  // Get the number of events in the tree
  long int nEvent = wcsimT->GetEntries();
  std::cout << "Number of events: "<< nEvent << std::endl;


  //  TTree  *wcsimT = (TTree*)f->Get("wcsimT");

  WCSimRootEvent *wcsimrootsuperevent = new WCSimRootEvent();
  //wcsimT->SetBranchAddress("wcsimrootevent_OD",&wcsimrootsuperevent);
  TBranch *branch = wcsimT->GetBranch("wcsimrootevent_OD");
  branch->SetAddress(&wcsimrootsuperevent);


  TTree *wcsimGeoT = (TTree*)f->Get("wcsimGeoT");
  //WCSimRootGeom* wcsimrootgeom = new WCSimRootGeom();
  WCSimRootGeom *wcsimrootgeom = 0;
  wcsimGeoT->SetBranchAddress("wcsimrootgeom",&wcsimrootgeom);
  std::cerr << "Geotree has " << wcsimGeoT->GetEntries() << " entries." << std::endl;
  //wcsimrootgeom->GetEntry(0);
  wcsimGeoT->GetEntry(0);

  // Force deletion to prevent memory leak when issuing multiple
  // calls to GetEvent()
  //wcsimT->GetBranch("wcsimrootevent_OD")->SetAutoDelete(kTRUE);

  // const long unsigned int nbEntries = wcsimT->GetEntries();
  const long unsigned int nbEntries = wcsimT->GetEntries();
  cout << "Nb of entries " << wcsimT->GetEntries() << endl;

  //////////////////////////////////////////
  // HISTOGRAMS DEFINITION /////////////////
  //////////////////////////////////////////

  const int nbBins = 50;
  const int nbPEMax = 5000;
  const int nbBinsByPMT = 25;
  const int nbPEMaxByPMT = 200;

  gStyle->SetOptStat(0);
  
  TH1D *hPEByEvtsByPMT = new TH1D("hPEByEvtsByPMT","RAW PE by Evts by PMT",
				  nbBinsByPMT,0,nbPEMaxByPMT);
  hPEByEvtsByPMT->GetXaxis()->SetTitle("raw PE");
  hPEByEvtsByPMT->SetLineColor(kBlue-4);
  hPEByEvtsByPMT->SetMarkerColor(kBlue-4);  
  TH1D *hPECollectedByEvtsByPMT = new TH1D("hPECollectedByEvtsByPMT","collected PE by Evts by PMT",
					   nbBinsByPMT,0,nbPEMaxByPMT);
  hPECollectedByEvtsByPMT->GetXaxis()->SetTitle("digi PE");
  hPECollectedByEvtsByPMT->SetLineColor(kRed-4);
  hPECollectedByEvtsByPMT->SetMarkerColor(kRed-4);

  //TH1F *digiHitsPerPMT = new TH1F("digiHitsPerPMT", "digiHitsPerPMT", nBinsHits, hitsMin, hitsMax );
  
  TH1D *hPEByEvts = new TH1D("hPEByEvts","Total RAW PE by Evts",nbBins,0,nbPEMax);
  hPEByEvts->GetXaxis()->SetTitle("raw PE");
  hPEByEvts->SetLineColor(kBlue+1);
  hPEByEvts->SetMarkerColor(kBlue+1);
  hPEByEvts->SetFillColor(kBlue+1);  
  TH1D *hPECollectedByEvts = new TH1D("hPECollectedByEvts","Total collected PE by Evts",nbBins,0,nbPEMax);
  hPECollectedByEvts->GetXaxis()->SetTitle("digi PE");
  hPECollectedByEvts->SetLineColor(kRed+1);
  hPECollectedByEvts->SetMarkerColor(kRed+1);
  hPECollectedByEvts->SetFillColor(kRed+1);

  TH1D *hNbTubesHit = new TH1D("hNbTubesHit","Nb of Tubes Hit",50,0,1000);

  //TH3D *BarrelHits = new TH3D("BarrelHits","Hits in barrel",740,-37,37,740,-37,37,800,-40,40);

  //TH2D *BarrelHits = new TH2D("BarrelHits","Hits in barrel",468,-11700,11700,120,-3000,3000);
  //TH2D *BarrelHits = new TH2D("BarrelHits","Hits in barrel",116,-5800,5800,60,-3000,3000);
  //TH2D *BarrelHits = new TH2D("BarrelHits","Hits in barrel",133,-11310,11310,33,-2800,2800);
  // Ideal for 8 inch PMT
  //TH2D *BarrelHits = new TH2D("BarrelHits","Hits in barrel",50,-12000,12000,28,-2800,2800);
  /*  TH2D *BarrelHits = new TH2D("BarrelHits","Hits in barrel",72,-12000,12000,28,-2800,2800);
  TH2D *TopCapHits = new TH2D("TopCapHits","Hits in top cap",22,-4000,4000,22,-4000,4000);
  TH2D *BottomCapHits = new TH2D("BottomCapHits","Hits in bottom cap",22,-4000,4000,22,-4000,4000);*/
  // Use this for 3 inch PMTs
  //TH2D *BarrelHits = new TH2D("BarrelHits","Hits in barrel",200,-12000,12000,35,-2800,2800);
  TH2D *BarrelHits = new TH2D("BarrelHits","Hits in barrel",129,-12000,12000,28,-2800,2800);
  TH2D *TopCapHits = new TH2D("TopCapHits","Hits in top cap",26,-4000,4000,26,-4000,4000);
  TH2D *BottomCapHits = new TH2D("BottomCapHits","Hits in bottom cap",26,-4000,4000,26,-4000,4000);
  //TH2D *BarrelHits = new TH2D("BarrelHits","Hits in barrel",1330,-11310,11310,330,-2800,2800);
  //TH2D *TopCapHits = new TH2D("TopCapHits","Hits in top cap",100,-3700,3700,100,-3700,3700);
  
  //TH2D *BottomCapHits = new TH2D("BottomCapHits","Hits in bottom cap",23,-4000,4000,23,-4000,4000);
  //TH2D *BottomCapHits = new TH2D("BottomCapHits","Hits in bottom cap",37,-3500,3500,35,-3500,3500);

  /*TH2Poly *BarrelPoly = new TH2Poly("BarrelHits","Hits in barrel",-12000,12000,-2800,2800);
  TH2Poly *TopCapPoly = new TH2Poly("TopCapHits","Hits in top cap",-4000,4000,-4000,4000);
  TH2Poly *BottomCapPoly = new TH2Poly("BottomCapHits","Hits in bottom cap",-4000,4000,-4000,4000);*/

  TH2Poly *ODpoly = new TH2Poly("ODHits","Charge in OD PMTs",-12000,12000,-11000,11000);
  TH2Poly *ODPEpoly = new TH2Poly("ODPE","Photons incident on OD PMTs",-12000,12000,-11000,11000);

  // Find a barrel pmt
  bool barrelOD = false;  // Boolean to see if a barrel OD pmt has been found
  bool capOD = false;  // Boolean to see if a cap OD pmt has been found
  int barrelPMTOD = -1; // PMT number of the barrel OD PMT
  int capPMTOD = -1; // PMT number of the cap OD PMT
  int pmtCount = 0; // Number used to count through all of the PMTs

  while (!barrelOD || !capOD ){ // Loop to look for barrel and cap PMTs to work out the radius and height respectively.

    if ( !barrelOD && (wcsimrootgeom->GetPMT(pmtCount).GetCylLoc() == 4 )  ) {barrelOD = true; barrelPMTOD = pmtCount; }
    if ( !capOD && (wcsimrootgeom->GetPMT(pmtCount).GetCylLoc() == 3 || wcsimrootgeom->GetPMT(pmtCount).GetCylLoc() == 5)  ) {capOD = true; capPMTOD = pmtCount; }
    pmtCount++;
    //pmtCount+= 10; // Can speed up this process by checking PMTs in multiples higher than 1
  }


  double num_OD_PMT = wcsimrootgeom->GetODWCNumPMT();
  std::cerr << num_OD_PMT << std::endl;
  double num_OD_PMTs_with_hits = 0;
  std::vector<double> coverage_i;
  double RadiusOD = sqrt( pow(wcsimrootgeom->GetPMT(barrelPMTOD).GetPosition(0),2) + pow(wcsimrootgeom->GetPMT(barrelPMTOD).GetPosition(1),2) );
  double HeightOD = 2*(abs(wcsimrootgeom->GetPMT(capPMTOD).GetPosition(2)));

  std::cerr << "RadiusOD " << RadiusOD << ", HeightOD " << HeightOD << std::endl;

  int MAXPMT = wcsimrootgeom->GetWCNumPMT(); // Get the maximum number of PMTs in the ID
  int MAXPMTA = wcsimrootgeom->GetODWCNumPMT(); // Get the maximum number of PMTs in the OD

  int num_barrel_PMTs = 0;
  
  for (int polybin = (MAXPMT - 1); polybin < (MAXPMT + num_OD_PMT); polybin++) {
    //
    double tube[3];
    int cylLoc = wcsimrootgeom->GetPMT(polybin).GetCylLoc();
    tube[0] = wcsimrootgeom->GetPMT(polybin).GetPosition(0);
    tube[1] = wcsimrootgeom->GetPMT(polybin).GetPosition(1);
    tube[2] = wcsimrootgeom->GetPMT(polybin).GetPosition(2);

    if ( cylLoc == 5){
      ODpoly->AddBin(tube[0] - WLSLength/2, tube[1] + RadiusOD + HeightOD/2 + 100 - WLSLength/2, tube[0] + WLSLength/2, tube[1] + RadiusOD + HeightOD/2 + 100 + WLSLength/2);
      ODPEpoly->AddBin(tube[0] - WLSLength/2, tube[1] + RadiusOD + HeightOD/2 + 100 - WLSLength/2, tube[0] + WLSLength/2, tube[1] + RadiusOD + HeightOD/2 + 100 + WLSLength/2);
    }
    //Bot OD
    else if ( cylLoc == 3){
      //std::cerr << "xpos " << tube[0] << ", ypos " << tube[1] << ", zpos " << tube[2] << std::endl;
      ODpoly->AddBin(tube[0] - WLSLength/2,-(HeightOD/2 +RadiusOD +tube[1] + 100) - WLSLength/2, tube[0] + WLSLength/2,-(HeightOD/2 +RadiusOD +tube[1] + 100) + WLSLength/2);
      ODPEpoly->AddBin(tube[0] - WLSLength/2,-(HeightOD/2 +RadiusOD +tube[1] + 100) - WLSLength/2, tube[0] + WLSLength/2,-(HeightOD/2 +RadiusOD +tube[1] + 100) + WLSLength/2);
    }
    //Barrel OD
    else {

      num_barrel_PMTs++;
      double pmtRad = sqrt((pow(tube[0], 2.0))+(pow(tube[1], 2.0)));
      double pmtTheta = std::atan(tube[1]/tube[0]);
      if (tube[0]>0) pmtTheta += (std::acos(-1)/2);
      else pmtTheta -= (std::acos(-1)/2);
      double pmtRadTheta = pmtRad*pmtTheta;
      
      //      double l = sqrt( pow((0 - tube[0]),2) + pow((-RadiusOD - tube[1]),2));
      //     double angle = 2*asin(l/(2*RadiusOD));
      //    double length = angle*RadiusOD ;
      //if (tube[0]<0) length *= -1;
      //ODpoly->AddBin(length - WLSLength/2, tube[2] - WLSLength/2, length + WLSLength/2, tube[2] + WLSLength/2);
      //ODPEpoly->AddBin(length - WLSLength/2, tube[2] - WLSLength/2, length + WLSLength/2, tube[2] + WLSLength/2);
      ODpoly->AddBin(pmtRadTheta - WLSLength/2, tube[2] - WLSLength/2, pmtRadTheta + WLSLength/2, tube[2] + WLSLength/2);
      ODPEpoly->AddBin(pmtRadTheta - WLSLength/2, tube[2] - WLSLength/2, pmtRadTheta + WLSLength/2, tube[2] + WLSLength/2);
    }

  }

  std::cerr << "num_barrel_PMTs " << num_barrel_PMTs << std::endl;
  
  TH1D *PE = new TH1D("PEmult","Photoelectron multiplicty", 16,-0.5,15.5);
  PE->SetXTitle("Photoelectrons");

  TH1D *PMT_hits = new TH1D("PMT_hits","Hits vs PMT detector number", 20000,-0.5,19999.5);
  PMT_hits->SetXTitle("PMT ID");
  PMT_hits->SetYTitle("Photoelectrons");
  PMT_hits->SetLineWidth(2);

  TH1D *check_PMT_hits = new TH1D("check_PMT_hits","check_Hits vs PMT detector number", 20000,-0.5,19999.5);
  check_PMT_hits->SetXTitle("PMT ID");
  check_PMT_hits->SetYTitle("Photoelectrons");
  check_PMT_hits->SetLineWidth(2);
  
  //PMT_hits->SetLabelSize(10);
  TH1D *PMT_hits100 = new TH1D("PMT_hits","Hits (>100) vs PMT detector number", 12000,-0.5,11999.5);

  TH2D *QvsTubeID = new TH2D("QvsTubesID","charge vs. PMTs hit", 20000,0,20000,100,0,100);
  QvsTubeID->SetXTitle("PMT ID");
  QvsTubeID->SetYTitle("Q");
  
  
  
  // END HISTOGRAMS DEFINITION /////////////
  //////////////////////////////////////////

  double total_charge = 0;
  double check_total_charge = 0;
  double trigger_Q = 0;
  for(long unsigned int iEntry = 0; iEntry < nbEntries; iEntry++){
    // Point to event iEntry inside WCSimTree
    wcsimT->GetEvent(iEntry);
    double num_hit_OD_PMTs = 0;
    /*    TH1D *check_PMT_hits = new TH1D("check_PMT_hits","check_Hits vs PMT detector number", 20000,-0.5,19999.5);
    check_PMT_hits->SetXTitle("PMT ID");
    check_PMT_hits->SetYTitle("Photoelectrons");
    check_PMT_hits->SetLineWidth(2);
    */
    // Nb of Trigger inside the event
    const unsigned int nbTriggers = wcsimrootsuperevent->GetNumberOfEvents();
    const unsigned int nbSubTriggers = wcsimrootsuperevent->GetNumberOfSubEvents();

    //cout << "iEntry : " << iEntry << endl;
    //cout << "nbTrig : " << nbTriggers << endl;
    //cout << "nbSubTrig : " << nbSubTriggers << endl;
    
    for(long unsigned int iTrig = 0; iTrig < nbTriggers; iTrig++){
      WCSimRootTrigger *wcsimrootevent = wcsimrootsuperevent->GetTrigger(iTrig);

      int max=wcsimrootevent->GetNcherenkovhits();
      double sumQ = wcsimrootevent->GetSumQ();
      trigger_Q += sumQ;
      for (int i = 0; i<max; i++){
	WCSimRootCherenkovHit *chit = (WCSimRootCherenkovHit*)wcsimrootevent->GetCherenkovHits()->At(i);
	int totalPE = chit->GetTotalPe(1);
	for (int j = 0; j < totalPE; j++) {
	  PMT_hits->Fill(chit->GetTubeID());

	  if (j>100) {
	    PMT_hits100->Fill(chit->GetTubeID());
	  }
	}
	//WCSimRootCherenkovHit has methods GetTubeId(), GetTotalPe(int)
	PE->Fill(chit->GetTotalPe(1));
	//if (chit->GetTubeID() == 3482) std::cerr << "3482 " << chit->GetTotalPe(1) << std::endl;
	//if (chit->GetTubeID() == 2157) std::cerr << "2157 " << chit->GetTotalPe(1) << std::endl;
	//if (chit->GetTotalPe(1) > 10) std::cerr << chit->GetTubeID() << std::endl;
	//QvsTubeID->Fill(chit->GetTubeID(),chit->GetTotalPe(1));
      }
      
      // RAW HITS
      int rawMax = wcsimrootevent->GetNcherenkovhits();
      int totRawPE = 0;
      for (int i = 0; i < rawMax; i++){
	WCSimRootCherenkovHit *chit = (WCSimRootCherenkovHit*)wcsimrootevent->GetCherenkovHits()->At(i);
	//int tubeId = chit->GetTubeID();
	int tubeId = MAXPMT + chit->GetTubeID() - 1;
	WCSimRootPMT pmt = wcsimrootgeom->GetPMT(tubeId);
	double pmtX = pmt.GetPosition(0);
	double pmtY = pmt.GetPosition(1);
	double pmtZ = pmt.GetPosition(2);
	double pmtRad = sqrt((pow(pmtX, 2.0))+(pow(pmtY, 2.0)));
	double pmtTheta = std::atan(pmtY/pmtX);
	//	double l = sqrt((pow((0-pmtX),2)) + pow((-(HeightOD/2)-pmtY),2));
	double l = sqrt((pow((0-pmtX),2)) + pow((-(RadiusOD/2)-pmtY),2));
	double angle = 2*asin(l/7200);
	double length = angle*3600;
	if(chit->GetTotalPe(1) != 0) {
	  hPEByEvtsByPMT->Fill(chit->GetTotalPe(1));
	  if (pmt.GetCylLoc() == 5) {
	    double ODpolyX = pmtX;
	    double ODpolyY = pmtY + RadiusOD + HeightOD/2 + 100;
	    ODPEpoly->Fill(ODpolyX,ODpolyY,chit->GetTotalPe(1));
	  }
	  else if (pmt.GetCylLoc() == 3) {
	    double ODpolyX = pmtX;
	    double ODpolyY = -(HeightOD/2 +RadiusOD + pmtY + 100);
	    ODPEpoly->Fill(ODpolyX,ODpolyY,chit->GetTotalPe(1));
	  }
	  else if (pmt.GetCylLoc() == 4) {
	    if (pmtX>0) pmtTheta += (std::acos(-1)/2);
	    else pmtTheta -= (std::acos(-1)/2);
	    double pmtRadTheta = pmtRad*pmtTheta;
	    ODPEpoly->Fill(pmtRadTheta,pmtZ,chit->GetTotalPe(1));
	  }
	}
	totRawPE+=chit->GetTotalPe(1);
      } // END FOR RAW HITS

      hNbTubesHit->Fill(rawMax);
      hPEByEvts->Fill(totRawPE);

      // DIGI HITS
      int digiMax = wcsimrootevent->GetNcherenkovdigihits();
      double totDigiPE = 0;
      for (int i = 0; i < digiMax; i++){
	WCSimRootCherenkovDigiHit *cDigiHit =
	  (WCSimRootCherenkovDigiHit*)wcsimrootevent->GetCherenkovDigiHits()->At(i);
	//WCSimRootChernkovDigiHit has methods GetTubeId(), GetT(), GetQ()
	WCSimRootCherenkovHitTime *cHitTime =
	  (WCSimRootCherenkovHitTime*)wcsimrootevent->GetCherenkovHitTimes()->At(i);
	//WCSimRootCherenkovHitTime has methods GetTubeId(), GetTruetime()
	//int tubeId = cDigiHit->GetTubeId();
	int tubeId = MAXPMT + cDigiHit->GetTubeId() - 1;
	WCSimRootPMT pmt = wcsimrootgeom->GetPMT(tubeId);
	double pmtX = pmt.GetPosition(0);
	double pmtY = pmt.GetPosition(1);
	double pmtZ = pmt.GetPosition(2);
	//std::cout << "CylLoc " << pmt.GetCylLoc() << std::endl;
	//std::cout << "Tube hit: " << tubeId << std::endl;
	//std::cout << "Tube position: (" << pmtX << "," << pmtY << "," << pmtZ << ")" << std::endl;
	//std::cout << "Charge deposited : " << cDigiHit->GetQ() << std::endl;
	double testcharge = cDigiHit->GetQ();
	if (testcharge > 0.0) check_PMT_hits->Fill(cDigiHit->GetTubeId());
	if (testcharge < 0) std::cerr << "WHY DO YOU HATE ME" << std::endl;
	QvsTubeID->Fill(cDigiHit->GetTubeId(),testcharge);
	//ODpoly->Fill(pmtX,pmtY,testcharge);
	//if (testcharge > 25) std::cerr << tubeId << " : " << testcharge << " x " << pmtX << " y " << pmtY << " z " << pmtZ << std::endl;
	double pmtRad = sqrt((pow(pmtX, 2.0))+(pow(pmtY, 2.0)));
	double pmtTheta = std::atan(pmtY/pmtX);
	double l = sqrt((pow((0-pmtX),2)) + pow((-3600-pmtY),2));
	double angle = 2*asin(l/7200);
	double length = angle*3600;
	//std::cerr << "cylloc " << pmt.GetCylLoc() << std::endl;
	if (pmt.GetCylLoc() == 5) {
	  TopCapHits->Fill(pmtX,pmtY,testcharge);
	  double ODpolyX = pmtX;
	  double ODpolyY = pmtY + RadiusOD + HeightOD/2 + 100;
	  //std::cerr << "CylLoc3 " << ODpolyX << "," << ODpolyY << ", charge " << testcharge << std::endl;
	  ODpoly->Fill(ODpolyX,ODpolyY,testcharge);
	  //std::cerr << "cyloc 3, z pos " << pmtZ << std::endl;
	  check_total_charge += testcharge;
	}
	else if (pmt.GetCylLoc() == 3) {
	  BottomCapHits->Fill(pmtX,pmtY,testcharge);
	  double ODpolyX = pmtX;
	  double ODpolyY = -(HeightOD/2 +RadiusOD + pmtY + 100);
	  //std::cerr << "CylLoc5 " << ODpolyX << "," << ODpolyY << ", charge " << testcharge << std::endl;
	  ODpoly->Fill(ODpolyX,ODpolyY,testcharge);
	  //std::cerr << "cyloc 5, z pos " << pmtZ << std::endl;
	  check_total_charge += testcharge;
	}
	else if (pmt.GetCylLoc() == 4) {
	  /*double tube[3];
	    tube[0] = wcsimrootgeom->GetPMT(tubeId).GetPosition(0);
	    tube[1] = wcsimrootgeom->GetPMT(tubeId).GetPosition(1);
	    tube[2] = wcsimrootgeom->GetPMT(tubeId).GetPosition(2);
	    double l = sqrt( pow((0 - tube[0]),2) + pow((-RadiusOD - tube[1]),2));
	    double angle = 2*asin(l/(2*RadiusOD));
	    double length = angle*RadiusOD ;
	    if (tube[0]<0) length *= -1;
	    double pmtRadTheta = angle*length;
	    BarrelHits->Fill(pmtRadTheta,pmtZ,testcharge);
	    ODpoly->Fill(pmtRadTheta,pmtZ,testcharge);*/
	  if (pmtX>0) pmtTheta += (std::acos(-1)/2);
	  else pmtTheta -= (std::acos(-1)/2);
	  double pmtRadTheta = pmtRad*pmtTheta;
	  //std::cerr << "pmtTheta " << pmtTheta << " pmtRad " << pmtRad << std::endl;
	  //std::cerr << "filling: " << pmtRadTheta << ","<< pmtZ << "," << testcharge << std::endl;
	  BarrelHits->Fill(pmtRadTheta,pmtZ,testcharge);
	  //std::cerr << "CylLoc4 " << pmtRadTheta << "," << pmtZ << ", charge " << testcharge << std::endl;
	  ODpoly->Fill(pmtRadTheta,pmtZ,testcharge);
	  //if (pmtX<0) length*= -1;
	  //std::cerr << "pmtRadTheta : " << pmtRadTheta << ", length : " << length << std::endl;
	  //BarrelHits->Fill(length,pmtZ,testcharge);
	  check_total_charge += testcharge;
	}
	if(cDigiHit->GetQ() != 0) {
	  hPECollectedByEvtsByPMT->Fill(cDigiHit->GetQ());
	  //check_PMT_hits->Fill(cDigiHit->GetTubeId());
	}
	totDigiPE+=cDigiHit->GetQ();
      } // END FOR DIGI HITS

      hPECollectedByEvts->Fill(totDigiPE);
      //      std::cerr << "digiMax " << digiMax << std::endl;
      total_charge +=totDigiPE;
    } // END FOR iTRIG
    /*  for (Int_t i = 0; i < check_PMT_hits->GetNbinsX(); i++) {
      if (check_PMT_hits->GetBinContent(i)) num_hit_OD_PMTs += 1;
    }

    double coverage = (num_hit_OD_PMTs/num_OD_PMT)*100;
    coverage_i.push_back(coverage);
    */
  } // END FOR iENTRY
  std::cerr << "total_charge " << total_charge << std::endl;
  std::cerr << "check_total_charge " << check_total_charge << std::endl;
  std::cerr << "trigger_Q " << trigger_Q << std::endl;
  /*  double sum_coverage = 0;
  double mean_coverage_1_hit = -999;
  for (int i = 0; i < coverage_i.size(); i++) {
    sum_coverage += coverage_i.at(i);
  }
  mean_coverage_1_hit = sum_coverage/coverage_i.size();
  std::cerr << "Average coverage 1 hit: " << mean_coverage_1_hit << std::endl;
  delete check_PMT_hits;*/
  for (Int_t i = 0; i < check_PMT_hits->GetNbinsX(); i++) {
    if (check_PMT_hits->GetBinContent(i)) num_OD_PMTs_with_hits += 1;
  }

  std::cerr << "num_OD_PMTs_with_hits " << (num_OD_PMTs_with_hits) << ", num_OD_PMT " << num_OD_PMT << std::endl;
  std::cerr << "PMT COVERAGE : " << (num_OD_PMTs_with_hits/num_OD_PMT)*100 << std::endl;
  //////////////////////////////////////////
  // DRAWING ///////////////////////////////
  //////////////////////////////////////////

  /*  TCanvas *c1;

  c1 = new TCanvas("cPE","cPE",800,600);
  c1->Divide(2,2);
  c1->cd(1);
  gPad->SetLogy();
  hPEByEvtsByPMT->Draw("HIST");
  c1->cd(2);
  hPEByEvts->Draw(""); hPEByEvts->Fit("gaus");
  c1->cd(3);
  gPad->SetLogy();
  hPECollectedByEvtsByPMT->Draw("HIST");
  c1->cd(4);  
  hPECollectedByEvts->Draw("HIST");

  TF1 *fit = hPEByEvts->GetFunction("gaus");

  //  c1 = new TCanvas("cBarrelHits","cBarrelHits",800,600);
  //BarrelHits->Draw("BOX");
  
  c1 = new TCanvas("cNbTubesHit","cNbTubesHit",800,600);
  hNbTubesHit->Draw("HIST");

  cout << "Mean nb of tubes hit by events : " << hNbTubesHit->GetMean()
       << " +- " << hNbTubesHit->GetRMS() << endl;
  cout << "Mean raw PE by events by PMT : " << hPEByEvtsByPMT->GetMean()
       << " +- " << hPEByEvtsByPMT->GetRMS() << endl;  
  cout << "Mean PE collected by events by PMT : " << hPECollectedByEvtsByPMT->GetMean()
       << " +- " << hPECollectedByEvtsByPMT->GetRMS() << endl;
  cout << "Mean raw PE by events : " << hPEByEvts->GetMean()
       << " +- " << hPEByEvts->GetRMS() << endl;  
  // cout << "Mean PE collected by events : " << hPECollectedByEvts->GetMean()
  //      << " +- " << hPECollectedByEvts->GetRMS() << endl;
  // cout << "FIT : " << fit->GetParameter(1)
  //      << " +- " << fit->GetParError(1) << endl;

  */
  // TFile *output=NULL;


  /* TCanvas *cp = new TCanvas("cp","cp",1000,1500);
  //  ODpoly->GetZaxis()->SetRangeUser(1.0,1000);
  gPad->SetLogz(1);
  ODpoly->Draw("COLZ 0");
  cp->SetRightMargin(0.15);
  //cp->SaveAs("3inch_1000_out_half_poly.pdf");*/
  gStyle->SetPalette(55);
  /*  TCanvas *cpe = new TCanvas("cpe","cpe",1500,1500);
  //  ODpoly->GetZaxis()->SetRangeUser(1.0,1000);
  gPad->SetLogz(1);
  ODPEpoly->Draw("COLZ 0");
  cpe->SetRightMargin(0.15);
  //cpe->SaveAs("/home/pidcott/Desktop/photons_400nm_bare_1000_028c.pdf");
  */
  TCanvas *cp = new TCanvas("cp","cp",900,900);
    ODpoly->GetZaxis()->SetRangeUser(1,300);
    ODpoly->GetYaxis()->SetRangeUser(3500, 10500);
  //  ODpoly->GetYaxis()->SetRangeUser(-10500,-3500);
  ODpoly->GetXaxis()->SetRangeUser(-3500,3500);
  gPad->SetLogz(1);
  cp->SetRightMargin(0.12);
  //ODpoly->Draw("COLZ 0 AH");
  ODpoly->Draw("COLZ 0");
  //  cp->SaveAs("/home/pidcott/OD_calib/calib_devel/WCSim/tyvek_100ev_struts_mvLI.pdf");

  /*TCanvas *c2 = new TCanvas("c2","c2",1200,300);
  //BarrelHits->Scale(0.2);
  BarrelHits->Draw("COLZ2");
  c2->SetRightMargin(0.15);
  //c2->SaveAs("testingdraw.C");
  c2->SaveAs("3inch_100_barrel_in_half.pdf");

  /*Int_t colours[] = {0, 1, 2, 3, 4, 5, 6};
  gStyle->SetPalette((sizeof(colours)/sizeof(Int_t)), colours);
  Double_t levels[] = {0.0, 200, 400, 600, 800, 1000, 1200, 1400};
  TopCapHits->SetContour((sizeof(levels)/sizeof(Double_t)));*/
  /* TCanvas *c3 = new TCanvas("c3","c3",800,750);
  //TopCapHits->Scale(0.2);
  TopCapHits->Draw("COLZ2");
  c3->SetRightMargin(0.15);
  c3->SaveAs("3inch_100_top_in_half.pdf");

  TCanvas *c4 = new TCanvas("c4","c4",800,750);
  //  BottomCapHits->Scale(0.2);
  BottomCapHits->Draw("COLZ2");
  c4->SetRightMargin(0.15);
  c4->SaveAs("3inch_100_bottom_in_half.pdf");*/
  // TCanvas *c3 = new TCanvas("c3","c3",1200,750);
  // QvsNTubesHit->Draw();

  // TCanvas *c4 = new TCanvas("c4","c4",1200,750);
  // QvsNHits->Draw();

  //  TCanvas *c5 = new TCanvas("c5","c5",1200,750);
  //QvsTubeID->Draw("COLZ2");

  //TCanvas *c6 = new TCanvas("c6","c6",1200,750);
  //PMT_hits->Scale(0.001);
  //  PMT_hits->Draw();
  //TCanvas *c8 = new TCanvas("c8","c8",1200,750);
  //check_PMT_hits->Draw();

  //TCanvas *c7 = new TCanvas("c7","c7",1200,750);
  //PMT_hits100->Draw();

  // char outputName[1000];
  // sprintf(outputName,"PROCESSED/PROCESSED_%s",bname);
  // output = new TFile(outputName,"recreate");
  // hPEByEvtsByPMT->Write();
  // hPEByEvts->Write();
  // hPECollectedByEvtsByPMT->Write();
  // hPECollectedByEvts->Write();

} // END MACRO
