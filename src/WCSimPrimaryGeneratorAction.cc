#include "WCSimPrimaryGeneratorAction.hh"
#include "WCSimDetectorConstruction.hh"
#include "WCSimPrimaryGeneratorMessenger.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4GeneralParticleSource.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4OpticalPhoton.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"
#include "Randomize.hh"
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <TFile.h>
#include <G4RandomDirection.hh>

#include "G4Navigator.hh"
#include "G4TransportationManager.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

using std::vector;
using std::string;
using std::fstream;

vector<string> tokenize( string separators, string input );

inline vector<string> readInLine(fstream& inFile, int lineSize, char* inBuf)
{
  // Read in line break it up into tokens
  inFile.getline(inBuf,lineSize);
  return tokenize(" $\r", inBuf);
}

inline float atof( const string& s ) {return std::atof( s.c_str() );}
inline int   atoi( const string& s ) {return std::atoi( s.c_str() );}

WCSimPrimaryGeneratorAction::WCSimPrimaryGeneratorAction(
					  WCSimDetectorConstruction* myDC)
  :myDetector(myDC), vectorFileName("")
{
  //T. Akiri: Initialize GPS to allow for the laser use 
  MyGPS = new G4GeneralParticleSource();

  // Initialize to zero
  mode = 0;
  vtxvol = 0;
  vtx = G4ThreeVector(0.,0.,0.);
  nuEnergy = 0.;
  _counterRock=0; // counter for generated in Rock
  _counterCublic=0; // counter generated
  
  //---Set defaults. Do once at beginning of session.
  
  G4int n_particle = 1;
  particleGun = new G4ParticleGun(n_particle);
  particleGun->SetParticleEnergy(1.0*GeV);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.0));

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  particleGun->
    SetParticleDefinition(particleTable->FindParticle(particleName="mu+"));

  particleGun->
    SetParticlePosition(G4ThreeVector(0.*m,0.*m,0.*m));
    
  messenger = new WCSimPrimaryGeneratorMessenger(this);
  useMulineEvt = true;
  useGunEvt    = false;
  useLaserEvt  = false;
  useGPSEvt    = false;
  useCalibration = false;
  useFullInjectors = false;
  useColBarrel = false;
  useColBottom = false;
  useColTop = false;
  useTimeTop = false;
  useSatTop = false;
  useSatBarrel = false;
  useIWCDFullInjectors = false;

  particleGunBarrel = new G4ParticleGun();
  particleGunTop = new G4ParticleGun();
  particleGunBottom = new G4ParticleGun();

// Create the relevant histograms to generate muons
  // according to SuperK flux extrapolated at HyperK site
  std::fstream inputFileCosmics;
  G4String vectorFileNameCosmics;
  altCosmics = 2*myDC->GetWCIDHeight();
  G4cout << "altCosmics : " << altCosmics << G4endl;
  if (inputFileCosmics.is_open())
    inputFileCosmics.close();

  vectorFileNameCosmics = "MuonFlux-HyperK-ThetaPhi.dat";
  inputFileCosmics.open(vectorFileNameCosmics, std::fstream::in);

  if (!inputFileCosmics.is_open()) {
    G4cout << "Muon Vector file " << vectorFileNameCosmics << " not found" << G4endl;
  } else {
    G4cout << "Muon Vector file " << vectorFileNameCosmics << " found" << G4endl;
    string line;
    vector<string> token(1);

    double binCos, binPhi;
    double cosThetaMean, cosThetaMin, cosThetaMax;
    double phiMean, phiMin, phiMax;
    double flux;
    double Emean;

    hFluxCosmics = new TH2D("hFluxCosmics","HK Flux", 180,0,360,100,0,1);
    hFluxCosmics->GetXaxis()->SetTitle("#phi (deg)");
    hFluxCosmics->GetYaxis()->SetTitle("cos #theta");
    hEmeanCosmics = new TH2D("hEmeanCosmics","HK Flux", 180,0,360,100,0,1);
    hEmeanCosmics->GetXaxis()->SetTitle("#phi (deg)");
    hEmeanCosmics->GetYaxis()->SetTitle("cos #theta");

    while ( getline(inputFileCosmics,line) ){
      token = tokenize(" $", line);

      binCos=(atof(token[0]));
      binPhi=(atof(token[1]));
      cosThetaMean=(atof(token[2]));
      cosThetaMin=(atof(token[3]));
      cosThetaMax=(atof(token[4]));
      phiMean=(atof(token[5]));
      phiMin=(atof(token[6]));
      phiMax=(atof(token[7]));
      flux=(atof(token[8]));
      Emean=(atof(token[9]));

      hFluxCosmics->SetBinContent(binPhi,binCos,flux);
      hEmeanCosmics->SetBinContent(binPhi,binCos,Emean);
    }

    TFile *file = new TFile("flux.root","RECREATE");
    hFluxCosmics->Write();
    hEmeanCosmics->Write();
    file->Close();

  }
}

WCSimPrimaryGeneratorAction::~WCSimPrimaryGeneratorAction()
{
  if (IsGeneratingVertexInRock()){
    G4cout << "Fraction of Rock volume is : " << G4endl;
      G4cout << " Random number generated in Rock / in Cublic = " 
             << _counterRock << "/" << _counterCublic 
             << " = " << _counterRock/(G4double)_counterCublic << G4endl;
  }
  inputFile.close();
  delete particleGun;
  delete MyGPS;   //T. Akiri: Delete the GPS variable
  delete messenger;
}

void WCSimPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{

  // We will need a particle table
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();

  // Temporary kludge to turn on/off vector text format 

  G4bool useNuanceTextFormat = true;


  // Do for every event

  if (useMulineEvt)
  { 

    if ( !inputFile.is_open() )
    {
      G4cout << "Set a vector file using the command /mygen/vecfile name"
	     << G4endl;
      exit(-1);
    }

    //
    // Documentation describing the nuance text format can be found here: 
    // http://neutrino.phy.duke.edu/nuance-format/
    //
    // The format must be strictly adhered to for it to be processed correctly.
    // The lines and their meanings from begin through info are fixed, and then
    // a variable number of tracks may follow.
    //
    if (useNuanceTextFormat)
      {
	const int lineSize=100;
	char      inBuf[lineSize];
	vector<string> token(1);
	
	token = readInLine(inputFile, lineSize, inBuf);
	  
        if (token.size() == 0) 
	  {
	    G4cout << "end of nuance vector file!" << G4endl;
	  }
	else if (token[0] != "begin")
	  {
	    G4cout << "unexpected line begins with " << token[0] << G4endl;
	  }
	else   // normal parsing begins here
	  {
	    // Read the nuance line (ignore value now)

	    token = readInLine(inputFile, lineSize, inBuf);
	    mode = atoi(token[1]);

	    // Read the Vertex line
	    token = readInLine(inputFile, lineSize, inBuf);
	    vtx = G4ThreeVector(atof(token[1])*cm,
				atof(token[2])*cm,
				atof(token[3])*cm);
	    
            // true : Generate vertex in Rock , false : Generate vertex in WC tank
            SetGenerateVertexInRock(false);

	    // Next we read the incoming neutrino and target
	    
	    // First, the neutrino line

	    token=readInLine(inputFile, lineSize, inBuf);
	    beampdg = atoi(token[1]);
	    beamenergy = atof(token[2])*MeV;
	    beamdir = G4ThreeVector(atof(token[3]),
				    atof(token[4]),
				    atof(token[5]));

	    // Now read the target line

	    token=readInLine(inputFile, lineSize, inBuf);
	    targetpdg = atoi(token[1]);
	    targetenergy = atof(token[2])*MeV;
	    targetdir = G4ThreeVector(atof(token[3]),
				      atof(token[4]),
				      atof(token[5]));

	    // Read the info line, basically a dummy
	    token=readInLine(inputFile, lineSize, inBuf);
	    G4cout << "Vector File Record Number " << token[2] << G4endl;
            vecRecNumber = atoi(token[2]);
	    
	    // Now read the outgoing particles
	    // These we will simulate.


	    while ( token=readInLine(inputFile, lineSize, inBuf),
		    token[0] == "track" )
	      {
		// We are only interested in the particles
		// that leave the nucleus, tagged by "0"


		if ( token[6] == "0")
		  {
		    G4int pdgid = atoi(token[1]);
		    G4double energy = atof(token[2])*MeV;
		    G4ThreeVector dir = G4ThreeVector(atof(token[3]),
						      atof(token[4]),
						      atof(token[5]));
		    particleGun->
		      SetParticleDefinition(particleTable->
					    FindParticle(pdgid));
		    G4double mass = 
		      particleGun->GetParticleDefinition()->GetPDGMass();

		    G4double ekin = energy - mass;

		    particleGun->SetParticleEnergy(ekin);
		    //G4cout << "Particle: " << pdgid << " KE: " << ekin << G4endl;
		    particleGun->SetParticlePosition(vtx);
		    particleGun->SetParticleMomentumDirection(dir);
		    particleGun->GeneratePrimaryVertex(anEvent);
		  }
	      }
	  }
      }
    else 
      {    // old muline format  
	inputFile >> nuEnergy >> energy >> xPos >> yPos >> zPos 
		  >> xDir >> yDir >> zDir;
	
	G4double random_z = ((myDetector->GetWaterTubePosition())
			     - .5*(myDetector->GetWaterTubeLength()) 
			     + 1.*m + 15.0*m*G4UniformRand())/m;
	zPos = random_z;
	G4ThreeVector vtx = G4ThreeVector(xPos, yPos, random_z);
	G4ThreeVector dir = G4ThreeVector(xDir,yDir,zDir);

	particleGun->SetParticleEnergy(energy*MeV);
	particleGun->SetParticlePosition(vtx);
	particleGun->SetParticleMomentumDirection(dir);
	particleGun->GeneratePrimaryVertex(anEvent);
      }
  }

  else if (useGunEvt)
  {      // manual gun operation
    particleGun->GeneratePrimaryVertex(anEvent);

    G4ThreeVector P  =anEvent->GetPrimaryVertex()->GetPrimary()->GetMomentum();
    G4ThreeVector vtx=anEvent->GetPrimaryVertex()->GetPosition();
    G4double m       =anEvent->GetPrimaryVertex()->GetPrimary()->GetMass();
    G4int pdg        =anEvent->GetPrimaryVertex()->GetPrimary()->GetPDGcode();

    G4ThreeVector dir  = P.unit();
    G4double E         = std::sqrt((P.dot(P))+(m*m));

//     particleGun->SetParticleEnergy(E);
//     particleGun->SetParticlePosition(vtx);
//     particleGun->SetParticleMomentumDirection(dir);

    SetVtx(vtx);
    SetBeamEnergy(E);
    SetBeamDir(dir);
    SetBeamPDG(pdg);
  }
  else if (useLaserEvt)
    {
      //T. Akiri: Create the GPS LASER event
      MyGPS->GeneratePrimaryVertex(anEvent);
      
      G4ThreeVector P   =anEvent->GetPrimaryVertex()->GetPrimary()->GetMomentum();
      G4ThreeVector vtx =anEvent->GetPrimaryVertex()->GetPosition();
      G4int pdg         =anEvent->GetPrimaryVertex()->GetPrimary()->GetPDGcode();
      
      G4ThreeVector dir  = P.unit();
      G4double E         = std::sqrt((P.dot(P)));
      
      SetVtx(vtx);
      SetBeamEnergy(E);
      SetBeamDir(dir);
      SetBeamPDG(pdg);
    }
  else if (useGPSEvt)
    {
      MyGPS->GeneratePrimaryVertex(anEvent);
      
      G4ThreeVector P   =anEvent->GetPrimaryVertex()->GetPrimary()->GetMomentum();
      G4ThreeVector vtx =anEvent->GetPrimaryVertex()->GetPosition();
      G4double m        =anEvent->GetPrimaryVertex()->GetPrimary()->GetMass();
      G4int pdg         =anEvent->GetPrimaryVertex()->GetPrimary()->GetPDGcode();
      
      G4ThreeVector dir  = P.unit();
      G4double E         = std::sqrt((P.dot(P))+(m*m));
      
      SetVtx(vtx);
      SetBeamEnergy(E);
      SetBeamDir(dir);
      SetBeamPDG(pdg);
    }
  else if(useCosmics){

    //////////////////
    // DEBUG PRINTS
    G4cout << G4endl;
    G4cout << "COSMYMYMATICS" << G4endl;
    G4cout << "#############" << G4endl;
    //////////////////

    double phiMuon, cosThetaMuon;
    energy = 0;
    while((int)(energy) == 0){
      hFluxCosmics->GetRandom2(phiMuon,cosThetaMuon);
      energy = hEmeanCosmics->GetBinContent(hFluxCosmics->GetBin(phiMuon,cosThetaMuon))*GeV;
    }
    

    G4ThreeVector dir(0,0,0);
    dir.setRThetaPhi(-1,acos(cosThetaMuon),phiMuon);
    G4ThreeVector vtx(0,0,0);
    vtx = -dir;
    vtx.setR(altCosmics);

    int pdgid = 13; // MUON
    particleGun->SetParticleDefinition(particleTable->FindParticle(pdgid));
    G4double mass =particleGun->GetParticleDefinition()->GetPDGMass();
    G4double ekin = energy - mass;

    //////////////////
    // DEBUG PRINTS
    G4cout << G4endl;
    G4cout << "Generated at position : " << vtx.getX()/m << "m "
           << vtx.getY()/m << "m "
           << vtx.getZ()/m << "m " << G4endl;
    G4cout << "phi : " << phiMuon << " cosTheta : " << cosThetaMuon << G4endl;
    G4cout << "E : " << energy/GeV << " GeV" << G4endl;
    G4cout << G4endl;
    //////////////////

    SetVtx(vtx);
    SetBeamEnergy(energy);
    SetBeamDir(dir);
    SetBeamPDG(pdgid);

    particleGun->SetParticleEnergy(ekin);
    particleGun->SetParticlePosition(vtx);
    particleGun->SetParticleMomentumDirection(dir);
    particleGun->GeneratePrimaryVertex(anEvent);

  }
  else if(useCalibration) {
    if(useFullInjectors) {

      //	std::cerr << "Running all injectors with " << c_particle << " photons and a half angle of " << ledtheta << " degrees" << std::endl;
      //------------------------------- Light array ------------------------     
      particleGunBarrel->SetParticleDefinition(G4OpticalPhoton::Definition());
      particleGunBarrel->SetNumberOfParticles(c_particle);
      particleGunTop->SetParticleDefinition(G4OpticalPhoton::Definition());
      particleGunTop->SetNumberOfParticles(c_particle);
      particleGunBottom->SetParticleDefinition(G4OpticalPhoton::Definition());
      particleGunBottom->SetNumberOfParticles(c_particle);

      G4double energytemp = 3.0996 * eV; // 400 nm
      double maxr = 3361*cm; // NEW HK GEOM facing out
      //double maxr = 3451*cm; // NEW HK GEOM facing in
      double rad1 = 800*cm;
      double rad2 = 1800*cm;
      double rad3 = 2800*cm;
      double xposb = -9999*cm;
      double yposb = -9999*cm;
      double zposb = -9999*cm;


      for(int ip = 0; ip < 5; ip++) {
	if (ip==0) zposb = 2700*cm;
	if (ip==1) zposb = 1350*cm;
	if (ip==2) zposb = 0*cm;
	if (ip==3) zposb = -1350*cm;
	if (ip==4) zposb = -2700*cm;
	for(int theta = 0; theta < 16; theta++) {
	  if (ip == 1 || ip == 3) {
	    xposb = maxr*(std::cos((theta+0.0)*(std::acos(-1))/8));
	    yposb = maxr*(std::sin((theta+0.0)*(std::acos(-1))/8));
	  }
	  else {
	    xposb = maxr*(std::cos((theta+0.5)*(std::acos(-1))/8));
	    yposb = maxr*(std::sin((theta+0.5)*(std::acos(-1))/8));
	  }
	  G4ThreeVector postempb(xposb,yposb,zposb);
	  particleGunBarrel->SetParticlePosition(postempb);
	  G4double rndm = G4UniformRand();
	  G4double rndm2 = G4UniformRand();
	  G4double costheta = 1- rndm*(1 - std::cos(ledtheta*deg));
	  G4double sintheta = std::sqrt(1. - costheta*costheta);
	  G4double phi = rndm2*2*(std::acos(-1));
	  G4double sinphi = std::sin(phi);
	  G4double cosphi = std::cos(phi);
	  G4double px = sintheta*cosphi;
	  G4double py = sintheta*sinphi;
	  G4double pz = costheta;
	  G4ThreeVector XPrime = G4ThreeVector(0.,0.,1.);
	  //G4ThreeVector YPrime = G4ThreeVector(-yposb,xposb,0.); // Facing in
	  G4ThreeVector YPrime = G4ThreeVector(yposb,-xposb,0.); // Facing out
	  //G4ThreeVector ZPrime = G4ThreeVector(-xposb,-yposb,0.); // Facing in
	  G4ThreeVector ZPrime = G4ThreeVector(xposb,yposb,0.); // Facing out
	  G4double finx = (px*XPrime.x()) + (py*YPrime.x()) + (pz*ZPrime.x());
	  G4double finy = (px*XPrime.y()) + (py*YPrime.y()) + (pz*ZPrime.y());
	  G4double finz = (px*XPrime.z()) + (py*YPrime.z()) + (pz*ZPrime.z());
	  //std::cerr << "finx " << finx << ", finy " << finy << ", finz " << finz << std::endl;
	  G4double time = G4RandFlat::shoot(20.0,40.0);
	  particleGunBarrel->SetParticleTime(time);
	  G4ThreeVector momtempb(finx,finy,finz);
	  particleGunBarrel->SetParticleMomentumDirection(momtempb);
	  particleGunBarrel->SetParticleEnergy(energytemp);
	  particleGunBarrel->SetParticlePolarization(G4RandomDirection());
	  particleGunBarrel->GeneratePrimaryVertex(anEvent);
	  SetVtx(postempb);
	  SetBeamEnergy(energytemp);
	  SetBeamDir(momtempb.unit());
	  SetBeamPDG(0);
	}
      }
      for (int i = 0; i < 3; i++) {
	//double zpos = 3650*cm; // top facing down
	double zpos = -3460*cm; // bottom facing down
	double xpos = 0.;
	double ypos = 0.;
	//G4ThreeVector postemp(0, 0, zpos);
	if (i==0) {
	  for (int angfrac1 = 0; angfrac1 < 3; angfrac1++) {
	    //	    if (angfrac1 % 2 == 0) continue;
	    double xpos = rad1*(std::cos((angfrac1*8 + 3)*(std::acos(-1)/12)));
	    double ypos = rad1*(std::sin((angfrac1*8 + 3)*(std::acos(-1)/12)));
	    G4ThreeVector postemp(xpos,ypos,zpos);
	    particleGunBottom->SetParticlePosition(postemp);
	    //particleGun->SetNumberOfParticles(5000);
	    //	    G4double energytemp = 2.505 * eV;
	    //G4double energytemp = 3.0996 * eV; // 400 nm
	    G4double rndm = G4UniformRand();
	    G4double rndm2 = G4UniformRand();
	    G4double costheta = 1- rndm*(1 - std::cos(ledtheta*deg));
	    G4double sintheta = std::sqrt(1. - costheta*costheta);
	    G4double phi = rndm2*2*(std::acos(-1));
	    G4double sinphi = std::sin(phi);
	    G4double cosphi = std::cos(phi);
	    G4double px = sintheta*cosphi;
	    G4double py = sintheta*sinphi;
	    G4double pz = -costheta; // Facing down
	    //G4double pz = -1;
	    G4ThreeVector momtemp(px,py,pz);
	    particleGunBottom->SetParticleMomentumDirection(momtemp);
	    particleGunBottom->SetParticleEnergy(energytemp);
	    particleGunBottom->SetParticlePolarization(G4RandomDirection());
	    particleGunBottom->GeneratePrimaryVertex(anEvent);

	    SetVtx(postemp);
	    SetBeamEnergy(energytemp);
	    SetBeamDir(momtemp.unit());
	    SetBeamPDG(0);
	  }
	}
	if (i==1) {
	  for (int angfrac2 = 0; angfrac2 < 6; angfrac2++) {
	    //	    if (angfrac2 % 4 == 0) continue;
	    double xpos = rad2*(std::cos((angfrac2*4 + 1)*(std::acos(-1)/12)));
	    double ypos = rad2*(std::sin((angfrac2*4 + 1)*(std::acos(-1)/12)));
	    G4ThreeVector postemp(xpos,ypos,zpos);
	    particleGunBottom->SetParticlePosition(postemp);
	    //particleGun->SetNumberOfParticles(5000);
	    //G4double energytemp = 2.505 * eV;
	    //G4double energytemp = 3.0996 * eV; // 400 nm
	    G4double rndm = G4UniformRand();
	    G4double rndm2 = G4UniformRand();
	    G4double costheta = 1- rndm*(1 - std::cos(ledtheta*deg));
	    G4double sintheta = std::sqrt(1. - costheta*costheta);
	    G4double phi = rndm2*2*(std::acos(-1));
	    G4double sinphi = std::sin(phi);
	    G4double cosphi = std::cos(phi);
	    G4double px = sintheta*cosphi;
	    G4double py = sintheta*sinphi;
	    G4double pz = -costheta;
	    //G4double pz = -1;
	    G4ThreeVector momtemp(px,py,pz);
	    particleGunBottom->SetParticleMomentumDirection(momtemp);
	    particleGunBottom->SetParticleEnergy(energytemp);
	    particleGunBottom->SetParticlePolarization(G4RandomDirection());
	    particleGunBottom->GeneratePrimaryVertex(anEvent);

	    SetVtx(postemp);
	    SetBeamEnergy(energytemp);
	    SetBeamDir(momtemp.unit());
	    SetBeamPDG(0);
	  }
	}
	if (i==2) {
	  for (int angfrac3 = 0; angfrac3 < 12; angfrac3++) {
	    //	    	      if (angfrac3 % 8 == 0) continue;
	    double xpos = rad3*(std::cos((angfrac3)*(std::acos(-1)/6)));
	    double ypos = rad3*(std::sin((angfrac3)*(std::acos(-1)/6)));
	    G4ThreeVector postemp(xpos,ypos,zpos);
	    particleGunBottom->SetParticlePosition(postemp);
	    //particleGun->SetNumberOfParticles(5000);
	    //G4double energytemp = 2.505 * eV;
	    //	    G4double energytemp = 3.0996 * eV; // 400 nm
	    G4double rndm = G4UniformRand();
	    G4double rndm2 = G4UniformRand();
	    G4double costheta = 1- rndm*(1 - std::cos(ledtheta*deg));
	    G4double sintheta = std::sqrt(1. - costheta*costheta);
	    G4double phi = rndm2*2*(std::acos(-1));
	    G4double sinphi = std::sin(phi);
	    G4double cosphi = std::cos(phi);
	    G4double px = sintheta*cosphi;
	    G4double py = sintheta*sinphi;
	    G4double pz = -costheta;
	    //G4double pz = -1;
	    G4ThreeVector momtemp(px,py,pz);
	    particleGunBottom->SetParticleMomentumDirection(momtemp);
	    particleGunBottom->SetParticleEnergy(energytemp);
	    particleGunBottom->SetParticlePolarization(G4RandomDirection());
	    particleGunBottom->GeneratePrimaryVertex(anEvent);

	    SetVtx(postemp);
	    SetBeamEnergy(energytemp);
	    SetBeamDir(momtemp.unit());
	    SetBeamPDG(0);
	  }
	}
      }
      for (int i = 0; i < 3; i++) {
	double zpos = 3460*cm; // top facing up
	//double zpos = 3650*cm; // top facing down
	double xpos = 0.;
	double ypos = 0.;
	//G4ThreeVector postemp(0, 0, zpos);
	if (i==0) {
	  for (int angfrac1 = 0; angfrac1 < 3; angfrac1++) {
	    //	    	      if (angfrac1 % 2 == 0) continue;
	    double xpos = rad1*(std::cos((angfrac1*8 + 3)*(std::acos(-1)/12)));
	    double ypos = rad1*(std::sin((angfrac1*8 + 3)*(std::acos(-1)/12)));
	    G4ThreeVector postemp(xpos,ypos,zpos);
	    particleGunTop->SetParticlePosition(postemp);
	    //particleGun->SetNumberOfParticles(5000);
	    //G4double energytemp = 2.505 * eV;
	    //G4double energytemp = 3.0996 * eV; // 400 nm
	    G4double rndm = G4UniformRand();
	    G4double rndm2 = G4UniformRand();
	    G4double costheta = 1- rndm*(1 - std::cos(ledtheta*deg));
	    G4double sintheta = std::sqrt(1. - costheta*costheta);
	    G4double phi = rndm2*2*(std::acos(-1));
	    G4double sinphi = std::sin(phi);
	    G4double cosphi = std::cos(phi);
	    G4double px = sintheta*cosphi;
	    G4double py = sintheta*sinphi;
	    G4double pz = costheta; // Facing up
	    //G4double pz = -costheta; // Facing down
	    G4ThreeVector momtemp(px,py,pz);
	    particleGunTop->SetParticleMomentumDirection(momtemp);
	    particleGunTop->SetParticleEnergy(energytemp);
	    particleGunTop->SetParticlePolarization(G4RandomDirection());
	    particleGunTop->GeneratePrimaryVertex(anEvent);

	    SetVtx(postemp);
	    SetBeamEnergy(energytemp);
	    SetBeamDir(momtemp.unit());
	    SetBeamPDG(0);
	  }
	}
	if (i==1) {
	  for (int angfrac2 = 0; angfrac2 < 6; angfrac2++) {
	    //	    	      if (angfrac2 % 4 == 0) continue;
	    double xpos = rad2*(std::cos((4*angfrac2 + 1)*(std::acos(-1)/12)));
	    double ypos = rad2*(std::sin((4*angfrac2 + 1)*(std::acos(-1)/12)));
	    G4ThreeVector postemp(xpos,ypos,zpos);
	    particleGunTop->SetParticlePosition(postemp);
	    //particleGun->SetNumberOfParticles(5000);
	    //G4double energytemp = 2.505 * eV;
	    //	    G4double energytemp = 3.0996 * eV; // 400 nm
	    G4double rndm = G4UniformRand();
	    G4double rndm2 = G4UniformRand();
	    G4double costheta = 1- rndm*(1 - std::cos(ledtheta*deg));
	    G4double sintheta = std::sqrt(1. - costheta*costheta);
	    G4double phi = rndm2*2*(std::acos(-1));
	    G4double sinphi = std::sin(phi);
	    G4double cosphi = std::cos(phi);
	    G4double px = sintheta*cosphi;
	    G4double py = sintheta*sinphi;
	    G4double pz = costheta;
	    //G4double pz = -1;
	    G4ThreeVector momtemp(px,py,pz);
	    particleGunTop->SetParticleMomentumDirection(momtemp);
	    particleGunTop->SetParticleEnergy(energytemp);
	    particleGunTop->SetParticlePolarization(G4RandomDirection());
	    particleGunTop->GeneratePrimaryVertex(anEvent);

	    SetVtx(postemp);
	    SetBeamEnergy(energytemp);
	    SetBeamDir(momtemp.unit());
	    SetBeamPDG(0);
	  }
	}
	if (i==2) {
	  for (int angfrac3 = 0; angfrac3 < 12; angfrac3++) {
	    //	    	      if (angfrac3 % 8 == 0) continue;
	    double xpos = rad3*(std::cos((angfrac3)*(std::acos(-1)/6)));
	    double ypos = rad3*(std::sin((angfrac3)*(std::acos(-1)/6)));
	    G4ThreeVector postemp(xpos,ypos,zpos);
	    particleGunTop->SetParticlePosition(postemp);
	    //particleGun->SetNumberOfParticles(5000);
	    //	    G4double energytemp = 2.505 * eV;
	    //G4double energytemp = 3.0996 * eV; // 400 nm
	    G4double rndm = G4UniformRand();
	    G4double rndm2 = G4UniformRand();
	    G4double costheta = 1- rndm*(1 - std::cos(ledtheta*deg));
	    G4double sintheta = std::sqrt(1. - costheta*costheta);
	    G4double phi = rndm2*2*(std::acos(-1));
	    G4double sinphi = std::sin(phi);
	    G4double cosphi = std::cos(phi);
	    G4double px = sintheta*cosphi;
	    G4double py = sintheta*sinphi;
	    G4double pz = costheta;
	    //G4double pz = -1;
	    G4ThreeVector momtemp(px,py,pz);
	    particleGunTop->SetParticleMomentumDirection(momtemp);
	    particleGunTop->SetParticleEnergy(energytemp);
	    particleGunTop->SetParticlePolarization(G4RandomDirection());
	    particleGunTop->GeneratePrimaryVertex(anEvent);

	    SetVtx(postemp);
	    SetBeamEnergy(energytemp);
	    SetBeamDir(momtemp.unit());
	    SetBeamPDG(0);
	  }
	}
      }
    } // ------------------------------ end of light array ------------------------------

    if(useSatTop) {

      particleGunTop->SetParticleDefinition(G4OpticalPhoton::Definition());
      particleGunTop->SetNumberOfParticles(c_particle);

      G4double energytemp = 3.0996 * eV; // 400 nm
      double rad1 = 800*cm;
      double zpos = 3460*cm; // top facing up
      //double zpos = 3650*cm; // top facing down
      double xpos = rad1*(std::cos(3*(std::acos(-1)/12)));
      double ypos = rad1*(std::sin(3*(std::acos(-1)/12)));
      G4ThreeVector postemp(xpos,ypos,zpos);
      particleGunTop->SetParticlePosition(postemp);
      G4double rndm = G4UniformRand();
      G4double rndm2 = G4UniformRand();
      G4double costheta = 1- rndm*(1 - std::cos(ledtheta*deg));
      G4double sintheta = std::sqrt(1. - costheta*costheta);
      G4double phi = rndm2*2*(std::acos(-1));
      G4double sinphi = std::sin(phi);
      G4double cosphi = std::cos(phi);
      G4double px = sintheta*cosphi;
      G4double py = sintheta*sinphi;
      G4double pz = costheta; // Facing up
      //G4double pz = -costheta; // Facing down
      G4double time = G4RandFlat::shoot(20.0,40.0);
      particleGunTop->SetParticleTime(time);
      G4ThreeVector momtemp(px,py,pz);
      std::cerr << "xpos " << xpos << ", ypos " << ypos << ", zpos " << zpos << std::endl;
      std::cerr << "finx " << px << ", finy " << py << ", finz " << pz << std::endl;
      particleGunTop->SetParticleMomentumDirection(momtemp);
      particleGunTop->SetParticleEnergy(energytemp);
      particleGunTop->SetParticlePolarization(G4RandomDirection());
      particleGunTop->GeneratePrimaryVertex(anEvent);

      SetVtx(postemp);
      SetBeamEnergy(energytemp);
      SetBeamDir(momtemp.unit());
      SetBeamPDG(0);
    }
    
    if(useSatBarrel) {
      
      particleGunBarrel->SetParticleDefinition(G4OpticalPhoton::Definition());
      particleGunBarrel->SetNumberOfParticles(c_particle);

      G4double energytemp = 3.0996 * eV; // 400 nm
      double maxr = 3361*cm; // NEW HK GEOM facing out
      //double maxr = 3451*cm; // NEW HK GEOM facing in
      double xposb = -9999*cm;
      double yposb = -9999*cm;
      double zposb = 0*cm;
      double theta = 0;
      xposb = maxr*(std::cos((theta+0.5)*(std::acos(-1))/8));
      yposb = maxr*(std::sin((theta+0.5)*(std::acos(-1))/8));
      G4ThreeVector postempb(xposb,yposb,zposb);
      particleGunBarrel->SetParticlePosition(postempb);
      G4double rndm = G4UniformRand();
      G4double rndm2 = G4UniformRand();
      G4double costheta = 1- rndm*(1 - std::cos(ledtheta*deg));
      G4double sintheta = std::sqrt(1. - costheta*costheta);
      G4double phi = rndm2*2*(std::acos(-1));
      G4double sinphi = std::sin(phi);
      G4double cosphi = std::cos(phi);
      G4double px = sintheta*cosphi;
      G4double py = sintheta*sinphi;
      G4double pz = costheta;
      G4ThreeVector XPrime = G4ThreeVector(0.,0.,1.);
      //G4ThreeVector YPrime = G4ThreeVector(-yposb,xposb,0.); // Facing in
      G4ThreeVector YPrime = G4ThreeVector(yposb,-xposb,0.); // Facing out
      //G4ThreeVector ZPrime = G4ThreeVector(-xposb,-yposb,0.); // Facing in
      G4ThreeVector ZPrime = G4ThreeVector(xposb,yposb,0.); // Facing out
      G4double finx = (px*XPrime.x()) + (py*YPrime.x()) + (pz*ZPrime.x());
      G4double finy = (px*XPrime.y()) + (py*YPrime.y()) + (pz*ZPrime.y());
      G4double finz = (px*XPrime.z()) + (py*YPrime.z()) + (pz*ZPrime.z());
      //std::cerr << "finx " << finx << ", finy " << finy << ", finz " << finz << std::endl;
      G4double time = G4RandFlat::shoot(20.0,40.0);
      particleGunBarrel->SetParticleTime(time);
      G4ThreeVector momtempb(finx,finy,finz);
      particleGunBarrel->SetParticleMomentumDirection(momtempb);
      particleGunBarrel->SetParticleEnergy(energytemp);
      particleGunBarrel->SetParticlePolarization(G4RandomDirection());
      particleGunBarrel->GeneratePrimaryVertex(anEvent);
      SetVtx(postempb);
      SetBeamEnergy(energytemp);
      SetBeamDir(momtempb.unit());
      SetBeamPDG(0);
    }

    else if(useColBarrel) {
      //std::cerr << "Running barrel with " << c_particle << " photons and a half angle of " << ledtheta << " degrees" << std::endl;
      // ------------------------- collimated barrel --------------------
      particleGunBarrel->SetParticleDefinition(G4OpticalPhoton::Definition());
      particleGunBarrel->SetNumberOfParticles(c_particle);
      double xposb = 3361*cm;
      double yposb = 0*cm;
      double zposb = 3460*cm;
      G4ThreeVector postemp(xposb,yposb,zposb);
      particleGunBarrel->SetParticlePosition(postemp);
      G4double energytemp = 3.0996 * eV; // 400 nm
      G4double rndm = G4UniformRand();
      G4double rndm2 = G4UniformRand();
      G4double costheta = 1- rndm*(1 - std::cos(ledtheta*deg));
      G4double sintheta = std::sqrt(1. - costheta*costheta);
      G4double phi = rndm2*2*(std::acos(-1));
      G4double sinphi = std::sin(phi);
      G4double cosphi = std::cos(phi);
      G4double px = sintheta*cosphi;
      G4double py = sintheta*sinphi;
      G4double pz = -costheta;
      G4ThreeVector momtemp(px,py,pz);
      particleGunBarrel->SetParticleMomentumDirection(momtemp);
      particleGunBarrel->SetParticleEnergy(energytemp);
      particleGunBarrel->SetParticlePolarization(G4RandomDirection());
      particleGunBarrel->GeneratePrimaryVertex(anEvent);
      //std::cerr << "Position x: " << particleGunBarrel->GetParticlePosition().x() << " y: " << particleGunBarrel->GetParticlePosition().y() << " z: " << particleGunBarrel->GetParticlePosition().z() << std::endl;
      //std::cerr << "Direction x: " << particleGunBarrel->GetParticleMomentumDirection().x() << " y: " << particleGunBarrel->GetParticleMomentumDirection().y() << " z: " << particleGunBarrel->GetParticleMomentumDirection().z() << std::endl;
      SetVtx(postemp);
      SetBeamEnergy(energytemp);
      SetBeamDir(momtemp.unit());
      SetBeamPDG(0);
      
      // ------------------------ end of collimated barrel -----------------

    }
    else if(useColBottom) {
      // ---------------- single collimated bottom cap pointing across -----

      //	std::cerr << "Running bottom with " << c_particle << " photons and a half angle of " << ledtheta << " degrees" << std::endl;
      particleGunBottom->SetParticleDefinition(G4OpticalPhoton::Definition());
      particleGunBottom->SetNumberOfParticles(c_particle);
      double xpos = 3361*cm;
      double ypos = 0*cm;
      double zpos = -3460*cm; // zpos of start of OD + PMT exposed height
      //double ledtheta = 2*deg;
      G4ThreeVector postemp(xpos,ypos,zpos);
      particleGunBottom->SetParticlePosition(postemp);
      G4double energytemp = 3.0996 * eV; // 400 nm
      G4double rndm = G4UniformRand();
      G4double rndm2 = G4UniformRand();
      G4double costheta = 1- rndm*(1 - std::cos(ledtheta*deg));
      G4double sintheta = std::sqrt(1. - costheta*costheta);
      G4double phi = rndm2*2*(std::acos(-1));
      G4double sinphi = std::sin(phi);
      G4double cosphi = std::cos(phi);
      G4double px = sintheta*cosphi;
      G4double py = sintheta*sinphi;
      G4double pz = costheta;
      G4ThreeVector XPrime = G4ThreeVector(0.,0.,1.);
      G4ThreeVector YPrime = G4ThreeVector(-ypos,xpos,0.);
      G4ThreeVector ZPrime = G4ThreeVector(-xpos,-ypos,0.);
      G4double finx = (px*XPrime.x()) + (py*YPrime.x()) + (pz*ZPrime.x());
      G4double finy = (px*XPrime.y()) + (py*YPrime.y()) + (pz*ZPrime.y());
      G4double finz = (px*XPrime.z()) + (py*YPrime.z()) + (pz*ZPrime.z());
      G4ThreeVector momtemp(finx,finy,finz);
      particleGunBottom->SetParticleMomentumDirection(momtemp);
      particleGunBottom->SetParticleEnergy(energytemp);
      particleGunBottom->SetParticlePolarization(G4RandomDirection());
      particleGunBottom->GeneratePrimaryVertex(anEvent);

      SetVtx(postemp);
      SetBeamEnergy(energytemp);
      SetBeamDir(momtemp.unit());
      SetBeamPDG(0);
      
      // ------------ end of single collimated bottom cap pointing across -----

    }
    else if(useColTop) {

      //	std::cerr << "Running top with " << c_particle << " photons and a half angle of " << ledtheta << " degrees" << std::endl;
      particleGunTop->SetParticleDefinition(G4OpticalPhoton::Definition());
      particleGunTop->SetNumberOfParticles(c_particle);
      double xpos = 3361*cm;
      double ypos = 0*cm;
      double zpos = 3460*cm; // zpos of start of OD + PMT exposed height
      //double ledtheta = 2*deg;
      G4ThreeVector postemp(xpos,ypos,zpos);
      particleGunTop->SetParticlePosition(postemp);
      G4double energytemp = 3.0996 * eV; // 400 nm
      G4double rndm = G4UniformRand();
      G4double rndm2 = G4UniformRand();
      G4double costheta = 1- rndm*(1 - std::cos(ledtheta*deg));
      G4double sintheta = std::sqrt(1. - costheta*costheta);
      G4double phi = rndm2*2*(std::acos(-1));
      G4double sinphi = std::sin(phi);
      G4double cosphi = std::cos(phi);
      G4double px = sintheta*cosphi;
      G4double py = sintheta*sinphi;
      G4double pz = costheta;
      G4ThreeVector XPrime = G4ThreeVector(0.,0.,1.);
      G4ThreeVector YPrime = G4ThreeVector(-ypos,xpos,0.);
      G4ThreeVector ZPrime = G4ThreeVector(-xpos,-ypos,0.);
      G4double finx = (px*XPrime.x()) + (py*YPrime.x()) + (pz*ZPrime.x());
      G4double finy = (px*XPrime.y()) + (py*YPrime.y()) + (pz*ZPrime.y());
      G4double finz = (px*XPrime.z()) + (py*YPrime.z()) + (pz*ZPrime.z());
      G4ThreeVector momtemp(finx,finy,finz);
      particleGunTop->SetParticleMomentumDirection(momtemp);
      particleGunTop->SetParticleEnergy(energytemp);
      particleGunTop->SetParticlePolarization(G4RandomDirection());
      particleGunTop->GeneratePrimaryVertex(anEvent);

      SetVtx(postemp);
      SetBeamEnergy(energytemp);
      SetBeamDir(momtemp.unit());
      SetBeamPDG(0);

    }
    else if(useTimeTop) {

      std::cerr << "using time top" << std::endl;
      particleGunTop->SetParticleDefinition(G4OpticalPhoton::Definition());
      particleGunTop->SetNumberOfParticles(c_particle);

      std::vector<double> gaussian_times;
      G4double rndm = G4UniformRand();
      G4double rndm2 = G4UniformRand();
      for (int time_vec = 0; time_vec < 20000; time_vec++) {
        G4double time = G4RandGauss::shoot(200.0,40.0);
        gaussian_times.push_back(time);
	if (time_vec % 1000 == 0) std::cerr << "processing times" << std::endl;
      }
      std::sort(gaussian_times.begin(),gaussian_times.end());
      for (int part_num = 0; part_num < gaussian_times.size(); part_num++) {
	if (part_num % 1000 == 0) std::cerr << "processing particles" << std::endl;
        G4double time = gaussian_times.at(part_num);
	double xpos = -1112.08*cm;
	double ypos = 327.084*cm;
	double zpos = 3600*cm;
	G4ThreeVector postemp(xpos,ypos,zpos);
	particleGunTop->SetParticlePosition(postemp);
	G4double energytemp = 3.0996 * eV; // 400 nm

	G4double costheta = 1- rndm*(1 - std::cos(ledtheta));
	G4double sintheta = std::sqrt(1. - costheta*costheta);
	G4double phi = rndm2*2*(std::acos(-1));
	G4double sinphi = std::sin(phi);
	G4double cosphi = std::cos(phi);
	G4double px = sintheta*cosphi;
	G4double py = sintheta*sinphi;
	G4double pz = -costheta;
	G4ThreeVector momtemp(px,py,pz);
	//particleGunTop->SetParticleTime(time);
	particleGunTop->SetParticleMomentumDirection(momtemp);
	particleGunTop->SetParticleEnergy(energytemp);
	particleGunTop->SetParticlePolarization(G4RandomDirection());
	particleGunTop->GeneratePrimaryVertex(anEvent);

	//	std::cerr << "Position x: " << particleGunTop->GetParticlePosition().x() << " y: " << particleGunTop->GetParticlePosition().y() << " z: " << particleGunTop->GetParticlePosition().z() << std::endl;
	//std::cerr << "Direction x: " << particleGunTop->GetParticleMomentumDirection().x() << " y: " << particleGunTop->GetParticleMomentumDirection().y() << " z: " << particleGunTop->GetParticleMomentumDirection().z() << std::endl;
	SetVtx(postemp);
	SetBeamEnergy(energytemp);
	SetBeamDir(momtemp.unit());
	SetBeamPDG(0);
      }
    } else if(useIWCDFullInjectors) {
      //      std::cerr << "Running IWCD injectors with " << c_particle << " photons and a half angle of " << ledtheta << " degrees" << std::endl;
      //------------------------------- Light array ------------------------     
      particleGunBarrel->SetParticleDefinition(G4OpticalPhoton::Definition());
      particleGunBarrel->SetNumberOfParticles(c_particle);
      particleGunTop->SetParticleDefinition(G4OpticalPhoton::Definition());
      particleGunTop->SetNumberOfParticles(c_particle);
      particleGunBottom->SetParticleDefinition(G4OpticalPhoton::Definition());
      particleGunBottom->SetNumberOfParticles(c_particle);
      //particleTable = FindParticle("opticalphoton")->DumpTable();

      double maxr = 335*cm; // IWCD GEOM
      //double maxr = 3400*cm;
      //double rad1 = 200*cm;
      double rad1 = 150*cm;
      //      double rad2 = 1800*cm;
      //double rad3 = 2800*cm;
      double xposb = -999*cm;
      double yposb = -999*cm;
      //double ledtheta = 40 * deg;

      //      for (int numpart = 0; numpart < 500; numpart++) {
      //for (int timestep = 0; timestep < 1000; timestep++) {
      for(int ip = 0; ip < 3; ip++) {
	//std::cerr << "looping through Z positions" << std::endl;
	//double zposb = (-1815 + ip*1210)*cm;
	//double zposb = (-1800 + ip*1200)*cm;
	double zposb = -9999999;
	if (ip==0) zposb = 250*cm;
	if (ip==1) zposb = 0*cm;
	if (ip==2) zposb = -250*cm;
	//double zposb = 1085*cm;
	//if (ip == 0 || ip == 2) continue;
	//for(int theta = 0; theta < 6; theta++) {
	for(int theta = 0; theta < 4; theta++) {
	  //  	  if (theta % 2 != 0) continue;
	  //MyGPS->GeneratePrimaryVertex( anEvent );
	  // CPFLAG These theta offsets are there because there seems to be some
	  // overlap in the geometry with the PMTs and injection positions without
	  // the +0.2
	  /*	  if (ip == 1) {
	    xposb = maxr*(std::cos((theta+0.2)*(std::acos(-1))/3));
	    yposb = maxr*(std::sin((theta+0.2)*(std::acos(-1))/3));
	    //	      if ((ip == 1 && theta == 0) || (ip == 3 && theta == 7)) continue;
	  }
	  else {
	    xposb = maxr*(std::cos((theta+0.7)*(std::acos(-1))/3));
	    yposb = maxr*(std::sin((theta+0.7)*(std::acos(-1))/3));
	    //if ((ip == 0 && theta == 13) || (ip == 2 && theta == 3) || (ip == 4 && theta == 11)) continue;
	    }*/
	  if (ip == 1) {
	    xposb = maxr*(std::cos((theta+0.2)*(std::acos(-1))/2));
	    yposb = maxr*(std::sin((theta+0.2)*(std::acos(-1))/2));
	    //	      if ((ip == 1 && theta == 0) || (ip == 3 && theta == 7)) continue;
	  }
	  else {
	    xposb = maxr*(std::cos((theta+0.7)*(std::acos(-1))/2));
	    yposb = maxr*(std::sin((theta+0.7)*(std::acos(-1))/2));
	    //if ((ip == 0 && theta == 13) || (ip == 2 && theta == 3) || (ip == 4 && theta == 11)) continue;
	  }

	  //	  std::cerr << "Position : " << xposb << "," << yposb << "," << zposb << std::endl;
	  //G4ThreeVector postemp(xpos, ypos, zpos);
	  G4ThreeVector postempb(xposb,yposb,zposb);
	  particleGunBarrel->SetParticlePosition(postempb);
	  //G4double energytempb = 2.505 * eV;
	  G4double energytempb = 3.0996 * eV; // 400 nm
	  //G4ThreeVector momtempb(-xposb, -yposb, 0);
	  G4double rndm = G4UniformRand();
	  G4double rndm2 = G4UniformRand();
	  G4double costheta = 1- rndm*(1 - std::cos(ledtheta*deg));
	  G4double sintheta = std::sqrt(1. - costheta*costheta);
	  G4double phi = rndm2*2*(std::acos(-1));
	  G4double sinphi = std::sin(phi);
	  G4double cosphi = std::cos(phi);
	  G4double px = sintheta*cosphi;
	  G4double py = sintheta*sinphi;
	  G4double pz = costheta;
	  G4ThreeVector XPrime = G4ThreeVector(0.,0.,1.);
	  //G4ThreeVector YPrime = G4ThreeVector(-yposb,xposb,0.); // Facing in
	  G4ThreeVector YPrime = G4ThreeVector(yposb,-xposb,0.); // Facing out
	  //G4ThreeVector ZPrime = G4ThreeVector(-xposb,-yposb,0.); // Facing in
	  G4ThreeVector ZPrime = G4ThreeVector(xposb,yposb,0.); // Facing out
	  //G4double pz = 0;
	  G4double finx = (px*XPrime.x()) + (py*YPrime.x()) + (pz*ZPrime.x());
	  G4double finy = (px*XPrime.y()) + (py*YPrime.y()) + (pz*ZPrime.y());
	  G4double finz = (px*XPrime.z()) + (py*YPrime.z()) + (pz*ZPrime.z());
	  //std::cerr << "finx " << finx << ", finy " << finy << ", finz " << finz << std::endl;
	  //G4ThreeVector momtempb(px,py,pz);
	  G4ThreeVector momtempb(finx,finy,finz);
	  particleGunBarrel->SetParticleMomentumDirection(momtempb);
	  particleGunBarrel->SetParticleEnergy(energytempb);
	  //particleGunBarrel->SetParticleTime(timestep);
	  //G4ThreeVector targetvector(-xposb,-yposb,0);
	  G4ThreeVector targetvector(xposb,yposb,0);
	  particleGunBarrel->SetParticlePolarization(G4RandomDirection());
	  particleGunBarrel->GeneratePrimaryVertex(anEvent);
	  SetVtx(postempb);
	  SetBeamEnergy(energytempb);
	  SetBeamDir(momtempb.unit());
	  SetBeamPDG(0);
	}
      }
      //}
      for (int i = 0; i < 3; i++) {
	//double zpos = 2990*cm; // top facing down
	// double zpos = -3555*cm; // bottom facing down
	double zpos = -435*cm; // bottom facing down
	double xpos = 0.;
	double ypos = 0.;
	//G4ThreeVector postemp(0, 0, zpos);
	if (i==0) {
	  //	  for (int angfrac1 = 0; angfrac1 < 4; angfrac1++) {
	  for (int angfrac1 = 0; angfrac1 < 2; angfrac1++) {
	    //if (angfrac1 % 2 == 0) continue;
	    //	    double xpos = rad1*(std::cos((angfrac1)*(std::acos(-1)/2)));
	    //double ypos = rad1*(std::sin((angfrac1)*(std::acos(-1)/2)));
	    double xpos = rad1*(std::cos((angfrac1)*(std::acos(-1))));
	    double ypos = rad1*(std::sin((angfrac1)*(std::acos(-1))));
	    G4ThreeVector postemp(xpos,ypos,zpos);
	    particleGunTop->SetParticlePosition(postemp);
	    //particleGun->SetNumberOfParticles(5000);
	    //	    G4double energytemp = 2.505 * eV;
	    G4double energytemp = 3.0996 * eV; // 400 nm
	    G4double rndm = G4UniformRand();
	    G4double rndm2 = G4UniformRand();
	    G4double costheta = 1- rndm*(1 - std::cos(ledtheta*deg));
	    G4double sintheta = std::sqrt(1. - costheta*costheta);
	    G4double phi = rndm2*2*(std::acos(-1));
	    G4double sinphi = std::sin(phi);
	    G4double cosphi = std::cos(phi);
	    G4double px = sintheta*cosphi;
	    G4double py = sintheta*sinphi;
	    G4double pz = -costheta; // Facing down
	    //G4double pz = -1;
	    G4ThreeVector momtemp(px,py,pz);
	    particleGunTop->SetParticleMomentumDirection(momtemp);
	    particleGunTop->SetParticleEnergy(energytemp);
	    particleGunTop->SetParticlePolarization(G4RandomDirection());
	    particleGunTop->GeneratePrimaryVertex(anEvent);

	    SetVtx(postemp);
	    SetBeamEnergy(energytemp);
	    SetBeamDir(momtemp.unit());
	    SetBeamPDG(0);
	  }
	}
	if (i==1) {
	  continue;
	}
	if (i==2) {
	  continue;
	}
      }
      for (int i = 0; i < 3; i++) {
	//double zpos = -2990*cm; // bottom facing up
	//double zpos = 3555*cm; // top facing up
	double zpos = 435*cm; // top facing up
	double xpos = 0.;
	double ypos = 0.;
	//G4ThreeVector postemp(0, 0, zpos);
	if (i==0) {
	  //for (int angfrac1 = 0; angfrac1 < 4; angfrac1++) {
	  for (int angfrac1 = 0; angfrac1 < 2; angfrac1++) {
	    //if (angfrac1 % 2 == 0) continue;
	    //	    double xpos = rad1*(std::cos((angfrac1)*(std::acos(-1)/2)));
	    //double ypos = rad1*(std::sin((angfrac1)*(std::acos(-1)/2)));
	    double xpos = rad1*(std::cos((angfrac1)*(std::acos(-1))));
	    double ypos = rad1*(std::sin((angfrac1)*(std::acos(-1))));
	    G4ThreeVector postemp(xpos,ypos,zpos);
	    particleGunBottom->SetParticlePosition(postemp);
	    //particleGun->SetNumberOfParticles(5000);
	    //G4double energytemp = 2.505 * eV;
	    G4double energytemp = 3.0996 * eV; // 400 nm
	    G4double rndm = G4UniformRand();
	    G4double rndm2 = G4UniformRand();
	    G4double costheta = 1- rndm*(1 - std::cos(ledtheta*deg));
	    G4double sintheta = std::sqrt(1. - costheta*costheta);
	    G4double phi = rndm2*2*(std::acos(-1));
	    G4double sinphi = std::sin(phi);
	    G4double cosphi = std::cos(phi);
	    G4double px = sintheta*cosphi;
	    G4double py = sintheta*sinphi;
	    G4double pz = costheta; // Facing up
	    G4ThreeVector momtemp(px,py,pz);
	    particleGunBottom->SetParticleMomentumDirection(momtemp);
	    particleGunBottom->SetParticleEnergy(energytemp);
	    particleGunBottom->SetParticlePolarization(G4RandomDirection());
	    particleGunBottom->GeneratePrimaryVertex(anEvent);

	    SetVtx(postemp);
	    SetBeamEnergy(energytemp);
	    SetBeamDir(momtemp.unit());
	    SetBeamPDG(0);
	  }
	}
	if (i==1) {
	  continue;
	}
	if (i==2) {
	  continue;
	}
      }
    } // ------------------------------ end of IWCD light array ------------------------------ 
  }
}

void WCSimPrimaryGeneratorAction::SaveOptionsToOutput(WCSimRootOptions * wcopt)
{
  if(useMulineEvt)
    wcopt->SetVectorFileName(vectorFileName);
  else
    wcopt->SetVectorFileName("");
  wcopt->SetGeneratorType(GetGeneratorTypeString());
}

G4String WCSimPrimaryGeneratorAction::GetGeneratorTypeString()
{
  if(useMulineEvt)
    return "muline";
  else if(useGunEvt)
    return "gun";
  else if(useGPSEvt)
    return "gps";
  else if(useLaserEvt)
    return "laser";
  else if(useCosmics)
    return "cosmics";
  else if(useCalibration)
    return "calibration";
  return "";
}

// Returns a vector with the tokens
vector<string> tokenize( string separators, string input ) 
{
  std::size_t startToken = 0, endToken; // Pointers to the token pos
  vector<string> tokens;  // Vector to keep the tokens
  
  if( separators.size() > 0 && input.size() > 0 ) 
    {
    
      while( startToken < input.size() )
	{
	  // Find the start of token
	  startToken = input.find_first_not_of( separators, startToken );
      
	  // If found...
	  if( startToken != input.npos ) 
	    {
	      // Find end of token
	      endToken = input.find_first_of( separators, startToken );
	      if( endToken == input.npos )
		// If there was no end of token, assign it to the end of string
		endToken = input.size();
        
	      // Extract token
	      tokens.push_back( input.substr( startToken, endToken - startToken ) );
        
	      // Update startToken
	      startToken = endToken;
	    }
	}
    }
  
  return tokens;
}

