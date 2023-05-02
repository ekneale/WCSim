#include "WCSimPrimaryGeneratorAction.hh"
#include "WCSimDetectorConstruction.hh"
#include "WCSimPrimaryGeneratorMessenger.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4GeneralParticleSource.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
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
  // Any line starting with # is ignored
	while(true)                                               
	{  
		if (inFile.getline(inBuf,lineSize))
		{
			if(inBuf[0]!='#')                                  
				return tokenize(" $\r", inBuf);
		}
		else
		{
		  if(inFile.fail())
		    G4cerr << "Failed to read line. Is the buffer size large enough?" << G4endl;
			vector<string> nullLine;                               
			return nullLine; 
		}
	}
}

inline double atof( const string& S ) {return std::atof( S.c_str() );}
inline int    atoi( const string& S ) {return std::atoi( S.c_str() );}

WCSimPrimaryGeneratorAction::WCSimPrimaryGeneratorAction(
					  WCSimDetectorConstruction* myDC)
  :myDetector(myDC), vectorFileName("")
{
  //T. Akiri: Initialize GPS to allow for the laser use 
  MyGPS = new G4GeneralParticleSource();

  // Initialize to zero
  mode[0] = 0;
  nvtxs = 0;
  for( Int_t u=0; u<MAX_N_VERTICES; u++){
    vtxsvol[u] = 0;
    vtxs[u] = G4ThreeVector(0.,0.,0.);
  }

  nuEnergy = 0.;
  _counterRock=0; // counter for generated in Rock
  _counterCublic=0; // counter generated
  
  //---Set defaults. Do once at beginning of session.
  
  G4int n_particle = 1;
  particleGun = new G4ParticleGun(n_particle);
  particleGun->SetParticleEnergy(1.0*GeV);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.0));

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
//  G4IonTable* ionTable = G4IonTable::GetIonTable();
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
  useSatBottom = false;
  useSatBarrel = false;
  useIWCDFullInjectors = false;
  useCosmics            = false;
  useRadioactiveEvt  	= false;
  useRadonEvt        	= false;

  particleGunBarrel = new G4ParticleGun();
  particleGunTop = new G4ParticleGun();
  particleGunBottom = new G4ParticleGun();

  // OD Detector radius and height
  OD_inner_radius = myDC->GetIDRadius() + 60*cm + 2*cm + 1*mm; // does not include plate or PMT exposure
  OD_inner_height = myDC->GetWCIDHeight() + 2*(60*cm + 2*cm + 1*mm); // does not include plate or PMT exposure

// Create the relevant histograms to generate muons
  // according to SuperK flux extrapolated at HyperK site
  altCosmics = 2*myDC->GetWCIDHeight();
  G4cout << "altCosmics : " << altCosmics << G4endl;
  if (inputCosmicsFile.is_open())
    inputCosmicsFile.close();


  inputCosmicsFile.open(cosmicsFileName, std::fstream::in);

  if (!inputCosmicsFile.is_open()) {
    G4cout << "Cosmics data file " << cosmicsFileName << " not found" << G4endl;
	exit(-1);
  } else {
    G4cout << "Cosmics data file " << cosmicsFileName << " found" << G4endl;
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

    while ( getline(inputCosmicsFile,line) ){
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

    TFile *file = new TFile("cosmicflux.root","RECREATE");
    hFluxCosmics->Write();
    hEmeanCosmics->Write();
    file->Close();

  }
  
  // Radioactive and Radon generator variables:
  radioactive_sources.clear();
  myRn222Generator	= 0;
  fRnScenario		= 0;
  fRnSymmetry		= 1;
  // Time units for vertices   
  fTimeUnit=CLHEP::nanosecond;
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
  inputCosmicsFile.close();
  delete particleGun;
  delete MyGPS;   //T. Akiri: Delete the GPS variable
  delete messenger;
}

void WCSimPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  // We will need a particle table
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4IonTable* ionTable = G4IonTable::GetIonTable();
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
	// The original documentation describing the nuance text format can be found here: 
    // http://neutrino.phy.duke.edu/nuance-format/
    //
    // Information specific to WCSim can be found in the file Nuance_MC_Format.txt in
    // the doc directory.
    // The format must be strictly adhered to for it to be processed correctly.
    // The lines and their meanings from begin through info are fixed, and then
    // a variable number of tracks may follow.
    //
    if (useNuanceTextFormat)
    {
    	const int lineSize=1000;
    	char      inBuf[lineSize];
    	vector<string> token(1);

    	token = readInLine(inputFile, lineSize, inBuf);

    	if (token.size() == 0 || token[0] == "stop" ) 
    	{
    		G4cout << "End of nuance vector file - run terminated..."<< G4endl;
    		G4RunManager::GetRunManager()-> AbortRun();
    	}
    	else if (token[0] != "begin" )
    	{
    		G4cout << "unexpected line begins with " << token[0] << " we were expecting \" begin \" "<<G4endl;
    	}
	else   // normal parsing begins here
	{
		// Read the nuance line 
        // should be nuance <value>
        // but could be just  
        // nuance 
		// if value is given set mode to equal it.

			token = readInLine(inputFile, lineSize, inBuf);
			int iVertex=0;
			while(token[0]=="nuance" && iVertex < MAX_N_VERTICES)
			{
				if(token.size()>1)
					mode[iVertex] = atoi(token[1]);
	            // Read the Vertex line
				token = readInLine(inputFile, lineSize, inBuf);
				vtxs[iVertex] = G4ThreeVector(atof(token[1])*cm,
					atof(token[2])*cm,
					atof(token[3])*cm);
				G4double VertexTime=atof(token[4])*fTimeUnit; 
				vertexTimes[iVertex]=VertexTime;
                // true : Generate vertex in Rock , false : Generate vertex in WC tank
				SetGenerateVertexInRock(false);
	            // Next we read the incoming neutrino and target
	            // First, the neutrino line
				token=readInLine(inputFile, lineSize, inBuf);
				beampdgs[iVertex] = atoi(token[1]);
				beamenergies[iVertex] = atof(token[2])*MeV;
				beamdirs[iVertex] = G4ThreeVector(atof(token[3]),
					atof(token[4]),
					atof(token[5]));
				SetBeamEnergy(beamenergies[iVertex]);         
                SetBeamDir(beamdirs[iVertex]);
	            // Now read the target line
				token=readInLine(inputFile, lineSize, inBuf);
				targetpdgs[iVertex] = atoi(token[1]);
				targetenergies[iVertex] = atof(token[2])*MeV;
				targetdirs[iVertex] = G4ThreeVector(atof(token[3]),
					atof(token[4]),
					atof(token[5]));
	            // Read the info line, basically a dummy
				token=readInLine(inputFile, lineSize, inBuf);
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
						G4double tempEnergy = atof(token[2])*MeV;
						G4ThreeVector dir = G4ThreeVector(atof(token[3]),
							atof(token[4]),
							atof(token[5]));
		                //must handle the case of an ion seperatly from other particles
		                //check PDG code if we have an ion.
		                //PDG code format for ions ±10LZZZAAAI
						char strPDG[11];
						char strA[10]={0};
						char strZ[10]={0};
						long int A=0,Z=0;
						if(abs(pdgid) >= 1000000000)
						{
			                //ion
							sprintf(strPDG,"%i",abs(pdgid));
							strncpy(strZ, &strPDG[3], 3);
							strncpy(strA, &strPDG[6], 3);
							strA[3]='\0';
							strZ[3]='\0';
							A=atoi(strA);
							Z=atoi(strZ);
							G4ParticleDefinition* ion;
							ion =  ionTable->GetIon(Z, A, 0.);
							particleGun->SetParticleDefinition(ion);
							particleGun->SetParticleCharge(0);
						}
						else {
		                    //not ion
							particleGun->
							SetParticleDefinition(particleTable->
								FindParticle(pdgid));
						}
						G4double mass = 
						particleGun->GetParticleDefinition()->GetPDGMass();

						G4double ekin = tempEnergy - mass;

						particleGun->SetParticleEnergy(ekin);
						particleGun->SetParticlePosition(vtxs[iVertex]);
						particleGun->SetParticleMomentumDirection(dir);
						particleGun->SetParticleTime(VertexTime);
						particleGun->GeneratePrimaryVertex(anEvent);

					}
				}
				iVertex++;
				if(iVertex > MAX_N_VERTICES)
					G4cout<<" CAN NOT DEAL WITH MORE THAN "<<MAX_N_VERTICES<<" VERTICES - TRUNCATING EVENT HERE "<<G4endl;
			}
			nvtxs=iVertex;
			SetNvtxs(nvtxs);

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

	    //To prevent occasional seg fault from an un assigned targetpdg 
    targetpdgs[0] = 2212; //ie. proton

    G4ThreeVector P  =anEvent->GetPrimaryVertex()->GetPrimary()->GetMomentum();
    G4ThreeVector vtx=anEvent->GetPrimaryVertex()->GetPosition();
    G4double mass    =anEvent->GetPrimaryVertex()->GetPrimary()->GetMass();
    G4int pdg        =anEvent->GetPrimaryVertex()->GetPrimary()->GetPDGcode();

    char strPDG[11];
    char strA[10]={0};
    char strZ[10]={0};
    
    
    long int A=0,Z=0;
    //		    A=strotl(strPDG,&str);
    if(abs(pdg) >= 1000000000)
      {
	//ion
	sprintf(strPDG,"%i",abs(pdg));
	strncpy(strZ, &strPDG[3], 3);
	strncpy(strA, &strPDG[6], 3);
	strA[3]='\0';
	strZ[3]='\0';
	A=atoi(strA);
	Z=atoi(strZ);

	G4ParticleDefinition* ion   = G4IonTable::GetIonTable()->GetIon(Z, A, 0);
	ion->SetPDGStable(false);
	ion->SetPDGLifeTime(0.);
	
	G4ParticleDefinition* ion2   = G4IonTable::GetIonTable()->GetIon(Z, A, 0);
	std::cout<<"ion2 "<<ion2->GetPDGLifeTime()<<"\n";
      }
    
    G4ThreeVector dir  = P.unit();
    G4double E         = std::sqrt((P.dot(P))+(mass*mass));

    SetVtx(vtx);
    SetBeamEnergy(E);
    SetBeamDir(dir);
    SetBeamPDG(pdg);
  }
  else if (useLaserEvt)
    {
      //T. Akiri: Create the GPS LASER event
      MyGPS->GeneratePrimaryVertex(anEvent);

	  //To prevent occasional seg fault from an un assigned targetpdg 
      targetpdgs[0] = 2212; //ie. proton
      
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
      G4double mass     =anEvent->GetPrimaryVertex()->GetPrimary()->GetMass();
      G4int pdg         =anEvent->GetPrimaryVertex()->GetPrimary()->GetPDGcode();
      
      G4ThreeVector dir  = P.unit();
      G4double E         = std::sqrt((P.dot(P))+(mass*mass));
      G4cout << " GPS primary vertex (" << vtx.x() << ", " << vtx.y() << ", " << vtx.z() << "), dir ("
	     << dir.x() << ", " << dir.y() << ", " << dir.z() << ") m " << m << " E "<< E << " pdg " << pdg << G4endl;
      
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

    double planck = 4.1357e-15; // ev.s
    double lightspeed = 299792458e9; // nm/s
    G4double energytemp = (planck*lightspeed/calib_wavelength) * eV;

    if(useFullInjectors) {

      //std::cerr << "Running all injectors with " << c_particle << " photons and a half angle of " << ledtheta << " degrees" << std::endl;
      //------------------------------- Light array ------------------------     
      particleGunBarrel->SetParticleDefinition(G4OpticalPhoton::Definition());
      particleGunBarrel->SetNumberOfParticles(c_particle);
      particleGunTop->SetParticleDefinition(G4OpticalPhoton::Definition());
      particleGunTop->SetNumberOfParticles(c_particle);
      particleGunBottom->SetParticleDefinition(G4OpticalPhoton::Definition());
      particleGunBottom->SetNumberOfParticles(c_particle);

      double maxr = 3317.01*cm;
      //double maxr = OD_inner_radius + 10*cm; // facing out/ For some reason getting the radius didn't work
      //double maxr = OD_inner_radius + 95*cm; // facing in/  hence the hard code above
      double rad1 = 800*cm;
      double rad2 = 1810*cm;// try adding 10 cm to miss struts
      double rad3 = 2800*cm;
      double xposb = -9999*cm;
      double yposb = -9999*cm;
      double zposb = -9999*cm;


      for(int ip = 0; ip < 5; ip++) {
	if (ip==0) zposb = 2700*cm;
	if (ip==1) zposb = 1350*cm;
	//	if (ip==2) zposb = 0*cm; // this position intersected PMTs
	if (ip==2) zposb = 60*cm;
	if (ip==3) zposb = -1350*cm;
	if (ip==4) zposb = -2700*cm;
	for(int theta = 0; theta < 16; theta++) {
	  if (ip == 1 || ip == 3) {
	    xposb = maxr*(std::cos((theta+0.0)*(std::acos(-1))/8));
	    yposb = maxr*(std::sin((theta+0.0)*(std::acos(-1))/8));
	    //std::cerr << "POS theta " << theta << " : " << xposb << ", " << yposb << std::endl;
	  }
	  else {
	    xposb = maxr*(std::cos((theta+0.5)*(std::acos(-1))/8));
	    yposb = maxr*(std::sin((theta+0.5)*(std::acos(-1))/8));
	  }
      /* // This is a 7 ring barrel option, if people decide they want
	 // as many vertical settings in the OD as the ID
	for(int ip = 0; ip < 7; ip++) {
	if (ip==0) zposb = 2900*cm;
	if (ip==1) zposb = 1950*cm;
	if (ip==2) zposb = 1000*cm;
	if (ip==3) zposb = 60*cm;
	if (ip==4) zposb = -900*cm;
	if (ip==5) zposb = -1900*cm;
	if (ip==6) zposb = -2900*cm;
	int theta_lim = 12;
	if (ip==2 || ip==4) theta_lim = 10;
	for(int theta = 0; theta < theta_lim; theta++) {
	  /*if (ip == 1 || ip == 3) {
	    xposb = maxr*(std::cos((theta+0.0)*(std::acos(-1))/(theta_lim/2)));
	    yposb = maxr*(std::sin((theta+0.0)*(std::acos(-1))/(theta_lim/2)));
	    //std::cerr << "POS theta " << theta << " : " << xposb << ", " << yposb << std::endl;
	  }
	  else {
	    xposb = maxr*(std::cos((theta+0.5)*(std::acos(-1))/(theta_lim/2)));
	    yposb = maxr*(std::sin((theta+0.5)*(std::acos(-1))/(theta_lim/2)));
	    //	  }*/
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
	double zpos = -3362.01*cm;
	//double zpos = -(OD_inner_height/2 + 10*cm); // bottom facing down
	//double zpos = (OD_inner_height/2 + 190*cm); // top facing down
	double xpos = 0.;
	double ypos = 0.;
	//G4ThreeVector postemp(0, 0, zpos);
	if (i==0) {
	  for (int angfrac1 = 0; angfrac1 < 3; angfrac1++) {
	    double xpos = rad1*(std::cos((angfrac1*8 + 3)*(std::acos(-1)/12)));
	    double ypos = rad1*(std::sin((angfrac1*8 + 3)*(std::acos(-1)/12)));
	    G4ThreeVector postemp(xpos,ypos,zpos);
	    particleGunBottom->SetParticlePosition(postemp);
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
	    double xpos = rad2*(std::cos((angfrac2*4 + 1)*(std::acos(-1)/12)));
	    double ypos = rad2*(std::sin((angfrac2*4 + 1)*(std::acos(-1)/12)));
	    G4ThreeVector postemp(xpos,ypos,zpos);
	    particleGunBottom->SetParticlePosition(postemp);
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
	    double xpos = rad3*(std::cos((angfrac3)*(std::acos(-1)/6)));
	    double ypos = rad3*(std::sin((angfrac3)*(std::acos(-1)/6)));
	    G4ThreeVector postemp(xpos,ypos,zpos);
	    particleGunBottom->SetParticlePosition(postemp);
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
	double zpos = 3362.01*cm;
	//double zpos = (OD_inner_height/2 + 10*cm); // top facing up
	double xpos = 0.;
	double ypos = 0.;
	if (i==0) {
	  for (int angfrac1 = 0; angfrac1 < 3; angfrac1++) {
	    double xpos = rad1*(std::cos((angfrac1*8 + 3)*(std::acos(-1)/12)));
	    double ypos = rad1*(std::sin((angfrac1*8 + 3)*(std::acos(-1)/12)));
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
	    double xpos = rad2*(std::cos((4*angfrac2 + 1)*(std::acos(-1)/12)));
	    double ypos = rad2*(std::sin((4*angfrac2 + 1)*(std::acos(-1)/12)));
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
	    G4double pz = costheta;
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
	    double xpos = rad3*(std::cos((angfrac3)*(std::acos(-1)/6)));
	    double ypos = rad3*(std::sin((angfrac3)*(std::acos(-1)/6)));
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
	    G4double pz = costheta;
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
    } // ------------------------------ end of diffuse light array ------------------------------

    // ------------- One single inner ring source for top cap saturation studies -----------------

    if(useSatTop) {

      particleGunTop->SetParticleDefinition(G4OpticalPhoton::Definition());
      particleGunTop->SetNumberOfParticles(c_particle);

      double rad1 = 800*cm;
      double zpos = 3362.01*cm; 
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

    // ------------- One single inner ring source for bottom cap saturation studies -----------------

    if(useSatBottom) {

      particleGunTop->SetParticleDefinition(G4OpticalPhoton::Definition());
      particleGunTop->SetNumberOfParticles(c_particle);

      double rad1 = 800*cm;
      double zpos = -3362.01*cm; 
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
      G4double pz = -costheta; // Facing down
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
    

    // ------------- One single mid-height source for barrel saturation studies -----------------    
    if(useSatBarrel) {
      
      particleGunBarrel->SetParticleDefinition(G4OpticalPhoton::Definition());
      particleGunBarrel->SetNumberOfParticles(c_particle);

      double temp_OD_inner_radius = 32.4;
      double maxr = 3317.01*cm; // facing out
      double xposb = -9999*cm;
      double yposb = -9999*cm;
      double zposb = 60*cm;
      double theta = 0.05;
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

    //--------------------------------- Barrel collimator -------------------------
    else if(useColBarrel) {
      std::cerr << "Running barrel with " << c_particle << " photons with a wavelength of " << calib_wavelength << " and a half angle of " << ledtheta << " degrees" << std::endl;
      // ------------------------- collimated barrel --------------------
      particleGunBarrel->SetParticleDefinition(G4OpticalPhoton::Definition());
      particleGunBarrel->SetNumberOfParticles(c_particle);
      double xposb = 3314.91*cm;//OD_inner_radius + 10*cm + 2.9*cm;
      double yposb = 0*cm;
      double zposb = -3352.01*cm;//-OD_inner_height/2;
      G4ThreeVector postemp(xposb,yposb,zposb);
      particleGunBarrel->SetParticlePosition(postemp);
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
      G4ThreeVector momtemp(px,py,pz);
      particleGunBarrel->SetParticleMomentumDirection(momtemp);
      particleGunBarrel->SetParticleEnergy(energytemp);
      particleGunBarrel->SetParticlePolarization(G4RandomDirection());
      particleGunBarrel->GeneratePrimaryVertex(anEvent);
      SetVtx(postemp);
      SetBeamEnergy(energytemp);
      SetBeamDir(momtemp.unit());
      SetBeamPDG(0);
      
      // ------------------------ end of collimated barrel -----------------

    }

    else if(useColBottom) {
      // ---------------- single collimated bottom cap pointing across -----
      std::cerr << "Running bottom with " << c_particle << " photons with a wavelength of " << calib_wavelength << " and a half angle of " << ledtheta << " degrees" << std::endl;

      particleGunBottom->SetParticleDefinition(G4OpticalPhoton::Definition());
      particleGunBottom->SetNumberOfParticles(c_particle);
      double xpos = 3302.01*cm;//OD_inner_radius;
      double ypos = 0*cm;
      double zpos = -3364.91*cm;//-(OD_inner_height/2 + 10*cm + 2.9*cm); // zpos of start of OD + 2.9cm PMT expose height + 10cm to clear PMTs
      G4ThreeVector postemp(xpos,ypos,zpos);
      particleGunBottom->SetParticlePosition(postemp);
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

    //---------------------- Top Cap Collimator
    else if(useColTop) {

      std::cerr << "Running top with " << c_particle << " photons with a wavelength of " << calib_wavelength << " and a half angle of " << ledtheta << " degrees" << std::endl;
      particleGunTop->SetParticleDefinition(G4OpticalPhoton::Definition());
      particleGunTop->SetNumberOfParticles(c_particle);
      double xpos = 3308.03*cm;//OD_inner_radius;
      double ypos = 0*cm;
      xpos = 3200*cm; // move in to avoid outer frame of support structure
      ypos = 0*cm;
      double zpos = 3364.91*cm;//(OD_inner_height/2 + 10*cm + 2.9*cm); // zpos of start of OD + 2.9cm PMT expose height + 10cm to clear PMTs
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
      //xpos 73.2285, ypos 73.2285, zpos -3348.93
      std::cerr << "using time top but now bottom and less about time" << std::endl;
      particleGunTop->SetParticleDefinition(G4OpticalPhoton::Definition());
      particleGunTop->SetNumberOfParticles(c_particle);

      G4double rndm = G4UniformRand();
      G4double rndm2 = G4UniformRand();
      double xpos = 73.2285*cm; // centre of PMT
      double ypos = 73.2285*cm; // centre of PMT
      double zpos = -3363.93*cm; // start 15 cm from centre(?) of PMT
      G4ThreeVector postemp(xpos,ypos,zpos);
      particleGunTop->SetParticlePosition(postemp);
      G4double energytemp = 3.0996 * eV; // 400 nm

      G4double costheta = 1- rndm*(1 - std::cos(0.1*deg));
      G4double sintheta = std::sqrt(1. - costheta*costheta);
      G4double phi = rndm2*2*(std::acos(-1));
      G4double sinphi = std::sin(phi);
      G4double cosphi = std::cos(phi);
      G4double px = sintheta*cosphi;
      G4double py = sintheta*sinphi;
      G4double pz = costheta; // face up at bottom PMT
      G4ThreeVector momtemp(px,py,pz);
      //particleGunTop->SetParticleTime(time);
      particleGunTop->SetParticleMomentumDirection(momtemp);
      particleGunTop->SetParticleEnergy(energytemp);
      particleGunTop->SetParticlePolarization(G4RandomDirection());
      particleGunTop->GeneratePrimaryVertex(anEvent);

      SetVtx(postemp);
      SetBeamEnergy(energytemp);
      SetBeamDir(momtemp.unit());
      SetBeamPDG(0);


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
	  particleGunBarrel->SetParticleEnergy(energytemp);
	  //particleGunBarrel->SetParticleTime(timestep);
	  //G4ThreeVector targetvector(-xposb,-yposb,0);
	  G4ThreeVector targetvector(xposb,yposb,0);
	  particleGunBarrel->SetParticlePolarization(G4RandomDirection());
	  particleGunBarrel->GeneratePrimaryVertex(anEvent);
	  SetVtx(postempb);
	  SetBeamEnergy(energytemp);
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
  else if (useRadioactiveEvt)
    {
      
      // initialize GPS properties
      MyGPS->ClearAll();
      
      MyGPS->SetMultipleVertex(true);
      
      std::vector<WCSimPmtInfo*> *pmts=NULL;
      
      std::vector<struct radioactive_source>::iterator it;
      
      for ( it = radioactive_sources.begin(); it != radioactive_sources.end(); it++ ){
      	G4String IsotopeName = it->IsotopeName;
      	G4String IsotopeLocation = it->IsotopeLocation;
      	G4double IsotopeActivity = it->IsotopeActivity;

      	double average= IsotopeActivity * GetRadioactiveTimeWindow();
      	if (IsotopeLocation.compareTo("PMT") == 0){
      		pmts = myDetector->Get_Pmts();
      		average *= pmts->size();
      	}
	  
	// random poisson number of vertices based on average
	int n_vertices = CLHEP::RandPoisson::shoot(average);

	//	n_vertices = 1; // qqq

	for(int u=0; u<n_vertices; u++){
	    
	  MyGPS->AddaSource(1.);
	    
	  MyGPS->SetCurrentSourceto(MyGPS->GetNumberofSource() - 1);

	  if (IsotopeName.compareTo("Tl208") == 0){
	    MyGPS->SetParticleDefinition(G4IonTable::GetIonTable()->GetIon( 81, 208, 0));
	  }else if (IsotopeName.compareTo("Bi214") == 0){
	    MyGPS->SetParticleDefinition(G4IonTable::GetIonTable()->GetIon( 83, 214, 0));
	  }else if (IsotopeName.compareTo("K40") == 0){
	    MyGPS->SetParticleDefinition(G4IonTable::GetIonTable()->GetIon( 19, 40, 0));
	  }else if (IsotopeName.compareTo("Rn220") == 0){
	    MyGPS->SetParticleDefinition(G4IonTable::GetIonTable()->GetIon( 86, 220, 0));
	  }else if (IsotopeName.compareTo("Po216") == 0){
	    MyGPS->SetParticleDefinition(G4IonTable::GetIonTable()->GetIon( 84, 216, 0));
	  }else if (IsotopeName.compareTo("Pb212") == 0){
	    MyGPS->SetParticleDefinition(G4IonTable::GetIonTable()->GetIon( 82, 212, 0));
	  }else if (IsotopeName.compareTo("Bi212") == 0){
	    MyGPS->SetParticleDefinition(G4IonTable::GetIonTable()->GetIon( 83, 212, 0));
	  }else if (IsotopeName.compareTo("Po212") == 0){
	    MyGPS->SetParticleDefinition(G4IonTable::GetIonTable()->GetIon( 84, 212, 0));
	  }else if (IsotopeName.compareTo("Rn222") == 0){
	    MyGPS->SetParticleDefinition(G4IonTable::GetIonTable()->GetIon( 86, 222, 0));
	  }else if (IsotopeName.compareTo("Po218") == 0){
	    MyGPS->SetParticleDefinition(G4IonTable::GetIonTable()->GetIon( 84, 218, 0));
	  }else if (IsotopeName.compareTo("At218") == 0){
	    MyGPS->SetParticleDefinition(G4IonTable::GetIonTable()->GetIon( 85, 218, 0));
	  }else if (IsotopeName.compareTo("Pb214") == 0){
	    MyGPS->SetParticleDefinition(G4IonTable::GetIonTable()->GetIon( 82, 214, 0));
	  }else if (IsotopeName.compareTo("Po214") == 0){
	    MyGPS->SetParticleDefinition(G4IonTable::GetIonTable()->GetIon( 84, 214, 0));
	  }else if (IsotopeName.compareTo("Tl210") == 0){
	    MyGPS->SetParticleDefinition(G4IonTable::GetIonTable()->GetIon( 81, 210, 0));
	  }else if (IsotopeName.compareTo("Pb210") == 0){
	    MyGPS->SetParticleDefinition(G4IonTable::GetIonTable()->GetIon( 82, 210, 0));
	  }else if (IsotopeName.compareTo("Bi210") == 0){
	    MyGPS->SetParticleDefinition(G4IonTable::GetIonTable()->GetIon( 83, 210, 0));
	  }else if (IsotopeName.compareTo("Po210") == 0){
	    MyGPS->SetParticleDefinition(G4IonTable::GetIonTable()->GetIon( 84, 210, 0));
	  }else if (IsotopeName.compareTo("Hg206") == 0){
	    MyGPS->SetParticleDefinition(G4IonTable::GetIonTable()->GetIon( 80, 206, 0));
	  }else if (IsotopeName.compareTo("Tl206") == 0){
	    MyGPS->SetParticleDefinition(G4IonTable::GetIonTable()->GetIon( 81, 206, 0));
	  }else if (IsotopeName.compareTo("Rn219") == 0){
	    MyGPS->SetParticleDefinition(G4IonTable::GetIonTable()->GetIon( 86, 219, 0));
	  }else if (IsotopeName.compareTo("Po215") == 0){
	    MyGPS->SetParticleDefinition(G4IonTable::GetIonTable()->GetIon( 84, 215, 0));
	  }else if (IsotopeName.compareTo("At215") == 0){
	    MyGPS->SetParticleDefinition(G4IonTable::GetIonTable()->GetIon( 85, 215, 0));
	  }else if (IsotopeName.compareTo("Pb211") == 0){
	    MyGPS->SetParticleDefinition(G4IonTable::GetIonTable()->GetIon( 82, 211, 0));
	  }else if (IsotopeName.compareTo("Bi211") == 0){
	    MyGPS->SetParticleDefinition(G4IonTable::GetIonTable()->GetIon( 83, 211, 0));
	  }else if (IsotopeName.compareTo("Po211") == 0){
	    MyGPS->SetParticleDefinition(G4IonTable::GetIonTable()->GetIon( 84, 211, 0));
	  }else if (IsotopeName.compareTo("Tl207") == 0){
	    MyGPS->SetParticleDefinition(G4IonTable::GetIonTable()->GetIon( 81, 207, 0));
	  }else if (IsotopeName.compareTo("Th232") == 0){
	    MyGPS->SetParticleDefinition(G4IonTable::GetIonTable()->GetIon( 90, 232, 0));
	  }else if (IsotopeName.compareTo("Ra228") == 0){
	    MyGPS->SetParticleDefinition(G4IonTable::GetIonTable()->GetIon( 88, 228, 0));
	  }else if (IsotopeName.compareTo("Ac228") == 0){
	    MyGPS->SetParticleDefinition(G4IonTable::GetIonTable()->GetIon( 89, 228, 0));
	  }else if (IsotopeName.compareTo("Th228") == 0){
	    MyGPS->SetParticleDefinition(G4IonTable::GetIonTable()->GetIon( 90, 228, 0));
	  }else if (IsotopeName.compareTo("Ra224") == 0){
	    MyGPS->SetParticleDefinition(G4IonTable::GetIonTable()->GetIon( 88, 224, 0));
	  }else if (IsotopeName.compareTo("U238") == 0){
	    MyGPS->SetParticleDefinition(G4IonTable::GetIonTable()->GetIon( 92, 238, 0));
	  }else if (IsotopeName.compareTo("Th234") == 0){
	    MyGPS->SetParticleDefinition(G4IonTable::GetIonTable()->GetIon( 90, 234, 0));
	  }else if (IsotopeName.compareTo("Pa234") == 0){
	    MyGPS->SetParticleDefinition(G4IonTable::GetIonTable()->GetIon( 91, 234, 0));
	  }else if (IsotopeName.compareTo("U234") == 0){
	    MyGPS->SetParticleDefinition(G4IonTable::GetIonTable()->GetIon( 92, 234, 0));
	  }else if (IsotopeName.compareTo("Th230") == 0){
	    MyGPS->SetParticleDefinition(G4IonTable::GetIonTable()->GetIon( 90, 230, 0));
	  }else if (IsotopeName.compareTo("Ra226") == 0){
	    MyGPS->SetParticleDefinition(G4IonTable::GetIonTable()->GetIon( 88, 226, 0));
	  }else if (IsotopeName.compareTo("U235") == 0){
	    MyGPS->SetParticleDefinition(G4IonTable::GetIonTable()->GetIon( 92, 235, 0));
	  }else if (IsotopeName.compareTo("Th231") == 0){
	    MyGPS->SetParticleDefinition(G4IonTable::GetIonTable()->GetIon( 90, 231, 0));
	  }else if (IsotopeName.compareTo("Pa231") == 0){
	    MyGPS->SetParticleDefinition(G4IonTable::GetIonTable()->GetIon( 91, 231, 0));
	  }else if (IsotopeName.compareTo("Ac227") == 0){
	    MyGPS->SetParticleDefinition(G4IonTable::GetIonTable()->GetIon( 89, 227, 0));
	  }else if (IsotopeName.compareTo("Th227") == 0){
	    MyGPS->SetParticleDefinition(G4IonTable::GetIonTable()->GetIon( 90, 227, 0));
	  }else if (IsotopeName.compareTo("Fr223") == 0){
	    MyGPS->SetParticleDefinition(G4IonTable::GetIonTable()->GetIon( 87, 223, 0));
	  }else if (IsotopeName.compareTo("Ra223") == 0){
	    MyGPS->SetParticleDefinition(G4IonTable::GetIonTable()->GetIon( 88, 223, 0));
	  }else if (IsotopeName.compareTo("At219") == 0){
	    MyGPS->SetParticleDefinition(G4IonTable::GetIonTable()->GetIon( 85, 219, 0));
	  }else if (IsotopeName.compareTo("Bi215") == 0){
	    MyGPS->SetParticleDefinition(G4IonTable::GetIonTable()->GetIon( 83, 215, 0));
	  }

	  if (IsotopeLocation.compareTo("water") == 0){
	    MyGPS->GetCurrentSource()->GetEneDist()->SetEnergyDisType("Mono");
	    MyGPS->GetCurrentSource()->GetEneDist()->SetMonoEnergy(0.);
	    MyGPS->GetCurrentSource()->GetPosDist()->SetPosDisType("Point");
	    MyGPS->GetCurrentSource()->GetPosDist()->SetCentreCoords(G4ThreeVector(0, 0, 0));
	    MyGPS->GetCurrentSource()->GetPosDist()->SetPosDisType("Volume");
	    MyGPS->GetCurrentSource()->GetPosDist()->SetPosDisShape("Cylinder");
	    G4String WCIDCollectionName = myDetector->GetIDCollectionName();
	    WCSimPMTObject *PMT = myDetector->GetPMTPointer(WCIDCollectionName);
	    MyGPS->GetCurrentSource()->GetPosDist()->SetRadius(myDetector->GetGeo_Dm(3)*CLHEP::cm - 2.*PMT->GetRadius());
	    MyGPS->GetCurrentSource()->GetPosDist()->SetHalfZ(myDetector->GetGeo_Dm(2)*CLHEP::cm/2. - 2.*PMT->GetRadius());
	    MyGPS->GetCurrentSource()->GetPosDist()->SetPosRot1(G4ThreeVector(1, 0, 0));
	    MyGPS->GetCurrentSource()->GetPosDist()->SetPosRot2(G4ThreeVector(0, 1, 0));

	  }
	  else if (IsotopeLocation.compareTo("PMT") == 0){
	    int npmts = pmts->size();
	    int random_pmt_id = CLHEP::RandFlat::shootInt(1,npmts);
	    WCSimPmtInfo* pmtinfo = (WCSimPmtInfo*)pmts->at( random_pmt_id - 1 );
	    G4ThreeVector random_pmt_center(pmtinfo->Get_transx()*CLHEP::cm, pmtinfo->Get_transy()*CLHEP::cm, pmtinfo->Get_transz()*CLHEP::cm);
	    double random_cos_theta = CLHEP::RandFlat::shoot(0., 1.);
	    double random_sin_theta = sqrt(1. - pow(random_cos_theta,2));
	    random_sin_theta *= (CLHEP::RandFlat::shootBit() == 0 ? -1 : 1);
	    double random_phi = CLHEP::RandFlat::shoot(0., 2.*CLHEP::pi*CLHEP::rad);
	    G4String WCIDCollectionName = myDetector->GetIDCollectionName();
	    WCSimPMTObject *PMT = myDetector->GetPMTPointer(WCIDCollectionName);
	    double PMT_radius = PMT->GetRadius();
	    double glassThickness = PMT->GetPMTGlassThickness();
	    double expose = PMT->GetExposeHeight();
	    double sphereRadius = (expose*expose+ PMT_radius*PMT_radius)/(2*expose);
	    double Rmin = sphereRadius-glassThickness;
	    double Rmax = sphereRadius;
	    double random_R = CLHEP::RandFlat::shoot(Rmin, Rmax);
	    G4ThreeVector orientation(pmtinfo->Get_orienx(), pmtinfo->Get_orieny(), pmtinfo->Get_orienz());
	    G4ThreeVector axis_1 = orientation.orthogonal();
	    G4ThreeVector axis_2 = orientation.cross(axis_1);
	    G4ThreeVector position = random_pmt_center + random_R*(orientation*random_cos_theta + axis_1*random_sin_theta*cos(random_phi) + axis_2*random_sin_theta*sin(random_phi));
	      
	    //G4cout << " random id " << random_pmt_id << " of " << npmts << " costheta " << random_cos_theta << " sintheta " << random_sin_theta << " phi " << random_phi << " WCIDCollectionName " << WCIDCollectionName << " PMT_radius " << PMT_radius << " expose " << expose << " sphereRadius " << sphereRadius << " Rmin " << Rmin << " Rmax " << Rmax << " random_R " << random_R << " orientation (" << orientation.x() << ", " << orientation.y() << ", " << orientation.z() << ") center (" << random_pmt_center.x() << ", " << random_pmt_center.y() << ", " << random_pmt_center.z() << ") position (" << position.x() << ", " << position.y() << ", " << position.z() << ") " << G4endl;
	      
	    MyGPS->GetCurrentSource()->GetEneDist()->SetEnergyDisType("Mono");
	    MyGPS->GetCurrentSource()->GetEneDist()->SetMonoEnergy(0.);
	    MyGPS->GetCurrentSource()->GetPosDist()->SetPosDisType("Point");
	    MyGPS->GetCurrentSource()->GetPosDist()->SetCentreCoords(position);
	  }
	    
	}

//	G4cout << " is " << IsotopeName << " of " << radioactive_sources.size() << " loc " << IsotopeLocation << " a " << IsotopeActivity << " nv " << n_vertices << G4endl;

      }

      G4int number_of_sources = MyGPS->GetNumberofSource();

      // this will generate several primary vertices
      MyGPS->GeneratePrimaryVertex(anEvent);

      SetNvtxs(number_of_sources);
      for( G4int u=0; u<number_of_sources; u++){
	targetpdgs[u] = 2212; //ie. proton 

      	G4ThreeVector P   =anEvent->GetPrimaryVertex(u)->GetPrimary()->GetMomentum();
      	G4ThreeVector vtx =anEvent->GetPrimaryVertex(u)->GetPosition();
      	G4int pdg         =anEvent->GetPrimaryVertex(u)->GetPrimary()->GetPDGcode();
      
      	//       G4ThreeVector dir  = P.unit();
      	G4double E         = std::sqrt((P.dot(P)));

//	G4cout << " vertex " << u << " of " << number_of_sources << " (" << vtx.x() << ", " << vtx.y() << ", " << vtx.z() << ")" << G4endl;

      	SetVtxs(u,vtx);
      	SetBeamEnergy(E,u);
      	//       SetBeamDir(dir);
      	SetBeamPDG(pdg,u);
      }

    }
  else if (useRadonEvt)
    { //G. Pronost: Add Radon (adaptation of Radioactive event)
    
      // Currently only one generator is possible
      // In order to have several, we need to find a solution for the fitting graphes (which are static currently)
      // Idea: array of fitting graphes? (each new generators having a specific ID)
      if ( !myRn222Generator ) {
      	myRn222Generator = new WCSimGenerator_Radioactivity(myDetector);
      	myRn222Generator->Configuration(fRnScenario);
      }
      
      //G4cout << " Generate radon events " << G4endl;
      // initialize GPS properties
      MyGPS->ClearAll();
      
      MyGPS->SetMultipleVertex(true);
      
      
      std::vector<struct radioactive_source>::iterator it;
      
      G4String IsotopeName = "Rn222";
      G4double IsotopeActivity = myRn222Generator->GetMeanActivity() * 1e-3; // mBq to Bq
      G4double iEventAvg = IsotopeActivity * GetRadioactiveTimeWindow();

      //G4cout << " Average " << iEventAvg << G4endl;
      // random poisson number of vertices based on average
      int n_vertices = CLHEP::RandPoisson::shoot(iEventAvg);

      if ( n_vertices < 1 ) {
      	 n_vertices = 1;
      }
      
      for(int u=0; u<n_vertices; u++){
	
	MyGPS->AddaSource(1.);	
	MyGPS->SetCurrentSourceto(MyGPS->GetNumberofSource() - 1);
	
	// Bi214 (source of electron in Rn222 decay chain, assumed to be in equilibrium)
	MyGPS->SetParticleDefinition(G4IonTable::GetIonTable()->GetIon( 83, 214, 0));
	
	// Get position (first position take few seconds to be produced, there after there is no trouble)
	//G4cout << "GetRandomVertex" << G4endl;
	G4ThreeVector position = myRn222Generator->GetRandomVertex(fRnSymmetry);
	//G4cout << "Done: " << position << G4endl;
	// energy 
	MyGPS->GetCurrentSource()->GetEneDist()->SetEnergyDisType("Mono");
	MyGPS->GetCurrentSource()->GetEneDist()->SetMonoEnergy(0.);
	    
	// position 
	MyGPS->GetCurrentSource()->GetPosDist()->SetPosDisType("Point");
	MyGPS->GetCurrentSource()->GetPosDist()->SetCentreCoords(position);

	//G4cout << u << " is " << IsotopeName << " loc " << position  << G4endl;

      }
      G4int number_of_sources = MyGPS->GetNumberofSource();

      // this will generate several primary vertices
      MyGPS->GeneratePrimaryVertex(anEvent);

      SetNvtxs(number_of_sources);
      for( G4int u=0; u<number_of_sources; u++){
	targetpdgs[u] = 2212; //ie. proton 

      	G4ThreeVector P   =anEvent->GetPrimaryVertex(u)->GetPrimary()->GetMomentum();
      	G4ThreeVector vtx =anEvent->GetPrimaryVertex(u)->GetPosition();
      	G4int pdg         =anEvent->GetPrimaryVertex(u)->GetPrimary()->GetPDGcode();
      
      	//       G4ThreeVector dir  = P.unit();
      	G4double E         = std::sqrt((P.dot(P)));
	
	//G4cout << " vertex " << u << " of " << number_of_sources << " (" << vtx.x() << ", " << vtx.y() << ", " << vtx.z() << ") with pdg: " << pdg << G4endl;

      	SetVtxs(u,vtx);
      	SetBeamEnergy(E,u);
      	//       SetBeamDir(dir);
      	SetBeamPDG(pdg,u);
      }

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

