#ifndef WCSimPrimaryGeneratorAction_h
#define WCSimPrimaryGeneratorAction_h

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

#include <fstream>

#include "WCSimRootOptions.hh"

#include "TH1D.h"
#include "TH2D.h"

class WCSimDetectorConstruction;
class G4ParticleGun;
class G4GeneralParticleSource;
class G4Event;
class WCSimPrimaryGeneratorMessenger;

class WCSimPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  WCSimPrimaryGeneratorAction(WCSimDetectorConstruction*);
  ~WCSimPrimaryGeneratorAction();
  
public:
  void GeneratePrimaries(G4Event* anEvent);

  // Gun, laser & gps setting calls these functions to fill jhfNtuple and Root tree
  void SetVtx(G4ThreeVector i)     { vtx = i; };
  void SetBeamEnergy(G4double i)   { beamenergy = i; };
  void SetBeamDir(G4ThreeVector i) { beamdir = i; };
  void SetBeamPDG(G4int i)         { beampdg = i; };

  // These go with jhfNtuple
  G4int GetVecRecNumber(){return vecRecNumber;}
  G4int GetMode() {return mode;};
  G4int GetVtxVol() {return vtxvol;};
  G4ThreeVector GetVtx() {return vtx;}
  G4int GetNpar() {return npar;};
  G4int GetBeamPDG() {return beampdg;};
  G4double GetBeamEnergy() {return beamenergy;};
  G4ThreeVector GetBeamDir() {return beamdir;};
  G4int GetTargetPDG() {return targetpdg;};
  G4double GetTargetEnergy() {return targetenergy;};
  G4ThreeVector GetTargetDir() {return targetdir;};

  // older ...
  G4double GetNuEnergy() {return nuEnergy;};
  G4double GetEnergy() {return energy;};
  G4double GetXPos() {return xPos;};
  G4double GetYPos() {return yPos;};
  G4double GetZPos() {return zPos;};
  G4double GetXDir() {return xDir;};
  G4double GetYDir() {return yDir;};
  G4double GetZDir() {return zDir;};

  G4String GetGeneratorTypeString();
  
  void SaveOptionsToOutput(WCSimRootOptions * wcopt);

private:
  WCSimDetectorConstruction*      myDetector;
  G4ParticleGun*                  particleGun;
  G4GeneralParticleSource*        MyGPS;  //T. Akiri: GPS to run Laser
  G4ParticleGun*                  particleGunBarrel;
  G4ParticleGun*                  particleGunBottom;
  G4ParticleGun*                  particleGunTop;
  WCSimPrimaryGeneratorMessenger* messenger;

  // Variables set by the messenger
  G4bool   useMulineEvt;
  G4bool   useGunEvt;
  G4bool   useLaserEvt;  //T. Akiri: Laser flag
  G4bool   useGPSEvt;
  G4bool   useCosmics;
  G4bool   useCalibration; //C. Pidcott: Calibrations flag
  G4bool   useFullInjectors;
  G4bool   useColBarrel;
  G4bool   useColBottom;
  G4bool   useColTop;
  G4bool   useSatTop;
  G4bool   useSatBarrel;
  G4bool   useTimeTop;
  G4bool   useIWCDFullInjectors;
  std::fstream inputFile;
  G4String vectorFileName;
  G4bool   GenerateVertexInRock;

  // These go with jhfNtuple
  G4int mode;
  G4int vtxvol;
  G4ThreeVector vtx;
  G4int npar;
  G4int beampdg, targetpdg;
  G4ThreeVector beamdir, targetdir;
  G4double beamenergy, targetenergy;
  G4int vecRecNumber;

  G4double nuEnergy;
  G4double energy;
  G4double xPos, yPos, zPos;
  G4double xDir, yDir, zDir;

  G4int    _counterRock; 
  G4int    _counterCublic; 

  // Use Histograms to generate cosmics
  TH2D *hFluxCosmics;
  TH2D *hEmeanCosmics;

  // Set cosmics altitude
  G4double altCosmics;

  //C. Pidcott: Calibrations variables
  G4int c_particle;
  G4double ledtheta;

 public:

  inline void SetMulineEvtGenerator(G4bool choice) { useMulineEvt = choice; }
  inline G4bool IsUsingMulineEvtGenerator() { return useMulineEvt; }

  inline void SetGunEvtGenerator(G4bool choice) { useGunEvt = choice; }
  inline G4bool IsUsingGunEvtGenerator()  { return useGunEvt; }

  //T. Akiri: Addition of function for the laser flag
  inline void SetLaserEvtGenerator(G4bool choice) { useLaserEvt = choice; }
  inline G4bool IsUsingLaserEvtGenerator()  { return useLaserEvt; }

  inline void SetGPSEvtGenerator(G4bool choice) { useGPSEvt = choice; }
  inline G4bool IsUsingGPSEvtGenerator()  { return useGPSEvt; }

  inline void SetCosmicsGenerator(G4bool choice) { useCosmics = choice; }
  inline G4bool IsUsingCosmicsGenerator()  { return useCosmics; }

  //C. Pidcott: Addition of function for the calibration flag
  inline void SetCalibrationGenerator(G4bool choice) { useCalibration = choice; }
  inline G4bool IsUsingCalibrationGenerator()  { return useCalibration; }

  inline void SetFullInjectors(G4bool choice) { useFullInjectors = choice; }
  inline G4bool IsUsingFullInjectors()  { return useFullInjectors; }

  inline void SetCollimatedBarrel(G4bool choice) { useColBarrel = choice; }
  inline G4bool IsUsingCollimatedBarrel()  { return useColBarrel; }
  
  inline void SetCollimatedBottom(G4bool choice) { useColBottom = choice; }
  inline G4bool IsUsingCollimatedBottom()  { return useColBottom; }

  inline void SetCollimatedTop(G4bool choice) { useColTop = choice; }
  inline G4bool IsUsingCollimatedTop()  { return useColTop; }

  inline void SetSaturationBarrel(G4bool choice) { useSatBarrel = choice; }
  inline G4bool IsUsingSaturationBarrel()  { return useSatBarrel; }

  inline void SetSaturationTop(G4bool choice) { useSatTop = choice; }
  inline G4bool IsUsingSaturationTop()  { return useSatTop; }

  inline void SetTimeTop(G4bool choice) { useTimeTop = choice; }
  inline G4bool IsUsingTimeTop()  { return useTimeTop; }

  inline void SetIWCDFullInjectors(G4bool choice) { useIWCDFullInjectors = choice; }
  inline G4bool IsUsingIWCDFullInjectors()  { return useIWCDFullInjectors; }

  inline void SetNumberCalibrationParticles(G4int choice) { c_particle = choice; }
  inline void SetCalibrationSourceHalfAngle(G4int choice) { ledtheta = choice; }

  inline void OpenVectorFile(G4String fileName) 
  {
    if ( inputFile.is_open() ) 
      inputFile.close();

    vectorFileName = fileName;
    inputFile.open(vectorFileName, std::fstream::in);

    if ( !inputFile.is_open() ) {
      G4cout << "Vector file " << vectorFileName << " not found" << G4endl;
      exit(-1);
    }
  }
  inline G4bool IsGeneratingVertexInRock() { return GenerateVertexInRock; }
  inline void SetGenerateVertexInRock(G4bool choice) { GenerateVertexInRock = choice; }

};

#endif


