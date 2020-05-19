#include "WCSimPrimaryGeneratorMessenger.hh"
#include "WCSimPrimaryGeneratorAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4ios.hh"

WCSimPrimaryGeneratorMessenger::WCSimPrimaryGeneratorMessenger(WCSimPrimaryGeneratorAction* pointerToAction)
:myAction(pointerToAction)
{
  mydetDirectory = new G4UIdirectory("/mygen/");
  mydetDirectory->SetGuidance("WCSim detector control commands.");

  genCmd = new G4UIcmdWithAString("/mygen/generator",this);
  genCmd->SetGuidance("Select primary generator.");
  //T. Akiri: Addition of laser
  genCmd->SetGuidance(" Available generators : muline, gun, laser, gps, cosmics");
  genCmd->SetParameterName("generator",true);
  genCmd->SetDefaultValue("muline");
  //T. Akiri: Addition of laser
  genCmd->SetCandidates("muline gun laser gps cosmics");
  //C. Pidcott: Addition of calibrations
  genCmd->SetCandidates("muline gun laser gps cosmics calibration");

  fileNameCmd = new G4UIcmdWithAString("/mygen/vecfile",this);
  fileNameCmd->SetGuidance("Select the file of vectors.");
  fileNameCmd->SetGuidance(" Enter the file name of the vector file");
  fileNameCmd->SetParameterName("fileName",true);
  fileNameCmd->SetDefaultValue("inputvectorfile");

  fileNameCmdCosmics = new G4UIcmdWithAString("/mygen/cosmicsfile",this);
  fileNameCmdCosmics->SetGuidance("Select the file of cosmics.");
  fileNameCmdCosmics->SetGuidance(" Enter the file name of the cosmics file");
  fileNameCmdCosmics->SetParameterName("fileName",true);
  fileNameCmdCosmics->SetDefaultValue("inputvectorfile");

    //C. Pidcott: Addition of calibration settings
  calibrationSource = new G4UIcmdWithAString("/mygen/calibrationsource",this);
  calibrationSource->SetGuidance("Select optical calibration source.");
  calibrationSource->SetGuidance(" Available sources : fullInjectors, colBarrel, colBottom, colTop, timeTop");
  calibrationSource->SetParameterName("calSource",true);
  calibrationSource->SetDefaultValue("fullInjectors");
  calibrationSource->SetCandidates("fullInjectors colBarrel colBottom colTop timeTop");

  calSourceNumParticles = new G4UIcmdWithAnInteger("/WCSim/calibrationsource/NumCalibParticles", this);
  calSourceNumParticles->SetGuidance("Set number of photons emitted per calibration source per event");
  calSourceNumParticles->SetParameterName("NumCalibParticles", true);
  calSourceNumParticles->SetDefaultValue(5000);
  calSourceNumParticles->SetRange("NumCalibParticles>0");

  calSourceHalfAngle = new G4UIcmdWithADouble("/WCSim/calibrationsource/CalibSourceHalfAngle", this);
  calSourceHalfAngle->SetGuidance("Set half angle of light source (degrees)");
  calSourceHalfAngle->SetParameterName("CalibHalfAngle",true);
  calSourceHalfAngle->SetDefaultValue(40);


}

WCSimPrimaryGeneratorMessenger::~WCSimPrimaryGeneratorMessenger()
{
  delete genCmd;
  delete mydetDirectory;
}

void WCSimPrimaryGeneratorMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{
  if( command==genCmd )
  {
    if (newValue == "muline")
    {
      myAction->SetMulineEvtGenerator(true);
      myAction->SetGunEvtGenerator(false);
      myAction->SetLaserEvtGenerator(false);
      myAction->SetGPSEvtGenerator(false);
      myAction->SetCosmicsGenerator(false);
    }
    else if ( newValue == "gun")
    {
      myAction->SetMulineEvtGenerator(false);
      myAction->SetGunEvtGenerator(true);
      myAction->SetLaserEvtGenerator(false);
      myAction->SetGPSEvtGenerator(false);
      myAction->SetCosmicsGenerator(false);
    }
    else if ( newValue == "laser")   //T. Akiri: Addition of laser
    {
      myAction->SetMulineEvtGenerator(false);
      myAction->SetGunEvtGenerator(false);
      myAction->SetLaserEvtGenerator(true);
      myAction->SetGPSEvtGenerator(false);
      myAction->SetCosmicsGenerator(false);
    }
    else if ( newValue == "gps")
    {
      myAction->SetMulineEvtGenerator(false);
      myAction->SetGunEvtGenerator(false);
      myAction->SetLaserEvtGenerator(false);
      myAction->SetGPSEvtGenerator(true);
      myAction->SetCosmicsGenerator(false);
    }
    else if ( newValue == "cosmics")
    {
      myAction->SetMulineEvtGenerator(false);
      myAction->SetGunEvtGenerator(false);
      myAction->SetLaserEvtGenerator(false);
      myAction->SetGPSEvtGenerator(false);
      myAction->SetCosmicsGenerator(true);
    }
    else if ( newValue == "calibration") //C. Pidcott: Addition of calibrations
    {
      myAction->SetMulineEvtGenerator(false);
      myAction->SetGunEvtGenerator(false);
      myAction->SetLaserEvtGenerator(false);
      myAction->SetGPSEvtGenerator(false);
      myAction->SetCosmicsGenerator(false);
      myAction->SetCalibrationGenerator(true);
    }
  }

  if( command == fileNameCmd || command == fileNameCmdCosmics )
  {
    myAction->OpenVectorFile(newValue);
    G4cout << "Input vector file set to " << newValue << G4endl;
  }

if( command==calibrationSource )
  {
    if (newValue == "fullInjectors")
    {
      myAction->SetFullInjectors(true);
      myAction->SetCollimatedBarrel(false);
      myAction->SetCollimatedBottom(false);
      myAction->SetCollimatedTop(false);
    }
    else if ( newValue == "colBarrel")
      {
	myAction->SetFullInjectors(false);
	myAction->SetCollimatedBarrel(true);
	myAction->SetCollimatedBottom(false);
	myAction->SetCollimatedTop(false);
      }
    else if ( newValue == "colBottom")
      {
	myAction->SetFullInjectors(false);
	myAction->SetCollimatedBarrel(false);
	myAction->SetCollimatedBottom(true);
	myAction->SetCollimatedTop(false);
      }
    else if ( newValue == "colTop")
      {
	myAction->SetFullInjectors(false);
	myAction->SetCollimatedBarrel(false);
	myAction->SetCollimatedBottom(false);
	myAction->SetCollimatedTop(true);
      }
    else if ( newValue == "timeTop")
      {
	myAction->SetFullInjectors(false);
	myAction->SetCollimatedBarrel(false);
	myAction->SetCollimatedBottom(false);
	myAction->SetCollimatedTop(false);
	myAction->SetTimeTop(true);
      }

  }

 if ( command == calSourceNumParticles )
   {
     myAction->SetNumberCalibrationParticles(calSourceNumParticles->GetNewIntValue(newValue));
   }

 if ( command == calSourceHalfAngle )
   {
     myAction->SetCalibrationSourceHalfAngle(calSourceHalfAngle->GetNewDoubleValue(newValue));
   }

}

G4String WCSimPrimaryGeneratorMessenger::GetCurrentValue(G4UIcommand* command)
{
  G4String cv;
  
  if( command==genCmd )
  {
    if(myAction->IsUsingMulineEvtGenerator())
      { cv = "muline"; }
    else if(myAction->IsUsingGunEvtGenerator())
      { cv = "gun"; }
    else if(myAction->IsUsingLaserEvtGenerator())
      { cv = "laser"; }   //T. Akiri: Addition of laser
    else if(myAction->IsUsingGPSEvtGenerator())
      { cv = "gps"; }
    else if(myAction->IsUsingCosmicsGenerator())
      { cv = "cosmics"; }
    else if(myAction->IsUsingCalibrationGenerator())
      { cv = "calibration"; } //C. Pidcott: Addition of calibrations
  }

  if( command==calibrationSource )
  {
    if(myAction->IsUsingFullInjectors())
      { cv = "fullInjectors"; }
    else if(myAction->IsUsingCollimatedBarrel())
      { cv = "colBarrel"; }
    else if(myAction->IsUsingCollimatedBottom())
      { cv = "colBottom"; }
    else if(myAction->IsUsingCollimatedTop())
      { cv = "colTop"; }
    else if(myAction->IsUsingTimeTop())
      { cv = "timeTop"; }
  }
  
  return cv;
}

