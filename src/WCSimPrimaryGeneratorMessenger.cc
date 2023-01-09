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
  genCmd->SetGuidance(" Available generators : muline, gun, laser, gps, cosmics, calibration, radioactive, radon");
  genCmd->SetParameterName("generator",true);
  genCmd->SetDefaultValue("muline");
  //T. Akiri: Addition of laser
  genCmd->SetCandidates("muline gun laser gps cosmics calibration radioactive radon");

  fileNameCmd = new G4UIcmdWithAString("/mygen/vecfile",this);
  fileNameCmd->SetGuidance("Select the file of vectors.");
  fileNameCmd->SetGuidance(" Enter the file name of the vector file");
  fileNameCmd->SetParameterName("fileName",true);
  fileNameCmd->SetDefaultValue("inputvectorfile");

  fileNameCmdCosmics = new G4UIcmdWithAString("/mygen/cosmicsfile",this);
  fileNameCmdCosmics->SetGuidance("Select the file of cosmics.");
  fileNameCmdCosmics->SetGuidance(" Enter the file name of the cosmics file");
  fileNameCmdCosmics->SetParameterName("fileName",true);

  //C. Pidcott: Addition of calibration settings
  calibrationSource = new G4UIcmdWithAString("/mygen/calibrationsource",this);
  calibrationSource->SetGuidance("Select optical calibration source.");
  calibrationSource->SetGuidance(" Available sources : fullInjectors, colBarrel, colBottom, colTop, satBarrel, satTop, satBottom, timeTop, IWCDfullInjectors");
  calibrationSource->SetParameterName("calSource",true);
  calibrationSource->SetDefaultValue("fullInjectors");
  calibrationSource->SetCandidates("fullInjectors colBarrel colBottom satBarrel satTop satBottom colTop timeTop IWCDfullInjectors");

  calSourceNumParticles = new G4UIcmdWithAnInteger("/WCSim/calibrationsource/NumCalibParticles", this);
  calSourceNumParticles->SetGuidance("Set number of photons emitted per calibration source per event");
  calSourceNumParticles->SetParameterName("NumCalibParticles", true);
  calSourceNumParticles->SetDefaultValue(5000);
  calSourceNumParticles->SetRange("NumCalibParticles>0");

  calSourceHalfAngle = new G4UIcmdWithADouble("/WCSim/calibrationsource/CalibSourceHalfAngle", this);
  calSourceHalfAngle->SetGuidance("Set half angle of light source (degrees)");
  calSourceHalfAngle->SetParameterName("CalibHalfAngle",true);
  calSourceHalfAngle->SetDefaultValue(40);

  calSourceWavelength = new G4UIcmdWithADouble("/WCSim/calibrationsource/CalibSourceWavelength", this);
  calSourceHalfAngle->SetGuidance("Set wavelength of light source (nm)");
  calSourceHalfAngle->SetParameterName("CalibWavelength",true);
  calSourceHalfAngle->SetDefaultValue(400);

  fileNameCmdCosmics->SetDefaultValue("data/MuonFlux-HyperK-ThetaPhi.dat");
  
  timeUnitCmd = new G4UIcmdWithAString("/mygen/time_unit",this);
  timeUnitCmd->SetGuidance("Define the units used for time in the input file.");
  timeUnitCmd->SetGuidance("Can be picosecond, ps, ns, nanosecond, ms, millisecond, s, sec or second");
  timeUnitCmd->SetGuidance("Default if not set is nanosecond");
  timeUnitCmd->SetParameterName("unit",true);
  timeUnitCmd->SetDefaultValue("ns");
  
  radioactive_time_window_Cmd = new G4UIcmdWithADouble("/mygen/radioactive_time_window",this);
  radioactive_time_window_Cmd->SetGuidance("Select time window for radioactivity");
  radioactive_time_window_Cmd->SetParameterName("radioactive_time_window",true);
  radioactive_time_window_Cmd->SetDefaultValue(0.);
  
  isotopeCmd = new G4UIcmdWithAString("/mygen/isotope",this);
  isotopeCmd->SetGuidance("Select properties of radioactive isotope");
  isotopeCmd->SetGuidance("[usage] /mygen/isotope ISOTOPE LOCATION ACTIVITY");
  isotopeCmd->SetGuidance("     ISOTOPE : Tl208, Bi214, K40");
  isotopeCmd->SetGuidance("     LOCATION : water PMT");
  isotopeCmd->SetGuidance("     ACTIVITY : (int) activity of isotope (Bq) ");
  G4UIparameter* param;
  param = new G4UIparameter("ISOTOPE",'s',true);
  param->SetDefaultValue("Tl208");
  isotopeCmd->SetParameter(param);
  param = new G4UIparameter("LOCATION",'s',true);
  param->SetDefaultValue("water");
  isotopeCmd->SetParameter(param);
  param = new G4UIparameter("ACTIVITY",'d',true);
  param->SetDefaultValue("0");
  isotopeCmd->SetParameter(param);
  
  radonScalingCmd = new G4UIcmdWithAString("/mygen/radon_scaling",this);
  radonScalingCmd->SetGuidance("Select scalling scenario, if scenario 0 is selected, Bi214 are generated uniformly");
  radonScalingCmd->SetGuidance("[usage] /mygen/radon SCENARIO ");
  radonScalingCmd->SetGuidance("     SCENARIO : 0, A, B");
  radonScalingCmd->SetCandidates("0 A B");
  param = new G4UIparameter("SCENARIO",'s',true);
  param->SetDefaultValue("A");
  radonScalingCmd->SetParameter(param);
  
  radonGeoSymCmd = new G4UIcmdWithAnInteger("/mygen/radon_symmetry",this);
  radonGeoSymCmd->SetGuidance("Select scalling scenario");
  radonGeoSymCmd->SetGuidance("[usage] /mygen/radon SCENARIO ");
  radonGeoSymCmd->SetGuidance("     SYMMETRY : 1 ... ");
  param = new G4UIparameter("SYMMETRY",'d',true);
  param->SetDefaultValue("1");
  radonScalingCmd->SetParameter(param);
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
      myAction->SetRadioactiveEvtGenerator(false);
      myAction->SetRadonEvtGenerator(false);
    }
    else if ( newValue == "gun")
    {
      myAction->SetMulineEvtGenerator(false);
      myAction->SetGunEvtGenerator(true);
      myAction->SetLaserEvtGenerator(false);
      myAction->SetGPSEvtGenerator(false);
      myAction->SetCosmicsGenerator(false);
      myAction->SetRadioactiveEvtGenerator(false);
      myAction->SetRadonEvtGenerator(false);
    }
    else if ( newValue == "laser")   //T. Akiri: Addition of laser
    {
      myAction->SetMulineEvtGenerator(false);
      myAction->SetGunEvtGenerator(false);
      myAction->SetLaserEvtGenerator(true);
      myAction->SetGPSEvtGenerator(false);
      myAction->SetCosmicsGenerator(false);
      myAction->SetRadioactiveEvtGenerator(false);
      myAction->SetRadonEvtGenerator(false);
    }
    else if ( newValue == "gps")
    {
      myAction->SetMulineEvtGenerator(false);
      myAction->SetGunEvtGenerator(false);
      myAction->SetLaserEvtGenerator(false);
      myAction->SetGPSEvtGenerator(true);
      myAction->SetCosmicsGenerator(false);
      myAction->SetRadioactiveEvtGenerator(false);
      myAction->SetRadonEvtGenerator(false);
    }
    else if ( newValue == "cosmics")
    {
      myAction->SetMulineEvtGenerator(false);
      myAction->SetGunEvtGenerator(false);
      myAction->SetLaserEvtGenerator(false);
      myAction->SetGPSEvtGenerator(false);
      myAction->SetCosmicsGenerator(true);
      myAction->SetRadioactiveEvtGenerator(false);
      myAction->SetRadonEvtGenerator(false);
    }
    else if ( newValue == "radioactive") //G. Pronost: Addition of Radioactivity (from F. Nova code)
    {
      myAction->SetMulineEvtGenerator(false);
      myAction->SetGunEvtGenerator(false);
      myAction->SetLaserEvtGenerator(false);
      myAction->SetGPSEvtGenerator(false);
      myAction->SetCosmicsGenerator(false);
      myAction->SetRadioactiveEvtGenerator(true);
      myAction->SetRadonEvtGenerator(false);
    }
    else if ( newValue == "radon" ) //G. Pronost: Addition of Radon generator (based on F. Nova's radioactive generator but dedicated to radioactive events in water)
    {
      myAction->SetMulineEvtGenerator(false);
      myAction->SetGunEvtGenerator(false);
      myAction->SetLaserEvtGenerator(false);
      myAction->SetGPSEvtGenerator(false);
      myAction->SetCosmicsGenerator(false);
      myAction->SetRadioactiveEvtGenerator(false);
      myAction->SetRadonEvtGenerator(true);
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

  if( command == fileNameCmd)
  {
    myAction->OpenVectorFile(newValue);
    G4cout << "Input vector file set to " << newValue << G4endl;
  }
  if(command == fileNameCmdCosmics )
  {
    myAction->OpenCosmicsFile(newValue);
    G4cout << "Input cosmics data file set to " << newValue << G4endl;
  }
  if( command==isotopeCmd )
  {
    IsotopeCommand(newValue);
  }

  if( command==radioactive_time_window_Cmd )
  {
    myAction->SetRadioactiveTimeWindow(StoD(newValue));
  }
  
  if ( command==radonScalingCmd ) 
  {
    RadonScalingCommand(newValue);
  }
  
  if ( command==radonGeoSymCmd ) 
  {
    myAction->SetRadonSymmetry(radonGeoSymCmd->GetNewIntValue(newValue));
  }

  if ( command==timeUnitCmd)
  {
    myAction->SetTimeUnit(newValue);
    G4cout << "Time unit set to " << newValue << G4endl;
  }

if( command==calibrationSource )
  {
    if (newValue == "fullInjectors")
    {
      myAction->SetFullInjectors(true);
      myAction->SetCollimatedBarrel(false);
      myAction->SetCollimatedBottom(false);
      myAction->SetCollimatedTop(false);
      myAction->SetSaturationBarrel(false);
      myAction->SetSaturationTop(false);
      myAction->SetSaturationBottom(false);
      myAction->SetTimeTop(false);
      myAction->SetIWCDFullInjectors(false);
    }
    else if ( newValue == "colBarrel")
      {
	myAction->SetFullInjectors(false);
	myAction->SetCollimatedBarrel(true);
	myAction->SetCollimatedBottom(false);
	myAction->SetCollimatedTop(false);
	myAction->SetSaturationBarrel(false);
        myAction->SetSaturationTop(false);
	myAction->SetSaturationBottom(false);
	myAction->SetTimeTop(false);
	myAction->SetIWCDFullInjectors(false);
      }
    else if ( newValue == "colBottom")
      {
	myAction->SetFullInjectors(false);
	myAction->SetCollimatedBarrel(false);
	myAction->SetCollimatedBottom(true);
	myAction->SetCollimatedTop(false);
	myAction->SetSaturationBarrel(false);
        myAction->SetSaturationTop(false);
	myAction->SetSaturationBottom(false);
        myAction->SetTimeTop(false);
	myAction->SetIWCDFullInjectors(false);
      }
    else if ( newValue == "colTop")
      {
	myAction->SetFullInjectors(false);
	myAction->SetCollimatedBarrel(false);
	myAction->SetCollimatedBottom(false);
	myAction->SetCollimatedTop(true);
	myAction->SetSaturationBarrel(false);
        myAction->SetSaturationTop(false);
	myAction->SetSaturationBottom(false);
	myAction->SetTimeTop(false);
	myAction->SetIWCDFullInjectors(false);
      }
    else if ( newValue == "satTop")
      {
	myAction->SetFullInjectors(false);
	myAction->SetCollimatedBarrel(false);
	myAction->SetCollimatedBottom(false);
	myAction->SetCollimatedTop(false);
	myAction->SetSaturationBarrel(false);
        myAction->SetSaturationTop(true);
	myAction->SetSaturationBottom(false);
	myAction->SetTimeTop(false);
	myAction->SetIWCDFullInjectors(false);
      }
    else if ( newValue == "satBarrel")
      {
	myAction->SetFullInjectors(false);
	myAction->SetCollimatedBarrel(false);
	myAction->SetCollimatedBottom(false);
	myAction->SetCollimatedTop(false);
	myAction->SetSaturationBarrel(true);
        myAction->SetSaturationTop(false);
	myAction->SetSaturationBottom(false);
	myAction->SetTimeTop(false);
	myAction->SetIWCDFullInjectors(false);
      }
    else if ( newValue == "satBottom")
      {
	myAction->SetFullInjectors(false);
	myAction->SetCollimatedBarrel(false);
	myAction->SetCollimatedBottom(false);
	myAction->SetCollimatedTop(false);
	myAction->SetSaturationBarrel(false);
        myAction->SetSaturationTop(false);
	myAction->SetSaturationBottom(true);
	myAction->SetTimeTop(false);
	myAction->SetIWCDFullInjectors(false);
      }
    else if ( newValue == "timeTop")
      {
	myAction->SetFullInjectors(false);
	myAction->SetCollimatedBarrel(false);
	myAction->SetCollimatedBottom(false);
	myAction->SetCollimatedTop(false);
	myAction->SetSaturationBarrel(false);
        myAction->SetSaturationTop(false);
	myAction->SetSaturationBottom(false);
	myAction->SetTimeTop(true);
	myAction->SetIWCDFullInjectors(false);
      }
    else if (newValue == "IWCDfullInjectors")
    {
      myAction->SetFullInjectors(false);
      myAction->SetCollimatedBarrel(false);
      myAction->SetCollimatedBottom(false);
      myAction->SetCollimatedTop(false);
      myAction->SetSaturationBarrel(false);
      myAction->SetSaturationTop(false);
      myAction->SetSaturationBottom(false);
      myAction->SetTimeTop(false);
      myAction->SetIWCDFullInjectors(true);
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

  if ( command == calSourceWavelength )
   {
     myAction->SetCalibrationSourceWavelength(calSourceWavelength->GetNewDoubleValue(newValue));
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
    else if(myAction->IsUsingRadioactiveEvtGenerator())
      { cv = "radioactive"; }
    else if(myAction->IsUsingRadonEvtGenerator())
      { cv = "radon"; }
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
    else if(myAction->IsUsingSaturationBarrel())
      { cv = "satBarrel"; }
    else if(myAction->IsUsingSaturationTop())
      { cv = "satTop"; }
    else if(myAction->IsUsingSaturationBottom())
      { cv = "satBottom"; }
    else if(myAction->IsUsingTimeTop())
      { cv = "timeTop"; }
    else if(myAction->IsUsingIWCDFullInjectors())
      { cv = "IWCDfullInjectors"; }
  }
  
  return cv;
}


void  WCSimPrimaryGeneratorMessenger::IsotopeCommand(G4String newValue)
{
  G4Tokenizer next( newValue );

  G4String isotope = next();
  G4String location = next();
  G4double activity = StoD(next());

  myAction->AddRadioactiveSource(isotope, location, activity);
}

void WCSimPrimaryGeneratorMessenger::RadonScalingCommand(G4String newValue)
{
  G4Tokenizer next( newValue );

  G4String scenario = next();
  G4int iScenario = 0;
   
  if ( scenario == "A" ) iScenario = 1; // Relative scaling with respect to full ID volume (Pessimistic)
  if ( scenario == "B" ) iScenario = 2; // Relative scaling with respect to fiducial volume
   
  myAction->SetRadonScenario(iScenario);
}
