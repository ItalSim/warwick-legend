//
// Created by moritz on 3/4/21.
//

#ifndef WARWICK_LEGEND_WLGDSTEPPINGACTION_HH
#define WARWICK_LEGEND_WLGDSTEPPINGACTION_HH

#include "G4GenericMessenger.hh"
#include "G4UserSteppingAction.hh"
#include "WLGDDetectorConstruction.hh"
#include "WLGDRunAction.hh"
#include "globals.hh"

class WLGDTrackingAction;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class WLGDSteppingAction : public G4UserSteppingAction
{
public:
  WLGDSteppingAction() = default;
  WLGDSteppingAction(WLGDEventAction* event, WLGDRunAction* run,
                     WLGDDetectorConstruction* det);
  virtual ~WLGDSteppingAction() = default;

  virtual void UserSteppingAction(const G4Step*);
  void         GetDepositionInfo(G4int answer);
  void         GetIndividualDepositionInfo(G4int answer);
  void         AllowForLongTimeEmissionReadout(G4int answer);
  void         DefineCommands();

private:
  WLGDRunAction*            fRunAction;
  WLGDEventAction*          fEventAction;
  WLGDDetectorConstruction* fDetectorConstruction;
  G4GenericMessenger*       fStepMessenger;
  G4int                     fDepositionInfo                  = 0;
  G4int                     fIndividualDepositionInfo        = 0;
  G4int                     fAllowForLongTimeEmissionReadout = 0;
  G4double                  X;
  G4double                  Y;
  G4double                  Z;
  G4float                   X1;
  G4float                   Y1;
  G4float                   Z1;
  G4float                   X2;
  G4float                   Y2;
  G4float                   Z2;
  G4double                  Time;
  G4double                  KineticEnergy;
  G4double                  EDep;
  G4int                     TrackID;
  G4int                     StepID;
  G4int                     EventID;
  G4int                     PID;
  G4String                  Process;
  G4String                  CreatorProcess;
  G4String                  Material;
  G4String                  Volume;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
