//
// Created by moritz on 3/4/21.
//


#include <iostream>

using namespace std;

#include "WLGDSteppingAction.hh"
#include "WLGDTrackingAction.hh"
#include "WLGDRunAction.hh"

#include "G4SystemOfUnits.hh"

#include "G4RunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WLGDSteppingAction::UserSteppingAction(const G4Step *aStep) {


  if(aStep->GetTrack()->GetParticleDefinition()->GetParticleName() == "neutron")
  {
    G4cout << "___________________________________________________________________________________________________________________"
           << G4endl;
    G4cout << "___________________________________________________________________________________________________________________"
           << G4endl;
    G4cout << "___________________________________________________________________________________________________________________"
           << G4endl;
    G4cout << "___________________________________________________________________________________________________________________"
           << G4endl;
    G4cout << "___________________________________________________________________________________________________________________"
           << G4endl;
    G4cout << "___________________________________________________________________________________________________________________"
           << G4endl;
    G4cout << "___________________________________________________________________________________________________________________"
           << G4endl;
    G4cout << "Before accessing the physVolume" << G4endl;
    auto physVol = aTrack->GetVolume();
    G4cout << "Pre: " << physVol->GetName() << " - " << aTrack->GetTrackID() << G4endl;
    if(aTrack->GetNextVolume())
    {
      auto physVol2 = aTrack->GetNextVolume();
      G4cout << "Post: " << physVol2->GetName() << " - " << aTrack->GetTrackID()
             << G4endl;
    }
    // if(aTrack->GetStep()->GetPreStepPoint()->GetPhysicalVolume()->GetName() ==
    // "Ge_phys")
  }

}
