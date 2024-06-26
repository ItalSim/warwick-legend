#ifndef WLGDTrackingAction_h
#define WLGDTrackingAction_h 1

#include "G4UserTrackingAction.hh"
#include "WLGDEventAction.hh"
#include "WLGDRunAction.hh"
#include <map>

class WLGDTrackingAction : public G4UserTrackingAction
{
public:
  WLGDTrackingAction();
  WLGDTrackingAction(WLGDEventAction* event, WLGDRunAction* run)
  {
    fEventAction = event;
    fRunAction   = run;
    WLGDTrackingAction();
  }
  virtual ~WLGDTrackingAction(){};

  virtual void PreUserTrackingAction(const G4Track*);
  virtual void PostUserTrackingAction(const G4Track*);

private:
  double           tmp_neutronXpos, tmp_neutronYpos, tmp_neutronZpos;
  double           tmp_neutronXmom, tmp_neutronYmom, tmp_neutronZmom;
  double           tmp_neutronTime;
  double           tmp_MuonXpos, tmp_MuonYpos, tmp_MuonZpos;
  double           tmp_MuonXmom, tmp_MuonYmom, tmp_MuonZmom;
  double           tmp_MuonEnergy;
  WLGDEventAction* fEventAction;
  WLGDRunAction*   fRunAction;

  G4double           X;
  G4double           Y;
  G4double           Z;
  G4double           px;
  G4double           py;
  G4double           pz;
  G4double           Time;
  G4double           KineticEnergy;
  G4double           Wavelength;
  G4int              TrackID;
  G4int              EventID;
  G4String           CreatorProcess;
  G4String           Material;
  G4String           Volume;
};

#endif
