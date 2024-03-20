#include "WLGDTrackingAction.hh"
// #include "WLGDTrackInformation.hh"
#include "WLGDTrajectory.hh"

#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4TrackingManager.hh"
#include "G4EventManager.hh"
#include "G4UnitsTable.hh"
#include <map>
#include "g4root.hh"

WLGDTrackingAction::WLGDTrackingAction() = default;

void WLGDTrackingAction::PreUserTrackingAction(const G4Track* aTrack)
{
  // Create trajectory for track if requested
  if(fpTrackingManager->GetStoreTrajectory() > 0)
    {
      fpTrackingManager->SetTrajectory(new WLGDTrajectory(aTrack));
    }
  
  // add Ge77 events to ListOfGe77
  if(aTrack->GetParticleDefinition()->GetAtomicMass() == 77 &&
     aTrack->GetParticleDefinition()->GetPDGCharge() == 32)
    {
      fEventAction->AddIDListOfGe77(aTrack->GetTrackID());
    }
  
  // Adding tracking of initial muons
  if(aTrack->GetParticleDefinition()->GetParticleName() == "mu-" ||
     aTrack->GetParticleDefinition()->GetParticleName() == "mu+")
    {
      auto tmp_vector = aTrack->GetVertexPosition();
      tmp_MuonXpos    = tmp_vector.getX() / m;
      tmp_MuonYpos    = tmp_vector.getY() / m;
      tmp_MuonZpos    = tmp_vector.getZ() / m;
      tmp_vector      = aTrack->GetMomentumDirection();
      tmp_MuonXmom    = tmp_vector.getX();
      tmp_MuonYmom    = tmp_vector.getY();
      tmp_MuonZmom    = tmp_vector.getZ();
      tmp_MuonZmom    = tmp_vector.getZ();
      tmp_MuonEnergy  = aTrack->GetKineticEnergy();
    }

  // Edit: 2021/03/30 by Moritz Neuberger
  // Adding tracking of neutrons being later captured by Ge-76 as well as general produced
  // in LAr
  if(aTrack->GetParticleDefinition()->GetParticleName() == "neutron")
    {
      {
	auto tmp_vector = aTrack->GetVertexPosition();
	tmp_neutronXpos = tmp_vector.getX() / m;
	tmp_neutronYpos = tmp_vector.getY() / m;
	tmp_neutronZpos = tmp_vector.getZ() / m;
	tmp_vector      = aTrack->GetMomentumDirection();
	tmp_neutronXmom = tmp_vector.getX();
	tmp_neutronYmom = tmp_vector.getY();
	tmp_neutronZmom = tmp_vector.getZ();
        
      }  // readout of neutron production info

      // if all neutron info is readout, sent to output
      if(fRunAction->getWriteOutAllNeutronInfoRoot() == 1)
	{
	  fEventAction->AddNeutronxLoc(tmp_neutronXpos);
	  fEventAction->AddNeutronyLoc(tmp_neutronYpos);
	  fEventAction->AddNeutronzLoc(tmp_neutronZpos);
	  fEventAction->AddNeutronxMom(tmp_neutronXmom);
	  fEventAction->AddNeutronyMom(tmp_neutronYmom);
	  fEventAction->AddNeutronzMom(tmp_neutronZmom);
	  fEventAction->AddNeutronTime(tmp_neutronTime);
	  fEventAction->AddNeutronEkin(aTrack->GetKineticEnergy() / eV);
	  fEventAction->AddNeutronID(aTrack->GetTrackID());
	  fEventAction->AddNeutronEventID(G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID());
	  G4String volumeName = aTrack->GetVolume()
	    ->GetLogicalVolume()
	    ->GetName();
	  G4int whichVolume = -9999;

	  if(volumeName == "Ge_log")
	    whichVolume = 1;
	  if(volumeName == "ULar_log")
	    whichVolume = 0;
	  if(volumeName == "Copper_log")
	    whichVolume = -1;
	  if(volumeName == "BoratedPET_Logical")
	    whichVolume = -2;
	  if(volumeName == "BoratedPET_Logical_Lid")
	    whichVolume = -3;
	  if(volumeName == "Lar_log")
	    whichVolume = -4;
	  if(volumeName == "Cinn_log")
	    whichVolume = -5;
	  if(volumeName == "Cvac_log")
	    whichVolume = -6;
	  if(volumeName == "Cout_log")
	    whichVolume = -7;
	  if(volumeName == "Lid_log")
	    whichVolume = -8;
	  if(volumeName == "Bot_log")
	    whichVolume = -9;
	  if(volumeName == "Water_log")
	    whichVolume = -10;    
	  fEventAction->AddNeutronVolume(whichVolume);
	}

      // initial value for furthest position of neutron away from center (for testing, can
      // be removed)
      fEventAction->SetMostOuterRadius(
				       sqrt(tmp_neutronXpos * tmp_neutronXpos + tmp_neutronYpos * tmp_neutronYpos));

      // increase number of neutrons produced
      if(fRunAction->getWriteOutGeneralNeutronInfo() == 1)
	fEventAction->IncreaseByOne_NeutronInEvent();

      // if production of neutrons in Ge76 nC is recorded, channel to output
      if(fRunAction->getWriteOutNeutronProductionInfo() == 1)
	{
	  fRunAction->increaseTotalNumberOfNeutronsInLAr();
	  fRunAction->addEventNumber(
				     G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID());
	  fRunAction->addCoordinatsToFile(tmp_neutronXpos, tmp_neutronYpos, tmp_neutronZpos);
	  fRunAction->addMomentumToFile(tmp_neutronXmom, tmp_neutronYmom, tmp_neutronZmom);
	  fRunAction->addEnergyToFile(aTrack->GetKineticEnergy() / eV);
	  fRunAction->addParentParticleType(
					    fEventAction->neutronProducerMap.find(aTrack->GetParentID())->second);
	}
    
    }//If particle is a neutron

  //If optical track-level output is enabled, fill the Optics ntuple

  if(fRunAction->getWriteOutOpticalData())
    {
      if(aTrack->GetDefinition()->GetPDGEncoding()==-22)
	{
	  //G4cout << G4endl << "XYZ" << G4endl;
	  X = aTrack->GetVertexPosition().x()/mm;
	  Y = aTrack->GetVertexPosition().y()/mm;
	  Z = aTrack->GetVertexPosition().z()/mm;
	  px = aTrack->GetMomentumDirection().x();
	  py = aTrack->GetMomentumDirection().y();
	  pz = aTrack->GetMomentumDirection().z();
	  //G4cout << G4endl << "Time,KE, Wave" << G4endl;
	  Time = aTrack->GetGlobalTime()/ns;
	  KineticEnergy = aTrack->GetVertexKineticEnergy()/eV;
	  Wavelength = 1239.8/KineticEnergy;//in nm
	  //G4cout << G4endl << "TrackID,EventID" << G4endl;
	  TrackID = aTrack->GetTrackID();
	  EventID = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
	  //G4cout << G4endl << "GetCreatorProcess" << G4endl;
	  if(aTrack->GetCreatorProcess())
	    CreatorProcess = aTrack->GetCreatorProcess()->GetProcessName();
	  else
	    CreatorProcess = "Primary";
	  //G4cout << G4endl << "Material" << G4endl;
	  //if(aTrack->GetMaterial())
	  //G4cout <<"Crasho?" << G4endl;
	  //Material      = aTrack->GetMaterial()->GetName();
	  //else
	  //Material    = "???";
	  //G4cout << G4endl << "Volume" << G4endl;
	  if(aTrack->GetVolume())
	    Volume        = aTrack->GetVolume()->GetName();

	  //Fill the ntuple

	  auto *am = G4AnalysisManager::Instance();
          am->FillNtupleDColumn(2,0,X);
          am->FillNtupleDColumn(2,1,Y);
          am->FillNtupleDColumn(2,2,Z);
          am->FillNtupleDColumn(2,3,Time);
          am->FillNtupleDColumn(2,4,KineticEnergy);
	  am->FillNtupleDColumn(2,5,Wavelength);
          am->FillNtupleIColumn(2,6,TrackID);
          am->FillNtupleIColumn(2,7,EventID);
          am->FillNtupleSColumn(2,8,CreatorProcess);
          am->FillNtupleSColumn(2,9,Volume);
	  am->FillNtupleDColumn(2,10,px);
	  am->FillNtupleDColumn(2,11,py);
	  am->FillNtupleDColumn(2,12,pz);
          am->AddNtupleRow(2);

	}      
    }
}//PreUserTrackingAction

void WLGDTrackingAction::PostUserTrackingAction(const G4Track* aTrack)
{
  // for tracking of particles creatd in Gd interactions
  if(fRunAction->getIndividualGdDepositionInfo())
    {
      if(fEventAction->GetIDListOfGdSiblingParticles().count(aTrack->GetParentID()))
	{
	  fEventAction->AddIDListOfGdSiblingParticles(aTrack->GetTrackID());
	}
    }

  
  // for tracking of sibling and secundary particles of Ge76 nC
  if(fEventAction->GetIDListOfGe77SiblingParticles().count(aTrack->GetParentID()) &&
     fRunAction->getIndividualGeDepositionInfo())
    {
      fEventAction->AddIDListOfGe77SiblingParticles(aTrack->GetTrackID());
    }
  if(fEventAction->GetIDListOfGe77().count(aTrack->GetParentID()) &&
     fRunAction->getIndividualGeDepositionInfo())
    {
      fEventAction->AddIDListOfGe77(aTrack->GetTrackID());
    }

  // write to output muon information

  
  if(aTrack->GetTrackID() == 1)
    {
      fEventAction->AddMuonxLoc(tmp_MuonXpos);
      fEventAction->AddMuonyLoc(tmp_MuonYpos);
      fEventAction->AddMuonzLoc(tmp_MuonZpos);
      fEventAction->AddMuonxMom(tmp_MuonXmom);
      fEventAction->AddMuonyMom(tmp_MuonYmom);
      fEventAction->AddMuonzMom(tmp_MuonZmom);
      fEventAction->AddMuonEnergy(tmp_MuonEnergy);
    }

  
  if(fRunAction->getWriteOutNeutronProductionInfo() == 1)
    {
      G4TrackVector* secondaries = fpTrackingManager->GimmeSecondaries();
      if(secondaries != nullptr)
	{
	  size_t nSeco = secondaries->size();
	  if(nSeco > 0)
	    {
	      for(size_t i = 0; i < nSeco; i++)
		{
		  // Edit: 2021/03/30 by Moritz Neuberger
		  // Adding map of parent particles that create neutrons used above
		  if((*secondaries)[i]->GetParticleDefinition()->GetParticleName() == "neutron")
		    {
		      fEventAction->neutronProducerMap.insert(
							      std::make_pair((int) aTrack->GetTrackID(),
									     (int) aTrack->GetParticleDefinition()->GetPDGEncoding()));
		    }
		}
	    }
	}
    }
  
  // For Ge77m IC readout
  if(aTrack->GetParticleDefinition()->GetPDGEncoding() == 1000320771)
    {
      fEventAction->SetisMetastable(1);
      int NumberOfSecundaries = aTrack->GetStep()->GetSecondaryInCurrentStep()->size();
      for(int i = 0; i < NumberOfSecundaries; i++)
	{
	  if(aTrack->GetStep()
             ->GetSecondaryInCurrentStep()
             ->at(i)
             ->GetParticleDefinition()
             ->GetParticleName() == "gamma" &&
	     abs(aTrack->GetStep()->GetSecondaryInCurrentStep()->at(i)->GetTotalEnergy() /
		 eV -
		 160e3) < 1e3)
	    fEventAction->SetisIC(1);
	}
    }

  // Edit: 2021/03/30 by Moritz Neuberger

  // Adding tracking of nC on different nuclei
  if(aTrack->GetParticleDefinition()->GetParticleName() == "neutron")
    {
      if(aTrack->GetStep()->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName() ==
	 "biasWrapper(nCapture)")  // altered name necessary due to biasing
	{
	  int NumberOfSecondaries = aTrack->GetStep()->GetSecondaryInCurrentStep()->size();
	  //G4cout << NumberOfSecondaries  << G4endl;
	  for(int i = 0; i < NumberOfSecondaries; i++)
	    {
	      if(aTrack->GetStep()
		 ->GetSecondaryInCurrentStep()
		 ->at(i)
		 ->GetParticleDefinition()
		 ->GetAtomicMass() == 77 &&
		 aTrack->GetStep()
		 ->GetSecondaryInCurrentStep()
		 ->at(i)
		 ->GetParticleDefinition()
		 ->GetPDGCharge() == 32)
		{
		  fEventAction->AddEkin(aTrack->GetStep()->GetPreStepPoint()->GetKineticEnergy() / eV);
		}
        
	      fEventAction->AddnCOther_timing(
					      aTrack->GetStep()->GetPostStepPoint()->GetGlobalTime() / s);
	      fEventAction->AddnCOther_x(
					 aTrack->GetStep()->GetPostStepPoint()->GetPosition().getX() / m);
	      fEventAction->AddnCOther_y(
					 aTrack->GetStep()->GetPostStepPoint()->GetPosition().getY() / m);
	      fEventAction->AddnCOther_z(
					 aTrack->GetStep()->GetPostStepPoint()->GetPosition().getZ() / m);
	      fEventAction->AddnCOther_A(aTrack->GetStep()
					 ->GetSecondaryInCurrentStep()
					 ->at(i)
					 ->GetParticleDefinition()
					 ->GetAtomicMass());
	      fEventAction->AddnCOther_ZC(aTrack->GetStep()
					  ->GetSecondaryInCurrentStep()
					  ->at(i)
					  ->GetParticleDefinition()
					  ->GetPDGCharge());
	      // fEventAction->AddnCOther_ID(aTrack->GetParentID());
	      fEventAction->AddnCOther_ID(aTrack->GetTrackID());
	      fEventAction->AddnCOther_EventID(G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID());

	      G4String volumeName = aTrack->GetVolume()
		->GetLogicalVolume()
		->GetName();
	      G4int whichVolume = -9999;

	      if(volumeName == "Ge_log")
		whichVolume = 1;
	      if(volumeName == "ULar_log")
		whichVolume = 0;
	      if(volumeName == "Copper_log")
		whichVolume = -1;
	      if(volumeName == "BoratedPET_Logical")
		whichVolume = -2;
	      if(volumeName == "BoratedPET_Logical_Lid")
		whichVolume = -3;
	      if(volumeName == "Lar_log")
		whichVolume = -4;
	      if(volumeName == "Cinn_log")
		whichVolume = -5;
	      if(volumeName == "Cvac_log")
		whichVolume = -6;
	      if(volumeName == "Cout_log")
		whichVolume = -7;
	      if(volumeName == "Lid_log")
		whichVolume = -8;
	      if(volumeName == "Bot_log")
		whichVolume = -9;
	      if(volumeName == "Water_log")
		whichVolume = -10;    
	      fEventAction->AddnCOther_Volume(whichVolume);
	    }
	}
    }//If neutron

  
  if(fEventAction->isAllProductions())
    {
    int NumberOfSecondaries = aTrack->GetStep()->GetSecondaryInCurrentStep()->size();
    for(int i = 0; i < NumberOfSecondaries; i++)
      {
	if(aTrack->GetStep()
	   ->GetSecondaryInCurrentStep()
	   ->at(i)
	   ->GetParticleDefinition()
	   ->GetAtomicMass() > 0){
	  fEventAction->Addprod_timing(
				       aTrack->GetStep()->GetPostStepPoint()->GetGlobalTime() / s);
	  fEventAction->Addprod_x(
				  aTrack->GetStep()->GetPostStepPoint()->GetPosition().getX() / m);
	  fEventAction->Addprod_y(
				  aTrack->GetStep()->GetPostStepPoint()->GetPosition().getY() / m);
	  fEventAction->Addprod_z(
				  aTrack->GetStep()->GetPostStepPoint()->GetPosition().getZ() / m);
	  fEventAction->Addprod_A(aTrack->GetStep()
				  ->GetSecondaryInCurrentStep()
				  ->at(i)
				  ->GetParticleDefinition()
				  ->GetAtomicMass());
	  fEventAction->Addprod_ZC(aTrack->GetStep()
				   ->GetSecondaryInCurrentStep()
				   ->at(i)
				   ->GetParticleDefinition()
				   ->GetPDGCharge());
	  fEventAction->Addprod_parentType(aTrack->GetParticleDefinition()
					   ->GetPDGEncoding());
	}
      }
    }//if(fEventAction->isAllProductions())
}
