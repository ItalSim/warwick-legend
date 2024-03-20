#include "WLGDRunAction.hh"
#include "WLGDEventAction.hh"
#include "g4root.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include <fstream>
#include <iostream>
using namespace std;

WLGDRunAction::WLGDRunAction(WLGDEventAction* eventAction, G4String name)
: G4UserRunAction()
, fEventAction(eventAction)
, fout(std::move(name))
{
  DefineCommands();

  // Create analysis manager
  auto analysisManager = G4AnalysisManager::Instance();

  // Create directories
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetNtupleMerging(true);
  //Warning: SetNtupleMerging can cause issues with outputs that are not thread-safe!

  // -- write out some branches only if the information is also stored


  // Creating ntuple with vector entries
  //
  analysisManager->CreateNtuple("Score", "Hits");//Index=0
  /*analysisManager->CreateNtupleIColumn("NGe77", fEventAction->GetNGe77());
  analysisManager->CreateNtupleIColumn("HitID", fEventAction->GetHitTID());
  analysisManager->CreateNtupleDColumn("Edep", fEventAction->GetHitEdep());
  analysisManager->CreateNtupleDColumn("Time", fEventAction->GetHitTime());
  analysisManager->CreateNtupleDColumn("Weight", fEventAction->GetHitWeight());
  analysisManager->CreateNtupleDColumn("Hitxloc", fEventAction->GetHitxLoc());
  analysisManager->CreateNtupleDColumn("Hityloc", fEventAction->GetHityLoc());
  analysisManager->CreateNtupleDColumn("Hitzloc", fEventAction->GetHitzLoc());

  analysisManager->CreateNtupleDColumn("Muonxloc", fEventAction->GetMuonxLoc());
  analysisManager->CreateNtupleDColumn("Muonyloc", fEventAction->GetMuonyLoc());
  analysisManager->CreateNtupleDColumn("Muonzloc", fEventAction->GetMuonzLoc());
  analysisManager->CreateNtupleDColumn("Muonxmom", fEventAction->GetMuonxMom());
  analysisManager->CreateNtupleDColumn("Muonymom", fEventAction->GetMuonyMom());
  analysisManager->CreateNtupleDColumn("Muonzmom", fEventAction->GetMuonzMom());
  analysisManager->CreateNtupleDColumn("MuonEnergy", fEventAction->GetMuonEnergy());

  // Edit: 2021/04/07 by Moritz Neuberger
  // Adding additional outputs to further investigate situations in which Ge-77 is
  // produced
  analysisManager->CreateNtupleDColumn("Ekin", fEventAction->GetHitEkin());
  analysisManager->CreateNtupleIColumn("ReentranceTube",
                                       fEventAction->GetReentranceTube());
  analysisManager->CreateNtupleIColumn("DetectorNumber",
                                       fEventAction->GetDetectorNumber());

  analysisManager->CreateNtupleDColumn("EdepWater_prompt",
                                       fEventAction->GetEdepWater_prompt());
  analysisManager->CreateNtupleDColumn("EdepWater_delayed",
                                       fEventAction->GetEdepWater_delayed());
  analysisManager->CreateNtupleIColumn("MunoVeto", fEventAction->GetMuonVeto_flag());
  analysisManager->CreateNtupleIColumn("isIC", fEventAction->GetisIC());
  analysisManager->CreateNtupleIColumn("isMetastable", fEventAction->GetisMetastable());
  */
  analysisManager->CreateNtupleDColumn("Neutronxloc", fEventAction->GetNeutronxLoc());
  analysisManager->CreateNtupleDColumn("Neutronyloc", fEventAction->GetNeutronyLoc());
  analysisManager->CreateNtupleDColumn("Neutronzloc", fEventAction->GetNeutronzLoc());
  analysisManager->CreateNtupleDColumn("Neutronxmom", fEventAction->GetNeutronxMom());
  analysisManager->CreateNtupleDColumn("Neutronymom", fEventAction->GetNeutronyMom());
  analysisManager->CreateNtupleDColumn("Neutronzmom", fEventAction->GetNeutronzMom());
  analysisManager->CreateNtupleDColumn("NeutronTime", fEventAction->GetNeutronTime());
  analysisManager->CreateNtupleDColumn("NeutronEkin", fEventAction->GetNeutronEkin());
  analysisManager->CreateNtupleDColumn("NeutronID", fEventAction->GetNeutronID());
  analysisManager->CreateNtupleDColumn("NeutronEventID", fEventAction->GetNeutronEventID());
  analysisManager->CreateNtupleDColumn("NeutronVolume", fEventAction->GetNeutronVolume());
  analysisManager->CreateNtupleIColumn("NeutronsInEvent",
                                       fEventAction->GetNumberOfNeutronsInEvent());
  analysisManager->CreateNtupleDColumn("NeutronsMostOuterRadius",
                                       fEventAction->GetNeutronsMostOuterRadius());

    //analysisManager->CreateNtupleDColumn("Neutronxtrack",
    //fEventAction->GetNeutronxTrack());
    //analysisManager->CreateNtupleDColumn("Neutronytrack",
    //fEventAction->GetNeutronyTrack());
    //analysisManager->CreateNtupleDColumn("Neutronztrack",
    //fEventAction->GetNeutronzTrack());
  /*
  analysisManager->CreateNtupleDColumn("LArEnergyDeposition",
                                       fEventAction->GetLArEnergyDeposition());
  analysisManager->CreateNtupleDColumn("GeEnergyDeposition",
                                       fEventAction->GetGeEnergyDeposition());
  analysisManager->CreateNtupleDColumn("LArEnergyDeposition_delayed",
                                       fEventAction->GetLArEnergyDeposition_delayed());
  analysisManager->CreateNtupleDColumn("GeEnergyDeposition_delayed",
                                       fEventAction->GetGeEnergyDeposition_delayed());
  if(fWriteOutAdvancedMultiplicity)
    analysisManager->CreateNtupleDColumn(
      "LArEnergyDeposition_delayed_long",
      fEventAction->GetLArEnergyDeposition_delayed_long());
  if(fWriteOutAdvancedMultiplicity)
    analysisManager->CreateNtupleDColumn(
      "GeEnergyDeposition_delayed_long",
      fEventAction->GetGeEnergyDeposition_delayed_long());
  if(fWriteOutAdvancedMultiplicity)
    analysisManager->CreateNtupleDColumn(
      "LArEnergyDeposition_after_delayed",
      fEventAction->GetLArEnergyDeposition_after_delayed());
  if(fWriteOutAdvancedMultiplicity)
    analysisManager->CreateNtupleDColumn(
      "GeEnergyDeposition_after_delayed",
      fEventAction->GetGeEnergyDeposition_after_delayed());

  analysisManager->CreateNtupleDColumn(
    "IndividualEnergyDeposition_Timing",
    fEventAction->GetIndividualEnergyDeposition_Timing());
  analysisManager->CreateNtupleDColumn(
    "IndividualEnergyDeposition_Energy",
    fEventAction->GetIndividualEnergyDeposition_Energy());
  analysisManager->CreateNtupleDColumn(
    "IndividualEnergyDeposition_Position_x",
    fEventAction->GetIndividualEnergyDeposition_Position_x());
  analysisManager->CreateNtupleDColumn(
    "IndividualEnergyDeposition_Position_y",
    fEventAction->GetIndividualEnergyDeposition_Position_y());
  analysisManager->CreateNtupleDColumn(
    "IndividualEnergyDeposition_Position_z",
    fEventAction->GetIndividualEnergyDeposition_Position_z());
  analysisManager->CreateNtupleIColumn(
    "IndividualEnergyDeposition_ReentranceTube",
    fEventAction->GetIndividualEnergyDeposition_ReentranceTube());
  analysisManager->CreateNtupleIColumn(
    "IndividualEnergyDeposition_Volume",
    fEventAction->GetIndividualEnergyDeposition_LArOrGe());
  analysisManager->CreateNtupleIColumn("IndividualEnergyDeposition_ID",
                                       fEventAction->GetIndividualEnergyDeposition_ID());
  analysisManager->CreateNtupleIColumn(
    "IndividualEnergyDeposition_Type",
    fEventAction->GetIndividualEnergyDeposition_Type());
  analysisManager->CreateNtupleIColumn(
    "IndividualEnergyDeposition_DetectorNumber",
    fEventAction->GetIndividualEnergyDeposition_DetectorNumber());

  analysisManager->CreateNtupleIColumn("Multiplicity_prompt",
                                       fEventAction->GetMultiplicity_prompt());
  analysisManager->CreateNtupleIColumn("Multiplicity_delayed",
                                       fEventAction->GetMultiplicity_delayed());
  if(fWriteOutAdvancedMultiplicity)
    analysisManager->CreateNtupleIColumn("Multiplicity_delayed_long",
                                         fEventAction->GetMultiplicity_delayed_long());
  analysisManager->CreateNtupleDColumn("EdepPerDetector_prompt",
                                       fEventAction->GetEdepPerDetector_prompt());
  analysisManager->CreateNtupleDColumn("EdepPerDetector_delayed",
                                       fEventAction->GetEdepPerDetector_delayed());
  if(fWriteOutAdvancedMultiplicity)
    analysisManager->CreateNtupleDColumn("EdepPerDetector_delayed_long",
                                         fEventAction->GetEdepPerDetector_delayed_long());
  analysisManager->CreateNtupleDColumn("CopyNDetector_prompt",
                                       fEventAction->GetNDetector_prompt());
  analysisManager->CreateNtupleDColumn("CopyNDetector_delayed",
                                       fEventAction->GetNDetector_delayed());
  if(fWriteOutAdvancedMultiplicity)
    analysisManager->CreateNtupleDColumn("CopyNDetector_delayed_long",
                                         fEventAction->GetNDetector_delayed_long());

  if(fWriteOutAdvancedMultiplicity)
  {
    //  analysisManager->CreateNtupleIColumn("Multiplicity_prompt_woGd",
    //  fEventAction->GetMultiplicity_prompt_woGd());
    //  analysisManager->CreateNtupleIColumn("Multiplicity_delayed_woGd",
    //  fEventAction->GetMultiplicity_delayed_woGd())
    analysisManager->CreateNtupleDColumn("EdepPerDetector_prompt_woGd",
                                         fEventAction->GetEdepPerDetector_prompt_woGd());
    analysisManager->CreateNtupleDColumn("EdepPerDetector_delayed_woGd",
                                         fEventAction->GetEdepPerDetector_delayed_woGd());
    analysisManager->CreateNtupleDColumn("CopyNDetector_prompt_woGd",
                                         fEventAction->GetNDetector_prompt_woGd());
    analysisManager->CreateNtupleDColumn("CopyNDetector_delayed_woGd",
                                         fEventAction->GetNDetector_delayed_woGd());

    //  analysisManager->CreateNtupleIColumn("Multiplicity_prompt_onlyGd",
    //  fEventAction->GetMultiplicity_prompt_onlyGd());
    //  analysisManager->CreateNtupleIColumn("Multiplicity_delayed_onlyGd",
    //  fEventAction->GetMultiplicity_delayed_onlyGd());
    analysisManager->CreateNtupleDColumn(
      "EdepPerDetector_prompt_onlyGd", fEventAction->GetEdepPerDetector_prompt_onlyGd());
    analysisManager->CreateNtupleDColumn(
      "EdepPerDetector_delayed_onlyGd",
      fEventAction->GetEdepPerDetector_delayed_onlyGd());
    analysisManager->CreateNtupleDColumn("CopyNDetector_prompt_onlyGd",
                                         fEventAction->GetNDetector_prompt_onlyGd());
    analysisManager->CreateNtupleDColumn("CopyNDetector_delayed_onlyGd",
                                         fEventAction->GetNDetector_delayed_onlyGd());
  }

  if(fIndividualGeDepositionInfo)
  {
    analysisManager->CreateNtupleDColumn("Ge77Siblings_timing",
                                         fEventAction->GetGe77Siblings_timing());
    analysisManager->CreateNtupleDColumn("Ge77Siblings_x",
                                         fEventAction->GetGe77Siblings_x());
    analysisManager->CreateNtupleDColumn("Ge77Siblings_y",
                                         fEventAction->GetGe77Siblings_y());
    analysisManager->CreateNtupleDColumn("Ge77Siblings_z",
                                         fEventAction->GetGe77Siblings_z());
    analysisManager->CreateNtupleDColumn("Ge77Siblings_edep",
                                         fEventAction->GetGe77Siblings_edep());
    analysisManager->CreateNtupleIColumn("Ge77Siblings_id",
                                         fEventAction->GetGe77Siblings_id());
    analysisManager->CreateNtupleIColumn("Ge77Siblings_type",
                                         fEventAction->GetGe77Siblings_type());
    analysisManager->CreateNtupleIColumn("Ge77Siblings_whichVolume",
                                         fEventAction->GetGe77Siblings_whichVolume());
  }
  //if(fNeutronCaptureSiblings){
    analysisManager->CreateNtupleIColumn("NeutronCaptureSiblings_type",
                                         fEventAction->GetNeutronCaptureSiblings_type());
    analysisManager->CreateNtupleDColumn("NeutronCaptureSiblings_edep",
                                         fEventAction->GetNeutronCaptureSiblings_edep());
  //y}
  if(fIndividualGdDepositionInfo)
  {
    analysisManager->CreateNtupleDColumn("GdSiblings_timing",
                                         fEventAction->GetGdSiblings_timing());
    analysisManager->CreateNtupleDColumn("GdSiblings_x", fEventAction->GetGdSiblings_x());
    analysisManager->CreateNtupleDColumn("GdSiblings_y", fEventAction->GetGdSiblings_y());
    analysisManager->CreateNtupleDColumn("GdSiblings_z", fEventAction->GetGdSiblings_z());
    analysisManager->CreateNtupleDColumn("GdSiblings_edep",
                                         fEventAction->GetGdSiblings_edep());
    analysisManager->CreateNtupleIColumn("GdSiblings_id",
                                         fEventAction->GetGdSiblings_id());
    analysisManager->CreateNtupleIColumn("GdSiblings_type",
                                         fEventAction->GetGdSiblings_type());
    analysisManager->CreateNtupleIColumn("GdSiblings_whichVolume",
                                         fEventAction->GetGdSiblings_whichVolume());
  }*/

  // analysisManager->CreateNtupleDColumn("nCAr_timing", fEventAction->GetnCAr_timing());
  // analysisManager->CreateNtupleDColumn("nCAr_x", fEventAction->GetnCAr_x());
  // analysisManager->CreateNtupleDColumn("nCAr_y", fEventAction->GetnCAr_y());
  // analysisManager->CreateNtupleDColumn("nCAr_z", fEventAction->GetnCAr_z());
  // analysisManager->CreateNtupleIColumn("nCAr_A", fEventAction->GetnCAr_A());

  // analysisManager->CreateNtupleDColumn("nCGd_timing", fEventAction->GetnCGd_timing());
  // analysisManager->CreateNtupleDColumn("nCGd_x", fEventAction->GetnCGd_x());
  // analysisManager->CreateNtupleDColumn("nCGd_y", fEventAction->GetnCGd_y());
  // analysisManager->CreateNtupleDColumn("nCGd_z", fEventAction->GetnCGd_z());
  // analysisManager->CreateNtupleIColumn("nCGd_A", fEventAction->GetnCGd_A());

  analysisManager->CreateNtupleDColumn("nCapture_timing",
                                       fEventAction->GetnCOther_timing());
  analysisManager->CreateNtupleDColumn("nCapture_x", fEventAction->GetnCOther_x());
  analysisManager->CreateNtupleDColumn("nCapture_y", fEventAction->GetnCOther_y());
  analysisManager->CreateNtupleDColumn("nCapture_z", fEventAction->GetnCOther_z());
  analysisManager->CreateNtupleIColumn("nCapture_A", fEventAction->GetnCOther_A());
  analysisManager->CreateNtupleIColumn("nCapture_Z", fEventAction->GetnCOther_ZC());
  analysisManager->CreateNtupleIColumn("nCapture_ID", fEventAction->GetnCOther_ID());
  analysisManager->CreateNtupleIColumn("nCapture_EventID", fEventAction->GetnCOther_EventID());
  analysisManager->CreateNtupleIColumn("nCapture_Volume", fEventAction->GetnCOther_Volume());

  analysisManager->CreateNtupleDColumn("prod_timing",
                                       fEventAction->Getprod_timing());
  analysisManager->CreateNtupleDColumn("prod_x", fEventAction->Getprod_x());
  analysisManager->CreateNtupleDColumn("prod_y", fEventAction->Getprod_y());
  analysisManager->CreateNtupleDColumn("prod_z", fEventAction->Getprod_z());
  analysisManager->CreateNtupleIColumn("prod_A", fEventAction->Getprod_A());
  analysisManager->CreateNtupleIColumn("prod_Z", fEventAction->Getprod_ZC());//prod_parentType
  analysisManager->CreateNtupleIColumn("prod_parentType", fEventAction->Getprod_parentType());//prod_parentType


  /*  if(fIndividualGeDepositionInfo)
  {
    analysisManager->CreateNtupleDColumn("Ge77mGammaEmission_timing",
                                         fEventAction->GetGe77mGammaEmission_timing());
    analysisManager->CreateNtupleDColumn("Ge77mGammaEmission_x",
                                         fEventAction->GetGe77mGammaEmission_x());
    analysisManager->CreateNtupleDColumn("Ge77mGammaEmission_y",
                                         fEventAction->GetGe77mGammaEmission_y());
    analysisManager->CreateNtupleDColumn("Ge77mGammaEmission_z",
                                         fEventAction->GetGe77mGammaEmission_z());
    analysisManager->CreateNtupleDColumn("Ge77mGammaEmission_edep",
                                         fEventAction->GetGe77mGammaEmission_edep());
    analysisManager->CreateNtupleIColumn("Ge77mGammaEmission_id",
                                         fEventAction->GetGe77mGammaEmission_id());
    analysisManager->CreateNtupleIColumn("Ge77mGammaEmission_type",
                                         fEventAction->GetGe77mGammaEmission_type());
    analysisManager->CreateNtupleIColumn(
      "Ge77mGammaEmission_whichVolume",
      fEventAction->GetGe77mGammaEmission_whichVolume());
    analysisManager->CreateNtupleIColumn("Ge77mGammaEmission_whichGe77",
                                         fEventAction->GetGe77mGammaEmission_whichGe77());
  }
  // GetnCAr_timing
  analysisManager->CreateNtupleIColumn("Trjpdg", fEventAction->GetTrjPDG());
  analysisManager->CreateNtupleIColumn("Trjentries", fEventAction->GetTrjEntries());
  analysisManager->CreateNtupleIColumn("VtxName", fEventAction->GetNameID());
  analysisManager->CreateNtupleDColumn("TrjXVtx", fEventAction->GetTrjXVtx());
  analysisManager->CreateNtupleDColumn("TrjYVtx", fEventAction->GetTrjYVtx());
  analysisManager->CreateNtupleDColumn("TrjZVtx", fEventAction->GetTrjZVtx());
  analysisManager->CreateNtupleDColumn("TrjXPos", fEventAction->GetTrjXPos());
  analysisManager->CreateNtupleDColumn("TrjYPos", fEventAction->GetTrjYPos());
  analysisManager->CreateNtupleDColumn("TrjZPos", fEventAction->GetTrjZPos());

  analysisManager->CreateNtupleDColumn("Muon_WLSR_intersect_x", fEventAction->Get_Muon_WLSR_intersect_x());
  analysisManager->CreateNtupleDColumn("Muon_WLSR_intersect_y", fEventAction->Get_Muon_WLSR_intersect_y());
  analysisManager->CreateNtupleDColumn("Muon_WLSR_intersect_z", fEventAction->Get_Muon_WLSR_intersect_z());
  analysisManager->CreateNtupleDColumn("Muon_WLSR_Edep", fEventAction->Get_Muon_WLSR_Edep());
  */
  analysisManager->FinishNtuple();

  //Even if they remain unused, we must create every ntuple
  //This is because, when using the analysis manager with multiple ntuples, the index of the ntuple is used when filling with data
  //So, if we don't create the ntuple which normally has index=1, the next ntuple (which should be index=2) will have index=1

  analysisManager->CreateNtuple("Steps","Step-level data");//Index=1
  
      analysisManager->CreateNtupleDColumn("X");             //0
      analysisManager->CreateNtupleDColumn("Y");             //1
      analysisManager->CreateNtupleDColumn("Z");             //2
      analysisManager->CreateNtupleDColumn("Time");          //3
      analysisManager->CreateNtupleDColumn("KineticEnergy"); //4
      analysisManager->CreateNtupleIColumn("TrackID");       //5
      analysisManager->CreateNtupleIColumn("ID");            //6
      analysisManager->CreateNtupleIColumn("EventID");       //7
      analysisManager->CreateNtupleIColumn("PID");           //8
      analysisManager->CreateNtupleSColumn("Process");       //9
      analysisManager->CreateNtupleSColumn("CreatorProcess");//10
      analysisManager->CreateNtupleSColumn("Material");      //11
      analysisManager->CreateNtupleSColumn("Volume");        //12
      analysisManager->CreateNtupleDColumn("EDep");          //13
  
      //analysisManager->CreateNtupleDColumn("X");             //0
      //analysisManager->CreateNtupleDColumn("Y");             //1
      //analysisManager->CreateNtupleDColumn("Z");             //2
      //analysisManager->CreateNtupleDColumn("Time");          //3
      //analysisManager->CreateNtupleIColumn("EventID");       //4
      //analysisManager->CreateNtupleDColumn("EDep");          //5

  analysisManager->FinishNtuple();

  
  analysisManager->CreateNtuple("OpticsTracks","Track-level optical photon data"); //Index=2
  
  analysisManager->CreateNtupleDColumn("X");             //0
  analysisManager->CreateNtupleDColumn("Y");             //1
  analysisManager->CreateNtupleDColumn("Z");             //2
  analysisManager->CreateNtupleDColumn("Time");          //3
  analysisManager->CreateNtupleDColumn("KineticEnergy"); //4
  analysisManager->CreateNtupleDColumn("Wavelength");    //5
  analysisManager->CreateNtupleIColumn("TrackID");       //6
  analysisManager->CreateNtupleIColumn("EventID");       //7
  analysisManager->CreateNtupleSColumn("CreatorProcess");//8
  analysisManager->CreateNtupleSColumn("Volume");        //9
  analysisManager->CreateNtupleDColumn("px");             //10
  analysisManager->CreateNtupleDColumn("py");             //11
  analysisManager->CreateNtupleDColumn("pz");             //12
  
  
  analysisManager->FinishNtuple();

  analysisManager->CreateNtuple("OpticalMapData","Data for creating the WLGD optical map");//Index=3
  
  analysisManager->CreateNtupleFColumn("X1");             //0
  analysisManager->CreateNtupleFColumn("Y1");             //1
  analysisManager->CreateNtupleFColumn("Z1");             //2
  analysisManager->CreateNtupleFColumn("X2");             //3
  analysisManager->CreateNtupleFColumn("Y2");             //4
  analysisManager->CreateNtupleFColumn("Z2");             //5
  
  analysisManager->FinishNtuple();

}

WLGDRunAction::~WLGDRunAction() { delete G4AnalysisManager::Instance(); }

void WLGDRunAction::BeginOfRunAction(const G4Run* /*run*/)
{
  //Start the clock
  runtimer = double(clock());
  
  // Get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();
  // Open an output file
  //
  analysisManager->OpenFile(fout);

  fNumberOfCrossingNeutrons   = 0;
  fTotalNumberOfNeutronsInLAr = 0;

  ofstream outputStream;
  if(fWriteOutGeneralNeutronInfo == 1)
  {
    fout2 = fout + "_NCrossing.txt";
    outputStream.open(fout2.c_str(), ios::trunc);
    outputStream.close();
  }

  if(fWriteOutNeutronProductionInfo == 1)
  {
    fout3 = fout + "_NCreationPosition.txt";
    outputStream_2.open(fout3.c_str(), ios::trunc);
  }
}

void WLGDRunAction::EndOfRunAction(const G4Run* /*run*/)
{
  // Get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();

  G4cout << "NumberOfNeutronCrossings: " << fNumberOfCrossingNeutrons << G4endl;
  G4cout << "TotalNumberOfNeutronInLAr: " << fTotalNumberOfNeutronsInLAr << G4endl;

  // save ntuple
  analysisManager->Write();
  analysisManager->CloseFile();

  // Edit: 2021/03/12 by Moritz Neuberger
  // Adding output for number of neutrons crossing the detectors and total produced in LAr
  if(fWriteOutGeneralNeutronInfo == 1)
  {
    ofstream outputStream;
    outputStream.open(fout2.c_str(), ios::app);
    outputStream << fNumberOfCrossingNeutrons << "   " << fTotalNumberOfNeutronsInLAr
                 << endl;
    outputStream.close();
  }

  // Edit: 2021/03/12 by Moritz Neuberger
  // Adding detail output for neutron production information
  G4cout << "N of neutrons total: " << vector_x_dir.size() << G4endl;
  if(fWriteOutNeutronProductionInfo == 1)
  {
    for(int i = 0; i < vector_x_dir.size(); i++)
    {
      outputStream_2 << vector_eventNumber[i] << " " << vector_x_dir[i] << " "
                     << vector_y_dir[i] << " " << vector_z_dir[i] << " "
                     << vector_x_mom[i] << " " << vector_y_mom[i] << " "
                     << vector_z_mom[i] << " " << vector_energy[i] << " "
                     << vector_parentParticleType[i] << endl;
    }
    outputStream_2.close();
  }
  runtimer = double(clock()) - runtimer;
  G4cout << G4endl << "Total processing time: " << runtimer/1000000 << " seconds." << G4endl << G4endl;
}

void WLGDRunAction::SetWriteOutNeutronProductionInfo(G4int answer)
{
  fWriteOutNeutronProductionInfo = answer;
}

void WLGDRunAction::SetWriteOutGeneralNeutronInfo(G4int answer)
{
  fWriteOutGeneralNeutronInfo = answer;
}

void WLGDRunAction::SetWriteOutAdvancedMultiplicity(G4int answer)
{
  fWriteOutAdvancedMultiplicity = answer;
}

void WLGDRunAction::SetWriteOutAllNeutronInfoRoot(G4int answer)
{
  fWriteOutAllNeutronInfoRoot = answer;
}

void WLGDRunAction::SetIndividualGeDepositionInfo(G4int answer)
{
  fIndividualGeDepositionInfo = answer;
}

void WLGDRunAction::SetIndividualGdDepositionInfo(G4int answer)
{
  fIndividualGdDepositionInfo = answer;
}

void WLGDRunAction::SetReadMuonCrossingWLSR(G4int answer)
{
  fReadMuonCrossingWLSR = answer;
}

void WLGDRunAction::SetNeutronCaptureSiblings(G4int answer)
{
  fNeutronCaptureSiblings = answer;
}

void WLGDRunAction::SetWriteOutOpticalData(G4int answer)
{
  fWriteOutOpticalData = answer;
}

void WLGDRunAction::SetWriteOutOpticalMapData(G4int answer)
{
  fWriteOutOpticalMapData = answer;
}

void WLGDRunAction::SetReduceStepsData(G4int answer)
{
  fReducedStepsData = answer;
}

void WLGDRunAction::SetWriteOutStepData(G4int answer)
{
  fWriteOutStepData = answer;
}


void WLGDRunAction::DefineCommands()
{
  // Define /WLGD/generator command directory using generic messenger class
  fMessenger = new G4GenericMessenger(this, "/WLGD/runaction/", "Run Action control");

  fMessenger
    ->DeclareMethod("WriteOutNeutronProductionInfo",
                    &WLGDRunAction::SetWriteOutNeutronProductionInfo)
    .SetGuidance("Set whether to write out Neutron Production Info")
    .SetGuidance("0 = without write out")
    .SetGuidance("1 = with write out")
    .SetCandidates("0 1")
    .SetDefaultValue("0");

  fMessenger
    ->DeclareMethod("WriteOutGeneralNeutronInfo",
                    &WLGDRunAction::SetWriteOutGeneralNeutronInfo)
    .SetGuidance("Set whether to write out General Neutron Info")
    .SetGuidance("0 = without write out")
    .SetGuidance("1 = with write out")
    .SetCandidates("0 1")
    .SetDefaultValue("0");

  fMessenger
    ->DeclareMethod("WriteOutAllNeutronInfoRoot",
                    &WLGDRunAction::SetWriteOutAllNeutronInfoRoot)
    .SetGuidance("Set whether to write out all Neutron Info into root")
    .SetGuidance("0 = without write out")
    .SetGuidance("1 = with write out")
    .SetCandidates("0 1")
    .SetDefaultValue("0");

  fMessenger
    ->DeclareMethod("WriteOutAdvancedMultiplicity",
                    &WLGDRunAction::SetWriteOutAdvancedMultiplicity)
    .SetGuidance("Set whether to write out General Neutron Info")
    .SetGuidance("0 = without write out")
    .SetGuidance("1 = with write out")
    .SetCandidates("0 1")
    .SetDefaultValue("0");

  fMessenger
    ->DeclareMethod("getIndividualGeDepositionInfo",
                    &WLGDRunAction::SetIndividualGeDepositionInfo)
    .SetGuidance(
      "Set whether to obtain individual energy deposition information inside Ge")
    .SetGuidance("0 = don't")
    .SetGuidance("1 = do")
    .SetCandidates("0 1")
    .SetDefaultValue("0");

  fMessenger
    ->DeclareMethod("getIndividualGdDepositionInfo",
                    &WLGDRunAction::SetIndividualGdDepositionInfo)
    .SetGuidance(
      "Set whether to obtain individual energy deposition information inside with Gd")
    .SetGuidance("0 = don't")
    .SetGuidance("1 = do")
    .SetCandidates("0 1")
    .SetDefaultValue("0");

  fMessenger
    ->DeclareMethod("getNeutronCaptureSiblings",
                    &WLGDRunAction::SetNeutronCaptureSiblings)
    .SetGuidance(
      "Set whether to obtain the information about the direct neutron capture siblings")
    .SetGuidance("0 = don't")
    .SetGuidance("1 = do")
    .SetCandidates("0 1")
    .SetDefaultValue("0");


  fMessenger
    ->DeclareMethod("readMuonCrossingWLSR",
                    &WLGDRunAction::SetReadMuonCrossingWLSR)
    .SetGuidance(
      "Set whether to write out the intersects of a muon crossing the WLSR plus the total energy deposited in it")
    .SetGuidance("0 = don't")
    .SetGuidance("1 = do")
    .SetCandidates("0 1")
    .SetDefaultValue("0");

  fMessenger
    ->DeclareMethod("WriteOpticalProductionData",
                    &WLGDRunAction::SetWriteOutOpticalData)
    .SetGuidance("Set whether to write out optical photon data at their time of production")
    .SetGuidance("0 = do not output this data")
    .SetGuidance("1 = output this data")
    .SetCandidates("0 1")
    .SetDefaultValue("0");

  fMessenger
    ->DeclareMethod("WriteOpticalMapData",
                    &WLGDRunAction::SetWriteOutOpticalMapData)
    .SetGuidance("Set whether to write out optical map data (advanced option; read the source code for more info)")
    .SetGuidance("0 = do not output this data")
    .SetGuidance("1 = output this data")
    .SetCandidates("0 1")
    .SetDefaultValue("0");

  fMessenger
    ->DeclareMethod("ReduceStepData",
                    &WLGDRunAction::SetReduceStepsData)
    .SetGuidance("Set whether to reduce data stored in the Steps tree (intended for optical studies)")
    .SetGuidance("0 = keep all info in the Steps tree")
    .SetGuidance("1 = only keep x,y,z,t,edep,event#")
    .SetCandidates("0 1")
    .SetDefaultValue("0");

  fMessenger
    ->DeclareMethod("WriteStepData",
                    &WLGDRunAction::SetWriteOutStepData)
    .SetGuidance("Set whether to write out step-level data within the cryostat (set cuts in WLGDSteppingAction.cc)")
    .SetGuidance("0 = do not output this data")
    .SetGuidance("1 = output this data")
    .SetCandidates("0 1")
    .SetDefaultValue("0");

}
