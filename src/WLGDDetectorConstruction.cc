#include "WLGDDetectorConstruction.hh"

#include <cmath>
#include <set>

#include "G4RunManager.hh"

#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4GDMLParser.hh"
#include "G4GeometryManager.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4Tubs.hh"

#include "G4Colour.hh"
#include "G4VisAttributes.hh"

#include "G4SDManager.hh"
#include "WLGDCrystalSD.hh"

#include "WLGDBiasMultiParticleChangeCrossSection.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "G4OpticalSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4MaterialPropertiesTable.hh"


WLGDDetectorConstruction::WLGDDetectorConstruction()
{
  DefineCommands();
  DefineMaterials();
}

WLGDDetectorConstruction::~WLGDDetectorConstruction()
{
  delete fDetectorMessenger;
  delete fBiasMessenger;
}

auto WLGDDetectorConstruction::Construct() -> G4VPhysicalVolume*
{

  // Clean up old geometry
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  if(fGeometryName == "baseline" || fGeometryName == "baseline_smaller" || fGeometryName == "baseline_large_reentrance_tube" || fGeometryName == "baseline_large_reentrance_tube_4m_cryo")
    return SetupBaseline();

  else if(fGeometryName == "hallA" || fGeometryName == "hallA_wo_ge" || fGeometryName == "hallA_only_WLSR")
    return SetupHallA();

  else
    return SetupAlternative();

}//Construct()



void WLGDDetectorConstruction::DefineMaterials()
{
  
  G4NistManager* nistManager = G4NistManager::Instance();
  nistManager->FindOrBuildMaterial("G4_Galactic");
  nistManager->FindOrBuildMaterial("G4_lAr");
  nistManager->FindOrBuildMaterial("G4_AIR");
  nistManager->FindOrBuildMaterial("G4_STAINLESS-STEEL");
  nistManager->FindOrBuildMaterial("G4_Cu");
  nistManager->FindOrBuildMaterial("G4_WATER");

  H    = new G4Element("Hydrogen", "H", 1., 1.00794 * g / mole);
  C    = new G4Element("Carbon", "C", 6., 12.011 * g / mole);
  N    = new G4Element("Nitrogen", "N", 7., 14.00 * g / mole);
  O    = new G4Element("Oxygen", "O", 8., 16.00 * g / mole);
  F    = new G4Element("Fluorine","F",9, 19.00 * g / mole);
  elS  = new G4Element("Sulfur", "S", 16., 32.066 * g / mole);
  Mg   = new G4Element("Magnesium", "Mg", 12., 24.31 * g / mole);
  Ca   = new G4Element("Calcium", "Ca", 20., 40.08 * g / mole);
  elGd = new G4Element("Gadolinium", "Gd", 64, 157.25 * g / mole);
  

  // Standard Rock definition, similar to Gran Sasso rock
  // with density from PDG report
  stdRock = new G4Material("StdRock", 2.65 * g / cm3, 4);
  stdRock->AddElement(O, 52.0 * perCent);
  stdRock->AddElement(Ca, 27.0 * perCent);
  stdRock->AddElement(C, 12.0 * perCent);
  stdRock->AddElement(Mg, 9.0 * perCent);

  puMat = new G4Material("polyurethane", 0.3 * g / cm3, 4);  // high density foam
  puMat->AddElement(H, 16);
  puMat->AddElement(O, 2);
  puMat->AddElement(C, 8);
  puMat->AddElement(N, 2);

  B10           = new G4Isotope("B10", 5, 10, 10 * g / mole);
  B11           = new G4Isotope("B11", 5, 11, 11 * g / mole);
  SpecialB      = new G4Element("SpecialB", "SpecialB", 2);
  G4double B_MassRatio   = 0.2;
  G4double B_NumberRatio = 1 / (1 + (1 - B_MassRatio) / B_MassRatio * 10. / 11.);
  SpecialB->AddIsotope(B10, B_NumberRatio);
  SpecialB->AddIsotope(B11, 1 - B_NumberRatio);

  // Borated Polyethylene according to GEM TN-92-172 TART Calculations of Neutron
  // Attenuation and Neutron-induced Photons on 5 % and 20 % Borated Polyethylene Slabs
  BoratedPET = new G4Material("BoratedPET", 0.95 * g / cm3, 4);  // high density foam
  BoratedPET->AddElement(H, 0.116);
  BoratedPET->AddElement(C, 0.612);
  BoratedPET->AddElement(SpecialB, 0.05);
  BoratedPET->AddElement(O, 0.222);

  // Estimated using the number of elements per chain elements  (C_5 O_2 H_8)
  PMMA = new G4Material("PMMA", 1.18 * g / cm3, 3);
  PMMA->AddElement(H, 0.08);
  PMMA->AddElement(C, 0.60);
  PMMA->AddElement(O, 0.32);


  //Doped PMMA shielding materials - added by CJ in May 2023
  //To calculate change in mass fraction, simply multiply all the current
  //mass fractions by (1 - doping fraction) 

  //Solid boron has multiple allotropes, between 2.35 and 2.52 g/cm^3 in density
  //For simplicity let's use a value of 2.4, it's not super impactful anyway

  //To calculate change in density is similarly easy
  //For PMMA, new density = (1 - doping fraction)*1.18 + (doping fraction) * 2.4
  
  PMMA1percentB = new G4Material("PMMA1percentB", 1.1922 * g / cm3, 4);
  PMMA1percentB->AddElement(H, 0.0792);
  PMMA1percentB->AddElement(C, 0.594 );
  PMMA1percentB->AddElement(O, 0.3168);
  PMMA1percentB->AddElement(SpecialB, 0.01);
  
  PMMA3percentB = new G4Material("PMMA3percentB", 1.2166 * g / cm3, 4);
  PMMA3percentB->AddElement(H, 0.0776);
  PMMA3percentB->AddElement(C, 0.582 );
  PMMA3percentB->AddElement(O, 0.3104);
  PMMA3percentB->AddElement(SpecialB, 0.03);
  
  PMMA5percentB = new G4Material("PMMA5percentB", 1.241 * g / cm3, 4);
  PMMA5percentB->AddElement(H, 0.076);
  PMMA5percentB->AddElement(C, 0.57 );
  PMMA5percentB->AddElement(O, 0.304);
  PMMA5percentB->AddElement(SpecialB, 0.05);
  
  PMMA7percentB = new G4Material("PMMA7percentB", 1.2654 * g / cm3, 4);
  PMMA7percentB->AddElement(H, 0.0744);
  PMMA7percentB->AddElement(C, 0.558 );
  PMMA7percentB->AddElement(O, 0.2976);
  PMMA7percentB->AddElement(SpecialB, 0.07);
  
  PMMA10percentB = new G4Material("PMMA10percentB", 1.302 * g / cm3, 4);
  PMMA10percentB->AddElement(H, 0.072);
  PMMA10percentB->AddElement(C, 0.54);
  PMMA10percentB->AddElement(O, 0.288);
  PMMA10percentB->AddElement(SpecialB, 0.1);


  //Things are just as simple for the Gd-doped PMMA, but adjusting the density is much more significant
  //I'm unsure whether the other study used pure Gd or Gd2O3
  //The first seems more convenient, but the second is more realistic
  //For simplicity and going with my gut feeling, I'll use pure Gd (density 7.9 g / cm3)
  
  PMMA1percentGd = new G4Material("PMMA1percentGd", 1.2472 * g / cm3, 4);
  PMMA1percentGd->AddElement(H, 0.0792);
  PMMA1percentGd->AddElement(C, 0.594 );
  PMMA1percentGd->AddElement(O, 0.3168);
  PMMA1percentGd->AddElement(elGd, 0.01);
  
  PMMA3percentGd = new G4Material("PMMA3percentGd", 1.3816 * g / cm3, 4);
  PMMA3percentGd->AddElement(H, 0.0776);
  PMMA3percentGd->AddElement(C, 0.582 );
  PMMA3percentGd->AddElement(O, 0.3104);
  PMMA3percentGd->AddElement(elGd, 0.03);
  
  PMMA5percentGd = new G4Material("PMMA5percentGd", 1.516 * g / cm3, 4);
  PMMA5percentGd->AddElement(H, 0.076);
  PMMA5percentGd->AddElement(C, 0.57 );
  PMMA5percentGd->AddElement(O, 0.304);
  PMMA5percentGd->AddElement(elGd, 0.05);
  
  PMMA7percentGd = new G4Material("PMMA7percentGd", 1.6504 * g / cm3, 4);
  PMMA7percentGd->AddElement(H, 0.0744);
  PMMA7percentGd->AddElement(C, 0.558 );
  PMMA7percentGd->AddElement(O, 0.2976);
  PMMA7percentGd->AddElement(elGd, 0.07);
  
  PMMA10percentGd = new G4Material("PMMA10percentGd", 1.852 * g / cm3, 4);
  PMMA10percentGd->AddElement(H, 0.072);
  PMMA10percentGd->AddElement(C, 0.54);
  PMMA10percentGd->AddElement(O, 0.288);
  PMMA10percentGd->AddElement(elGd, 0.1);

  
  //The situation is a little more complicated for the Poly-Gd
  //It's easier for us to just make a separate material definition for the Poly-Gd, then add it to the PMMA

  //I assume this is the article which is referenced: https://www.sciencedirect.com/science/article/abs/pii/S1002072117301370
  //I don't have free access to it, however. We'll play it by ear using the definition from another source:
  //https://pubchem.ncbi.nlm.nih.gov/compound/Gadolinium-methacrylate
  //Chemical formula: C12H18GdO6
  //I couldn't find any information on the density of this substance
  //Since this is basically three PMMA polymers with gadolinum jammed between, I'm gonna assume the density is slightly higher
  //Essentially, 3 PMMAs (formula C5H8O2) have a CH3 stripped and replaced with one Gd bond and an H
  //(C5H8O2 - CH2) x 3 = = C12H18O6 + Gd = C12H18GdO6

  //Total molecular weight = 415.5
  //Hydrogen content: (1.0078 x 18)/415.5 = 0.0436
  //Carbon content: (12.011 x 12)/415.5 = 0.3469
  //Oxygen content: (15.999 x 6)/415.5 = .2310
  //Gadolinium content: 157.25/415.5 = 0.3785

  //Density is just a guess
  PolyGd = new G4Material("PolyGd", 1.25 * g / cm3, 4);
  PolyGd->AddElement(H,    0.0436);
  PolyGd->AddElement(C,    0.3469);
  PolyGd->AddElement(O,    0.2310);
  PolyGd->AddElement(elGd, 0.3785);
  
  //Then make a compound material
  PMMA038percentPolyGd = new G4Material("PMMA038percentPolyGd", 1.1803 * g / cm3, 2);
  PMMA038percentPolyGd->AddMaterial(PMMA,   0.9962);
  PMMA038percentPolyGd->AddMaterial(PolyGd, 0.0038);

  PMMA191percentPolyGd = new G4Material("PMMA191percentPolyGd", 1.1813 * g / cm3, 2);
  PMMA191percentPolyGd->AddMaterial(PMMA,   0.9809);
  PMMA191percentPolyGd->AddMaterial(PolyGd, 0.0191);

  PMMA381percentPolyGd = new G4Material("PMMA381percentPolyGd", 1.1827 * g / cm3, 2);
  PMMA381percentPolyGd->AddMaterial(PMMA,   0.9619);
  PMMA381percentPolyGd->AddMaterial(PolyGd, 0.0381);

  
  // Estimated using the number of elements per molecule (C_2 H_4)
  PolyEthylene = new G4Material("PolyEthylene", 0.95 * g / cm3, 2);
  PolyEthylene->AddElement(H, 0.142);
  PolyEthylene->AddElement(C, 0.857);

  G4double density = 3.01 * g / cm3;  // https://www.sigmaaldrich.com/catalog/product/aldrich/203300
                                      // @room temp
  gadoliniumSulfate = new G4Material("GadoliniumSulfate", density, 3);  // Gd2(SO4)3
  gadoliniumSulfate->AddElement(elGd, 2);
  gadoliniumSulfate->AddElement(elS, 3);
  gadoliniumSulfate->AddElement(O, 12);

  purewater = G4Material::GetMaterial(
    "G4_WATER");  // EDIT: changed water to purewater & use it to create "special" water
  water = new G4Material("GdLoadedWater", 1.000000 * g / cm3, 2);
  water->AddMaterial(purewater, 1. - 0.002);
  water->AddMaterial(gadoliniumSulfate, 0.002);

  // enriched Germanium from isotopes
  Ge_74 = new G4Isotope("Ge74", 32, 74, 74.0 * g / mole);
  Ge_76 = new G4Isotope("Ge76", 32, 76, 76.0 * g / mole);

  Ge_70 = new G4Isotope("Ge70",  32, 70, 69.92*g/mole);
  Ge_72 = new G4Isotope("Ge72",  32, 72, 71.92*g/mole);
  Ge_73 = new G4Isotope("Ge73",  32, 73, 73.0*g/mole);

  eGe = new G4Element("enriched Germanium", "enrGe", 2);
  eGe->AddIsotope(Ge_76, 87. * perCent);
  eGe->AddIsotope(Ge_74, 13. * perCent);
  density      = 5.323 * g / cm3;

  if(fMaGeMaterial)
    {
      density = 5.56 * g / cm3;
      eGe = new G4Element("enriched Germanium", "enrGe", 5);  
      eGe->AddIsotope(Ge_70,0.0*perCent);
      eGe->AddIsotope(Ge_72,0.1*perCent);
      eGe->AddIsotope(Ge_73,0.2*perCent);
      eGe->AddIsotope(Ge_74,13.1*perCent);
      eGe->AddIsotope(Ge_76,86.6*perCent);
    }


  roiMat = new G4Material("enrGe", density, 1);
  roiMat->AddElement(eGe, 1);

  
  // Edit: 2020/02/17 by Moritz Neuberger
  // Added new def. of LAr in order to add doping with Xe and He-3

  G4double dLAr  = 1.393 * g / cm3;
  G4double dLXe  = 3.02 * g / cm3;
  G4double dLHe3 = 0.059 * g / cm3;

  G4double fArConc = 1 - (fXeConc + fHe3Conc);

  G4double dComb = 1 / ((fArConc / dLAr) + (fXeConc / dLXe) + (fHe3Conc / dLHe3));

  if(fXeConc || fHe3Conc)//Only output nonzero concentrations
    {
      G4cout << "___________________________________________" << G4endl;
      G4cout << "Mass ratios of cryostat:" << G4endl;
      G4cout << "LAr:   " << fArConc << G4endl;
      G4cout << "LXe:   " << fXeConc << G4endl;
      G4cout << "LHe3:   " << fHe3Conc << G4endl;
      G4cout << "dComb:   " << dComb << G4endl;
      G4cout << "___________________________________________" << G4endl;
    }
  
  //auto* eLAr = new G4Element("LAr", "Ar", 18., 39.95 * g / mole);
  larMat = G4Material::GetMaterial("G4_lAr");
  elXe   = new G4Element("LXe", "Xe", 54., 131.29 * g / mole);
  eHe3   = new G4Element("He3", "He3", 1);
  iHe3   = new G4Isotope("He3", 2, 3);
  eHe3->AddIsotope(iHe3, 1);

  if(fMaGeMaterial)
    {
      density = 1.396  * g / cm3; 
      larEl = new G4Element("LAr", "Ar", 18., 39.95 * g / mole);
      larMat = new G4Material("LAr", density, 1);
      larMat->AddElement(larEl, 1);
    }

  CombinedArXeHe3 = new G4Material("CombinedArXeHe3", dComb, 3, kStateLiquid, 87. * kelvin);
  CombinedArXeHe3->AddMaterial(larMat, fArConc);
  CombinedArXeHe3->AddElement(eHe3, fHe3Conc);
  CombinedArXeHe3->AddElement(elXe, fXeConc);

  //G4cout << CombinedArXeHe3 << G4endl;

  PEN = new G4Material("PEN", 1.35 * g / cm3, 3, kStateSolid);
  PEN->AddElement(C,14);
  PEN->AddElement(H,10);
  PEN->AddElement(O,4);

  polystyrene = new G4Material("polystyrene", 1.050 * g / cm3, 2, kStateSolid);
  polystyrene->AddElement(C,8);
  polystyrene->AddElement(H,8);

  TPB = new G4Material("TPB", 1.08 * g / cm3, 2, kStateSolid);
  TPB->AddElement(C,28);
  TPB->AddElement(H,22);

  tetratex = new G4Material("tetratex", 0.35 * g / cm3, 2, kStateSolid);
  tetratex->AddElement(F,0.76);
  tetratex->AddElement(C,0.24);

  
}//DefineMaterials()



void WLGDDetectorConstruction::ConstructSDandField()
{

  G4SDManager::GetSDMpointer()->SetVerboseLevel(1);

  // Only need to construct the (per-thread) SD once
  if(!fSD.Get())
  {
    G4String       crystalSDname = "CrystalSD";
    WLGDCrystalSD* aCrystalSD =
      new WLGDCrystalSD(crystalSDname, "CrystalHitsCollection", fGeometryName);
    fSD.Put(aCrystalSD);

    // Also only add it once to the SD manager!
    G4SDManager::GetSDMpointer()->AddNewDetector(fSD.Get());
    SetSensitiveDetector("Ge_log", fSD.Get());

    
    // ----------------------------------------------
    // -- operator creation and attachment to volume:
    // ----------------------------------------------
    G4LogicalVolumeStore* volumeStore = G4LogicalVolumeStore::GetInstance();
    /*
        // -- Attach neutron XS biasing to Germanium -> enhance nCapture
        auto *biasnXS = new WLGDBiasMultiParticleChangeCrossSection();
        biasnXS->SetNeutronFactor(fNeutronBias);
        biasnXS->SetMuonFactor(fMuonBias);
        biasnXS->SetNeutronYieldFactor(fNeutronYieldBias);
        G4cout << " >>> Detector: set neutron bias to " << fNeutronBias << G4endl;
        biasnXS->AddParticle("neutron");
        G4LogicalVolume *logicGe = volumeStore->GetVolume("Ge_log");
        biasnXS->AttachTo(logicGe);
        G4LogicalVolume *logicLar = volumeStore->GetVolume("Lar_log");
        biasnXS->AttachTo(logicLar);
        G4LogicalVolume *logicTank = volumeStore->GetVolume("Tank_log");
        biasnXS->AttachTo(logicTank);
        G4LogicalVolume *logicCavern = volumeStore->GetVolume("Cavern_log");
        biasnXS->AttachTo(logicCavern);
        G4LogicalVolume *logicHall = volumeStore->GetVolume("Hall_log");
        biasnXS->AttachTo(logicHall);
    */

    
    // -- Attach neutron XS biasing to Germanium -> enhance nCapture
    auto* biasnXS = new WLGDBiasMultiParticleChangeCrossSection();
    biasnXS->SetNeutronFactor(fNeutronBias);
    biasnXS->SetMuonFactor(fMuonBias);
    G4cout << " >>> Detector: set neutron bias to " << fNeutronBias << G4endl;
    biasnXS->AddParticle("neutron");
    G4LogicalVolume* logicGe = volumeStore->GetVolume("Ge_log");
    biasnXS->AttachTo(logicGe);

    // -- Attach muon XS biasing to all required volumes consistently
    auto* biasmuXS = new WLGDBiasMultiParticleChangeCrossSection();
    biasmuXS->SetNeutronFactor(fNeutronBias);
    biasmuXS->SetMuonFactor(fMuonBias);
    biasmuXS->SetNeutronYieldFactor(fNeutronYieldBias);
    G4cout << " >>> Detector: set muon bias to " << fMuonBias << G4endl;
    biasmuXS->AddParticle("mu-");

    if(fNeutronYieldBias != 1)
      {
	biasmuXS->AddParticle("neutron");
	biasmuXS->AddParticle("pi+");
	biasmuXS->AddParticle("pi-");
	biasmuXS->AddParticle("gamma");
	biasmuXS->AddParticle("kaon-");
	biasmuXS->AddParticle("proton");
      }

    
    // G4LogicalVolume* logicGe = volumeStore->GetVolume("Ge_log");
    // biasmuXS->AttachTo(logicGe);
    G4LogicalVolume* logicCavern = volumeStore->GetVolume("Cavern_log");
    biasmuXS->AttachTo(logicCavern);
    G4LogicalVolume* logicHall = volumeStore->GetVolume("Hall_log");
    biasmuXS->AttachTo(logicHall);
    G4LogicalVolume* logicTank = volumeStore->GetVolume("Tank_log");
    biasmuXS->AttachTo(logicTank);
    G4LogicalVolume* logicLar = volumeStore->GetVolume("Lar_log");
    biasmuXS->AttachTo(logicLar);

    
    // non hallA have these volumes
    if(fGeometryName != "hallA")
      {
	G4LogicalVolume* logicCu = volumeStore->GetVolume("Copper_log");
	biasmuXS->AttachTo(logicCu);
	G4LogicalVolume* logicULar = volumeStore->GetVolume("ULar_log");
	biasmuXS->AttachTo(logicULar);
      }

    
    // Baseline also has a water volume and cryostat
    if(fGeometryName == "baseline" || fGeometryName == "hallA" || fGeometryName == "hallA_wo_ge" || fGeometryName == "hallA_only_WLSR" ||
       fGeometryName == "baseline_smaller" || fGeometryName == "baseline_large_reentrance_tube" || fGeometryName == "baseline_large_reentrance_tube_4m_cryo")
      {
	G4LogicalVolume* logicWater = volumeStore->GetVolume("Water_log");
	biasmuXS->AttachTo(logicWater);
	G4LogicalVolume* logicCout = volumeStore->GetVolume("Cout_log");
	biasmuXS->AttachTo(logicCout);
	G4LogicalVolume* logicCinn = volumeStore->GetVolume("Cinn_log");
	biasmuXS->AttachTo(logicCinn);
	G4LogicalVolume* logicCLid = volumeStore->GetVolume("Lid_log");
	biasmuXS->AttachTo(logicCLid);
	G4LogicalVolume* logicCBot = volumeStore->GetVolume("Bot_log");
	biasmuXS->AttachTo(logicCBot);
      }
    
    // Alternative has the membrane and insulator
    else if(fGeometryName == "alternative")
      {
	G4LogicalVolume* logicPu = volumeStore->GetVolume("Pu_log");
	biasmuXS->AttachTo(logicPu);
	G4LogicalVolume* logicMembrane = volumeStore->GetVolume("Membrane_log");
	biasmuXS->AttachTo(logicMembrane);
      }
  }//if(!fSD.Get())
  
  //else
      //G4cout << " >>> fSD has entry. Repeated call." << G4endl;

}//ConstructSDandField()



auto WLGDDetectorConstruction::SetupAlternative() -> G4VPhysicalVolume*
{
  // Get materials
  worldMaterial = G4Material::GetMaterial("G4_Galactic");
  //auto* larMat        = G4Material::GetMaterial("G4_lAr");
  airMat        = G4Material::GetMaterial("G4_AIR");
  steelMat      = G4Material::GetMaterial("G4_STAINLESS-STEEL");
  copperMat     = G4Material::GetMaterial("G4_Cu");
  stdRock       = G4Material::GetMaterial("StdRock");
  puMat         = G4Material::GetMaterial("polyurethane");
  roiMat        = G4Material::GetMaterial("enrGe");
  larMat_alt    = G4Material::GetMaterial("CombinedArXeHe3");

  if(fXeConc != 0 || fHe3Conc != 0)
    larMat = larMat_alt;


  // size parameter, unit [cm]

  // cavern
  G4double stone     = 100.0;  // Hall wall thickness 1 m
  G4double hallhside = 850.0;  // Hall cube side 17 m

  // cryostat
  G4double tankhside  = 650;   // cryostat cube side 13 m
  G4double outerwall  = 1.2;   // outer SS wall thickness
  G4double insulation = 80.0;  // polyurethane foam
  G4double innerwall  = 0.12;  // inner SS membrane

  // copper tubes with Germanium ROI
  G4double copper    = 0.35;    // tube thickness 3.5 mm
  G4double curad     = 40.0;    // copper tube diam 80 cm
  G4double cuhheight = 334.34;  // copper tube height 7 m inside cryostat
  G4double cushift   = 234.34;  // shift cu tube inside cryostat to top
  G4double ringrad   = 100.0;   // cu tube placement ring radius

  // Ge cylinder for 250 kg at 5.32 g/cm3
  G4double roiradius = 30.0;  // string radius curad - Ge radius - gap

  G4double gerad          = 4.0;                    // Ge radius
  G4double gehheight      = 5.0;                    // full height 10 cm
  G4double gegap          = 3.0;                    // gap between Ge 3cm
  G4double layerthickness = gegap + 2 * gehheight;  // 13 cm total

  G4int    nofLayers      = 8;   // 8 Ge + 7 gaps = 1010 mm string height
  G4int    nofStrings     = 12;  // 12 strings  of 8 Ge each

  // total
  G4double offset    = hallhside - tankhside;  // shift cavern floor to keep detector centre at origin
  G4double worldside = hallhside + stone + offset + 0.1;  // larger than rest
  G4double larside   = tankhside - outerwall - insulation - innerwall;  // cube side of LAr volume

  fvertexZ = (worldside - stone - 0.1) * cm;  // max vertex height
  fmaxrad  = hallhside * cm;                  // max vertex circle radius


  // Volumes for this geometry

  //
  // World
  //
  auto* worldSolid = new G4Box("World", worldside * cm, worldside * cm, worldside * cm);
  auto* fWorldLogical  = new G4LogicalVolume(worldSolid, worldMaterial, "World_log");
  auto* fWorldPhysical = new G4PVPlacement(nullptr, G4ThreeVector(), fWorldLogical,
                                           "World_phys", nullptr, false, 0);

  //
  // Cavern
  //
  auto* cavernSolid    = new G4Box("Cavern", (hallhside + stone) * cm,
                                   (hallhside + stone) * cm, (hallhside + stone) * cm);
  auto* fCavernLogical = new G4LogicalVolume(cavernSolid, stdRock, "Cavern_log");
  auto* fCavernPhysical =
    new G4PVPlacement(nullptr, G4ThreeVector(0., 0., offset * cm), fCavernLogical,
                      "Cavern_phys", fWorldLogical, false, 0);

  //
  // Hall
  //
  auto* hallSolid = new G4Box("Cavern", hallhside * cm, hallhside * cm, hallhside * cm);
  auto* fHallLogical = new G4LogicalVolume(hallSolid, airMat, "Hall_log");
  auto* fHallPhysical =
    new G4PVPlacement(nullptr, G4ThreeVector(0., 0., -stone * cm), fHallLogical,
                      "Hall_phys", fCavernLogical, false, 0, true);

  //
  // Tank
  //
  auto* tankSolid    = new G4Box("Tank", tankhside * cm, tankhside * cm, tankhside * cm);
  auto* fTankLogical = new G4LogicalVolume(tankSolid, steelMat, "Tank_log");
  auto* fTankPhysical =
    new G4PVPlacement(nullptr, G4ThreeVector(0., 0., -stone * cm), fTankLogical,
                      "Tank_phys", fHallLogical, false, 0, true);

  //
  // Insulator
  //
  auto* puSolid     = new G4Box("Insulator", (tankhside - outerwall) * cm,
                                (tankhside - outerwall) * cm, (tankhside - outerwall) * cm);
  auto* fPuLogical  = new G4LogicalVolume(puSolid, puMat, "Pu_log");
  auto* fPuPhysical = new G4PVPlacement(nullptr, G4ThreeVector(), fPuLogical, "Pu_phys",
                                        fTankLogical, false, 0, true);

  //
  // Membrane
  //
  auto* membraneSolid = new G4Box("Membrane", (tankhside - outerwall - insulation) * cm,
                                  (tankhside - outerwall - insulation) * cm,
                                  (tankhside - outerwall - insulation) * cm);
  auto* fMembraneLogical = new G4LogicalVolume(membraneSolid, steelMat, "Membrane_log");
  auto* fMembranePhysical =
    new G4PVPlacement(nullptr, G4ThreeVector(), fMembraneLogical, "Membrane_phys",
                      fPuLogical, false, 0, true);

  //
  // LAr
  //
  auto* larSolid     = new G4Box("LAr", larside * cm, larside * cm, larside * cm);
  auto* fLarLogical  = new G4LogicalVolume(larSolid, larMat, "Lar_log");
  auto* fLarPhysical = new G4PVPlacement(nullptr, G4ThreeVector(), fLarLogical,
                                         "Lar_phys", fMembraneLogical, false, 0, true);

  //
  // copper tubes, hollow cylinder shell
  //
  auto* copperSolid = new G4Tubs("Copper", (curad - copper) * cm, curad * cm,
                                 cuhheight * cm, 0.0, CLHEP::twopi);

  //
  // ULAr bath, solid cylinder
  //
  auto* ularSolid = new G4Tubs("ULar", 0.0 * cm, (curad - copper) * cm, cuhheight * cm,
                               0.0, CLHEP::twopi);

  // tower; logical volumes
  auto* fCopperLogical = new G4LogicalVolume(copperSolid, copperMat, "Copper_log");
  auto* fUlarLogical   = new G4LogicalVolume(ularSolid, larMat, "ULar_log");

  //
  // Germanium, solid cylinder
  //
  // layers in tower
  auto* layerSolid = new G4Tubs("LayerSolid", 0.0 * cm, gerad * cm,
                                (gehheight + gegap / 2.0) * cm, 0.0, CLHEP::twopi);

  auto* fLayerLogical = new G4LogicalVolume(layerSolid, larMat, "Layer_log");

  
  // fill one layer
  auto* geSolid =
    new G4Tubs("ROI", 0.0 * cm, gerad * cm, gehheight * cm, 0.0, CLHEP::twopi);

  auto* fGeLogical = new G4LogicalVolume(geSolid, roiMat, "Ge_log");
  new G4PVPlacement(nullptr, G4ThreeVector(0.0, 0.0, -gegap / 2.0 * cm), fGeLogical,
                    "Ge_phys", fLayerLogical, false, 0, true);

  auto* gapSolid =
    new G4Tubs("Gap", 0.0 * cm, gerad * cm, gegap / 2.0 * cm, 0.0, CLHEP::twopi);

  auto* fGapLogical = new G4LogicalVolume(gapSolid, larMat, "Gap_log");
  new G4PVPlacement(nullptr, G4ThreeVector(0.0, 0.0, gehheight * cm), fGapLogical,
                    "Gap_phys", fLayerLogical, false, 0, true);

  
  // place layers as mother volume with unique copy number
  G4double step = (gehheight + gegap / 2) * cm;
  G4double xpos;
  G4double ypos;
  G4double angle = CLHEP::twopi / nofStrings;

  // layer logical into ULarlogical
  for(G4int j = 0; j < nofStrings; j++)
  {
    xpos = roiradius * cm * std::cos(j * angle);
    ypos = roiradius * cm * std::sin(j * angle);
    for(G4int i = 0; i < nofLayers; i++)
    {
      new G4PVPlacement(
        nullptr,
        G4ThreeVector(xpos, ypos,
                      -step + (nofLayers / 2 * layerthickness - i * layerthickness) * cm),
        fLayerLogical, "Layer_phys", fUlarLogical, false, i + j * nofLayers, true);
    }
  }

  
  // placements
  new G4PVPlacement(nullptr, G4ThreeVector(ringrad * cm, 0., cushift * cm),
                    fCopperLogical, "Copper_phys", fLarLogical, false, 0, true);

  new G4PVPlacement(nullptr, G4ThreeVector(ringrad * cm, 0., cushift * cm), fUlarLogical,
                    "ULar_phys", fLarLogical, false, 0, true);

  // tower 2
  new G4PVPlacement(nullptr, G4ThreeVector(0., ringrad * cm, cushift * cm),
                    fCopperLogical, "Copper_phys2", fLarLogical, false, 1, true);

  new G4PVPlacement(nullptr, G4ThreeVector(0., ringrad * cm, cushift * cm), fUlarLogical,
                    "ULar_phys2", fLarLogical, false, 1, true);

  // tower 3
  new G4PVPlacement(nullptr, G4ThreeVector(-ringrad * cm, 0., cushift * cm),
                    fCopperLogical, "Copper_phys3", fLarLogical, false, 2, true);

  new G4PVPlacement(nullptr, G4ThreeVector(-ringrad * cm, 0., cushift * cm), fUlarLogical,
                    "ULar_phys3", fLarLogical, false, 2, true);

  // tower 4
  new G4PVPlacement(nullptr, G4ThreeVector(0., -ringrad * cm, cushift * cm),
                    fCopperLogical, "Copper_phys4", fLarLogical, false, 3, true);

  new G4PVPlacement(nullptr, G4ThreeVector(0., -ringrad * cm, cushift * cm), fUlarLogical,
                    "ULar_phys4", fLarLogical, false, 3, true);

  //
  // Visualization attributes
  //
  fWorldLogical->SetVisAttributes(G4VisAttributes::GetInvisible());

  auto* redVisAtt = new G4VisAttributes(G4Colour::Red());
  redVisAtt->SetVisibility(true);
  auto* greyVisAtt = new G4VisAttributes(G4Colour::Grey());
  greyVisAtt->SetVisibility(true);
  auto* greenVisAtt = new G4VisAttributes(G4Colour::Green());
  greenVisAtt->SetVisibility(true);
  auto* blueVisAtt = new G4VisAttributes(G4Colour::Blue());
  blueVisAtt->SetVisibility(true);

  fCavernLogical->SetVisAttributes(redVisAtt);
  fHallLogical->SetVisAttributes(greyVisAtt);
  fTankLogical->SetVisAttributes(blueVisAtt);
  fPuLogical->SetVisAttributes(greyVisAtt);
  fMembraneLogical->SetVisAttributes(blueVisAtt);
  fLarLogical->SetVisAttributes(greyVisAtt);
  fCopperLogical->SetVisAttributes(greenVisAtt);
  fUlarLogical->SetVisAttributes(greyVisAtt);
  fGeLogical->SetVisAttributes(redVisAtt);

  
  return fWorldPhysical;

}//SetupAlternative()



auto WLGDDetectorConstruction::SetupBaseline() -> G4VPhysicalVolume*
{

  // Get materials
  worldMaterial = G4Material::GetMaterial("G4_Galactic");
  // auto* larMat        = G4Material::GetMaterial("G4_lAr");
  airMat   = G4Material::GetMaterial("G4_AIR");
  waterMat = G4Material::GetMaterial("G4_WATER");
  if(fWithGdWater == 1)
    waterMat = water;
  if(fWithWoWater == 1)
    waterMat = airMat;
  steelMat      = G4Material::GetMaterial("G4_STAINLESS-STEEL");
  copperMat     = G4Material::GetMaterial("G4_Cu");
  stdRock       = G4Material::GetMaterial("StdRock");
  roiMat        = G4Material::GetMaterial("enrGe");
  BoratedPETMat = G4Material::GetMaterial("BoratedPET");
  if(fSetMaterial == "PolyEthylene")
    BoratedPETMat = G4Material::GetMaterial("PolyEthylene");
  if(fSetMaterial == "PMMA")
    BoratedPETMat = G4Material::GetMaterial("PMMA");

  if(fSetMaterial == "PMMA1percentB")
    BoratedPETMat = G4Material::GetMaterial("PMMA1percentB");
  if(fSetMaterial == "PMMA3percentB")
    BoratedPETMat = G4Material::GetMaterial("PMMA3percentB");
  if(fSetMaterial == "PMMA5percentB")
    BoratedPETMat = G4Material::GetMaterial("PMMA5percentB");
  if(fSetMaterial == "PMMA7percentB")
    BoratedPETMat = G4Material::GetMaterial("PMMA7percentB");
  if(fSetMaterial == "PMMA10percentB")
    BoratedPETMat = G4Material::GetMaterial("PMMA10percentB");

  if(fSetMaterial == "PMMA1percentGd")
    BoratedPETMat = G4Material::GetMaterial("PMMA1percentGd");
  if(fSetMaterial == "PMMA3percentGd")
    BoratedPETMat = G4Material::GetMaterial("PMMA3percentGd");
  if(fSetMaterial == "PMMA5percentGd")
    BoratedPETMat = G4Material::GetMaterial("PMMA5percentGd");
  if(fSetMaterial == "PMMA7percentGd")
    BoratedPETMat = G4Material::GetMaterial("PMMA7percentGd");
  if(fSetMaterial == "PMMA10percentGd")
    BoratedPETMat = G4Material::GetMaterial("PMMA10percentGd");

  if(fSetMaterial == "PMMA038percentPolyGd")
    BoratedPETMat = G4Material::GetMaterial("PMMA038percentPolyGd");
  if(fSetMaterial == "PMMA191percentPolyGd")
    BoratedPETMat = G4Material::GetMaterial("PMMA191percentPolyGd");
  if(fSetMaterial == "PMMA381percentPolyGd")
    BoratedPETMat = G4Material::GetMaterial("PMMA381percentPolyGd");


  //auto* larMat /*_alt*/ = G4Material::GetMaterial("CombinedArXeHe3");
  // if(fXeConc != 0 || fHe3Conc != 0) BoratedPET
  
  if(!((fXeConc + fHe3Conc) == 0))
    larMat = CombinedArXeHe3;
  //larMat = G4Material::GetMaterial("G4_lAr");

  // G4cout << larMat << G4endl;

  
  // Edit: 2021/03/30 by Moritz Neuberger
  // Adjusted size of stone (1m->5m) s.t. MUSUN cuboid lies inside.
  // Also adjusted relative geometry relations s.t. they are independent of the stone
  // size.

  // constants
  // size parameter, unit [cm]
  G4double offset = 200.0;  // shift cavern floor to keep detector centre at origin
  G4double offset_2 =
    100.0;  // shift s.t. cavern, hall and tank are in line for different stone sizes
  G4double offset_3 = 100;  // 69.5; // shift to get to the baseline lowered position
  if(fDetectorPosition == "original")
    offset_3 = 0;
  if(fGeometryName == "baseline_smaller"  || fGeometryName == "baseline_large_reentrance_tube_4m_cryo")
    offset_3 = 37.5;

  // cavern
  G4double stone       = 500.0;  // Hall wall thickness 5 m
  G4double hallrad     = 600.0;  // Hall cylinder diam 12 m
  G4double hallhheight = 850.0;  // Hall cylinder height 17 m
  // water tank
  G4double tankwalltop = 0.6;  // water tank thickness at top 6 mm
  G4double tankwallbot = 0.8;  // water tank thickness at bottom 8 mm
  G4double tankrad     = 550;  // water tank diam 11 m
  G4double tankhheight = 650;  // water tank height 13 m
  // cryostat
  G4double cryowall   = 3.0;                   // cryostat wall thickness from GERDA
  G4double vacgap     = 1.0;                   // vacuum gap between walls
  G4double cryrad     = fCryostatOuterRadius;  // 350.0;  // cryostat diam 7 m
  G4double cryhheight = fCryostatHeight;       // 350.0;  // cryostat height 7 m

  if(fGeometryName == "baseline_large_reentrance_tube")
    vacgap     = 50.0;  

  if(fGeometryName == "baseline_smaller" || fGeometryName == "baseline_large_reentrance_tube_4m_cryo")
    {
      cryrad     = 200.0;  // cryostat diam 4 m
      cryhheight = 225.0;  // cryostat height 4.5 m
    }
  
  // Borated PET tubes around copper tubes
  G4double BoratedPETouterrad = 5.0;  // tube thickness 5 cm
  // copper tubes with Germanium ROI
  G4double copper    = 0.35;  // tube thickness 3.5 mm
  G4double curad     = 40.0;  // copper tube diam 80 cm
  G4double cuhheight = (400 - (350 - fCryostatHeight)) / 2.;  // 200.0;  // copper tube height 4 m inside cryostat
  G4double cushift = fCryostatHeight - cuhheight;  // 150.0;  // shift cu tube inside cryostat to top
  if(fGeometryName == "baseline_smaller" || fGeometryName == "baseline_large_reentrance_tube_4m_cryo")
    {
      cuhheight = 137.5;  // cupper height 2.25 m
      cushift   = 87.5;   // shift
    }
  if(fGeometryName == "baseline_large_reentrance_tube" || fGeometryName == "baseline_large_reentrance_tube_4m_cryo")
      curad     = 95.0;  


  G4double ringrad = 100.0;  // cu tube placement ring radius
  // Ge cylinder 2.67 kg at 5.32 g/cm3
  G4double roiradius = 30.0;  // string radius curad - Ge radius - gap
  // total mass 1026.86 kg in 4 towers, each with 8 Ge stacked in 12 strings
  G4double gerad          = 4.0;                    // Ge radius
  G4double gehheight      = 5.0;                    // full height 10 cm
  G4double gegap          = 3.0;                    // gap between Ge 3cm
  G4double layerthickness = gegap + 2 * gehheight;  // 13 cm total
  //    G4int nofLayers = 8;   // 8 Ge + 7 gaps = 1010 mm string height
  //    G4int nofStrings = 12;  // 12 strings  of 8 Ge each
  G4int nofLayers  = 7;   // 8 Ge + 7 gaps = 1010 mm string height
  G4int nofStrings = 14;  // 12 strings  of 8 Ge each

  fvertexZ = (hallhheight + offset) * cm;
  fmaxrad  = hallrad * cm;

  // Volumes for this geometry

  //
  // World
  //
  auto* worldSolid =
    new G4Tubs("World", 0.0 * cm, (hallrad + stone + 0.1) * cm,
               (hallhheight + stone + offset + 0.1) * cm, 0.0, CLHEP::twopi);
  auto* fWorldLogical  = new G4LogicalVolume(worldSolid, worldMaterial, "World_log");
  auto* fWorldPhysical = new G4PVPlacement(nullptr, G4ThreeVector(), fWorldLogical,
                                           "World_phys", nullptr, false, 0);

  //
  // Cavern
  //
  auto* cavernSolid    = new G4Tubs("Cavern", 0.0 * cm, (hallrad + stone) * cm,
                                    (hallhheight + stone) * cm, 0.0, CLHEP::twopi);
  auto* fCavernLogical = new G4LogicalVolume(cavernSolid, stdRock, "Cavern_log");
  auto* fCavernPhysical =
    new G4PVPlacement(nullptr, G4ThreeVector(0., 0., offset * cm), fCavernLogical,
                      "Cavern_phys", fWorldLogical, false, 0);

  //
  // Hall
  //
  auto* hallSolid =
    new G4Tubs("Hall", 0.0 * cm, hallrad * cm, hallhheight * cm, 0.0, CLHEP::twopi);
  auto* fHallLogical = new G4LogicalVolume(hallSolid, airMat, "Hall_log");
  auto* fHallPhysical =
    new G4PVPlacement(nullptr, G4ThreeVector(0., 0., /*-stone*/ -offset_2 * cm),
                      fHallLogical, "Hall_phys", fCavernLogical, false, 0, true);

  //
  // Tank
  //
  auto* tankSolid =
    new G4Cons("Tank", 0.0 * cm, (tankrad + tankwallbot) * cm, 0.0 * cm,
               (tankrad + tankwalltop) * cm, tankhheight * cm, 0.0, CLHEP::twopi);
  auto* fTankLogical  = new G4LogicalVolume(tankSolid, steelMat, "Tank_log");
  auto* fTankPhysical = new G4PVPlacement(
    nullptr,
    G4ThreeVector(0., 0., /*-stone -offset_2*/ -(hallhheight - tankhheight) * cm),
    fTankLogical, "Tank_phys", fHallLogical, false, 0, true);

  //
  // Water
  //
  auto* waterSolid     = new G4Tubs("Water", 0.0 * cm, tankrad * cm,
                                    (tankhheight - tankwallbot) * cm, 0.0, CLHEP::twopi);
  auto* fWaterLogical  = new G4LogicalVolume(waterSolid, waterMat, "Water_log");
  auto* fWaterPhysical = new G4PVPlacement(nullptr, G4ThreeVector(), fWaterLogical,
                                           "Water_phys", fTankLogical, false, 0, true);

  //
  // outer cryostat
  //
  auto* coutSolid =
    new G4Tubs("Cout", 0.0 * cm, cryrad * cm, cryhheight * cm, 0.0, CLHEP::twopi);
  auto* fCoutLogical  = new G4LogicalVolume(coutSolid, steelMat, "Cout_log");
  auto* fCoutPhysical = new G4PVPlacement(nullptr, G4ThreeVector(), fCoutLogical,
                                          "Cout_phys", fWaterLogical, false, 0, true);

  //
  // vacuum gap
  //
  auto* cvacSolid     = new G4Tubs("Cvac", 0.0 * cm, (cryrad - cryowall) * cm,
                                   (cryhheight - cryowall) * cm, 0.0, CLHEP::twopi);
  auto* fCvacLogical  = new G4LogicalVolume(cvacSolid, worldMaterial, "Cvac_log");
  auto* fCvacPhysical = new G4PVPlacement(nullptr, G4ThreeVector(), fCvacLogical,
                                          "Cvac_phys", fCoutLogical, false, 0, true);

  //
  // inner cryostat
  //
  auto* cinnSolid     = new G4Tubs("Cinn", 0.0 * cm, (cryrad - cryowall - vacgap) * cm,
                                   (cryhheight - cryowall - vacgap) * cm, 0.0, CLHEP::twopi);
  auto* fCinnLogical  = new G4LogicalVolume(cinnSolid, steelMat, "Cinn_log");
  auto* fCinnPhysical = new G4PVPlacement(nullptr, G4ThreeVector(), fCinnLogical,
                                          "Cinn_phys", fCvacLogical, false, 0, true);

  //
  // LAr bath
  //
  auto* larSolid     = new G4Tubs("LAr", 0.0 * cm, (cryrad - 2 * cryowall - vacgap) * cm,
                                  (cryhheight - 2 * cryowall - vacgap) * cm, 0.0, CLHEP::twopi);
  auto* fLarLogical  = new G4LogicalVolume(larSolid, larMat, "Lar_log");
  auto* fLarPhysical = new G4PVPlacement(nullptr, G4ThreeVector(), fLarLogical,
                                         "Lar_phys", fCinnLogical, false, 0, true);

  //
  // cryostat Lid
  //
  auto* lidSolid =
    new G4Tubs("Lid", 0.0 * cm, cryrad * cm, cryowall / 2.0 * cm, 0.0, CLHEP::twopi);
  auto* fLidLogical = new G4LogicalVolume(lidSolid, steelMat, "Lid_log");
  auto* fLidPhysical =
    new G4PVPlacement(nullptr, G4ThreeVector(0., 0., (cryhheight + cryowall / 2.0) * cm),
                      fLidLogical, "Lid_phys", fWaterLogical, false, 0, true);
  auto* fBotLogical = new G4LogicalVolume(lidSolid, steelMat, "Bot_log");
  auto* fBotPhysical =
    new G4PVPlacement(nullptr, G4ThreeVector(0., 0., -(cryhheight + cryowall / 2.0) * cm),
                      fBotLogical, "Bot_phys", fWaterLogical, false, 0, true);

  //
  // copper tubes, hollow cylinder shell
  //

  // Box variables
  G4double b_width  = fBoratedTurbineWidth / 2. * cm;   // 2.5 * cm;
  G4double b_length = fBoratedTurbineLength / 2. * cm;  // 0.25 * m;
  G4double b_height = fBoratedTurbineHeight / 2. * cm;  // fCryostatHeight * cm - 0.5 * m;

 /* auto* boratedPETSolid_Tube =
    new G4Tubs("BoratedPET", curad * cm, (BoratedPETouterrad + curad) * cm,
               cuhheight * cm, 0.0, CLHEP::twopi);
 */

  auto* boratedPETSolid_Tube =
    new G4Tubs("BoratedPET", curad * cm, (2*b_width + curad) * cm,
               cuhheight * cm, 0.0, CLHEP::twopi);

  auto* boratedPETSolid_Box = new G4Box("BoratedPET", b_length, b_width, b_height);

  //
  // copper tubes, hollow cylinder shell
  auto* copperSolid = new G4Tubs("Copper", (curad - copper) * cm, curad * cm,
                                 cuhheight * cm, 0.0, CLHEP::twopi);
  //
  // ULAr bath, solid cylinder
  auto* ularSolid = new G4Tubs("ULar", 0.0 * cm, (curad - copper) * cm,
                               (cuhheight - copper) * cm, 0.0, CLHEP::twopi);

  // tower; logical volumes
  auto* fBoratedPETLogical_Tube =
    new G4LogicalVolume(boratedPETSolid_Tube, BoratedPETMat, "BoratedPET_Logical");
  auto* fBoratedPETLogical_Box =
    new G4LogicalVolume(boratedPETSolid_Box, BoratedPETMat, "BoratedPET_Logical");
  auto* fCopperLogical = new G4LogicalVolume(copperSolid, copperMat, "Copper_log");
  auto* fUlarLogical   = new G4LogicalVolume(ularSolid, larMat, "ULar_log");

  //
  // Germanium, solid cylinder
  //
  // layers in tower
  auto* layerSolid = new G4Tubs("LayerSolid", 0.0 * cm, gerad * cm,
                                (gehheight + gegap / 2.0) * cm, 0.0, CLHEP::twopi);

  auto* fLayerLogical = new G4LogicalVolume(layerSolid, larMat, "Layer_log");

  // fill one layer
  auto* geSolid =
    new G4Tubs("ROI", 0.0 * cm, gerad * cm, gehheight * cm, 0.0, CLHEP::twopi);

  auto* fGeLogical = new G4LogicalVolume(geSolid, roiMat, "Ge_log");
  new G4PVPlacement(nullptr, G4ThreeVector(0.0, 0.0, -gegap / 2.0 * cm), fGeLogical,
                    "Ge_phys", fLayerLogical, false, 0, true);

  auto* gapSolid =
    new G4Tubs("Gap", 0.0 * cm, gerad * cm, gegap / 2.0 * cm, 0.0, CLHEP::twopi);

  auto* fGapLogical = new G4LogicalVolume(gapSolid, larMat, "Gap_log");
  new G4PVPlacement(nullptr, G4ThreeVector(0.0, 0.0, gehheight * cm), fGapLogical,
                    "Gap_phys", fLayerLogical, false, 0, true);

  // place layers as mother volume with unique copy number
  G4double step = (gehheight + gegap / 2) * cm;
  G4double xpos;
  G4double ypos;
  G4double angle = CLHEP::twopi / nofStrings;

  // layer logical into ULarlogical
  if(fGeometryName == "baseline_large_reentrance_tube" || fGeometryName == "baseline_large_reentrance_tube_4m_cryo")
    {
      
      G4double length = 16.75 * cm;
      
      G4double vec_main_x = 1/2.;
      G4double vec_main_y = sqrt(1 - vec_main_x*vec_main_x);
      G4double vec_left_x = -1/2.;
      G4double vec_left_y = sqrt(1 - vec_left_x*vec_left_x);
      G4double vec_right_x = 1;
      G4double vec_right_y = 0;
      
      int N = 0;

      //Construct crystal array
      for(G4int i = -4; i <= 4; i++)
	{
	  
	  for(G4int j = 0; j <= 4 && j <= 4 - i; j++)
	    {
	      
	      // leave out calibration ports
	      if( ( i == -3 && j == 0 ) || ( i == 0 && j == 0 ) || ( i == 3 && j == 0) || ( i == -3 && j == 3) || ( i == 0 && j == 3) ) continue;
	      
	      xpos = length * ( vec_main_x * i + vec_left_x * j);
	      ypos = length * ( vec_main_y * i + vec_left_y * j) ;
	      
	      G4cout << "i: " << i << " j: " << j << " | x: " << xpos << " - y: " << ypos << G4endl;
	  
	      for(G4int k = 0; k < nofLayers; k++)
		{
		  // Cube coordinates
		  int x_coor = i;
		  int y_coor = j;
		  int z_coor = -(i + j);
		  int coordinate = N;//(x_coor<0)*1e6 + abs(x_coor)*1e5 + (y_coor<0)*1e4 + abs(y_coor)*1e3 + (z_coor<0)*1e2 + abs(z_coor)*1e1 + k;
		  G4cout << "coordinate: " << coordinate << G4endl;
		  new G4PVPlacement(nullptr,G4ThreeVector(xpos, ypos,
				-step + (nofLayers / 2 * layerthickness - k * layerthickness) * cm - offset_3 * cm),
				fLayerLogical, "Layer_phys", fUlarLogical, false, coordinate, true);
		  N++;
		}
	    }//for(G4int j = 0;...
	  
	  
	  for(G4int j = 1; j <= 4 && j <= 4 - i; j++)
	    {
	      
	      // leave out calibration ports
	      if( ( i == -3 && j == 3 ) || ( i == 0 && j == 3 ) ) continue;
	      
	      xpos = length * ( vec_main_x * i + vec_right_x * j);
	      ypos = length * ( vec_main_y * i + vec_right_y * j) ;
	      
	      for(G4int k = 0; k < nofLayers; k++)
		{
		  // Cube coordinates
		  int x_coor = i + j;
		  int y_coor = -j;
		  int z_coor = -i;
		  int coordinate = N;//(x_coor<0)*1e6 + abs(x_coor)*1e5 + (y_coor<0)*1e4 + abs(y_coor)*1e3 + (z_coor<0)*1e2 + abs(z_coor)*1e1 + k;
		  G4cout << "coordinate: " << coordinate << G4endl;
		  new G4PVPlacement(
				    nullptr,G4ThreeVector(xpos, ypos,
					      -step + (nofLayers / 2 * layerthickness - k * layerthickness) * cm - offset_3 * cm),
				              fLayerLogical, "Layer_phys", fUlarLogical, false, coordinate, true);
		N++;
	      }

	  }//for(G4int j = 1; ...
	  
	}//for(G4int i = -4; ...
      
    }//if(fGeometryName == "baseline_large_reentrance_tube" || ...

  
  else
    {
      for(G4int j = 0; j < nofStrings; j++)
	{
	  xpos = roiradius * cm * std::cos(j * angle);
	  ypos = roiradius * cm * std::sin(j * angle);
	  
	  for(G4int i = 0; i < nofLayers; i++)
	    {
	      new G4PVPlacement(nullptr, G4ThreeVector(xpos, ypos,
                        -step + (nofLayers / 2 * layerthickness - i * layerthickness) * cm - offset_3 * cm),
			fLayerLogical, "Layer_phys", fUlarLogical, false, i + j * nofLayers, true);
	    }
	}
    }

  
  if(fGeometryName == "baseline")
    {
      // placements
      if(fWithBoratedPET == 1)
	new G4PVPlacement(nullptr, G4ThreeVector(ringrad * cm, 0., cushift * cm),
			  fBoratedPETLogical_Tube, "BoratedPET_phys", fLarLogical, false, 0,
			  true);
      
      if(fWithOutCupperTubes == 0)
	new G4PVPlacement(nullptr, G4ThreeVector(ringrad * cm, 0., cushift * cm),
			  fCopperLogical, "Copper_phys", fLarLogical, false, 0, true);
      
      new G4PVPlacement(nullptr, G4ThreeVector(ringrad * cm, 0., cushift * cm),
			fUlarLogical, "ULar_phys", fLarLogical, false, 0, true);
      
      // tower 2
      if(fWithBoratedPET == 1)
	new G4PVPlacement(nullptr, G4ThreeVector(0., ringrad * cm, cushift * cm),
			  fBoratedPETLogical_Tube, "BoratedPET_phys2", fLarLogical, false,
			  1, true);
      
      if(fWithOutCupperTubes == 0)
	new G4PVPlacement(nullptr, G4ThreeVector(0., ringrad * cm, cushift * cm),
			  fCopperLogical, "Copper_phys2", fLarLogical, false, 1, true);
      
      new G4PVPlacement(nullptr, G4ThreeVector(0., ringrad * cm, cushift * cm),
			fUlarLogical, "ULar_phys2", fLarLogical, false, 1, true);
      
      // tower 3
      if(fWithBoratedPET == 1)
	new G4PVPlacement(nullptr, G4ThreeVector(-ringrad * cm, 0., cushift * cm),
			  fBoratedPETLogical_Tube, "BoratedPET_phys3", fLarLogical, false,
			  2, true);
      
      if(fWithOutCupperTubes == 0)
	new G4PVPlacement(nullptr, G4ThreeVector(-ringrad * cm, 0., cushift * cm),
			  fCopperLogical, "Copper_phys3", fLarLogical, false, 2, true);
      
      new G4PVPlacement(nullptr, G4ThreeVector(-ringrad * cm, 0., cushift * cm),
			fUlarLogical, "ULar_phys3", fLarLogical, false, 2, true);
      
      // tower 4
      if(fWithBoratedPET == 1)
	new G4PVPlacement(nullptr, G4ThreeVector(0., -ringrad * cm, cushift * cm),
			  fBoratedPETLogical_Tube, "BoratedPET_phys4", fLarLogical, false,
			  3, true);
      
      if(fWithOutCupperTubes == 0)
	new G4PVPlacement(nullptr, G4ThreeVector(0., -ringrad * cm, cushift * cm),
			  fCopperLogical, "Copper_phys4", fLarLogical, false, 3, true);
      
      new G4PVPlacement(nullptr, G4ThreeVector(0., -ringrad * cm, cushift * cm),
			fUlarLogical, "ULar_phys4", fLarLogical, false, 3, true);
    }//if(fGeometryName == "baseline")

  
  if(fGeometryName == "baseline_smaller" || fGeometryName == "baseline_large_reentrance_tube" || fGeometryName == "baseline_large_reentrance_tube_4m_cryo")
    {
      // placements
      if(fWithBoratedPET == 1)
	new G4PVPlacement(nullptr, G4ThreeVector(0., 0., cushift * cm),
			  fBoratedPETLogical_Tube, "BoratedPET_phys", fLarLogical, false, 0,
			  true);
      
      if(fWithOutCupperTubes == 0)
	new G4PVPlacement(nullptr, G4ThreeVector(0., 0., cushift * cm), fCopperLogical,
			  "Copper_phys", fLarLogical, false, 0, true);
      
      new G4PVPlacement(nullptr, G4ThreeVector(0., 0., cushift * cm), fUlarLogical,
			"ULar_phys", fLarLogical, false, 0, true);
    }


#define whichGeometry 1

#if whichGeometry == 0
  if(fWithBoratedPET == 2)
  {
    int               NPanels        = 28;
    double            radiusOfPanels = 2 * m;
    double            anglePanel     = 360 / 28. * deg;
    G4double          zpos           = 0 * cm;
    G4RotationMatrix* rotMat;
    for(G4int j = 0; j < NPanels; j++)
      {
	xpos   = radiusOfPanels * std::cos(j * anglePanel);
	ypos   = radiusOfPanels * std::sin(j * anglePanel);
	rotMat = new G4RotationMatrix;
	rotMat->rotateZ(-(j + 1) * anglePanel + 90 * deg);
	new G4PVPlacement(rotMat, G4ThreeVector(xpos, ypos, zpos), fBoratedPETLogical_Box,
			  "BoratedPET_phys", fLarLogical, false, j, true);
      }
  }
#endif

#if whichGeometry == 1
  if(fWithBoratedPET == 2 || fWithBoratedPET == 4)
    {
      G4double densityOfBPE   = 0.95;
      double   radiusOfPanels = fBoratedTurbineRadius * cm;
      double   constantAngle  = fBoratedTurbineAngle * deg;  // 45 * deg;
      int      NPanels;
      
      if(fBoratedTurbineNPanels == 0)
	NPanels = ceil(2 * 3.14159265 * radiusOfPanels / cm /
		       (0.95 * b_length / cm * cos(constantAngle)));
      else
	NPanels = fBoratedTurbineNPanels;
      
      fNPanels          = NPanels;
      double anglePanel = 360. / NPanels * deg;
      
      G4double totalVolume =
	NPanels * 2 * b_length / cm * 2 * b_width / cm * 2 * b_height / cm;
      
      G4cout << "Total Mass of B-PE: " << totalVolume * densityOfBPE << G4endl;
      
      
      G4double          zpos = 0 * cm;
      G4RotationMatrix* rotMat;
      
      for(G4int j = 0; j < NPanels; j++)
	{
	  xpos   = radiusOfPanels * std::cos(j * anglePanel);
	  ypos   = radiusOfPanels * std::sin(j * anglePanel);
	  rotMat = new G4RotationMatrix;
	  rotMat->rotateZ(-j * anglePanel + 90 * deg + constantAngle);
	  new G4PVPlacement(rotMat, G4ThreeVector(xpos, ypos, zpos), fBoratedPETLogical_Box,
                        "BoratedPET_phys", fLarLogical, false, j, true);
	}

      
      if(fWithBoratedPET == 4)
	{
	  boratedPETSolid_Tube =
	    new G4Tubs("BoratedPET", 0, (fBoratedTurbineRadius * cm + b_width * 2), b_width,
		       0.0, CLHEP::twopi);
	  fBoratedPETLogical_Tube =
	    new G4LogicalVolume(boratedPETSolid_Tube, BoratedPETMat, "BoratedPET_Logical");
	  
	  new G4PVPlacement(
			    nullptr, G4ThreeVector(0, 0, fBoratedTurbinezPosition * cm - b_height - b_width),
			    fBoratedPETLogical_Tube, "BoratedPET_phys", fLarLogical, false, 0, true);
	  new G4PVPlacement(
			    nullptr, G4ThreeVector(0, 0, fBoratedTurbinezPosition * cm + b_height + b_width),
			    fBoratedPETLogical_Tube, "BoratedPET_phys", fLarLogical, false, 0, true);
	}
    }
#endif

  
  if(fWithBoratedPET == 3)
    {
      G4double densityOfBPE   = 0.95;
      double   radiusOfPanels = fBoratedTurbineRadius * cm;
      
      boratedPETSolid_Tube =
	new G4Tubs("BoratedPET", fBoratedTurbineRadius * cm,
		   (fBoratedTurbineRadius * cm + b_width * 2), b_height, 0.0, CLHEP::twopi);
      fBoratedPETLogical_Tube =
	new G4LogicalVolume(boratedPETSolid_Tube, BoratedPETMat, "BoratedPET_Logical");
      
      G4cout << "Total Mass of B-PE: "
	     << 3.141592653589 * b_height / cm *
	     (pow(fBoratedTurbineRadius + b_width / cm * 2, 2) -
	     pow(fBoratedTurbineRadius, 2)) *
	     densityOfBPE
	     << G4endl;

      new G4PVPlacement(nullptr, G4ThreeVector(0, 0, fBoratedTurbinezPosition * cm),
			fBoratedPETLogical_Tube, "BoratedPET_phys", fLarLogical, false, 0,
			true);
      
      boratedPETSolid_Tube =
	new G4Tubs("BoratedPET", 0, (fBoratedTurbineRadius * cm + b_width * 2), b_width,
		   0.0, CLHEP::twopi);
      fBoratedPETLogical_Tube =
	new G4LogicalVolume(boratedPETSolid_Tube, BoratedPETMat, "BoratedPET_Logical_Lid");
      
      new G4PVPlacement(
			nullptr, G4ThreeVector(0, 0, fBoratedTurbinezPosition * cm - b_height - b_width),
			fBoratedPETLogical_Tube, "BoratedPET_phys", fLarLogical, false, 0, true);
      new G4PVPlacement(
			nullptr, G4ThreeVector(0, 0, fBoratedTurbinezPosition * cm + b_height + b_width),
			fBoratedPETLogical_Tube, "BoratedPET_phys", fLarLogical, false, 0, true);
    }//if(fWithBoratedPET == 3)

  if(fWithBoratedPET == 5)
    {
      //For use with the LEGEND-1000 single reentrance tube design only
      //Shield is made of a hollow n-sided prism, with a hole in the top for the reentrant tube
      //The size of the hole is fixed, but everything else is set by the user, including n


      //For each segment, make 3 basic trapezoids, flat on 4 sides
      //The user inputs a 'radius' for the shield dimensions, but the polygon shield doesn't have a radius...
      //Any regular polygon can be inscribed inside of a circle, so the radius is used as the radius of an inscribing circle

      //Formula for an inscribed regular polygon's side compared to the inscribed circle's radius:
      //l = 2*r*sin(pi/n)

      //All in meters
      G4double shieldheight = fBoratedTurbineHeight/100;//Half-height
      G4double shieldradius = fBoratedTurbineRadius/100;
      G4double shieldthickness = fBoratedTurbineWidth/200;//Half-thickness
      //shieldnsides is set by the user

      //G4double shieldwedgeouterlength = shieldradius*std::sin(CLHEP::pi/shieldnsides);//Half-length
      //G4double shieldwedgeinnerlength = (shieldradius-shieldthickness)*std::sin(CLHEP::pi/shieldnsides);//Half-length
      G4double shieldwedgeouterlength = shieldradius*std::tan(CLHEP::pi/shieldnsides);//Half-length
      G4double shieldwedgeinnerlength = (shieldradius-shieldthickness)*std::tan(CLHEP::pi/shieldnsides);//Half-length

      //Piece of the shield's wall
      G4Trd *shieldwedge = new G4Trd("shieldwedge",shieldheight*m,shieldheight*m,shieldwedgeouterlength*m,shieldwedgeinnerlength*m,shieldthickness*m);

      G4LogicalVolume *shieldwedgelogical =
        new G4LogicalVolume(shieldwedge, BoratedPETMat, "Shieldwedge_log");


      //Piece of the shield's bottom, which should be the same thickness as the wall, and one end of the trapezoid tapers off to 0, making it a triangle
      G4Trd *shieldbottompiece = new G4Trd("shieldbottompiece",shieldthickness*m,shieldthickness*m,shieldwedgeinnerlength*m,0*m,(shieldradius-shieldthickness)/2*m);

      G4LogicalVolume *shieldbottompiecelogical =
        new G4LogicalVolume(shieldbottompiece, BoratedPETMat, "Shieldbottompiece_log");


      //The top is like the bottom, but the bottom pieces all combine to form a solid floor
      //The top stops 1 meter short of the center because of the hole for the reentrance tube
      //Since the trapezoid needs to be 1m shorter, we need to calculate the width of the shortest side, and we need to move its radial position back by 500 cm

      //G4double shieldtopinnerlength = 1*std::sin(CLHEP::pi/shieldnsides);
      G4double shieldtopinnerlength = 1*std::tan(CLHEP::pi/shieldnsides);

      G4Trd *shieldtoppiece = new G4Trd("shieldtoppiece",shieldthickness*m,shieldthickness*m,shieldwedgeinnerlength*m,shieldtopinnerlength*m,(shieldradius-shieldthickness-1)/2*m);

      G4LogicalVolume *shieldtoppiecelogical =
        new G4LogicalVolume(shieldtoppiece, BoratedPETMat, "Shieldtoppiece_log");


      //Each piece needs to be rotated when it is placed
      G4RotationMatrix *rotme;
      

      for(int i = 0; i < shieldnsides; i++)
        {
          rotme  = new G4RotationMatrix(-(360/shieldnsides)*i*deg,90*deg,90*deg);

          //wall
          new G4PVPlacement(
                            rotme, G4ThreeVector(std::sin(2*i*CLHEP::pi/shieldnsides)*(shieldradius-shieldthickness)*m,
                                                 std::cos(2*i*CLHEP::pi/shieldnsides)*(shieldradius-shieldthickness)*m,
                                                 fBoratedTurbinezPosition * cm - b_height - b_width),
                            shieldwedgelogical, "Shield_phys", fLarLogical, false, 0, true);

          //bottom
          new G4PVPlacement(
                            rotme, G4ThreeVector(std::sin(2*i*CLHEP::pi/shieldnsides)*(shieldradius-shieldthickness)/2*m,
                                                 std::cos(2*i*CLHEP::pi/shieldnsides)*(shieldradius-shieldthickness)/2*m,
                                                 fBoratedTurbinezPosition * cm - b_height - b_width - shieldheight*m),
                            shieldbottompiecelogical, "Shield_phys", fLarLogical, false, 0, true);
          //top
          new G4PVPlacement(
                            rotme, G4ThreeVector(std::sin(2*i*CLHEP::pi/shieldnsides)*(shieldradius-shieldthickness+1)/2*m,
                                                 std::cos(2*i*CLHEP::pi/shieldnsides)*(shieldradius-shieldthickness+1)/2*m,
                                                 fBoratedTurbinezPosition * cm - b_height
                                                 - b_width + shieldheight*m),
                            shieldtoppiecelogical, "Shield_phys", fLarLogical, false, 0, true);

        }


    }//if(fWithBoratedPET == 5)

  
  //
  // Visualization attributes
  //
  fWorldLogical->SetVisAttributes(G4VisAttributes::GetInvisible());

  G4Color testColor(0., 109 / 225., 119 / 225.);
  auto*   testVisAtt = new G4VisAttributes(testColor);
  testVisAtt->SetVisibility(true);

  G4Color testColor2(131 / 255., 197 / 225., 190 / 225.);
  auto*   testVisAtt2 = new G4VisAttributes(testColor2);
  testVisAtt2->SetVisibility(true);

  G4Color testColor3(226 / 255., 149 / 225., 120 / 225.);
  auto*   testVisAtt3 = new G4VisAttributes(testColor3);
  testVisAtt3->SetVisibility(true);

  G4Color testColor4(255 / 255., 221 / 225., 210 / 225.);
  auto*   testVisAtt4 = new G4VisAttributes(testColor4);
  testVisAtt4->SetVisibility(true);


  G4Color testColor_copper(1., 190 / 225., 102 / 225.); //255, 190, 182
  auto*   testVisAtt_copper = new G4VisAttributes(testColor_copper);
  testVisAtt_copper->SetVisibility(true);

  G4Color testColor_LAr(153 / 255., 193 / 225., 151 / 225.); //153, 193, 151
  auto*   testVisAtt_LAr = new G4VisAttributes(testColor_LAr);
  testVisAtt_LAr->SetVisibility(true);

  G4Color testColor_Ge(193 / 255., 193 / 225., 193 / 225.); //193, 193, 193
  auto*   testVisAtt_Ge = new G4VisAttributes(testColor_Ge);
  testVisAtt_Ge->SetVisibility(true);

  G4Color testColor_water(170 / 255., 191 / 225., 219 / 225.); //170, 191, 219
  auto*   testVisAtt_water = new G4VisAttributes(testColor_water);
  testVisAtt_water->SetVisibility(true);


  auto* redVisAtt = new G4VisAttributes(G4Colour::Red());
  redVisAtt->SetVisibility(true);
  auto* whiteVisAtt = new G4VisAttributes(G4Colour::White());
  whiteVisAtt->SetVisibility(true);
  auto* orangeVisAtt = new G4VisAttributes(G4Colour::Brown());
  orangeVisAtt->SetVisibility(true);

  auto* greyVisAtt = new G4VisAttributes(G4Colour::Grey());
  greyVisAtt->SetVisibility(true);
  auto* greenVisAtt = new G4VisAttributes(G4Colour::Green());
  greenVisAtt->SetVisibility(true);
  auto* blueVisAtt = new G4VisAttributes(G4Colour::Blue());
  blueVisAtt->SetVisibility(true);

  fCavernLogical->SetVisAttributes(greyVisAtt);
  fHallLogical->SetVisAttributes(whiteVisAtt);
  fTankLogical->SetVisAttributes(greyVisAtt);
  //fWaterLogical->SetVisAttributes(testVisAtt);
  fWaterLogical->SetVisAttributes(testVisAtt_water);
  fLarLogical->SetVisAttributes(testVisAtt2);
  fCoutLogical->SetVisAttributes(greyVisAtt);
  fCvacLogical->SetVisAttributes(greyVisAtt);
  fCinnLogical->SetVisAttributes(greyVisAtt);
  fLidLogical->SetVisAttributes(greyVisAtt);
  fBotLogical->SetVisAttributes(greyVisAtt);
  //fCopperLogical->SetVisAttributes(testVisAtt4);
  fCopperLogical->SetVisAttributes(testVisAtt_copper);
  //fUlarLogical->SetVisAttributes(testVisAtt2);
  fUlarLogical->SetVisAttributes(testVisAtt_LAr);
  //fGapLogical->SetVisAttributes(testVisAtt2);
  fGapLogical->SetVisAttributes(testVisAtt_LAr);
  //fGeLogical->SetVisAttributes(testVisAtt3);
  fGeLogical->SetVisAttributes(testVisAtt_Ge);
  fBoratedPETLogical_Tube->SetVisAttributes(testVisAtt4);
  fBoratedPETLogical_Box->SetVisAttributes(testVisAtt4);

  SetupOpticalProperties();
  
  return fWorldPhysical;

}//SetupBaseline()



auto WLGDDetectorConstruction::SetupHallA() -> G4VPhysicalVolume*
{

  // Full copy of baseline set up but smaller as a Gerda mock-up.

  // Get materials
  worldMaterial = G4Material::GetMaterial("G4_Galactic");
  //auto* larMat        = G4Material::GetMaterial("G4_lAr");
  airMat        = G4Material::GetMaterial("G4_AIR");
  waterMat      = G4Material::GetMaterial("G4_WATER");
  if(fWithGdWater == 1)
    waterMat = water;
  if(fWithWoWater == 1)
    waterMat = airMat;
  steelMat   = G4Material::GetMaterial("G4_STAINLESS-STEEL");
  copperMat  = G4Material::GetMaterial("G4_Cu");
  stdRock    = G4Material::GetMaterial("StdRock");
  roiMat     = G4Material::GetMaterial("enrGe");
  larMat_alt = G4Material::GetMaterial("CombinedArXeHe3");

  if(fXeConc != 0 || fHe3Conc != 0)
    larMat = larMat_alt;

  // constants
  // size parameter, unit [cm]
  G4double offset = 250.0;  // shift cavern floor to keep detector centre at origin
  // cavern
  G4double stone       = 400.0;  // Hall wall thickness 1 m
  G4double hallrad     = 700.0;  // Hall cylinder diam 16 m
  G4double hallhheight = 650.0;  // Hall cylinder height 13 m
  // water tank
  G4double tankwalltop = 0.6;  // water tank thickness at top 6 mm
  G4double tankwallbot = 0.8;  // water tank thickness at bottom 8 mm
  G4double tankrad     = 500;  // water tank diam 10 m
  G4double tankhheight = 400;  // water tank height 8 m
  // cryostat
  G4double cryowall   = 1.5;    // cryostat wall thickness from GERDA
  G4double vacgap     = 1.0;    // vacuum gap between walls
  G4double cryrad     = 200.0;  // cryostat diam 4 m
  G4double cryhheight = 225.0;  // cryostat height 4.5 m
  // Ge cylinder for 35.6 kg at 5.32 g/cm3
  G4double roiradius      = 15.0;                     // detector region diam 30 cm
  G4double roihalfheight  = 20.0;                     // detector region height 40 cm
  G4double gerad          = 3.7;                      // BEGe radius
  G4double gehheight      = 1.5;                      // full height 3 cm
  G4double begegap        = 2.0;                      // gap between BEGe 2cm
  G4double layerthickness = begegap + 2 * gehheight;  // 5 cm total
  G4int    nofLayers      = 8;

  fvertexZ = (hallhheight + offset) * cm;
  fmaxrad  = hallrad * cm;  // 8 m radius

  // Volumes for this geometry

  //
  // World
  //
  auto* worldSolid =
    new G4Tubs("World", 0.0 * cm, (hallrad + stone + 0.1) * cm,
               (hallhheight + stone + offset + 0.1) * cm, 0.0, CLHEP::twopi);
  auto* fWorldLogical  = new G4LogicalVolume(worldSolid, worldMaterial, "World_log");
  auto* fWorldPhysical = new G4PVPlacement(nullptr, G4ThreeVector(), fWorldLogical,
                                           "World_phys", nullptr, false, 0);

  //
  // Cavern
  //
  auto* cavernSolid    = new G4Tubs("Cavern", 0.0 * cm, (hallrad + stone) * cm,
                                    (hallhheight + stone) * cm, 0.0, CLHEP::twopi);
  auto* fCavernLogical = new G4LogicalVolume(cavernSolid, stdRock, "Cavern_log");
  auto* fCavernPhysical =
    new G4PVPlacement(nullptr, G4ThreeVector(0., 0., offset * cm), fCavernLogical,
                      "Cavern_phys", fWorldLogical, false, 0);

  //
  // Hall
  //
  auto* hallSolid =
    new G4Tubs("Hall", 0.0 * cm, hallrad * cm, hallhheight * cm, 0.0, CLHEP::twopi);
  auto* fHallLogical  = new G4LogicalVolume(hallSolid, airMat, "Hall_log");
  auto* fHallPhysical = new G4PVPlacement(nullptr, G4ThreeVector(), fHallLogical,
                                          "Hall_phys", fCavernLogical, false, 0, true);

  //
  // Tank
  //
  auto* tankSolid =
    new G4Cons("Tank", 0.0 * cm, (tankrad + tankwallbot) * cm, 0.0 * cm,
               (tankrad + tankwalltop) * cm, tankhheight * cm, 0.0, CLHEP::twopi);
  auto* fTankLogical = new G4LogicalVolume(tankSolid, steelMat, "Tank_log");
  auto* fTankPhysical =
    new G4PVPlacement(nullptr, G4ThreeVector(0., 0., -offset * cm), fTankLogical,
                      "Tank_phys", fHallLogical, false, 0, true);

  //
  // Water
  //
  auto* waterSolid     = new G4Tubs("Water", 0.0 * cm, tankrad * cm,
                                    (tankhheight - tankwallbot) * cm, 0.0, CLHEP::twopi);
  auto* fWaterLogical  = new G4LogicalVolume(waterSolid, waterMat, "Water_log");
  auto* fWaterPhysical = new G4PVPlacement(nullptr, G4ThreeVector(), fWaterLogical,
                                           "Water_phys", fTankLogical, false, 0, true);

  //
  // outer cryostat
  //
  auto* coutSolid =
    new G4Tubs("Cout", 0.0 * cm, cryrad * cm, cryhheight * cm, 0.0, CLHEP::twopi);
  auto* fCoutLogical  = new G4LogicalVolume(coutSolid, steelMat, "Cout_log");
  auto* fCoutPhysical = new G4PVPlacement(nullptr, G4ThreeVector(), fCoutLogical,
                                          "Cout_phys", fWaterLogical, false, 0, true);

  //
  // vacuum gap
  //
  auto* cvacSolid     = new G4Tubs("Cvac", 0.0 * cm, (cryrad - cryowall) * cm,
                                   (cryhheight - cryowall) * cm, 0.0, CLHEP::twopi);
  auto* fCvacLogical  = new G4LogicalVolume(cvacSolid, worldMaterial, "Cvac_log");
  auto* fCvacPhysical = new G4PVPlacement(nullptr, G4ThreeVector(), fCvacLogical,
                                          "Cvac_phys", fCoutLogical, false, 0, true);

  //
  // inner cryostat
  //
  auto* cinnSolid     = new G4Tubs("Cinn", 0.0 * cm, (cryrad - cryowall - vacgap) * cm,
                                   (cryhheight - cryowall - vacgap) * cm, 0.0, CLHEP::twopi);
  auto* fCinnLogical  = new G4LogicalVolume(cinnSolid, steelMat, "Cinn_log");
  auto* fCinnPhysical = new G4PVPlacement(nullptr, G4ThreeVector(), fCinnLogical,
                                          "Cinn_phys", fCvacLogical, false, 0, true);

  //
  // LAr bath
  //
  auto* larSolid     = new G4Tubs("LAr", 0.0 * cm, (cryrad - 2 * cryowall - vacgap) * cm,
                                  (cryhheight - 2 * cryowall - vacgap) * cm, 0.0, CLHEP::twopi);
  auto* fLarLogical  = new G4LogicalVolume(larSolid, larMat, "Lar_log");
  auto* fLarPhysical = new G4PVPlacement(nullptr, G4ThreeVector(), fLarLogical,
                                         "Lar_phys", fCinnLogical, false, 0, true);

/*
  //
  // cryostat Lid
  //
  auto* lidSolid =
    new G4Tubs("Lid", 0.0 * cm, cryrad * cm, cryowall / 2.0 * cm, 0.0, CLHEP::twopi);
  auto* fLidLogical = new G4LogicalVolume(lidSolid, steelMat, "Lid_log");
  auto* fLidPhysical =
    new G4PVPlacement(nullptr, G4ThreeVector(0., 0., (cryhheight + cryowall / 2.0) * cm),
                      fLidLogical, "Lid_phys", fWaterLogical, false, 0, true);
  auto* fBotLogical = new G4LogicalVolume(lidSolid, steelMat, "Bot_log");
  auto* fBotPhysical =
    new G4PVPlacement(nullptr, G4ThreeVector(0., 0., -(cryhheight + cryowall / 2.0) * cm),
                      fBotLogical, "Bot_phys", fWaterLogical, false, 0, true);
*/

  //
  // String tower
  //
  // auto* towerSolid =
  //   new G4Tubs("String", 0.0 * cm, gerad * cm, roihalfheight * cm, 0.0, CLHEP::twopi);

  // auto* fTowerLogical  = new G4LogicalVolume(towerSolid, larMat, "Tower_log");

  // layers in tower
  auto* layerSolid = new G4Tubs("LayerSolid", 0.0 * cm, gerad * cm,
                                (gehheight + begegap / 2.0) * cm, 0.0, CLHEP::twopi);

  auto* fLayerLogical = new G4LogicalVolume(layerSolid, larMat, "Layer_log");

  // fill one layer
  auto* geSolid =
    new G4Tubs("ROI", 0.0 * cm, gerad * cm, gehheight * cm, 0.0, CLHEP::twopi);

  auto* fGeLogical = new G4LogicalVolume(geSolid, roiMat, "Ge_log");
  new G4PVPlacement(nullptr, G4ThreeVector(0.0, 0.0, -begegap / 2.0 * cm), fGeLogical,
                    "Ge_phys", fLayerLogical, false, 0, true);

  auto* gapSolid =
    new G4Tubs("Gap", 0.0 * cm, gerad * cm, begegap / 2.0 * cm, 0.0, CLHEP::twopi);

  auto* fGapLogical = new G4LogicalVolume(gapSolid, larMat, "Gap_log");
  new G4PVPlacement(nullptr, G4ThreeVector(0.0, 0.0, gehheight * cm), fGapLogical,
                    "Gap_phys", fLayerLogical, false, 0, true);

  // place layers as mother volume with unique copy number
  G4double step = (gehheight + begegap / 2) * cm;
  G4double xpos;
  G4double ypos;
  G4double angle = CLHEP::twopi / 6.0;

  
  if(fGeometryName != "hallA_wo_ge" && fGeometryName != "hallA_only_WLSR")
    {
      for(G4int j = 0; j < 6; j++)
	{
	  xpos = roiradius * cm * std::cos(j * angle);
	  ypos = roiradius * cm * std::sin(j * angle);
	  for(G4int i = 0; i < nofLayers; i++)
	    {
	      new G4PVPlacement(nullptr, G4ThreeVector(xpos, ypos,
                        -step + (nofLayers / 2 * layerthickness - i * layerthickness) * cm),
                        fLayerLogical, "Layer_phys", fLarLogical, false, i + j * nofLayers, true);
	    }
	}
    }
  

  // WLSR volume - its a workaround to access the volume of the WLSR w/o implementing it properly
  if(fGeometryName == "hallA_only_WLSR")
    {
      // cavern
      G4double WLSR_r       = 70.0;  
      G4double WLSR_h       = 300.0/2.; 
      
      auto WLSR_LAr_solid = new G4Tubs("WLSR_LAr_solid", 0.0 * cm, WLSR_r * cm, WLSR_h * cm, 0.0, CLHEP::twopi);  
      auto* fWLSR_LAr_logical = new G4LogicalVolume(WLSR_LAr_solid, larMat, "WLSR_LAr_logical");
      new G4PVPlacement(nullptr, G4ThreeVector(), fWLSR_LAr_logical,
			"WLSR_LAr_physical", fLarLogical, false, 0, true);
    }

  
  //
  // Visualization attributes
  //
  fWorldLogical->SetVisAttributes(G4VisAttributes::GetInvisible());

  auto* redVisAtt = new G4VisAttributes(G4Colour::Red());
  redVisAtt->SetVisibility(true);
  auto* greyVisAtt = new G4VisAttributes(G4Colour::Grey());
  greyVisAtt->SetVisibility(true);
  auto* greenVisAtt = new G4VisAttributes(G4Colour::Green());
  greenVisAtt->SetVisibility(true);
  auto* blueVisAtt = new G4VisAttributes(G4Colour::Blue());
  blueVisAtt->SetVisibility(true);

  fCavernLogical->SetVisAttributes(redVisAtt);
  fHallLogical->SetVisAttributes(greyVisAtt);
  fTankLogical->SetVisAttributes(greenVisAtt);
  fWaterLogical->SetVisAttributes(greyVisAtt);
  fLarLogical->SetVisAttributes(greyVisAtt);
  fCoutLogical->SetVisAttributes(blueVisAtt);
  fCvacLogical->SetVisAttributes(greyVisAtt);
  fCinnLogical->SetVisAttributes(blueVisAtt);
  //fLidLogical->SetVisAttributes(blueVisAtt);
  //fBotLogical->SetVisAttributes(blueVisAtt);
  fLayerLogical->SetVisAttributes(G4VisAttributes::Invisible);
  fGapLogical->SetVisAttributes(greyVisAtt);
  fGeLogical->SetVisAttributes(redVisAtt);

  return fWorldPhysical;

}//SetupHallA()

void WLGDDetectorConstruction::SetupOpticalProperties(void)
{
  // - Gives the LAr scintillation properties (now)
  // - Gives other materials within the cryostat optical properties (later)
  // - Makes the code slow as hell if applied at the inappropriate time
  //Most of this is adapted from some of the MaGe source code

  G4cout << G4endl << "Optical physics enabled" << G4endl << G4endl;
  
  using namespace CLHEP;

  //Some constants which will be useful throughout

  static const G4double LambdaE = twopi *1.973269602e-16 * m * GeV;//Planck constant, for energy-wavelength conversions
  const G4int           ArEntries = 69;
  static const G4double temperature = 88.5*kelvin;


  //Single-value properties, which can be taken and applied directly
  G4double scint_yield = 23.6*eV;  // Nominal energy to produce a photon
  G4double photon_yield = 1.0*MeV/scint_yield;
  G4double tau_s = 5.95*ns;
  G4double tau_l = 1590.*ns;
  G4double fano = 0.11;
  G4double yieldratio = 0.3;

  G4double EnergyMax = LambdaE/(115*nm);
  G4double EnergyMin = LambdaE/(650*nm);
  G4double dE = (EnergyMax-EnergyMin)/(ArEntries - 1);

  G4double EnergyMinScint = LambdaE/(136*nm);
  G4double dES = (EnergyMax-EnergyMinScint)/(ArEntries - 1);

  //Ideally, the optical properties of a material can be described by a smoothly varying function of the photon energy
  //Unfortunately, this is not how things are done in Geant4 at the moment
  //Instead, the user generates pairs of arrays; the first array has a list of photon energies,
  //and the second one has some property of the photon at that energy
  //In the case where the arrays have an infinite number of entries, we recover the smooth function
  //The dE and dES variables represent the spacing between each array entry (dES is for scintillation)
  //ArEntries is obviously the number of entries in these arrays

  G4double LArEnergyArray[ArEntries];
  G4double LArRIArray[ArEntries];//Refractive index
  G4double LArRaylArray[ArEntries];//Rayleigh scattering length
  G4double LArAbsArray[ArEntries];//Absorption length

  G4double LArScintEnergyArray[ArEntries];
  G4double LArScintArray[ArEntries];


  //Two energy ranges - one for the scintillation spectrum (narrow), one for the optical properties (broad)
  G4double E = EnergyMin - dE;//To offset the first iteration of the next loop
  G4double EScint = EnergyMinScint - dES;

  for(int i = 0; i < ArEntries; i++)
    {

      E+=dE;//Starts at EnergyMin
      EScint+=dES;//Starts at EnergyMinScint

      LArEnergyArray[i] = E;
      LArScintEnergyArray[i] = EScint;

      //Functions are all dependent on wavelength
      LArRIArray[i] = LArRefIndex((LambdaE/E));
      LArRaylArray[i] = LArRayLength((LambdaE/E),temperature);
      LArAbsArray[i] = LArAbsLength((LambdaE/E));

      LArScintArray[i] = LArScintSpec((LambdaE/EScint));

    }//for(int i

  //Finally, with all the optical properties calculated, create a G4 object to store them and assign them to the LAr
  
  G4MaterialPropertiesTable* ArMPT = new G4MaterialPropertiesTable();
  
  //Add array properties
  ArMPT->AddProperty("RINDEX",   LArEnergyArray,LArRIArray,ArEntries);
  ArMPT->AddProperty("RAYLEIGH", LArEnergyArray,LArRaylArray,ArEntries);
  ArMPT->AddProperty("ABSLENGTH",LArEnergyArray,LArAbsArray,ArEntries);
  ArMPT->AddProperty("FASTCOMPONENT",LArScintEnergyArray,LArScintArray,ArEntries);
  ArMPT->AddProperty("SLOWCOMPONENT",LArScintEnergyArray,LArScintArray,ArEntries);

  //Add single-value properties
  ArMPT->AddConstProperty("SCINTILLATIONYIELD",photon_yield);
  ArMPT->AddConstProperty("FASTTIMECONSTANT",tau_s);
  ArMPT->AddConstProperty("SLOWTIMECONSTANT",tau_l);
  ArMPT->AddConstProperty("YIELDRATIO",yieldratio);
  ArMPT->AddConstProperty("RESOLUTIONSCALE",fano);


  larMat->SetMaterialPropertiesTable(ArMPT);
  larMat->GetIonisation()->SetBirksConstant(5.1748e-4*cm/MeV);



  //For most other materials, we can simply set the absorption length to be arbitrarily small
  //For the fibers, we need to be more careful
  //Also, for all materials within the cryostat, we have to define the reflection properties


  //Ge optical properties

  const static int GeEntries = 423;
  
  G4double GeWavelengthArray[GeEntries] = {110,279,280,281,282,283,284,285,286,287,288,289,290,291,292,293,294,295,296,297,298,299,300,301,302,303,304,305,306,307,308,309,310,311,312,313,314,315,316,317,318,319,320,321,322,323,324,325,326,327,328,329,330,331,332,333,334,335,336,337,338,339,340,341,342,343,344,345,346,347,348,349,350,351,352,353,354,355,356,357,358,359,360,361,362,363,364,365,366,367,368,369,370,371,372,373,374,375,376,377,378,379,380,381,382,383,384,385,386,387,388,389,390,391,392,393,394,395,396,397,398,399,400,401,402,403,404,405,406,407,408,409,410,411,412,413,414,415,416,417,418,419,420,421,422,423,424,425,426,427,428,429,430,431,432,433,434,435,436,437,438,439,440,441,442,443,444,445,446,447,448,449,450,451,452,453,454,455,456,457,458,459,460,461,462,463,464,465,466,467,468,469,470,471,472,473,474,475,476,477,478,479,480,481,482,483,484,485,486,487,488,489,490,491,492,493,494,495,496,497,498,499,500,501,502,503,504,505,506,507,508,509,510,511,512,513,514,515,516,517,518,519,520,521,522,523,524,525,526,527,528,529,530,531,532,533,534,535,536,537,538,539,540,541,542,543,544,545,546,547,548,549,550,551,552,553,554,555,556,557,558,559,560,561,562,563,564,565,566,567,568,569,570,571,572,573,574,575,576,577,578,579,580,581,582,583,584,585,586,587,588,589,590,591,592,593,594,595,596,597,598,599,600,601,602,603,604,605,606,607,608,609,610,611,612,613,614,615,616,617,618,619,620,621,622,623,624,625,626,627,628,629,630,631,632,633,634,635,636,637,638,639,640,641,642,643,644,645,646,647,648,649,650,651,652,653,654,655,656,657,658,659,660,661,662,663,664,665,666,667,668,669,670,671,672,673,674,675,676,677,678,679,680,681,682,683,684,685,686,687,688,689,690,691,692,693,694,695,696,697,698,699,700};//Wavelengths used to determine the appropriate energy for each reflectivity value

  G4double GeReflArray[GeEntries] = {0.65,0.65,0.5869,0.5842,0.5822,0.574,0.5694,0.5647,0.5575,0.5498,0.542,0.5341,0.5239,0.517,0.511,0.504,0.4973,0.4928,0.487,0.4831,0.4763,0.4743,0.4688,0.4653,0.4615,0.4586,0.4556,0.453,0.4507,0.4477,0.4441,0.4419,0.4387,0.4372,0.4344,0.4295,0.429,0.4261,0.4239,0.4262,0.4236,0.4226,0.4175,0.4167,0.4128,0.4141,0.4119,0.4107,0.4088,0.4086,0.4062,0.4078,0.4018,0.4041,0.4042,0.4012,0.4008,0.3997,0.4018,0.3976,0.3943,0.4,0.3961,0.3924,0.3981,0.3953,0.3914,0.3927,0.3921,0.3893,0.3902,0.3887,0.3891,0.389,0.3872,0.3887,0.3885,0.3863,0.3861,0.3872,0.3859,0.3848,0.3849,0.3824,0.3847,0.3845,0.3835,0.3822,0.3828,0.3821,0.3804,0.381,0.3793,0.3806,0.3792,0.3779,0.3777,0.3783,0.3764,0.376,0.3755,0.3752,0.3749,0.3734,0.3728,0.3734,0.3709,0.3689,0.3709,0.369,0.3704,0.3695,0.3678,0.3669,0.365,0.3662,0.3637,0.3636,0.3642,0.3609,0.3624,0.3612,0.3599,0.3596,0.3578,0.3575,0.3587,0.3563,0.3564,0.3559,0.3561,0.3548,0.3527,0.3534,0.353,0.3525,0.3524,0.3523,0.3512,0.3507,0.3514,0.3522,0.3505,0.3503,0.3502,0.3507,0.3499,0.3509,0.3506,0.3504,0.352,0.351,0.3507,0.3502,0.3511,0.352,0.3525,0.3523,0.3527,0.3523,0.3517,0.3518,0.352,0.3521,0.3532,0.3534,0.3538,0.355,0.3557,0.3556,0.3563,0.3562,0.3563,0.3574,0.3577,0.3577,0.3581,0.358,0.3589,0.3588,0.3595,0.3606,0.3606,0.3608,0.3613,0.3616,0.3618,0.3624,0.3635,0.363,0.3646,0.3644,0.3651,0.3647,0.3661,0.3664,0.3672,0.3674,0.3678,0.3686,0.369,0.3694,0.3707,0.3706,0.3717,0.3717,0.372,0.3727,0.3734,0.3736,0.3743,0.3755,0.3756,0.3762,0.3771,0.377,0.3781,0.3784,0.3793,0.38,0.3807,0.3809,0.3813,0.3829,0.3823,0.3839,0.3849,0.3848,0.3857,0.3862,0.3871,0.3873,0.3883,0.3891,0.3899,0.3911,0.3909,0.3923,0.393,0.3931,0.3945,0.3951,0.3955,0.3964,0.3969,0.3981,0.3985,0.3991,0.3998,0.4002,0.4005,0.402,0.402,0.4022,0.4034,0.4037,0.4035,0.4043,0.4042,0.405,0.4051,0.4046,0.4052,0.4049,0.4047,0.4039,0.4044,0.4043,0.4037,0.4033,0.4037,0.4034,0.4035,0.4035,0.403,0.4029,0.4034,0.4026,0.4028,0.4028,0.4027,0.403,0.4031,0.4031,0.4033,0.4045,0.4037,0.4044,0.4042,0.4039,0.4052,0.4054,0.4062,0.4062,0.407,0.4075,0.4083,0.4085,0.4094,0.4101,0.4105,0.4114,0.4115,0.4118,0.4126,0.4127,0.4134,0.4133,0.4134,0.4135,0.4132,0.413,0.4136,0.4129,0.4124,0.4119,0.4113,0.4108,0.4103,0.4098,0.4089,0.4084,0.4081,0.4076,0.4058,0.4053,0.4048,0.4044,0.4033,0.4026,0.4017,0.4008,0.3997,0.3993,0.3982,0.3975,0.3967,0.3964,0.3957,0.3948,0.394,0.3933,0.3919,0.3916,0.3908,0.3903,0.3894,0.3886,0.3879,0.3876,0.3869,0.3861,0.3855,0.3847,0.3841,0.3832,0.3827,0.3822,0.3817,0.3812,0.3809,0.3795,0.3793,0.3788,0.3786,0.3781,0.3772,0.3766,0.3763,0.3758,0.3752,0.3743,0.3741,0.3739,0.3732,0.3724,0.3722,0.3715,0.3714,0.3705,0.3701,0.3698,0.3698,0.3692,0.3683,0.3685,0.3678,0.3665,0.3671,0.3667,0.3651,0.3654,0.3648,0.3645,0.3642,0.3636,0.3639,0.3628,0.3628,0.3621,0.3615,0.3615,0.3612,0.3606,0.3606,0.3603,0.3601,0.3592,0.3588,0.3585,0.3585,0.3576,0.3576,0.3573,0.3567,0.3575,0.3565,0.3559,0.3558,0.3552,0.3546,0.3551,0.3537};//Reflectivity

  G4double GeEnergyArray[GeEntries];
  G4double GeAbsArray[GeEntries];//Absorption length

  for(int i = 0; i < GeEntries; i++)
    {
      GeEnergyArray[i] = (LambdaE/(GeWavelengthArray[i]*nm));      
      GeAbsArray[i] = 1*nm;//Arbitrarily small
    }

  G4MaterialPropertiesTable* GeMPT = new G4MaterialPropertiesTable();

  //Add array properties
  GeMPT->AddProperty("REFLECTIVITY",GeEnergyArray,GeReflArray,GeEntries);
  GeMPT->AddProperty("ABSLENGTH",   GeEnergyArray,GeAbsArray ,GeEntries);
  //Refractive index will be calculated automatically later, using the surface definitions
  roiMat = G4Material::GetMaterial("enrGe");  
  roiMat->SetMaterialPropertiesTable(GeMPT);


  //Cu optical properties
  const static int CuEntries = 424;

  G4double CuWavelengthArray[CuEntries] = {110,150,180,280,281,282,283,284,285,286,287,288,289,290,291,292,293,294,295,296,297,298,299,300,301,302,303,304,305,306,307,308,309,310,311,312,313,314,315,316,317,318,319,320,321,322,323,324,325,326,327,328,329,330,331,332,333,334,335,336,337,338,339,340,341,342,343,344,345,346,347,348,349,350,351,352,353,354,355,356,357,358,359,360,361,362,363,364,365,366,367,368,369,370,371,372,373,374,375,376,377,378,379,380,381,382,383,384,385,386,387,388,389,390,391,392,393,394,395,396,397,398,399,400,401,402,403,404,405,406,407,408,409,410,411,412,413,414,415,416,417,418,419,420,421,422,423,424,425,426,427,428,429,430,431,432,433,434,435,436,437,438,439,440,441,442,443,444,445,446,447,448,449,450,451,452,453,454,455,456,457,458,459,460,461,462,463,464,465,466,467,468,469,470,471,472,473,474,475,476,477,478,479,480,481,482,483,484,485,486,487,488,489,490,491,492,493,494,495,496,497,498,499,500,501,502,503,504,505,506,507,508,509,510,511,512,513,514,515,516,517,518,519,520,521,522,523,524,525,526,527,528,529,530,531,532,533,534,535,536,537,538,539,540,541,542,543,544,545,546,547,548,549,550,551,552,553,554,555,556,557,558,559,560,561,562,563,564,565,566,567,568,569,570,571,572,573,574,575,576,577,578,579,580,581,582,583,584,585,586,587,588,589,590,591,592,593,594,595,596,597,598,599,600,601,602,603,604,605,606,607,608,609,610,611,612,613,614,615,616,617,618,619,620,621,622,623,624,625,626,627,628,629,630,631,632,633,634,635,636,637,638,639,640,641,642,643,644,645,646,647,648,649,650,651,652,653,654,655,656,657,658,659,660,661,662,663,664,665,666,667,668,669,670,671,672,673,674,675,676,677,678,679,680,681,682,683,684,685,686,687,688,689,690,691,692,693,694,695,696,697,698,699,700};


  G4double CuReflArray[CuEntries] = {0.1,0.2,0.3,0.2698,0.2687,0.2688,0.2681,0.2649,0.2639,0.2637,0.2641,0.2611,0.2602,0.2598,0.2608,0.2596,0.2599,0.2628,0.2643,0.2624,0.2634,0.2663,0.2671,0.2679,0.2721,0.2723,0.2744,0.2744,0.2771,0.2785,0.2812,0.2815,0.2829,0.2829,0.2863,0.2863,0.2887,0.2875,0.2877,0.2878,0.2906,0.2916,0.2925,0.2944,0.2968,0.2954,0.2958,0.2963,0.2971,0.2996,0.298,0.2999,0.3051,0.302,0.3044,0.3035,0.3059,0.3072,0.3063,0.3067,0.3077,0.3076,0.31,0.3147,0.3114,0.3102,0.312,0.3162,0.3198,0.32,0.3193,0.3212,0.3203,0.3212,0.3219,0.324,0.3244,0.3258,0.3267,0.3271,0.329,0.3307,0.3314,0.3324,0.3337,0.3343,0.3359,0.3383,0.3386,0.3407,0.3422,0.3429,0.3448,0.3442,0.3468,0.3473,0.3483,0.3504,0.352,0.3523,0.355,0.3562,0.3569,0.359,0.3602,0.3593,0.3613,0.3629,0.3637,0.3642,0.3667,0.3655,0.3683,0.3704,0.3712,0.3726,0.3744,0.3752,0.3763,0.3777,0.3786,0.3802,0.3805,0.3821,0.3839,0.3855,0.3871,0.3877,0.3894,0.3907,0.3908,0.3929,0.3932,0.3954,0.3967,0.3984,0.3995,0.4005,0.4019,0.4034,0.4042,0.4065,0.4079,0.4096,0.4095,0.4121,0.4128,0.4135,0.4159,0.4177,0.4178,0.4192,0.4211,0.4222,0.4226,0.4247,0.4255,0.4268,0.4282,0.4291,0.4307,0.4317,0.433,0.4336,0.4348,0.4359,0.4372,0.4385,0.4393,0.4398,0.4421,0.443,0.4442,0.4448,0.4467,0.4464,0.448,0.4493,0.4504,0.4509,0.4521,0.4525,0.4544,0.4551,0.456,0.4567,0.4576,0.4583,0.4592,0.4605,0.4611,0.4616,0.462,0.4639,0.4642,0.4656,0.4657,0.4669,0.4667,0.468,0.4687,0.4693,0.4704,0.4706,0.471,0.4715,0.472,0.4731,0.4738,0.4744,0.4753,0.4758,0.476,0.4762,0.477,0.4777,0.4784,0.479,0.4792,0.48,0.4804,0.4807,0.4814,0.4817,0.4818,0.4814,0.4836,0.4837,0.4844,0.4846,0.4838,0.4846,0.4852,0.4858,0.4856,0.4863,0.4863,0.4868,0.4876,0.4882,0.4885,0.489,0.4898,0.4898,0.4907,0.4908,0.4919,0.4913,0.4925,0.4929,0.4937,0.4943,0.4942,0.4954,0.4962,0.4968,0.4972,0.498,0.5,0.5,0.5014,0.5022,0.5029,0.5046,0.5056,0.5066,0.5082,0.5094,0.5118,0.5138,0.5154,0.5174,0.5195,0.5222,0.5246,0.5276,0.5308,0.5342,0.5381,0.5415,0.546,0.5495,0.5543,0.5604,0.566,0.5725,0.5767,0.5837,0.5895,0.5961,0.6039,0.6113,0.6196,0.6266,0.6346,0.6423,0.6512,0.6599,0.6682,0.6763,0.6844,0.6936,0.701,0.709,0.7164,0.7255,0.7332,0.7408,0.7484,0.7556,0.7621,0.7683,0.7754,0.7805,0.7876,0.7936,0.7985,0.8036,0.8097,0.814,0.8183,0.8228,0.827,0.8308,0.834,0.8376,0.8416,0.8445,0.8469,0.8514,0.8534,0.8556,0.8576,0.8606,0.8625,0.8644,0.8671,0.8686,0.8709,0.8725,0.8733,0.8756,0.8777,0.8789,0.8805,0.881,0.8826,0.8842,0.8859,0.8872,0.8885,0.889,0.8904,0.8916,0.8928,0.8935,0.8946,0.8961,0.897,0.8974,0.8978,0.8992,0.8997,0.9009,0.901,0.9023,0.903,0.9039,0.9044,0.9057,0.9054,0.9069,0.9068,0.9086,0.9086,0.9099,0.9107,0.9114,0.9115,0.9116,0.9127,0.9122,0.9132,0.914,0.9151,0.9147,0.9158,0.917,0.9172,0.9167,0.9175,0.9176,0.9181,0.9188,0.9187,0.9187,0.9204,0.9205,0.9215,0.9204,0.9224,0.9226,0.9226,0.9213,0.9241,0.9234,0.9234,0.9242,0.925,0.925,0.9253,0.926,0.9268,0.9255,0.9265,0.9287,0.9282,0.9283,0.928,0.9294,0.9309,0.9294,0.9298,0.9303};

  G4double CuEnergyArray[CuEntries];
  G4double CuAbsArray[CuEntries];//Absorption length

  for(int i = 0; i < CuEntries; i++)
    {
      CuEnergyArray[i] = (LambdaE/(CuWavelengthArray[i]*nm));      
      CuAbsArray[i] = 1*nm;//Arbitrarily small
    }

  G4MaterialPropertiesTable* CuMPT = new G4MaterialPropertiesTable();

  //Add array properties
  CuMPT->AddProperty("REFLECTIVITY",CuEnergyArray,CuReflArray,CuEntries);
  CuMPT->AddProperty("ABSLENGTH",   CuEnergyArray,CuAbsArray ,CuEntries);

  copperMat = G4Material::GetMaterial("G4_Cu");
  copperMat->SetMaterialPropertiesTable(CuMPT);
  

  //For the steel, let Geant4 calculate the reflectivity automatically using the real and imaginary refractive indices

  const static int SSEntries = 28;

  G4double SSWavelengthArray[SSEntries] = {103.32,107.812,112.713,118.08,123.984,130.51,137.76,145.864,154.98,165.312,177.12,190.745,206.64,225.426,247.968,261.019,275.52,291.728,309.96,330.625,354.241,381.49,413.281,450.852,495.937,551.041,619.921,708.481};
  G4double SSRealRIArray[SSEntries] = {0.8523,0.8791,0.9105,0.9448,0.9811,1.0174,1.0508,1.078,1.0952,1.0979,1.0834,1.0524,1.0114,0.9786,0.9848,1.0157,1.0709,1.1477,1.2262,1.2657,1.2321,1.1422,1.0547,1.0197,1.0666,1.2276,1.5753,2.3096};
  G4double SSImagRIArray[SSEntries] = {1.1434,1.1966,1.2494,1.3003,1.3485,1.3943,1.4374,1.4792,1.523,1.5749,1.6440,1.7437,1.8905,2.1009,2.3822,2.5445,2.7112,2.8678,2.9991,3.1106,3.2521,3.4954,3.8792,4.4072,5.0863,5.9541,7.0964,8.6751};
  G4double SSAbsArray[SSEntries];
  G4double SSEnergyArray[SSEntries];

  for(int i = 0; i < SSEntries; i++)
    {
      SSEnergyArray[i] = (LambdaE/(SSWavelengthArray[i]*nm));
      SSAbsArray[i] = 1*nm;
    }

  G4MaterialPropertiesTable* SSMPT = new G4MaterialPropertiesTable();

  //Add array properties
  SSMPT->AddProperty("REALRINDEX",     SSEnergyArray,SSRealRIArray,SSEntries);
  SSMPT->AddProperty("IMAGINARYRINDEX",SSEnergyArray,SSImagRIArray,SSEntries);
  SSMPT->AddProperty("ABSLENGTH",      SSEnergyArray,SSAbsArray,SSEntries);
  steelMat = G4Material::GetMaterial("G4_STAINLESS-STEEL");
  steelMat->SetMaterialPropertiesTable(SSMPT);

  
  //For the Tetratex, according to Luigi's notes, the reflectivity might be a bit off due to the TPB coating
  //I think this should be a pretty minor effect, though...

  const static int TTEntries = 7;

  G4double TTWavelengthArray[TTEntries] = {250,260,400,560,680,780,800};
  G4double TTReflArray[TTEntries] = {0.975,0.973,0.95,0.925,0.9,0.875,0.865};
  
  //Most of the boundaries are between the LAr and something else
  //LAr to Ge, LAr to Cu, LAr to steel, LAr to moderator, LAr to fibers
  //As a result, these will be defined first
  //The remainder of the boundaries are mostly inside the fiber
  //The fiber needs a carefully defined, complex structure, since this is the most important area for optical sims
  //The fibers will be the last volumes to receive the optical treatment


  
}//SetupOpticalProperties


double WLGDDetectorConstruction::LArEpsilon(double lambda)
{

  // Calculates the dielectric constant of LAr from the Bideau-Sellmeier formula.
  // See : A. Bideau-Mehu et al., "Measurement of refractive indices of Ne, Ar,
  // Kr and Xe ...", J. Quant. Spectrosc. Radiat. Transfer, Vol. 25 (1981), 395

  G4double epsilon;
  if (lambda < 110*nanometer) return 1.0e4; // lambda MUST be > 110.0 nm
  epsilon = lambda / micrometer; // switch to micrometers
  epsilon = 1.0 / (epsilon * epsilon); // 1 / (lambda)^2
  epsilon = 1.2055e-2 * ( 0.2075 / (91.012 - epsilon) +
                          0.0415 / (87.892 - epsilon) +
                          4.3330 / (214.02 - epsilon) );
  epsilon *= (8./12.); // Bideau-Sellmeier -> Clausius-Mossotti
  G4double LArRho = 1.396*g/cm3;
  G4double GArRho = 1.66e-03*g/cm3;
  epsilon *= (LArRho / GArRho); // density correction (Ar gas -> LAr liquid)
  if (epsilon < 0.0 || epsilon > 0.999999) return 4.0e6;
  epsilon = (1.0 + 2.0 * epsilon) / (1.0 - epsilon); // solve Clausius-Mossotti
  return epsilon;

}//LarEpsilon



double WLGDDetectorConstruction::LArRefIndex(double LambdaE)
{
  return sqrt(LArEpsilon(LambdaE)); // square root of dielectric constant
}//LarRefIndex



double WLGDDetectorConstruction::LArRayLength(double LambdaE, double temperature)
{
  // Calculates the Rayleigh scattering length using equations given in
  // G. M. Seidel at al., "Rayleigh scattering in rare-gas liquids",
  // arXiv:hep-ex/0111054 v2 22 Apr 2002

G4double dyne = 1.0e-5*newton;
  static const G4double LArKT = 2.18e-10 * cm2/dyne; // LAr isothermal compressibility
  static const G4double k = 1.380658e-23 * joule/kelvin; // the Boltzmann constant
  G4double h;
  h = LArEpsilon(LambdaE);
  if (h < 1.00000001) h = 1.00000001; // just a precaution

  h = (h - 1.0) * (h + 2.0); // the "dielectric constant" dependance
  h *= h; // take the square
  h *= LArKT * temperature * k; // compressibility * temp * Boltzmann constant
  h /= LambdaE * LambdaE * LambdaE * LambdaE; // (lambda)^4
  h *= 9.18704494231105429; // (2 * Pi / 3)^3

  if ( h < (1.0 / (10.0 * km)) ) h = 1.0 / (10.0 * km); // just a precaution
  if ( h > (1.0 / (0.1 * nanometer)) ) h = 1.0 / (0.1 * nanometer); // just a precaution
  return ( 1.0 / h );
}//LarRayLength



double WLGDDetectorConstruction::LArAbsLength(double LambdaE)
{

  //For now, just implement absorption length as a stepwise function
  //For VUV photons (under 200 nm wavelength), the attenuation length is 60 cm
  //For visible spectrum photons, the LAr is (ideally) transparent
  //This second assumption depends on the purity of the LAr inside the tank
  //Could possibly implement a more realistic value based on L200 measurements
  //For L1000, situation is more complicated - two types of LAr (atm and UG)

  G4double LArAbsL = 0;
  G4double LArAbsVUV = 60*cm;
  G4double LArAbsVis = 1000*m;


  if ((LambdaE/nm) < 200.0)
    LArAbsL = LArAbsVUV;
  else
    LArAbsL = LArAbsVis;

  return LArAbsL;

}//LarAbsLength



double WLGDDetectorConstruction::LArScintSpec(double LambdaES)
{

  G4double waveL;
  waveL =exp(-0.5*((LambdaES-128.0)/(2.929))*((LambdaES-128.0)/(2.929)));
  return waveL;

}//LarScintSpec


void WLGDDetectorConstruction::SetPositionOfDetectors(G4String name)
{

  std::set<G4String> knownGeometries = { "baseline", "original" };
  
  if(knownGeometries.count(name) == 0)
  {
    G4Exception("WLGDDetectorConstruction::SetGeometry", "WLGD0001", JustWarning,
                ("Invalid geometry setup name '" + name + "'").c_str());
    return;
  }

  fDetectorPosition = name;
  // Reinit wiping out stores
  G4RunManager::GetRunManager()->ReinitializeGeometry();

}



void WLGDDetectorConstruction::SetGeometry(const G4String& name)
{

  std::set<G4String> knownGeometries = { "baseline", "baseline_smaller", "baseline_large_reentrance_tube", "alternative",
                                         "hallA", "hallA_wo_ge", "hallA_only_WLSR", "baseline_large_reentrance_tube_4m_cryo" };
  if(knownGeometries.count(name) == 0)
    {
      G4Exception("WLGDDetectorConstruction::SetGeometry", "WLGD0001", JustWarning,
		  ("Invalid geometry setup name '" + name + "'").c_str());
      return;
    }

  fGeometryName = name;
  // Reinit wiping out stores
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}



void WLGDDetectorConstruction::ExportGeometry(const G4String& file)
{

  G4GDMLParser parser;
  parser.Write(file);
}



void WLGDDetectorConstruction::SetNeutronBiasFactor(G4double nf) { fNeutronBias = nf; }

void WLGDDetectorConstruction::SetMuonBiasFactor(G4double mf) { fMuonBias = mf; }

void WLGDDetectorConstruction::SetNeutronYieldBias(G4double nf) { fNeutronYieldBias = nf; }


// Additional settings for adjusting the detector geometry

// changing the concentration of Xe in the LAr
void WLGDDetectorConstruction::SetXeConc(G4double nf)
{
  fXeConc = nf * 1e-3;
  WLGDDetectorConstruction::DefineMaterials();
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

// changing the concentration of He-3 in the LAr
void WLGDDetectorConstruction::SetHe3Conc(G4double nf)
{
  fHe3Conc = nf * 1e-3;
  WLGDDetectorConstruction::DefineMaterials();
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

// changing the radius of the cryostat
void WLGDDetectorConstruction::SetOuterCryostatRadius(G4double rad)
{
  fCryostatOuterRadius = rad;
  WLGDDetectorConstruction::DefineMaterials();
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

// changing the height of the cryostat
void WLGDDetectorConstruction::SetCryostatHeight(G4double height)
{
  fCryostatHeight = height;
  WLGDDetectorConstruction::DefineMaterials();
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

// option to remove the cupper tubes to see what impact it has on the production rate
void WLGDDetectorConstruction::SetWithoutCupperTubes(G4int answer)
{
  fWithOutCupperTubes = answer;
  WLGDDetectorConstruction::DefineMaterials();
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

// option to include borated PE in the setup (1: tubes around the re-entrance tubes, 2:
// turbine structure)
void WLGDDetectorConstruction::SetNeutronModerator(G4int answer)
{
  fWithBoratedPET = answer;
  WLGDDetectorConstruction::DefineMaterials();
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

// option to include borated PE in the setup (1: tubes around the re-entrance tubes, 2:
// turbine structure)
void WLGDDetectorConstruction::SetMaterial(G4String answer)
{
  fSetMaterial = answer;
  WLGDDetectorConstruction::DefineMaterials();
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

// option to set the radius of the turbine structure
void WLGDDetectorConstruction::SetTurbineAndTubeRadius(G4double radius)
{
  fBoratedTurbineRadius = radius;
}

// option to set the length of the turbine structure
void WLGDDetectorConstruction::SetTurbineAndTubeLength(G4double length)
{
  fBoratedTurbineLength = length;
}

// option to set the width of the turbine structure
void WLGDDetectorConstruction::SetTurbineAndTubeWidth(G4double width)
{
  fBoratedTurbineWidth = width;
}

// option to set the angle of the turbine structure
void WLGDDetectorConstruction::SetTurbineAndTubeAngle(G4double deg)
{
  fBoratedTurbineAngle = deg;
}

// option to set the height of the turbine structure
void WLGDDetectorConstruction::SetTurbineAndTubeHeight(G4double height)
{
  fBoratedTurbineHeight = height;
}

// option to set the zPosition of the turbine structure
void WLGDDetectorConstruction::SetTurbineAndTubezPosition(G4double zPosition)
{
  fBoratedTurbinezPosition = zPosition;
}

// option to set the zPosition of the turbine structure
void WLGDDetectorConstruction::SetTurbineAndTubeNPanels(G4double nPanels)
{
  fBoratedTurbineNPanels = nPanels;
}

// option to change the water to gadolinium weighted water in the water tank
void WLGDDetectorConstruction::SetGdWater(G4int answer)
{
  fWithGdWater = answer;
  WLGDDetectorConstruction::DefineMaterials();
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

void WLGDDetectorConstruction::SetWoWater(G4int answer)
{
  fWithWoWater = answer;
  WLGDDetectorConstruction::DefineMaterials();
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

void WLGDDetectorConstruction::SetMaGeMaterial(G4int answer)
{
  fMaGeMaterial = answer;
  WLGDDetectorConstruction::DefineMaterials();
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

void WLGDDetectorConstruction::DefineCommands()
{

  // Define geometry command directory using generic messenger class
  fDetectorMessenger = new G4GenericMessenger(this, "/WLGD/detector/",
                                              "Commands for controlling detector setup");

  // switch command
  fDetectorMessenger
    ->DeclareMethod("setPositionOfDetectors",
                    &WLGDDetectorConstruction::SetPositionOfDetectors)
    .SetGuidance("Set geometry model of cavern and detector")
    .SetGuidance("original = original Warwick-Legend Position")
    .SetGuidance("baseline = lowered position")
    .SetCandidates("original baseline")
    .SetStates(G4State_PreInit)
    .SetToBeBroadcasted(false);

  // switch command
  fDetectorMessenger->DeclareMethod("setGeometry", &WLGDDetectorConstruction::SetGeometry)
    .SetGuidance("Set geometry model of cavern and detector")
    .SetGuidance("baseline = NEEDS DESCRIPTION")
    .SetGuidance("baseline_smaller = Gerda cryostat with only one module")
    .SetGuidance("baseline_large_reentrance_tube = large single re-entrance tube")
    .SetGuidance("baseline_large_reentrance_tube_4m_cryo")
    .SetGuidance("alternative = NEEDS DESCRIPTION")
    .SetGuidance("hallA = Gerda mock-up for validation.")
    .SetGuidance("hallA_wo_ge = Gerda mock-up w/o HPGe")
    .SetGuidance("hallA_only_WLSR = Gerda w/o detectors but with a workaround WLSR volume")
    .SetCandidates("baseline baseline_smaller baseline_large_reentrance_tube alternative hallA hallA_wo_ge hallA_only_WLSR baseline_large_reentrance_tube_4m_cryo")
    .SetStates(G4State_PreInit)
    .SetToBeBroadcasted(false);

  // GDML Export
  fDetectorMessenger
    ->DeclareMethod("exportGeometry", &WLGDDetectorConstruction::ExportGeometry)
    .SetGuidance("Export current geometry to a GDML file")
    .SetParameterName("filename", false)
    .SetDefaultValue("wlgd.gdml")
    .SetStates(G4State_Idle)
    .SetToBeBroadcasted(false);

  // Edit: 2021/02/17 by Moritz Neuberger
  // Adding options to adjust the concentration of Xe and He-3 in the LAr to test change
  // in Ge-77 production

  // changing the concentration of Xe in the LAr
  fDetectorMessenger->DeclareMethod("XeConc", &WLGDDetectorConstruction::SetXeConc)
    .SetGuidance("Set concentration of Xe in the LAr [mg/g]")
    .SetDefaultValue("0.0")
    .SetStates(G4State_PreInit)
    .SetToBeBroadcasted(false);

  // changing the concentration of He-3 in the LAr
  fDetectorMessenger->DeclareMethod("He3Conc", &WLGDDetectorConstruction::SetHe3Conc)
    .SetGuidance("Set concentration of He3 in the LAr [mg/g]")
    .SetDefaultValue("0.0")
    .SetStates(G4State_PreInit)
    .SetToBeBroadcasted(false);

  // Edit: 2021/03/30 by Moritz Neuberger
  // Mainly for the baseline geometry

  // changing the radius if the cryostat
  fDetectorMessenger
    ->DeclareMethod("Cryostat_Radius_Outer",
                    &WLGDDetectorConstruction::SetOuterCryostatRadius)
    .SetGuidance("Set the outer radius of the cryostat [cm]")
    .SetDefaultValue("350.0")
    .SetStates(G4State_PreInit)
    .SetToBeBroadcasted(false);

  // changing the height of the cryostat
  fDetectorMessenger
    ->DeclareMethod("Cryostat_Height", &WLGDDetectorConstruction::SetCryostatHeight)
    .SetGuidance("Set the height of the cryostat [cm]")
    .SetDefaultValue("350.0")
    .SetStates(G4State_PreInit)
    .SetToBeBroadcasted(false);

  // option to remove the cupper tubes to see what impact it has on the production rate
  fDetectorMessenger
    ->DeclareMethod("Without_Cupper_Tubes",
                    &WLGDDetectorConstruction::SetWithoutCupperTubes)
    .SetGuidance("Set whether to include cupper tubes or not")
    .SetGuidance("0 = with cupper tubes")
    .SetGuidance("1 = without cupper tubes")
    .SetCandidates("0 1")
    .SetDefaultValue("0");

  // option to include borated PE in the setup (1: tubes around the re-entrance tubes, 2:
  // trubine structure)
  fDetectorMessenger
    ->DeclareMethod("With_NeutronModerators",
                    &WLGDDetectorConstruction::SetNeutronModerator)
    .SetGuidance("Set whether to include Neutron Moderators or not")
    .SetGuidance("0 = without Neutron Moderators")
    .SetGuidance("1 = with Neutron Moderators around Re-Entrance tubes")
    .SetGuidance("2 = with Neutron Moderators in turbine mode")
    .SetGuidance("3 = with Neutron Moderators in large tub")
    .SetGuidance("4 = with Neutron Moderators in turbine mode with lids")
    .SetCandidates("0 1 2 3 4")
    .SetDefaultValue("0");

  // option to include borated PE in the setup (1: tubes around the re-entrance tubes, 2:
  // trubine structure)
  fDetectorMessenger
    ->DeclareMethod("Which_Material", &WLGDDetectorConstruction::SetMaterial)
    .SetGuidance("Set which material should be used instead of Boarated PE")
    .SetGuidance("BoratedPE = normal case")
    .SetGuidance("PolyEthylene = without Boron")
    .SetGuidance("PMMA = instead using PMMA")
    .SetCandidates("BoratedPE PolyEthylene PMMA PMMA1percentB PMMA3percentB PMMA5percentB PMMA7percentB PMMA10percentB PMMA1percentGd PMMA3percentGd PMMA5percentGd PMMA7percentGd PMMA10percentGd PMMA038percentPolyGd PMMA191percentPolyGd PMMA381percentPolyGd")
    .SetDefaultValue("BoratedPE");

  // option to set the radius of the turbine structure
  fDetectorMessenger
    ->DeclareMethod("TurbineAndTube_Radius",
                    &WLGDDetectorConstruction::SetTurbineAndTubeRadius)
    .SetGuidance("Set the radius on which the borated PE pannels are aligned on [cm]")
    .SetDefaultValue("200.0")
    .SetToBeBroadcasted(false);

  // option to set the radius of the turbine structure
  fDetectorMessenger
    ->DeclareMethod("TurbineAndTube_Length",
                    &WLGDDetectorConstruction::SetTurbineAndTubeLength)
    .SetGuidance("Set the length of a panel of borated PE [cm]")
    .SetDefaultValue("50.0")
    .SetToBeBroadcasted(false);

  // option to set the radius of the turbine structure
  fDetectorMessenger
    ->DeclareMethod("TurbineAndTube_Angle",
                    &WLGDDetectorConstruction::SetTurbineAndTubeAngle)
    .SetGuidance("Set the angle on which the borated PE pannels are aligned on [deg]")
    .SetDefaultValue("45.0")
    .SetToBeBroadcasted(false);

  // option to set the radius of the turbine structure
  fDetectorMessenger
    ->DeclareMethod("TurbineAndTube_Width",
                    &WLGDDetectorConstruction::SetTurbineAndTubeWidth)
    .SetGuidance("Set the width of the borated PE pannels [cm]")
    .SetDefaultValue("5.0")
    .SetToBeBroadcasted(false);

  // option to set the radius of the turbine structure
  fDetectorMessenger
    ->DeclareMethod("TurbineAndTube_Height",
                    &WLGDDetectorConstruction::SetTurbineAndTubeHeight)
    .SetGuidance("Set the height of the borated PE pannels [cm] (this is total height)")
    .SetDefaultValue("400")
    .SetToBeBroadcasted(false);

  // option to set the radius of the turbine structure
  fDetectorMessenger
    ->DeclareMethod("TurbineAndTube_zPosition",
                    &WLGDDetectorConstruction::SetTurbineAndTubezPosition)
    .SetGuidance("Set the zPosition of the borated PE pannels [cm]")
    .SetDefaultValue("0")
    .SetToBeBroadcasted(false);

  // option to set the number of panels of the turbine structure
  fDetectorMessenger
    ->DeclareMethod("TurbineAndTube_NPanels",
                    &WLGDDetectorConstruction::SetTurbineAndTubeNPanels)
    .SetGuidance("Set the number of panels of the borated PE pannels [cm]")
    .SetDefaultValue("0")
    .SetToBeBroadcasted(false);

  // option to change the water to gadolinium weighted water in the water tank
  fDetectorMessenger
    ->DeclareMethod("With_Gd_Water", &WLGDDetectorConstruction::SetGdWater)
    .SetGuidance("Set whether to include Gd water or not")
    .SetGuidance("0 = without Gd water")
    .SetGuidance("1 = with Gd water")
    .SetCandidates("0 1")
    .SetDefaultValue("0");
  // option to change the water to gadolinium weighted water in the water tank
  fDetectorMessenger
    ->DeclareMethod("Without_Water", &WLGDDetectorConstruction::SetWoWater)
    .SetGuidance("Set whether to include water or not")
    .SetGuidance("0 = with water")
    .SetGuidance("1 = without water")
    .SetCandidates("0 1")
    .SetDefaultValue("0");


  // Define bias operator command directory using generic messenger class
  fBiasMessenger =
    new G4GenericMessenger(this, "/WLGD/bias/", "Commands for controlling bias factors");

  // switch commands
  fBiasMessenger
    ->DeclareMethod("setNeutronBias", &WLGDDetectorConstruction::SetNeutronBiasFactor)
    .SetGuidance("Set Bias factor for neutron capture process.")
    .SetDefaultValue("1.0")
    .SetStates(G4State_PreInit)
    .SetToBeBroadcasted(false);
  fBiasMessenger
    ->DeclareMethod("setMuonBias", &WLGDDetectorConstruction::SetMuonBiasFactor)
    .SetGuidance("Set Bias factor for muon nuclear process.")
    .SetDefaultValue("1.0")
    .SetStates(G4State_PreInit)
    .SetToBeBroadcasted(false);
  fBiasMessenger
    ->DeclareMethod("setNeutronYieldBias", &WLGDDetectorConstruction::SetNeutronYieldBias)
    .SetGuidance("Set Bias factor for neutron yield process.")
    .SetDefaultValue("1.0")
    .SetStates(G4State_PreInit)
    .SetToBeBroadcasted(false);

  // Define bias operator command directory using generic messenger class
  fMaterialMessenger =
    new G4GenericMessenger(this, "/WLGD/material/", "Commands for controlling the material choice");  
    
  fMaterialMessenger
    ->DeclareMethod("setMaGeMaterial", &WLGDDetectorConstruction::SetMaGeMaterial)
    .SetGuidance("Set LAr and Ge material to MaGe settings")
    .SetDefaultValue("0")
    .SetStates(G4State_PreInit)
    .SetToBeBroadcasted(false);
}//DefineCommands()
