// -- this is as far as I can see a copy of the tutorial how to use the Bias Change Cross-Section from the G4 tutorials
// -- for a better explanation go there

#include "WLGDBiasMultiParticleChangeCrossSection.hh"


WLGDBiasMultiParticleChangeCrossSection::WLGDBiasMultiParticleChangeCrossSection()
: G4VBiasingOperator("TestManyExponentialTransform")
{}

WLGDBiasMultiParticleChangeCrossSection::~WLGDBiasMultiParticleChangeCrossSection()
{
  for(auto& it : fBOptrForParticle)
  {
    delete it.second;
  }
}

void WLGDBiasMultiParticleChangeCrossSection::AddParticle(const G4String& particleName)
{
  const G4ParticleDefinition* particle =
    G4ParticleTable::GetParticleTable()->FindParticle(particleName);

  if(particle == nullptr)
  {
    G4ExceptionDescription ed;
    ed << "Particle `" << particleName << "' not found !" << G4endl;
    G4Exception("WLGDBiasMultiParticleChangeCrossSection::AddParticle(...)", "exWLGD.02",
                JustWarning, ed);
    return;
  }

  WLGDBiasChangeCrossSection* optr = new WLGDBiasChangeCrossSection(particleName);
  optr->SetNeutronFactor(fNeutronBias);
  optr->SetMuonFactor(fMuonBias);
  optr->SetNeutronYieldFactor(fNeutronYieldBias);
  G4cout << " >>> MultiBias: set neutron and muon factors to " << fNeutronBias << ", "
         << fMuonBias << ", " << fNeutronYieldBias << G4endl;
  fParticlesToBias.push_back(particle);
  fBOptrForParticle[particle] = optr;
}

G4VBiasingOperation*
WLGDBiasMultiParticleChangeCrossSection::ProposeOccurenceBiasingOperation(
  const G4Track* track, const G4BiasingProcessInterface* callingProcess){
  if(fCurrentOperator != nullptr)
  {
    return fCurrentOperator->GetProposedOccurenceBiasingOperation(track, callingProcess);
  }

  return nullptr;
}

void WLGDBiasMultiParticleChangeCrossSection::StartTracking(const G4Track* track){
  // -- fetch the underneath biasing operator, if any, for the current particle type:
  const G4ParticleDefinition* definition = track->GetParticleDefinition();
  auto                        it         = fBOptrForParticle.find(definition);
  fCurrentOperator                       = nullptr;
  if(it != fBOptrForParticle.end())
  {
    fCurrentOperator = (*it).second;
  }
}

void WLGDBiasMultiParticleChangeCrossSection::OperationApplied(
  const G4BiasingProcessInterface* callingProcess, G4BiasingAppliedCase biasingCase,
  G4VBiasingOperation* occurenceOperationApplied, G4double weightForOccurenceInteraction,
  G4VBiasingOperation*     finalStateOperationApplied,
  const G4VParticleChange* particleChangeProduced){
  // -- inform the underneath biasing operator that a biased interaction occured:
  if(fCurrentOperator != nullptr)
  {
    fCurrentOperator->ReportOperationApplied(
      callingProcess, biasingCase, occurenceOperationApplied,
      weightForOccurenceInteraction, finalStateOperationApplied, particleChangeProduced);
  }
}
