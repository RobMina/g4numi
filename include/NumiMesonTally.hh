//----------------------------------------------------------------------
// NumiMesonTally.hh -- g4numi glue for the LDM neutral-meson tally.
//
// The tally itself is ldm::MesonTally (include/MesonTally.hh, src/MesonTally.cc), vendored unchanged from
// beamline-ldm-tools geant4/MesonTally. This file only adapts it to g4numi: UI directory /NuMI (commands
// /NuMI/mesons/Create, FileName, PDGs), run number padded to 4 digits like the other g4numi outputs, weight and
// tgen from NumiTrackInformation. Hooks: NumiRunAction (constructor, Begin/EndOfRunAction),
// NumiStackingAction::ClassifyNewTrack (first statement, before any kill), and g4numi.cc (physics-list name).
//
// Unlike the legacy `meson` tree in the zpNtuple file (NumiAnalysis::FillMeson, filled at track start in
// NumiTrackingAction::PreUserTrackingAction), this tally also records the mesons the stacking action kills
// (pz < 0, kinetic energy below KillTrackingThreshold, importance-weight roulette). Each entry carries the
// parent's importance weight, which NumiTrackingAction::PostUserTrackingAction copies onto every secondary
// before G4EventManager stacks them.
//----------------------------------------------------------------------
#ifndef NumiMesonTally_h
#define NumiMesonTally_h 1

#include <string>

#include "MesonTally.hh"

namespace NumiMesonTally {
ldm::MesonTally* Instance();
void SetPhysicsListName(const std::string& name);
void EndRun(int nEvents);
}

#endif
