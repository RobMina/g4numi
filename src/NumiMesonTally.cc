//----------------------------------------------------------------------
// NumiMesonTally.cc -- see include/NumiMesonTally.hh
//----------------------------------------------------------------------
#include "NumiMesonTally.hh"

#include "G4Track.hh"

#include "NumiTrackInformation.hh"

// Version string written to the tally's meta tree; the build may pass -DG4NUMI_COMMIT="\"<sha>\"".
#ifndef G4NUMI_COMMIT
#define G4NUMI_COMMIT "g4numi (commit not recorded at build)"
#endif

namespace NumiMesonTally {

namespace {
std::string& PhysicsListName()
{
   static std::string name;
   return name;
}
}  // namespace

ldm::MesonTally* Instance()
{
   // Secondaries carry the parent's importance weight and tgen = parent + 1 (NumiTrackingAction).
   // The primary has no NumiTrackInformation yet when it is stacked (PreUserTrackingAction adds it).
   static ldm::MesonTally* tally = new ldm::MesonTally(
      "/NuMI",
      [](const G4Track* t) {
         const auto* info = dynamic_cast<const NumiTrackInformation*>(t->GetUserInformation());
         if (info != nullptr) return ldm::MesonTally::Info{info->GetNImpWt(), info->GetTgen()};
         return ldm::MesonTally::Info{1.0, t->GetParentID() == 0 ? 1 : -1};
      },
      4);
   return tally;
}

void SetPhysicsListName(const std::string& name) { PhysicsListName() = name; }

void EndRun(int nEvents) { Instance()->EndRun(nEvents, PhysicsListName(), G4NUMI_COMMIT); }

}  // namespace NumiMesonTally
