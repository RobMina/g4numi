//----------------------------------------------------------------------------------------------------------------
// MesonTally.hh — portable light-dark-matter meson tally for Geant4 beamline simulations
// Part of beamline-ldm-tools (MIT). Vendor MesonTally.hh / MesonTally.cc into the application unchanged.
//
// Records every particle of the selected PDG codes (default pi0, eta, eta', omega) when it is stacked: call
// Consider(track) first thing in the application's G4UserStackingAction::ClassifyNewTrack, before any kill
// (kinetic-energy threshold, importance-weight Russian roulette, EM kill). Each entry carries the weight and
// generation the application assigns the new track (normally the parent's importance weight), so the weighted
// list is unbiased whatever the kills do afterwards.
//
// Output <stem>_<runID>.mesons.root (run number zero-padded to `runDigits`):
//   meson : px py pz E [GeV], x y z [cm, global], t [ns], weight, pdg, tgen, evt, trk, ppdg, pE [GeV, parent's
//           total energy entering the step that made the meson; -1 if unknown],
//           proc, vol, mat (ids into names). The first 11 match the g4numi `meson` tree.
//   meta  : pots (events in the run), run, nmeson, physlist, commit
//   names : id, kind (0 creator process, 1 logical volume, 2 material), name
//
// UI commands under <uiDir>/mesons/: Create (off by default), FileName, PDGs.
//
// Integration (see README.md):
//   static ldm::MesonTally tally("/LBNE", [](const G4Track* t) {           // construct before the macro runs
//       auto* i = dynamic_cast<const MyTrackInformation*>(t->GetUserInformation());
//       return ldm::MesonTally::Info{i ? i->GetNImpWt() : 1.0, i ? i->GetTgen() : -1}; });
//   BeginOfRunAction: tally.BeginRun(run->GetRunID());
//   EndOfRunAction:   tally.EndRun(run->GetNumberOfEvent(), physicsListName, versionString);
//   ClassifyNewTrack: tally.Consider(track);                                // before any kill
//----------------------------------------------------------------------------------------------------------------
#ifndef LDM_MesonTally_h
#define LDM_MesonTally_h 1

#include <functional>
#include <map>
#include <set>
#include <string>
#include <utility>

#include "G4UImessenger.hh"
#include "globals.hh"

class G4Track;
class G4Navigator;
class G4UIdirectory;
class G4UIcmdWithABool;
class G4UIcmdWithAString;
class TFile;
class TTree;

namespace ldm {

class MesonTally : public G4UImessenger
{
public:
   struct Info { double weight; int tgen; };
   using InfoFn = std::function<Info(const G4Track*)>;

   MesonTally(const std::string& uiDir, InfoFn info, int runDigits = 5);
   ~MesonTally() override;

   void BeginRun(int runID);
   void EndRun(int nEvents, const std::string& physicsList = "", const std::string& commit = "");
   void Consider(const G4Track* aTrack);

   bool IsEnabled() const { return fEnabled; }
   void SetEnabled(bool on) { fEnabled = on; }
   void SetFileStem(const std::string& stem) { fFileName = stem; }
   void SetPDGs(const std::string& list);

   void SetNewValue(G4UIcommand* command, G4String newValue) override;
   G4String GetCurrentValue(G4UIcommand* command) override;

private:
   int NameId(int kind, const std::string& name);

   InfoFn fInfo;
   int fRunDigits;
   bool fEnabled;
   std::string fFileName;
   std::set<int> fPDGs;

   TFile* fFile;
   TTree* fMeson;
   G4Navigator* fNavigator;
   int fRun;
   long fNMeson;
   std::map<std::pair<int, std::string>, int> fNames;

   double fPx, fPy, fPz, fE, fX, fY, fZ, fT, fWeight, fPE;
   int fPdg, fTgen, fEvt, fTrk, fPPdg, fProc, fVol, fMat;

   G4UIdirectory* fDir;
   G4UIcmdWithABool* fCreateCmd;
   G4UIcmdWithAString* fFileNameCmd;
   G4UIcmdWithAString* fPDGsCmd;
};

}  // namespace ldm

#endif
