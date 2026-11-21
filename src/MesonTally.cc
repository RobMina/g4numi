//----------------------------------------------------------------------------------------------------------------
// MesonTally.cc — see MesonTally.hh. Part of beamline-ldm-tools (MIT).
//----------------------------------------------------------------------------------------------------------------
#include "MesonTally.hh"

#include <cstdio>
#include <iostream>
#include <sstream>

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4Navigator.hh"
#include "G4ParticleDefinition.hh"
#include "G4RunManager.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4TrackingManager.hh"
#include "G4TransportationManager.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIdirectory.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VProcess.hh"

#include "TDirectory.h"
#include "TFile.h"
#include "TTree.h"

namespace ldm {

MesonTally::MesonTally(const std::string& uiDir, InfoFn info, int runDigits)
   : fInfo(std::move(info)), fRunDigits(runDigits), fEnabled(false), fFileName("mesons"), fFile(nullptr),
     fMeson(nullptr), fNavigator(nullptr), fRun(0), fNMeson(0),
     fPx(0), fPy(0), fPz(0), fE(0), fX(0), fY(0), fZ(0), fT(0), fWeight(0), fPE(0),
     fPdg(0), fTgen(0), fEvt(0), fTrk(0), fPPdg(0), fProc(0), fVol(0), fMat(0)
{
   SetPDGs("111 221 331 223");
   std::string dir = uiDir;
   if (dir.empty() || dir.back() != '/') dir += "/";
   dir += "mesons/";

   fDir = new G4UIdirectory(dir.c_str());
   fDir->SetGuidance("Neutral-meson tally for light-dark-matter studies (written at stacking, before kills)");

   fCreateCmd = new G4UIcmdWithABool((dir + "Create").c_str(), this);
   fCreateCmd->SetGuidance("Write <FileName>_<runID>.mesons.root with every selected meson (true/false)");
   fCreateCmd->SetParameterName("create", true);
   fCreateCmd->SetDefaultValue(true);
   fCreateCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

   fFileNameCmd = new G4UIcmdWithAString((dir + "FileName").c_str(), this);
   fFileNameCmd->SetGuidance("Output file name stem (default: mesons)");
   fFileNameCmd->SetParameterName("fileName", false);
   fFileNameCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

   fPDGsCmd = new G4UIcmdWithAString((dir + "PDGs").c_str(), this);
   fPDGsCmd->SetGuidance("Space-separated PDG codes to record (default: \"111 221 331 223\")");
   fPDGsCmd->SetParameterName("pdgs", false);
   fPDGsCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}

MesonTally::~MesonTally()
{
   delete fCreateCmd;
   delete fFileNameCmd;
   delete fPDGsCmd;
   delete fDir;
   delete fNavigator;
}

void MesonTally::SetPDGs(const std::string& list)
{
   fPDGs.clear();
   std::istringstream in(list);
   int code;
   while (in >> code) fPDGs.insert(code);
}

void MesonTally::SetNewValue(G4UIcommand* command, G4String newValue)
{
   if (command == fCreateCmd) fEnabled = fCreateCmd->GetNewBoolValue(newValue);
   else if (command == fFileNameCmd) fFileName = newValue;
   else if (command == fPDGsCmd) SetPDGs(newValue);
}

G4String MesonTally::GetCurrentValue(G4UIcommand* command)
{
   if (command == fCreateCmd) return fCreateCmd->ConvertToString(fEnabled);
   if (command == fFileNameCmd) return G4String(fFileName);
   if (command == fPDGsCmd) {
      std::ostringstream out;
      for (int p : fPDGs) out << p << " ";
      return G4String(out.str());
   }
   return G4String("");
}

int MesonTally::NameId(int kind, const std::string& name)
{
   auto key = std::make_pair(kind, name);
   auto it = fNames.find(key);
   if (it != fNames.end()) return it->second;
   int id = static_cast<int>(fNames.size());
   fNames[key] = id;
   return id;
}

void MesonTally::BeginRun(int runID)
{
   if (!fEnabled) return;
   fRun = runID;
   fNMeson = 0;
   fNames.clear();

   char fileName[1024];
   snprintf(fileName, sizeof(fileName), "%s_%0*d.mesons.root", fFileName.c_str(), fRunDigits, runID);
   TDirectory::TContext keepDir;   // leave the application's current ROOT directory alone
   fFile = new TFile(fileName, "RECREATE", "neutral meson tally (beamline-ldm-tools MesonTally)");
   if (fFile == nullptr || fFile->IsZombie()) {
      G4Exception("ldm::MesonTally::BeginRun", "MesonTally001", FatalException,
                  ("cannot open " + std::string(fileName)).c_str());
   }
   fMeson = new TTree("meson", "neutral mesons at stacking (before kills)");
   fMeson->Branch("px", &fPx, "px/D");
   fMeson->Branch("py", &fPy, "py/D");
   fMeson->Branch("pz", &fPz, "pz/D");
   fMeson->Branch("E", &fE, "E/D");
   fMeson->Branch("x", &fX, "x/D");
   fMeson->Branch("y", &fY, "y/D");
   fMeson->Branch("z", &fZ, "z/D");
   fMeson->Branch("t", &fT, "t/D");
   fMeson->Branch("weight", &fWeight, "weight/D");
   fMeson->Branch("pdg", &fPdg, "pdg/I");
   fMeson->Branch("tgen", &fTgen, "tgen/I");
   fMeson->Branch("evt", &fEvt, "evt/I");
   fMeson->Branch("trk", &fTrk, "trk/I");
   fMeson->Branch("ppdg", &fPPdg, "ppdg/I");
   fMeson->Branch("pE", &fPE, "pE/D");
   fMeson->Branch("proc", &fProc, "proc/I");
   fMeson->Branch("vol", &fVol, "vol/I");
   fMeson->Branch("mat", &fMat, "mat/I");

   if (fNavigator == nullptr) {
      fNavigator = new G4Navigator();   // private navigator: never move the tracking navigator
      fNavigator->SetWorldVolume(
         G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking()->GetWorldVolume());
   }
   std::cout << "ldm::MesonTally: writing " << fileName << ", PDGs";
   for (int p : fPDGs) std::cout << " " << p;
   std::cout << std::endl;
}

void MesonTally::Consider(const G4Track* aTrack)
{
   if (!fEnabled || fMeson == nullptr) return;
   const int pdg = aTrack->GetDefinition()->GetPDGEncoding();
   if (fPDGs.count(pdg) == 0) return;
   if (aTrack->GetCurrentStepNumber() != 0) return;   // re-stacked (suspended) track: already counted

   const G4ThreeVector& p = aTrack->GetMomentum();
   const G4ThreeVector& x = aTrack->GetPosition();
   fPx = p.x() / GeV; fPy = p.y() / GeV; fPz = p.z() / GeV;
   fE = aTrack->GetTotalEnergy() / GeV;
   fX = x.x() / cm; fY = x.y() / cm; fZ = x.z() / cm;
   fT = aTrack->GetGlobalTime() / ns;
   fPdg = pdg;
   fTrk = aTrack->GetTrackID();

   const Info info = fInfo ? fInfo(aTrack) : Info{1.0, -1};
   fWeight = info.weight;
   fTgen = info.tgen;

   const G4Event* evt = G4RunManager::GetRunManager()->GetCurrentEvent();
   fEvt = evt ? evt->GetEventID() : -1;

   // Secondaries are stacked right after their parent's tracking ends, while the tracking manager still holds
   // the parent (G4EventManager::DoProcessing deletes it only after StackTracks).
   // By then the parent's own energy is its post-interaction value (its rest mass if the interaction killed it), so
   // pE is the parent's total energy at the pre-step point of its last step: the energy entering the step whose
   // end created this meson (it omits only that step's continuous loss). Mesons made at an earlier step of a
   // parent that survived its interaction get pE = -1, since that step is no longer available.
   fPPdg = 0;
   fPE = -1.;
   if (aTrack->GetParentID() > 0) {
      const G4Track* parent = G4EventManager::GetEventManager()->GetTrackingManager()->GetTrack();
      if (parent != nullptr && parent->GetTrackID() == aTrack->GetParentID()) {
         fPPdg = parent->GetDefinition()->GetPDGEncoding();
         const G4Step* last = parent->GetStep();
         if (last != nullptr && last->GetPreStepPoint() != nullptr && last->GetPostStepPoint() != nullptr &&
             (last->GetPostStepPoint()->GetPosition() - aTrack->GetPosition()).mag() < 1.0 * um) {
            fPE = last->GetPreStepPoint()->GetTotalEnergy() / GeV;
         }
      }
   }

   const G4VProcess* creator = aTrack->GetCreatorProcess();
   fProc = NameId(0, creator ? std::string(creator->GetProcessName()) : std::string("primary"));

   G4VPhysicalVolume* pv = fNavigator->LocateGlobalPointAndSetup(x, nullptr, false, true);
   std::string vol("outside"), mat("none");
   if (pv != nullptr && pv->GetLogicalVolume() != nullptr) {
      vol = pv->GetLogicalVolume()->GetName();
      if (pv->GetLogicalVolume()->GetMaterial() != nullptr) mat = pv->GetLogicalVolume()->GetMaterial()->GetName();
   }
   fVol = NameId(1, vol);
   fMat = NameId(2, mat);

   fMeson->Fill();
   ++fNMeson;
}

void MesonTally::EndRun(int nEvents, const std::string& physicsList, const std::string& commit)
{
   if (!fEnabled || fFile == nullptr) return;
   TDirectory::TContext keepDir(fFile);

   int pots = nEvents, run = fRun;
   long long nmeson = fNMeson;
   char physlist[256], version[256];
   snprintf(physlist, sizeof(physlist), "%s", physicsList.empty() ? "unknown" : physicsList.c_str());
   snprintf(version, sizeof(version), "%s", commit.empty() ? "unknown" : commit.c_str());
   TTree meta("meta", "meson tally run metadata");
   meta.Branch("pots", &pots, "pots/I");
   meta.Branch("run", &run, "run/I");
   meta.Branch("nmeson", &nmeson, "nmeson/L");
   meta.Branch("physlist", physlist, "physlist/C");
   meta.Branch("commit", version, "commit/C");
   meta.Fill();

   int id = 0, kind = 0;
   char name[256];
   TTree names("names", "id legend for meson proc / vol / mat (kind 0 / 1 / 2)");
   names.Branch("id", &id, "id/I");
   names.Branch("kind", &kind, "kind/I");
   names.Branch("name", name, "name/C");
   for (const auto& kv : fNames) {
      kind = kv.first.first;
      snprintf(name, sizeof(name), "%s", kv.first.second.c_str());
      id = kv.second;
      names.Fill();
   }

   fMeson->Write();
   meta.Write();
   names.Write();
   fFile->Close();
   delete fFile;   // also deletes fMeson
   fFile = nullptr;
   fMeson = nullptr;
   std::cout << "ldm::MesonTally: " << nmeson << " mesons from " << pots << " POT written" << std::endl;
}

}  // namespace ldm
