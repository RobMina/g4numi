//
// NumiAnalysis.hh
//
// Modified Jul 2005 by A. MArino to make data_t and hadmmtuple_t classes

#ifndef NUMIANALYSIS_HH
#define NUMIANALYSIS_HH

#include "globals.hh"
#include "NumiTrajectory.hh"
#include "NumiDataInput.hh"

//root
#include "TSystem.h"

//G4
#include "G4ios.hh"
#include "G4TrajectoryContainer.hh"

#include <map>
#include <vector>

class G4ParticleDefinition;
class G4Step;
class G4StepPoint;
class G4Material;
class G4VPhysicalVolume;
class TFile;
class TTree;
class G4Track;
class G4VTrajectory;
class data_t;
class hadmmtuple_t;
class draytupleMIB_t;
class draytupleSPB_t;
class absbkgtuple_t;
class zptuple_t;
class target_exit_t;

//New for DK2NU:///
#include "dk2nu/tree/dkmeta.h"
#include "dk2nu/tree/dk2nu.h"
//class bsim::Location;
//class bsim::DkMeta;
//class bsim::NuRay;
//class bsim::Decay;
//class bsim::Ancestor;
//class bsim::TgtExit;
//class bsim::Traj;
//class bsim::Dk2Nu;
///////////////////

class NumiAnalysis
{
public:

   NumiAnalysis();
   ~NumiAnalysis();
   
   void book();
   void finish();
   void finishPL();
   void FillMeta(); //new function to fill dk2meta
   void FillNeutrinoNtuple(const G4Track& track,
                           const std::vector<G4VTrajectory*>& nuHistory);
   G4int GetNucleus(G4String nucl_name); //new function to get target element
   void FillHadmmNtuple(const G4Track& track,Int_t mm_num,Int_t cellNum);
   void FillHadmmNtuple();
   void FillBXDRAW(const G4Step* aStep);
   void FillAlcEdepInfo(const G4Track& track, G4int alc);
   void FillEdep(G4double edep, const G4ParticleDefinition* particleDefinition,
                 const G4int alc, const G4int IntExt, const G4double zpos,
                 Int_t cellNum, const Int_t trackID, const G4double weight);
   void FillAbsorberBkgrdNtuple(const G4Track& track);
   void WriteHadmmNtuple();
   
   void FillZpNtuple();
   void FillZpNtuple(const G4Track& track, Int_t zpnum);
   void WriteZpNtuple();
   
   void FillTarNtuple();
   void FillTarNtuple(const G4Track& track);
   void WriteTarNtuple();
    
   void FillMeson(const G4Track& track);

   // Mesons crossing the outer surface of the target assembly, i.e. the
   // flux that actually leaves the target.  Same variables as the meson
   // tree, but recorded on exit instead of at creation, so the two can be
   // compared directly (see dm_production_output/g4numi/comparisons/mipp).
   void FillTgtExit(const G4Track& track);

   // Provenance bookkeeping.  RecordTrack() is called for *every* track at the
   // start of its tracking, so that when a meson shows up we can walk its
   // ancestry back to the beam proton without keeping trajectories around.
   // Cheap: one small struct per track, cleared at the start of each event.
   void RecordTrack(const G4Track& track);
   void FillMesonProvenance(const G4Track& track);

   // EM cascade seeds (/NuMI/output/outputEMSeeds).  A seed is a gamma or
   // e+- whose creator is not itself an EM-cascade process (see
   // IsEMCascadeProcess), i.e. it comes from a pi0/eta decay, hadronic
   // de-excitation, a delta ray of a hadron or muon, a Dalitz pair, ...
   //  - RecordTrackPdg: every tracked particle, so a seed can report its
   //    parent's PDG code (called from PreUserTrackingAction).
   //  - FillEMSeedAtCreation: e+- seeds, recorded from NumiStackingAction at
   //    birth (they are then killed).
   //  - CacheEMSeedCreation / FillEMSeedAtInteraction: photon seeds, recorded
   //    by NumiSteppingAction at their first interaction (then killed).
   static G4bool IsEMCascadeProcess(const G4String& procName);
   void RecordTrackPdg(const G4Track& track);
   void FillEMSeedAtCreation(const G4Track& track);
   void CacheEMSeedCreation(const G4Track& track, const G4StepPoint* pre);
   void FillEMSeedAtInteraction(const G4Track& track, const G4Step* step,
                                const G4String& procName);
   // pi0/eta that NumiStackingAction kills before they reach the meson tree
   // (diagnostic; summary printed at end of run when emseed is on).
   void CountKilledMeson(const G4Track& track, G4int reason);

   NumiTrajectory* GetParentTrajectory(G4int parentID);
   static NumiAnalysis* getInstance();
   
   void SetAlcEdepFlag(G4bool AlcEdep);
   void SetMuInAlcFlag(G4bool MuInAlc);
   
   void SetCount(G4int count);
   G4int GetCount();
   void SetEntry(G4int entry);
   G4int GetEntry();
   G4double GetDistanceInVolume(NumiTrajectory* wanted_traj, G4String wanted_vol);
private:
   std::string GetOFileName(std::string ifilename);
   
private:
   static NumiAnalysis* instance;
   
   G4double x;
   G4double y;
   G4double z;

   G4double noProtons;
   char asciiFileName[50], nuNtupleFileName[50], hadmmNtupleFileName[50], tarNtupleFileName[50];
   char zpNtupleFileName[50], bxdrawFileName[50];
   std::map<int, int> code;
   
   NumiDataInput* NumiData;
   
   TFile* hadmmNtuple;
   TFile* absbkgNtuple;
   TFile* nuNtuple;
   TFile* zpNtuple;
   TFile* tarNtuple;
   
   TTree* tree;
   TTree* hadmmtree;
   TTree* absbkgtree;
   TTree* zptree;
   TTree* zpmetatree;
   TTree* tartree;
   TTree* mesontree;
   Double_t meson_px;
   Double_t meson_py;
   Double_t meson_pz;
   Double_t meson_E;
   Double_t meson_x;
   Double_t meson_y;
   Double_t meson_z;
   Double_t meson_t;
   Double_t meson_weight;
   Int_t meson_pdg;
   Int_t meson_tgen;

   // --- meson provenance (filled when NumiData->createMesonProvenance) -------
   // this meson
   Int_t    meson_proc;      // creator process id (see the provnames tree)
   Int_t    meson_vol;       // volume id at the production vertex
   Int_t    meson_intgt;     // 1 if produced inside the target mother volume
   // its parent
   Int_t    meson_ppdg;
   Int_t    meson_pproc;
   Int_t    meson_pvol;
   Int_t    meson_pintgt;
   Double_t meson_pE;        // parent total energy at the parent's own birth
   Double_t meson_pvx, meson_pvy, meson_pvz;
   // the ancestor that carried the energy out of the target: the last one born
   // inside the target whose next-in-chain daughter was born outside it.
   // meson_cflag: 0 = meson itself born in the target (nothing left it),
   //              1 = carrier is that last in-target ancestor,
   //              2 = no ancestor was born in the target, so the carrier is the
   //                  beam proton itself and meson_i* is where it first
   //                  interacted along this chain.
   Int_t    meson_cflag;
   Int_t    meson_cpdg;
   Int_t    meson_cproc;
   Int_t    meson_cgen;      // generations between the carrier and this meson
   Double_t meson_cE;
   Double_t meson_cvx, meson_cvy, meson_cvz;
   Double_t meson_ivx, meson_ivy, meson_ivz;  // where the chain left the target

   // id <-> name legend for meson_proc / meson_vol, written as its own tree
   TTree* provnamestree;
   Int_t  provname_id;
   Int_t  provname_kind;     // 0 = process, 1 = volume
   char   provname_name[80];
   std::map<G4String, G4int> fProcIds;
   std::map<G4String, G4int> fVolIds;
   G4int ProvId(std::map<G4String, G4int>& ids, const G4String& name);

   struct TrackProv {
      G4int    pdg;
      G4int    parentID;
      G4int    proc;
      G4int    vol;
      G4int    tgen;
      G4bool   inTgt;
      G4double E;
      G4double x, y, z;
   };
   std::map<G4int, TrackProv> fTrackProv;
   G4bool fTgtMotherResolved;
   const G4VPhysicalVolume* fTgtMotherPV;

   // --- emseed tree ----------------------------------------------------------
   TTree* emseedtree;
   Int_t    es_evt, es_trackID, es_parentID, es_pdg, es_parent_pdg, es_cproc, es_tgen;
   Double_t es_nimpwt;
   // at creation
   Double_t es_E, es_px, es_py, es_pz, es_x, es_y, es_z, es_t;
   Int_t    es_vol, es_mat;
   // at the first interaction (e+- seeds: copy of the creation point, iproc 0)
   Int_t    es_iproc, es_ivol, es_imat;
   Double_t es_ix, es_iy, es_iz, es_it, es_iE;
   TTree* emseednamestree;
   Int_t  esname_id;
   Int_t  esname_kind;       // 0 = process, 1 = volume, 2 = material
   char   esname_name[80];
   std::map<G4String, G4int> fESProcIds;
   std::map<G4String, G4int> fESVolIds;
   std::map<G4String, G4int> fESMatIds;
   std::vector<G4int> fESTrackPdg;   // trackID -> PDG, current event
   G4int fESTrackPdgEvt;
   G4int ESVolId(const G4VPhysicalVolume* pv);
   G4int ESMatId(const G4Material* mat);
   void  FillEMSeedCommon(const G4Track& track);
   // cache of the photon currently being tracked (creation point)
   G4int    fESCacheTrackID, fESCacheEvt;
   Double_t fESCacheE, fESCachePx, fESCachePy, fESCachePz;
   Double_t fESCacheX, fESCacheY, fESCacheZ, fESCacheT;
   G4int    fESCacheVol, fESCacheMat;
   // counters
   G4long   fESNPhoton, fESNLepton, fESNNoCache;
   // killed pi0/eta: [species 0 pi0, 1 eta][reason 0 KE<thr, 1 pz<0]
   G4long   fKilledN[2][2];
   G4double fKilledW[2][2];
   G4double fKilledWE[2][2];

   TTree* tgtexittree;
   Double_t tgtexit_px;
   Double_t tgtexit_py;
   Double_t tgtexit_pz;
   Double_t tgtexit_E;
   Double_t tgtexit_x;
   Double_t tgtexit_y;
   Double_t tgtexit_z;
   Double_t tgtexit_t;
   Double_t tgtexit_weight;
   Int_t tgtexit_pdg;
   Int_t tgtexit_tgen;

   data_t *g4data;
   hadmmtuple_t *g4hmmdata;
   draytupleMIB_t *g4draydataMIB;
   draytupleSPB_t *g4draydataSPB;
   absbkgtuple_t *g4absbkgdata;
   zptuple_t* g4zpdata; 
   target_exit_t *g4tardata;
   
   G4int fcount;
   G4int fentry;
   G4int fZpNPots;
   G4int fZpNGeneral;
   G4int fZpNMuon;
   
   std::vector<bool> fAlcEdep_called;
   std::vector<bool> fMuInAlc;
   
   //New fr DK2NU:////////////////////////
  std::vector<bsim::Location> vec_loc;
  bsim::DkMeta*  this_meta;
  TTree* meta;

  bsim::Dk2Nu* this_dk2nu;
  bsim::Decay this_decay; 
  bsim::TgtExit this_tgtexit; 
  std::vector<bsim::Traj> vec_traj;
  std::vector<bsim::NuRay> vec_nuray;   
  std::vector<bsim::Ancestor> vec_ancestor;
  std::vector<int> vec_int;
  std::vector<double> vec_dbl;
  //////////////////////////////////////////
  // September 2017 Paul Lebrun.. 
  // Study of the inclusive  numu spectrum for various magnetic field 
  // To speed the I/O and processing on the farm, we include here a fixed size
  // weighted histogram 
  double energyBinSimpleHistoMinerva; // depend on the LE vs ME configuration.. 
  std::vector<double> MinervaNuMuHisto;
  std::vector<double> MinervaNuMuBarHisto;
  std::vector<double> NovaFarNuMuHisto;
  std::vector<double> NovaNearNuMuHisto;
  std::vector<double> NovaFarNuMuBarHisto;
  std::vector<double> NovaNearNuMuBarHisto;
};

#endif 
