 //----------------------------------------------------------------------
// NumiAnalysis.cc
//
// $Id: NumiAnalysis.cc,v 1.26.4.28 2020/09/13 20:35:59 laliaga Exp $
//----------------------------------------------------------------------

#include <vector>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <stdlib.h>
#include <typeinfo>
#include <cstring>
#if defined(__clang__) || defined(__GNUC__)
  // C++ name mangling
  #include <cxxabi.h>
#endif

//Root
#include "TSystem.h"        // ROOT head file for a generic interface to the OS
#include "TStopwatch.h"     // ROOT head file for stopwatch for real and cpu time
#include "TFile.h"
#include "TTree.h"

//GEANT4
#include "globals.hh"
#include "G4ios.hh"
#include "G4Track.hh"
#include "G4SteppingManager.hh"
#include "G4ThreeVector.hh"
#include "G4TrajectoryContainer.hh"
#include "G4RunManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4Navigator.hh"
#include "G4TransportationManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VTouchable.hh"
#include "G4VProcess.hh"
#include "G4Run.hh"

#include "G4VUserPhysicsList.hh"
#include "G4Version.hh"
#include "G4Material.hh"
#include "G4LogicalVolume.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4Event.hh"


#include "CLHEP/Units/PhysicalConstants.h"

//g4numi
#include "dataProducts/data_t.hh"
#include "dataProducts/hadmmtuple_t.hh"
#include "dataProducts/draytupleMIB_t.hh"
#include "dataProducts/draytupleSPB_t.hh"
#include "dataProducts/absbkgtuple_t.hh"
#include "dataProducts/zptuple_t.hh" // for raytracing
#include "dataProducts/target_exit_t.hh"// for hadron production studies
#include "NumiParticleCode.hh"
#include "NumiAnalysis.hh"
#include "NumiTrackInformation.hh"
#include "NumiPrimaryGeneratorAction.hh"
#include "NumiDataInput.hh"
#include "NumiNuWeight.hh"
#include "NumiTrajectory.hh"
#include "NumiRunManager.hh"
#include "NumiSteppingAction.hh"

//NEW For DK2NU
#include <string>
#include <sstream>
#include <iostream>

using namespace std;

NumiAnalysis* NumiAnalysis::instance = 0;

//------------------------------------------------------------------------------------
NumiAnalysis::NumiAnalysis()
   :g4hmmdata(0),
    g4draydataMIB(0),
    g4draydataSPB(0),
    g4absbkgdata(0),
    energyBinSimpleHistoMinerva(-1.0e10)
// The above histo must be of the same size..
{
  //
  //  Simple histogram for quick study of the flux spectrum, Paul Lebrun, September/October 2017.
  // To turn on, simply uncomment the bin size.
  //
//  energyBinSimpleHistoMinerva = 0.1;
  if (energyBinSimpleHistoMinerva > 0.) {
    MinervaNuMuHisto = std::vector<double>(250,0.);
    MinervaNuMuBarHisto = std::vector<double>(250,0.);
    NovaNearNuMuHisto = std::vector<double>(100,0.);
    NovaNearNuMuBarHisto = std::vector<double>(100,0.);
    NovaFarNuMuHisto = std::vector<double>(100,0.);
    NovaFarNuMuBarHisto = std::vector<double>(100,0.);
  }

  NumiData = NumiDataInput::GetNumiDataInput();

  if(NumiData->GetDebugLevel() > 0)
  {
     std::cout << "NumiAnalysis Constructor Called." << std::endl;
  }

#ifdef G4ANALYSIS_USE
#endif


   G4cout << "NumiAnalysis" << G4endl;

  g4data = new data_t();
  g4zpdata = new zptuple_t();
  g4tardata = new target_exit_t();

  //New for DK2NU:
  this_meta  = new bsim::DkMeta();
  this_dk2nu = new bsim::Dk2Nu();

  fcount = 0;
  fentry = 0;
  for(int i = 0; i < 3; ++i)
  {
     fAlcEdep_called.push_back(false);
     fMuInAlc.push_back(false);
  }

  code[-13]   = 10;
  code[13]    = 11;
  code[111]   = 23;
  code[211]   = 13;
  code[-211]  = 14;
  code[130]   = 12;
  code[321]   = 15;
  code[-321]  = 16;
  code[2112]  = 8;
  code[2212]  = 1;
  code[-2212] = 2;
  code[310]   = 19;
  code[3122]  = 17;
  code[3222]  = 21;
  code[3212]  = 22;
  code[3112]  = 20;
  mesontree = 0;
  tgtexittree = 0;
  provnamestree = 0;
  emseedtree = 0;
  emseednamestree = 0;
  fESTrackPdgEvt = -1;
  fESCacheTrackID = -1;
  fESCacheEvt = -1;
  fESNPhoton = fESNLepton = fESNNoCache = 0;
  for(int i = 0; i < 2; ++i) for(int j = 0; j < 2; ++j){
     fKilledN[i][j] = 0; fKilledW[i][j] = 0.; fKilledWE[i][j] = 0.;
  }
  fTgtMotherResolved = false;
  fTgtMotherPV = 0;
}
//------------------------------------------------------------------------------------
NumiAnalysis::~NumiAnalysis()
{
   if(NumiData->GetDebugLevel() > 0)
   {
      std::cout << "NumiAnalysis Destructor Called." << std::endl;
   }

#ifdef G4ANALYSIS_USE
  // delete things
#endif
}
//------------------------------------------------------------------------------------
NumiAnalysis* NumiAnalysis::getInstance()
{
  if (instance == 0) instance = new NumiAnalysis;
  return instance;
}
//------------------------------------------------------------------------------------
void NumiAnalysis::book()
{
  //
  // Set the Minerva flux binning depending on the horn configuration. (le vs me )
  //
  if (NumiData->GetHornConfig() == G4String("le")) {
    std::cerr << " Low energy config, Setting energyBinSimpleHistoMinerva to 0.02.. " << std::endl;
    energyBinSimpleHistoMinerva = 0.02;
  }

   if(NumiData->GetDebugLevel() == 10) { G4cout << "NumiAnalysis::book() called." << G4endl;}

  G4RunManager* pRunManager = G4RunManager::GetRunManager();
  if (NumiData->createNuNtuple){
    sprintf(nuNtupleFileName,"%s_%04d%s.root",(NumiData->nuNtupleName).c_str(),pRunManager->GetCurrentRun()->GetRunID(), (NumiData->geometry).c_str());
    nuNtuple = new TFile(nuNtupleFileName,"RECREATE","root ntuple");
    G4cout << "Creating neutrino ntuple: "<<nuNtupleFileName<<G4endl;
    /*
    tree = new TTree("nudata","g4numi Neutrino ntuple");
    tree->Branch("data","data_t",&g4data,32000,1);
    */
    //NEW DK2NU:
    tree = new TTree("dk2nuTree","g4numi Neutrino ntuple");
    tree->Branch("dk2nu","bsim::Dk2Nu",&this_dk2nu,32000,99);
    meta = new TTree("dkmetaTree","Meta info for neutrino run");
    meta->Branch("dkmeta","bsim::DkMeta",&this_meta,32000,99);

  }

  if (NumiData->createTarNtuple){
      sprintf(tarNtupleFileName,"%s_%04d%s.root",(NumiData->tarNtupleName).c_str(),pRunManager->GetCurrentRun()->GetRunID(), (NumiData->geometry).c_str());
      tarNtuple = new TFile(tarNtupleFileName,"RECREATE","root ntuple");
    if(tarNtuple){
      G4cout << "Creating target ntuple: "<<tarNtupleFileName<<G4endl;
      tarNtuple->cd();
      tartree = new TTree("tardata","g4numi Neutrino ntuple");
      //tartree->BranchOld("exit_tar","exit_tar_t",&g4tardata,32000,1);
      tartree->Branch("target_exit","target_exit_t",&g4tardata,32000,99);
    }
    else{
      G4Exception("NumiAnalysis","NumiAnalysis",FatalException,"Something went wrong with tarNtuple");
    }
  }

  if(NumiData->GetSimAbsBkg() && NumiData->GetCreateAbsBkgNtuple())
  {
     std::string FileNameStr;
     if((NumiData->GetAbsBkgNtupleName()).empty())
     {
        G4cout << "NumiAnalysis::book() - Output file name not set" << G4endl;
        FileNameStr = NumiData->GetAbsBkgNtupleDir() + "/AbsBkg.root";
     }
     else FileNameStr = NumiData->GetAbsBkgNtupleDir() + "/" + NumiData->GetAbsBkgNtupleName();

     G4cout << "Creating absorber backgrounds ntuple: "<< FileNameStr.c_str()<<G4endl;

     g4absbkgdata = new absbkgtuple_t();
     absbkgNtuple = new TFile(FileNameStr.c_str(), "RECREATE","abs bkg ntuple");
     if(!absbkgNtuple) G4cout << "Failed to create file" << G4endl;
     absbkgtree = new TTree("absbkg","g4numi Absorber Backgrounds ntuple");
     absbkgtree->Branch("AbsBkgData","absbkgtuple_t",&g4absbkgdata,32000,1);
     g4absbkgdata->Clear();
  }

//////////////////////////////////////////////////////////////////
  /*
  if (NumiData->createHadmmNtuple && NumiData->useMuonInput)
  {
     std::string hadmmNtupleFileNameStr;
     if((NumiData->GetHadmmNtupleName()).empty())
        hadmmNtupleFileNameStr = GetOFileName(NumiData->GetExtNtupleFileName());
     else hadmmNtupleFileNameStr = NumiData->GetHadmmNtupleDir() + "/" + NumiData->GetHadmmNtupleName();

     G4cout << "Creating hadron and muon monitors ntuple: "<<hadmmNtupleFileNameStr.c_str()<<G4endl;

     g4hmmdata = new hadmmtuple_t();
     //
     //Create ntuple for mu just getting into mons
     //

     hadmmNtuple = new TFile(hadmmNtupleFileNameStr.c_str(), "RECREATE","hadmm ntuple");
     if(!hadmmNtuple) G4cout << "Failed to create file" << G4endl;
     hadmmtree = new TTree("hadmm","g4numi Hadron and muon monitor ntuple");
     hadmmtree->Branch("hadmmdata","hadmmtuple_t",&g4hmmdata,32000,1);
     g4hmmdata->Clear();


  }
  */
  /////////////////////////////////////////////////

  if (NumiData->createHadmmNtuple && NumiData->useMuonInput)
  {
     std::string hadmmNtupleFileNameStr;
     if((NumiData->GetHadmmNtupleName()).empty())
        hadmmNtupleFileNameStr = GetOFileName(NumiData->GetExtNtupleFileName());
     else hadmmNtupleFileNameStr = NumiData->GetHadmmNtupleDir() + "/" + NumiData->GetHadmmNtupleName();

     G4cout << "Creating hadron and muon monitors ntuple: "<<hadmmNtupleFileNameStr.c_str()<<G4endl;

     if(!NumiData->GetSimDRays())
     {
        g4hmmdata = new hadmmtuple_t();
        //
        //Create ntuple for mu just getting into mons
        //

        hadmmNtuple = new TFile(hadmmNtupleFileNameStr.c_str(), "RECREATE","hadmm ntuple");
        if(!hadmmNtuple) G4cout << "Failed to create file" << G4endl;
        hadmmtree = new TTree("hadmm","g4numi Hadron and muon monitor ntuple");
        hadmmtree->Branch("hadmmdata","hadmmtuple_t",&g4hmmdata,32000,1);
        g4hmmdata->Clear();

     }
     else if(NumiData->GetSimDRays())
     {
        g4draydataMIB = new draytupleMIB_t();

        hadmmNtuple = new TFile(hadmmNtupleFileNameStr.c_str(), "RECREATE","dray ntuple");
        if(!hadmmNtuple) G4cout << "Failed to create file" << G4endl;
        hadmmtree = new TTree("dRay","g4numi Hadron and muon monitor ntuple");
        hadmmtree->Branch("dRaydata","draytupleMIB_t",&g4draydataMIB,32000,1);
        g4draydataMIB->Clear();
     }
  }
  else if (NumiData->createHadmmNtuple && !(NumiData->useMuonInput))
  {
     std::string hadmmNtupleFileNameStr;
     if((NumiData->GetHadmmNtupleName()).empty())
     {
        stringstream path;
        path << NumiData->GetHadmmNtupleDir()
             << "/hadmmNtuple_P"
             << (int)NumiData->GetMuonBeamMomentum()
             << "_"
             << pRunManager->GetCurrentRun()->GetRunID()
             << ".root";

        hadmmNtupleFileNameStr = path.str();
     }
     else hadmmNtupleFileNameStr = NumiData->GetHadmmNtupleDir() + "/" + NumiData->GetHadmmNtupleName();

     G4cout << "Creating hadron and muon monitors ntuple: "<<hadmmNtupleFileNameStr.c_str()<<G4endl;

     if(NumiData->GetSimDRays())
     {
        g4draydataSPB = new draytupleSPB_t();

        hadmmNtuple = new TFile(hadmmNtupleFileNameStr.c_str(), "RECREATE","dray ntuple");
        if(!hadmmNtuple) G4cout << "Failed to create file" << G4endl;
        hadmmtree = new TTree("dRay","g4numi Hadron and muon monitor ntuple");
        hadmmtree->Branch("dRaydata","draytupleSPB_t",&g4draydataSPB,32000,1);
        g4draydataSPB->Clear();
     }
  }
///////////////////////////////////////////////////////

  if (NumiData->createASCII) {
    G4RunManager* pRunManager = G4RunManager::GetRunManager();
    sprintf(asciiFileName,"%s%s%04d%s",(NumiData->asciiName).c_str(),"_",pRunManager->GetCurrentRun()->GetRunID(),".txt");
    G4cout << "Creating ASCII output file : "<<asciiFileName<<G4endl;
    std::ofstream asciiFile(asciiFileName);
  }

  if (NumiData->createBXDRAW) {
    sprintf(bxdrawFileName,"%s%s%03d%s",(NumiData->bxdrawName).c_str(),"_",pRunManager->GetCurrentRun()->GetRunID(),".dat");
    G4cout << "Creating BXDRAW output file : "<< bxdrawFileName <<G4endl;
    std::ofstream bxdrawFile(bxdrawFileName);
  }
  if(NumiData->createZpNtuple || NumiData->raytracing){
    sprintf(zpNtupleFileName,"%s_%04d%s.root",(NumiData->zpNtupleName).c_str(),pRunManager->GetCurrentRun()->GetRunID(),(NumiData->geometry).c_str());
    G4cout<<"Creating validation and muon tracking ntuple:"<<zpNtupleFileName<<G4endl;
    zpNtuple = new TFile(zpNtupleFileName, "RECREATE","zptrack ntuple");
    if(zpNtuple){
      zpNtuple->cd();
      zptree = new TTree("zp","g4numi Tracking particles through ZPoints");
      zptree->Branch("zpdata", "zptuple_t", &g4zpdata,32000,1);
      zpmetatree = new TTree("zpmeta","zpNtuple run metadata");
      zpmetatree->Branch("pots", &fZpNPots, "pots/I");
      zpmetatree->Branch("n_general_zplanes", &fZpNGeneral, "n_general_zplanes/I");
      zpmetatree->Branch("n_muon_zplanes", &fZpNMuon, "n_muon_zplanes/I");
      
      if(NumiData->raytracing){
      mesontree = new TTree("meson", "g4numi Generated Mesons");
      mesontree->Branch("px", &meson_px, "px/D");
      mesontree->Branch("py", &meson_py, "py/D");
      mesontree->Branch("pz", &meson_pz, "pz/D");
      mesontree->Branch("E", &meson_E, "E/D");
      mesontree->Branch("x", &meson_x, "x/D");
      mesontree->Branch("y", &meson_y, "y/D");
      mesontree->Branch("z", &meson_z, "z/D");
      mesontree->Branch("t", &meson_t, "t/D");
      mesontree->Branch("weight", &meson_weight, "weight/D");
      mesontree->Branch("pdg", &meson_pdg, "pdg/I");
      mesontree->Branch("tgen", &meson_tgen, "tgen/I");

      // Ancestry of each meson.  Off by default; the branches only appear when
      // /NuMI/output/outputMesonProvenance is set, so existing readers of the
      // meson tree are unaffected.
      if(NumiData->createMesonProvenance){
         mesontree->Branch("proc",   &meson_proc,   "proc/I");
         mesontree->Branch("vol",    &meson_vol,    "vol/I");
         mesontree->Branch("intgt",  &meson_intgt,  "intgt/I");
         mesontree->Branch("ppdg",   &meson_ppdg,   "ppdg/I");
         mesontree->Branch("pproc",  &meson_pproc,  "pproc/I");
         mesontree->Branch("pvol",   &meson_pvol,   "pvol/I");
         mesontree->Branch("pintgt", &meson_pintgt, "pintgt/I");
         mesontree->Branch("pE",     &meson_pE,     "pE/D");
         mesontree->Branch("pvx",    &meson_pvx,    "pvx/D");
         mesontree->Branch("pvy",    &meson_pvy,    "pvy/D");
         mesontree->Branch("pvz",    &meson_pvz,    "pvz/D");
         mesontree->Branch("cflag",  &meson_cflag,  "cflag/I");
         mesontree->Branch("cpdg",   &meson_cpdg,   "cpdg/I");
         mesontree->Branch("cproc",  &meson_cproc,  "cproc/I");
         mesontree->Branch("cgen",   &meson_cgen,   "cgen/I");
         mesontree->Branch("cE",     &meson_cE,     "cE/D");
         mesontree->Branch("cvx",    &meson_cvx,    "cvx/D");
         mesontree->Branch("cvy",    &meson_cvy,    "cvy/D");
         mesontree->Branch("cvz",    &meson_cvz,    "cvz/D");
         mesontree->Branch("ivx",    &meson_ivx,    "ivx/D");
         mesontree->Branch("ivy",    &meson_ivy,    "ivy/D");
         mesontree->Branch("ivz",    &meson_ivz,    "ivz/D");

         provnamestree = new TTree("provnames", "meson provenance id legend");
         provnamestree->Branch("id",   &provname_id,   "id/I");
         provnamestree->Branch("kind", &provname_kind, "kind/I");
         provnamestree->Branch("name", provname_name,  "name[80]/C");
      }

      if(NumiData->createTgtExitTree){
      tgtexittree = new TTree("tgtexit", "g4numi Mesons Leaving the Target");
      tgtexittree->Branch("px", &tgtexit_px, "px/D");
      tgtexittree->Branch("py", &tgtexit_py, "py/D");
      tgtexittree->Branch("pz", &tgtexit_pz, "pz/D");
      tgtexittree->Branch("E", &tgtexit_E, "E/D");
      tgtexittree->Branch("x", &tgtexit_x, "x/D");
      tgtexittree->Branch("y", &tgtexit_y, "y/D");
      tgtexittree->Branch("z", &tgtexit_z, "z/D");
      tgtexittree->Branch("t", &tgtexit_t, "t/D");
      tgtexittree->Branch("weight", &tgtexit_weight, "weight/D");
      tgtexittree->Branch("pdg", &tgtexit_pdg, "pdg/I");
      tgtexittree->Branch("tgen", &tgtexit_tgen, "tgen/I");
      }
      }

      // EM cascade seeds.  Off by default: the tree only exists when
      // /NuMI/output/outputEMSeeds is set, so the file is otherwise unchanged.
      // Not compatible with the delta-ray / muon-beam modes, which track
      // photons and e+- themselves.
      if(NumiData->createEMSeeds){
         if(NumiData->GetSimDRays() || NumiData->useMuonBeam){
            G4cout << "NumiAnalysis: outputEMSeeds ignored, simDRays or "
                   << "useMuonBeam is on" << G4endl;
            NumiData->OutputEMSeeds(false);
         } else {
            emseedtree = new TTree("emseed", "g4numi EM cascade seeds (gamma, e+-)");
            emseedtree->Branch("evt",        &es_evt,        "evt/I");
            emseedtree->Branch("trackID",    &es_trackID,    "trackID/I");
            emseedtree->Branch("parentID",   &es_parentID,   "parentID/I");
            emseedtree->Branch("pdg",        &es_pdg,        "pdg/I");
            emseedtree->Branch("parent_pdg", &es_parent_pdg, "parent_pdg/I");
            emseedtree->Branch("cproc",      &es_cproc,      "cproc/I");
            emseedtree->Branch("nimpwt",     &es_nimpwt,     "nimpwt/D");
            emseedtree->Branch("tgen",       &es_tgen,       "tgen/I");
            emseedtree->Branch("E",   &es_E,   "E/D");
            emseedtree->Branch("px",  &es_px,  "px/D");
            emseedtree->Branch("py",  &es_py,  "py/D");
            emseedtree->Branch("pz",  &es_pz,  "pz/D");
            emseedtree->Branch("x",   &es_x,   "x/D");
            emseedtree->Branch("y",   &es_y,   "y/D");
            emseedtree->Branch("z",   &es_z,   "z/D");
            emseedtree->Branch("t",   &es_t,   "t/D");
            emseedtree->Branch("vol", &es_vol, "vol/I");
            emseedtree->Branch("mat", &es_mat, "mat/I");
            emseedtree->Branch("iproc", &es_iproc, "iproc/I");
            emseedtree->Branch("ix",    &es_ix,    "ix/D");
            emseedtree->Branch("iy",    &es_iy,    "iy/D");
            emseedtree->Branch("iz",    &es_iz,    "iz/D");
            emseedtree->Branch("it",    &es_it,    "it/D");
            emseedtree->Branch("iE",    &es_iE,    "iE/D");
            emseedtree->Branch("ivol",  &es_ivol,  "ivol/I");
            emseedtree->Branch("imat",  &es_imat,  "imat/I");

            emseednamestree = new TTree("emseednames", "emseed id legend");
            emseednamestree->Branch("id",   &esname_id,   "id/I");
            emseednamestree->Branch("kind", &esname_kind, "kind/I");
            emseednamestree->Branch("name", esname_name,  "name[80]/C");
            G4cout << "NumiAnalysis: emseed tree on, threshold E >= "
                   << NumiData->EMSeedMinEnergy/CLHEP::MeV << " MeV" << G4endl;
         }
      }
    }
    else{
      G4Exception("NumiAnalysis","NumiAnalysis",FatalException,"Something went wrong with zpNtuple");
    }
    G4cout<<" End of if statement"<<G4endl;
  }

  //book histograms
}

//------------------------------------------------------------------------------------
std::string NumiAnalysis::GetOFileName(std::string ifilename)
{
   std::string filestart = ifilename;

   string::size_type slash_pos = ifilename.find_last_of('/');
   if(slash_pos == string::npos)
   {
      cerr << "MakeMuNtp::GetOFileName() - bad filename" << endl;
      return string("temp.root");
   }

   ifilename.erase(0, slash_pos+1);

   filestart.erase(slash_pos+1);

   string::size_type loc_underscore = ifilename.find("_", 0);
   if(loc_underscore == string::npos)
   {
      cerr << "MakeMuNtp::GetOFileName() - bad filename" << endl;
      return string("temp.root");
   }

   ifilename.erase(0, loc_underscore);

   string::size_type dot_pos = ifilename.find('.');
   if(dot_pos == string::npos)
   {
      cerr << "MakeMuNtp::GetOFileName() - can't find \".\" " << endl;
      return string("temp.root");
   }

   string filePrefix = "";
   ///////////////////////////////////////////////////////////////////////////////////////

   if(NumiData->GetSimDRays()) filePrefix = "dRayNtuple";
   else filePrefix = "hadmmNtuple";

///////////////////////////////////////////////////////////////////////////////////////
   /*
   filePrefix = "hadmmNtuple";
   */

   const G4int    NInputPart  = NumiData->GetNInputPart();
        if(NInputPart == 1)   ifilename.insert(dot_pos, "a");
   else if(NInputPart == 2)   ifilename.insert(dot_pos, "b");
   else if(NInputPart == 3)   ifilename.insert(dot_pos, "c");
   else if(NInputPart == 4)   ifilename.insert(dot_pos, "d");
   else if(NInputPart == 5)   ifilename.insert(dot_pos, "e");
   else if(NInputPart == 6)   ifilename.insert(dot_pos, "f");
   else if(NInputPart == 7)   ifilename.insert(dot_pos, "g");
   else if(NInputPart == 8)   ifilename.insert(dot_pos, "h");
   else if(NInputPart == 9)   ifilename.insert(dot_pos, "i");
   else if(NInputPart == 10)  ifilename.insert(dot_pos, "j");
   else if(NInputPart == 11)  ifilename.insert(dot_pos, "k");
   else if(NInputPart == 12)  ifilename.insert(dot_pos, "l");
   else if(NInputPart == 13)  ifilename.insert(dot_pos, "m");
   else if(NInputPart == 14)  ifilename.insert(dot_pos, "n");
   else if(NInputPart == 15)  ifilename.insert(dot_pos, "o");
   else if(NInputPart == 16)  ifilename.insert(dot_pos, "p");
   else if(NInputPart == 17)  ifilename.insert(dot_pos, "q");
   else if(NInputPart == 18)  ifilename.insert(dot_pos, "r");
   else if(NInputPart == 19)  ifilename.insert(dot_pos, "s");
   else if(NInputPart == 20)  ifilename.insert(dot_pos, "t");


   stringstream ofilename;
   ofilename << NumiData->GetHadmmNtupleDir() << "/" << filePrefix << ifilename;


   return ofilename.str();
}

//------------------------------------------------------------------------------------
void NumiAnalysis::finishPL()
{
  if (energyBinSimpleHistoMinerva < 0.) return;
   G4cout << "NumiAnalysis::finishPL() called." << G4endl;
   std::cerr << "NumiAnalysis::finishPL() called." << G4endl;

  // Study of the inclusive  numu spectrum for various magnetic field
  // To speed the I/O and processing on the farm, we include here a fixed size
  // weighted histogram
  std::string aDir("./");
  const char* aDirGG = getenv("_CONDOR_SCRATCH_DIR");
  if (aDirGG != 0) aDir = std::string(aDirGG);
  {
    std::string aNamePL(aDir); aNamePL += std::string("/FluxHistoPL.txt");
    std::ofstream fOutPLHisto(aNamePL.c_str());
    fOutPLHisto << " k e numu numub " << std::endl;
    G4RunManager* pRunManager = G4RunManager::GetRunManager();
     int NPots = pRunManager->GetCurrentRun()->GetNumberOfEventToBeProcessed();
    for (size_t k=0; k != MinervaNuMuHisto.size();  k++) {
       fOutPLHisto << " " << k << " " <<  (energyBinSimpleHistoMinerva/2. + energyBinSimpleHistoMinerva*k)
                << " " << MinervaNuMuHisto[k]/NPots << " "
                << MinervaNuMuBarHisto[k]/NPots << std::endl;
    }
    fOutPLHisto.close();
  }
  {
    std::string aNamePL(aDir); aNamePL += std::string("/FluxHistoNovaPL.txt");
    std::ofstream fOutPLHisto(aNamePL.c_str());
    fOutPLHisto << " k e numuN numubN numuF numubF " << std::endl;
    G4RunManager* pRunManager = G4RunManager::GetRunManager();
    int NPots = pRunManager->GetCurrentRun()->GetNumberOfEventToBeProcessed();
    for (size_t k=0; k != NovaNearNuMuHisto.size();  k++) {
      fOutPLHisto << " " << k << " " <<  (0.05 + 0.1*k)
                << " " << NovaNearNuMuHisto[k]/NPots << " "  << NovaNearNuMuBarHisto[k]/NPots
                << " " << NovaFarNuMuHisto[k]/NPots << " "
                << NovaFarNuMuBarHisto[k]/NPots << std::endl;
    }
    fOutPLHisto.close();
  }
}
void NumiAnalysis::finish() {

   G4cout << "NumiAnalysis::finish() called." << G4endl;
   std::cerr << "NumiAnalysis::finish() called." << G4endl;

  if (NumiData->createNuNtuple){
    nuNtuple->cd();
    meta->Write(); //write dkmeta
    tree->Write();
    nuNtuple->Close();
    delete nuNtuple;
  }

  if (NumiData->createHadmmNtuple){
    hadmmNtuple->cd();
    hadmmtree->Write();
    hadmmNtuple->Close();
    delete hadmmNtuple;
  }

  if (NumiData->GetCreateAbsBkgNtuple()){
    absbkgNtuple->cd();
    absbkgtree->Write();
    absbkgNtuple->Close();
    delete absbkgNtuple;
  }

  if (NumiData->createZpNtuple || NumiData->raytracing){
    zpNtuple->cd();
    zptree->Write();
    fZpNPots = G4RunManager::GetRunManager()->GetCurrentRun()->GetNumberOfEventToBeProcessed();
    fZpNGeneral = NumiData->GetNGeneralZPlanes();
    fZpNMuon = NumiData->GetNMuonZPlanes();
    zpmetatree->Fill();
    zpmetatree->Write();
    if (mesontree) mesontree->Write();
    if (tgtexittree) tgtexittree->Write();
    if (provnamestree){
       // dump the id legend once, at end of run, from the maps built on the fly
       for(int kind = 0; kind < 2; ++kind){
          const std::map<G4String, G4int>& ids = (kind == 0) ? fProcIds : fVolIds;
          for(std::map<G4String, G4int>::const_iterator it = ids.begin();
              it != ids.end(); ++it){
             provname_id   = it->second;
             provname_kind = kind;
             strncpy(provname_name, it->first.c_str(), sizeof(provname_name) - 1);
             provname_name[sizeof(provname_name) - 1] = 0;
             provnamestree->Fill();
             G4cout << "provname " << (kind == 0 ? "proc " : "vol  ")
                    << it->second << " " << it->first << G4endl;
          }
       }
       provnamestree->Write();
    }
    if (emseedtree){
       emseedtree->Write();
       // legend: id 0 of each kind is "none"/unknown
       for(int kind = 0; kind < 3; ++kind){
          const std::map<G4String, G4int>& ids =
             (kind == 0) ? fESProcIds : (kind == 1) ? fESVolIds : fESMatIds;
          esname_id = 0;  esname_kind = kind;
          strncpy(esname_name, kind == 0 ? "none" : "unknown", sizeof(esname_name) - 1);
          esname_name[sizeof(esname_name) - 1] = 0;
          emseednamestree->Fill();
          for(std::map<G4String, G4int>::const_iterator it = ids.begin();
              it != ids.end(); ++it){
             esname_id   = it->second;
             esname_kind = kind;
             strncpy(esname_name, it->first.c_str(), sizeof(esname_name) - 1);
             esname_name[sizeof(esname_name) - 1] = 0;
             emseednamestree->Fill();
          }
       }
       emseednamestree->Write();

       // run summary: seed counts and the pi0/eta the stacking cuts drop
       TTree* esmeta = new TTree("emseedmeta", "emseed run summary");
       Int_t    m_pots = fZpNPots;
       Double_t m_emin = NumiData->EMSeedMinEnergy/CLHEP::GeV;
       Double_t m_kthr = NumiData->GetKillTracking() ? NumiData->GetKillTrackingThreshold()/CLHEP::GeV : 0.;
       Long64_t m_nph = fESNPhoton, m_nlep = fESNLepton;
       // order: pi0 KE<thr, pi0 pz<0, eta KE<thr, eta pz<0
       Long64_t m_kn[4];
       Double_t m_kw[4], m_kwe[4];
       for(int i = 0; i < 2; ++i) for(int j = 0; j < 2; ++j){
          m_kn[2*i+j] = fKilledN[i][j]; m_kw[2*i+j] = fKilledW[i][j]; m_kwe[2*i+j] = fKilledWE[i][j];
       }
       esmeta->Branch("pots",      &m_pots, "pots/I");
       esmeta->Branch("emin",      &m_emin, "emin/D");
       esmeta->Branch("kill_thr",  &m_kthr, "kill_thr/D");
       esmeta->Branch("n_photon",  &m_nph,  "n_photon/L");
       esmeta->Branch("n_lepton",  &m_nlep, "n_lepton/L");
       esmeta->Branch("killed_n",  m_kn,    "killed_n[4]/L");
       esmeta->Branch("killed_w",  m_kw,    "killed_w[4]/D");
       esmeta->Branch("killed_wE", m_kwe,   "killed_wE[4]/D");
       esmeta->Fill();
       esmeta->Write();

       G4cout << "EMSeed summary: " << fESNPhoton << " photon seeds, "
              << fESNLepton << " e+- seeds over " << fZpNPots << " POT ("
              << (fZpNPots > 0 ? double(fESNPhoton + fESNLepton)/fZpNPots : 0.)
              << " per POT); photons without a creation cache: " << fESNNoCache << G4endl;
       const char* sp[2] = {"pi0", "eta"};
       const char* rs[2] = {"KE<KillTrackingThreshold", "pz<0"};
       for(int i = 0; i < 2; ++i) for(int j = 0; j < 2; ++j){
          G4cout << "EMSeed killed " << sp[i] << " " << rs[j] << ": n=" << fKilledN[i][j]
                 << " sum(w)/POT=" << (fZpNPots > 0 ? fKilledW[i][j]/fZpNPots : 0.)
                 << " sum(w*E)/POT=" << (fZpNPots > 0 ? fKilledWE[i][j]/fZpNPots : 0.)
                 << " GeV" << G4endl;
       }
    }
    zpNtuple->Close();
    delete zpNtuple;
  }

  if (NumiData->createTarNtuple){
    tarNtuple->cd();
    tartree->Write();
    tarNtuple->Close();
    delete tarNtuple;
  }
}

void NumiAnalysis::FillMeta(){


  //This class fills dk2meta:
  G4RunManager* pRunManager = G4RunManager::GetRunManager();
  int NPots = pRunManager->GetCurrentRun()->GetNumberOfEventToBeProcessed();

  G4String namentp = (NumiData->nuNtupleName);
  G4int pos_last = namentp.rfind('_');
  namentp.erase(0,pos_last+1);
  istringstream buffer(namentp);
  int valjob;
  buffer >> valjob;
  //if we are running interactevely (or more than 1000 jobs), we will write
  //job ID = -1
  if(valjob<0 || valjob>=1000)valjob=-1;

  NumiPrimaryGeneratorAction *NPGA = (NumiPrimaryGeneratorAction*)(pRunManager)->GetUserPrimaryGeneratorAction();

  this_meta->job  = valjob;
  this_meta->pots = NPots;

  //We are going to do this manually for now:
  //Perhaps should be part of NumiDataInput??

#ifdef USEMODGEANT4
  this_meta->beamsim = "g4numi_v6_MODGEANT4";
#else
  this_meta->beamsim = "g4numi_v6";
#endif

  // try to get the version of Geant4
  /// static const G4String G4Version = "$Name: geant4-10-04-patch-02b $";
  G4String version = G4Version.substr(7,string::npos);
  // get the physics list as used
  G4String physListName = NumiData->GetPhysicsListName();

//#ifdef THE_HARD_WAY
//  const G4VUserPhysicsList* pl = pRunManager->GetUserPhysicsList();
//  const std::type_info& ti = typeid(*pl);
//  G4String mangled = ti.name();
//#if defined(HAVE_CXA_DEMANGLE) || defined(__clang__) || defined(__GNUC__)
//  int status;
//  char *realname = abi::__cxa_demangle(ti.name(),0,0,&status);
//  G4String demangled = realname;
//  free(realname);
//#else
//  G4String demangled = mangled;  // will have to do
//#endif
//  G4String physListName = demanangled;
//#endif

  //this_meta->physics = "geant4_9_2_p03_FTFP_BERT1.0";
  this_meta->physics = version + " " + physListName;
  this_meta->physcuts = "nofillyet";

  G4String hornC = NumiData->GetBeamConfig();
  G4String tgtC  = NumiData->GetBeamConfig();
  G4int confsize = (NumiData->GetBeamConfig()).length();
  // G4String::remove(pos,N) erased N+pos characters, not N -- keep that.
  hornC.erase(0,6);
  tgtC.erase(6,confsize+6);
  G4String playlist  = NumiData->GetPlaylist();

  G4bool isHe = NumiData->HeInDecayPipe;
  std::string decayVolMat = "vacuum";
  if(isHe)decayVolMat = "helium";
  this_meta->tgtcfg = std::string(tgtC+"_"+playlist);
  this_meta->horncfg = std::string(hornC);
  this_meta->dkvolcfg = decayVolMat;
  //////////////////////////////////////////

  this_meta->beam0x = NumiData->beamPosition[0]/CLHEP::cm;
  this_meta->beam0y = NumiData->beamPosition[1]/CLHEP::cm;
  this_meta->beam0z = NumiData->beamPosition[2]/CLHEP::cm;
  this_meta->beamhwidth = NumiData->beamSigmaX/CLHEP::cm;
  this_meta->beamvwidth = NumiData->beamSigmaY/CLHEP::cm;

  G4ThreeVector protonMomentum = NPGA->GetProtonMomentum();
  this_meta->beamdxdz = protonMomentum[0]/protonMomentum[2];
  this_meta->beamdydz = protonMomentum[1]/protonMomentum[2];

  G4double mm2cm = 0.1;
  vec_loc.push_back(bsim::Location(mm2cm*0.0,mm2cm*0.0,mm2cm*0.0,"random"));
  for (G4int ii=0; ii<NumiData->nNear; ++ii){
    G4double xxdet = mm2cm*(NumiData->xdet_near[ii]);
    G4double yydet = mm2cm*(NumiData->ydet_near[ii]);
    G4double zzdet = mm2cm*(NumiData->zdet_near[ii]);
    vec_loc.push_back(bsim::Location(xxdet,yydet,zzdet,NumiData->det_near_name[ii]));
  }
  for (G4int ii=0; ii<NumiData->nFar; ++ii){
    G4double xxdet = mm2cm*(NumiData->xdet_far[ii]);
    G4double yydet = mm2cm*(NumiData->ydet_far[ii]);
    G4double zzdet = mm2cm*(NumiData->zdet_far[ii]);
    vec_loc.push_back(bsim::Location(xxdet,yydet,zzdet,NumiData->det_far_name[ii]));
  }

  this_meta->location = vec_loc;

  for(int ii=0;ii<NumiData->nVintTot;ii++){
    (this_meta->vintnames).push_back(NumiData->VolVintName[ii]);
  }
  for(int ii=0;ii<NumiData->nVdblTot;ii++){
    (this_meta->vdblnames).push_back(NumiData->VolVdblName[ii]);
  }
  if (NumiData->createNuNtuple) meta->Fill(); // If we don't have the Ntuple properly open, won't work.

}

//------------------------------------------------------------------------------------
void NumiAnalysis::FillHadmmNtuple(const G4Track& track, Int_t hmm_num, Int_t cellNum)
{
   if(NumiData->GetDebugLevel() == 10) { G4cout << "NumiAnalysis::FillHadmmNtuple() called." << G4endl;}

   if (!NumiData->createHadmmNtuple) return;

   G4ParticleDefinition* particleDefinition = track.GetDefinition();

   if ( track.GetTrackID() != 1 ||
        (particleDefinition != G4MuonPlus::MuonPlusDefinition()
         && particleDefinition != G4MuonMinus::MuonMinusDefinition()) )
   {
      G4cout << "******** Entry " << fentry
           << ": Problem in NumiAnalysis::FillHadmmNtuple() - The particle is either"
           << " not a muon or the not the orignial muon. Particle = "
           << particleDefinition -> GetParticleName() << " the track number = "
           << track.GetTrackID() << G4endl;
      return;
   }

   if(g4hmmdata)
   {
      if(!(NumiData->useMuonInput))
      {
         //g4hmmdata->hmmenergy = track.GetTotalEnergy();

         if(hmm_num == 4)
         {
            /*g4hmmdata->hmmxpos = track.GetPosition()[0];
              g4hmmdata->hmmypos = track.GetPosition()[1];
              g4hmmdata->hmmzpos = track.GetPosition()[2];
              g4hmmdata->hmmpx = track.GetMomentum()[0];
              g4hmmdata->hmmpy = track.GetMomentum()[1];
              g4hmmdata->hmmpz = track.GetMomentum()[2];
            */
         }
      }
      if(hmm_num != 4)
      {
         if(fMuInAlc[hmm_num] == true) return;

         if( !(NumiData->GetMaterialSigma() < 0) && !(NumiData->GetMaterialSigma() > 0))
         {
            g4hmmdata->mmxpos[hmm_num] = track.GetPosition()[0];
            g4hmmdata->mmypos[hmm_num] = track.GetPosition()[1];
            //g4hmmdata->mmzpos[hmm_num] = track.GetPosition()[2];
         }

         g4hmmdata->mmpx[hmm_num]   = track.GetMomentum()[0];
         g4hmmdata->mmpy[hmm_num]   = track.GetMomentum()[1];
         g4hmmdata->mmpz[hmm_num]   = track.GetMomentum()[2];
         g4hmmdata->cell[hmm_num]   = cellNum;

         /*if(hmm_num == 0)
           G4cout << GetEntry() << ": mmpz = " << g4hmmdata->mmpz[hmm_num]
           << " mmxpos = " << g4hmmdata->mmxpos[hmm_num]
           << " mmcell = " << g4hmmdata->cell[hmm_num] << G4endl;
         */

         fMuInAlc[hmm_num] = true;
      }
   }
   else if(g4draydataMIB)
   {
      if(hmm_num != 4)
      {
         if(fMuInAlc[hmm_num] == true) return;

         g4draydataMIB->mmxpos[hmm_num] = track.GetPosition()[0];
         g4draydataMIB->mmpx[hmm_num]   = track.GetMomentum()[0];
         g4draydataMIB->mmypos[hmm_num] = track.GetPosition()[1];
         g4draydataMIB->mmpy[hmm_num]   = track.GetMomentum()[1];
         //g4draydataMIB->mmzpos[hmm_num] = track.GetPosition()[2];
         g4draydataMIB->mmpz[hmm_num]   = track.GetMomentum()[2];
         g4draydataMIB->cell[hmm_num]   = cellNum;

         fMuInAlc[hmm_num] = true;
      }
   }
   else if(g4draydataSPB)
   {
      if(hmm_num != 4)
      {
         if(fMuInAlc[hmm_num] == true) return;

         g4draydataSPB->mmxpos[hmm_num] = track.GetPosition()[0];
         g4draydataSPB->mmpx[hmm_num]   = track.GetMomentum()[0];
         g4draydataSPB->mmypos[hmm_num] = track.GetPosition()[1];
         g4draydataSPB->mmpy[hmm_num]   = track.GetMomentum()[1];
         //g4draydataSPB->mmzpos[hmm_num] = track.GetPosition()[2];
         g4draydataSPB->mmpz[hmm_num]   = track.GetMomentum()[2];
         g4draydataSPB->cell[hmm_num]   = cellNum;

         fMuInAlc[hmm_num] = true;
      }
   }

}

//------------------------------------------------------------------------------------
void NumiAnalysis::FillEdep(G4double edep, const G4ParticleDefinition* particleDefinition,
                            const G4int alc, const G4int IntExt, const G4double /*zpos*/, Int_t cellNum,
                            const Int_t trackID, const G4double weight)
{

   if(NumiData->GetDebugLevel() == 10) { G4cout << "NumiAnalysis::FillEdep() called." << G4endl;}

   edep = edep/CLHEP::eV;

   if(g4draydataMIB)
   {
      if ( particleDefinition == G4MuonPlus::MuonPlusDefinition()
           || particleDefinition == G4MuonMinus::MuonMinusDefinition() )
      {
         g4draydataMIB -> SetMuEdep(alc, cellNum, edep);
      }
      else if ( particleDefinition == G4Electron::ElectronDefinition()
                || particleDefinition == G4Positron::PositronDefinition())
      {
         if(IntExt == 0) // ext
         {
            g4draydataMIB -> SetExtEdep(alc, cellNum, edep, trackID, weight);
            //g4draydataMIB->zpos_edep[alc] = zpos;
         }
         else if(IntExt == 1) // int
         {
            g4draydataMIB -> SetIntEdep(alc, cellNum, edep, trackID, weight);
            //g4draydataMIB->zpos_edep[alc] = zpos;
         }
         else
            G4cerr << "NumiAnalysis::FillEdep - Problem: Invalid IntExt flag " << IntExt << G4endl;
      }
      else
         G4cerr << "NumiAnalysis::FillEdep - Problem: Invalid particle type " << particleDefinition -> GetParticleName() << G4endl;
   }
   else if(g4draydataSPB)
   {
      if ( particleDefinition == G4MuonPlus::MuonPlusDefinition()
           || particleDefinition == G4MuonMinus::MuonMinusDefinition() )
      {
         g4draydataSPB->mu_edep[alc] += edep;
      }
      else if ( particleDefinition == G4Electron::ElectronDefinition()
                || particleDefinition == G4Positron::PositronDefinition() )
      {
         if(IntExt == 0) // ext
         {
            g4draydataSPB->ext_edep[alc] += edep;
         }
         else if(IntExt == 1) // int
         {
            g4draydataSPB->int_edep[alc] += edep;
         }
         else
            G4cerr << "NumiAnalysis::FillEdep - Problem: Invalid IntExt flag " << IntExt << G4endl;
      }
      else
         G4cerr << "NumiAnalysis::FillEdep - Problem: Invalid particle type " << particleDefinition -> GetParticleName() << G4endl;
   }

}





//
//record mu momentum and position at 0.5m into the rock upstream
//of the monitors to apply delta ray corrections
//
//------------------------------------------------------------------------------------
void NumiAnalysis::FillAlcEdepInfo(const G4Track& track, G4int alc)
{
   if(NumiData->GetDebugLevel() == 10) { G4cout << "NumiAnalysis::FillAlcEdepInfo() called." << G4endl;}

  if(fAlcEdep_called[alc] == true) return;

  //G4cout << "fentry = " << fentry
  //<< " Calling NumiAnalysis::FillAlcEdepInfo() for alc "
  //<< alc << " zpos = " << track.GetPosition()[2] << G4endl;

  //g4hmmdata->mmxpos_Edep[alc] = track.GetPosition()[0];
  g4hmmdata->mmpx_Edep[alc]   = track.GetMomentum()[0];
  //g4hmmdata->mmypos_Edep[alc] = track.GetPosition()[1];
  g4hmmdata->mmpy_Edep[alc]   = track.GetMomentum()[1];
  //g4hmmdata->mmzpos_Edep[alc] = track.GetPosition()[2];
  g4hmmdata->mmpz_Edep[alc]   = track.GetMomentum()[2];

  fAlcEdep_called[alc] = true;

}


//
//this function only gets called if useMuonInput = true
//
//------------------------------------------------------------------------------------
void NumiAnalysis::FillHadmmNtuple()
{
   if(NumiData->GetDebugLevel() == 10) { G4cout << "NumiAnalysis::FillHadmmNtuple() (with no args) called." << G4endl;}

  if (!NumiData->createHadmmNtuple) return;

  G4RunManager* pRunManager = G4RunManager::GetRunManager();
  NumiPrimaryGeneratorAction *NPGA = (NumiPrimaryGeneratorAction*)(pRunManager)->GetUserPrimaryGeneratorAction();

  if(g4hmmdata)
  {

     if( !(NumiData->GetMaterialSigma() < 0) && !(NumiData->GetMaterialSigma() > 0))
     {
        G4ThreeVector particlePosition = NPGA->GetParticlePosition();
        g4hmmdata->muvx = particlePosition[0];
        g4hmmdata->muvy = particlePosition[1];
        //g4hmmdata->muvz = particlePosition[2];

        G4ThreeVector tparentPosition = NPGA->GetMuTParentPosition();
        g4hmmdata->tvx = tparentPosition[0];
        g4hmmdata->tvy = tparentPosition[1];
        g4hmmdata->tvz = tparentPosition[2];

        G4ThreeVector parentMomentum = NPGA->GetMuParentMomentum();
        g4hmmdata->pdpx = parentMomentum[0];
        g4hmmdata->pdpy = parentMomentum[1];
        g4hmmdata->pdpz = parentMomentum[2];

        G4ThreeVector parentPosition = NPGA->GetMuParentPosition();
        g4hmmdata->pdvx = parentPosition[0];
        g4hmmdata->pdvy = parentPosition[1];
        g4hmmdata->pdvz = parentPosition[2];

        G4ThreeVector parentProdPosition = NPGA->GetMuParentProdPosition();
        g4hmmdata->ppvx = parentProdPosition[0];
        g4hmmdata->ppvy = parentProdPosition[1];
        g4hmmdata->ppvz = parentProdPosition[2];

        g4hmmdata->pptype = NPGA->GetMuParentType();
        g4hmmdata->ppmedium = NPGA->GetMuParentProdMedium();
        g4hmmdata->pgen = NPGA->GetMuParentGen();
        g4hmmdata->ptype = NPGA->GetMuType();

        g4hmmdata->evtno = NPGA->GetEvtno();
     }

     G4ThreeVector particleMomentum = NPGA->GetParticleMomentum();
     g4hmmdata->mupx = particleMomentum[0];
     g4hmmdata->mupy = particleMomentum[1];
     g4hmmdata->mupz = particleMomentum[2];

     G4ThreeVector tparentMomentum = NPGA->GetMuTParentMomentum();
     g4hmmdata->tpx = tparentMomentum[0];
     g4hmmdata->tpy = tparentMomentum[1];
     g4hmmdata->tpz = tparentMomentum[2];

     g4hmmdata->muweight = NPGA->GetMuWeight();
     g4hmmdata->tpptype = NPGA->GetTParentType();
     g4hmmdata->nimpwt = NPGA->GetImpWeight();

     //g4hmmdata->run = pRunManager->GetCurrentRun()->GetRunID();
     //g4hmmdata->mtgthsig = NumiData->beamSigmaX/CLHEP::cm;
     //g4hmmdata->mtgtvsig = NumiData->beamSigmaY/CLHEP::cm;
     //g4hmmdata->mtgthpos = NumiData->beamPosition[0]/CLHEP::cm;
     //g4hmmdata->mtgtvpos = NumiData->beamPosition[1]/CLHEP::cm;
     //g4hmmdata->evtno = pRunManager->GetCurrentEvent()->GetEventID();



     //g4hmmdata->hmmenergy = -81000;
     //g4hmmdata->hmmxpos = -81000;
     //g4hmmdata->hmmypos = -81000;
     //g4hmmdata->hmmzpos = -81000;
     //g4hmmdata->hmmpx = -81000;
     //g4hmmdata->hmmpy = -81000;
     //g4hmmdata->hmmpz = -81000;

     for(Int_t i=0;i<3;++i)
     {
        if( !(NumiData->GetMaterialSigma() < 0) && !(NumiData->GetMaterialSigma() > 0))
        {
           g4hmmdata->mmxpos[i] = -99999.;
           g4hmmdata->mmypos[i] = -99999.;
           //g4hmmdata->mmzpos[i] = -99999.;
        }

        g4hmmdata->mmpx[i]   = -99999.;
        g4hmmdata->mmpy[i]   = -99999.;
        g4hmmdata->mmpz[i]   = -99999.;
        g4hmmdata->cell[i]   = -999;
     }
  }
  else if(g4draydataMIB)
  {
     G4ThreeVector particlePosition = NPGA->GetParticlePosition();
     g4draydataMIB->muvx = particlePosition[0];
     g4draydataMIB->muvy = particlePosition[1];
     //g4draydataMIB->muvz = particlePosition[2];

     G4ThreeVector particleMomentum = NPGA->GetParticleMomentum();
     g4draydataMIB->mupx = particleMomentum[0];
     g4draydataMIB->mupy = particleMomentum[1];
     g4draydataMIB->mupz = particleMomentum[2];

     g4draydataMIB->muweight = NPGA->GetMuWeight();
     g4draydataMIB->nimpwt = NPGA->GetImpWeight();
     g4draydataMIB->ptype = NPGA->GetMuType();

  }
  else if(g4draydataSPB)
  {
     G4ThreeVector particlePosition = NPGA->GetParticlePosition();
     g4draydataSPB->muvx = particlePosition[0];
     g4draydataSPB->muvy = particlePosition[1];
     g4draydataSPB->muvz = particlePosition[2];

     G4ThreeVector particleMomentum = NPGA->GetParticleMomentum();
     g4draydataSPB->mupx = particleMomentum[0];
     g4draydataSPB->mupy = particleMomentum[1];
     g4draydataSPB->mupz = particleMomentum[2];

     g4draydataSPB->ptype = NPGA->GetMuType();

  }


}

//------------------------------------------------------------------------------------
void NumiAnalysis::FillAbsorberBkgrdNtuple(const G4Track& track)
{
   if(NumiData->GetDebugLevel() == 10) { G4cout << "NumiAnalysis::FillAbsorberBkgrdNtuple() called." << G4endl;}

   G4ParticleDefinition * particleType = track.GetDefinition();
   G4String partType = particleType->GetParticleType();
   g4absbkgdata -> ptype = NumiParticleCode::AsInt(NumiParticleCode::StringToEnum(particleType->GetParticleName()));

   g4absbkgdata->x    = track.GetPosition()[0]/CLHEP::cm;
   g4absbkgdata->px   = track.GetMomentum()[0]/CLHEP::GeV;
   g4absbkgdata->y    = track.GetPosition()[1]/CLHEP::cm;
   g4absbkgdata->py   = track.GetMomentum()[1]/CLHEP::GeV;
   g4absbkgdata->z    = track.GetPosition()[2]/CLHEP::cm;
   g4absbkgdata->pz   = track.GetMomentum()[2]/CLHEP::GeV;

   g4absbkgdata->KE = track.GetKineticEnergy()/CLHEP::GeV;

   NumiTrackInformation* info = (NumiTrackInformation*)(track.GetUserInformation());
   g4absbkgdata->impwt = info->GetNImpWt();  // Importance weight
   g4absbkgdata->tgen = info->GetTgen()-1;

   int tgtz = -999;
   if(NumiData->TargetZ0/CLHEP::cm == -35.)
      tgtz = 0;
   else if(NumiData->TargetZ0/CLHEP::cm == -45.)
      tgtz = 10;
   else if(NumiData->TargetZ0/CLHEP::cm == -135.)
      tgtz = 100;
   else if(NumiData->TargetZ0/CLHEP::cm == -185.)
      tgtz = 150;
   else if(NumiData->TargetZ0/CLHEP::cm == -285.)
      tgtz = 250;

   g4absbkgdata->tgtz = tgtz;
   g4absbkgdata->ihorn = NumiData->HornCurrent/CLHEP::ampere/1000.;

   absbkgtree->Fill();
   g4absbkgdata->Clear();
}




//------------------------------------------------------------------------------------
void NumiAnalysis::SetAlcEdepFlag(G4bool AlcEdep)
{
  fAlcEdep_called[0] = AlcEdep;
  fAlcEdep_called[1] = AlcEdep;
  fAlcEdep_called[2] = AlcEdep;
}
//------------------------------------------------------------------------------------
void NumiAnalysis::SetMuInAlcFlag(G4bool MuInAlc)
{
  fMuInAlc[0] = MuInAlc;
  fMuInAlc[1] = MuInAlc;
  fMuInAlc[2] = MuInAlc;
}


//------------------------------------------------------------------------------------
void NumiAnalysis::SetCount(G4int count)
{
  fcount = count;
}
//------------------------------------------------------------------------------------
G4int NumiAnalysis::GetCount()
{
  return fcount;
}
//------------------------------------------------------------------------------------
void NumiAnalysis::SetEntry(G4int entry)
{
  fentry = entry;
}
//------------------------------------------------------------------------------------
G4int NumiAnalysis::GetEntry()
{
  return fentry;
}

//------------------------------------------------------------------------------------
void NumiAnalysis::WriteHadmmNtuple()
{
   if(NumiData->GetDebugLevel() == 10) { G4cout << "NumiAnalysis::WriteHadmmNtuple() called." << G4endl;}

   if (!(NumiData->createHadmmNtuple)) return;

   if(g4draydataMIB) g4draydataMIB->ClearTrackIDVectors();

   hadmmtree->Fill();

   if(g4hmmdata) g4hmmdata->Clear();
   if(g4draydataMIB) g4draydataMIB->Clear();
   if(g4draydataSPB) g4draydataSPB->Clear();


}

//------------------------------------------------------------------------------------
void NumiAnalysis::FillNeutrinoNtuple(const G4Track& track, const std::vector<G4VTrajectory*>& history)
{

   if(NumiData->GetDebugLevel() == 10) { G4cout << "NumiAnalysis::FillNeutrinoNtuple() called." << G4endl;}

  if (!NumiData->createNuNtuple) return;

  G4String myHorn1Config =  NumiData->GetHornConfig();
  const NumiRunManager* theRunManager = reinterpret_cast<const NumiRunManager*>(NumiRunManager::GetRunManager());
  const NumiSteppingAction* theSteppingAction =
     reinterpret_cast<const NumiSteppingAction*>(theRunManager->GetUserSteppingAction());

  //Neutrino vertex position and momentum
  G4ThreeVector pos = track.GetPosition()/CLHEP::mm;
  x = pos.x();
  y = pos.y();
  z = pos.z();
  G4ThreeVector NuMomentum = track.GetMomentum();
  G4int parentID = track.GetParentID();

  if (parentID == 0) return; //I have to make some changes so that neutrinos in fluka/mars ntuples don't crash

  NumiTrajectory* NuParentTrack = GetParentTrajectory(parentID);
  G4int point_no = NuParentTrack->GetPointEntries();
  G4ThreeVector ParentMomentumFinal = NuParentTrack->GetMomentum(point_no-1);
  G4ThreeVector vertex_r = (NuParentTrack->GetPoint(point_no-1)->GetPosition()/CLHEP::m)*CLHEP::m; //Should be the same as Neutrino vertex
  G4String parent_name = NuParentTrack->GetParticleName();
  G4double Parent_mass = NuParentTrack->GetMass();
  G4double gamma = sqrt(ParentMomentumFinal*ParentMomentumFinal+Parent_mass*Parent_mass)/Parent_mass;
  G4double Parent_energy = gamma*Parent_mass;
  G4ThreeVector beta_vec = ParentMomentumFinal/Parent_energy;
  G4double partial = gamma*(beta_vec*NuMomentum);

  G4double enuzr = gamma*(track.GetTotalEnergy())-partial; //neutrino energy in parent rest frame

  //fill histograms, ntuples,...
  G4RunManager* pRunManager = G4RunManager::GetRunManager();
  NumiPrimaryGeneratorAction *NPGA = (NumiPrimaryGeneratorAction*)(pRunManager)->GetUserPrimaryGeneratorAction();

  g4data->run = pRunManager->GetCurrentRun()->GetRunID();
  g4data->evtno = pRunManager->GetCurrentEvent()->GetEventID();
  g4data->beamHWidth = NumiData->beamSigmaX/CLHEP::cm;
  g4data->beamVWidth = NumiData->beamSigmaY/CLHEP::cm;
  g4data->beamX = NumiData->beamPosition[0]/CLHEP::cm;
  g4data->beamY = NumiData->beamPosition[1]/CLHEP::cm;

  G4int particleID = track.GetParentID();
  G4ThreeVector protonOrigin = NPGA->GetProtonOrigin();
  g4data->protonX = protonOrigin[0];
  g4data->protonY = protonOrigin[1];
  g4data->protonZ = protonOrigin[2];

  G4ThreeVector protonMomentum = NPGA->GetProtonMomentum();
  g4data->protonPx = protonMomentum[0];
  g4data->protonPy = protonMomentum[1];
  g4data->protonPz = protonMomentum[2];

  g4data->nuTarZ = NumiData->TargetZ0;
  g4data->hornCurrent = NumiData->HornCurrent/CLHEP::ampere/1000.;

  // Random decay - these neutrinos rarely hit any of the detectors
  g4data->Ndxdz = NuMomentum[0]/NuMomentum[2];
  g4data->Ndydz = NuMomentum[1]/NuMomentum[2];
  g4data->Npz = NuMomentum[2]/CLHEP::GeV;
  g4data->Nenergy = track.GetTotalEnergy()/CLHEP::GeV;

   //other info
  // Neutrino origin:
  // 3 From muon decay
  // 1 From particle from target
  // 2 From scraping
  //check if nu is from muon decay or from a particle from target, otherwise Norig = 2
  G4int Norig = 2;
  if ((parent_name=="mu+") || (parent_name=="mu-")) Norig = 3;
  G4String firstvolname = NuParentTrack->GetPreStepVolumeName(0);
  if (firstvolname.contains("Baffle") || firstvolname.contains("TGT")) Norig = 1;

  g4data->Norig = Norig;
  g4data->Ndecay = NuParentTrack->GetDecayCode();

  G4ParticleDefinition * particleType = track.GetDefinition();
  G4int ntype = NumiParticleCode::AsInt(NumiParticleCode::StringToEnum(particleType->GetParticleName()));
  g4data->Ntype = ntype;
  g4data->Vx = x/CLHEP::cm;
  g4data->Vy = y/CLHEP::cm;
  g4data->Vz = z/CLHEP::cm;
  g4data->pdPx = ParentMomentumFinal[0]/CLHEP::GeV;
  g4data->pdPy = ParentMomentumFinal[1]/CLHEP::GeV;
  g4data->pdPz = ParentMomentumFinal[2]/CLHEP::GeV;

  G4ThreeVector ParentMomentumProduction = NuParentTrack->GetMomentum(0);
  g4data->ppdxdz = ParentMomentumProduction[0]/ParentMomentumProduction[2];
  g4data->ppdydz = ParentMomentumProduction[1]/ParentMomentumProduction[2];
  g4data->pppz = ParentMomentumProduction[2]/CLHEP::GeV;

  G4double parentp = sqrt(ParentMomentumProduction*ParentMomentumProduction);

  g4data->ppenergy = sqrt((parentp*parentp+Parent_mass*Parent_mass))/CLHEP::GeV;

  g4data->ppmedium = 0.; //this is still empty

  g4data->ptype = NumiParticleCode::AsInt(NumiParticleCode::StringToEnum(parent_name));

  G4ThreeVector production_vertex = (NuParentTrack->GetPoint(0)->GetPosition()/CLHEP::m)*CLHEP::m;
  g4data->ppvx = production_vertex[0]/CLHEP::cm;
  g4data->ppvy = production_vertex[1]/CLHEP::cm;
  g4data->ppvz = production_vertex[2]/CLHEP::cm;

  //if nu parent is a muon then find muon parent info
  if ((parent_name=="mu+" || parent_name=="mu-") && NuParentTrack->GetParentID()!=0)
    {
      G4int mupar = NuParentTrack->GetParentID();
      NumiTrajectory* MuParentTrack = GetParentTrajectory(mupar);
      G4int nopoint_mupar = MuParentTrack->GetPointEntries();
      G4ThreeVector muparp = MuParentTrack->GetMomentum(nopoint_mupar-1);
      G4double muparm = MuParentTrack->GetMass();
      g4data->muparpx = muparp[0]/CLHEP::GeV; // vector of hadron parent of muon
      g4data->muparpy = muparp[1]/CLHEP::GeV; //
      g4data->muparpz = muparp[2]/CLHEP::GeV;
      g4data->mupare = (sqrt(muparp*muparp+muparm*muparm))/CLHEP::GeV;
    }
  else
    {
      g4data->muparpx = -999999.;
      g4data->muparpy = -999999.;
      g4data->muparpz = -999999.;
      g4data->mupare = -999999.;
    }

  g4data->Necm = enuzr/CLHEP::GeV; // Neutrino energy in parent rest frame
  NumiTrackInformation* info = (NumiTrackInformation*)(track.GetUserInformation());
  g4data->Nimpwt = info->GetNImpWt();  // Importance weight
  g4data->tgen = info->GetTgen()-1;

  g4data->xpoint = 0.;  // x, y, z of parent at user selected vol
  g4data->xpoint = 0.;
  g4data->xpoint = 0.;

  /*
        tgen is is the "generation" number
        of the particle that makes it out of the target. Beam protons have
        tgen=1, any particle produced by a p-C interaction would have tgen=2,
        particles produced from interactions of those products have tgen=3 etc.
        etc. until the cascade exiting the target core.
  */

  //if not using external ntuple then need to find the particle that exited the target

  G4int tar_pdg = 0; //auxiliar variable to convert to PDG code
  G4int tar_trackId = -1;
  if(!NumiData->useFlukaInput && !NumiData->useMarsInput)
    {
        G4bool findTarget = false;
        G4ThreeVector ParticleMomentum = G4ThreeVector(-999999,-999999,-999999);
        G4ThreeVector ParticlePosition = G4ThreeVector(-999999,-999999,-999999);
        NumiTrajectory* PParentTrack = GetParentTrajectory(track.GetParentID());
        G4int tptype = NumiParticleCode::AsInt(NumiParticleCode::StringToEnum(PParentTrack->GetParticleName()));
        particleID = PParentTrack->GetTrackID();

        while (!findTarget && particleID!=1){
            G4int numberOfPoints = PParentTrack->GetPointEntries();
            for (G4int ii=0;ii<numberOfPoints-1; ++ii){
                G4String lastVolName = PParentTrack->GetPreStepVolumeName(ii);
                G4String nextVolName = PParentTrack->GetPreStepVolumeName(ii+1);
                if (theSteppingAction->EscapingTarget(lastVolName, nextVolName)) {
                    // tv_ and tp_ are equal to position and
                    // momentum of the particle exiting the target (actually shell around target)
                    ParticleMomentum = PParentTrack->GetMomentum(ii);
                    ParticlePosition = PParentTrack->GetPoint(ii)->GetPosition();
                    tptype = NumiParticleCode::AsInt(NumiParticleCode::StringToEnum(PParentTrack->GetParticleName()));
                    tar_pdg = PParentTrack->GetPDGEncoding();
                    tar_trackId = PParentTrack->GetTrackID();
                    findTarget = true;
//                  std::cerr << " NumiAnalysis::FillNeutrinoNtuple Find Target is true Position "
//                            << ParticlePosition << " momentum " << ParticleMomentum << std::endl;
                  }
            }
            PParentTrack = GetParentTrajectory(PParentTrack->GetParentID());
            particleID = PParentTrack->GetTrackID();
        }
        g4data->tvx = ParticlePosition[0]/CLHEP::cm;
        g4data->tvy = ParticlePosition[1]/CLHEP::cm;
        g4data->tvz = ParticlePosition[2]/CLHEP::cm;
        g4data->tpx = ParticleMomentum[0]/CLHEP::GeV;
        g4data->tpy = ParticleMomentum[1]/CLHEP::GeV;
        g4data->tpz = ParticleMomentum[2]/CLHEP::GeV;
        g4data->tptype = tptype;
    }
  else          // using external ntuple, so set these to whatever comes from that ntuple
  {
      G4ThreeVector ParticlePosition = NPGA->GetParticlePosition();
      g4data->tvx = ParticlePosition[0]/CLHEP::cm;
      g4data->tvy = ParticlePosition[1]/CLHEP::cm;
      g4data->tvz = ParticlePosition[2]/CLHEP::cm;
      G4ThreeVector ParticleMomentum = NPGA->GetParticleMomentum();
      g4data->tpx = ParticleMomentum[0]/CLHEP::GeV;
      g4data->tpy = ParticleMomentum[1]/CLHEP::GeV;
      g4data->tpz = ParticleMomentum[2]/CLHEP::GeV;
      g4data->tptype = NPGA->GetParticleType();
  }

      //set all trk_ & trkp_ to -999999
  for (G4int ii=0; ii<10; ++ii){
      g4data->trkx[ii] = -999999.;
      g4data->trky[ii] = -999999.;
      g4data->trkz[ii] = -999999.;
      g4data->trkpx[ii] = -999999.;
      g4data->trkpy[ii] = -999999.;
      g4data->trkpz[ii] = -999999.;
  }
      // Neutrino data at different points
      // need neutrino parent info to be filled in g4data by this point

  for (G4int ii=0; ii<NumiData->nNear; ++ii){          // near detector
      g4data->NdxdzNear[ii] = (x-NumiData->xdet_near[ii])/(z-NumiData->zdet_near[ii]);
      g4data->NdydzNear[ii] = (y-NumiData->ydet_near[ii])/(z-NumiData->zdet_near[ii]);

      NumiNuWeight nuwgh;
      G4double nu_wght;
      G4double nu_energy;
      std::vector<double> r_det;
      r_det.push_back(NumiData->xdet_near[ii]/CLHEP::cm);
      r_det.push_back(NumiData->ydet_near[ii]/CLHEP::cm);
      r_det.push_back(NumiData->zdet_near[ii]/CLHEP::cm);
      nuwgh.GetWeight(g4data, r_det,nu_wght,nu_energy);
      g4data->NenergyN[ii] = nu_energy; //in GeV
      g4data->NWtNear[ii]  = nu_wght;
  }

  for (G4int ii=0; ii<NumiData->nFar; ++ii){         // far detector
      g4data->NdxdzFar[ii] = (x-NumiData->xdet_far[ii])/(z-NumiData->zdet_far[ii]);
      g4data->NdydzFar[ii] = (y-NumiData->ydet_far[ii])/(z-NumiData->zdet_far[ii]);

      NumiNuWeight nuwgh;
      G4double nu_wght;
      G4double nu_energy;
      std::vector<double> r_det;
      r_det.push_back(NumiData->xdet_far[ii]/CLHEP::cm);
      r_det.push_back(NumiData->ydet_far[ii]/CLHEP::cm);
      r_det.push_back(NumiData->zdet_far[ii]/CLHEP::cm);
      nuwgh.GetWeight(g4data, r_det,nu_wght,nu_energy);
      g4data->NenergyF[ii] = nu_energy; //in GeV
      g4data->NWtFar[ii]   = nu_wght;

  }

      //if(parentID!=0){
  G4ThreeVector ParentMomentum;
  G4ThreeVector ParentPosition;

  G4String sconf = NumiData->GetBeamConfig();
  G4bool is_me   = sconf.contains("me") || sconf.contains("ME");
  //std::cout<<"=> (check) Beam Configuration: "<< (NumiData->GetBeamConfig()) <<" "<< is_me <<std::endl;

  for (G4int ii=0; ii<point_no-1; ++ii){
      ParentMomentum = NuParentTrack->GetMomentum(ii);
      ParentPosition = (NuParentTrack->GetPoint(ii)->GetPosition()/CLHEP::m)*CLHEP::m;

      G4String postvolname = "";
      G4String prevolname = NuParentTrack->GetPreStepVolumeName(ii);
      if (ii<point_no-2) postvolname = NuParentTrack->GetPreStepVolumeName(ii+1);

          // parent created inside target
      if ((prevolname.contains("TGT")||prevolname.contains("Budal")) && ii==0){
          g4data->trkx[0] = ParentPosition[0]/CLHEP::cm;
          g4data->trky[0] = ParentPosition[1]/CLHEP::cm;
          g4data->trkz[0] = ParentPosition[2]/CLHEP::cm;
          g4data->trkpx[0] = ParentMomentum[0]/CLHEP::GeV;
          g4data->trkpy[0] = ParentMomentum[1]/CLHEP::GeV;
          g4data->trkpz[0] = ParentMomentum[2]/CLHEP::GeV;}
          //parent at exits target
      if (theSteppingAction->EscapingTarget(prevolname, postvolname)) {
            g4data->trkx[1] = ParentPosition[0]/CLHEP::cm;
            g4data->trky[1] = ParentPosition[1]/CLHEP::cm;
            g4data->trkz[1] = ParentPosition[2]/CLHEP::cm;
            g4data->trkpx[1] = ParentMomentum[0]/CLHEP::GeV;
            g4data->trkpy[1] = ParentMomentum[1]/CLHEP::GeV;
            g4data->trkpz[1] = ParentMomentum[2]/CLHEP::GeV;
//          std::cerr << " NumiAnalysis::FillNeutrinoNtuple filling point 1 ... " << ParentPosition << std::endl;
         }

          //enter horn1
      if (prevolname.contains("TGAR") && postvolname.contains("Horn1")){
          g4data->trkx[2] = ParentPosition[0]/CLHEP::cm;
          g4data->trky[2] = ParentPosition[1]/CLHEP::cm;
          g4data->trkz[2] = ParentPosition[2]/CLHEP::cm;
          g4data->trkpx[2] = ParentMomentum[0]/CLHEP::GeV;
          g4data->trkpy[2] = ParentMomentum[1]/CLHEP::GeV;
          g4data->trkpz[2] = ParentMomentum[2]/CLHEP::GeV;
      }
          //exit horn1
      if (prevolname.contains("Horn1") && postvolname.contains("TGAR")){
          g4data->trkx[3] = ParentPosition[0]/CLHEP::cm;
          g4data->trky[3] = ParentPosition[1]/CLHEP::cm;
          g4data->trkz[3] = ParentPosition[2]/CLHEP::cm;
          g4data->trkpx[3] = ParentMomentum[0]/CLHEP::GeV;
          g4data->trkpy[3] = ParentMomentum[1]/CLHEP::GeV;
          g4data->trkpz[3] = ParentMomentum[2]/CLHEP::GeV;
      }
          //enter horn2
      if (prevolname.contains("TGAR") && postvolname.contains("Horn2")){
          g4data->trkx[4] = ParentPosition[0]/CLHEP::cm;
          g4data->trky[4] = ParentPosition[1]/CLHEP::cm;
          g4data->trkz[4] = ParentPosition[2]/CLHEP::cm;
          g4data->trkpx[4] = ParentMomentum[0]/CLHEP::GeV;
          g4data->trkpy[4] = ParentMomentum[1]/CLHEP::GeV;
          g4data->trkpz[4] = ParentMomentum[2]/CLHEP::GeV;
      }
          //exit horn2
    if (prevolname.contains("Horn2") && postvolname.contains("TGAR")){
        g4data->trkx[5] = ParentPosition[0]/CLHEP::cm;
        g4data->trky[5] = ParentPosition[1]/CLHEP::cm;
        g4data->trkz[5] = ParentPosition[2]/CLHEP::cm;
        g4data->trkpx[5] = ParentMomentum[0]/CLHEP::GeV;
        g4data->trkpy[5] = ParentMomentum[1]/CLHEP::GeV;
        g4data->trkpz[5] = ParentMomentum[2]/CLHEP::GeV;
    }
        //enter decay pipe
    if (prevolname.contains("DVOL") && (postvolname.contains("UpWn"))){
        g4data->trkx[6] = ParentPosition[0]/CLHEP::cm;
        g4data->trky[6] = ParentPosition[1]/CLHEP::cm;
        g4data->trkz[6] = ParentPosition[2]/CLHEP::cm;
        g4data->trkpx[6] = ParentMomentum[0]/CLHEP::GeV;
        g4data->trkpy[6] = ParentMomentum[1]/CLHEP::GeV;
        g4data->trkpz[6] = ParentMomentum[2]/CLHEP::GeV;}



    // check if the particle passes through the neck of the horn
    // if yes then set the trk_ to +999999
    // need to make this work for arbitrary horn position!!
    //LE:
    if ((ParentPosition[2]>0.&&ParentPosition[2]<3.*CLHEP::m)&&  // horn 1 position 0-3m
        (sqrt(ParentPosition[0]*ParentPosition[0]+ParentPosition[1]*ParentPosition[1])<5.*CLHEP::cm)&&
        (g4data->trkx[2]==-999999. || g4data->trkx[2]==999999.))
    {
        g4data->trkx[2] = 999999.;
        g4data->trky[2] = 999999.;
        g4data->trkz[2] = 999999.;
    }
    //Different Horn2 positions between LE and ME (Leo, Sept 13, 2020)
    if(!is_me){
      if ((ParentPosition[2]>10.*CLHEP::m&&ParentPosition[2]<13.*CLHEP::m)&&  //horn 2 position 10-13m
          (sqrt(ParentPosition[0]*ParentPosition[0]+ParentPosition[1]*ParentPosition[1])<5.*CLHEP::cm)&&
          (g4data->trkx[4]==-999999. || g4data->trkx[4]==999999.))
        {
          g4data->trkx[4] = 999999.;
          g4data->trky[4] = 999999.;
          g4data->trkz[4] = 999999.;
        }
    }
    //ME
    else if(is_me){
      if ((ParentPosition[2]>19.*CLHEP::m&&ParentPosition[2]<22.*CLHEP::m)&&  //horn 2 position 19-22m
          (sqrt(ParentPosition[0]*ParentPosition[0]+ParentPosition[1]*ParentPosition[1])<5.*CLHEP::cm)&&
          (g4data->trkx[4]==-999999. || g4data->trkx[4]==999999.))
        {
          g4data->trkx[4] = 999999.;
          g4data->trky[4] = 999999.;
          g4data->trkz[4] = 999999.;
        }
    }
  }

  ParentMomentum = NuParentTrack->GetMomentum(point_no-1);
  ParentPosition = (NuParentTrack->GetPoint(point_no-1)->GetPosition()/CLHEP::m)*CLHEP::m;
  g4data->trkx[7] = ParentPosition[0]/CLHEP::cm;
  g4data->trky[7] = ParentPosition[1]/CLHEP::cm;
  g4data->trkz[7] = ParentPosition[2]/CLHEP::cm;
  g4data->trkpx[7] = ParentMomentum[0]/CLHEP::GeV;
  g4data->trkpy[7] = ParentMomentum[1]/CLHEP::GeV;
  g4data->trkpz[7] = ParentMomentum[2]/CLHEP::GeV;
      // }

      // Reset the neutrino history ntuple variables
  g4data->ntrajectory = -1;
  g4data->overflow = false;
  for (std::size_t index = 0; index < maxGen; ++index) {
      g4data->pdg[index] = 0;
      g4data->trackId[index] = -1;
      g4data->parentId[index] = -1;
      g4data->startx[index]  = 999.999;
      g4data->starty[index]  = 999.999;
      g4data->startz[index]  = 999.999;
      g4data->startpx[index] = 999.999;
      g4data->startpy[index] = 999.999;
      g4data->startpz[index] = 999.999;

      g4data->stopx[index]  = 999.999;
      g4data->stopy[index]  = 999.999;
      g4data->stopz[index]  = 999.999;
      g4data->stoppx[index] = 999.999;
      g4data->stoppy[index] = 999.999;
      g4data->stoppz[index] = 999.999;

      g4data->pprodpx[index] = 999.999;
      g4data->pprodpy[index] = 999.999;
      g4data->pprodpz[index] = 999.999;

      g4data->proc[index] = "NotDefined";
      g4data->ivol[index] = "NotDefined";
      g4data->fvol[index] = "NotDefined";
  }


      // Save neutrino history to the ntuple. If the depth of the history is more than the
      // allocated memory in the arrays, then the variable 'overflow' will be set to 'TRUE'.
      // Entries with overflow set to true should be excluded from analysis. Currently at
      // maxGen = 12, this happens approximately every 10,000 neutrinos. If this is unacceptable,
      // the maxGen value should be increased to 13 or higher.
  g4data->ntrajectory = history.size();
  if (history.size()>maxGen) g4data->overflow = true;
  int idx_tar_in_chain = -1;

      // Create a temporary vector so that it can be reversed.
  std::vector<G4VTrajectory*> tmpHistory(history);
  std::reverse(tmpHistory.begin(),tmpHistory.end());
  for (std::size_t index = 0; index < std::min(maxGen,history.size()); ++index) {
      NumiTrajectory* traj = dynamic_cast<NumiTrajectory*>(tmpHistory.at(index));
      g4data->pdg[index] = traj->GetPDGEncoding();
      g4data->trackId[index] = traj->GetTrackID();
      g4data->parentId[index] = traj->GetParentID();
      g4data->startx[index]  = traj->GetPoint(0)->GetPosition().x();
      g4data->starty[index]  = traj->GetPoint(0)->GetPosition().y();
      g4data->startz[index]  = traj->GetPoint(0)->GetPosition().z();
      g4data->startpx[index] = traj->GetMomentum(0).x();
      g4data->startpy[index] = traj->GetMomentum(0).y();
      g4data->startpz[index] = traj->GetMomentum(0).z();

      const int lastPoint = traj->GetPointEntries()-1;
      g4data->stopx[index]  = traj->GetPoint(lastPoint)->GetPosition().x();
      g4data->stopy[index]  = traj->GetPoint(lastPoint)->GetPosition().y();
      g4data->stopz[index]  = traj->GetPoint(lastPoint)->GetPosition().z();
      g4data->stoppx[index] = traj->GetMomentum(lastPoint).x();
      g4data->stoppy[index] = traj->GetMomentum(lastPoint).y();
      g4data->stoppz[index] = traj->GetMomentum(lastPoint).z();

      g4data->pprodpx[index] = traj->GetParentMomentumAtThisProduction().x();
      g4data->pprodpy[index] = traj->GetParentMomentumAtThisProduction().y();
      g4data->pprodpz[index] = traj->GetParentMomentumAtThisProduction().z();

      g4data->proc[index] = traj->GetProcessName();
      g4data->ivol[index] = traj->GetPreStepVolumeName(0);
      g4data->fvol[index] = traj->GetPreStepVolumeName(lastPoint);

      //Searching for the index of the hadron that the target:
      //This is saved later in dk2nu:
      if(g4data->trackId[index] == tar_trackId){
        idx_tar_in_chain = int(index);
      }

  }

  //For now, I write it by hand but I need to figure out how to do it automatically
  //this values were extracted from NumiMaterials and checking the volumes definition
  //of IC, DPIP and DVOL. the fact_X is mass number over density ( gr/mole / g /cm3 )
  //and looking at parent(0), grand-parent(1), great-gran-parent(2)
  // (Leo Aliaga. Feb18, 2015)

  const int Nanc = NumiData->nGenAbs;
//  std::cerr << " NumiAnalysis fill ntuple,  Nanc " << Nanc << std::endl;

  G4double fact_Al = (26.98) / (2.7);
  G4double fact_Fe = (55.85) / (7.869990);
  G4double fact_He = (4.003) / (0.000145);
  G4double dist_IC1[3];
  G4double dist_IC2[3];
  G4double dist_DPIP[3];
  G4double dist_DVOL[3];

  if (Nanc > 3) {
    std::cerr << " NumiAnalysis::FillNeutrinoNtuple, Nanc greter than 3, Overwrite!... Quit here " << std::endl;
    exit(2);
  }

  //initializing:
  for(int ii=0; ii<Nanc;ii++){
    dist_IC1[ii]  = -1.0*fact_Al;
    dist_IC2[ii]  = -1.0*fact_Al;
    dist_DPIP[ii] = -1.0*fact_Fe;
    dist_DVOL[ii] = -1.0*fact_He;
  }
  //Getting values:
  NumiTrajectory* tmp_traj;
  for(int ii=0;ii<Nanc;ii++){
    if(history.size()<=3 && ii==2)continue;
    tmp_traj = dynamic_cast<NumiTrajectory*>(tmpHistory.at(history.size()-(ii+2)));
    //If the Horn1 is alternate (See doc 10573), the name convention of the Horn1 IC volumes changes.
    //Ths part of the code checks if the horn1 is alternate. If not, it assumes it is the nominal (Horn 1
    //model we used in g4numi v5).
    double tmp_dist_h1 = 0;
    if( NumiData->GetHorn1IsAlternate() ){
      tmp_dist_h1 += GetDistanceInVolume(tmp_traj,"Horn1IOTransCont0_P");
      tmp_dist_h1 += GetDistanceInVolume(tmp_traj,"Horn1IOTransCont1_P");
      tmp_dist_h1 += GetDistanceInVolume(tmp_traj,"Horn1IOTransCont2_P");
      tmp_dist_h1 += GetDistanceInVolume(tmp_traj,"Horn1IOTransCont3_P");

      tmp_dist_h1 += GetDistanceInVolume(tmp_traj,"Horn1UpstrSubSect0WeldUpstr_P");

      tmp_dist_h1 += GetDistanceInVolume(tmp_traj,"Horn1UpstrSubSect0_P");
      tmp_dist_h1 += GetDistanceInVolume(tmp_traj,"Horn1UpstrSubSect1_P");
      tmp_dist_h1 += GetDistanceInVolume(tmp_traj,"Horn1UpstrSubSect2_P");

      tmp_dist_h1 += GetDistanceInVolume(tmp_traj,"Horn1UpstrSubSect1Weld0_P");

      tmp_dist_h1 += GetDistanceInVolume(tmp_traj,"Horn1ToNeckPartM0SubSect0_P");
      tmp_dist_h1 += GetDistanceInVolume(tmp_traj,"Horn1ToNeckPartM0SubSect1_P");
      tmp_dist_h1 += GetDistanceInVolume(tmp_traj,"Horn1ToNeckPartM0SubSect2_P");
      tmp_dist_h1 += GetDistanceInVolume(tmp_traj,"Horn1ToNeckPartM0SubSect3_P");

      tmp_dist_h1 += GetDistanceInVolume(tmp_traj,"Horn1Neck_P");

      tmp_dist_h1 += GetDistanceInVolume(tmp_traj,"Horn1DownstrPart1SubSect0_P");
      tmp_dist_h1 += GetDistanceInVolume(tmp_traj,"Horn1DownstrPart1SubSect1_P");
      tmp_dist_h1 += GetDistanceInVolume(tmp_traj,"Horn1DownstrPart1SubSect2_P");
      tmp_dist_h1 += GetDistanceInVolume(tmp_traj,"Horn1DownstrPart1SubSect3_P");
      tmp_dist_h1 += GetDistanceInVolume(tmp_traj,"Horn1DownstrPart1SubSect4_P");
      tmp_dist_h1 += GetDistanceInVolume(tmp_traj,"Horn1DownstrPart1SubSect5_P");

      tmp_dist_h1 += GetDistanceInVolume(tmp_traj,"Horn1DownstrPart1Weld1_P");
    }
    else{
      tmp_dist_h1 = GetDistanceInVolume(tmp_traj,"PHorn1IC");
    }
    dist_IC1[ii]  = tmp_dist_h1;
    dist_IC2[ii]  = GetDistanceInVolume(tmp_traj,"PHorn2IC");
    dist_DPIP[ii] = GetDistanceInVolume(tmp_traj,"DPIP");
    dist_DVOL[ii] = GetDistanceInVolume(tmp_traj,"DVOL");

 }


  /////////////////////////////////////////////////////////////////////////
  //In this approach, I am not deleting any g4data variables.
  //Just copying or accomodating values for dk2nu:

  //For conversions:
  G4double mm2cm   = 0.1;
  G4double MeV2GeV = 0.001;

  //cleanning everything:
  vec_traj.clear();
  vec_nuray.clear();
  vec_ancestor.clear();
  vec_int.clear();
 vec_dbl.clear();

  //1) Ancestries:

#ifdef USEMODGEANT4

  MinervaElementInter*  ei =  MinervaElementInter::getInstance();
  std::map<G4int,G4String> ei_map =  ei->GetInfo();
  std::map<G4int,G4String>::iterator ei_it;
  std::map<G4int,G4ThreeVector> pi_map =  ei->GetP();
  std::map<G4int,G4ThreeVector>::iterator pi_it;

#endif

  std::size_t ntrajectory = history.size();
  for (std::size_t index = 0; index < ntrajectory; ++index) {
      NumiTrajectory* traj = dynamic_cast<NumiTrajectory*>(tmpHistory.at(index));

      bsim::Ancestor tmp_ancestor;

      G4double startx = mm2cm*(traj->GetPoint(0)->GetPosition().x());
      G4double starty = mm2cm*(traj->GetPoint(0)->GetPosition().y());
      G4double startz = mm2cm*(traj->GetPoint(0)->GetPosition().z());
      G4double localt = traj->GetTimeStart();
      tmp_ancestor.SetStartXYZT(startx,starty,startz,localt);

      G4double startpx = MeV2GeV*(traj->GetMomentum(0).x());
      G4double startpy = MeV2GeV*(traj->GetMomentum(0).y());
      G4double startpz = MeV2GeV*(traj->GetMomentum(0).z());
      tmp_ancestor.SetStartP(startpx,startpy,startpz);

      tmp_ancestor.pdg = traj->GetPDGEncoding();

      const int lastPoint = traj->GetPointEntries()-1;

      G4double stoppx  = MeV2GeV*(traj->GetMomentum(lastPoint).x());
      G4double stoppy  = MeV2GeV*(traj->GetMomentum(lastPoint).y());
      G4double stoppz  = MeV2GeV*(traj->GetMomentum(lastPoint).z());
      tmp_ancestor.SetStopP(stoppx,stoppy,stoppz);

      G4double pprodpx = MeV2GeV*(traj->GetParentMomentumAtThisProduction().x());
      G4double pprodpy = MeV2GeV*(traj->GetParentMomentumAtThisProduction().y());
      G4double pprodpz = MeV2GeV*(traj->GetParentMomentumAtThisProduction().z());

      tmp_ancestor.proc    = traj->GetProcessName();
      tmp_ancestor.ivol    = traj->GetPreStepVolumeName(0);
      tmp_ancestor.imat    = "NotFillYet";
      G4int elemID = 0;

      const G4ThreeVector pol = traj->GetPolarization();
      tmp_ancestor.polx = pol.x();
      tmp_ancestor.poly = pol.y();
      tmp_ancestor.polz = pol.z();


#ifdef USEMODGEANT4

      //filling inter. elem. and momuntum info:
      ei_it = ei_map.find(traj->GetTrackID());
      pi_it = pi_map.find(traj->GetTrackID());

      if(ei_it != ei_map.end()){
        elemID = GetNucleus(ei_it->second);
        // if(elemID == 0) G4cout<< "NO FOUND IN THE TABLE "<< ei_it->second << G4endl;
      }
      else{
        //the final neutrino does not interact:
        elemID = 0;
      }

      if(pi_it != pi_map.end()){
        G4ThreeVector mom_inter = pi_it->second;
        pprodpx  = mom_inter.x()CLHEP::*meV2GeV;
        pprodpy  = mom_inter.y()CLHEP::*meV2GeV;
        pprodpz  = mom_inter.z()CLHEP::*meV2GeV;
      }

#else
      elemID = traj->GetPDGNucleus();

      G4ThreeVector pprodp = traj->GetParentMomentumAtThisProduction();
      if ( pprodp.z() == 0 ) {

        //// amit bashyal 2021-11-17 "fix"????
        //G4int np = traj->GetPointEntries();
        //if (traj->GetMomentum(np-2)[2] > 1e-100 ) pprodp = traj->GetMomentum(np-2);
        //else {
        //  G4ThreeVector pLast = traj->GetMomentum(np-1);
        //  pprodp = pLast;
        //}

        // alternative to amit's "fix"
        // ParentMomentumAtThisProductions z was exactly zero
        // ... look for a non-zero case of GetMomentum
        // starting at the end, and moving backwards
        G4int np = traj->GetPointEntries();
        G4int ip;
        for (ip=np-1; ip>=0; --ip) {
          pprodp = traj->GetMomentum(ip);
          if ( pprodp.z() > 1e-100 ) break;
        }

      }
      pprodpx = MeV2GeV*pprodp.x();
      pprodpy = MeV2GeV*pprodp.y();
      pprodpz = MeV2GeV*pprodp.z();

#endif

      tmp_ancestor.nucleus = elemID;
      tmp_ancestor.SetPProdP(pprodpx,pprodpy,pprodpz);

      vec_ancestor.push_back(tmp_ancestor);
  }
  /////

  //2) TGT Exit:

  this_tgtexit.tvx = g4data->tvx;
  this_tgtexit.tvy = g4data->tvy;
  this_tgtexit.tvz = g4data->tvz;
  this_tgtexit.tpx = g4data->tpx;
  this_tgtexit.tpy = g4data->tpy;
  this_tgtexit.tpz = g4data->tpz;
  this_tgtexit.tptype = tar_pdg;
  this_tgtexit.tgen   =g4data->tgen;

  //2.1)Storing integer values:

  //-) Storing the index of tar
  vec_int.push_back(idx_tar_in_chain);

  //-) We need to see how to do it better. We already store the playlsit in the meta tree
  // but GENIE is not propagating dkmeta (just dk2nu) and then this playlist ID is going
  // to know the playlist used.
  G4String playlist  = NumiData->GetPlaylist();
  int plID = -1;
  if(playlist=="downstream")plID = 0;
  if(playlist=="minerva1")  plID = 1;
  if(playlist=="minerva2")  plID = 2;
  if(playlist=="minerva3")  plID = 3;
  if(playlist=="minerva4")  plID = 4;
  if(playlist=="minerva5")  plID = 5;
  if(playlist=="minerva6")  plID = 6;
  if(playlist=="minerva7")  plID = 7;
  if(playlist=="minerva8")  plID = 8;
  if(playlist=="minerva9")  plID = 9;
  if(playlist=="minerva10") plID = 10;
  if(playlist=="minerva11") plID = 11;
  if(playlist=="minerva12") plID = 12;
  if(playlist=="minerva13") plID = 13;
  vec_int.push_back(plID);

  //2.2)Storing the distance by density of the particle in IC:
  for(int ii=0;ii<Nanc;ii++)vec_dbl.push_back(dist_IC1[ii]/fact_Al);
  for(int ii=0;ii<Nanc;ii++)vec_dbl.push_back(dist_IC2[ii]/fact_Al);
  for(int ii=0;ii<Nanc;ii++)vec_dbl.push_back(dist_DPIP[ii]/fact_Fe);
  for(int ii=0;ii<Nanc;ii++)vec_dbl.push_back(dist_DVOL[ii]/fact_He);

  /////


  //3) Trajectories (ray tracing):
  for (G4int ii=0; ii<10; ++ii){
    bsim::Traj tmp_traj;
    tmp_traj.trkx  = g4data->trkx[ii];
    tmp_traj.trky  = g4data->trky[ii];
    tmp_traj.trkz  = g4data->trkz[ii];
    tmp_traj.trkpx = g4data->trkpx[ii];
    tmp_traj.trkpy = g4data->trkpy[ii];
    tmp_traj.trkpz = g4data->trkpz[ii];
    vec_traj.push_back(tmp_traj);
  }
  ///////

  //4) Decay:
  this_decay.norig  = Norig;
  this_decay.ndecay = NuParentTrack->GetDecayCode();
  this_decay.ntype = particleType->GetPDGEncoding();
  this_decay.vx = x/CLHEP::cm;
  this_decay.vy = y/CLHEP::cm;
  this_decay.vz = z/CLHEP::cm;
  this_decay.pdpx     = g4data->pdPx;
  this_decay.pdpy     = g4data->pdPy;
  this_decay.pdpz     = g4data->pdPz;
  this_decay.ppdxdz   = ParentMomentumProduction[0]/ParentMomentumProduction[2];
  this_decay.ppdydz   = ParentMomentumProduction[1]/ParentMomentumProduction[2];
  this_decay.pppz     = ParentMomentumProduction[2]/CLHEP::GeV;
  this_decay.ppenergy = g4data->ppenergy;
  this_decay.ppmedium = int(g4data->ppmedium);
  this_decay.ptype    = NuParentTrack->GetPDGEncoding();
  this_decay.muparpx  = g4data->muparpx;
  this_decay.muparpy  = g4data->muparpy;
  this_decay.muparpz  = g4data->muparpz;
  this_decay.mupare   = g4data->mupare;
  this_decay.necm     = g4data->Necm;
  this_decay.nimpwt   = g4data->Nimpwt;
  //////

  //5) NuRay:

  // Random decay:
  G4double RdecPx = NuMomentum[0]/CLHEP::GeV;
  G4double RdecPy = NuMomentum[1]/CLHEP::GeV;
  G4double RdecPz = NuMomentum[2]/CLHEP::GeV;
  G4double RdecE  = track.GetTotalEnergy()/CLHEP::GeV;
  bsim::NuRay tmp_nuray_random(RdecPx,RdecPy,RdecPz,RdecE,1.0);
  vec_nuray.push_back(tmp_nuray_random);

  //calculating again...
  //I will optimize this.

  //for near detectors:
  for(G4int ii=0; ii<NumiData->nNear; ++ii){
    NumiNuWeight nuwgh;
    G4double nu_wght;
    G4double nu_energy;
    std::vector<double> r_det;
    r_det.push_back(NumiData->xdet_near[ii]/CLHEP::cm);
    r_det.push_back(NumiData->ydet_near[ii]/CLHEP::cm);
    r_det.push_back(NumiData->zdet_near[ii]/CLHEP::cm);
    nuwgh.GetWeight(g4data, r_det,nu_wght,nu_energy);
    G4double mom_nu[3];
    double rad = sqrt((r_det[0]- g4data->Vx)*(r_det[0]- g4data->Vx) +
                      (r_det[1]- g4data->Vy)*(r_det[1]- g4data->Vy) +
                      (r_det[2]- g4data->Vz)*(r_det[2]- g4data->Vz));
    mom_nu[0] = (r_det[0]- g4data->Vx) * nu_energy / rad;
    mom_nu[1] = (r_det[1]- g4data->Vy) * nu_energy / rad;
    mom_nu[2] = (r_det[2]- g4data->Vz) * nu_energy / rad;
    bsim::NuRay tmp_nuray(mom_nu[0],mom_nu[1],mom_nu[2],nu_energy,nu_wght);
    vec_nuray.push_back(tmp_nuray);
    if ((energyBinSimpleHistoMinerva > 0.) &&
        (ii == 1) && (nu_energy > 0.)) { // 2nd detector is centered on axis. We
      size_t iBin = nu_energy/energyBinSimpleHistoMinerva;
      if ((iBin < MinervaNuMuHisto.size()) &&  (particleType->GetPDGEncoding() == 14))
        MinervaNuMuHisto[iBin] += nu_wght*g4data->Nimpwt;
      if ((iBin < MinervaNuMuBarHisto.size()) &&  (particleType->GetPDGEncoding() == -14))
        MinervaNuMuBarHisto[iBin] += nu_wght*g4data->Nimpwt;
    }
    if ((energyBinSimpleHistoMinerva > 0.) &&
        (ii == 2) && (nu_energy > 0.)) { // Nova Near Detector off axis..
      size_t iBin = nu_energy/energyBinSimpleHistoMinerva;
      if ((iBin < NovaNearNuMuHisto.size()) &&  (particleType->GetPDGEncoding() == 14))
        NovaNearNuMuHisto[iBin] += nu_wght*g4data->Nimpwt;
      if ((iBin < NovaNearNuMuBarHisto.size()) &&  (particleType->GetPDGEncoding() == -14))
        NovaNearNuMuBarHisto[iBin] += nu_wght*g4data->Nimpwt;
    }
  }
  //for far detectors:
  for(G4int ii=0; ii<NumiData->nFar; ++ii){
      NumiNuWeight nuwgh;
      G4double nu_wght;
      G4double nu_energy;
      std::vector<double> r_det;
      r_det.push_back(NumiData->xdet_far[ii]/CLHEP::cm);
      r_det.push_back(NumiData->ydet_far[ii]/CLHEP::cm);
      r_det.push_back(NumiData->zdet_far[ii]/CLHEP::cm);
      nuwgh.GetWeight(g4data, r_det,nu_wght,nu_energy);
      double rad = sqrt((r_det[0]- g4data->Vx)*(r_det[0]- g4data->Vx) +
                        (r_det[1]- g4data->Vy)*(r_det[1]- g4data->Vy) +
                        (r_det[2]- g4data->Vz)*(r_det[2]- g4data->Vz));
      G4double mom_nu[3];
      mom_nu[0] = (r_det[0]- g4data->Vx) * nu_energy / rad;
      mom_nu[1] = (r_det[1]- g4data->Vy) * nu_energy / rad;
      mom_nu[2] = (r_det[2]- g4data->Vz) * nu_energy / rad;
      bsim::NuRay tmp_nuray(mom_nu[0],mom_nu[1],mom_nu[2],nu_energy,nu_wght);
      vec_nuray.push_back(tmp_nuray);
      if ((energyBinSimpleHistoMinerva > 0.) &&
          (ii == 1) && (nu_energy > 0.)) { // Nova Near Detector off axis..
        size_t iBin = nu_energy/0.1;
        if ((iBin < NovaFarNuMuHisto.size()) &&  (particleType->GetPDGEncoding() == 14))
          NovaFarNuMuHisto[iBin] += nu_wght*g4data->Nimpwt;
        if ((iBin < NovaFarNuMuBarHisto.size()) &&  (particleType->GetPDGEncoding() == -14))
          NovaFarNuMuBarHisto[iBin] += nu_wght*g4data->Nimpwt;
      }
  }
  ////
  //6) Others:

  //calculating the job number:
  G4String namentp = (NumiData->nuNtupleName);
  G4int namesize = (NumiData->nuNtupleName).rfind('_');
  namentp.erase(0,namesize+1);
  istringstream buffer(namentp);
  int valjob;
  buffer >> valjob;
  this_dk2nu->job = valjob;

  this_dk2nu->potnum = g4data->evtno;
  this_dk2nu->ppvx = g4data->ppvx;
  this_dk2nu->ppvy = g4data->ppvy;
  this_dk2nu->ppvz = g4data->ppvz;

  //Placing in dk2nu:
  this_dk2nu->ancestor = vec_ancestor;
  this_dk2nu->tgtexit = this_tgtexit;
  this_dk2nu->traj = vec_traj;
  this_dk2nu->decay = this_decay;
  this_dk2nu->nuray = vec_nuray;
  this_dk2nu->vint = vec_int;  //Storing the vint
  this_dk2nu->vdbl = vec_dbl; //Storing the vdbl

  /////////////////////////////////////////////////////////////////////////


  tree->Fill();


      // Write to file
  if (NumiData->createASCII)
  {
      std::ofstream asciiFile(asciiFileName, std::ios::app);
      if(asciiFile.is_open())
      {
          asciiFile << g4data->Ntype<< " " << g4data->Nenergy << " " << g4data->NenergyN[0] << " " << g4data->NWtNear[0];
          asciiFile << " " << g4data->NenergyF[0] << " " << g4data->NWtFar[0] <<" "<<g4data->Nimpwt<< G4endl;
          asciiFile.close();
      }
  }


}

//------------------------------------------------------------------------------------
void NumiAnalysis::FillBXDRAW(const G4Step* aStep) {

  if (!NumiData->createBXDRAW) return;
  // Based on Fortran code in mgdraw.f in *_fluka/for/
  G4int JTRACK, LTRACK;
  G4int TMREG, TNEWREG;
  G4int NCASE;
  G4double WTRACK;
  G4double ETRACK, AM;
  G4double XSCO, YSCO, RSCO, ZSCO;
  G4double CXTRCK, CYTRCK, CZTRCK;

  G4Track * theTrack = aStep->GetTrack();
  G4StepPoint* startPoint = aStep->GetPreStepPoint();
  G4StepPoint* endPoint = aStep->GetPostStepPoint();


  G4VPhysicalVolume *Start = startPoint->GetPhysicalVolume();
  G4VPhysicalVolume *End   = endPoint->GetPhysicalVolume();

  if (!End || !Start) // Left simulated region
        return;

  G4String StartName = Start->GetName();
  G4String EndName   = End->GetName();

  TMREG = 0;
  TNEWREG = 0;
  WTRACK = theTrack->GetWeight();

  if (StartName.contains("Horn1Box"))
        TMREG = 1;             // Horn 1
  else if (StartName.contains("Horn2Box"))
        TMREG = 2;             // Horn 2
  else if (StartName.contains("DVOL") )
        TMREG = 3;             // Decay Pipe - DVOL only

  if (EndName.contains("Horn1Box"))
        TNEWREG = 1;             // Horn 1
  else if (EndName.contains("Horn2Box"))
        TNEWREG = 2;             // Horn 2
  else if (EndName.contains("DVOL") )
        TNEWREG = 3;             // Decay Pipe - DVOL only

  ZSCO = theTrack->GetPosition().z()/CLHEP::cm;

  if ( TMREG != TNEWREG )  {                 // Change Regions
        if ( TMREG == 1 || TMREG == 2 ||       // Exit horns
                 (ZSCO > 4569  && ZSCO < 4571) || // Enter Decay Pipe
                 (ZSCO > 72237 && ZSCO < 72238)   // Exit Decay Pipe
         ) {

          ETRACK = theTrack->GetTotalEnergy()/CLHEP::GeV;
          G4ParticleDefinition * particleType = theTrack->GetDefinition();
          JTRACK = code[particleType->GetPDGEncoding()];
          if (JTRACK == 8 || JTRACK == 9)
                return;
          if (ETRACK > 0.5 ) { // 0.5 GeV Cut off
                NumiTrackInformation* info=(NumiTrackInformation*)(theTrack->GetUserInformation());
                if (info) {
                  LTRACK = info->GetTgen();
                }
                else {
                  LTRACK = 0;
                }

                NCASE = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();

                AM = particleType->GetPDGMass()/CLHEP::GeV;

                XSCO = theTrack->GetPosition().x()/CLHEP::cm;
                YSCO = theTrack->GetPosition().y()/CLHEP::cm;

                RSCO = sqrt(XSCO*XSCO + YSCO*YSCO);

                CXTRCK = theTrack->GetMomentumDirection().x();
                CYTRCK = theTrack->GetMomentumDirection().y();
                CZTRCK = theTrack->GetMomentumDirection().z();


                std::ofstream bxdrawFile(bxdrawFileName, std::ios::app);
                bxdrawFile << JTRACK << "  " << TMREG  << "  " << TNEWREG<< "  "  // Code, regions
                                   << NCASE  << "  " << WTRACK << "  "                    // Event, Track Weight
                                   << LTRACK << "  " << ETRACK << "  " << AM     << "  "  // Gen, E, M
                                   << XSCO   << "  " << YSCO   << "  " << RSCO   << "  " << ZSCO << "  "  // crossing point
                                   << CXTRCK << "  " << CYTRCK << "  " << CZTRCK << "  "  // direction cosines
                                   << G4endl;
                bxdrawFile.close();
          }
        }
  }
}

//------------------------------------------------------------------------------------
NumiTrajectory* NumiAnalysis::GetParentTrajectory(G4int parentID)
{

   if(NumiData->GetDebugLevel() == 10) { G4cout << "NumiAnalysis::GetParentTrajectory() called." << G4endl;}

  G4TrajectoryContainer* container =
    G4RunManager::GetRunManager()->GetCurrentEvent()->GetTrajectoryContainer();
  if(container==0) return 0;

  TrajectoryVector* vect = container->GetVector();
  G4VTrajectory* tr;
  G4int ii = 0;
  while (ii<G4int(vect->size())){
    tr = (*vect)[ii];
    NumiTrajectory* tr1 = (NumiTrajectory*)(tr);
    if(tr1->GetTrackID() == parentID) return tr1;
    ++ii;
  }

  return 0;
}
//=============================================================================
//-----------------------------------------------------------------------------
void NumiAnalysis::FillZpNtuple(const G4Track& track,Int_t zpnum)
{//---Jasmine added
  if(!NumiData->raytracing && !NumiData->createZpNtuple) return ;
  G4RunManager* pRunManager = G4RunManager::GetRunManager();

  /*  if(g4zpdata->evtno!=  pRunManager->GetCurrentEvent()->GetEventID())
        {
        g4zpdata->run=pRunManager->GetCurrentRun()->GetRunID();
        g4zpdata->evtno = pRunManager->GetCurrentEvent()->GetEventID();
        WriteZpNtuple();//store info and initialize
        }*/

  if(zpnum<(Int_t)NumiData->Zpoint.size()){

    G4ParticleDefinition* particleDefinition = track.GetDefinition();

    g4zpdata->xposatz= track.GetPosition()[0]/CLHEP::cm;
    g4zpdata->yposatz= track.GetPosition()[1]/CLHEP::cm;
    g4zpdata->zposatz= track.GetPosition()[2]/CLHEP::cm;
    g4zpdata->xmomatz= track.GetMomentum()[0]/CLHEP::GeV;
    g4zpdata->ymomatz= track.GetMomentum()[1]/CLHEP::GeV;
    g4zpdata->zmomatz= track.GetMomentum()[2]/CLHEP::GeV;
    g4zpdata->matilen= track.GetMaterial()->GetNuclearInterLength();
    g4zpdata->field= track.GetStepLength()/CLHEP::cm;
    g4zpdata->pathlength= track.GetTrackLength()/CLHEP::cm;
    g4zpdata->ptypeatz= track.GetParentID();

    NumiTrackInformation* info = (NumiTrackInformation*)(track.GetUserInformation());
    if (info) {
      g4zpdata->nimpwt = info->GetNImpWt();
      g4zpdata->tgen   = info->GetTgen();
    }

    g4zpdata->pidtype= NumiParticleCode::AsInt(NumiParticleCode::StringToEnum(particleDefinition->GetParticleName()));
    g4zpdata->zpoint=NumiData->Zpoint[zpnum]/CLHEP::cm ;
    if(zpnum < (Int_t)NumiData->Zplane_modes.size()){
      g4zpdata->zplane_mode = NumiData->Zplane_modes[zpnum];
      g4zpdata->zplane_index = zpnum;
    }

    g4zpdata->run=pRunManager->GetCurrentRun()->GetRunID();
    g4zpdata->evtno = pRunManager->GetCurrentEvent()->GetEventID();
    WriteZpNtuple();//store info and initialize

  }
}
void NumiAnalysis::WriteZpNtuple(){//Jasmine added

  zptree->Fill();
    g4zpdata->xposatz=-10000;
    g4zpdata->yposatz=-10000;
    g4zpdata->zposatz=-10000;
    g4zpdata->xmomatz=-10000;
    g4zpdata->ymomatz=-10000;
    g4zpdata->zmomatz=-10000;
    g4zpdata->matilen=-10000;
    g4zpdata->field =-10000;
    g4zpdata->pathlength=-10000;
    g4zpdata->zpoint=-10000;
    g4zpdata->ptypeatz=-10;
    g4zpdata->pidtype=-10;
    g4zpdata->nimpwt=-10000;
    g4zpdata->tgen=-10;
    g4zpdata->zplane_mode=-10;
    g4zpdata->zplane_index=-10;
    g4zpdata->run=-100;
    g4zpdata->evtno=-100;
}

//------------------------------------------------------------------------------------
G4int NumiAnalysis::ProvId(std::map<G4String, G4int>& ids, const G4String& name)
{
   std::map<G4String, G4int>::const_iterator it = ids.find(name);
   if(it != ids.end()) return it->second;
   G4int id = (G4int)ids.size() + 1;   // 0 is reserved for "unknown"
   ids[name] = id;
   return id;
}

//------------------------------------------------------------------------------------
// Remember where every track came from, so the ancestry of a meson can be
// reconstructed without keeping the (much heavier) trajectory objects.  Called
// from NumiTrackingAction::PreUserTrackingAction for every track, before the
// meson branch, so a meson's own entry is always in the map by the time
// FillMeson runs.
void NumiAnalysis::RecordTrack(const G4Track& track)
{
   if(!NumiData->createMesonProvenance) return;

   // one primary proton per event, always tracked first: use it as the reset
   if(track.GetTrackID() == 1) fTrackProv.clear();

   if(!fTgtMotherResolved){
      fTgtMotherResolved = true;
      G4PhysicalVolumeStore* pvStore = G4PhysicalVolumeStore::GetInstance();
      for(size_t ipv = 0; ipv < pvStore->size(); ++ipv){
         const G4String& pvName = (*pvStore)[ipv]->GetName();
         if(pvName == "TargetMotherVol" || pvName == "pvTargetMother"){
            fTgtMotherPV = (*pvStore)[ipv];
            break;
         }
      }
      G4cout << "NumiAnalysis: meson provenance on, target mother volume "
             << (fTgtMotherPV ? fTgtMotherPV->GetName() : G4String("NOT FOUND"))
             << G4endl;
   }

   TrackProv prov;
   prov.pdg      = track.GetDefinition()->GetPDGEncoding();
   prov.parentID = track.GetParentID();
   prov.tgen     = 1;
   prov.E        = track.GetTotalEnergy()/CLHEP::GeV;
   prov.x        = track.GetPosition()[0]/CLHEP::cm;
   prov.y        = track.GetPosition()[1]/CLHEP::cm;
   prov.z        = track.GetPosition()[2]/CLHEP::cm;

   const G4VProcess* creator = track.GetCreatorProcess();
   prov.proc = ProvId(fProcIds, creator ? creator->GetProcessName()
                                        : G4String("primary"));

   const G4VPhysicalVolume* pv = track.GetVolume();
   prov.vol = pv ? ProvId(fVolIds, pv->GetName()) : 0;

   // "inside the target" means anywhere inside the target mother volume --
   // graphite, casing and cooling water alike.  The vertex position is stored
   // too, so a tighter graphite-only cut can still be made offline.
   prov.inTgt = false;
   if(fTgtMotherPV){
      const G4VTouchable* touch = track.GetTouchable();
      if(touch){
         for(G4int d = 0; d <= touch->GetHistoryDepth(); ++d){
            if(touch->GetVolume(d) == fTgtMotherPV){ prov.inTgt = true; break; }
         }
      }
   }

   NumiTrackInformation* info = (NumiTrackInformation*)(track.GetUserInformation());
   if(info) prov.tgen = info->GetTgen();

   fTrackProv[track.GetTrackID()] = prov;
}

//------------------------------------------------------------------------------------
void NumiAnalysis::FillMeson(const G4Track& track)
{
  if(!NumiData->raytracing) return ;
  if(!mesontree) return;

  G4ParticleDefinition* particleDefinition = track.GetDefinition();

  meson_px = track.GetMomentum()[0]/CLHEP::GeV;
  meson_py = track.GetMomentum()[1]/CLHEP::GeV;
  meson_pz = track.GetMomentum()[2]/CLHEP::GeV;
  meson_E  = track.GetTotalEnergy()/CLHEP::GeV;
  meson_x  = track.GetPosition()[0]/CLHEP::cm;
  meson_y  = track.GetPosition()[1]/CLHEP::cm;
  meson_z  = track.GetPosition()[2]/CLHEP::cm;
  meson_t  = track.GetGlobalTime()/CLHEP::ns;

  meson_weight = 1.0;
  meson_tgen   = 1;

  NumiTrackInformation* info = (NumiTrackInformation*)(track.GetUserInformation());
  if (info) {
     meson_weight = info->GetNImpWt();
     meson_tgen   = info->GetTgen();
  }
  meson_pdg = particleDefinition->GetPDGEncoding();

  if(NumiData->createMesonProvenance) FillMesonProvenance(track);

  mesontree->Fill();
}

//------------------------------------------------------------------------------------
// Walk this meson's ancestry back towards the beam proton and record (a) its
// own creator process and production volume, (b) the same for its immediate
// parent, and (c) the "carrier": the last ancestor that was born inside the
// target and whose next-in-chain daughter was born outside it, i.e. the
// particle that actually carried the energy out of the target.  If no ancestor
// was born in the target at all, the carrier is the beam proton itself and the
// chain never had anything to do with target production -- that is the
// signature of a proton that made it through the target and showered
// downstream.
void NumiAnalysis::FillMesonProvenance(const G4Track& track)
{
  meson_proc = meson_vol = meson_intgt = 0;
  meson_ppdg = meson_pproc = meson_pvol = meson_pintgt = 0;
  meson_pE = meson_pvx = meson_pvy = meson_pvz = -1e6;
  meson_cflag = -1;
  meson_cpdg = meson_cproc = meson_cgen = 0;
  meson_cE = meson_cvx = meson_cvy = meson_cvz = -1e6;
  meson_ivx = meson_ivy = meson_ivz = -1e6;

  std::map<G4int, TrackProv>::const_iterator it =
     fTrackProv.find(track.GetTrackID());
  if(it == fTrackProv.end()) return;      // should not happen

  const TrackProv self = it->second;
  meson_proc  = self.proc;
  meson_vol   = self.vol;
  meson_intgt = self.inTgt ? 1 : 0;

  // the meson was made in the target: nothing carried it out
  if(self.inTgt){
     meson_cflag = 0;
     meson_ivx = self.x;  meson_ivy = self.y;  meson_ivz = self.z;
  }

  TrackProv child = self;
  G4bool haveParent = false;
  G4int  guard = 0;
  while(child.parentID != 0 && ++guard < 1000){
     std::map<G4int, TrackProv>::const_iterator ip = fTrackProv.find(child.parentID);
     if(ip == fTrackProv.end()) break;
     const TrackProv& parent = ip->second;

     if(!haveParent){
        haveParent = true;
        meson_ppdg   = parent.pdg;
        meson_pproc  = parent.proc;
        meson_pvol   = parent.vol;
        meson_pintgt = parent.inTgt ? 1 : 0;
        meson_pE     = parent.E;
        meson_pvx    = parent.x;
        meson_pvy    = parent.y;
        meson_pvz    = parent.z;
     }

     // first target -> outside step going up the chain is the last one that
     // happened in time, which is the one we want
     if(meson_cflag < 0 && parent.inTgt && !child.inTgt){
        meson_cflag = 1;
        meson_cpdg  = parent.pdg;
        meson_cproc = parent.proc;
        meson_cE    = parent.E;
        meson_cvx   = parent.x;
        meson_cvy   = parent.y;
        meson_cvz   = parent.z;
        meson_cgen  = self.tgen - parent.tgen;
        meson_ivx   = child.x;
        meson_ivy   = child.y;
        meson_ivz   = child.z;
     }

     if(parent.parentID == 0){
        // reached the beam proton without ever passing through the target
        if(meson_cflag < 0){
           meson_cflag = 2;
           meson_cpdg  = parent.pdg;
           meson_cproc = parent.proc;
           meson_cE    = parent.E;
           meson_cvx   = parent.x;
           meson_cvy   = parent.y;
           meson_cvz   = parent.z;
           meson_cgen  = self.tgen - parent.tgen;
           meson_ivx   = child.x;   // where the proton first interacted
           meson_ivy   = child.y;
           meson_ivz   = child.z;
        }
        break;
     }
     child = parent;
  }
}

//------------------------------------------------------------------------------------
void NumiAnalysis::FillTgtExit(const G4Track& track)
{
  if(!tgtexittree) return;

  G4ParticleDefinition* particleDefinition = track.GetDefinition();

  tgtexit_px = track.GetMomentum()[0]/CLHEP::GeV;
  tgtexit_py = track.GetMomentum()[1]/CLHEP::GeV;
  tgtexit_pz = track.GetMomentum()[2]/CLHEP::GeV;
  tgtexit_E  = track.GetTotalEnergy()/CLHEP::GeV;
  tgtexit_x  = track.GetPosition()[0]/CLHEP::cm;
  tgtexit_y  = track.GetPosition()[1]/CLHEP::cm;
  tgtexit_z  = track.GetPosition()[2]/CLHEP::cm;
  tgtexit_t  = track.GetGlobalTime()/CLHEP::ns;

  tgtexit_weight = 1.0;
  tgtexit_tgen   = 1;

  NumiTrackInformation* info = (NumiTrackInformation*)(track.GetUserInformation());
  if (info) {
     tgtexit_weight = info->GetNImpWt();
     tgtexit_tgen   = info->GetTgen();
  }
  tgtexit_pdg = particleDefinition->GetPDGEncoding();

  tgtexittree->Fill();
}

//------------------------------------------------------------------------------------
// EM cascade seeds.
//
// Processes whose secondaries are part of an EM cascade that is already being
// followed from its seed.  Any gamma/e+- made by something else starts a new
// cascade.  GammaGeneralProc is the Geant4 11 wrapper of conv/compt/phot.
G4bool NumiAnalysis::IsEMCascadeProcess(const G4String& procName)
{
   static const char* const kCascade[] = {
      "conv", "compt", "phot", "eBrem", "eIoni", "annihil", "msc",
      "Rayl", "ePairProd", "CoulombScat", "eCoulombScat", "GammaGeneralProc"
   };
   for(size_t i = 0; i < sizeof(kCascade)/sizeof(kCascade[0]); ++i)
      if(procName == kCascade[i]) return true;
   return false;
}

G4int NumiAnalysis::ESVolId(const G4VPhysicalVolume* pv)
{
   return pv ? ProvId(fESVolIds, pv->GetName()) : ProvId(fESVolIds, "OutOfWorld");
}

G4int NumiAnalysis::ESMatId(const G4Material* mat)
{
   return mat ? ProvId(fESMatIds, mat->GetName()) : 0;
}

void NumiAnalysis::RecordTrackPdg(const G4Track& track)
{
   if(!emseedtree) return;
   G4int evt = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
   if(evt != fESTrackPdgEvt){ fESTrackPdg.clear(); fESTrackPdgEvt = evt; }
   G4int id = track.GetTrackID();
   if(id < 0) return;
   if((size_t)id >= fESTrackPdg.size()) fESTrackPdg.resize(2*id + 64, 0);
   fESTrackPdg[id] = track.GetDefinition()->GetPDGEncoding();
}

void NumiAnalysis::FillEMSeedCommon(const G4Track& track)
{
   es_evt      = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
   es_trackID  = track.GetTrackID();
   es_parentID = track.GetParentID();
   es_pdg      = track.GetDefinition()->GetPDGEncoding();
   es_parent_pdg = 0;
   if(es_evt == fESTrackPdgEvt && es_parentID > 0 &&
      (size_t)es_parentID < fESTrackPdg.size())
      es_parent_pdg = fESTrackPdg[es_parentID];
   const G4VProcess* creator = track.GetCreatorProcess();
   es_cproc = ProvId(fESProcIds, creator ? creator->GetProcessName() : G4String("primary"));
   es_nimpwt = 1.0;
   es_tgen   = 1;
   NumiTrackInformation* info = (NumiTrackInformation*)(track.GetUserInformation());
   if(info){ es_nimpwt = info->GetNImpWt(); es_tgen = info->GetTgen(); }
}

// e+- seed: creation point only.  The "interaction" block repeats the
// creation point with iproc = 0, since an electron starts showering at birth.
void NumiAnalysis::FillEMSeedAtCreation(const G4Track& track)
{
   if(!emseedtree) return;
   FillEMSeedCommon(track);
   es_E  = track.GetTotalEnergy()/CLHEP::GeV;
   es_px = track.GetMomentum()[0]/CLHEP::GeV;
   es_py = track.GetMomentum()[1]/CLHEP::GeV;
   es_pz = track.GetMomentum()[2]/CLHEP::GeV;
   es_x  = track.GetPosition()[0]/CLHEP::cm;
   es_y  = track.GetPosition()[1]/CLHEP::cm;
   es_z  = track.GetPosition()[2]/CLHEP::cm;
   es_t  = track.GetGlobalTime()/CLHEP::ns;
   const G4VPhysicalVolume* pv = track.GetVolume();
   es_vol = ESVolId(pv);
   es_mat = (pv && pv->GetLogicalVolume()) ? ESMatId(pv->GetLogicalVolume()->GetMaterial()) : 0;
   es_iproc = 0;
   es_ix = es_x; es_iy = es_y; es_iz = es_z; es_it = es_t; es_iE = es_E;
   es_ivol = es_vol; es_imat = es_mat;
   emseedtree->Fill();
   ++fESNLepton;
}

// photon seed, first step: remember where it was born
void NumiAnalysis::CacheEMSeedCreation(const G4Track& track, const G4StepPoint* pre)
{
   if(!emseedtree || !pre) return;
   fESCacheTrackID = track.GetTrackID();
   fESCacheEvt = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
   fESCacheE  = pre->GetTotalEnergy()/CLHEP::GeV;
   fESCachePx = pre->GetMomentum()[0]/CLHEP::GeV;
   fESCachePy = pre->GetMomentum()[1]/CLHEP::GeV;
   fESCachePz = pre->GetMomentum()[2]/CLHEP::GeV;
   fESCacheX  = pre->GetPosition()[0]/CLHEP::cm;
   fESCacheY  = pre->GetPosition()[1]/CLHEP::cm;
   fESCacheZ  = pre->GetPosition()[2]/CLHEP::cm;
   fESCacheT  = pre->GetGlobalTime()/CLHEP::ns;
   fESCacheVol = ESVolId(pre->GetPhysicalVolume());
   fESCacheMat = ESMatId(pre->GetMaterial());
}

// photon seed, first interaction (or leaving the world)
void NumiAnalysis::FillEMSeedAtInteraction(const G4Track& track, const G4Step* step,
                                           const G4String& procName)
{
   if(!emseedtree) return;
   FillEMSeedCommon(track);
   if(fESCacheTrackID == es_trackID && fESCacheEvt == es_evt){
      es_E = fESCacheE; es_px = fESCachePx; es_py = fESCachePy; es_pz = fESCachePz;
      es_x = fESCacheX; es_y = fESCacheY; es_z = fESCacheZ; es_t = fESCacheT;
      es_vol = fESCacheVol; es_mat = fESCacheMat;
   } else {
      // should not happen: fall back on the vertex stored in the track
      ++fESNNoCache;
      G4double ke = track.GetVertexKineticEnergy();
      const G4ThreeVector& d = track.GetVertexMomentumDirection();
      es_E  = ke/CLHEP::GeV;
      es_px = ke*d[0]/CLHEP::GeV; es_py = ke*d[1]/CLHEP::GeV; es_pz = ke*d[2]/CLHEP::GeV;
      es_x  = track.GetVertexPosition()[0]/CLHEP::cm;
      es_y  = track.GetVertexPosition()[1]/CLHEP::cm;
      es_z  = track.GetVertexPosition()[2]/CLHEP::cm;
      es_t  = (track.GetGlobalTime() - track.GetLocalTime())/CLHEP::ns;
      es_vol = 0;
      const G4LogicalVolume* lv = track.GetLogicalVolumeAtVertex();
      es_mat = lv ? ESMatId(lv->GetMaterial()) : 0;
   }
   const G4StepPoint* pre  = step->GetPreStepPoint();
   const G4StepPoint* post = step->GetPostStepPoint();
   es_iproc = ProvId(fESProcIds, procName);
   es_ix = post->GetPosition()[0]/CLHEP::cm;
   es_iy = post->GetPosition()[1]/CLHEP::cm;
   es_iz = post->GetPosition()[2]/CLHEP::cm;
   es_it = post->GetGlobalTime()/CLHEP::ns;
   es_iE = pre->GetTotalEnergy()/CLHEP::GeV;
   es_ivol = ESVolId(pre->GetPhysicalVolume());
   es_imat = ESMatId(pre->GetMaterial());
   emseedtree->Fill();
   ++fESNPhoton;
}

// reason 0: KE below KillTrackingThreshold, 1: pz < 0.  The weight is the
// parent's importance weight, which is what the track would have carried on
// average had it survived the importance sampling.
void NumiAnalysis::CountKilledMeson(const G4Track& track, G4int reason)
{
   G4int pdg = track.GetDefinition()->GetPDGEncoding();
   G4int isp = (pdg == 111) ? 0 : (pdg == 221) ? 1 : -1;
   if(isp < 0 || reason < 0 || reason > 1) return;
   G4double w = 1.0;
   NumiTrackInformation* info = (NumiTrackInformation*)(track.GetUserInformation());
   if(info) w = info->GetNImpWt();
   fKilledN[isp][reason]  += 1;
   fKilledW[isp][reason]  += w;
   fKilledWE[isp][reason] += w*track.GetTotalEnergy()/CLHEP::GeV;
}

//------------------------------------------------------------------------------------
void NumiAnalysis::FillTarNtuple(const G4Track& track)
{//---Melissa added
  if(!NumiData->createTarNtuple) return ;

  G4RunManager* pRunManager = G4RunManager::GetRunManager();

    G4ParticleDefinition* particleDefinition = track.GetDefinition();

    NumiTrackInformation* info = (NumiTrackInformation*)(track.GetUserInformation());
    g4tardata->impwt = info->GetNImpWt();  // Importance weight

    g4tardata->tvx= track.GetPosition()[0]/CLHEP::cm;
    g4tardata->tvy= track.GetPosition()[1]/CLHEP::cm;
    g4tardata->tvz= track.GetPosition()[2]/CLHEP::cm;
    g4tardata->tpx= track.GetMomentum()[0]/CLHEP::GeV;
    g4tardata->tpy= track.GetMomentum()[1]/CLHEP::GeV;
    g4tardata->tpz= track.GetMomentum()[2]/CLHEP::GeV;

    g4tardata->tptype= NumiParticleCode::AsInt(NumiParticleCode::StringToEnum(particleDefinition->GetParticleName()));

    g4tardata->run=pRunManager->GetCurrentRun()->GetRunID();
    g4tardata->evtno = pRunManager->GetCurrentEvent()->GetEventID();
    WriteTarNtuple();//store info and initialize
}

//------------------------------------------------------------------------------------
void NumiAnalysis::WriteTarNtuple(){//Melissa added

  tartree->Fill();
    g4tardata->impwt=-10000;
    g4tardata->tvx=-10000;
    g4tardata->tvy=-10000;
    g4tardata->tvz=-10000;
    g4tardata->tpx=-10000;
    g4tardata->tpy=-10000;
    g4tardata->tpz=-10000;
    g4tardata->tptype=-10;
    g4tardata->run=-100;
    g4tardata->evtno=-100;
}

//----------------------------
G4int NumiAnalysis::GetNucleus(G4String nucl_name){

  //The new dk2nu store the nucleus element as integer.
  //We have the element name from our custumized geant4.
  //Form now, I created this map will all elements in g4numi:
  //I am contructing the pdg code based on the convention:
  //+-10LZZZAAAI, where L is the number of s quarks,
  //I is the isomer id (both set to zero).
  G4int id = 0;

  if(nucl_name == "Hydrogen")        return  1000010010;
  if(nucl_name == "Helium")          return  1000020040;
  if(nucl_name == "Carbon")          return  1000060120;
  if(nucl_name == "Nitrogen")        return  1000070140;
  if(nucl_name == "Oxygen")          return  1000080160;
  if(nucl_name == "Natrium")         return  1000110230;
  if(nucl_name == "Magnesium")       return  1000120240;
  if(nucl_name == "Aluminum")        return  1000130270;
  if(nucl_name == "Silicon")         return  1000140280;
  if(nucl_name == "Phosphorus")      return  1000150310;
  if(nucl_name == "Sulfur")          return  1000160320;
  if(nucl_name == "Argon")           return  1000180400;
  if(nucl_name == "Potassium")       return  1000190390;
  if(nucl_name == "Calcium")         return  1000200400;
  if(nucl_name == "Titanium")        return  1000220480;
  if(nucl_name == "Chromium")        return  1000240520;
  if(nucl_name == "Manganese")       return  1000250550;
  if(nucl_name == "Iron")            return  1000260560;
  if(nucl_name == "Nickel")          return  1000280590;
  if(nucl_name == "Copper")          return  1000290640;
  if(nucl_name == "Gallium")         return  1000310700;
  if(nucl_name == "Mercury")         return  1000802010;
  if(nucl_name == "Sodium")          return  1000110230;
  if(nucl_name == "Phospho")         return  1000150310;
  if(nucl_name == "Berylliu")        return  1000040090;
  if(nucl_name == "Berillium")       return  1000040090;
  if(nucl_name == "Target")          return  1000060120;
  if(nucl_name == "Lead")            return  1000822070;
  if(nucl_name == "SecMonitorHelium")return  1000020040;
  if(nucl_name == "TargetHelium")    return  1000020040;
  if(nucl_name == "HeGas")           return  1000020040;

  return id;

}

//----------------------------
G4double NumiAnalysis::GetDistanceInVolume(NumiTrajectory* wanted_traj, G4String wanted_vol){
  double dist_vol = 0;
  if(wanted_traj==0)return -1.;

  G4ThreeVector ParticlePos;
  G4int npoints = GetParentTrajectory(wanted_traj->GetTrackID())->GetPointEntries();

  G4ThreeVector tmp_ipos,tmp_fpos;
  G4ThreeVector tmp_disp;

  G4double tmp_dist = 0.0;
  G4bool enter_vol = false;
  G4bool exit_vol  = false;
  for(G4int ii=0; ii<npoints; ++ii){
    ParticlePos = (wanted_traj->GetPoint(ii)->GetPosition()/CLHEP::m)*CLHEP::m;
    G4String postvol = "";
    G4String prevol  = wanted_traj->GetPreStepVolumeName(ii);
    if(ii<npoints-1) postvol = wanted_traj->GetPreStepVolumeName(ii+1);

    G4bool vol_in  = ( (prevol != wanted_vol) && (postvol == wanted_vol) ) || ( ii==0 && prevol== wanted_vol);
    G4bool vol_out = (prevol == wanted_vol) && (postvol != wanted_vol);
    if(vol_in){
      enter_vol = true;
      exit_vol  = false;
      tmp_ipos = G4ThreeVector(ParticlePos[0]/CLHEP::cm,ParticlePos[1]/CLHEP::cm,ParticlePos[2]/CLHEP::cm);
      tmp_dist = 0.0;
    }
    if(enter_vol && !exit_vol){
      tmp_fpos = G4ThreeVector(ParticlePos[0]/CLHEP::cm,ParticlePos[1]/CLHEP::cm,ParticlePos[2]/CLHEP::cm);
      tmp_disp = tmp_fpos - tmp_ipos;
      tmp_dist += tmp_disp.mag();
      tmp_ipos = tmp_fpos;
    }
    if(enter_vol && vol_out){
      tmp_fpos  = G4ThreeVector(ParticlePos[0]/CLHEP::cm,ParticlePos[1]/CLHEP::cm,ParticlePos[2]/CLHEP::cm);
      tmp_disp  = tmp_fpos - tmp_ipos;
      tmp_dist += tmp_disp.mag();
      tmp_ipos  = tmp_fpos;
      exit_vol  = true;
      enter_vol = false;
      dist_vol += tmp_dist;
    }
  }

  return dist_vol;

}
