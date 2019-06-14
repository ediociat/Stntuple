///////////////////////////////////////////////////////////////////////////////
// Class StntupleMaker : fills Stntuple (P.Murat)
// ------------------------------------------
// order of the data blocks is essential - they are filled according to the
// order in which they are declared...
//
///////////////////////////////////////////////////////////////////////////////

#ifdef __GNUG__
#pragma implementation
#endif

#include <string>
#include <cstdio>

#include "Stntuple/obj/TStnEvent.hh"

#include "fhiclcpp/ParameterSet.h"

#include <assert.h>
#include <iostream>
#include <iomanip>

#include "TNamed.h"
#include "TH1.h"
#include "TString.h"
#include "TProfile.h"
#include "TFolder.h"
#include "TSystem.h"

#include "Stntuple/obj/TStnNode.hh"
#include "Stntuple/obj/TStnErrorLogger.hh"
#include "Stntuple/obj/TStnTrackBlock.hh"
#include "Stntuple/obj/TStnHelixBlock.hh"
#include "Stntuple/obj/TStrawDataBlock.hh"
#include "Stntuple/obj/TCalDataBlock.hh"
#include "Stntuple/obj/TSimpBlock.hh"
#include "Stntuple/obj/TGenpBlock.hh"
#include "Stntuple/obj/TStepPointMCBlock.hh"
#include "Stntuple/obj/TStnHeaderBlock.hh"
//-----------------------------------------------------------------------------
// nothing new: link-wise, TModule depends on TAnaDump
//-----------------------------------------------------------------------------
#include "Stntuple/mod/TAnaDump.hh"
#include "Stntuple/mod/StntupleModule.hh"
#include "Stntuple/mod/StntupleGlobals.hh"

#include "Stntuple/mod/InitCrvPulseBlock.hh"
#include "Stntuple/mod/InitCrvClusterBlock.hh"
#include "Stntuple/mod/InitSimpBlock.hh"
#include "Stntuple/mod/InitGenpBlock.hh"

#include "Stntuple/mod/InitStntupleDataBlocks.hh"
#include "Stntuple/mod/StntupleUtilities.hh"

#include "Stntuple/base/TNamedHandle.hh"
#include "Stntuple/alg/TStntuple.hh"

#include "Mu2eUtilities/inc/SimParticleTimeOffset.hh"
#include "TrkReco/inc/DoubletAmbigResolver.hh"
#include "TrkDiag/inc/KalDiag.hh"
#include "MCDataProducts/inc/GenId.hh"
#include "RecoDataProducts/inc/HelixSeed.hh"

using namespace std; 

// ClassImp(StntupleMaker)

static const char rcsid[] = "$Name:  $";
// stntuple_get_version is autogenerated, see Stntuple/scripts/create_print_header_routine.sh
void stntuple_get_version(char* ver, char* test);
namespace mu2e {
class StntupleMaker : public StntupleModule {
//------------------------------------------------------------------------------
//  data members
//------------------------------------------------------------------------------
protected:
					// process name, default - PROD
  std::string              fProcessName;
//-----------------------------------------------------------------------------
// switches for individual branches
//-----------------------------------------------------------------------------
  int                      fMakeCalData;
  int                      fMakeClusters;
  int                      fMakeGenp;
  int                      fMakePid;
  int                      fMakeSimp;         // 0:dont store, 1:all;
  int                      fMakeStepPointMC;
  int                      fMakeStrawData;
  int                      fMakeTracks;
  int                      fMakeTrackStrawHits;
  int                      fMakeTimeClusters;
  int                      fMakeHelices;
  int                      fMakeTrackSeeds;
  int                      fMakeTrigger;
  int                      fMakeCrvPulses;
  int                      fMakeCrvClusters;
//-----------------------------------------------------------------------------
// module parameters
//-----------------------------------------------------------------------------
  string                   fGenpCollTag;

  string                   fSimpCollTag;
  string                   fStrawHitsCollTag;
  string                   fStrawDigiMCCollTag;

  string                   fCrvRecoPulseCollTag;            //
  string                   fCrvCoincidenceCollTag;          //
  string                   fCrvCoincidenceClusterCollTag;   //

  string                   fVDHitsCollTag;                  // hits on virtual detectors (StepPointMCCollection)

  vector<string>           fTimeClusterBlockName;
  vector<string>           fTimeClusterCollTag;

  vector<string>           fHelixBlockName;
  vector<string>           fHelixCollTag;

  vector<string>           fTrackSeedBlockName;
  vector<string>           fTrackSeedCollTag;

  vector<string>           fTrackBlockName;
  vector<string>           fTrackCollTag;

  vector<string>           fTrackTsBlockName;  // for each track block, the tag  of the corresponding TrackSeed coll
  vector<string>           fTrackTsCollTag;    // for each track block, the name of the corresponding TrackSeed block

  vector<string>           fTciCollTag;        // collection produced by TrackCaloIntersection module
  vector<string>           fTcmCollTag;        // collection produced by TrackCaloMatching     module
  vector<string>           fTrkQualCollTag;    // collection produced by TrackQuality          module

  vector<string>           fPidBlockName;
  vector<string>           fPidCollTag;

  vector<string>           fTrackSHBlockName;

  vector<string>           fSpmcBlockName;
  vector<string>           fSpmcCollTag;

  string                   fCaloCrystalHitMaker;
  string                   fCaloClusterMaker;
//-----------------------------------------------------------------------------
// initialization of various data blocks
//-----------------------------------------------------------------------------
  StntupleInitCrvPulseBlock*   fInitCrvPulseBlock;
  StntupleInitCrvClusterBlock* fInitCrvClusterBlock;
  StntupleInitSimpBlock*       fInitSimpBlock;
  StntupleInitGenpBlock*       fInitGenpBlock;
//-----------------------------------------------------------------------------
// cut-type parameters
//-----------------------------------------------------------------------------
  GenId                    fGenId       ;  // generated process ID
  
  double                   fMinTActive  ;  // start of the active window
  double                   fMinECrystal ;  // 
  double                   fMinSimpEnergy; // min energy of a particle stored in SIMP block

  string                   fCutHelixSeedCollTag; // helix collection to cut on
  int                      fMinNHelices    ; // min number of helices (for cosmics)

  TNamed*                  fVersion;

  TNamedHandle*            fDarHandle;
  TNamedHandle*            fKalDiagHandle;
  TNamedHandle*            fTimeOffsetMapsHandle;
  //  TNamedHandle*            fMinSimpEnergyHandle;

  DoubletAmbigResolver*    fDar;
  KalDiag*                 fKalDiag;
  SimParticleTimeOffset*   fTimeOffsets;
//------------------------------------------------------------------------------
// function members
//------------------------------------------------------------------------------
public:
					// constructors and destructor

  StntupleMaker(fhicl::ParameterSet const& pset);

  ~StntupleMaker();
//-----------------------------------------------------------------------------
// functions of the module
//-----------------------------------------------------------------------------
  void GetDefTrackCollName(char* Name);

					// ****** setters

//-----------------------------------------------------------------------------
// overloaded virtual functions of EDFilter
//-----------------------------------------------------------------------------
  virtual bool beginRun(art::Run& ARun);
  virtual bool endRun  (art::Run& ARun);
  virtual void beginJob();
  virtual void endJob  ();
  virtual bool filter  (AbsEvent& event);

  //  ClassDef(StntupleMaker,0)
};



//------------------------------------------------------------------------------
// constructors
//------------------------------------------------------------------------------
StntupleMaker::StntupleMaker(fhicl::ParameterSet const& PSet): 
  StntupleModule             (PSet,                   "StntupleMaker"         )
  , fProcessName             (PSet.get<string>        ("processName"         ))

  , fMakeCalData             (PSet.get<int>           ("makeCalData"         ))
  , fMakeClusters            (PSet.get<int>           ("makeClusters"        ))
  , fMakeGenp                (PSet.get<int>           ("makeGenp"            ))
  , fMakePid                 (PSet.get<int>           ("makePid"             ))
  , fMakeSimp                (PSet.get<int>           ("makeSimp"            ))
  , fMakeStepPointMC         (PSet.get<int>           ("makeStepPointMC"     ))
  , fMakeStrawData           (PSet.get<int>           ("makeStrawData"       ))
  , fMakeTracks              (PSet.get<int>           ("makeTracks"          ))
  , fMakeTrackStrawHits      (PSet.get<int>           ("makeTrackStrawHits"  ))
  , fMakeTimeClusters        (PSet.get<int>           ("makeTimeClusters"    ))
  , fMakeHelices             (PSet.get<int>           ("makeHelices"         ))
  , fMakeTrackSeeds          (PSet.get<int>           ("makeTrackSeeds"      ))
  , fMakeTrigger             (PSet.get<int>           ("makeTrigger"         ))
  , fMakeCrvPulses           (PSet.get<int>           ("makeCrvPulses"       ))
  , fMakeCrvClusters         (PSet.get<int>           ("makeCrvClusters"     ))
  
  , fGenpCollTag             (PSet.get<string>        ("genpCollTag"         ))
  , fSimpCollTag             (PSet.get<string>        ("simpCollTag"         ))
  , fStrawHitsCollTag        (PSet.get<string>        ("strawHitsCollTag"    ))
  , fStrawDigiMCCollTag      (PSet.get<string>        ("strawDigiMCCollTag"  ))

  , fCrvRecoPulseCollTag         (PSet.get<string>    ("crvRecoPulseCollTag"         ))
  , fCrvCoincidenceCollTag       (PSet.get<string>    ("crvCoincidenceCollTag"       ))
  , fCrvCoincidenceClusterCollTag(PSet.get<string>    ("crvCoincidenceClusterCollTag"))

  , fVDHitsCollTag           (PSet.get<string>        ("vdHitsCollTag"       ))
  , fTimeClusterBlockName    (PSet.get<vector<string>>("timeClusterBlockName"))
  , fTimeClusterCollTag      (PSet.get<vector<string>>("timeClusterCollTag"  ))
  , fHelixBlockName          (PSet.get<vector<string>>("helixBlockName"      ))
  , fHelixCollTag            (PSet.get<vector<string>>("helixCollTag"        ))

  , fTrackSeedBlockName      (PSet.get<vector<string>>("trackSeedBlockName"  ))
  , fTrackSeedCollTag        (PSet.get<vector<string>>("trackSeedCollTag"    ))

  , fTrackBlockName          (PSet.get<vector<string>>("trackBlockName"      ))
  , fTrackCollTag            (PSet.get<vector<string>>("trackCollTag"        ))
  , fTrackTsBlockName        (PSet.get<vector<string>>("trackTsBlockName"    ))
  , fTrackTsCollTag          (PSet.get<vector<string>>("trackTsCollTag"      ))
  , fTciCollTag              (PSet.get<vector<string>>("tciCollTag"          ))
  , fTcmCollTag              (PSet.get<vector<string>>("tcmCollTag"          ))
  , fTrkQualCollTag          (PSet.get<vector<string>>("trkQualCollTag"      ))
  , fPidBlockName            (PSet.get<vector<string>>("pidBlockName"        ))
  , fPidCollTag              (PSet.get<vector<string>>("pidCollTag"          ))
  , fTrackSHBlockName        (PSet.get<vector<string>>("trackSHBlockName"    ))
  
  , fSpmcBlockName           (PSet.get<vector<string>>("spmcBlockName"       ))
  , fSpmcCollTag             (PSet.get<vector<string>>("spmcCollTag"         ))

  , fCaloCrystalHitMaker     (PSet.get<string>        ("caloCrystalHitsMaker"))
  , fCaloClusterMaker        (PSet.get<string>        ("caloClusterMaker"    ))

  , fGenId(GenId::findByName (PSet.get<std::string>   ("genId"              )))

  , fMinTActive              (PSet.get<double>        ("minTActive"          ))
  , fMinECrystal             (PSet.get<double>        ("minECrystal"         ))
  , fMinSimpEnergy           (PSet.get<double>        ("minSimpEnergy"       ))
  , fCutHelixSeedCollTag     (PSet.get<string>        ("cutHelixSeedCollTag" ))
  , fMinNHelices             (PSet.get<int>           ("minNHelices"         ))
{

  char  ver[20], text[200];
  stntuple_get_version(ver,text);

  fVersion      = new TNamed(ver,text);
  TModule::fFolder->Add(fVersion);

  fInitCrvPulseBlock   = nullptr;
  fInitCrvClusterBlock = nullptr;
  fInitSimpBlock       = nullptr;
  fInitGenpBlock       = nullptr;
//-----------------------------------------------------------------------------
// fTimeOffsets is owned by the TAnaDump singleton
//-----------------------------------------------------------------------------
  fTimeOffsets          = TModule::fDump->TimeOffsets();

  fTimeOffsetMapsHandle = new TNamedHandle("TimeOffsetMapsHandle",fTimeOffsets);
  fDar                  = new DoubletAmbigResolver (PSet.get<fhicl::ParameterSet>("DoubletAmbigResolver"),0.,0,0);
  fDarHandle            = new TNamedHandle("DarHandle",fDar);
  fKalDiag              = new KalDiag     (PSet.get<fhicl::ParameterSet>("KalDiag",fhicl::ParameterSet()));
  fKalDiagHandle        = new TNamedHandle("KalDiagHandle"      ,fKalDiag);

  fFolder->Add(fDarHandle);
  fFolder->Add(fKalDiagHandle);
  fFolder->Add(fTimeOffsetMapsHandle);
}


//------------------------------------------------------------------------------
StntupleMaker::~StntupleMaker() {
  delete fDar;
  delete fDarHandle;
  delete fVersion;

  if (fInitCrvPulseBlock  ) delete fInitCrvPulseBlock;
  if (fInitCrvClusterBlock) delete fInitCrvClusterBlock;
  if (fInitSimpBlock      ) delete fInitSimpBlock;
}


//------------------------------------------------------------------------------
bool StntupleMaker::beginRun(art::Run& aRun) {

  static int first_begin_run = 1;

  THistModule::beforeBeginRun(aRun);

  if (first_begin_run) {
//-----------------------------------------------------------------------------
// if we runnning stnmaker_prod.exe, save revision of the TCL file in STNTUPLE
//-----------------------------------------------------------------------------
    first_begin_run = 0;
    const char* c = gSystem->Getenv("STNMAKER_PROD_TCL");
    if (c) TModule::fFolder->Add(new TNamed("STNMAKER_PROD_TCL",c));
    else   TModule::fFolder->Add(new TNamed("STNMAKER_PROD_TCL","unknown"));
  }

  THistModule::afterBeginRun(aRun);

  return 1;
}

//------------------------------------------------------------------------------
bool StntupleMaker::endRun(art::Run& aRun ) {
  THistModule::beforeEndRun(aRun);
  THistModule::afterEndRun (aRun);
  return 1;
}


//------------------------------------------------------------------------------
void StntupleMaker::endJob() {

  THistModule::beforeEndJob();
  THistModule::afterEndJob ();

}

//------------------------------------------------------------------------------
void StntupleMaker::beginJob() {

  int split_mode, compression_level, buffer_size;

  THistModule::beforeBeginJob();

  // create data blocks and branches

  fgStntupleFolder->Add(new TNamed("ProcessName"     ,fProcessName));

					// for the moment do it by hands...
					// create default branches to go into 
					// STNTUPLE

  split_mode        = THistModule::SplitLevel();
  compression_level = THistModule::CompressionLevel();
  buffer_size       = THistModule::BufferSize();
//-----------------------------------------------------------------------------
// header block is always there
//-----------------------------------------------------------------------------
  TStnHeaderBlock* header = (TStnHeaderBlock*) Event()->GetDataBlock("HeaderBlock");
  header->AddCollName("mu2e::StrawHitCollection",fStrawHitsCollTag.data());
//-----------------------------------------------------------------------------
// calorimeter hit data
// this is not RAW hit data yet...
//-----------------------------------------------------------------------------
  if (fMakeCalData) {
    TStnDataBlock* cal_data;

    cal_data = AddDataBlock("CalDataBlock","TCalDataBlock",
			    StntupleInitMu2eCalDataBlock,
			    buffer_size,
			    split_mode,
			    compression_level);
    if (cal_data) {
      cal_data->AddCollName("mu2e::CaloCrystalHitCollection",fCaloCrystalHitMaker.data());
    }
  }
//-----------------------------------------------------------------------------
// CRV
//-----------------------------------------------------------------------------
  if (fMakeCrvPulses) {

    fInitCrvPulseBlock = new StntupleInitCrvPulseBlock();
    fInitCrvPulseBlock->SetCrvRecoPulseCollTag(fCrvRecoPulseCollTag);
    fInitCrvPulseBlock->SetCrvCoincidenceCollTag(fCrvCoincidenceCollTag);
    fInitCrvPulseBlock->SetCrvCoincidenceClusterCollTag(fCrvCoincidenceClusterCollTag);

    AddDataBlock("CrvPulseBlock","TCrvPulseBlock",
		 fInitCrvPulseBlock,
		 buffer_size,
		 split_mode,
		 compression_level);
  }

  if (fMakeCrvClusters) {

    fInitCrvClusterBlock = new StntupleInitCrvClusterBlock();
    fInitCrvClusterBlock->SetCrvRecoPulseCollTag(fCrvRecoPulseCollTag);
    fInitCrvClusterBlock->SetCrvCoincidenceClusterCollTag(fCrvCoincidenceClusterCollTag);

    AddDataBlock("CrvClusterBlock","TCrvClusterBlock",
		 fInitCrvClusterBlock,
		 buffer_size,
		 split_mode,
		 compression_level);
  }
//-----------------------------------------------------------------------------
// straw hit data
//-----------------------------------------------------------------------------
  if (fMakeStrawData) {
    TStnDataBlock* db = AddDataBlock("StrawDataBlock","TStrawDataBlock",
				     StntupleInitMu2eStrawDataBlock,
				     buffer_size,
				     split_mode,
				     compression_level);
    if (db) {
      db->AddCollName("mu2e::StrawHitCollection"   ,fStrawHitsCollTag.data());
      db->AddCollName("mu2e::StrawDigiMCCollection",fStrawDigiMCCollTag.data());
    }
  }
//--------------------------------------------------------------------------------
// time clusters
//--------------------------------------------------------------------------------
  if (fMakeTimeClusters) {
    int nblocks = fTimeClusterBlockName.size();

    for (int i=0; i<nblocks; i++) {
      const char* block_name = fTimeClusterBlockName[i].data();
      TStnDataBlock* db   = AddDataBlock(block_name, 
					 "TStnTimeClusterBlock",
					 StntupleInitMu2eTimeClusterBlock,
					 buffer_size,
					 split_mode,
					 compression_level);
      if (db) {
	db->AddCollName("mu2e::HelixSeedCollection"  ,fHelixCollTag[i].data()      );
 	db->AddCollName("mu2e::TimeClusterCollection",fTimeClusterCollTag[i].data());
	//      SetResolveLinksMethod(block_name,StntupleInitMu2eTimeClusterBlockLinks);
      }
    }
  }
//--------------------------------------------------------------------------------
// helix data
//--------------------------------------------------------------------------------
  if (fMakeHelices) {
    int nb = fHelixBlockName.size();

    for (int i=0; i<nb; i++) {
      const char*    block_name = fHelixBlockName[i].data();
      TStnDataBlock* db         = AddDataBlock(block_name, 
					       "TStnHelixBlock",
					       StntupleInitMu2eHelixBlock,
					       buffer_size,
					       split_mode,
					       compression_level);
      if (db) {
	db->AddCollName("mu2e::HelixSeedCollection"  ,fHelixCollTag[i].data()      );
 	db->AddCollName("mu2e::TimeClusterCollection",fTimeClusterCollTag[i].data());
	db->AddCollName("mu2e::StrawDigiMCCollection",fStrawDigiMCCollTag.data()   );
	//      SetResolveLinksMethod(block_name,StntupleInitMu2eHelixBlockLinks);
     }
    }
  }
//--------------------------------------------------------------------------------
// trackSeed data
//--------------------------------------------------------------------------------
  if (fMakeTrackSeeds) {
    int nb = fTrackSeedBlockName.size();

    for (int i=0; i<nb; i++) {
      const char* block_name = fTrackSeedBlockName[i].data();
      TStnDataBlock* db      = AddDataBlock(block_name, 
					    "TStnTrackSeedBlock",
					    StntupleInitMu2eTrackSeedBlock,
					    buffer_size,
					    split_mode,
					    compression_level);
      if (db) {
	db->AddCollName("mu2e::HelixSeedCollection"  ,fHelixCollTag[i].data()    );
	db->AddCollName("mu2e::KalSeedCollection"    ,fTrackSeedCollTag[i].data());
	db->AddCollName("HelixBlockName"             ,fHelixBlockName[i].data()  );
	db->AddCollName("mu2e::StrawDigiMCCollection",fStrawDigiMCCollTag.data() );
	db->AddCollName("mu2e::ComboHitCollection"   ,fStrawHitsCollTag.data()   );
	SetResolveLinksMethod(block_name,StntupleInitMu2eTrackSeedBlockLinks);	
      }
    }
  }
//-----------------------------------------------------------------------------
// track straw hits
//-----------------------------------------------------------------------------
  if (fMakeTrackStrawHits) {
    int nblocks = fTrackSHBlockName.size();

    for (int i=0; i<nblocks; i++) {
      const char* block_name = fTrackSHBlockName[i].data();
      TStnDataBlock  *block = AddDataBlock(block_name,
					   "TTrackStrawHitBlock",
					   StntupleInitMu2eTrackStrawHitBlock,
					   buffer_size,
					   split_mode,
					   compression_level);
      if (block) {
	block->AddCollName("mu2e::KalRepCollection"              ,fTrackCollTag[i].data ());
	block->AddCollName("mu2e::StrawHitCollection"            ,fStrawHitsCollTag.data());
	//      SetResolveLinksMethod(block_name,StntupleInitMu2eTrackBlockLinks);
      }
    }
  }
//-----------------------------------------------------------------------------
// track branches: for ROOT v3 to use streamers one has to specify split=-1
//-----------------------------------------------------------------------------
  if (fMakeTracks) {
    int nblocks = fTrackBlockName.size();

    for (int i=0; i<nblocks; i++) {
      const char* block_name = fTrackBlockName[i].data();
      TStnDataBlock *track_data = AddDataBlock(block_name,
					       "TStnTrackBlock",
					       StntupleInitMu2eTrackBlock,
					       buffer_size,
					       split_mode,
					       compression_level);
//-----------------------------------------------------------------------------
// each track points back to its seed
// if nshortblocks != 0, for each track we store an index to that tracks's seed
//-----------------------------------------------------------------------------
      if (track_data) {
	track_data->AddCollName("mu2e::KalRepCollection"              ,fTrackCollTag[i].data()    );
	track_data->AddCollName("mu2e::ComboHitCollection"            ,fStrawHitsCollTag.data()   );
	track_data->AddCollName("mu2e::StrawDigiMCCollection"         ,fStrawDigiMCCollTag.data() );
	track_data->AddCollName("mu2e::TrkCaloIntersectCollection"    ,fTciCollTag [i].data()     );
	track_data->AddCollName("mu2e::CaloClusterCollection"         ,fCaloClusterMaker.data()   );
	track_data->AddCollName("mu2e::TrackClusterMatchCollection"   ,fTcmCollTag[i].data()      );
	track_data->AddCollName("mu2e::TrkQualCollection"             ,fTrkQualCollTag[i].data()  );
	track_data->AddCollName("mu2e::PIDProductCollection"          ,fPidCollTag[i].data()      );
	track_data->AddCollName("mu2e::StepPointMCCollection"         ,fVDHitsCollTag.data()      );
	track_data->AddCollName("DarHandle"                           ,GetName()                  ,"DarHandle"    );
	track_data->AddCollName("KalDiagHandle"                       ,GetName()                  ,"KalDiagHandle");
	if (fTrackTsBlockName.size() > 0) {
	  track_data->AddCollName("TrackTsBlockName"                    ,fTrackTsBlockName[i].data());
	  track_data->AddCollName("TrackTsCollTag"                      ,fTrackTsCollTag  [i].data());
	}

	SetResolveLinksMethod(block_name,StntupleInitMu2eTrackBlockLinks);
      }
    }
  }
//-----------------------------------------------------------------------------
// calorimeter clusters 
//-----------------------------------------------------------------------------
  if (fMakeClusters) {
    TStnDataBlock* db = AddDataBlock("ClusterBlock",
				     "TStnClusterBlock",
				     StntupleInitMu2eClusterBlock,
				     buffer_size,
				     split_mode,
				     compression_level);
    if (db) {
      db->AddCollName("mu2e::CaloClusterCollection",fCaloClusterMaker.data());
      SetResolveLinksMethod("ClusterBlock",StntupleInitMu2eClusterBlockLinks);
    }
  }
//-----------------------------------------------------------------------------
// PID - one PID block per track block
//-----------------------------------------------------------------------------
  if (fMakePid) {
    int nblocks = fTrackBlockName.size();

    for (int i=0; i<nblocks; i++) {
      const char* block_name  = fPidBlockName[i].data();
      TStnDataBlock* db = AddDataBlock(block_name,
				       "TStnPidBlock",
				       StntupleInitMu2ePidBlock,
				       buffer_size,
				       split_mode,
				       compression_level);
      if (db) {
	db->AddCollName("mu2e::AvikPIDNewProductCollection",fPidCollTag[i].data());
	SetResolveLinksMethod(block_name,StntupleInitMu2ePidBlockLinks);
      }
    }
  }
//-----------------------------------------------------------------------------
// generator particles 
//-----------------------------------------------------------------------------
  if (fMakeGenp) {
    fInitGenpBlock = new StntupleInitGenpBlock();
    fInitGenpBlock->SetGenpCollTag(fGenpCollTag);
    fInitGenpBlock->SetGenProcessID (fGenId.id());

    AddDataBlock("GenpBlock","TGenpBlock",fInitGenpBlock,buffer_size,split_mode,compression_level);
  }
//-----------------------------------------------------------------------------
// simulated particles 
//-----------------------------------------------------------------------------
  if (fMakeSimp) {
    fInitSimpBlock = new StntupleInitSimpBlock();
    fInitSimpBlock->SetSimpCollTag(fSimpCollTag);
    fInitSimpBlock->SetStrawHitsCollTag(fStrawHitsCollTag);
    fInitSimpBlock->SetStrawDigiMCCollTag(fStrawDigiMCCollTag);
    fInitSimpBlock->SetVDHitsCollTag(fVDHitsCollTag);
    fInitSimpBlock->SetMinSimpEnergy(fMinSimpEnergy);
    fInitSimpBlock->SetGenProcessID (fGenId.id());

    AddDataBlock("SimpBlock","TSimpBlock",fInitSimpBlock,buffer_size,split_mode,compression_level);
  }
//-----------------------------------------------------------------------------
// StepPointMC collections - could be several
//-----------------------------------------------------------------------------
  if (fMakeStepPointMC) {
    int nblocks = fSpmcBlockName.size();

    for (int i=0; i<nblocks; i++) {
      const char* block_name = fSpmcBlockName[i].data();
      TStnDataBlock* db = AddDataBlock(block_name,
				       "TStepPointMCBlock",
				       StntupleInitMu2eStepPointMCBlock,
				       buffer_size,
				       split_mode,
				       compression_level);
      if (db) {
	((TStepPointMCBlock*) db)->SetGenProcessID(fGenId.id());

	db->AddCollName("mu2e::StepPointMCCollection",fSpmcCollTag[i].data());
	db->AddCollName("TimeOffsetMapsHandle"       ,GetName()            ,"TimeOffsetMapsHandle");
      }
    }
    //      SetResolveLinksMethod(block_name,StntupleInitMu2eTrackBlockLinks);
  }

  THistModule::afterEndJob();
}

//_____________________________________________________________________________
bool StntupleMaker::filter(AbsEvent& AnEvent) {

  // when execution comes here al the registered data blocks are already
  // initialized with the event data. Left: variables in the data blocks
  // which depend on the variable defined in other blocks, like track number
  // for a muon or an electron - the idea is that these are defined during the
  // 2nd loop in FillStntupleModule, where ResolveLinks methods are called
  // for each data block

//-----------------------------------------------------------------------------
// connect to the error reporting facility
//-----------------------------------------------------------------------------
//  TStnErrorLogger* logger = Event()->GetErrorLogger();
//   logger->Connect("Report(Int_t, const char*)",
// 		  "StntupleModule",
// 		  this,
// 		  "LogError(const char*)");
//-----------------------------------------------------------------------------
// disconnect from the error reporting signal and return back to AC++
//-----------------------------------------------------------------------------
//   logger->Disconnect("Report(Int_t,const char*)",
// 		     this,"LogError(Int_t,const char*)");

  bool passed(1);

  if (fMinNHelices > 0) {
    auto hH = AnEvent.getValidHandle<mu2e::HelixSeedCollection>(fCutHelixSeedCollTag);
    int  nh = hH->size();
    passed  = (nh >= fMinNHelices);
  }

  return passed;
}


//------------------------------------------------------------------------------
void StntupleMaker::GetDefTrackCollName(char* Name) {
  // put in a working kludge first

  strcpy(Name,"default");
}


// //_____________________________________________________________________________
// int StntupleMaker::InitCalDataBlock(TStnDataBlock* Block) {
//   int mode = 0;
//   AbsEvent* event = AbsEnv::instance()->theEvent();
//   return StntupleInitMu2eCalDataBlock(Block,event,mode);
// }

// //_____________________________________________________________________________
// int StntupleMaker::InitHeaderBlock(TStnDataBlock* Block) {
//   int mode = 0;
//   AbsEvent* event = AbsEnv::instance()->theEvent();
//   return StntupleInitMu2eHeaderBlock(Block,event,mode);
// }

// //_____________________________________________________________________________
// int StntupleMaker::InitTrackBlock(TStnDataBlock* Block) {
//   int mode = 0;
//   AbsEvent* event = AbsEnv::instance()->theEvent();
//   return StntupleInitMu2eTrackBlock(Block,event,mode);
// }

//_____________________________________________________________________________
// int StntupleMaker::InitTriggerBlock(TStnDataBlock* Block) {
//   int mode = 0;
//   AbsEvent* event = AbsEnv::instance()->theEvent();
//   return StntupleInitMu2eTriggerBlock(Block,event,mode);
// }

} // end namespace mu2e

using mu2e::StntupleMaker;

DEFINE_ART_MODULE(StntupleMaker);
