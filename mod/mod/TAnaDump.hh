//
#ifndef __murat_inc_TAnaDump_hh__
#define __murat_inc_TAnaDump_hh__

#include "TObject.h"
#include "TObjArray.h"
#include "TString.h"
#include "TGraph.h"
#include "TMarker.h"
#include "TCanvas.h"
#include "TEllipse.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"

#ifndef __CINT__

#include "art/Framework/Principal/Event.h"

#include "MCDataProducts/inc/CaloClusterMCTruthAssn.hh"
#include "MCDataProducts/inc/CaloHitMCTruthAssn.hh"
#include "MCDataProducts/inc/StrawDigiMCCollection.hh"

#include "RecoDataProducts/inc/ComboHit.hh"
#include "RecoDataProducts/inc/HelixHit.hh"
#include "RecoDataProducts/inc/StrawHitPositionCollection.hh"
#include "DataProducts/inc/XYZVec.hh"

#include "TrkReco/inc/TrkPrintUtils.hh"

#else

namespace art {
  class Event;
}
#endif

namespace mu2e {

#ifndef MCDataProducts_StrawDigiMC_hh
  class StrawDigiMCCollection;
#endif

  class StrawHit;
  class StrawHitMCTruth;
  class CaloCluster;
  class CaloProtoCluster;
  class TrkToCaloExtrapol;
  class StepPointMC;
  class GenParticle;
  class SimParticle;
  class TimeCluster;
  class KalSeed;
  class ComboHit;
  class HelixSeed;
  class TrackClusterMatch;
  class TrkCaloHit;
  class TrkStrawHit;
  class SimParticleTimeOffset;
  class TrkPrintUtils;
}

class KalRep;

class TAnaDump : public TObject {
public:

  const art::Event*               fEvent;
  TObjArray*                      fListOfObjects;
  TString                         fFlagBgrHitsModuleLabel;
  TString                         fStrawDigiMCCollTag;
  mu2e::SimParticleTimeOffset*    fTimeOffsets;
  const mu2e::StrawDigiMCCollection*    _mcdigis;
  double                          fTmp[100];  // for testing

  mu2e::TrkPrintUtils*            _printUtils;

private:

  TAnaDump(const fhicl::ParameterSet* Pset = NULL);
  ~TAnaDump();
  
  class  Cleaner {
  public: 
    Cleaner();
    ~Cleaner();
  };
  
  friend class Cleaner;

  static TAnaDump*  fgInstance;
public:
//-----------------------------------------------------------------------------
// TAnaDump gets initialized by the first TModule instantiated
//-----------------------------------------------------------------------------
  static TAnaDump*  Instance(const fhicl::ParameterSet* PSet = NULL);
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  mu2e::SimParticleTimeOffset*   TimeOffsets() { return fTimeOffsets; }
//-----------------------------------------------------------------------------
// other methods
//-----------------------------------------------------------------------------
  void   AddObject      (const char* Name, void* Object);
  void*  FindNamedObject(const char* Name);

  void SetEvent(art::Event& Evt) { fEvent = &Evt; }

  void SetFlagBgrHitsModuleLabel(const char* Tag) { fFlagBgrHitsModuleLabel = Tag; }
  void SetStrawDigiMCCollTag    (const char* Tag) { fStrawDigiMCCollTag     = Tag; }

  double evalWeight(const mu2e::ComboHit* Hit   ,
		    XYZVec& StrawDir ,
		    XYZVec& HelCenter, 
		    double             Radius   ,
		    int                WeightMode,
		    fhicl::ParameterSet const& Pset);

  void printEventHeader();

  void printCaloCrystalHits (const char* ModuleLabel, 
			     const char* ProductName = "", 
			     const char* ProcessName = ""); 

  void printCaloDigiCollection(const char* ModuleLabel, 
			       const char* ProductName = "", 
			       const char* ProcessName = ""); 

  void printCaloRecoDigiCollection(const char* ModuleLabel, 
				   const char* ProductName = "", 
				   const char* ProcessName = "");

  void printCaloHits        (const char* ModuleLabel, 
			     const char* ProductName = "", 
			     const char* ProcessName = ""); 

  void printDiskCalorimeter();

  void printCaloCluster          (const mu2e::CaloCluster* Cluster ,
				  const char*              Opt = "",
				  const mu2e::CaloHitMCTruthAssns* CaloHitTruth=NULL);
  
  void printCaloClusterCollection (const char* ModuleLabel, 
				   const char* ProductName= "",
				   const char* ProcessName= "",
				   double      Emin=50.,
				   int         hitOpt=0,
				   const char* MCModuleLabel= "");
  
  void printCaloProtoCluster      (const mu2e::CaloProtoCluster* Clu     ,
				   const char*                   Opt = "");
  
  void printCaloProtoClusterCollection (const char* ModuleLabel, 
					const char* ProductName,
					const char* ProcessName);
  
  void printTrackSeed          (const mu2e::KalSeed* TrkSeed, 
				const char* Opt = ""    ,
				const char* ModuleLabelStrawHit="makeSH");

  void printTrackSeedCollection(const char* ModuleLabel             , 
				const char* ProductName         = "", 
				const char* ProcessName         = "",
				int         hitOpt              = 0,
				const char* ModuleLabelStrawHit ="makeSH"); 

  void printHelixSeed          (const mu2e::HelixSeed*         Helix   , 
				//				const char* HelixSeedCollTag           ,
				const char* StrawHitCollTag            ,  // usually - "makeSH"
				const char* StrawDigiCollTag = "makeSD",
				const char* Opt              = ""      );

  void printHelixSeedCollection(const char* HelixSeedCollTag           ,  // always needed
				int         PrintHits        = 0       ,
				const char* StrawHitCollTag  = "makeSH",  // usually, "makeSH"
				const char* StrawDigiCollTag = "makeSD");

  void printTrkCaloHit(const KalRep* Krep, mu2e::TrkCaloHit* CaloHit);

  void printKalRep(const KalRep* Krep, const char* Opt = "", const char* Prefix = "");

  void printKalRepCollection(const char* ModuleLabel     , 
			     const char* ProductName = "", 
			     const char* ProcessName = "",
			     int         hitOpt      = 0); 
  
  void printTrkToCaloExtrapol           (const mu2e::TrkToCaloExtrapol*extrk,
					 const char* Opt = "");

  void printTrkToCaloExtrapolCollection (const char* ModuleLabel, 
					 const char* ProductName = "", 
					 const char* ProcessName = "");

  void printTimeCluster   (const mu2e::TimeCluster* TimePeak, const char* Opt = "", 
			   const mu2e::ComboHitCollection* ChColl=0,
			   const char*StrawDigiMCModuleLabel = "makeSD");

  void printTimeClusterCollection(const char* ModuleLabel        , 
				  const char* ComboHitModuleLabel,
				  const char* ProductName = ""   , 
				  const char* ProcessName = ""   ,
				  int         PrintHits   = 0    ,
				  const char* StrawDigiMCModuleLabel = "makeSD");

//   void printCaloCrystalHit(const CaloCrystalHit* Hit, const char* Opt = "");
//   void printCaloHit       (const CaloHit*        Hit, const char* Opt = "");
  void printGenParticle   (const mu2e::GenParticle*    P  , const char* Opt = "");

  void printGenParticleCollections();

  void printSimParticle   (const mu2e::SimParticle*    P  , const char* Opt = "", const void* PrintData = nullptr);

  void printSimParticleCollection(const char* ModuleLabel     , 
				  const char* ProductName = "", 
				  const char* ProcessName = "");

  void printStrawHit      (const mu2e::StrawHit*    Hit, 
			   const mu2e::StepPointMC* Step,
			   const char*              Opt   = "", 
			   int                      INit  = -1,
			   int                      Flags = -1);
  
  void printComboHit      (const mu2e::ComboHit*    Hit, 
			   const mu2e::StepPointMC* Step,
			   const char*              Opt   = "", 
			   int                      INit  = -1,
			   int                      Flags = -1);
  
  void printComboHitCollection (const char* StrawHitCollTag, 
				const char* FlagBgrHitsCollTag,  // 'FlagBgrHits:StrawHits' or 'FlagBgrHits:ComboHits'
				const char* StrawDigiMCCollTag,  // "makeSD" or "compressDigiMCs"
				const char* ProcessName = "",
				double TMin = -1.e6,
				double TMax =  1.e6);
 
  void printHelixHit      (const mu2e::HelixHit*    HelHit,
			   const mu2e::ComboHit*    Hit, 
			   const mu2e::StepPointMC* Step,
			   const char*              Opt   = "", 
			   int                      INit  = -1,
			   int                      Flags = -1);
  
  void printStrawHitCollection (const char* ModuleLabel, 
				const char* ProductName = "", 
				const char* ProcessName = "",
				double TMin = -1.e6,
				double TMax =  1.e6);

  void printStrawHitMCTruth      (const mu2e::StrawHitMCTruth* Hit, const char* Opt = "");

  void printStrawHitMCTruthCollection (const char* ModuleLabel, 
				       const char* ProductName = "", 
				       const char* ProcessName = "");

  // pass the detector name to know what to print for different detectors
  // tested for Detector = 'tracker', 'calorimeter'
  void printStepPointMC(const mu2e::StepPointMC* Step, const char* Detector, const char* Opt = "");

  void printStepPointMCCollection (const char* ModuleLabel     , 
				   const char* ProductName = "", 
				   const char* ProcessName = "");

  void printStrawHitPosition(const mu2e::StrawHitPosition* Pos, const char* Opt = "");

  void printStrawHitPositionCollection (const char* ModuleLabel     , 
					const char* ProductName = "", 
					const char* ProcessName = "");

  void printTrackClusterMatch(const mu2e::TrackClusterMatch* TcMatch, const char* Option);

  void printTrackClusterMatchCollection (const char* ModuleLabel     , 
					 const char* ProductName = "", 
					 const char* ProcessName = "");
			      
					// refit track dropping hits away > NSig sigma (0.1)
  void  refitTrack(void* Trk, double NSig);

  void  Test_000(const KalRep* Krep, mu2e::TrkStrawHit* Hit);

  ClassDef(TAnaDump,0)
};


#endif
