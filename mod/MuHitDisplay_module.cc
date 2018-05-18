///////////////////////////////////////////////////////////////////////////////
// A half-interactive 2D event display. 
//
// $Id: MuHitDisplay_module.cc,v 1.6 2014/09/20 17:54:06 murat Exp $
// $Author: murat $
// $Date: 2014/09/20 17:54:06 $
//
// Contact person:  Pavel Murat, Gianantonio Pezzulo
//
// What this event display shows: 2D view of the straw hits, tracks, and calorimeter clusters
//
// Straw hit display mode: 
// -----------------------
// displayBackgroundHits : false - respect hit flags set by FlagStrawHits_module (default)
//                         the timing hit flag check is commented out - I didn't have a MC file
//                         in hands to check
//                       : true  - display all hits
// useStereoHits         : true : displayed hit position errors are defined by the StrawHitPositionCollection
//                       : false: semi-random assignment (sigr=5mm, sigp = strawTimeDivisionErr/2.)
//
// small black triangles: MC truth on the trajectory
//
// red   hits: hits on a track reconstructed as a downstream-moving conversion electron within the time window
//             size of the timing window can be redefined via  talk-to (input .fcl file)
// blue  hits: hits of the conversion electron track, which fall outside the time window
// black hits: hits produced by anything, but the conversion electron
//
// a hit is displayed by a cross with the radial error bar of 5mm and the resolution 
// along the straw of sigma(time division)/2, to better guide the eye
//
// green circles - contour of the disk-based calorimeter
// clusters on the 1st disk are shown in red, on the second disk - in pink
// 
// there are few other features which we need to document
//
// .fcl file to use: Analyses/test/hitDisplay.fcl
///////////////////////////////////////////////////////////////////////////////

// C++ includes.
#include <map>
#include <iostream>
#include <string>

// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Principal/Selector.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Mu2e includes.
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/TrackerCalibrations.hh"

#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"

#include "TTrackerGeom/inc/TTracker.hh"
// #include "CalorimeterGeom/inc/VaneCalorimeter.hh"
#include "CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"
//#include "CaloCluster/inc/CaloClusterUtilities.hh"
#include "Mu2eUtilities/inc/SimParticlesWithHits.hh"
#include "Mu2eUtilities/inc/SortedStepPoints.hh"
#include "Mu2eUtilities/inc/TrackTool.hh"

#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/StrawDigiMCCollection.hh"

#include "DataProducts/inc/StrawEnd.hh"

#include "BTrk/TrkBase/HelixParams.hh"
// #include "BTrk/TrkBase/TrkHotList.hh"
#include "BTrk/KalmanTrack/KalHit.hh"
#include "RecoDataProducts/inc/TrkFitDirection.hh"
#include "BTrk/TrkBase/TrkParticle.hh"

#include "RecoDataProducts/inc/CaloCrystalHit.hh"
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "RecoDataProducts/inc/CaloClusterCollection.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/XYZVec.hh"
#include "RecoDataProducts/inc/StrawHitFlagCollection.hh"

#include "RecoDataProducts/inc/CrvRecoPulseCollection.hh"
#include "RecoDataProducts/inc/CrvRecoPulse.hh"
#include "RecoDataProducts/inc/TimeCluster.hh"

#include "BTrkData/inc/TrkStrawHit.hh"
#include "RecoDataProducts/inc/KalRepPtrCollection.hh"


#include "TrkReco/inc/DoubletAmbigResolver.hh"
#include "TrkDiag/inc/KalDiag.hh"

// ROOT includes
#include "TApplication.h"
#include "TArc.h"
#include "TArrow.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TLine.h"
#include "TBox.h"
#include "TMarker.h"
#include "TEllipse.h"
#include "TText.h"
#include "TNtuple.h"
#include "TDatabasePDG.h"

#include "Stntuple/base/TNamedHandle.hh"
#include "Stntuple/gui/gui/TStnVisManager.hh"
#include "Stntuple/gui/THeaderVisNode.hh"
#include "Stntuple/gui/TCalVisNode.hh"
#include "Stntuple/gui/TCrvVisNode.hh"
#include "Stntuple/gui/TTrkVisNode.hh"
#include "Stntuple/gui/TMcTruthVisNode.hh"

#include "Stntuple/alg/TStnTrackID.hh"

#include "Stntuple/mod/TAnaDump.hh"
#include "Stntuple/mod/THistModule.hh"

#include "Stntuple/obj/TStnHeaderBlock.hh"

#include "Stntuple/obj/TStnTrackBlock.hh"
#include "Stntuple/obj/TStnClusterBlock.hh"

#include "CalPatRec/inc/McUtilsToolBase.hh"

using namespace std;
using CLHEP::Hep3Vector;

namespace mu2e {

  class MuHitDisplay : public THistModule {
  private:
//-----------------------------------------------------------------------------
// Input parameters: Module labels 
//-----------------------------------------------------------------------------
    std::string        _moduleLabel;	             // this module label
    std::string        _processName;
    std::string        _generatorModuleLabel;
    std::string        _g4ModuleLabel;
    std::string        _caloClusterModuleLabel;
    std::string	       _crvRecoPulsesModuleLabel;

    std::string        _producerName;

    std::string        _makeStrawHitModuleLabel;
    std::string        _makeComboHitModuleLabel;
    std::string        _makeStrawDigiModuleLabel;
    // std::string        fStrawHitPosMaker;
    std::string        fStrawHitFlagMaker;

    std::string        fTrkRecoModuleLabel;
    std::string        fCrystalHitMaker;
    std::string        fTrkExtrapol;
    std::string        fTrkCalMatch;
    std::string        fPidModuleLabel;

    TrkFitDirection    fTrkDirection;
    TrkParticle        fParticleHypo;


    int                fGeneratorID;
    string             _trackerStepPoints;
    string             fDirectionAndParticle;

    double             _minEnergyDep;
    double             _timeWindow;

    mu2e::StrawHitFlag fGoodHitMask;
    mu2e::StrawHitFlag fBadHitMask;
    size_t             _minHits;


    bool               fDisplayBackgroundHits;
    bool               fPrintHits;
    int                fUseStereoHits;
//-----------------------------------------------------------------------------
// Control for CRV-specific viewing
//-----------------------------------------------------------------------------
    bool				_showCRVOnly;
    bool				_showTracks;
    bool				foundCRV;
    bool				foundTrkr;
    bool				foundTrkr_StrawHitColl;
    bool				foundTrkr_StrawDigiMCColl;
    bool				foundTrkr_StrawHitPosColl;
    bool				foundTrkr_StrawHitFlagColl;
    bool				foundCalo_CrystalHitColl;
    bool				foundCalo_ClusterColl;
    bool				foundCalo;

    fhicl::ParameterSet                 _vmConfig;
//-----------------------------------------------------------------------------
// end of input parameters
// Options to control the display
// hit flag bits which should be ON and OFF
//-----------------------------------------------------------------------------
    TApplication*                               fApplication;
    TCanvas*                                    fCanvas;

    const mu2e::Calorimeter*                    fCal;              //

    const mu2e::GenParticleCollection*          _genParticleColl;         // 

    const mu2e::StrawHitCollection*             fStrawHitColl;     // 
    const mu2e::ComboHitCollection*             fShComboHitColl;     // 
    const mu2e::ComboHitCollection*             fComboHitColl;     // 

    // const mu2e::StrawHitPositionCollection*     fStrawHitPosColl;  //
    const mu2e::StrawHitFlagCollection*         fStrawHitFlagColl; //
    const mu2e::StrawDigiMCCollection*         _strawDigiMCColl; //

    const mu2e::CaloCrystalHitCollection*       fListOfCrystalHits;//
    const mu2e::CaloClusterCollection*          fListOfClusters;   //

    const mu2e::PtrStepPointMCVectorCollection* _hits_mcptr;        //

    const mu2e::StepPointMCCollection*          _stepPointMCColl;  //

    const mu2e::TimeClusterCollection*          fTimeClusterColl;  //

    mu2e::CrvRecoPulseCollection*		fCrvPulseColl_Right;
    mu2e::CrvRecoPulseCollection*		fCrvPulseColl_Left;
    mu2e::CrvRecoPulseCollection*		fCrvPulseColl_TopDS;
    mu2e::CrvRecoPulseCollection*		fCrvPulseColl_TopTS;
    mu2e::CrvRecoPulseCollection*		fCrvPulseColl_Dwnstrm;
    mu2e::CrvRecoPulseCollection*		fCrvPulseColl_Upstrm;
		
    const mu2e::KalRepPtrCollection*            _kalRepPtrColl;

    mu2e::SimParticlesWithHits*                 _simParticlesWithHits; // 
    
    std::unique_ptr<McUtilsToolBase>            _mcUtils;


    int                    fNClusters;
    // 4 hypotheses: dem, uep, dmm, ump
    int                    fNTracks[4];
		
    TMarker*               fMarker;

    TStnVisManager*        fVisManager;
    TStnTrackBlock*        fTrackBlock;
    TStnClusterBlock*      fClusterBlock;
    TStnHeaderBlock*       fHeaderBlock;

    TStnTrackID*           fTrackID;

    double                 fTrackerP;   	// CE momentum on entrance to the straw tracker
    double                 fTrackerPt;	// CE Pt on entrance to the straw tracker

    const TimeCluster*     fTimeCluster;

    const TTracker*        fTracker;    // straw tracker geometry

    TNamedHandle*          fDarHandle;
    TNamedHandle*          fKalDiagHandle;

    DoubletAmbigResolver*  fDar;
    KalDiag*               fKalDiag;

  public:
    explicit MuHitDisplay(fhicl::ParameterSet const& pset);
    virtual ~MuHitDisplay();

    int      getData(const art::Event* Evt);
    void     Init(art::Event* Evt);
    void     InitVisManager();
    int		 getCRVSection(int shieldNumber);

    void     printCaloCluster(const CaloCluster* Cl, const char* Opt);
    //-----------------------------------------------------------------------------
    // overloaded virtual methods of the base class
    //-----------------------------------------------------------------------------
    virtual void     beginJob();
    virtual bool     beginRun(art::Run& aRun);
    virtual bool     filter(art::Event& Evt);
  };


  //-----------------------------------------------------------------------------
  MuHitDisplay::MuHitDisplay(fhicl::ParameterSet const& pset) :
    THistModule(pset, "MuHitDisplay"),
    _moduleLabel(pset.get<std::string>("module_label")),
    _processName(pset.get<string>("processName", "")),
    _generatorModuleLabel(pset.get<std::string>("generatorModuleLabel")),
    _g4ModuleLabel(pset.get<std::string>("g4ModuleLabel")),
    _caloClusterModuleLabel(pset.get<std::string>("caloClusterModuleLabel")),
    _crvRecoPulsesModuleLabel(pset.get<std::string>("crvRecoPulsesModuleLabel")),
				
    _producerName(""),

    _makeStrawHitModuleLabel(pset.get<std::string>("strawHitMakerModuleLabel")),
    _makeComboHitModuleLabel(pset.get<std::string>("comboHitMakerModuleLabel")),
    _makeStrawDigiModuleLabel(pset.get<std::string>("strawDigiMakerModuleLabel")),
    // fStrawHitPosMaker(pset.get<std::string>("strawHitPosMakerModuleLabel")),
    fStrawHitFlagMaker(pset.get<std::string>("strawHitFlagMakerModuleLabel")),

    fTrkRecoModuleLabel(pset.get<std::string>("trkRecoModuleLabel")),
    fCrystalHitMaker(pset.get<std::string>("caloCrystalHitsMaker")),
    fTrkExtrapol(pset.get<std::string>("trkExtrapol")),
    fTrkCalMatch(pset.get<std::string>("trkCalMatch")),
    fPidModuleLabel(pset.get<std::string>("pidModuleLabel")),


    fTrkDirection((TrkFitDirection::FitDirection)(pset.get<int>("fitDirection"))),
    fParticleHypo((TrkParticle::type)            (pset.get<int>("fitParticle"))),

    fGeneratorID(pset.get<int>("generatorID", GenId::conversionGun)),
    _trackerStepPoints(pset.get<std::string>("trackerStepPoints")),
    _minEnergyDep(pset.get<double>("minEnergyDep", 0)),
    _timeWindow(pset.get<double>("timeWindow", 1.e6)),
    fGoodHitMask(pset.get<std::vector<std::string> >("goodHitMask")),
    fBadHitMask(pset.get<std::vector<std::string> >("badHitMask")),
    _minHits(pset.get<unsigned>("minHits")),

    fDisplayBackgroundHits(pset.get<bool>("displayBackgroundHits", false)),
    fPrintHits(pset.get<bool>("printHits", false)),
    fUseStereoHits(pset.get<bool>("useStereoHits", false)),
    _showCRVOnly(pset.get<bool>("showCRVOnly", false)),
    _showTracks (pset.get<bool>("showTracks")),


    _vmConfig(pset.get<fhicl::ParameterSet>("visManager", fhicl::ParameterSet()))
  {

    fApplication = 0;
    fCanvas = 0;
    fCal = 0;

    fTrackBlock = new TStnTrackBlock();
    fClusterBlock = new TStnClusterBlock();
    fHeaderBlock = new TStnHeaderBlock();

    fVisManager = TStnVisManager::Instance();

    _simParticlesWithHits = NULL;
    
    _mcUtils = std::make_unique<McUtilsToolBase>();
    
    fTrackID = new TStnTrackID();

    fDirectionAndParticle = fTrkDirection.name() + fParticleHypo.name();

    foundCRV  = false;
    foundTrkr = false;
    foundCalo = false;

    fCrvPulseColl_Right   = new CrvRecoPulseCollection();
    fCrvPulseColl_Left    = new CrvRecoPulseCollection();
    fCrvPulseColl_TopDS   = new CrvRecoPulseCollection();
    fCrvPulseColl_TopTS   = new CrvRecoPulseCollection();
    fCrvPulseColl_Dwnstrm = new CrvRecoPulseCollection();
    fCrvPulseColl_Upstrm  = new CrvRecoPulseCollection();

    fDar           = new DoubletAmbigResolver (pset.get<fhicl::ParameterSet>("DoubletAmbigResolver"),0.,0,0);
    fDarHandle     = new TNamedHandle("DarHandle",fDar);

    fKalDiag       = new KalDiag(pset.get<fhicl::ParameterSet>("KalDiag",fhicl::ParameterSet()));
    fKalDiagHandle = new TNamedHandle("KalDiagHandle",fKalDiag);

    fFolder->Add(fDarHandle);
    fFolder->Add(fKalDiagHandle);
  }

  //-----------------------------------------------------------------------------
  MuHitDisplay::~MuHitDisplay() {

    if (fApplication) delete fApplication;

    delete fCrvPulseColl_Right;
    delete fCrvPulseColl_Left;
    delete fCrvPulseColl_TopDS;
    delete fCrvPulseColl_TopTS;
    delete fCrvPulseColl_Dwnstrm;
    delete fCrvPulseColl_Upstrm;

    delete fKalDiag;
    delete fDar;
  }


  //-----------------------------------------------------------------------------
  void MuHitDisplay::beginJob() {

    //     const char oname[] = "MuHitDisplay::beginJob";
    int    tmp_argc(0);
    char** tmp_argv(0);

    if (!gApplication) {
      fApplication = new TApplication("MuHitDisplay_module", &tmp_argc, tmp_argv);
    }

    // Create a canvas with a guaranteed unique name; the module label is unique within a job.
    TString name  = "canvas_"     + _moduleLabel;
    TString title = "Canvas for " + _moduleLabel;

    fCanvas = new TCanvas(name, title, 800, 800);
    fMarker = new TMarker(0, 0, 20);
    fMarker->SetMarkerSize(0.3);
//-----------------------------------------------------------------------------
// make sure needed plugin libraries are loaded
//-----------------------------------------------------------------------------
    // std::string lib = gSystem->Getenv("MU2E_BASE_RELEASE");
    // lib += "/lib/libmu2e_CalPatRec_utils.so";
    // if (! gInterpreter->IsLoaded(lib.data())) gSystem->Load(lib.data());
//-----------------------------------------------------------------------------
// define collection names to be used for initialization
//-----------------------------------------------------------------------------
//    const char* charDirectionAndParticle = fDirectionAndParticle.c_str();


    fClusterBlock->AddCollName("mu2e::CaloClusterCollection"       ,_caloClusterModuleLabel.data()    ,"");

    fTrackBlock->AddCollName("mu2e::KalRepCollection"              ,fTrkRecoModuleLabel.data()        ,"");
    fTrackBlock->AddCollName("mu2e::StrawHitCollection"            ,_makeStrawHitModuleLabel.data()   ,"");
    fTrackBlock->AddCollName("mu2e::StrawDigiMCCollection"         ,_makeStrawDigiModuleLabel.data()   ,"");
    fTrackBlock->AddCollName("mu2e::PtrStepPointMCVectorCollection",_makeStrawDigiModuleLabel.data()   ,"");
    fTrackBlock->AddCollName("mu2e::TrkCaloIntersectCollection"    ,fTrkExtrapol.data()               ,"");
    fTrackBlock->AddCollName("mu2e::CaloClusterCollection"         ,_caloClusterModuleLabel.data()    ,"");
    fTrackBlock->AddCollName("mu2e::TrackClusterMatchCollection"   ,fTrkCalMatch.data()               ,"");
    fTrackBlock->AddCollName("mu2e::PIDProductCollection"          ,fPidModuleLabel.data()            ,"");
    fTrackBlock->AddCollName("mu2e::StepPointMCCollection"         ,_g4ModuleLabel.data()             ,"");
    fTrackBlock->AddCollName("DarHandle"                           ,GetName()                         ,"DarHandle");
    fTrackBlock->AddCollName("KalDiagHandle"                       ,GetName()                         ,"KalDiagHandle");

    TAnaDump::Instance()->AddObject("MuHitDisplay::TrackBlock"  , fTrackBlock);
    TAnaDump::Instance()->AddObject("MuHitDisplay::ClusterBlock", fClusterBlock);
  }

  //-----------------------------------------------------------------------------
  bool MuHitDisplay::beginRun(art::Run& Run) {
    mu2e::GeomHandle<mu2e::TTracker> handle;
    fTracker = handle.get();
    return true;
  }


  //-----------------------------------------------------------------------------
  // initialize the visualization manager
  // TStnVisManager is responsible for deleting all nodes created here
  //-----------------------------------------------------------------------------
  void MuHitDisplay::InitVisManager() {
    const char oname [] = "MuHitDisplay::InitVisManager";

    fVisManager->SetTitleNode(new THeaderVisNode("HeaderVisNode", fHeaderBlock));

    if (_showCRVOnly) {
      TCrvVisNode      *crv_node[6];
      
      crv_node[0] = new TCrvVisNode("CrvVisNode#0", 0);
      crv_node[0]->SetRecoPulsesCollection(&fCrvPulseColl_Right);
      fVisManager->AddNode(crv_node[0]);
      crv_node[1] = new TCrvVisNode("CrvVisNode#1", 1);
      crv_node[1]->SetRecoPulsesCollection(&fCrvPulseColl_Left);
      fVisManager->AddNode(crv_node[1]);
      crv_node[2] = new TCrvVisNode("CrvVisNode#2", 2);
      crv_node[2]->SetRecoPulsesCollection(&fCrvPulseColl_TopDS);
      fVisManager->AddNode(crv_node[2]);
      crv_node[3] = new TCrvVisNode("CrvVisNode#3", 3);
      crv_node[3]->SetRecoPulsesCollection(&fCrvPulseColl_Dwnstrm);
      fVisManager->AddNode(crv_node[3]);
      crv_node[4] = new TCrvVisNode("CrvVisNode#4", 4);
      crv_node[4]->SetRecoPulsesCollection(&fCrvPulseColl_Upstrm);
      fVisManager->AddNode(crv_node[4]);
      crv_node[5] = new TCrvVisNode("CrvVisNode#5", 8);
      crv_node[5]->SetRecoPulsesCollection(&fCrvPulseColl_TopTS);
      fVisManager->AddNode(crv_node[5]);
    }
    else {
      TCalVisNode      *cal_node[2];
      TTrkVisNode      *trk_node;
      TMcTruthVisNode  *mctr_node;

      const mu2e::DiskCalorimeter  *dc;
//-----------------------------------------------------------------------------
// parse the configuration module parameters
//-----------------------------------------------------------------------------
      int debug_level = _vmConfig.get<int>("debugLevel");
      fVisManager->SetDebugLevel(debug_level);

      int display_straw_digi_mc = _vmConfig.get<int>("displayStrawDigiMC");
      fVisManager->SetDisplayStrawDigiMC(display_straw_digi_mc);
//-----------------------------------------------------------------------------
// do the geometry
//-----------------------------------------------------------------------------
      art::ServiceHandle<mu2e::GeometryService> geom;

      if (geom->hasElement<mu2e::DiskCalorimeter>()) {
	mu2e::GeomHandle<mu2e::DiskCalorimeter> dc_handle;
	dc = dc_handle.get();
//-----------------------------------------------------------------------------
// TCalVisNode: calorimeter, CaloCrystalHits and CaloCrystals
//-----------------------------------------------------------------------------
	cal_node[0] = new TCalVisNode("CalVisNode#0", &dc->disk(0), 0);
	cal_node[0]->SetListOfClusters(&fListOfClusters);
	cal_node[0]->SetListOfCrystalHits(&fListOfCrystalHits);
	cal_node[0]->SetTimeClusterColl(&fTimeClusterColl);
	fVisManager->AddNode(cal_node[0]);

	cal_node[1] = new TCalVisNode("CalVisNode#1", &dc->disk(1), 1);
	cal_node[1]->SetListOfClusters(&fListOfClusters);
	cal_node[1]->SetListOfCrystalHits(&fListOfCrystalHits);
	cal_node[1]->SetTimeClusterColl(&fTimeClusterColl);
	fVisManager->AddNode(cal_node[1]);
      }
      else {
	printf("[%s] >>> ERROR : DiskCalorimeter is not defined\n", oname);
	dc = 0;
      }
//-----------------------------------------------------------------------------
// TrkVisNode: tracker, tracks and straw hits
//-----------------------------------------------------------------------------
      trk_node = new TTrkVisNode("TrkVisNode", fTracker, NULL);
      trk_node->SetStrawHitColl(&fStrawHitColl);
      trk_node->SetStrawHitPosColl(&fShComboHitColl);//fStrawHitPosColl);
      trk_node->SetStrawHitFlagColl(&fStrawHitFlagColl);
      trk_node->SetTimeClusterColl(&fTimeClusterColl);
      trk_node->SetKalRepPtrColl(&_kalRepPtrColl);
      trk_node->SetMcPtrColl(&_hits_mcptr);
      trk_node->SetStrawDigiMCColl(&_strawDigiMCColl);
      fVisManager->AddNode(trk_node);
//-----------------------------------------------------------------------------
// MC truth: StepPointMC's and something else - not entirely sure
//-----------------------------------------------------------------------------
      mctr_node = new TMcTruthVisNode("McTruthVisNode");
      mctr_node->SetListOfHitsMcPtr(&_hits_mcptr);
      mctr_node->SetStepPointMCCollection(&_stepPointMCColl);
      mctr_node->SetSimParticlesWithHits(&_simParticlesWithHits);
      mctr_node->SetGenpColl(&_genParticleColl);
      fVisManager->AddNode(mctr_node);
    }
  }

//-----------------------------------------------------------------------------
// get data from the event record
//-----------------------------------------------------------------------------
  int MuHitDisplay::getData(const art::Event* Evt) {
    const char* oname = "MuHitDisplay::getData";
//-----------------------------------------------------------------------------
//  CRV pulse information
//-----------------------------------------------------------------------------
    art::Handle<CrvRecoPulseCollection> pulsesHandle;
    Evt->getByLabel(_crvRecoPulsesModuleLabel, pulsesHandle);
    
    if (pulsesHandle.isValid()) {
      foundCRV = true;
      printf(">>> [%s] MSG: CrvRecoPulsesCollection by %s, found. CONTINUE\n", oname, _crvRecoPulsesModuleLabel.data());
      const mu2e::CrvRecoPulseCollection* fCrvPulseColl = (CrvRecoPulseCollection*) pulsesHandle.product();
      
      // Clear the map pointers in preperation to (re)fill them with new information
      fCrvPulseColl_Right->clear();
      fCrvPulseColl_Left->clear();
      fCrvPulseColl_TopDS->clear();
      fCrvPulseColl_TopTS->clear();
      fCrvPulseColl_Dwnstrm->clear();
      fCrvPulseColl_Upstrm->clear();
      printf(">>> Section-collections cleared\n");
      
      mu2e::GeomHandle<mu2e::CosmicRayShield> CRS;
      
      // Loop over the RecoPulses in the collection and sort each pulse/bar into its appropriate section
      // for (mu2e::CrvRecoPulseCollection::const_iterator icrpc = fCrvPulseColl->begin(), ecrpc = fCrvPulseColl->end(); icrpc != ecrpc; ++icrpc) {
      int    crvCollSize = fCrvPulseColl->size();
      for (int ic=0; ic < crvCollSize; ++ic) {
	mu2e::CrvRecoPulse                   icprc       = fCrvPulseColl->at(ic);
	const mu2e::CRSScintillatorBarIndex &CRVBarIndex = icprc.GetScintillatorBarIndex();
	
	switch (getCRVSection(CRS->getBar(CRVBarIndex).id().getShieldNumber()))
	  {
	  case 0:
	    fCrvPulseColl_Right->  push_back(icprc);
	    break;
	  case 1:
	    fCrvPulseColl_Left->   push_back(icprc);
	    break;
	  case 2:
	    fCrvPulseColl_TopDS->  push_back(icprc);
	    break;
	  case 3:
	    fCrvPulseColl_Dwnstrm->push_back(icprc);
	    break;
	  case 4:
	    fCrvPulseColl_Upstrm-> push_back(icprc);
	    break;
	  case 5:
	  case 6:
	  case 7:
	    break;
	  case 8:
	    fCrvPulseColl_TopTS->push_back(icprc);
	    break;
	  }
      }
								
      printf(">>> Section-collections filled\n");
    }
    else {
      printf(">>> [%s] WARNING: CrvRecoPulsesCollection by %s is missing. CONTINUE.\n",
	     oname, _crvRecoPulsesModuleLabel.data());
    }
    
    if (_showCRVOnly) { //If only displaying the CRV, skip everything else
      printf("!!!!! ONLY SHOWING THE CRV\n");
      return 0;
    }
    else {
//-----------------------------------------------------------------------------
//  MC truth - gen particles
//-----------------------------------------------------------------------------
      art::Handle<GenParticleCollection> gensHandle;
      Evt->getByLabel(_generatorModuleLabel, gensHandle);

      if (gensHandle.isValid()) _genParticleColl = gensHandle.product();
      else {
	_genParticleColl = 0;
	printf(">>> [%s] WARNING: GenParticleCollection by %s is missing. CONTINUE\n",
	       oname, _generatorModuleLabel.data());
      }

      art::Handle<PtrStepPointMCVectorCollection> mcptrHandle;
      Evt->getByLabel(_makeStrawDigiModuleLabel, "", mcptrHandle);


      if (mcptrHandle.isValid()) {
	_hits_mcptr = mcptrHandle.product();
	if (_hits_mcptr->size() <= 0) {
	  printf(">>> [%s] WARNING: PtrStepPointMCVectorCollection by %s has zero length. CONTINUE\n",
		 oname, _makeStrawDigiModuleLabel.data());
	}
      }
      else {
	printf(">>> [%s] ERROR: PtrStepPointMCVectorCollection by %s is missing. BAIL OUT\n",
	       oname, _makeStrawDigiModuleLabel.data());
	return -1;
      }

      art::Handle<StepPointMCCollection> stepsHandle;
      art::Selector getTrackerSteps(art::ProductInstanceNameSelector(_trackerStepPoints) &&
				    art::ProcessNameSelector(_processName) &&
				    art::ModuleLabelSelector(_g4ModuleLabel));
      Evt->get(getTrackerSteps, stepsHandle);

      if (stepsHandle.isValid()) _stepPointMCColl = stepsHandle.product();

      else                       _stepPointMCColl = NULL;
//-----------------------------------------------------------------------------
//  straw hit information
//-----------------------------------------------------------------------------
      art::Handle<StrawHitCollection> shH;
      Evt->getByLabel(_makeStrawHitModuleLabel, shH);


      if (shH.isValid()) fStrawHitColl = shH.product();
      else {
	printf(">>> [%s] ERROR: StrawHitCollection by %s is missing. BAIL OUT\n",
	       oname, _makeStrawHitModuleLabel.data());
	return -1;
      }

      art::Handle<mu2e::ComboHitCollection> schH;
      Evt->getByLabel(_makeStrawHitModuleLabel, schH);


      if (shH.isValid()) fShComboHitColl = schH.product();
      else {
	printf(">>> [%s] ERROR: ComboHitCollection by %s is missing. BAIL OUT\n",
	       oname, _makeStrawHitModuleLabel.data());
	return -1;
      }

      art::Handle<mu2e::ComboHitCollection> chH;
      Evt->getByLabel(_makeComboHitModuleLabel, chH);


      if (chH.isValid()) fComboHitColl = chH.product();
      else {
	printf(">>> [%s] ERROR: ComboHitCollection by %s is missing. BAIL OUT\n",
	       oname, _makeStrawHitModuleLabel.data());
	return -1;
      }

      art::Handle<StrawDigiMCCollection> handle;
      art::Selector sel_straw_digi_mc(art::ProductInstanceNameSelector("") &&
				      art::ProcessNameSelector(_processName) &&
				      art::ModuleLabelSelector(_makeStrawDigiModuleLabel));
      Evt->get(sel_straw_digi_mc, handle);
      if (handle.isValid()) {
	_strawDigiMCColl = handle.product();
	if (_strawDigiMCColl->size() <= 0) {
	  printf(">>> [%s] WARNING:StrawDigiMCCollection by %s has zero length. CONTINUE\n",
		 oname, _makeStrawDigiModuleLabel.data());
	}
      }
      else {
	printf(">>> [%s] ERROR: mu2e::StrawDigiMCCollection by %s is missing. BAIL OUT\n",
	       oname, _makeStrawDigiModuleLabel.data());
	return -1;
      }

      // art::Handle<mu2e::StrawHitPositionCollection> shpH;
      // Evt->getByLabel(fStrawHitPosMaker, shpH);


      // if (shpH.isValid()) fStrawHitPosColl = shpH.product();
      // else {
      // 	printf(">>> [%s] ERROR:StrawHitPositionCollection by %s is missing. BAIL OUT\n",
      // 	       oname, fStrawHitPosMaker.data());
      // 	return -1;
      // }

      art::Handle<mu2e::StrawHitFlagCollection> shfH;
      Evt->getByLabel(fStrawHitFlagMaker, shfH);


      if (shfH.isValid()) fStrawHitFlagColl = shfH.product();
      else {
	printf(">>> [%s] ERROR: StrawHitFlagCollection by %s is missing. BAIL OUT\n",
	       oname, fStrawHitFlagMaker.data());
	return -1;
      }
//-----------------------------------------------------------------------------
// calorimeter crystal hit data
//-----------------------------------------------------------------------------
      art::Handle<CaloCrystalHitCollection> ccHandle;
      Evt->getByLabel(fCrystalHitMaker.data(), ccHandle);


      if (ccHandle.isValid()) {
	fListOfCrystalHits = (CaloCrystalHitCollection*) ccHandle.product();
      }
      else {
	fListOfCrystalHits = NULL;
	printf(">>> [%s] ERROR: CaloCrystalHitCollection by %s is missing. BAIL OUT\n",
	       oname, fCrystalHitMaker.data());
      }
//-----------------------------------------------------------------------------
// calorimeter cluster data
//-----------------------------------------------------------------------------
      art::Handle<CaloClusterCollection> calo_cluster_handle;
      Evt->getByLabel(_caloClusterModuleLabel, _producerName, calo_cluster_handle);

      if (calo_cluster_handle.isValid()) {
	fListOfClusters = calo_cluster_handle.product();
      }
      else {
	fListOfClusters = NULL;
	printf(">>> [%s] ERROR: CaloClusterCollection by %s is missing. BAIL OUT\n",
	       oname, _caloClusterModuleLabel.data());
      }
//-----------------------------------------------------------------------------
// timepeaks 
//-----------------------------------------------------------------------------
      fTimeClusterColl = NULL;
      fTimeCluster = NULL;

      art::Handle<TimeClusterCollection> tpch;
      //      const char* charDirectionAndParticle = fDirectionAndParticle.c_str();

      if (_showTracks){
	Evt->getByLabel("CalTimePeakFinder", "", tpch);
      } 
      else {
//-----------------------------------------------------------------------------
// not sure I understand the clause
//-----------------------------------------------------------------------------
	Evt->getByLabel("CalTimePeakFinder", tpch);
      }
      if (tpch.isValid()) {
	fTimeClusterColl = tpch.product();
//-----------------------------------------------------------------------------
// find the right time peak to display - display the first one with the track
//-----------------------------------------------------------------------------
	const TimeCluster* tp;
	int ipeak = -1;
	if (fTimeClusterColl != NULL) {
	  int ntp = fTimeClusterColl->size();
	  for (int i = 0; i<ntp; i++) {
	    tp = &fTimeClusterColl->at(i);
	    // if (tp->CprIndex() >= 0) {
	    fTimeCluster = tp;
	    ipeak = i;
	    break;
	    // }
	  }
	}
	fVisManager->SetTimeCluster(ipeak);
      }
//-----------------------------------------------------------------------------
// tracking data - downstream moving electrons
//-----------------------------------------------------------------------------
      art::Handle<KalRepPtrCollection> krepHandle;
      //      Evt->getByLabel(fTrkRecoModuleLabel.data(), charDirectionAndParticle, krepHandle);
      Evt->getByLabel(fTrkRecoModuleLabel.data(), "", krepHandle);

      fNTracks[0] = 0;
      _kalRepPtrColl = NULL;
      if (krepHandle.isValid()) {
	_kalRepPtrColl = krepHandle.product();
	fNTracks[0] = _kalRepPtrColl->size();
      }
    }
    return 0;
  }

//-----------------------------------------------------------------------------
  void MuHitDisplay::Init(art::Event* Evt) {
    //    TStnCluster*    cluster;
    int             id_word, ntrk;
    TStnTrack*      track;

//-----------------------------------------------------------------------------
// initialize tracks and determine track quality
//-----------------------------------------------------------------------------
    ntrk        = 0;

    if (_showTracks){
      StntupleInitMu2eTrackBlock(fTrackBlock, Evt, 0);

      fNTracks[0] = fTrackBlock->NTracks();
      ntrk        = fNTracks[0];
    }

    for (int i = 0; i<ntrk; i++) {
      track          = fTrackBlock->Track(i);
      id_word        = fTrackID->IDWord(track);
      track->fIDWord = id_word;
    }
//-----------------------------------------------------------------------------
// initialize clusters
//-----------------------------------------------------------------------------
    StntupleInitMu2eClusterBlock(fClusterBlock, Evt, 0);
    fNClusters = fClusterBlock->NClusters();
//-----------------------------------------------------------------------------
// find momentum of the G4 particle in front of the tracker
//-----------------------------------------------------------------------------
    art::Handle<mu2e::StepPointMCCollection> handle;
    const mu2e::StepPointMCCollection*       coll(0);

    art::Selector  selector(art::ProductInstanceNameSelector("virtualdetector") &&
			    art::ProcessNameSelector("") &&
			    art::ModuleLabelSelector(_g4ModuleLabel.data()));
    Evt->get(selector, handle);

    if (handle.isValid()) coll = handle.product();
    else {
      printf(">>> ERROR in MuHitDisplay::Init: failed to locate StepPointMCCollection for virtual detectors");
      printf(". BAIL OUT. \n");
      return;
    }

    int nsteps = coll->size();

    const mu2e::StepPointMC* step;

    const CLHEP::Hep3Vector* mom;

    double px, py, pt;

    fTrackerP = -1.;
    fTrackerPt = -1.;

    for (int i = 0; i<nsteps; i++) {
      step = &coll->at(i);
      if (step->volumeId() == 13) {
	mom = &step->momentum();
	fDump->printStepPointMC(step, "banner+data");

	px = mom->x();
	py = mom->x();

	pt = sqrt(px*px + py*py);

	fTrackerP = mom->mag();
	fTrackerPt = pt;

	break;
      }
    }
  }


//-----------------------------------------------------------------------------
  bool MuHitDisplay::filter(art::Event& Evt) {
    const char* oname = "MuHitDisplay::filter";

    static int          firstCall(1);
    TText               t;
    TEllipse*           e;
    TGraph*             graph(0);

    TrackTool           *tt(0), *tg(0), *tmid(0);
    const GenParticle*  gen_signal(0);
    int                 pdg_id, ndisks, rc;

    TObjArray           list_of_ellipses;
    int                 n_displayed_hits, color(1), intime;
    // size_t              nmc;

    const int           module_color[2] = { kRed, kMagenta };

    vector<double>      xStep, yStep;

    const KalRep*       trk;
    const CaloCluster   *cl;
    double              xl, yl, event_time;
    float               q;

    TString             opt;
    Hep3Vector          pos1, mom1, pos2, mom2;

    printf("[%s] RUN: %10i EVENT: %10i\n", oname, Evt.run(), Evt.event());

    if (_showCRVOnly) {
	if (firstCall == 1) {
	  firstCall = 0;
	  InitVisManager();
	}
	fDump->SetEvent(Evt);

	rc = getData(&Evt);

	printf("Line 856 \n");

	fCanvas->cd(0);
	fCanvas->Clear();
//-----------------------------------------------------------------------------
// Draw the frame
//-----------------------------------------------------------------------------
	double plotLimits(850.);
	fCanvas->DrawFrame(-plotLimits, -plotLimits, plotLimits, plotLimits);

	fCanvas->Modified();
	fCanvas->Update();
	fVisManager->SetEvent(Evt);
	printf("Line 869: Passed SetEvent \n");
	fVisManager->DisplayEvent();
	// go into interactive mode, till '.q'
	TModule::filter(Evt);
//-----------------------------------------------------------------------------
// memory cleanup
//-----------------------------------------------------------------------------
	if (graph) delete graph;

	return true;
      }
//-----------------------------------------------------------------------------
// Geometry of the tracker envelope.
//-----------------------------------------------------------------------------
    TubsParams envelope(fTracker->getInnerTrackerEnvelopeParams());

    // Tracker calibration object.
    mu2e::ConditionsHandle<mu2e::TrackerCalibrations> trackCal("ignored");
//-----------------------------------------------------------------------------
// init VisManager - failed to do it in beginJob - what is the right place for doing it?
//-----------------------------------------------------------------------------
    if (firstCall == 1) {
      firstCall = 0;
      InitVisManager();
    }
    fDump->SetEvent(Evt);
//-----------------------------------------------------------------------------
// get event data and initialize data blocks
//-----------------------------------------------------------------------------
    rc = getData(&Evt);

    if (rc < 0) {
      printf(">>> [%s] ERROR: not all data products present, BAIL OUT\n", oname);
      return true;
    }

    Init(&Evt);
//-----------------------------------------------------------------------------
// Construct an object that ties together all of the simulated particle and hit info.
//-----------------------------------------------------------------------------
    art::ServiceHandle<mu2e::GeometryService> geom;

    if (_simParticlesWithHits) delete _simParticlesWithHits;

    _simParticlesWithHits = new SimParticlesWithHits(Evt,
						     _g4ModuleLabel,
						     _makeStrawHitModuleLabel,
						     _trackerStepPoints,
						     _minEnergyDep,
						     _minHits);
//-----------------------------------------------------------------------------
// is he finding the first particle?
//-----------------------------------------------------------------------------
    SimParticleCollection::key_type key(1);

    SimParticleInfo const *info(_simParticlesWithHits->findOrNull(key));
    const StepPointMC     *firstStep;

    if (!info) {
      printf(">>> [%s] WARNING: SimParticleInfo is missing for event %10i. CONTINUE.\n", oname, Evt.event());
      firstStep = 0;
    }
    else {
      firstStep = &info->firstStepPointMCinTracker();
    }

    static TLine*  line(0);
    static TArc*   arc(0);
    static TArc*   arccalo(0);
    static TArrow* arrow(0);
    static TBox*   box(0);

    if (line == 0) {
      line = new TLine;
      arc = new TArc;
      arccalo = new TArc;
      arrow = new TArrow;
      box = new TBox;

      box->SetFillStyle(3002);
      box->SetFillColor(4);
      box->SetLineColor(4);
    }
//-----------------------------------------------------------------------------
//  DISPLAY: 1. Draw the frame
//-----------------------------------------------------------------------------
    fCanvas->cd(0);
    fCanvas->Clear();

    double plotLimits(850.);
    fCanvas->DrawFrame(-plotLimits, -plotLimits, plotLimits, plotLimits);

    t.SetText(-800., 900., Form("[%s] RUN: %10i EVENT: %10i NTRACKS: %4i NCLUSTERS: %4i",
				oname, Evt.run(), Evt.event(), fNTracks[0], fNClusters));
    t.SetTextSize(0.02);
    t.Draw();
//-----------------------------------------------------------------------------
// Draw the inner and outer boundaries of the tracker
//-----------------------------------------------------------------------------
    arc->SetFillStyle(0);
    arc->SetLineColor(kBlack);

    arc->DrawArc(0., 0., envelope.outerRadius());
    arc->DrawArc(0., 0., envelope.innerRadius());

    if (_stepPointMCColl != NULL) {

      SortedStepPoints sortedSteps(key, *_stepPointMCColl);

      const StepPointMC* midStep = &sortedSteps.middleByZ();

      if ((midStep == 0) || (firstStep == 0)) {
	printf("[%s] **** ERROR : firstStep = %8p midstep = %8p, BAIL OUT\n",
	       oname, firstStep, midStep);
      }

      if (firstStep) {
	pos1 = firstStep->position();
	mom1 = firstStep->momentum();
      }
      else {
	pos1.set(1., 1., 1.);
	mom1.set(1., 1., 1.);
      }

      if (midStep) {
	pos2 = midStep->position();
	mom2 = midStep->momentum();
      }
      else {
	pos2.set(1., 1., 1.);
	mom2.set(1., 1., 1.);
      }
//-----------------------------------------------------------------------------
// The generated signal particle.
//-----------------------------------------------------------------------------
      // this is not necessarily correct (for mixed events) !

      if (_genParticleColl) {

	gen_signal = &_genParticleColl->at(0);
	pdg_id = gen_signal->pdgId();
	// ROOT returns charge in units of dbar-quark charge

	q = TDatabasePDG::Instance()->GetParticle(pdg_id)->Charge() / 3.;
	tt = new TrackTool(pdg_id, q, pos1, mom1, 1., Hep3Vector());
	tg = new TrackTool(pdg_id, q, gen_signal->position(), gen_signal->momentum(), 1., Hep3Vector());
	tmid = new TrackTool(pdg_id, q, pos2, mom2, 1., Hep3Vector());

	int npt = _stepPointMCColl->size();
	for (int ipt = 0; ipt<npt; ++ipt){
	  StepPointMC const& step = _stepPointMCColl->at(ipt);
	  if (step.totalEDep() > _minEnergyDep) {
	    xStep.push_back(step.position().x());
	    yStep.push_back(step.position().y());
	  }
	}
	arc->SetFillStyle(0);
	arccalo->SetFillStyle(0);
	
	arc->SetLineColor(kRed);
	arc->DrawArc(tt->xc(), tt->yc(), tt->rho());
	arc->SetLineColor(kMagenta);
	arc->DrawArc(tmid->xc(), tmid->yc(), tmid->rho());
	arc->SetLineColor(kRed);
      }
    }
//-----------------------------------------------------------------------------
// Draw disks: the inner and outer calorimeter boundaries
//-----------------------------------------------------------------------------
    double rmin, rmax;
    mu2e::GeomHandle<mu2e::DiskCalorimeter> dc;

    ndisks = dc->nDisk();

    for (int id = 0; id<ndisks; id++) {
      rmin = dc->disk(id).innerRadius();
      rmax = dc->disk(id).outerRadius();
      arccalo->SetLineColor(kGreen + id);
      arccalo->DrawArc(0., 0., rmax);
      arccalo->DrawArc(0., 0., rmin);
    }
//-----------------------------------------------------------------------------
// print reconstructed tracks - all 4 hypotheses for comparison...
//-----------------------------------------------------------------------------
    printf(">>> [%s] NTRACKS = %4i, NCLUSTERS = %4i\n", oname, fNTracks[0], fNClusters);

    opt = "data";
    if (fPrintHits) opt += "+hits";

    if (fNTracks[0] > 0) {

      double  d0, om, r, phi0, x0, y0;

      TAnaDump::Instance()->printKalRep(0, "banner"," I alg ");
      for (int i = 0; i<fNTracks[0]; i++) {
	trk = _kalRepPtrColl->at(i).get();
	printf("%2i dem ", i);
	TAnaDump::Instance()->printKalRep(trk, opt);
//-----------------------------------------------------------------------------
// also display the reconstructed track, use s=0
//-----------------------------------------------------------------------------
	//	HelixParams hel = trk->helix(0);
	d0 = trk->helix(0.).d0();
	om = trk->helix(0.).omega();
	r = fabs(1. / om);
	phi0 = trk->helix(0.).phi0();

	x0 = -(1 / om + d0)*sin(phi0);
	y0 = (1 / om + d0)*cos(phi0);
	// 	printf("[MuHitDispla::printHelixParams] d0 = %5.3f r = %5.3f phi0 = %5.3f x0 = %5.3f y0 = %5.3f\n",
	// 	       d0, r, phi0, x0, y0);

	e = new TEllipse(x0, y0, r);
	e->SetFillStyle(3001);		// make it transparent

	e->SetLineColor(kBlue - 7);
	list_of_ellipses.Add(e);
	e->Draw();
      }
    }
//-----------------------------------------------------------------------------
// Calorimeter Clusters: print and prepare to display
//-----------------------------------------------------------------------------
    event_time = -1.;
    if (fNClusters > 0) {
      TAnaDump::Instance()->printCaloCluster(0, "banner");

      for (int i = 0; i<fNClusters; i++) {
	cl = &fListOfClusters->at(i);

	TAnaDump::Instance()->printCaloCluster(cl, opt);
//-----------------------------------------------------------------------------
// event time defined by the first, most energetic cluster
// no more division into split/non-split clusters
//-----------------------------------------------------------------------------
	if (i == 0) event_time = cl->time() + 15.;
      }

      for (int i = 0; i<fNClusters; i++) {
	cl = &fListOfClusters->at(i);
//-----------------------------------------------------------------------------
// display only clusters with E > 5 MeV
//-----------------------------------------------------------------------------
	if (cl->energyDep() > 5.) {

	  xl = cl->cog3Vector().x(); // +3904.;
	  yl = cl->cog3Vector().y();

	  e = new TEllipse(xl, yl, 50.*cl->energyDep() / 100.);
	  e->SetFillStyle(3001);

	  if (geom->hasElement<mu2e::DiskCalorimeter>()) color = module_color[cl->diskId()];

	  e->SetFillColor(color);
	  e->SetLineColor(color);

	  list_of_ellipses.Add(e);

	  e->Draw();
	}
      }
    }
//-----------------------------------------------------------------------------
// Loop over straw hits. If flagBackgroundHits = true, filter them out
//-----------------------------------------------------------------------------
    int                          ihit, n_combo_hits, display_hit, index_sh(0);
    bool                         isFromConversion;
    double                       sigv, /*vnorm, v,*/ sigr;
    const ComboHit*              hit;
    const ComboHit*              sh;
    // const StrawHitPosition*      hitpos;
    const StrawHitFlag*          hit_id_word;
    // const Straw*                 straw;
    // const PtrStepPointMCVector*  mcptr;
    XYZVec                       vx0, vx1, vx2, w;

    n_displayed_hits = 0;

    if (fTimeCluster != NULL) n_combo_hits = fTimeCluster->nhits();// NHits();
    else                      n_combo_hits = fComboHitColl->size();// fStrawHitColl->size();

    for (int ih = 0; ih<n_combo_hits; ++ih) {

      if (fTimeCluster != NULL) ihit = fTimeCluster->hits().at(ih);//HitIndex(ih);
      else                      ihit = ih;

      hit         = &fComboHitColl->at(ihit);//fStrawHitColl->at(ihit);
      hit_id_word = &fStrawHitFlagColl->at(ihit);

      int    nsh = hit->nStrawHits();
      
      for (int ish=0; ish<nsh; ++ish){
	// hitpos      = &fStrawHitPosColl->at(ihit);
	index_sh    = hit->index(ish);
	sh          = &fShComboHitColl->at(index_sh);

	display_hit = 1;

	if (display_hit && (_hits_mcptr->size() > 0)) {

	  mu2e::GenId gen_id(fGeneratorID);

	  std::vector<StrawDigiIndex> shids;
	  fShComboHitColl->fillStrawDigiIndices(Evt,ish,shids);
	  mu2e::StrawDigiMC     const& mcd1  = _strawDigiMCColl->at(shids[0]);
	  art::Ptr<StepPointMC> const& spmcp = mcd1.stepPointMC(StrawEnd::cal);
	  const mu2e::SimParticle*     sim   = spmcp->simParticle().get();//_mcUtils->getSimParticle(_hits_mcptr,shids[0]);//ish);
	
	  // mcptr = &_hits_mcptr->at(ihit);

	  // Get the straw information:
	  // straw = &fTracker->getStraw(hit->strawIndex());

	  w = sh->wdir();//Geom::toXYZVec(straw->getDirection());

	  isFromConversion = false;

	  // nmc = mcptr->size();
	  // // for (size_t j = 0; j<nmc; ++j){
	  //   const art::Ptr<mu2e::StepPointMC>& sptr = mcptr->at(j);
	  //   const mu2e::StepPointMC* step = sptr.get();

	  //   //	  SimParticleCollection::key_type trackId(step->trackId());
	  //   art::Ptr<SimParticle> const& simptr = step->simParticle();
	  //   const SimParticle* sim = simptr.operator ->();
	  if (sim == NULL) {
	    printf(">>> ERROR: %s sim == NULL\n", oname);
	  }
	  else {
	    if (sim->fromGenerator()){
	      GenParticle* gen = (GenParticle*) &(*sim->genParticle());
	      //	      if ( gen->generatorId() == gen_id /*GenId::conversionGun*/ ){
	      if (gen == gen_signal) {
		isFromConversion = true;
	      }
	    }
	  }
	  // }
	  //-----------------------------------------------------------------------------
	  // Position along the wire comes from time division - but we are only displaying
	  // StrawHitPositions
	  //-----------------------------------------------------------------------------
	  // v = trackCal->TimeDiffToDistance(straw->index(), hit->dt());
	  // vnorm = v / straw->getHalfLength();
	  // if (fUseStereoHits) {
	    //-----------------------------------------------------------------------------
	    // new default, hit position errors come from StrawHitPositionCollection
	    //-----------------------------------------------------------------------------
	  sigv = sh->posRes(ComboHit::wire);//hitpos->posRes(StrawHitPosition::wire);
	  sigr = sh->posRes(ComboHit::trans);//hitpos->posRes(StrawHitPosition::trans);
	  // }
	  // else {
	  //   //-----------------------------------------------------------------------------
	  //   // old default, draw semi-random errors
	  //   //-----------------------------------------------------------------------------
	  //   sigv = trackCal->TimeDivisionResolution(sh->sid().asUint16(), vnorm) / 2.; // P.Murat
	  //   sigr = 5.; // in mm
	  // }
				
	  vx0 = sh->pos();//hitpos->pos();

	  vx1 = vx0 + sigv*w;
	  vx2 = vx0 - sigv*w;

	  intime = fabs(sh->time() - event_time) < _timeWindow;
	  //-----------------------------------------------------------------------------
	  // hit colors:
	  // -----------
	  // CE hits : red (in time), blue (out-of-time)
	  // the rest: cyan+3 (if any bad bit), otherwise - black 
	  //-----------------------------------------------------------------------------
	  int width(1);
	  if (isFromConversion) {
	    if (intime) color = kRed-3;
	    else        color = kBlue;
	    width = 2;
	  }
	  else {
	    if (hit_id_word->hasAnyProperty(fBadHitMask)) color = kCyan+3;
	    else                                          color = kBlack;
	  }
	  //-----------------------------------------------------------------------------
	  // not sure why the logic is that complicated
	  //-----------------------------------------------------------------------------
	  if ((isFromConversion) || (intime)) {
	    line->SetLineColor(color);
	    line->SetLineWidth(width);
	    line->DrawLine(vx1.x(), vx1.y(), vx2.x(), vx2.y());

	    line->DrawLine(vx0.x() + sigr*w.y(), vx0.y() - sigr*w.x(),
			   vx0.x() - sigr*w.y(), vx0.y() + sigr*w.x());

	    n_displayed_hits++;
	  }
	}
      }//end loop over the strawHits associated to the ComboHit
    }
//-----------------------------------------------------------------------------
// Draw StepPointMCCollection.
//-----------------------------------------------------------------------------
    if (xStep.size() > 0) {
      graph = new TGraph(xStep.size(), &xStep[0], &yStep[0]);
      graph->SetMarkerStyle(kFullTriangleUp);
      graph->SetMarkerColor(kRed-3);
      graph->SetMarkerSize(0.8);
      graph->Draw("PSAME");
    }
//-----------------------------------------------------------------------------
// red marker and arrow: first point on the track.
//-----------------------------------------------------------------------------
    double xf1 = pos1.x();
    double yf1 = pos1.y();
    TGraph genPoint(1, &xf1, &yf1);
    genPoint.SetMarkerColor(kRed-3);
    genPoint.SetMarkerSize(1.0);
    genPoint.SetMarkerStyle(kFullCircle);
    genPoint.Draw("PSAME");

    double arrowLength(200.);
    double xf2 = xf1 + arrowLength*mom1.x() / mom1.perp();
    double yf2 = yf1 + arrowLength*mom1.y() / mom1.perp();
    arrow->SetLineColor(kRed-3);
    arrow->DrawArrow(xf1, yf1, xf2, yf2, 0.01, ">");

    if (_genParticleColl && (_stepPointMCColl != NULL)) {
      double d0x = tt->d0x();
      double d0y = tt->d0y();
      double d0x2 = tt->d0x() + arrowLength*tt->u0();
      double d0y2 = tt->d0y() + arrowLength*tt->v0();

      arrow->SetLineColor(kBlue);
      arrow->DrawArrow(d0x, d0y, d0x2, d0y2, 0.01, ">");
    }
//-----------------------------------------------------------------------------
// blue cross - origin.
//-----------------------------------------------------------------------------
    double xo = 0.;
    double yo = 0.;
    TGraph genPointo(1, &xo, &yo);
    genPointo.SetMarkerColor(kBlue);
    genPointo.SetMarkerSize(2.5);
    genPointo.SetMarkerStyle(kPlus);
    genPointo.Draw("PSAME");

    fCanvas->Modified();
    fCanvas->Update();

    printf("N(combo hits) = %5i, N(displayed straw hits): %5i\n", n_combo_hits, n_displayed_hits);

    fVisManager->SetEvent(Evt);
    fVisManager->DisplayEvent();
//-----------------------------------------------------------------------------
// go into interactive mode, till '.q' is pressed
//-----------------------------------------------------------------------------
    TModule::filter(Evt);
//-----------------------------------------------------------------------------
// memory cleanup
//-----------------------------------------------------------------------------
    list_of_ellipses.Delete();

    if (tt) {
      delete tt;
      delete tg;
      delete tmid;
    }

    if (graph) delete graph;

    return true;
  } 

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
  int MuHitDisplay::getCRVSection(int shieldNumber) {
    int CRVSection = -1;
    switch (shieldNumber) {
      case 0:
      case 1:
      case 2:
      case 3:
      case 4:
      case 5:  CRVSection = 0; break;  //R
      case 6:
      case 7:
      case 8:  CRVSection = 1; break;  //L
      case 9:  CRVSection = 8; break;  //TS T
      case 10:
      case 11:
      case 12: CRVSection = 2; break;  //T
      case 13: CRVSection = 3; break;  //D
      case 14: CRVSection = 4; break;  //U
      case 15: CRVSection = 5; break;  //CU
      case 16: CRVSection = 6; break;  //CD
      case 17: CRVSection = 7; break;  //CT
      }
    return CRVSection;
  }

}


using mu2e::MuHitDisplay;
DEFINE_ART_MODULE(MuHitDisplay);
