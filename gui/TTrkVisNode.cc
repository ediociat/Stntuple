///////////////////////////////////////////////////////////////////////////////
// May 04 2013 P.Murat
// 
// in 'XY' mode draw calorimeter clusters as circles with different colors 
// in 'Cal' mode draw every detail...
///////////////////////////////////////////////////////////////////////////////
#include "TVirtualX.h"
#include "TPad.h"
#include "TStyle.h"
#include "TVector3.h"
#include "TLine.h"
#include "TArc.h"
#include "TArrow.h"
#include "TBox.h"
// #include "art/Framework/Principal/Event.h"
// #include "art/Framework/Principal/Handle.h"

#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "TrackerConditions/inc/StrawResponse.hh"

#include "Stntuple/gui/TEvdTrack.hh"
#include "Stntuple/gui/TTrkVisNode.hh"
#include "Stntuple/gui/TEvdStraw.hh"
#include "Stntuple/gui/TEvdStrawHit.hh"
#include "Stntuple/gui/TEvdTrkStrawHit.hh"
#include "Stntuple/gui/TEvdStation.hh"
// #include "Stntuple/gui/TEvdFace.hh"
#include "Stntuple/gui/TEvdPanel.hh"
#include "Stntuple/gui/TEvdPlane.hh"
#include "Stntuple/gui/TEvdStrawTracker.hh"
#include "Stntuple/gui/TStnVisManager.hh"

#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitPositionCollection.hh"

#include "RecoDataProducts/inc/XYZVec.hh"

// #include "Stntuple/mod/TAnaDump.hh"

ClassImp(TTrkVisNode)

//_____________________________________________________________________________
TTrkVisNode::TTrkVisNode() : TVisNode("") {
}

//_____________________________________________________________________________
TTrkVisNode::TTrkVisNode(const char* name, const mu2e::TTracker* Tracker, TStnTrackBlock* TrackBlock): 
  TVisNode(name) {
  fTracker    = new TEvdStrawTracker(Tracker);
  fTrackBlock = TrackBlock;

  fArc        = new TArc;
  fEventTime  = 0;
  fTimeWindow = 1.e6;

  fListOfStrawHits = new TObjArray();
  fTimeCluster        = NULL;
  fUseStereoHits   = 0;

  fListOfTracks    = new TObjArray();
}

//_____________________________________________________________________________
TTrkVisNode::~TTrkVisNode() {
  delete fArc;
  
  delete fListOfStrawHits;

  fListOfTracks->Delete();
  delete fListOfTracks;
}


//_____________________________________________________________________________
void TTrkVisNode::Paint(Option_t* Option) {
  //

  TStnVisManager* vm = TStnVisManager::Instance();

  const char* view = vm->GetCurrentView();

  if      (strstr(view,"trkxy") != 0) PaintXY  (Option);
  else if (strstr(view,"trkrz") != 0) PaintRZ  (Option);
  else if (strstr(view,"cal"  ) != 0) PaintCal (Option);
  else {
				// what is the default?
    printf(Form("[TTrkVisNode::Paint] >>> ERROR: Unknown option %s",Option));
  }
}

//-----------------------------------------------------------------------------
int TTrkVisNode::InitEvent() {

  const char* oname = "TTrkVisNode::InitEvent";

  mu2e::GeomHandle<mu2e::TTracker> ttHandle;
  const mu2e::TTracker* tracker = ttHandle.get();

  // Tracker calibration object.

  // mu2e::ConditionsHandle<mu2e::TrackerCalibrations> trackCal("ignored");
  mu2e::ConditionsHandle<mu2e::StrawResponse> srep = mu2e::ConditionsHandle<mu2e::StrawResponse>("ignored");

  fListOfStrawHits->Delete();

  fListOfTracks->Delete();

  const mu2e::ComboHit              *hit;
  //  const XYZVec                      *hit_pos;
  //  const mu2e::StrawHitFlag          *hit_id_word;
  const mu2e::PtrStepPointMCVector  *mcptr;
  const mu2e::StrawDigiMC           *hit_digi_mc;

  TEvdStrawHit                      *evd_straw_hit; 
  const CLHEP::Hep3Vector           /**mid,*/ *w; 
  const mu2e::Straw                 *straw; 

  int                               n_straw_hits, display_hit, color, nl, ns; // , ipeak, ihit;
  bool                              isFromConversion, intime;
  size_t                            nmc;
  double                            sigw, /*vnorm, v,*/ sigr; 
  CLHEP::Hep3Vector                 vx0, vx1, vx2;
//-----------------------------------------------------------------------------
// first, clear the cached hit information from the previous event
//-----------------------------------------------------------------------------
  TEvdStation*   station;
  TEvdPlane*     plane;
  TEvdPanel*     panel;
  //  TEvdFace*      face;

  int            nst, nplanes, npanels/*, isec*/; 

  nst = tracker->nStations();
  for (int ist=0; ist<nst; ist++) {
    station = fTracker->Station(ist);
    nplanes = station->NPlanes();
    for (int iplane=0; iplane<nplanes; iplane++) {
      plane = station->Plane(iplane);
      npanels = plane->NPanels();
      for (int ipanel=0; ipanel<npanels; ipanel++) {
	panel = plane->Panel(ipanel);
	nl    = panel->NLayers();
	for (int il=0; il<nl; il++) {
	  ns = panel->NStraws(il);
	  for (int is=0; is<ns; is++) {
	    panel->Straw(il,is)->Clear();
	  }
	}
      }
    }
  }

  fListOfTracks->Delete();
//-----------------------------------------------------------------------------
// display hits corresponding to a given time peak, or all hits, 
// if the time peak is not found
//-----------------------------------------------------------------------------
  TEvdStraw* evd_straw;
  n_straw_hits = (*fStrawHitColl)->size();

  for (int ihit=0; ihit<n_straw_hits; ihit++ ) {

    hit         = &(*fStrawHitColl)    ->at(ihit);
    //    hit_pos     = &hit->pos();
    //    hit_pos     = &(*fStrawHitPosColl) ->at(ihit);
    //    hit_id_word = &(*fStrawHitFlagColl)->at(ihit);

    if ((*fStrawDigiMCColl)->size() > 0) hit_digi_mc = &(*fStrawDigiMCColl)->at(ihit);
    else                                 hit_digi_mc = NULL; // normally, should not be happening, but it does

    straw       = &tracker->getStraw(hit->strawId());//strawIndex());
    display_hit = 1;

    if (! display_hit)                                      continue;
//-----------------------------------------------------------------------------
// deal with MC information - later
//-----------------------------------------------------------------------------
    if ((*fStrawDigiMCColl)->size() > 0) mcptr = &(*fMcPtrColl)->at(ihit); // this seems to be wrong
    else                                 mcptr = NULL;
    // Get the straw information:

    //    mid   = &straw->getMidPoint();
    w     = &straw->getDirection();

    isFromConversion = false;

    if (mcptr != NULL) nmc = mcptr->size();
    else               nmc = 0;

    for (size_t j=0; j<nmc; ++j ){
      const mu2e::StepPointMC& step = *mcptr->at(j);
      art::Ptr<mu2e::SimParticle> const& simptr = step.simParticle();
      mu2e::SimParticleCollection::key_type trackId(step.trackId());
      const mu2e::SimParticle* sim  = simptr.get();
      if (sim == NULL) {
	printf(">>> ERROR: %s sim == NULL\n",oname);
      }
      else {
	if ( sim->fromGenerator() ){
	  mu2e::GenParticle* gen = (mu2e::GenParticle*) &(*sim->genParticle());
	  //	    if ( gen->generatorId() == mu2e::GenId::conversionGun ){
	  if ( gen->generatorId() == mu2e::GenId::StoppedParticleReactionGun ){
	    isFromConversion = true;
	    break;
	  }
	}
      }
    }
//-----------------------------------------------------------------------------
// Position along wire displayed are StrawHitPosition's
//-----------------------------------------------------------------------------	
    // v     = trackCal->TimeDiffToDistance( straw->index(), hit->dt() );
    // vnorm = v/straw->getHalfLength();

    if (fUseStereoHits) {
//-----------------------------------------------------------------------------
// new default, hit position errors come from StrawHitPositionCollection
//-----------------------------------------------------------------------------
      sigw  = hit->wireRes(); 
      sigr  = hit->transRes(); 
    }
    else {
//-----------------------------------------------------------------------------
// old default, draw semi-random errors
//-----------------------------------------------------------------------------
      sigw  = hit->wireRes()/2.; // P.Murat
      sigr  = 5.; // in mm
    }
	
    intime = fabs(hit->time()-fEventTime) < fTimeWindow;
	
    if ( isFromConversion ) {
      if (intime) color = kRed;
      else        color = kBlue;
    }
    else          color = kBlack;
//-----------------------------------------------------------------------------
// add a pointer to the hit to the straw 
//-----------------------------------------------------------------------------
    int mask = 0;
    if (intime          ) mask |= TEvdStrawHit::kInTimeBit;
    if (isFromConversion) mask |= TEvdStrawHit::kConversionBit;
    
    int ist, ipl, ippl, /*ifc,*/ ipn, il, is;

    ipl = straw->id().getPlane(); // plane number here runs from 0 to 2*NStations-1

    //    int idd = straw->id().getDevice();	// really it is a plane number in a throughout enumeration

    ist  = ipl / 2 ; // *** should become straw->id().getStation();
    ippl = ipl % 2 ; // plane number within the station

    //    ipl = idd % 2 ; // straw->id().getSector();

    ipn = straw->id().getPanel();

//     ifc  = isec%2;
//     ipn  = isec/2;

    il   = straw->id().getLayer();
    is   = straw->id().getStraw();

    // from now on, the straw number in a layer is is' = is/2
      
//     ifc = -1; //## FIXME
//     ipn = -1;

    evd_straw     = fTracker->Station(ist)->Plane(ippl)->Panel(ipn)->Straw(il,is/2);

    evd_straw_hit = new TEvdStrawHit(hit,
				     evd_straw,
				     hit_digi_mc,
				     hit->pos().x(),
				     hit->pos().y(),
				     hit->pos().z(),
				     w->x(),w->y(),
				     sigw,sigr,
				     mask,color);

    evd_straw->AddHit(evd_straw_hit);

    fListOfStrawHits->Add(evd_straw_hit);
  }
//-----------------------------------------------------------------------------
// hit MC truth unformation from StepPointMC's
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// now initialize tracks
//-----------------------------------------------------------------------------
  TEvdTrack                *trk;
  const KalRep             *krep;  
  const mu2e::TrkStrawHit  *track_hit;

  int ntrk = 0;

  if ( (*fKalRepPtrColl) != 0 )ntrk = (*fKalRepPtrColl)->size();
  
  for (int i=0; i<ntrk; i++) {
    krep = (*fKalRepPtrColl)->at(i).get();
    trk  = new TEvdTrack(i,krep);
//-----------------------------------------------------------------------------
// add hits
//-----------------------------------------------------------------------------
    const TrkHitVector* hits = &krep->hitVector();
    for (auto it=hits->begin(); it!=hits->end(); it++) {
      track_hit = (const mu2e::TrkStrawHit*) (*it);
      TEvdTrkStrawHit* h = new TEvdTrkStrawHit(track_hit);
      trk->AddHit(h);
    }

    fListOfTracks->Add(trk);
  }
  
  return 0;
}

//-----------------------------------------------------------------------------
// draw reconstructed tracks and straw hits
//-----------------------------------------------------------------------------
void TTrkVisNode::PaintXY(Option_t* Option) {

  double        time;
  int           station, display_hit, ntrk(0);
  TEvdStrawHit  *hit;

  const mu2e::ComboHit   *straw_hit;
  const mu2e::Straw      *straw; 

  //  const char* view = TVisManager::Instance()->GetCurrentView();

  mu2e::GeomHandle<mu2e::TTracker> ttHandle;
  const mu2e::TTracker* tracker = ttHandle.get();

  TStnVisManager* vm = TStnVisManager::Instance();

  int ipeak = vm->TimeCluster();

  if (ipeak >= 0) {
    if ((*fTimeClusterColl) != NULL) {
      int ntp = (*fTimeClusterColl)->size();
      if (ipeak < ntp) fTimeCluster = &(*fTimeClusterColl)->at(ipeak);
      else             fTimeCluster = NULL;
    }
  }

  double tMin = fTimeCluster->t0().t0() - 30;//FIXME!
  double tMax = fTimeCluster->t0().t0() + 20;//FIXME!
  
  int nhits = fListOfStrawHits->GetEntries();
  for (int i=0; i<nhits; i++) {
    hit       = (TEvdStrawHit*) fListOfStrawHits->At(i);
    straw_hit = hit->StrawHit();
    straw     = &tracker->getStraw(straw_hit->strawId());//strawIndex());
    station   = straw->id().getStation();
    time      = straw_hit->time();

    if ((station >= vm->MinStation()) && (station <= vm->MaxStation())) { 
      display_hit = 1;
      if (fTimeCluster != NULL) {
	if ((time < tMin) || (time > tMax)) display_hit = 0;
      }
      if (display_hit) {
	hit->Paint(Option);
      }
    }
  }

//-----------------------------------------------------------------------------
// now - tracks
//-----------------------------------------------------------------------------
  TEvdTrack* evd_trk;
  //  TAnaDump::Instance()->printKalRep(0,"banner");

  if ( (fListOfTracks) != 0 )  ntrk = fListOfTracks->GetEntriesFast();

  for (int i=0; i<ntrk; i++ ) {
    evd_trk = (TEvdTrack*) fListOfTracks->At(i);
    evd_trk->Paint(Option);
  }

//-----------------------------------------------------------------------------
// seedfit, if requested - not implemented yet
//-----------------------------------------------------------------------------
//   TAnaDump::Instance()->printKalRep(0,"banner");

//   ntrk = fListOfTracks->GetEntriesFast();
//   for (int i=0; i<ntrk; i++ ) {
//     evd_trk = (TEvdTrack*) fListOfTracks->At(i);
//     evd_trk->Paint(Option);
//   }

  gPad->Modified();
}

//-----------------------------------------------------------------------------
void TTrkVisNode::PaintRZ(Option_t* Option) {
  int             ntrk(0), nhits;
  TEvdTrack*      evd_trk;

  //  TStnVisManager* vm = TStnVisManager::Instance();

  fTracker->PaintRZ(Option);
//-----------------------------------------------------------------------------
// do not draw all straw hits - just redraw straws in different color instead
//-----------------------------------------------------------------------------
//   int nhits = fListOfStrawHits->GetEntries();
//   for (int i=0; i<nhits; i++) {
//     hit       = (TEvdStrawHit*) fListOfStrawHits->At(i);

//     if ((station >= vm->MinStation()) && (station <= vm->MaxStation())) continue;
//     if ((time    <  vm->TMin()      ) || (time     > vm->TMax()      )) continue; 

//     hit->Paint(Option);
//   }
//-----------------------------------------------------------------------------
// display tracks and track hits
//-----------------------------------------------------------------------------
  if (fListOfTracks != 0)  ntrk = fListOfTracks->GetEntriesFast();

  for (int i=0; i<ntrk; i++ ) {
    evd_trk = (TEvdTrack*) fListOfTracks->At(i);
    evd_trk->Paint(Option);

    nhits = evd_trk->NHits();
    for (int ih=0; ih<nhits; ih++) {
      TEvdTrkStrawHit* hit = evd_trk->Hit(ih);
      hit->PaintRZ(Option);
    }
  }

  gPad->Modified();
}

//_____________________________________________________________________________
void TTrkVisNode::PaintCal(Option_t* option) {
  // so far do nothing
 
}

//_____________________________________________________________________________
Int_t TTrkVisNode::DistancetoPrimitive(Int_t px, Int_t py) {
  return 9999;
}

//_____________________________________________________________________________
Int_t TTrkVisNode::DistancetoPrimitiveXY(Int_t px, Int_t py) {

  Int_t dist = 9999;

  static TVector3 global;
  //  static TVector3 local;

  //  Double_t    dx1, dx2, dy1, dy2, dx_min, dy_min, dr;

  global.SetXYZ(gPad->AbsPixeltoX(px),gPad->AbsPixeltoY(py),0);

  return dist;
}

//_____________________________________________________________________________
Int_t TTrkVisNode::DistancetoPrimitiveRZ(Int_t px, Int_t py) {
  return 9999;
}

