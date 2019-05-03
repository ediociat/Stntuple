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

#include "DataProducts/inc/StrawId.hh"
#include "DataProducts/inc/XYZVec.hh"

// #include "Stntuple/mod/TAnaDump.hh"

ClassImp(TTrkVisNode)

//_____________________________________________________________________________
TTrkVisNode::TTrkVisNode() : TVisNode("") {
}

//_____________________________________________________________________________
TTrkVisNode::TTrkVisNode(const char* name, const mu2e::Tracker* Tracker, TStnTrackBlock* TrackBlock): 
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

  //  const char* oname = "TTrkVisNode::InitEvent";

  mu2e::GeomHandle<mu2e::Tracker> ttHandle;
  const mu2e::Tracker* tracker = ttHandle.get();

  // Tracker calibration object.
  //  mu2e::ConditionsHandle<mu2e::StrawResponse> srep = mu2e::ConditionsHandle<mu2e::StrawResponse>("ignored");

  fListOfStrawHits->Delete();

  fListOfTracks->Delete();

  const mu2e::ComboHit              *hit;
  const mu2e::StrawDigiMC           *hit_digi_mc;

  TEvdStrawHit                      *evd_straw_hit; 
  const CLHEP::Hep3Vector           /**mid,*/ *w; 
  const mu2e::Straw                 *straw; 

  int                               n_straw_hits, display_hit, color, nl, ns; // , ipeak, ihit;
  bool                              isFromConversion, intime;
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

  nst = mu2e::StrawId::_nstations; // tracker->nStations();
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
  n_straw_hits = (*fComboHitColl)->size();

  for (int ihit=0; ihit<n_straw_hits; ihit++ ) {

    hit         = &(*fComboHitColl)->at(ihit);

    if ((*fStrawDigiMCColl)->size() > 0) hit_digi_mc = &(*fStrawDigiMCColl)->at(ihit);
    else                                 hit_digi_mc = NULL; // normally, should not be happening, but it does

    straw       = &tracker->getStraw(hit->strawId());
    display_hit = 1;

    if (! display_hit)                                      continue;
//-----------------------------------------------------------------------------
// deal with MC information - later
//-----------------------------------------------------------------------------
    const mu2e::StrawDigiMC             *mcdigi(0);
    if ((*fStrawDigiMCColl)->size() > 0) mcdigi = &(*fStrawDigiMCColl)->at(ihit); // this seems to be wrong
    // Get the straw information:

    w     = &straw->getDirection();

    isFromConversion = false;

    const mu2e::StepPointMC   *stmc;
    if (mcdigi->wireEndTime(mu2e::StrawEnd::cal) < mcdigi->wireEndTime(mu2e::StrawEnd::hv)) {
      stmc = mcdigi->stepPointMC(mu2e::StrawEnd::cal).get();
    }
    else {
      stmc = mcdigi->stepPointMC(mu2e::StrawEnd::hv ).get();
    }

    const mu2e::SimParticle* sim = &(*stmc->simParticle());
    
    if ( sim->fromGenerator() ){
      mu2e::GenParticle* gen = (mu2e::GenParticle*) &(*sim->genParticle());
      //	    if ( gen->generatorId() == mu2e::GenId::conversionGun ){
      if ( gen->generatorId() == mu2e::GenId::StoppedParticleReactionGun ){
	isFromConversion = true;
      }
    }

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

    ipl  = straw->id().getPlane(); // plane number here runs from 0 to 2*NStations-1
    ist  = straw->id().getStation();

    ippl = ipl % 2 ; // plane number within the station

    ipn = straw->id().getPanel();

    il   = straw->id().getLayer();
    is   = straw->id().getStraw();

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
  int           station, ntrk(0);
  TEvdStrawHit  *hit;

  const mu2e::ComboHit   *straw_hit;
  const mu2e::Straw      *straw; 

  //  const char* view = TVisManager::Instance()->GetCurrentView();

  mu2e::GeomHandle<mu2e::Tracker> ttHandle;
  const mu2e::Tracker* tracker = ttHandle.get();

  TStnVisManager* vm = TStnVisManager::Instance();

  int ipeak = vm->TimeCluster();

  if (ipeak >= 0) {
    if ((*fTimeClusterColl) != NULL) {
      int ntp = (*fTimeClusterColl)->size();
      if (ipeak < ntp) fTimeCluster = &(*fTimeClusterColl)->at(ipeak);
      else             fTimeCluster = NULL;
    }
  }

  double tmin(0), tmax(2000.);

  if (fTimeCluster) {
    tmin = fTimeCluster->t0().t0() - 30;//FIXME!
    tmax = fTimeCluster->t0().t0() + 20;//FIXME!
  }
  
  int nhits = fListOfStrawHits->GetEntries();
  for (int i=0; i<nhits; i++) {
    hit       = (TEvdStrawHit*) fListOfStrawHits->At(i);
    straw_hit = hit->StrawHit();
    straw     = &tracker->getStraw(straw_hit->strawId());//strawIndex());
    station   = straw->id().getStation();
    time      = straw_hit->time();

    if ((station >= vm->MinStation()) && (station <= vm->MaxStation())) { 
      if ((time >= tmin) && (time <= tmax)) {
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

