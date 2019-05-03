//-----------------------------------------------------------------------------
//  Dec 26 2000 P.Murat: initialization of the STNTUPLE track block
//  2014-06-23: remove vane support
//-----------------------------------------------------------------------------
#include <cstdio>
#include "TROOT.h"
#include "TFolder.h"
#include "TLorentzVector.h"
#include "TVector2.h"
					// Mu2e 
#include "Stntuple/obj/TStnTrack.hh"
#include "Stntuple/obj/TStnTrackBlock.hh"
#include "Stntuple/obj/TStnEvent.hh"
#include "Stntuple/obj/TStnHelix.hh"
#include "Stntuple/obj/TStnHelixBlock.hh"

#include "Stntuple/obj/TStnTrackSeed.hh"
#include "Stntuple/obj/TStnTrackSeedBlock.hh"

#include "art/Framework/Principal/Selector.h"
#include "art/Framework/Principal/Handle.h"

#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/VirtualDetector.hh"
#include "GeometryService/inc/DetectorSystem.hh"

//#include "TrackerConditions/inc/Types.hh"
#include "TrackerGeom/inc/Tracker.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "CalorimeterGeom/inc/DiskCalorimeter.hh"
// #include "CalorimeterGeom/inc/VaneCalorimeter.hh"

#include "RecoDataProducts/inc/HelixSeed.hh"
#include "RecoDataProducts/inc/KalSeed.hh"
#include "RecoDataProducts/inc/KalRepPtrCollection.hh"
#include "RecoDataProducts/inc/KalRepCollection.hh"
#include "BTrkData/inc/TrkStrawHit.hh"
#include "BTrkData/inc/Doublet.hh"
#include "RecoDataProducts/inc/TrkQual.hh"
#include "TrkReco/inc/DoubletAmbigResolver.hh"

#include "RecoDataProducts/inc/TrkCaloIntersectCollection.hh"
#include "RecoDataProducts/inc/TrackClusterMatch.hh"

#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/StrawDigiMCCollection.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "DataProducts/inc/VirtualDetectorId.hh"

#include "RecoDataProducts/inc/StrawDigi.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "RecoDataProducts/inc/CaloHitCollection.hh"
#include "RecoDataProducts/inc/CaloClusterCollection.hh"
#include "RecoDataProducts/inc/PIDProductCollection.hh"
#include "RecoDataProducts/inc/AlgorithmIDCollection.hh"

					          // BaBar 
#include "BTrk/ProbTools/ChisqConsistency.hh"
#include "BTrk/BbrGeom/BbrVectorErr.hh"
#include "BTrk/BbrGeom/TrkLineTraj.hh"
#include "BTrk/TrkBase/TrkPoca.hh"
#include "BTrk/KalmanTrack/KalHit.hh"

#include "Stntuple/mod/THistModule.hh"
#include "Stntuple/base/TNamedHandle.hh"

#include "TrkDiag/inc/KalDiag.hh"

namespace {

  struct ZMap_t {
    int    fMap[44][6][2];		// 44 means "by plane"
    double fZ  [88];
  };

//-----------------------------------------------------------------------------
// Map[iplane][ipanel][il]: index of the layer in Z-ordered sequence
//-----------------------------------------------------------------------------
  void InitTrackerZMap(const mu2e::Tracker* Tracker, ZMap_t* Map) {
    int      ix, loc;
    double   z0, z1;

    int nplanes = Tracker->nPlanes();

    for (int ipl=0; ipl<nplanes; ipl++) {
      for (int isec=0; isec<6; isec++) {
	ix  = isec % 2;
	loc = 2*ipl+ix;
	Map->fMap[ipl][isec][0] = loc;
	Map->fMap[ipl][isec][1] = loc;
      }
					// form the list of Z-coordinates
      const mu2e::Straw *s0, *s1;

      s0 = &Tracker->getPlane(ipl).getPanel(0).getStraw(0);
      s1 = &Tracker->getPlane(ipl).getPanel(0).getStraw(1);
      z0 = s0->getMidPoint().z();
      z1 = s1->getMidPoint().z();

      Map->fZ[2*ipl] = (z0+z1)/2.;

      s0 = &Tracker->getPlane(ipl).getPanel(1).getStraw(0);
      s1 = &Tracker->getPlane(ipl).getPanel(1).getStraw(1);
      z0 = s0->getMidPoint().z();
      z1 = s1->getMidPoint().z();

      Map->fZ[2*ipl+1] = (z0+z1)/2.;
    }
  }

//-----------------------------------------------------------------------------
// for a given Z find closest Z-layer, returns 'Plane'
// 'Offset' is a 'face number' within the Plane
//-----------------------------------------------------------------------------
  void get_station(const mu2e::Tracker* Tracker, ZMap_t* Map, double Z, int* Plane, int* Offset) {

    double dz, dz_min(1.e10);
    int    iface(-1);
					// looks that Device == Plane 
    int nplanes = Tracker->nPlanes();
					// a plane has 2 "faces", 2 layers in each
    int nfaces  = 2*nplanes;

    for (int i=0; i<nfaces; i++) {
      dz = Map->fZ[i]-Z;
      if (fabs(dz) < dz_min) {
	iface  = i;
	dz_min = fabs(dz);
      }
    }

    *Plane   = iface / 2;
    *Offset  = iface % 2;
  }


//-----------------------------------------------------------------------------
// extrapolate track to a given Z
//-----------------------------------------------------------------------------
  double s_at_given_z(const KalRep* Krep, double Z) {

    double  ds(10.), s0, s1, s2, z0, z1, z2, dzds, sz, sz1, z01;

    s1     = Krep->firstHit()->kalHit()->hit()->fltLen();
    s2     = Krep->lastHit ()->kalHit()->hit()->fltLen();
    z1     = Krep->position(s1).z();
    z2     = Krep->position(s2).z();

    dzds   = (z2-z1)/(s2-s1);
//-----------------------------------------------------------------------------
// iterate once, choose the closest point
//-----------------------------------------------------------------------------
    if (fabs(Z-z1) > fabs(Z-z2)) {
      z0 = z2;
      s0 = s2;
    }
    else {
      z0 = z1;
      s0 = s1;
    }

    sz    = s0+(Z-z0)/dzds;

    z0     = Krep->position(sz).z();     // z0 has to be close to Z(TT_FrontPA)
    z01    = Krep->position(sz+ds).z();

    dzds   = (z01-z0)/ds;
    sz1    = sz+(Z-z0)/dzds;	          // should be good enough

    return sz1;
  }

};


//-----------------------------------------------------------------------------
Int_t StntupleInitMu2eTrackBlock  (TStnDataBlock* Block, AbsEvent* AnEvent, Int_t Mode) {
  const char*               oname = {"InitMu2eTrackBlock"};
//-----------------------------------------------------------------------------
// cached pointers, owned by the StntupleMaker_module
//-----------------------------------------------------------------------------
  static int                          initialized(0);
  static mu2e::DoubletAmbigResolver* _dar;
  static mu2e::KalDiag*              _kalDiag;
  
  int                       ntrk, ev_number, rn_number;
  TStnTrack*                track;
  TStnTrackBlock            *data(0);   

  static  ZMap_t            zmap;

  const mu2e::Tracker*     tracker;

  mu2e::AlgorithmIDCollection*             list_of_algs               (0);
  const mu2e::KalRepPtrCollection*         list_of_kreps              (0);
  const mu2e::TrkQualCollection*           list_of_trkqual            (0);
  const mu2e::StrawDigiMCCollection*       list_of_mc_straw_hits      (0);
  const mu2e::ComboHitCollection*          list_of_combo_hits         (0);
  const mu2e::TrkCaloIntersectCollection*  list_of_extrapolated_tracks(0);
  const mu2e::PIDProductCollection*        list_of_pidp               (0);

  char   algs_module_label[100], algs_description[100];
  char   krep_module_label[100], krep_description[100];
  char   trex_module_label[100], trex_description[100];
  char   trcl_module_label[100], trcl_description[100];

  char   trkq_coll_tag[100];

  char   cmbh_module_label[100], cmbh_description[100];
  char   sdmc_module_label[100], sdmc_description[100];
  char   stmc_module_label[100], stmc_description[100];
  char   calo_module_label[100], calo_description[100];
  char   pidp_module_label[100], pidp_description[100];
  char   spmc_module_label[100], spmc_description[100];

  char   module_name      [100], dar_name        [100], kaldiag_name[100];

  ev_number = AnEvent->event();
  rn_number = AnEvent->run();

  if (Block->Initialized(ev_number,rn_number)) return 0;

  mu2e::GeomHandle<mu2e::Tracker> ttHandle;
  tracker = ttHandle.get();

  data = (TStnTrackBlock*) Block;
  data->Clear();

  if (initialized == 0) {
    initialized = 1;

    InitTrackerZMap(tracker,&zmap);

    data->GetModuleLabel("DarHandle",module_name);
    data->GetDescription("DarHandle",dar_name   );

    THistModule*  m;
    TNamedHandle* nh;

    m    = static_cast<THistModule*>  (THistModule::GetListOfModules()->FindObject(module_name));
    nh   = static_cast<TNamedHandle*> (m->GetFolder()->FindObject(dar_name));
    _dar = static_cast<mu2e::DoubletAmbigResolver*> (nh->Object());

    data->GetModuleLabel("KalDiagHandle",module_name );
    data->GetDescription("KalDiagHandle",kaldiag_name);

    m        = static_cast<THistModule*>  (THistModule::GetListOfModules()->FindObject(module_name));
    nh       = static_cast<TNamedHandle*> (m->GetFolder()->FindObject(kaldiag_name));
    _kalDiag = static_cast<mu2e::KalDiag*> (nh->Object());
  }

  data->GetModuleLabel("mu2e::AlgorithmIDCollection",algs_module_label);
  data->GetDescription("mu2e::AlgorithmIDCollection",algs_description );

  data->GetModuleLabel("mu2e::KalRepCollection",krep_module_label);
  data->GetDescription("mu2e::KalRepCollection",krep_description );

  data->GetModuleLabel("mu2e::TrkCaloIntersectCollection",trex_module_label);
  data->GetDescription("mu2e::TrkCaloIntersectCollection",trex_description );

  data->GetModuleLabel("mu2e::TrackClusterMatchCollection",trcl_module_label);
  data->GetDescription("mu2e::TrackClusterMatchCollection",trcl_description );

  data->GetModuleLabel("mu2e::ComboHitCollection",cmbh_module_label);
  data->GetDescription("mu2e::ComboHitCollection",cmbh_description );
  
  data->GetModuleLabel("mu2e::PtrStepPointMCVectorCollection",stmc_module_label);
  data->GetDescription("mu2e::PtrStepPointMCVectorCollection",stmc_description );

  data->GetModuleLabel("mu2e::StrawDigiMCCollection",sdmc_module_label);
  data->GetDescription("mu2e::StrawDigiMCCollection",sdmc_description );

  data->GetModuleLabel("mu2e::CaloClusterCollection",calo_module_label);
  data->GetDescription("mu2e::CaloClusterCollection",calo_description );

  data->GetModuleLabel("mu2e::PIDProductCollection",pidp_module_label);
  data->GetDescription("mu2e::PIDProductCollection",pidp_description );

  data->GetModuleLabel("mu2e::StepPointMCCollection",spmc_module_label);
  data->GetDescription("mu2e::StepPointMCCollection",spmc_description );

  data->GetModuleLabel("mu2e::TrkQualCollection",trkq_coll_tag);

  art::Handle<mu2e::AlgorithmIDCollection> algsHandle;
  strcpy(algs_module_label, krep_module_label);
  strcpy(algs_description , krep_description);
  if (algs_module_label[0] != 0) {
    if (algs_description[0] == 0) AnEvent->getByLabel(algs_module_label,algsHandle);
    else                          AnEvent->getByLabel(algs_module_label,algs_description, algsHandle);
    if (algsHandle.isValid()) list_of_algs = (mu2e::AlgorithmIDCollection*) algsHandle.product();
  }

  art::Handle<mu2e::KalRepPtrCollection> krepsHandle;
  if (krep_module_label[0] != 0) {
    AnEvent->getByLabel(krep_module_label,krepsHandle);
    if (krepsHandle.isValid())    list_of_kreps = krepsHandle.product();
  }

  art::Handle<mu2e::TrkQualCollection> trkQualHandle;
  if (trkq_coll_tag[0] != 0) {
    AnEvent->getByLabel(trkq_coll_tag,trkQualHandle);
    if (trkQualHandle.isValid()) list_of_trkqual = trkQualHandle.product();
  }


  art::Handle<mu2e::PtrStepPointMCVectorCollection> mcptrHandle;
  if (stmc_module_label[0] != 0) {
    if (stmc_description[0] == 0) AnEvent->getByLabel(stmc_module_label,mcptrHandle);
    else                       AnEvent->getByLabel(stmc_module_label,stmc_description, mcptrHandle);
  }

  art::Handle<mu2e::ComboHitCollection> shHandle;
  if (cmbh_module_label[0] != 0) {
    if (cmbh_description[0] == 0) AnEvent->getByLabel(cmbh_module_label,shHandle);
    else                          AnEvent->getByLabel(cmbh_module_label,cmbh_description,shHandle);
    if (shHandle.isValid()) list_of_combo_hits = shHandle.product();
  }

  art::Handle<mu2e::StrawDigiMCCollection> sdmcHandle;
  if (sdmc_module_label[0] != 0) {
    AnEvent->getByLabel(sdmc_module_label,sdmcHandle);
    if (sdmcHandle.isValid()) list_of_mc_straw_hits = sdmcHandle.product();
  }

  art::Handle<mu2e::TrkCaloIntersectCollection>  texHandle;
  if (trex_module_label[0] != 0) {
    AnEvent->getByLabel(trex_module_label,texHandle);
    if (texHandle.isValid()) list_of_extrapolated_tracks = texHandle.product();
  }

  art::Handle<mu2e::TrackClusterMatchCollection>  tcmH;
  if (trcl_module_label[0] != 0) {
    if (trcl_description[0] == 0) AnEvent->getByLabel(trcl_module_label,tcmH);
    else                          AnEvent->getByLabel(trcl_module_label,trcl_description,tcmH);
  }

  art::Handle<mu2e::PIDProductCollection>  pidpHandle;
  if (pidp_module_label[0] != 0) {
    AnEvent->getByLabel(pidp_module_label,pidp_description,pidpHandle);
    if (pidpHandle.isValid()) list_of_pidp = pidpHandle.product();
  }

  art::ServiceHandle<mu2e::GeometryService>   geom;
  mu2e::GeomHandle<mu2e::DetectorSystem>      ds;
  mu2e::GeomHandle<mu2e::VirtualDetector>     vdet;

  const mu2e::AlgorithmID*  alg_id;
  int                       mask;

  const mu2e::Calorimeter* bc(NULL);
  
  if (geom->hasElement<mu2e::DiskCalorimeter>() ) {
    mu2e::GeomHandle<mu2e::DiskCalorimeter> h;
    bc = (const mu2e::Calorimeter*) h.get();
  }

  ntrk = list_of_kreps->size();
  
  // int xxx(0);

  for (int itrk=0; itrk<ntrk; itrk++) {
    track          = data->NewTrack();
    art::Handle<mu2e::KalRepCollection> handle;
    art::Ptr<KalRep> const& ptr = list_of_kreps->at(itrk);
    AnEvent->get(ptr.id(), handle);
    // fhicl::ParameterSet const& pset = handle.provenance()->parameterSet();
    // string module_type = pset.get<std::string>("module_type");
    // if      (module_type == "CalTrkFit"  ) xxx =  1;
    // else if (module_type == "KalFinalFit") xxx =  0;
    // else                                   xxx = -1;

//-----------------------------------------------------------------------------
// track-only-based particle ID, initialization ahs already happened in the constructor
//-----------------------------------------------------------------------------
    if (list_of_pidp) {
      const mu2e::PIDProduct* pidp = &list_of_pidp->at(itrk);
      track->fEleLogLHDeDx = pidp->GetLogEProb();
      track->fMuoLogLHDeDx = pidp->GetLogMProb();
      track->fRSlope       = pidp->GetResidualsSlope();
      track->fRSlopeErr    = pidp->GetResidualsSlopeError();
    }

    const KalRep* krep = ptr.get();
    track->fKalRep[0] = (KalRep*) krep;
    mask = (0x0001 << 16) | 0x0000;

    if (list_of_algs) {
      alg_id = &list_of_algs->at(itrk);
      mask   = alg_id->BestID() | (alg_id->AlgMask() << 16);

      // if (xxx != alg_id->BestID()) { 
      // 	printf (" *** InitTrackBlock ERROR: we are in alg_id trouble: xxx = %2i best = %i\n",xxx,alg_id->BestID());
      // }
    }

    track->fAlgorithmID = mask;
//-----------------------------------------------------------------------------
// in all cases define momentum at lowest Z - ideally, at the tracker entrance
// 'entlen' - trajectory length, corresponding to the first point in Z (?) 
//-----------------------------------------------------------------------------
    double  h1_fltlen, hn_fltlen, entlen;

    h1_fltlen      = krep->firstHit()->kalHit()->hit()->fltLen();
    hn_fltlen      = krep->lastHit ()->kalHit()->hit()->fltLen();
    entlen         = std::min(h1_fltlen,hn_fltlen);

    CLHEP::Hep3Vector fitmom = krep->momentum(entlen);
    HepPoint   pos           = krep->position(entlen);

    track->fX1 = pos.x();
    track->fY1 = pos.y();
    track->fZ1 = pos.z();

    double px, py, pz;
    px = fitmom.x();
    py = fitmom.y();
    pz = fitmom.z();
//-----------------------------------------------------------------------------
// track parameters in the first point
//-----------------------------------------------------------------------------
    track->Momentum()->SetXYZM(px,py,pz,0.511);
    track->fP         = fitmom.mag();
    track->fPt        = track->Momentum()->Pt();
    track->fChi2      = krep->chisq();
    track->fFitCons   = krep->chisqConsistency().consistency();
    track->fT0        = krep->t0().t0();
    track->fT0Err     = krep->t0().t0Err();
//-----------------------------------------------------------------------------
// momentum error in the first point
//-----------------------------------------------------------------------------
    CLHEP::Hep3Vector momdir = fitmom.unit();
    BbrVectorErr      momerr = krep->momentumErr(entlen);
    
    CLHEP::HepVector momvec(3);
    for (int i=0; i<3; i++) momvec[i] = momdir[i];

    track->fFitMomErr = sqrt(momerr.covMatrix().similarity(momvec));
//-----------------------------------------------------------------------------
// determine, approximately, 'sz0' - flight length corresponding to the 
// virtual detector at the tracker front
//-----------------------------------------------------------------------------
    Hep3Vector tfront = ds->toDetector(vdet->getGlobal(mu2e::VirtualDetectorId::TT_FrontPA));
    double     zfront = tfront.z();
    double     sz0    = s_at_given_z(krep,zfront);
//-----------------------------------------------------------------------------
// fP0 - track momentum value at Z(TT_FrontPA)
//-----------------------------------------------------------------------------
    CLHEP::Hep3Vector fitmom2 = krep->momentum(sz0);
    track->fP0 = fitmom2.mag();
//-----------------------------------------------------------------------------
// helical parameters at Z(TT_FrontPA)
//-----------------------------------------------------------------------------
    HelixParams helx  = krep->helix(sz0);
    track->fC0        = helx.omega(); // old
    track->fD0        = helx.d0();
    track->fZ0        = helx.z0();
    track->fPhi0      = helx.phi0();
    track->fTanDip    = helx.tanDip(); // old

    track->fCharge    = krep->charge();
//-----------------------------------------------------------------------------
// virtual detector at the tracker exit: Time at Z(TT_Back)
//-----------------------------------------------------------------------------
    Hep3Vector vd_tt_back = ds->toDetector(vdet->getGlobal(mu2e::VirtualDetectorId::TT_Back));
    double     zback      = vd_tt_back.z();
    double     szb        = s_at_given_z(krep,zback);
    double     tback      = krep->arrivalTime(szb);
					// rename later
    track->fTBack         = tback;
//-----------------------------------------------------------------------------
// fP2 : track momentum at Z(TT_Back), just for fun, should not be used for anything
//-----------------------------------------------------------------------------
    track->fP2 = krep->momentum(szb).mag();
//-----------------------------------------------------------------------------
    const mu2e::TrkStrawHit  *hit, *closest_hit(NULL);
    const TrkHitVector*       krep_hits = &krep->hitVector();
    
    for (int j=0; j<40; j++) {
      track->fNHPerStation[j] = 0;
    }
    
    int     loc, ipart, nhits, n_straw_hits, found, ntrkhits(0), nhitsambig0(0); // , pdg_code;
    int     id(-1),  npart(0), part_pdg_code[100], part_nh[100], part_id[100];
    int     nwrong = 0;
    double  mcdoca;

    const mu2e::ComboHit      *s_hit0;
    const mu2e::ComboHit      *s_hit; 
    const mu2e::SimParticle   *sim(NULL); 
    const mu2e::StepPointMC   *stmc[2];
    const mu2e::Straw         *straw;

    n_straw_hits = list_of_combo_hits->size();

    if (n_straw_hits <= 0) {
      printf(">>> ERROR in StntupleInitMu2eTrackBlock: ComboHitCollection by module %s is empty, NHITS = %i\n",
	     cmbh_module_label,n_straw_hits);
    }
    else {
      s_hit0 = &list_of_combo_hits->at(0);

      for (auto it=krep_hits->begin(); it!=krep_hits->end(); it++) {
	++ntrkhits;
	hit   = dynamic_cast<const mu2e::TrkStrawHit*> (*it);
//-----------------------------------------------------------------------------
// skip calorimeter hit
//-----------------------------------------------------------------------------
	if (! hit) continue;
	straw = &hit->straw();
//-----------------------------------------------------------------------------
// the rest makes sense only for active hits
//-----------------------------------------------------------------------------
	if (hit->isActive()) {
	  s_hit = &hit->comboHit();
	  loc   = s_hit-s_hit0;
	  if ((loc >= 0) && (loc < n_straw_hits)) {
	    if ((list_of_mc_straw_hits != NULL) && (list_of_mc_straw_hits->size() > 0)) {

	      const mu2e::StrawDigiMC* mcdigi = &list_of_mc_straw_hits->at(loc);

	      if (mcdigi->wireEndTime(mu2e::StrawEnd::cal) < mcdigi->wireEndTime(mu2e::StrawEnd::hv)) {
		stmc[0] = mcdigi->stepPointMC(mu2e::StrawEnd::cal).get();
	      }
	      else {
		stmc[0] = mcdigi->stepPointMC(mu2e::StrawEnd::hv ).get();
	      }

	      //-----------------------------------------------------------------------------
	      // count number of active hits with R > 200 um and misassigned drift signs
	      //-----------------------------------------------------------------------------
	      if (hit->driftRadius() > 0.2) {
		const CLHEP::Hep3Vector* v1 = &straw->getMidPoint();
		HepPoint p1(v1->x(),v1->y(),v1->z());
	      
		const CLHEP::Hep3Vector* v2 = &stmc[0]->position();
		HepPoint    p2(v2->x(),v2->y(),v2->z());
	      
		TrkLineTraj trstraw(p1,straw->getDirection()  ,0.,0.);
		TrkLineTraj trstep (p2,stmc[0]->momentum().unit(),0.,0.);
	      
		TrkPoca poca(trstep, 0., trstraw, 0.);
	      
		mcdoca = poca.doca();
		//-----------------------------------------------------------------------------
		// if mcdoca and hit->_iamb have different signs, the hit drift direction 
		// has wrong sign
		//-----------------------------------------------------------------------------
		if (hit->ambig()*mcdoca < 0) nwrong += 1;

		if (hit->ambig()       == 0) ++nhitsambig0;
	      }
	    
	      sim = &(*stmc[0]->simParticle());
	      if (sim != NULL) id = sim->id().asInt();
	      else {
		printf(">>> ERROR in %s : sim is NULL, set PDG_CODE to -1\n",oname);
		id = -1;
	      }

	      found = 0;
	      for (int ip=0; ip<npart; ip++) {
		if (id == part_id[ip]) {
		  found        = 1;
		  part_nh[ip] += 1;
		  break;
		}
	      }
	    
	      if (found == 0) {
		part_id      [npart] = id;
		part_pdg_code[npart] = sim->pdgId();
		part_nh      [npart] = 1;
		npart               += 1;
	      }
	    }
	  }
	  else {
	    printf(">>> ERROR in StntupleInitMu2eTrackBlock: wrong hit collection used, loc = %10i, n_straw_hits = %10i\n",
		   loc,n_straw_hits);
	  }

	  const mu2e::StrawId& straw_id = straw->id();
	
	  int ist = straw_id.getStation();
	
	  track->fNHPerStation[ist] += 1;
	
	  int pan = straw_id.getPanel();
	  int lay = straw_id.getLayer();
	  int bit = zmap.fMap[ist][pan][lay];

	  track->fHitMask.SetBit(bit,1);
	}
      }
    }
//-----------------------------------------------------------------------------
// Dave's variables calculated by KalDiag
//-----------------------------------------------------------------------------
    _kalDiag->kalDiag(krep,false);
//-----------------------------------------------------------------------------
// total number of hits associated with the trackand the number of bend sites
//-----------------------------------------------------------------------------
    track->fNHits     = ntrkhits | (_kalDiag->_trkinfo._nbend << 16);
    track->fNMatSites = _kalDiag->_trkinfo._nmat | (_kalDiag->_trkinfo._nmatactive << 16);

    if (list_of_trkqual) track->fTrkQual = list_of_trkqual->at(itrk).MVAOutput();//_kalDiag->_trkinfo._trkqual;
    else                 track->fTrkQual = -1.e6;
//-----------------------------------------------------------------------------
// defined bit-packed fNActive word
//-----------------------------------------------------------------------------
    track->fNActive   = krep->nActive() | (nwrong << 16);
    
    mu2e::Doublet*                     d;
    mu2e::DoubletAmbigResolver::Data_t r;

    int   nd, nd_tot(0), nd_os(0), nd_ss(0), ns;
    vector<mu2e::Doublet> list_of_doublets;

    _dar->findDoublets(krep,&list_of_doublets);

    nd = list_of_doublets.size();
//-----------------------------------------------------------------------------
// counting only 2+ hit doublets
//-----------------------------------------------------------------------------
    int nad(0);  // number of doublets with all hits active

    for (int i=0; i<nd; i++) {
      d  = &list_of_doublets.at(i);
      ns = d->fNStrawHits;
					
      if (ns > 1) { 
	nd_tot += 1;
	if (d->isSameSign()) nd_ss += 1;
	else                 nd_os += 1;

	int active = 1;
	for (int is=0; is<ns; is++) {
	  if (!d->fHit[is]->isActive()) {
	    active = 0;
	    break;
	  }
	}

	if (active == 1) {
	  nad += 1;
	}
      }
    }

  track->fNDoublets = nd_os | (nd_ss << 8) | (nhitsambig0 << 16) | (nad << 24);
//-----------------------------------------------------------------------------
// given track parameters, build the expected hit mask
//-----------------------------------------------------------------------------
  double z, closest_z(-1.e6), dz, zw, dz_min, s, s0;
    int    iplane, offset;
    int    nz(88);

    for (int iz=0; iz<nz; iz++) {
      z = zmap.fZ[iz];
					// find the track hit closest to that Z
      dz_min = 1.e10;
      for(auto it=krep_hits->begin(); it!=krep_hits->end(); it++) {

	hit = dynamic_cast<const mu2e::TrkStrawHit*> (*it);
	if (! hit) continue;

	s_hit = &hit->comboHit();
	loc   = s_hit-s_hit0;
	zw    = hit->straw().getMidPoint().z();
	dz    = z-zw;

	if (fabs(dz) < dz_min) {
	  dz_min      = fabs(dz);
	  closest_hit = hit;
	  closest_z   = zw;
	}
      }
//-----------------------------------------------------------------------------
// found closest hit and the extrapolation length, then extrapolate track
//-----------------------------------------------------------------------------
      s0  = closest_hit->fltLen();
      s   = (z-track->fZ0)/(closest_z-track->fZ0)*s0;

      HepPoint    pz    = krep->position(s);

      get_station(tracker,&zmap,z,&iplane,&offset);

      const mu2e::Panel*  panel0 = NULL;
      const mu2e::Panel*  panel;
      const mu2e::Plane*  plane;
      double              min_dphi(1.e10);
      double              dphi, nx, ny, wx, wy, wrho, rho, phi0;
      CLHEP::Hep3Vector   w0mid;
      int                 ipanel;

      plane = &tracker->getPlane(iplane);

      for (int i=0; i<3; i++) {
	ipanel = 2*i+offset;		// panel number
					// check if point pz(x0,y0) overlaps with the segment iseg
					// expected mask is set to zero
	panel = &plane->getPanel(ipanel);
	w0mid = panel->straw0MidPoint();
					// also calculate pho wrt the sector
	wx    = w0mid.x();
	wy    = w0mid.y();

	phi0  = w0mid.phi();            // make sure I understand the range

	wrho  = sqrt(wx*wx+wy*wy);

	nx    = wx/wrho;
	ny    = wy/wrho;

	rho = nx*pz.x()+ny*pz.y();

	dphi = TVector2::Phi_mpi_pi(phi0-pz.phi());
	if ((abs(dphi) < min_dphi) && (rho > wrho)) {
	  min_dphi = fabs(dphi);
	  panel0   = panel;
	}
      }
//-----------------------------------------------------------------------------
// OK, closest segment found, set the expected bit..
//-----------------------------------------------------------------------------
      if (panel0 != NULL) {
	track->fExpectedHitMask.SetBit(iz,1);
      }
    }
//-----------------------------------------------------------------------------
// identify track with the particle which produced most hits
//-----------------------------------------------------------------------------
    ipart = 0;
    nhits = part_nh[0];

    for (int ip=1; ip<npart; ip++) {
      if (part_nh[ip] > nhits) {
	nhits = part_nh[ip];
	ipart = ip;
      }
    }

    track->fPdgCode     = part_pdg_code[ipart];
    track->fPartID      = part_id      [ipart];
    track->fNGoodMcHits = nhits;
//-----------------------------------------------------------------------------
// particle parameters at virtual detectors
//-----------------------------------------------------------------------------
    mu2e::GeomHandle<mu2e::VirtualDetector> vdg;

    double t_front(1.e6), t_stout(1.e6);

    track->fPFront = -1.;
    track->fPStOut = -1.;

    if (vdg->nDet() > 0) {
      art::Handle<mu2e::StepPointMCCollection> vdhits;
      AnEvent->getByLabel(spmc_module_label,spmc_description,vdhits);
      if (!vdhits.isValid()) {
	char warning[100];
	sprintf(warning,"WARNING: StepPointMCCollection %s:virtualdetector not found\n",spmc_module_label);
	mf::LogWarning(oname) << warning;
      }
      else {
	int nvdhits = vdhits->size();
	for (int i=0; i<nvdhits; i++) {
	  const mu2e::StepPointMC* hit = &(*vdhits)[i];
	  
	  //int vdid = hit.volumeId();
	  mu2e::VirtualDetectorId vdid(hit->volumeId());

	  if (vdid.id() == mu2e::VirtualDetectorId::ST_Out) {

	    //	    TAnaDump::Instance()->printStepPointMC(hit,"");

	    art::Ptr<mu2e::SimParticle> const& simptr = hit->simParticle();
	    const mu2e::SimParticle* sim  = simptr.operator ->();

	    if (sim == NULL) {
	      char warning[100];
	      printf(">>> ERROR: %s sim == NULL\n",oname);
	      sprintf(warning,"WARNING: SimParticle for step %i = NULL\n",i);
	      mf::LogWarning(oname) << warning;
	    }
	    int sim_id = sim->id().asInt();
	    if ((sim_id == track->fPartID)  && (hit->time() <  t_stout)) {
	      track->fPStOut = hit->momentum().mag();
	      t_stout        = hit->time();
	    }
	  }
	  else if (vdid.isTrackerFront()) {

	    //	    TAnaDump::Instance()->printStepPointMC(hit,"");

	    art::Ptr<mu2e::SimParticle> const& simptr = hit->simParticle();
	    const mu2e::SimParticle* sim  = simptr.operator ->();

	    if (sim == NULL) {
	      printf(">>> ERROR: %s sim == NULL\n",oname);
	    }
	    int sim_id = sim->id().asInt();
	    if ((sim_id == track->fPartID) && (hit->time() < t_front)) {
	      track->fPFront = hit->momentum().mag();
	      t_front        = hit->time();
	    }
	  }
	}
      }
    }
//-----------------------------------------------------------------------------
// number of MC hits produced by the mother particle
//-----------------------------------------------------------------------------
    const mu2e::StepPointMC* step;

    track->fNMcStrawHits = 0;

    if (list_of_mc_straw_hits->size() > 0) {
      for (int i=0; i<n_straw_hits; i++) {
	const mu2e::StrawDigiMC* mcdigi = &list_of_mc_straw_hits->at(i);

	if (mcdigi->wireEndTime(mu2e::StrawEnd::cal) < mcdigi->wireEndTime(mu2e::StrawEnd::hv)) {
	  step = mcdigi->stepPointMC(mu2e::StrawEnd::cal).get();
	}
	else {
	  step = mcdigi->stepPointMC(mu2e::StrawEnd::hv ).get();
	}

	art::Ptr<mu2e::SimParticle> const& simptr = step->simParticle(); 
	art::Ptr<mu2e::SimParticle> mother        = simptr;
	while(mother->hasParent())  mother        = mother->parent();
	const mu2e::SimParticle*    sim           = mother.operator ->();

	int sim_id = sim->id().asInt();

	if (sim_id == track->fPartID) {
	  track->fNMcStrawHits += 1;
	}
      }
    }
//-----------------------------------------------------------------------------
// consider half-ready case when can't use the extrapolator yet; turn it off softly
//-----------------------------------------------------------------------------
    const mu2e::TrkCaloIntersect*  extrk;
    //    const KalRep *                 krep;
    CLHEP::Hep3Vector                     x1, x2;
    HepPoint                       extr_point;
    CLHEP::Hep3Vector                     extr_mom;
    int                            iv, next;
    TStnTrack::InterData_t*        vint;

    if (list_of_extrapolated_tracks != NULL) next = list_of_extrapolated_tracks->size();
    else                                     next = 0;
  
    for (int iext=0; iext<next; iext++) {
      extrk = &list_of_extrapolated_tracks->at(iext);
      krep  = extrk->trk().get();
      if (krep == track->fKalRep[0]) {
	if (track->fExtrk == 0) {
	  track->fExtrk = (const mu2e::TrkToCaloExtrapol*) extrk;
	}
	if (bc) { 
//-----------------------------------------------------------------------------
// store coordinates of the best intersection in a plane
//-----------------------------------------------------------------------------
	  iv   = extrk->diskId();
	  vint = &(track->fDisk[iv]);

	  if (vint->fID == -1) {
	    vint->fID        = iv;
	    vint->fExtrk     = (const mu2e::TrkToCaloExtrapol*) extrk;
	    vint->fChi2Match = 1.e6;
	  }
	  else {
	    printf("Run:Event: %i:%i %s : ADDITIONAL EXTR POINT for track %i on vane = %i\n", 
		   rn_number,ev_number,oname,itrk,iv);
	  }
	}
      }
    }
//-----------------------------------------------------------------------------
// now loop over track-cluster matches and find the right ones to associate 
// with the track
//-----------------------------------------------------------------------------
    unsigned int nm (0);

    const mu2e::TrackClusterMatchCollection* tcmcoll;

    if (tcmH.isValid()) {
      tcmcoll = tcmH.product();
      nm      = tcmcoll->size();
    }

    const mu2e::TrackClusterMatch* tcm;

    double best_chi2_match(1.e6);

    for (size_t im=0; im<nm; im++) {
      tcm   = &tcmcoll->at(im);
      extrk = tcm->textrapol();
      krep  = extrk->trk().get();
      if (krep == track->fKalRep[0]) {
	const mu2e::CaloCluster* cl = tcm->caloCluster();
	iv   = cl->diskId();
	vint = &track->fDisk[iv];
	if (bc == 0) {
	  printf(">>> ERROR: %s VANE calorimeter is not defined \n",oname);
	  continue;
	}

	x1   = bc->geomUtil().mu2eToDisk(iv,cl->cog3Vector());//toSectionFrame(iv,cl->cog3Vector());

	if ((track->fClosestCaloCluster == NULL) || (tcm->chi2() < best_chi2_match )) {
//-----------------------------------------------------------------------------
// if closest cluster has not been defined or the energy of the new one is higher
// depending on the calorimeter geometry choice either DX or DZ is meaningful
//-----------------------------------------------------------------------------
	  track->fClosestCaloCluster = cl;
	  track->fExtrk              = (const mu2e::TrkToCaloExtrapol*) extrk;
	  best_chi2_match            = tcm->chi2();
	}

	vint->fXTrk  = tcm->xtrk();
	vint->fYTrk  = tcm->ytrk();
	vint->fZTrk  = tcm->ztrk();
	vint->fTime  = tcm->ttrk();
	
	vint->fNxTrk = tcm->nx();
	vint->fNyTrk = tcm->ny();
	vint->fNzTrk = tcm->nz();

	if (vint->fCluster == 0) {
	  vint->fCluster      = cl;
	  vint->fClusterIndex = tcm->icl();
	  vint->fEnergy       = cl->energyDep();
	  vint->fXCl          = x1.x();
	  vint->fYCl          = x1.y();
	  vint->fZCl          = x1.z();
	  vint->fDt           = tcm->dt();
	  vint->fDx           = tcm->dx();
	  vint->fDy           = tcm->dy();
	  vint->fDz           = tcm->dz();
	  vint->fDu           = tcm->du();
	  vint->fDv           = tcm->dv();
	  vint->fChi2Match    = tcm->chi2();
	  vint->fChi2Time     = tcm->chi2_time();
	  vint->fIntDepth     = tcm->int_depth();
	  vint->fPath         = tcm->ds();
	  vint->fDr           = tcm->dr();
	  vint->fSInt         = tcm->sint();
	}
	else {
	  printf("%s : ADDITIONAL MATCH for track %i on vane = %i\n", oname,itrk,iv);
	}
      }
    }
//-----------------------------------------------------------------------------
// find intersections to use for electron ID, 
// in this version both VMinS and VMaxEp are the same
//-----------------------------------------------------------------------------
    double                    min_chi2_match(1.e6);
    TStnTrack::InterData_t*   v;

    track->fVMinS  = 0;
    track->fVMaxEp = 0;

    int ndisks = bc->nDisk();

    for (int iv=0; iv<ndisks; iv++) {
      v = &track->fDisk[iv];
      if (v->fID >= 0) {
	if (v->fCluster) {
	  if (v->fChi2Match < min_chi2_match) {
	    track->fVMaxEp = v;
	    track->fVMinS  = v;
	    min_chi2_match = v->fChi2Match;
	  }
	}
      }
    }
//-----------------------------------------------------------------------------
// define E/P by the first intersection, if it exists, the second in the 
// high-occupancy environment may be unreliable
//----------------------------------------------------
    track->fClusterE = -track->fP;
    track->fDt       = -1.e12;
    track->fDx       = -1.e12;
    track->fDy       = -1.e12;
    track->fDz       = -1.e12;

    if (track->fVMinS != 0) {
      if (track->fVMinS->fCluster) {
	track->fClusterE = track->fVMinS->fCluster->energyDep();
	track->fDx       = track->fVMinS->fDx;
	track->fDy       = track->fVMinS->fDy;
	track->fDz       = track->fVMinS->fDz;
	track->fDt       = track->fVMinS->fDt;
      }
      else {
//-----------------------------------------------------------------------------
// intersection with minimal S doesn't have a cluster, check MaxP
//-----------------------------------------------------------------------------
	if (track->fVMaxEp != 0) {
	  if (track->fVMaxEp->fCluster) {
	    track->fClusterE = track->fVMaxEp->fCluster->energyDep();
	    track->fDx       = track->fVMaxEp->fDx;
	    track->fDy       = track->fVMaxEp->fDy;
	    track->fDz       = track->fVMaxEp->fDz;
	    track->fDt       = track->fVMaxEp->fDt;
	  }
	}
      }
    }

    track->fEp = track->fClusterE/track->fP2;
  }
					// on return set event and run numbers
					// to mark block as initialized
  data->f_RunNumber   = rn_number;
  data->f_EventNumber = ev_number;

  return 0;
}

//-----------------------------------------------------------------------------
// 2015-04-02: this routine is not finished yet
//-----------------------------------------------------------------------------
Int_t StntupleInitMu2eTrackBlockLinks(TStnDataBlock* Block, AbsEvent* AnEvent, int Mode) 
{
  // Mu2e version, do nothing

  int    ev_number, rn_number;
  char   calo_module_label[100], calo_description[100];

  //  const mu2e::CaloClusterCollection*  list_of_clusters(0);

  ev_number = AnEvent->event();
  rn_number = AnEvent->run();

  if (! Block->Initialized(ev_number,rn_number)) return -1;
//-----------------------------------------------------------------------------
// do not do initialize links for 2nd time
//-----------------------------------------------------------------------------
  if (Block->LinksInitialized()) return 0;

  TStnTrackBlock* tb = (TStnTrackBlock*) Block;
  //  TStnEvent* ev      = tb->GetEvent();

  tb->GetModuleLabel("mu2e::CaloClusterCollection",calo_module_label);
  tb->GetDescription("mu2e::CaloClusterCollection",calo_description );

  art::Handle<mu2e::CaloClusterCollection> cch;
  if (calo_module_label[0] != 0) {
    if (calo_description[0] == 0) AnEvent->getByLabel(calo_module_label,cch);
    else                          AnEvent->getByLabel(calo_module_label,calo_description,cch);
  }
//-----------------------------------------------------------------------------
// TStnTrackBlock is already initialized
//-----------------------------------------------------------------------------
  TStnEvent*           ev;
  char                 krep_module_label[100];
  tb->GetModuleLabel("mu2e::KalRepCollection",krep_module_label);
  
  TStnTrackSeedBlock   *tsb;
  //  TStnTrackSeed*       ftrkseed, *trkseed;

  TStnTrack*           trk;

  //  mu2e::KalSeed*   kseed, *fkseed;

  char                 kseed_block_name[100], tsc_label[100], tsc_description[100] ;

  ev     = Block->GetEvent();

  tb->GetModuleLabel("TrackTsBlockName",kseed_block_name);
//-----------------------------------------------------------------------------
// 'TrackTs' collection stores, for each track, its KalSeed
//-----------------------------------------------------------------------------
  tb->GetModuleLabel("TrackTsCollTag"  ,tsc_label);
  tb->GetDescription("TrackTsCollTag"  ,tsc_description);

  art::Handle<mu2e::KalSeedCollection>  tsch;
  mu2e::KalSeedCollection*              list_of_trackSeeds;

  if (tsc_description[0] == 0) AnEvent->getByLabel(tsc_label,tsch);
  else                         AnEvent->getByLabel(tsc_label,tsc_description,tsch);
  list_of_trackSeeds = (mu2e::KalSeedCollection*) tsch.product();

  tsb    = (TStnTrackSeedBlock*) ev->GetDataBlock(kseed_block_name);

  int    ntrk     = tb->NTracks();
  int    ntrkseed = tsb->NTrackSeeds();

  const mu2e::KalSeed *ts, *tss;

  for (int i=0; i<ntrk; i++) {
    trk      = tb->Track(i);
    ts       = &list_of_trackSeeds->at(i); // seed corresponding to track # i

    int  loc(-1);
    for (int j=0; j<ntrkseed; ++j) {
      tss = tsb->TrackSeed(j)->fTrackSeed;
      if (ts == tss) {
	loc = j;
	break;
      }
    }
    
    if (loc < 0) {
      printf(">>> ERROR: %s track %i -> no TrackSeed associated\n", krep_module_label, i);
      continue;
    }
    
    trk->SetTrackSeedIndex(loc);
  }
//-----------------------------------------------------------------------------
// mark links as initialized
//-----------------------------------------------------------------------------
  Block->SetLinksInitialized();

  return 0;
}

