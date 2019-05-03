//

#include "mod/TAnaDump.hh"
#include "TROOT.h"
#include "TVector2.h"


#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Selector.h"

#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"

#include "TrackerGeom/inc/Tracker.hh"
//  #include "CalorimeterGeom/inc/VaneCalorimeter.hh"
#include "CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"

#include "RecoDataProducts/inc/CaloCluster.hh"
#include "RecoDataProducts/inc/CaloClusterCollection.hh"
#include "RecoDataProducts/inc/CaloProtoCluster.hh"
#include "RecoDataProducts/inc/CaloProtoClusterCollection.hh"
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "RecoDataProducts/inc/CaloCrystalHit.hh"
#include "RecoDataProducts/inc/CaloHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitFlagCollection.hh"
#include "RecoDataProducts/inc/KalSeed.hh"
#include "RecoDataProducts/inc/HelixSeed.hh"
#include "RecoDataProducts/inc/ComboHit.hh"
#include "RecoDataProducts/inc/TimeCluster.hh"

#include "RecoDataProducts/inc/CaloDigi.hh"
#include "RecoDataProducts/inc/CaloDigiCollection.hh"
#include "RecoDataProducts/inc/CaloRecoDigi.hh"
#include "RecoDataProducts/inc/CaloRecoDigiCollection.hh"

#include "MCDataProducts/inc/GenParticle.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/SimParticle.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/StepPointMC.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/StrawHitMCTruth.hh"
#include "MCDataProducts/inc/StrawHitMCTruthCollection.hh"
#include "MCDataProducts/inc/PhysicalVolumeInfo.hh"
#include "MCDataProducts/inc/PhysicalVolumeInfoCollection.hh"
#include "CaloMC/inc/ClusterContentMC.hh"
#include "CaloMC/inc/CrystalContentMC.hh"

#include "BTrkData/inc/TrkCaloHit.hh"
#include "BTrkData/inc/TrkStrawHit.hh"
#include "RecoDataProducts/inc/KalRepPtrCollection.hh"

#include "RecoDataProducts/inc/TrkToCaloExtrapol.hh"
#include "RecoDataProducts/inc/TrkCaloIntersectCollection.hh"
#include "RecoDataProducts/inc/TrackClusterMatch.hh"
#include "TrkDiag/inc/TrkMCTools.hh"

#include "Mu2eUtilities/inc/LsqSums4.hh"

#include "Stntuple/base/TNamedHandle.hh"

#include "Mu2eUtilities/inc/SimParticleTimeOffset.hh"
#include "TrkDiag/inc/TrkMCTools.hh"


//BaBar includes
#include "BTrk/BbrGeom/TrkLineTraj.hh"
#include "BTrk/TrkBase/TrkPoca.hh"
#include "BTrk/KalmanTrack/KalHit.hh"
#include "BTrk/KalmanTrack/KalRep.hh"
#include "BTrk/TrkBase/HelixParams.hh"
#include "BTrk/ProbTools/ChisqConsistency.hh"

using namespace std;

ClassImp(TAnaDump)

TAnaDump* TAnaDump::fgInstance = 0;

void TAnaDump::evalHelixInfo(const mu2e::HelixSeed*         Helix,
			     int   &NLoops,
			     int   &NHitsLoopFailed){
  const mu2e::ComboHit*     hit(0);

  //init
  NLoops          = 0;
  NHitsLoopFailed = 0;
  
  int           nMinHitsLoop(3), nHitsLoop(0), nHitsLoopChecked(0);
  float         meanHitRadialDist(0), z_first_hit(0), z_last_hit(0), counter(0); 
  bool          isFirst(true);
  float         half_pitch  =  M_PI*fabs(Helix->_helix._lambda);
  float         dz_min_toll = 600.;
  unsigned      nhits = Helix->_hhits.size();

  for (unsigned f=0; f<nhits; ++f){
    hit = &Helix->_hhits[f];
    if (hit->_flag.hasAnyProperty(mu2e::StrawHitFlag::outlier))     continue;
      
    meanHitRadialDist += sqrtf(hit->pos().x()*hit->pos().x() + hit->pos().y()*hit->pos().y());
    ++counter;
    float z = hit->pos().z();
    if (isFirst){
      z_first_hit = z;
      z_last_hit  = z;
      nHitsLoop   = 1;
      isFirst     = false;
    }else {
      float    dz_last_hit  = z - z_last_hit;
      float    dz_first_hit = z - z_first_hit;

      if ( ( dz_first_hit < half_pitch) && ( dz_last_hit < dz_min_toll)){
	++nHitsLoop;
	z_last_hit        = z;
      } else {
	if (nHitsLoop >= nMinHitsLoop) {
	  ++NLoops;
	  nHitsLoopChecked +=  nHitsLoop;
	}
	nHitsLoop = 0;

	if ( (z - z_last_hit) >= half_pitch){
	  //reset the hits-per-loop counter
	  nHitsLoop = 1;
	      
	  //re-set the values of the first and last z of the hits within the loop
	  z_first_hit = z;
	  z_last_hit  = z;
	}
      }
    }
	
  }//end loop over the hits

  if (counter > 0) meanHitRadialDist /= counter;
  if (nHitsLoop >= nMinHitsLoop) {
    ++NLoops;
    nHitsLoopChecked +=  nHitsLoop;
  }
      
  NHitsLoopFailed   =  nHitsLoopChecked;
}

  

double TAnaDump::evalWeight(const mu2e::ComboHit* Hit   ,
			    XYZVec&   StrawDir ,
			    XYZVec&   HelCenter, 
			    double             Radius   ,
			    int                WeightMode,
			    fhicl::ParameterSet const& Pset) {//WeightMode = 1 is for XY chi2 , WeightMode = 0 is for Phi-z chi2
  
  // double    rs(2.5);   // straw radius, mm
  // double    ew(30.0);  // assumed resolution along the wire, mm
  double    transErr = 5./sqrt(12.);
  //scale the error based on the number of the strawHits that are within teh ComboHit
  if (Hit->nStrawHits() > 1) transErr *= 1.5;
  double    transErr2 = transErr*transErr;

  double x  = Hit->pos().x();
  double y  = Hit->pos().y();
  double dx = x-HelCenter.x();
  double dy = y-HelCenter.y();
  
  double costh  = (dx*StrawDir.x()+dy*StrawDir.y())/sqrt(dx*dx+dy*dy);
  double costh2 = costh*costh;
  double sinth2 = 1-costh2;
  
  double wt(0), wtXY(1), wtPhiZ(1);

  //  fhicl::ParameterSet const& pset = helix_handle.provenance()->parameterSet();
  string                module     = Pset.get<string>("module_type");
  
  if ( module == "CalHelixFinder"){
    fhicl::ParameterSet const& psetHelFit = Pset.get<fhicl::ParameterSet>("HelixFinderAlg", fhicl::ParameterSet());
    wtXY   = psetHelFit.get<double>("weightXY");
    wtPhiZ = psetHelFit.get<double>("weightZPhi");
  }
                                            //scale the weight for having chi2/ndof distribution peaking at 1
  if ( WeightMode == 1){//XY-Fit
    // double e2     = ew*ew*sinth2+rs*rs*costh*costh;
    double e2     = Hit->wireErr2()*sinth2+transErr2*costh2;
    wt  = 1./e2;
    wt *= wtXY;
  } else if (WeightMode ==0 ){//Phi-Z Fit
    // double e2     = ew*ew*costh*costh+rs*rs*sinth2;
    double e2     = Hit->wireErr2()*costh2+transErr2*sinth2;
    wt     = Radius*Radius/e2;
    wt    *= wtPhiZ;
  }
  
  return wt;
}

//-----------------------------------------------------------------------------
TAnaDump::TAnaDump(const fhicl::ParameterSet* PSet) {

  fEvent                  = nullptr;
  fListOfObjects          = new TObjArray();
  fFlagBgrHitsModuleLabel = "FlagBkgHits";
  fStrawDigiMCCollTag     = "compressDigiMCs";

  if (PSet) {
    fhicl::ParameterSet to_maps = PSet->get<fhicl::ParameterSet>("timeOffsetMaps" );
    fTimeOffsets = new mu2e::SimParticleTimeOffset(to_maps);
  }
  else {
    fTimeOffsets = NULL;
  }

  _printUtils = new mu2e::TrkPrintUtils(PSet->get<fhicl::ParameterSet>("printUtils",fhicl::ParameterSet()));
}

//------------------------------------------------------------------------------
TAnaDump* TAnaDump::Instance(const fhicl::ParameterSet* PSet) {
  static TAnaDump::Cleaner cleaner;

  if  (! fgInstance) fgInstance  = new TAnaDump(PSet);
  return fgInstance;
}


//______________________________________________________________________________
TAnaDump::~TAnaDump() {
  fListOfObjects->Delete();
  delete fListOfObjects;
  delete fTimeOffsets;
  delete _printUtils;
}

//------------------------------------------------------------------------------
TAnaDump::Cleaner::Cleaner() {
}

//------------------------------------------------------------------------------
  TAnaDump::Cleaner::~Cleaner() {
    if (TAnaDump::fgInstance) {
      delete TAnaDump::fgInstance;
      TAnaDump::fgInstance = 0;
    }
  }


//-----------------------------------------------------------------------------
void TAnaDump::AddObject(const char* Name, void* Object) {
  TNamedHandle* h = new TNamedHandle(Name,Object);
  fListOfObjects->Add(h);
}

//-----------------------------------------------------------------------------
void* TAnaDump::FindNamedObject(const char* Name) {
  void* o(NULL);

  TNamedHandle* h = (TNamedHandle*) fListOfObjects->FindObject(Name);
  if (h != NULL) {
    o = h->Object();
  }
  return o;
}

//-----------------------------------------------------------------------------
// print position of the cluster in the tracker system
//-----------------------------------------------------------------------------
void TAnaDump::printCaloCluster(const mu2e::CaloCluster* Cl, 
				const char* Opt, 
				const mu2e::CaloHitMCTruthAssns* CaloHitTruth) {
  int row, col;
  TString opt = Opt;

  art::ServiceHandle<mu2e::GeometryService> geom;
  mu2e::GeomHandle  <mu2e::Calorimeter>     cal;
  Hep3Vector        gpos, tpos;

  if ((opt == "") || (opt.Index("banner") >= 0)) {
    printf("-----------------------------------------------------------------------------------------------");
    printf("-------------------------------\n");
    printf(" Row Col Address        Disk Parent  NC   Energy   Time       X(loc)     Y(loc)   Z(loc)");
    printf("        X          Y          Z\n");
    printf("-----------------------------------------------------------------------------------------------");
    printf("-------------------------------\n");
  }
 
  if ((opt == "") || (opt.Index("data") >= 0)) {
    row = -1; // Cl->cogRow();
    col = -1; // Cl->cogColumn();
    
    const mu2e::CaloCluster::CaloCrystalHitPtrVector caloClusterHits = Cl->caloCrystalHitsPtrVector();
    int nh = caloClusterHits.size();

    if ((row < 0) || (row > 9999)) row = -99;
    if ((col < 0) || (col > 9999)) col = -99;
//-----------------------------------------------------------------------------
// transform cluster coordinates to the tracker coordiante system
//-----------------------------------------------------------------------------
    gpos = cal->geomUtil().diskToMu2e(Cl->diskId(),Cl->cog3Vector());
    tpos = cal->geomUtil().mu2eToTracker(gpos);

    printf(" %3i %3i %-16p %2i %6i %3i %8.3f %8.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f\n",
	   row, col,
	   Cl,
	   Cl->diskId(),
	   -999, 
	   nh,
	   Cl->energyDep(),
	   Cl->time(),
	   Cl->cog3Vector().x(),
	   Cl->cog3Vector().y(),
	   Cl->cog3Vector().z(),
	   tpos.x(),
	   tpos.y(),
	   tpos.z() 
	   );
  }
  
  if (opt.Index("hits") >= 0) {

    const mu2e::Crystal           *cr;
    const CLHEP::Hep3Vector       *pos;
    // int  iz, ir;
//-----------------------------------------------------------------------------
// print individual crystals in local vane coordinate system
//-----------------------------------------------------------------------------
    const mu2e::CaloCluster::CaloCrystalHitPtrVector caloClusterHits = Cl->caloCrystalHitsPtrVector();
    int nh = caloClusterHits.size();
    
    mu2e::Calorimeter const & calo = *(mu2e::GeomHandle<mu2e::Calorimeter>());

    printf("-----------------------------------------------------------------------------------------------");
    printf("-------------------------------\n");
    printf("    Id       time      Gen-Code     ID   PDG  PDG(M)      energy       X(loc)     Y(loc)   Z(loc)    energy-MC     nSimP\n");
    printf("-----------------------------------------------------------------------------------------------");
    printf("-------------------------------\n");
   
    for (int i=0; i<nh; i++) {
      const mu2e::CaloCrystalHit* hit = &(*caloClusterHits.at(i));
      int id = hit->id();
      
      cr = &cal->crystal(id);

      pos = &cr->localPosition();
  
      mu2e::CrystalContentMC contentMC(calo, *CaloHitTruth, *hit);
      double                 simMaxEdep(0);
      int                    simPDGId(9999), simPDGM(9999), simCreationCode(9999), /*simTime(9999),*/ simID(9999);

      //search the simParticle that gave the largest energy contribution
      int  nSimPart(0);
      for (const auto& contentMap : contentMC.simContentMap() )	{	       
	art::Ptr<mu2e::SimParticle> sim  = contentMap.first;
	mu2e::CaloContentSim        data = contentMap.second;
               
	auto parent(sim);
	while ( parent->hasParent()) parent = parent->parent();               
               
	if (data.edep() > simMaxEdep){
	  simMaxEdep = data.edep();
	    
	  simID            = sim->id().asInt();
	  simPDGId         = sim->pdgId();
	  simCreationCode  = parent->creationCode();
	  // simTime          = data.time();
	  simPDGM          = parent->pdgId(); 
	}
	++nSimPart;
      }
      
      printf("%6i   %10.3f %8i %8i %5i %5i    %10.3f   %10.3f %10.3f %10.3f %10.3f %8i\n",
	     id,
	     hit->time(),
	     simCreationCode, simID,simPDGId,simPDGM,
	     // iz,ir,
	     hit->energyDep(),
	     pos->x(),
	     pos->y(),
	     pos->z(),
	     simMaxEdep,
	     nSimPart
	     );
    }
  }
}


//-----------------------------------------------------------------------------
void TAnaDump::printCaloClusterCollection(const char* ModuleLabel, 
					  const char* ProductName,
					  const char* ProcessName,
					  double      Emin,
					  int         HitOpt,
					  const char* MCModuleLabel) {

  printf(">>>> ModuleLabel = %s\n",ModuleLabel);

  //data about hits in the calorimeter crystals

  art::Handle<mu2e::CaloClusterCollection> handle;
  const mu2e::CaloClusterCollection* caloCluster;

  art::Handle<mu2e::CaloHitMCTruthAssns> caloHitTruthHandle;
  const mu2e::CaloHitMCTruthAssns* caloHitTruth(0);
  
  if (ProductName[0] != 0) {
    art::Selector  selector(art::ProductInstanceNameSelector(ProductName) &&
			    art::ProcessNameSelector(ProcessName)         && 
			    art::ModuleLabelSelector(ModuleLabel)            );
    fEvent->get(selector,handle);
    art::Selector  selectorMC(art::ProductInstanceNameSelector(ProductName) &&
			      art::ProcessNameSelector(ProcessName)         && 
			      art::ModuleLabelSelector(MCModuleLabel)            );
   
    fEvent->get(selectorMC,caloHitTruthHandle);
  }
  else {
    art::Selector  selector(art::ProcessNameSelector(ProcessName)         && 
			    art::ModuleLabelSelector(ModuleLabel)            );
    fEvent->get(selector,handle);
    art::Selector  selectorMC(art::ProcessNameSelector(ProcessName)         && 
			      art::ModuleLabelSelector(MCModuleLabel)            );
    fEvent->get(selectorMC,caloHitTruthHandle);
  }
//-----------------------------------------------------------------------------
// make sure collection exists
//-----------------------------------------------------------------------------
  if (! handle.isValid()) {
    printf("TAnaDump::printCaloClusterCollection: no CaloClusterCollection ");
    printf("for module %s and ProductName=%s found, BAIL OUT\n",
	   ModuleLabel,ProductName);
    return;
  }

  caloCluster = handle.product();

  int ncl = caloCluster->size();

  if (caloHitTruthHandle.isValid()) caloHitTruth = caloHitTruthHandle.product();

  const mu2e::CaloCluster* cl;


  int banner_printed = 0;
  for (int i=0; i<ncl; i++) {
    cl = &caloCluster->at(i);
    if (cl->energyDep() < Emin )    continue;
    
    if ( (banner_printed == 0) && (HitOpt==0) ) {
      printCaloCluster(cl, "banner");
      banner_printed = 1;
    }
    if(HitOpt==0) {
      printCaloCluster(cl,"data");
    } else{
      printCaloCluster(cl,"banner+data+hits",caloHitTruth);
    }
  }
 
}


//-----------------------------------------------------------------------------
void TAnaDump::printCaloProtoCluster(const mu2e::CaloProtoCluster* Cluster, const char* Opt) {

  TString opt = Opt;

  int section_id(-1), iz, ir;

  const mu2e::Calorimeter       * cal(NULL);
  const mu2e::Crystal           *cr;
  const CLHEP::Hep3Vector       *pos;

  art::ServiceHandle<mu2e::GeometryService> geom;
  mu2e::GeomHandle  <mu2e::Calorimeter>     cg;

  cal = cg.get();

  if ((opt == "") || (opt == "banner")) {
    printf("-----------------------------------------------------------------------------------------------------\n");
    printf("       Address  SectionID  IsSplit  NC    Time    Energy      \n");
    printf("-----------------------------------------------------------------------------------------------------\n");
  }
 
  const mu2e::CaloProtoCluster::CaloCrystalHitPtrVector caloClusterHits = Cluster->caloCrystalHitsPtrVector();
  int nh = caloClusterHits.size();

  if ((opt == "") || (opt.Index("data") >= 0)) {

    printf("%16p  %3i %5i %5i %10.3f %10.3f\n",
	   Cluster,
	   section_id,
	   nh,
	   Cluster->isSplit(),
	   Cluster->time(),
	   Cluster->energyDep()
	   ); 
  }
  
  if (opt.Index("hits") >= 0) {
//-----------------------------------------------------------------------------
// print individual crystals in local vane coordinate system
//-----------------------------------------------------------------------------
    for (int i=0; i<nh; i++) {
      const mu2e::CaloCrystalHit* hit = &(*caloClusterHits.at(i));
      int id = hit->id();
      
      cr  = &cal->crystal(id);
      pos = &cr->localPosition();
      iz  = -1;
      ir  = -1;
      
      printf("%6i     %10.3f %5i %5i %8.3f %10.3f %10.3f %10.3f %10.3f\n",
	     id,
	     hit->time(),
	     iz,ir,
	     hit->energyDep(),
	     pos->x(),
	     pos->y(),
	     pos->z(),
	     hit->energyDepTot()
	     );
    }
  }
}


//-----------------------------------------------------------------------------
void TAnaDump::printCaloProtoClusterCollection(const char* ModuleLabel, 
					       const char* ProductName,
					       const char* ProcessName) {

  art::Handle<mu2e::CaloProtoClusterCollection> handle;
  const mu2e::CaloProtoClusterCollection       *coll;
  const mu2e::CaloProtoCluster                 *cluster;

  int banner_printed(0), nclusters;

  if (ProductName[0] != 0) {
    art::Selector  selector(art::ProductInstanceNameSelector(ProductName) &&
			    art::ProcessNameSelector(ProcessName)         && 
			    art::ModuleLabelSelector(ModuleLabel)            );
    fEvent->get(selector,handle);
  }
  else {
    art::Selector  selector(art::ProcessNameSelector(ProcessName)         && 
			    art::ModuleLabelSelector(ModuleLabel)            );
    fEvent->get(selector,handle);
  }
//-----------------------------------------------------------------------------
// make sure collection exists
//-----------------------------------------------------------------------------
  if (! handle.isValid()) {
    printf("TAnaDump::printCaloProtoClusterCollection: no CaloProtoClusterCollection ");
    printf("for module %s and ProductName=%s found, BAIL OUT\n",
	   ModuleLabel,ProductName);
    return;
  }

  coll      = handle.product();
  nclusters = coll->size();

  for (int i=0; i<nclusters; i++) {
    cluster = &coll->at(i);
    if (banner_printed == 0) {
      printCaloProtoCluster(cluster, "banner");
      banner_printed = 1;
    }
    printCaloProtoCluster(cluster,"data");
  }
}

//-----------------------------------------------------------------------------
void TAnaDump::printTrackSeed(const mu2e::KalSeed* TrkSeed        ,
			      const char*          Opt            , 
			      const char*          StrawHitCollTag) {
  TString opt = Opt;
  
  if ((opt == "") || (opt == "banner")) {
    printf("-------------------------------------------------------------");
    printf("---------------------------------------------------------------------------------\n");
    printf("  TrkID       Address    N      P      pT      T0       T0err    fmin   fmax");
    printf("      D0       Z0      Phi0    TanDip    radius    caloEnergy     chi2      FitCon\n");
    printf("-------------------------------------------------------------");
    printf("---------------------------------------------------------------------------------\n");
  }
 
  art::InputTag sdmccT(fStrawDigiMCCollTag.Data());
  art::Handle<mu2e::StrawDigiMCCollection> sdmccH;
  const mu2e::StrawDigiMCCollection* sdmcc(nullptr);

  fEvent->getByLabel(sdmccT,sdmccH);
  if (sdmccH.isValid()) sdmcc = sdmccH.product();
  else {
    printf("TAnaDump::printTrackSeed: no StrawDigiMCCollection with tag=%s, BAIL OUT\n", fStrawDigiMCCollTag.Data());
    return;
  }

  if ((opt == "") || (opt.Index("data") >= 0)) {
    int    nhits   = TrkSeed->hits().size();
    
    double t0     = TrkSeed->t0()._t0;
    double t0err  = TrkSeed->t0()._t0err;

    for (const mu2e::KalSegment& kalSeg : TrkSeed->segments()) {

      double d0     = kalSeg.helix().d0();
      double z0     = kalSeg.helix().z0();
      double phi0   = kalSeg.helix().phi0();
      double tandip = kalSeg.helix().tanDip();
      double mm2MeV = 3/10.;
      double mom    = kalSeg.mom();
      double pt     = mom*cos(atan(tandip));
      double radius = pt/mm2MeV;


      const mu2e::CaloCluster*cluster = TrkSeed->caloCluster().get();
      double clusterEnergy(-1);
      if (cluster != 0) clusterEnergy = cluster->energyDep();
      printf("%5i %16p %3i %8.3f %8.5f %7.3f %7.3f %8.3f %8.3f",
	     -1,
	     TrkSeed,
	     nhits,
	     mom, pt, t0, t0err,kalSeg.fmin(),kalSeg.fmax() );

      float chi2    = TrkSeed->chisquared()/double(nhits - 5.);
      float fitCons = TrkSeed->fitConsistency();
      
      printf(" %8.3f %8.3f %8.3f %8.4f %10.4f %10.3f %12.3f %12.3e\n",
	     d0,z0,phi0,tandip,radius,clusterEnergy,chi2,fitCons);
    }
  }

  if ((opt == "") || (opt.Index("hits") >= 0) ){
    int nsh = TrkSeed->hits().size();

    const mu2e::ComboHit* hit(0), *hit_0(0);
    const mu2e::StepPointMC* step(0);

    art::Handle<mu2e::ComboHitCollection> shcolH;
    const mu2e::ComboHitCollection*       shcol(0);
    fEvent->getByLabel<mu2e::ComboHitCollection>(StrawHitCollTag,shcolH);
    if (shcolH.isValid()) shcol = shcolH.product();
    else {
      printf("ERROR in TAnaDump::printTrackSeed: no ComboHitCollection with tag=%s, BAIL OUT\n",StrawHitCollTag);
      return;
    }

    hit_0    = &shcol->at(0);
    int      loc(-1);
    int banner_printed(0);
    for (int i=0; i<nsh; ++i){
      int  hitIndex  = int(TrkSeed->hits().at(i).index());
      hit            = &shcol->at(hitIndex);
      loc            = hit - hit_0;
      if (sdmcc) {
	const mu2e::StrawDigiMC* sdmc = &sdmcc->at(loc);
	if (sdmc->wireEndTime(mu2e::StrawEnd::cal) < sdmc->wireEndTime(mu2e::StrawEnd::hv)) {
	  step = sdmc->stepPointMC(mu2e::StrawEnd::cal).get();
	}
	else {
	  step = sdmc->stepPointMC(mu2e::StrawEnd::hv ).get();
	}
      }

      if (banner_printed == 0){
	// printStrawHit(hit, step, "banner", -1, 0);
	printComboHit(hit, step, "banner", -1, 0);
	banner_printed = 1;
      } else {
	//	printStrawHit(hit, step, "data"  , -1, 0);
	printComboHit(hit, step, "data"  , -1, 0);
      }
     
    }
    
    
  }
}

//-----------------------------------------------------------------------------

void TAnaDump::printTrackSeedCollection(const char* ModuleLabel, 
					const char* ProductName, 
					const char* ProcessName,
					int         hitOpt     ,
					const char* StrawHitCollTag){
  

  art::Handle<mu2e::KalSeedCollection> trkseedsHandle;

  if (ProductName[0] != 0) {
    art::Selector  selector(art::ProductInstanceNameSelector(ProductName) &&
			    art::ProcessNameSelector(ProcessName)         && 
			    art::ModuleLabelSelector(ModuleLabel)            );
    fEvent->get(selector,trkseedsHandle);
  }
  else {
    art::Selector  selector(art::ProcessNameSelector(ProcessName)         && 
			    art::ModuleLabelSelector(ModuleLabel)            );
    fEvent->get(selector,trkseedsHandle);
  }
//-----------------------------------------------------------------------------
// make sure collection exists
//-----------------------------------------------------------------------------
  if (! trkseedsHandle.isValid()) {
    printf("TAnaDump::printTrackSeedCollection: no KalSeedCollection for module %s, BAIL OUT\n",
	   ModuleLabel);
    return;
  }
  
  art::Handle<mu2e::StrawDigiMCCollection> sdmccH;

  fEvent->getByLabel(fStrawDigiMCCollTag.Data(),sdmccH);
  if (sdmccH.isValid()) _mcdigis = sdmccH.product();
  else {
    printf(">>> ERROR in TAnaDump::printTrackSeedCollection: failed to locate StrawDigiMCCollection by %s\n",
	   fStrawDigiMCCollTag.Data());
    _mcdigis = nullptr;
  }

  mu2e::KalSeedCollection*  list_of_trackSeeds;
  list_of_trackSeeds    = (mu2e::KalSeedCollection*) &(*trkseedsHandle);

  int ntrkseeds = list_of_trackSeeds->size();

  const mu2e::KalSeed *trkseed;

  int banner_printed = 0;
  for (int i=0; i<ntrkseeds; i++) {
    trkseed = &list_of_trackSeeds->at(i);
    if (banner_printed == 0) {
      printTrackSeed(trkseed,"banner");
      banner_printed = 1;
    }
    printTrackSeed(trkseed,"data");
    if(hitOpt>0) printTrackSeed(trkseed,"hits",StrawHitCollTag);

  }
}

//-----------------------------------------------------------------------------
void TAnaDump::printHelixSeed(const mu2e::HelixSeed* Helix          , 
			      const char*            StrawHitCollTag, 
			      const char*            StrawDigiMCCollTag,
			      const char*            Opt                ) {
  TString opt(Opt);
  
  if ((opt == "") || (opt == "banner")) {
    printf("------------------------------------------------------------------");
    printf("--------------------------------------------------------------------------------------\n");
    printf("  HelID   Address    N nL nCln     P        pT      T0     T0err  ");
    printf("    D0      FZ0      X0       Y0    Lambda    radius   ECal   chi2XY  chi2ZPhi    flag\n");
    printf("------------------------------------------------------------------");
    printf("--------------------------------------------------------------------------------------\n");
  }
 
  if ((opt == "") || (opt.Index("data") >= 0)) {
    const mu2e::StrawDigiMCCollection* mcdigis(0);
    art::Handle<mu2e::StrawDigiMCCollection> mcdH;
    fEvent->getByLabel(StrawDigiMCCollTag, mcdH);
    if (mcdH.isValid()) mcdigis = mcdH.product();
    else {
      printf("ERROR in TAnaDump::printHelixSeed : no StrawDigiMCCollection tag=%s, BAIL OUT\n",StrawDigiMCCollTag);
      return;
    }

    const mu2e::RobustHelix*robustHel = &Helix->helix();

    if ((opt == "") || (opt.Index("data") >= 0)) {
      int    nhits  = Helix->hits().size();
      int    flag   = *((int*) &Helix->status()); 

      double t0     = Helix->t0()._t0;
      double t0err  = Helix->t0()._t0err;

      double fz0    = robustHel->fz0();
      //    double phi0   = robustHel->phi0();
      double radius = robustHel->radius();
      double d0     = robustHel->rcent() - radius;

      double lambda = robustHel->lambda();
      //      double tandip = lambda/radius;

      double mm2MeV = 3/10.;
      double mom    = robustHel->momentum()*mm2MeV;
      double pt     = radius*mm2MeV;

      double x0     = robustHel->centerx();
      double y0     = robustHel->centery();
    
      int    nLoops(0), nhitsLoopChecked(0);
      evalHelixInfo(Helix, nLoops, nhitsLoopChecked);

      const mu2e::CaloCluster*cluster = Helix->caloCluster().get();
      double clusterEnergy(-1);
      if (cluster != 0) clusterEnergy = cluster->energyDep();
      printf("%5i %12p %3i %2i %4i %8.3f %8.3f %7.3f %7.3f",
	     -1,
	     Helix,
	     nhits,
	     nLoops, 
	     nhitsLoopChecked,
	     mom, pt, t0, t0err );

      float chi2xy   = robustHel->chi2dXY();
      float chi2zphi = robustHel->chi2dZPhi();

      printf(" %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %7.3f %8.3f %8.3f %08x\n",
	     d0,fz0,x0,y0,lambda,radius,clusterEnergy,chi2xy,chi2zphi, flag);
    }

    if ((opt == "") || (opt.Index("hits") >= 0) ){
      int nsh = Helix->hits().size();

      art::Handle<mu2e::ComboHitCollection> shcHandle;
      const mu2e::ComboHitCollection*       shcol;

      const char* ProductName = "";
      const char* ProcessName = "";
//-----------------------------------------------------------------------------
// get ComboHitCollection with single straw hits (makeSH)
//-----------------------------------------------------------------------------
      if (ProductName[0] != 0) {
	art::Selector  selector(art::ProductInstanceNameSelector(ProductName) &&
				art::ProcessNameSelector(ProcessName)         && 
				art::ModuleLabelSelector(StrawHitCollTag)        );
	fEvent->get(selector, shcHandle);
      }
      else {
	art::Selector  selector(art::ProcessNameSelector(ProcessName    ) && 
				art::ModuleLabelSelector(StrawHitCollTag)    );
	fEvent->get(selector, shcHandle);
      }
    
      shcol = shcHandle.product();

      int banner_printed(0);
      for (int i=0; i<nsh; ++i){
	const mu2e::ComboHit*  helHit = &Helix->hits().at(i);
	int  hitIndex                 = helHit->index(0);     // index of the first straw hit
      
	vector<StrawDigiIndex> shids;
	Helix->hits().fillStrawDigiIndices(*(fEvent),i,shids);

	const mu2e::ComboHit* hit     = &shcol->at(hitIndex);
 
	mu2e::StrawDigiMC const&     mcdigi = mcdigis->at(shids[0]);
	art::Ptr<mu2e::StepPointMC>  spmcp;
	mu2e::TrkMCTools::stepPoint(spmcp,mcdigi);

	const mu2e::StepPointMC* Step = spmcp.get();
    
	if (banner_printed == 0) {
	  printHelixHit(helHit, hit, Step, "banner", -1, 0);
	  banner_printed = 1;
	} 
	printHelixHit(helHit, hit, Step, "data", -1, 0);
      }
    }
  }
}

//-----------------------------------------------------------------------------
void TAnaDump::printHelixSeedCollection(const char* HelixSeedCollTag, 
					int         PrintHits       ,
					const char* StrawHitCollTag ,
					const char* StrawDigiCollTag) {
  
  const mu2e::HelixSeedCollection*       list_of_helixSeeds;
  art::Handle<mu2e::HelixSeedCollection> hsH;

  fEvent->getByLabel<mu2e::HelixSeedCollection>(HelixSeedCollTag,hsH);

  if (hsH.isValid()) list_of_helixSeeds = hsH.product();
  else {
    printf("ERROR: cant find HelixSeedCollection tag=%s, avalable collections are:\n",HelixSeedCollTag);

    vector<art::Handle<mu2e::HelixSeedCollection>> list;
    const art::Handle<mu2e::HelixSeedCollection>*  handle;
    const art::Provenance*                         prov;

    fEvent->getManyByType(list);

    for (auto it = list.begin(); it != list.end(); it++) {
      handle = it.operator -> ();
      if (handle->isValid()) {
	prov = handle->provenance();
	
	printf("moduleLabel: %-20s, productInstanceName: %-20s, processName:= %-30s\n" ,
	       prov->moduleLabel().data(),
	       prov->productInstanceName().data(),
	       prov->processName().data()
	       );
      }
    }
    return;
  }

  int nhelices = list_of_helixSeeds->size();

  const mu2e::HelixSeed *helix;

  int banner_printed = 0;

  char popt[20];
  strcpy(popt,"data");
  if (PrintHits > 0) strcat(popt,"+hits");

  for (int i=0; i<nhelices; i++) {
    helix = &list_of_helixSeeds->at(i);
    if (banner_printed == 0) {
      printHelixSeed(helix,StrawHitCollTag,StrawDigiCollTag,"banner");
      banner_printed = 1;
    }

    printHelixSeed(helix,StrawHitCollTag,StrawDigiCollTag,popt);
  }
}

// //-----------------------------------------------------------------------------
// void TAnaDump::plotTrackSeed(int Index, const char* ModuleLabelTrkSeeds, const char* ModuleLabelHitPos){
//-----------------------------------------------------------------------------

  
//   const char* ProductName = "";
//   const char* ProcessName = "";


//   art::Handle<mu2e::StrawHitPositionCollection> spcHandle;
//   const mu2e::StrawHitPositionCollection*       shpcol;

//   if (ProductName[0] != 0) {
//     art::Selector  selector(art::ProductInstanceNameSelector(ProductName) &&
// 			    art::ProcessNameSelector(ProcessName)         && 
// 			    art::ModuleLabelSelector(ModuleLabelHitPos)            );
//     fEvent->get(selector, spcHandle);
//   }
//   else {
//     art::Selector  selector(art::ProcessNameSelector(ProcessName)         && 
// 			    art::ModuleLabelSelector(ModuleLabelHitPos)            );
//     fEvent->get(selector, spcHandle);
//   }

           //                 //now get the strawhitposcollection
//   if (ProductName[0] != 0) {
//     art::Selector  selector(art::ProductInstanceNameSelector(ProductName) &&
// 			    art::ProcessNameSelector(ProcessName)         && 
// 			    art::ModuleLabelSelector(ModuleLabelHitPos)            );
//     fEvent->get(selector, spcHandle);
//   }
//   else {
//     art::Selector  selector(art::ProcessNameSelector(ProcessName)         && 
//// 			    art::ModuleLabelSelector(ModuleLabelHitPos)            );
//     fEvent->get(selector, spcHandle);
//}

//   shpcol = spcHandle.product();
  

//   art::Handle<mu2e::TrackSeedCollection> trkseedsHandle;
//   const mu2e::TrackSeedCollection       *trkseeds(0);


//   if (ProductName[0] != 0) {
//     art::Selector  selector(art::ProductInstanceNameSelector(ProductName) &&
// 			    art::ProcessNameSelector(ProcessName)         && 
// 			    art::ModuleLabelSelector(ModuleLabelTrkSeeds)            );
//     fEvent->get(selector,trkseedsHandle);
//   }
//   else {
//     art::Selector  selector(art::ProcessNameSelector(ProcessName)         && 
// 			    art::ModuleLabelSelector(ModuleLabelTrkSeeds)            );
//     fEvent->get(selector,trkseedsHandle);
//   }

//   if (trkseedsHandle.isValid())  trkseeds = trkseedsHandle.product();

//   if (trkseeds == 0) {
//     printf("no TrackSeedCollection FOUND\n");
//     return;
//   }
//   printf("---------------------------------------------------------------------------------");
//   printf("---------------------------------------------------------------------------------\n");
//   printf("  TrkID       Address    N      P      pT      T0     T0err  ");
//   printf("      D0       Z0      Phi0   TanDip      DfDz     X0       Y0      radius    caloEnergy     chi2XY     chi2ZPhi\n");
//   printf("---------------------------------------------------------------------------------");
//   printf("---------------------------------------------------------------------------------\n");
  
//   const mu2e::TrackSeed* TrkSeed = &trkseeds->at(Index);

//   int    nhits   = TrkSeed->_selectedTrackerHits.size();

  
//   double t0     = TrkSeed->t0();
//   double t0err  = TrkSeed->errt0();
  
//   double d0     = TrkSeed->d0();
//   double z0     = TrkSeed->z0();
//   double phi0   = TrkSeed->phi0();
//   double radius = 1./fabs(TrkSeed->omega());
//   double tandip = TrkSeed->tanDip();
  
//   double mm2MeV = 3/10.;
//   double mom    = radius*mm2MeV/cos(atan(tandip));
//   double pt     = radius*mm2MeV;
  
//   double dfdz   = 1./(radius*tandip);
  
//   const mu2e::CaloCluster*cluster = TrkSeed->_caloCluster.get();
  
//   const mu2e::DiskCalorimeter* cal;
//   art::ServiceHandle<mu2e::GeometryService> geom;
    
//   if (geom->hasElement<mu2e::DiskCalorimeter>() ) {
//     mu2e::GeomHandle<mu2e::DiskCalorimeter> dc;
//     cal = dc.operator->();
//   }
//   else {
//     printf(">>> ERROR: disk calorimeter not found.\n");
//     return;
//   }
 
//   int               idisk   = cluster->sectionId();
//   CLHEP::Hep3Vector cp_mu2e = cal->fromSectionFrame(idisk, cluster->cog3Vector());
//   CLHEP::Hep3Vector cp_st   = cal->toTrackerFrame(cp_mu2e);

//   double x0   = -(radius + d0)*sin(phi0);
//   double y0   =  (radius + d0)*cos(phi0);

//   phi0        = -z0*dfdz + phi0 - M_PI/2.;

//   printf("%5i %16p %3i %8.3f %8.5f %7.3f %7.3f",
// 	 -1,
// 	 TrkSeed,
// 	 nhits,
// 	 mom, pt, t0, t0err );
  
//   printf(" %8.3f %8.3f %8.3f %8.4f %10.4f  %10.4f %10.4f %10.4f %8.3f %10.3f %8.3f\n",
// 	 d0, z0, phi0, tandip, dfdz, x0, y0, radius, cluster->energyDep(), TrkSeed->chi2XY(), TrkSeed->chi2ZPhi());


//   //fill Z and Phi of the straw hits
//   const mu2e::StrawHitPosition   *pos(0);
//   const mu2e::HelixVal           *helix = &TrkSeed->_fullTrkSeed;

  //fill Z and Phi of the straw hits
  //  const mu2e::HelixVal           *helix = &TrkSeed->_helix;
  
//   int    hitIndex(0);
//   int    np = nhits + 1;
//   double x[np], y[np], z[np];
  
//   for (int i=0; i<nhits; ++i){
//     hitIndex  = helix->_selectedTrackerHitsIdx.at(i)._index;
    
//     pos       = &shpcol->at(hitIndex);
//     x[i]      = pos->pos().x();  
//     y[i]      = pos->pos().y();  
//     z[i]      = pos->pos().z();  
//   } 

//   //fill calo-cluster position
//   x[nhits]    = cluster->cog3Vector().x();
//   y[nhits]    = cluster->cog3Vector().y();
//   z[nhits]    = cp_st.z();
  

//   // TGraph* gr_xy = new TGraph(np,x,y);
//   // TGraph* gr_yz = new TGraph(np,z,y);

//   if (c_plot_hits_xy != 0) delete c_plot_hits_xy;
//   c_plot_hits_xy = new TCanvas("htrakseedview","htrakseedview",1200,600);

//   c_plot_hits_xy->Divide(2,1);

//   c_plot_hits_xy->cd(1);

//   if (h2_xy != 0) delete h2_xy;
//   h2_xy = new TH2F("h2_xy","XY View",140,-700,700,140,-700,700);
//   h2_xy->SetStats(0);
//   h2_xy->Draw();

//   TMarker* m;

//   int color;

//   for (int i=0; i<np; i++) {
//     m = new TMarker(x[i],y[i],2);
//     if (i < np - 1) color = kBlack;
//     else            color = kRed;
//     m->SetMarkerSize(1.5);
//     m->SetMarkerColor(color);
//     m->Draw();
//   }

//   if (e != 0) delete e;
//   e = new TEllipse(x0, y0, radius);

//   e->SetFillStyle(0);
//   e->Draw();


//   c_plot_hits_xy->cd(2);

//   if (h2_yz != 0) delete h2_yz;
//   h2_yz = new TH2F("h2_yz","YZ VIEW",2500,-2500, 2500, 140,-700,700);
//   h2_yz->SetStats(0);
//   h2_yz->Draw();


//   for (int i=0; i<np; i++) {
//     m = new TMarker(z[i],y[i],2);
//     if (i < np - 1) color = kBlack;
//     else            color = kRed;
//     m->SetMarkerColor(color);
//     m->SetMarkerSize(1.5);
//     m->Draw();
//   }

//   if (yf != 0) delete yf;
//   yf = new TF1("yf","[0]*TMath::Sin([1]+x*[2]) + [3]",-1600., 3000.);

//   //correct phi0...

 
//   yf->SetParameters(radius, phi0, dfdz, y0);
//   yf->Draw("same");
// }
//-----------------------------------------------------------------------------
void TAnaDump::printTimeCluster(const mu2e::TimeCluster* TPeak, const char* Opt, const mu2e::ComboHitCollection* ChColl, const char*StrawDigiMCModuleLabel) {

  const mu2e::ComboHit*      hit;
  int                        flags;
  //  const mu2e::StepPointMC*   step(NULL);
  // art::Handle<mu2e::PtrStepPointMCVectorCollection> mcptrHandleStraw;
  // fEvent->getByLabel("makeSD","",mcptrHandleStraw);
  // mu2e::PtrStepPointMCVectorCollection const* hits_mcptrStraw = mcptrHandleStraw.product();
  const mu2e::StrawDigiMCCollection* mcdigis(0);
  art::Handle<mu2e::StrawDigiMCCollection> mcdH;
  fEvent->getByLabel(StrawDigiMCModuleLabel, mcdH);
  if (mcdH.isValid()) mcdigis = mcdH.product();
  
  TString opt = Opt;
  opt.ToLower();

  if ((opt == "") || (opt.Index("banner") >= 0)) {
    printf("-------------------------------------------------------------------\n");
    printf("    Energy       X          Y          Z        T0       NCH   NSH \n");
    printf("-------------------------------------------------------------------\n");
  }
  double caloClusterEnergy(-1);
  if (TPeak->caloCluster().get()!=0)  caloClusterEnergy = TPeak->caloCluster()->energyDep();

  if ((opt == "") || (opt.Index("data") >= 0)) {

    int  nhits, nsh(0), loc;
    nhits = TPeak->nhits();

    for (int i=0; i<nhits; i++) {
      loc  = TPeak->hits().at(i);
      hit  = &(ChColl->at(loc));
      nsh += hit->nStrawHits();
    }

    printf("%10.3f %10.3f %10.3f %10.3f %10.3f %5i %5i\n",
	   caloClusterEnergy, 
	   TPeak->position().x(),
	   TPeak->position().y(),
	   TPeak->position().z(),
	   TPeak->t0().t0(),
	   nhits           ,
	   nsh             );

    if (opt.Index("hits") >= 0) {
//-----------------------------------------------------------------------------
// print straw hits in the list
//-----------------------------------------------------------------------------
      printComboHit(0,0,"banner",0,0);

      int  nhits, loc;
      nhits = TPeak->nhits();

      for (int i=0; i<nhits; i++) {
	loc         = int(TPeak->hits().at(i));
	hit         = &(ChColl->at(loc));
	vector<StrawDigiIndex> shids;
	ChColl->fillStrawDigiIndices(*(fEvent),loc,shids);
	mu2e::StrawDigiMC const&     mcdigi = mcdigis->at(shids[0]);
	art::Ptr<mu2e::StepPointMC>  spmcp;
	mu2e::TrkMCTools::stepPoint(spmcp,mcdigi);
	const mu2e::StepPointMC* step = spmcp.get();
	flags = 0; // *((int*) &TPeak->_shfcol->at(loc));

	printComboHit(hit,step,"data",i,flags);
      }
    }
  }
}





//-----------------------------------------------------------------------------
void TAnaDump::printTimeClusterCollection(const char* ModuleLabel         , 
					  const char* ComboHitModuleLabel,
					  const char* ProductName         , 
					  const char* ProcessName         ,
					  int         hitOpt              ,
					  const char* StrawDigiMCModuleLabel) {

  art::Handle<mu2e::TimeClusterCollection>  handle;
  const mu2e::TimeClusterCollection*        coll(0);
  const mu2e::TimeCluster*                  tp(0);

  art::Handle<mu2e::ComboHitCollection>     chHandle;
  const mu2e::ComboHitCollection*           chColl(0);
 

  if (ProductName[0] != 0) {
    art::Selector  selector(art::ProductInstanceNameSelector(ProductName) &&
			    art::ProcessNameSelector        (ProcessName) && 
			    art::ModuleLabelSelector        (ModuleLabel)    );
    fEvent->get(selector, handle);
  }
  else {
    art::Selector  selector(art::ProcessNameSelector(ProcessName)         && 
			    art::ModuleLabelSelector(ModuleLabel)            );
    fEvent->get(selector, handle);
  }

  if (handle.isValid()) coll = (mu2e::TimeClusterCollection* )handle.product();
  else {
    printf(">>> ERROR in TAnaDump::printSimParticleCollection: failed to locate collection");
    printf(". BAIL OUT. \n");
    return;
  }

  art::Selector  selectorCombo(art::ProcessNameSelector(ProcessName)         && 
			       art::ModuleLabelSelector(ComboHitModuleLabel)  );
  fEvent->get(selectorCombo, chHandle);
  if (chHandle.isValid()) chColl = (mu2e::ComboHitCollection* )chHandle.product();

  int banner_printed(0);

  for ( mu2e::TimeClusterCollection::const_iterator j=coll->begin(); j != coll->end(); ++j ){
    tp = &(*j);

    if (banner_printed == 0) {
      printTimeCluster(tp,"banner");
      banner_printed = 1;
    }

    if (hitOpt > 0) printTimeCluster(tp,"data+hits",chColl,StrawDigiMCModuleLabel);
    else            printTimeCluster(tp,"data",chColl,StrawDigiMCModuleLabel);

  }
}

//-----------------------------------------------------------------------------
void TAnaDump::printEventHeader() {

  printf(" Run / Subrun / Event : %10i / %10i / %10i\n",
	 fEvent->run(),
	 fEvent->subRun(),
	 fEvent->event());
}

//-----------------------------------------------------------------------------
void TAnaDump::printTrkCaloHit(const KalRep* Krep, mu2e::TrkCaloHit* CaloHit){
  double    len  = CaloHit->fltLen();
  HepPoint  plen = Krep->position(len);
  
  printf("%3i %5i 0x%08x %1i %9.3f %8.3f %8.3f %9.3f %8.3f %7.3f",
	 -1,//++i,
	 0, //straw->index().asInt(), 
	 CaloHit->hitFlag(),
	 //	     hit->isUsable(),
	 CaloHit->isActive(),
	 len,
	 //	     hit->hitRms(),
	 plen.x(),plen.y(),plen.z(),
	 CaloHit->time(), -1.//sh->dt()
	 );

  printf(" %2i %2i %2i %2i",
	 -1,//straw->id().getPlane(),
	 -1,//straw->id().getPanel(),
	 -1,//straw->id().getLayer(),
	 -1//straw->id().getStraw()
	 );

  printf(" %8.3f",CaloHit->hitT0().t0());
  
  double res, sigres;
  CaloHit->resid(res, sigres, true);

  CLHEP::Hep3Vector  pos;
  CaloHit->hitPosition(pos);

  printf("%8.3f %8.3f %9.3f %7.3f %7.3f",
	 pos.x(),
	 pos.y(),
	 pos.z(),
	 res,
	 sigres
	 );
      
  printf("   %6.3f", -1.);//CaloHit->driftRadius());
  
	  

  printf("  %7.3f",-1.);

  double exterr = CaloHit->temperature();//*CaloHit->driftVelocity();

  printf(" %6.3f %6.3f %6.3f %6.3f %6.3f",		 
	 -1.,//CaloHit->totalErr(),
	 CaloHit->hitErr(),
	 CaloHit->hitT0().t0Err(),
	 -1.,//,CaloHit->penaltyErr(),
	 exterr
	 );
  //-----------------------------------------------------------------------------
  // test: calculated residual in fTmp[0]
  //-----------------------------------------------------------------------------
  //       Test_000(Krep,hit);
  //       printf(" %7.3f",fTmp[0]);

  printf("\n");

}

//-----------------------------------------------------------------------------
// ""       : banner+track parameters (default)
// "banner" : banner only
// "data"   : track parameters only
// "hits"   : hits
//-----------------------------------------------------------------------------
void TAnaDump::printKalRep(const KalRep* Krep, const char* Opt, const char* Prefix) {

  //  TString opt = Opt;

  _printUtils->printTrack(fEvent,Krep,Opt,Prefix);
}

//-----------------------------------------------------------------------------
void TAnaDump::printKalRepCollection(const char* ModuleLabel, 
				     const char* ProductName,
				     const char* ProcessName,
				     int         hitOpt     ) {

  art::Handle<mu2e::KalRepPtrCollection> krepsHandle;

  if (ProductName[0] != 0) {
    art::Selector  selector(art::ProductInstanceNameSelector(ProductName) &&
			    art::ProcessNameSelector(ProcessName)         && 
			    art::ModuleLabelSelector(ModuleLabel)            );
    fEvent->get(selector,krepsHandle);
  }
  else {
    art::Selector  selector(art::ProcessNameSelector(ProcessName)         && 
			    art::ModuleLabelSelector(ModuleLabel)            );
    fEvent->get(selector,krepsHandle);
  }
//-----------------------------------------------------------------------------
// make sure collection exists
//-----------------------------------------------------------------------------
  if (! krepsHandle.isValid()) {
    printf("TAnaDump::printKalRepCollection: no KalRepPtrCollection for module %s, BAIL OUT\n",
	   ModuleLabel);
    return;
  }

  _mcdigis = &fEvent->getByLabel<mu2e::StrawDigiMCCollection>(fStrawDigiMCCollTag.Data());
  if (_mcdigis == nullptr) {
    printf(">>> ERROR in TAnaDump::printKalRepCollection: failed to locate StepPointMCCollection:: by %s\n",fStrawDigiMCCollTag.Data());
  }

  int ntrk = krepsHandle->size();

  const KalRep *trk;

  int banner_printed = 0;
  for (int i=0; i<ntrk; i++) {
    art::Ptr<KalRep> kptr = krepsHandle->at(i);
    //    fEvent->get(kptr.id(), krepsHandle);
    fhicl::ParameterSet const& pset = krepsHandle.provenance()->parameterSet();
    string module_type = pset.get<string>("module_type");
 
    trk = kptr.get();
    if (banner_printed == 0) {
      printKalRep(trk,"banner",module_type.data());
      banner_printed = 1;
    }
    printKalRep(trk,"data",module_type.data());
    if (hitOpt > 0) printKalRep(trk,"hits");
  }
 
}


//-----------------------------------------------------------------------------
void TAnaDump::printGenParticle(const mu2e::GenParticle* P, const char* Opt) {

  TString opt = Opt;
  
  if ((opt == "") || (opt == "banner")) {
    printf("------------------------------------------------------------------------------------\n");
    printf("Index                 generator     PDG      Time      Momentum       Pt       CosTh\n");
    printf("------------------------------------------------------------------------------------\n");
  }
  
  if ((opt == "") || (opt == "data")) {
    int    gen_code   = P->generatorId().id();
    string gen_name   = P->generatorId().name();
    int    pdg_code   = P->pdgId();
    double time       = P->time();
    
    double mom   = P->momentum().vect().mag();
    double pt    = P->momentum().vect().perp();
    double costh = P->momentum().vect().cosTheta();
    
    printf("%5i %2i:%-26s %3i %10.3f %10.3f %10.3f %10.3f\n",
	   -1,gen_code,gen_name.data(),pdg_code,time,mom,pt,costh);
  }
}

//-----------------------------------------------------------------------------
// there could be multiple collections in the event
//-----------------------------------------------------------------------------
void TAnaDump::printGenParticleCollections() {
  
  vector<art::Handle<mu2e::GenParticleCollection>> list_of_gp;

  const mu2e::GenParticleCollection*        coll(0);
  const mu2e::GenParticle*        genp(0);

  const art::Provenance* prov;

  //  art::Selector  selector(art::ProductInstanceNameSelector("mu2e::GenParticleCollection"));
  art::Selector  selector(art::ProductInstanceNameSelector(""));

  fEvent->getMany(selector, list_of_gp);

  const art::Handle<mu2e::GenParticleCollection>* handle;

  int banner_printed;
  for (vector<art::Handle<mu2e::GenParticleCollection>> ::const_iterator it = list_of_gp.begin();
       it != list_of_gp.end(); it++) {
    handle = it.operator -> ();

    if (handle->isValid()) {
      coll = handle->product();
      prov = handle->provenance();

      printf("moduleLabel = %-20s, producedClassname = %-30s, productInstanceName = %-20s\n",
	     prov->moduleLabel().data(),
	     prov->producedClassName().data(),
	     prov->productInstanceName().data());

      banner_printed = 0;
      for (vector<mu2e::GenParticle>::const_iterator ip = coll->begin();
	   ip != coll->end(); ip++) {
	genp = ip.operator -> ();
	if (banner_printed == 0) {
	  printGenParticle(genp,"banner");
	  banner_printed = 1;
	}
	printGenParticle(genp,"data");
      }

      
    }
    else {
      printf(">>> ERROR in TAnaDump::printStepPointMCCollection: failed to locate collection");
      printf(". BAIL OUT. \n");
      return;
    }
  }
}


// //-----------------------------------------------------------------------------
//   void TAnaDump::printCaloHit(const CaloHit* Hit, const char* Opt) {
//     //    int row, col;
//     TString opt = Opt;

//     if ((opt == "") || (opt == "banner")) {
//       printf("--------------------------------------\n");
//       printf("RID      Time   Energy                \n");
//       printf("--------------------------------------\n");
//     }
    
//     if ((opt == "") || (opt == "data")) {
//       printf("%7i  %10.3f %10.3f \n",
// 	     Hit->id(),
// 	     Hit->time(),
// 	     Hit->energyDep()); 
//     }
//   }


//-----------------------------------------------------------------------------
void TAnaDump::printCaloHits(const char* ModuleLabel, 
			     const char* ProductName, 
			     const char* ProcessName) {

  printf(">>>> ModuleLabel = %s\n",ModuleLabel);

  //data about hits in the calorimeter crystals

  art::Selector  selector(art::ProductInstanceNameSelector(ProductName) &&
			  art::ProcessNameSelector(ProcessName)         && 
			  art::ModuleLabelSelector(ModuleLabel)            );

  art::Handle<mu2e::CaloHitCollection> caloHitsHandle;

  fEvent->get(selector,caloHitsHandle);

  const mu2e::CaloHitCollection* caloHits;

  caloHits = caloHitsHandle.operator ->();

  int nhits = caloHits->size();

  const mu2e::CaloHit* hit;

  printf("--------------------------------------\n");
  printf("RID      Time   Energy                \n");
  printf("--------------------------------------\n");

  for (int ic=0; ic<nhits; ic++) {
    hit  = &caloHits->at(ic);
    printf("%7i  %10.3f %10.3f \n",
	   hit->id(),
	   hit->time(),
	   hit->energyDep()); 
  }
}


//-----------------------------------------------------------------------------
void TAnaDump::printDiskCalorimeter() {
  const mu2e::DiskCalorimeter* cal;
  const mu2e::Disk* disk;
  
  art::ServiceHandle<mu2e::GeometryService> geom;
    
  if (geom->hasElement<mu2e::DiskCalorimeter>() ) {
    mu2e::GeomHandle<mu2e::DiskCalorimeter> dc;
    cal = dc.operator->();
  }
  else {
    printf(">>> ERROR: disk calorimeter not found.\n");
    return;
  }

  int nd = cal->nDisk();
  printf(" ndisks = %i\n", nd);
  printf(" crystal size  : %10.3f\n", cal->caloInfo().getDouble("crystalXYLength"));
  printf(" crystal length: %10.3f\n", cal->caloInfo().getDouble("crystalZLength"));

  for (int i=0; i<nd; i++) {
    disk = &cal->disk(i);
    printf(" ---- disk # %i\n",i);
    printf(" Rin  : %10.3f  Rout : %10.3f\n", disk->innerRadius(),disk->outerRadius());
    printf(" X : %12.3f Y : %12.3f Z : %12.3f\n",
	   disk->geomInfo().origin().x(),
	   disk->geomInfo().origin().y(),
	   disk->geomInfo().origin().z());
    // printf(" Xsize : %10.3f Ysize : %10.3f Zsize : %10.3f\n", 
    // 	   disk->size().x(),
    // 	   disk->size().y(),
    // 	   disk->size().z()
    // 	   );
  }
}


//-----------------------------------------------------------------------------
void TAnaDump::printCaloCrystalHits(const char* ModuleLabel, 
				    const char* ProductName,
				    const char* ProcessName) {

  art::Selector  selector(art::ProductInstanceNameSelector(ProductName) &&
			  art::ProcessNameSelector(ProcessName)         && 
			  art::ModuleLabelSelector(ModuleLabel)            );

  art::Handle<mu2e::CaloCrystalHitCollection> caloCrystalHitsHandle;

  fEvent->get(selector,caloCrystalHitsHandle);

  const mu2e::CaloCrystalHitCollection* caloCrystalHits;

  caloCrystalHits = caloCrystalHitsHandle.operator->();

  int nhits = caloCrystalHits->size();

  const mu2e::CaloCrystalHit* hit;

  printf("----------------------------------------------------------------\n");
  printf("CrystalID      Time   Energy    EnergyTot  NRoids               \n");
  printf("----------------------------------------------------------------\n");

  for (int ic=0; ic<nhits; ic++) {
    hit  = &caloCrystalHits->at(ic);

    printf("%7i  %10.3f %10.3f %10.3f %5i\n",
	   hit->id(),
	   hit->time(),
	   hit->energyDep(),
	   hit->energyDepTot(),
	   hit->nROId());
  }
}

//-----------------------------------------------------------------------------
void TAnaDump::printCaloDigiCollection(const char* ModuleLabel, 
				       const char* ProductName,
				       const char* ProcessName) {

  art::Selector  selector(art::ProductInstanceNameSelector(ProductName) &&
			  art::ProcessNameSelector(ProcessName)         && 
			  art::ModuleLabelSelector(ModuleLabel)            );

  art::Handle<mu2e::CaloDigiCollection> calodigisHandle;

  fEvent->get(selector,calodigisHandle);

  const mu2e::CaloDigiCollection* calodigis;

  calodigis = calodigisHandle.operator->();

  int nhits = calodigis->size();

  const mu2e::CaloDigi* hit;

  printf("----------------------------------------------------------------\n");
  printf("ReadoutID      Time      NSamples               \n");
  printf("----------------------------------------------------------------\n");

  for (int ic=0; ic<nhits; ic++) {
    hit  = &calodigis->at(ic);
    int pulse_size =  hit->waveform().size();

    printf("%7i  %5i %5i\n",
	   hit->roId(),
	   hit->t0(),
	   pulse_size);
  }
}



//-----------------------------------------------------------------------------
void TAnaDump::printCaloRecoDigiCollection(const char* ModuleLabel, 
				       const char* ProductName,
				       const char* ProcessName) {

  art::Selector  selector(art::ProductInstanceNameSelector(ProductName) &&
			  art::ProcessNameSelector(ProcessName)         && 
			  art::ModuleLabelSelector(ModuleLabel)            );

  art::Handle<mu2e::CaloRecoDigiCollection> recocalodigisHandle;

  fEvent->get(selector,recocalodigisHandle);

  const mu2e::CaloRecoDigiCollection* recocalodigis;

  recocalodigis = recocalodigisHandle.operator->();

  int nhits = recocalodigis->size();

  const mu2e::CaloRecoDigi* hit;

  printf("-----------------------------------------------------------------------------------\n");
  printf("ReadoutID      Time      Time-Chi2     Energy     Amplitude      PSD               \n");
  printf("-----------------------------------------------------------------------------------\n");

  for (int ic=0; ic<nhits; ic++) {
    hit  = &recocalodigis->at(ic);

    printf("%7i  %10.3f   %10.3f   %10.3f   %10.3f   %10.3f\n",
	   hit->ROid(),
	   hit->time(),
	   hit->chi2(), 
	   hit->energyDep(),
	   -1., //hit->amplitude(),
	   -1.);//hit->psd());
  }
}

//------------------------------------------------------------------
// void TAnaDump::printTrackClusterLink(const char* ModuleLabel, 
// 			     const char* ProductName,
// 			     const char* ProcessName) {

//   printf(">>>> ModuleLabel = %s\n",ModuleLabel);

//   //data about hits in the calorimeter crystals

//   art::Handle<mu2e::TrackClusterLink> trkCluLinkHandle;
//   const mu2e::TrackClusterLink* trkCluLink;

//   if (ProductName[0] != 0) {
//     art::Selector  selector(art::ProductInstanceNameSelector(ProductName) &&
// 			    art::ProcessNameSelector(ProcessName)         && 
// 			    art::ModuleLabelSelector(ModuleLabel)            );
//     fEvent->get(selector, trkCluLinkHandle);
//   }
//   else {
//     art::Selector  selector(art::ProcessNameSelector(ProcessName)         && 
// 			    art::ModuleLabelSelector(ModuleLabel)            );
//     fEvent->get(selector, trkCluLinkHandle);
//   }
  
//   trkCluLink = trkCluLinkHandle.operator ->();

//   int nhits = trkCluLink->size();

//   const mu2e::CaloCluster* clu;
//   //const KalRep *trk;

//   int banner_printed = 0;
//   const mu2e::TrkToCaloExtrapol* extrk;
  


//   for (int i=0; i<nhits; i++) {
//     clu = &(*(trkCluLink->at(i).second.get()) );
//     extrk = &(*(trkCluLink->at(i).first));
//     if (banner_printed == 0) {
//       printCaloCluster(clu, "banner");
//       //banner_printed = 1;
//     }
//     printCaloCluster(clu,"data");
    
//     if (banner_printed == 0) {
//       printTrkToCaloExtrapol(extrk,"banner");
//       //printKalRep(trk,"banner");
//       banner_printed = 1;
//     }
//     // printKalRep(trk,"data");
//     printTrkToCaloExtrapol(extrk,"data");
//   }
  
// }

////////////////////////////////////////////////////////////////////////////////

void TAnaDump::printTrkToCaloExtrapol(const mu2e::TrkToCaloExtrapol* trkToCalo,
				      const char* Opt) {
 TString opt = Opt;

  if ((opt == "") || (opt == "banner")) {
    printf("-------------------------------------------------------------------------------------------------------\n");
    printf("sectionId      Time     ExtPath     Ds       FitCon      t0          X           Y        Z          Mom  \n");
    printf("-------------------------------------------------------------------------------------------------------\n");
  }
  
  if ((opt == "") || (opt.Index("data") >= 0)) {

    double ds = trkToCalo->pathLengthExit()-trkToCalo->pathLengthEntrance();
  
    printf("%6i %10.3f %10.3f %8.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f \n",
	   trkToCalo->diskId(),
	   trkToCalo->time(),
	   trkToCalo->pathLengthEntrance(),
	   ds,
	   trkToCalo->fitConsistency(),
	   trkToCalo->t0(),
	   trkToCalo->entrancePosition().x(),
	   trkToCalo->entrancePosition().y(),
	   trkToCalo->entrancePosition().z(),
	   trkToCalo->momentum().mag() );
  }
  
}

////////////////////////////////////////////////////////////////////////////////

void TAnaDump::printTrkToCaloExtrapolCollection(const char* ModuleLabel, 
						const char* ProductName,
						const char* ProcessName) {

  printf(">>>> ModuleLabel = %s\n",ModuleLabel);

  //data about hits in the calorimeter crystals

  art::Handle<mu2e::TrkToCaloExtrapolCollection> trkToCaloHandle;
  const mu2e::TrkToCaloExtrapolCollection* trkToCalo;

  if (ProductName[0] != 0) {
    art::Selector  selector(art::ProductInstanceNameSelector(ProductName) &&
			    art::ProcessNameSelector(ProcessName)         && 
			    art::ModuleLabelSelector(ModuleLabel)            );
    fEvent->get(selector, trkToCaloHandle);
  }
  else {
    art::Selector  selector(art::ProcessNameSelector(ProcessName)         && 
			    art::ModuleLabelSelector(ModuleLabel)            );
    fEvent->get(selector, trkToCaloHandle);
  }

  trkToCalo = trkToCaloHandle.operator ->();

  int nhits = trkToCalo->size();

  const mu2e::TrkToCaloExtrapol* hit;
  
  int banner_printed = 0;
  for (int i=0; i<nhits; i++) {
    hit = &trkToCalo->at(i);
    if (banner_printed == 0) {
      printTrkToCaloExtrapol(hit, "banner");
      banner_printed = 1;
    }
    printTrkToCaloExtrapol(hit,"data");
  }
  
}

////////////////////////////////////////////////////////////////////////////////


//-----------------------------------------------------------------------------
void TAnaDump::printStrawHit(const mu2e::StrawHit* Hit, const mu2e::StepPointMC* Step, const char* Opt, int IHit, int Flags) {
    TString opt = Opt;
    opt.ToLower();

    if ((opt == "") || (opt.Index("banner") >= 0)) {
      printf("-----------------------------------------------------------------------------------");
      printf("------------------------------------------------------------\n");
      printf("    I   SID    Flags  Plane   Panel  Layer   Straw     Time          dt       eDep ");
      printf("           PDG       PDG(M)   Generator      SimpID      p  \n");
      printf("-----------------------------------------------------------------------------------");
      printf("------------------------------------------------------------\n");
    }

    if (opt == "banner") return;

    mu2e::GeomHandle<mu2e::Tracker> ttH;
    const mu2e::Tracker* tracker = ttH.get();

    const mu2e::Straw* straw = &tracker->getStraw(Hit->strawId());

    const mu2e::SimParticle* sim(0);
    
    int      pdg_id(-1), mother_pdg_id(-1), generator_id(-1), simp_id(-1);
    double   mc_mom(-1.);

    mu2e::GenId gen_id;

    if (Step) {
      art::Ptr<mu2e::SimParticle> const& simptr = Step->simParticle(); 
      art::Ptr<mu2e::SimParticle> mother = simptr;

      while(mother->hasParent()) mother = mother->parent();

      sim = mother.operator ->();

      pdg_id        = simptr->pdgId();
      mother_pdg_id = sim->pdgId();

      if (simptr->fromGenerator()) generator_id = simptr->genParticle()->generatorId().id();
      else                         generator_id = -1;

      simp_id       = simptr->id().asInt();
      mc_mom        = Step->momentum().mag();
    }
    
    if ((opt == "") || (opt == "data")) {
      if (IHit  >= 0) printf("%5i " ,IHit);
      else            printf("    ");

      printf("%5i",Hit->strawId().asUint16());

      if (Flags >= 0) printf(" %08x",Flags);
      else            printf("        ");
      printf("  %5i  %5i   %5i   %5i   %8.3f   %8.3f   %9.6f   %10i   %10i  %10i  %10i %8.3f\n",
	     straw->id().getPlane(),
	     straw->id().getPanel(),
	     straw->id().getLayer(),
	     straw->id().getStraw(),
	     Hit->time(),
	     Hit->dt(),
	     Hit->energyDep(),
	     pdg_id,
	     mother_pdg_id,
	     generator_id,
	     simp_id,
	     mc_mom);
    }
  }


//-----------------------------------------------------------------------------
void TAnaDump::printComboHit(const mu2e::ComboHit* Hit, const mu2e::StepPointMC* Step, const char* Opt, int IHit, int Flags) {
    TString opt = Opt;
    opt.ToLower();

    if ((opt == "") || (opt.Index("banner") >= 0)) {
      printf("--------------------------------------------------------------------------------------------");
      printf("------------------------------------------------------------------------------------------------\n");
      printf("    I nsh   SID    Flags  Pln:Pnl:Lay:Str      X        Y        Z        Time         eDep ");
      printf("   End  DrTime    TRes   WDist     WRes         PDG       PDG(M)   Generator         ID       p  \n");
      printf("--------------------------------------------------------------------------------------------");
      printf("------------------------------------------------------------------------------------------------\n");
    }

    if (opt == "banner") return;

    const mu2e::SimParticle * sim (0);
    
    int      pdg_id(-1), mother_pdg_id(-1), generator_id(-1), sim_id(-1);
    double   mc_mom(-1.);

    mu2e::GenId gen_id;

    if (Step) {
      art::Ptr<mu2e::SimParticle> const& simptr = Step->simParticle(); 
      art::Ptr<mu2e::SimParticle> mother        = simptr;

      while(mother->hasParent()) mother = mother->parent();

      sim           = mother.operator ->();

      pdg_id        = simptr->pdgId();
      mother_pdg_id = sim->pdgId();

      if (simptr->fromGenerator()) generator_id = simptr->genParticle()->generatorId().id();
      else                         generator_id = -1;

      sim_id        = simptr->id().asInt();
      mc_mom        = Step->momentum().mag();
    }
    
    if ((opt == "") || (opt == "data")) {
      if (IHit  >= 0) printf("%5i " ,IHit);
      else            printf("    ");

      printf("%3i ",Hit->nStrawHits());

      printf("%5u",Hit->strawId().asUint16());

      if (Flags >= 0) printf(" %08x",Flags);
      else            printf("        ");

      printf("  %3i %3i %3i %3i  %8.3f %8.3f %9.3f  %8.3f    %9.6f  %3i %7.2f %7.2f %8.3f %8.3f %10i   %10i  %10i  %10i %8.3f\n",
	     Hit->strawId().plane(),
	     Hit->strawId().panel(),
	     Hit->strawId().layer(),
	     Hit->strawId().straw(),
	     Hit->pos().x(),Hit->pos().y(),Hit->pos().z(),
	     Hit->time(),
	     Hit->energyDep(),

	     (int) Hit->driftEnd(),
	     Hit->driftTime(),
	     Hit->transRes(),
	     Hit->wireDist(),
	     Hit->wireRes(),

	     pdg_id,
	     mother_pdg_id,
	     generator_id,
	     sim_id,
	     mc_mom);
    }
  }


//-----------------------------------------------------------------------------
void TAnaDump::printHelixHit(const mu2e::HelixHit*    HelHit, const mu2e::ComboHit* Hit, 
			     const mu2e::StepPointMC* Step, const char* Opt, int IHit, int Flags) {
    TString opt = Opt;
    opt.ToLower();

    if ((opt == "") || (opt.Index("banner") >= 0)) {
      printf("----------------------------------------------------------------------------------------");
      printf("------------------------------------------------------------------\n");
      printf("   I NSH  SHID   Flags  Pl Pn L  S      x         y         z       phi    Time    eDep ");
      printf("       PDG     PDG(M)  GenID      SimID      p       pT         pZ\n");
      printf("----------------------------------------------------------------------------------------");
      printf("------------------------------------------------------------------\n");
    }

    if (opt == "banner") return;

    mu2e::GeomHandle<mu2e::Tracker> ttHandle;
    const mu2e::Tracker* tracker = ttHandle.get();

    const mu2e::Straw* straw;

    straw = &tracker->getStraw(Hit->strawId());//ndex());

    const mu2e::SimParticle * sim (0);
    
    int      pdg_id(-1), mother_pdg_id(-1), generator_id(-1), sim_id(-1);
    double   mc_mom(-1.), mc_pT(-1.), mc_pZ(0.);
    double   x(0), y(0), z(0), phi(0);

    XYZVec shPos = HelHit->pos();
    x      = shPos.x();
    y      = shPos.y();
    z      = shPos.z();    
    phi    = HelHit->helixPhi();
    
    mu2e::GenId gen_id;

    if (Step) {
      art::Ptr<mu2e::SimParticle> const& simptr = Step->simParticle(); 
      art::Ptr<mu2e::SimParticle> mother = simptr;

      while(mother->hasParent()) mother = mother->parent();

      sim = mother.operator ->();

      pdg_id        = simptr->pdgId();
      mother_pdg_id = sim->pdgId();

      if (simptr->fromGenerator()) generator_id = simptr->genParticle()->generatorId().id();
      else                         generator_id = -1;

      sim_id        = simptr->id().asInt();
      mc_mom        = Step->momentum().mag();
      mc_pT         = sqrt(Step->momentum().x()*Step->momentum().x() + Step->momentum().y()*Step->momentum().y());
      mc_pZ         = Step->momentum().z();
    }
    
    if ((opt == "") || (opt == "data")) {
      if (IHit  >= 0) printf("%5i " ,IHit);
      else            printf("    ");

      printf(" %3i ",HelHit->nStrawHits());

      printf("%5i",Hit->strawId().asUint16());

      if (Flags >= 0) printf(" %08x",Flags);
      else            printf("        ");
      printf(" %2i %2i %1i %2i  %8.3f  %8.3f %9.3f %6.3f %8.3f %6.3f %10i %10i %6i %10i %8.3f %8.3f %8.3f\n",
	     straw->id().getPlane(),
	     straw->id().getPanel(),
	     straw->id().getLayer(),
	     straw->id().getStraw(),
	     x, y, z, phi,
	     Hit->time(),
	     Hit->energyDep(),
	     pdg_id,
	     mother_pdg_id,
	     generator_id,
	     sim_id,
	     mc_mom,
	     mc_pT,
	     mc_pZ);
    }
  }


//-----------------------------------------------------------------------------
void TAnaDump::printStrawHitCollection(const char* ModuleLabel, 
				       const char* ProductName,
				       const char* ProcessName,
				       double TMin, double TMax) {

  const char* oname = "TAnaDump::printStrawHitCollection";

  art::Handle<mu2e::StrawHitCollection> shcHandle;
  const mu2e::StrawHitCollection*       shc;

  // art::Handle<mu2e::StrawHitFlagCollection> shflagH;
  // const mu2e::StrawHitFlagCollection*       shfcol;
  

  art::Handle<mu2e::ComboHitCollection> chcHandle;
  const mu2e::ComboHitCollection*       chcol;
//-----------------------------------------------------------------------------
// get straw hits
//-----------------------------------------------------------------------------
  if (ProductName[0] != 0) {
    art::Selector  selector(art::ProductInstanceNameSelector(ProductName) &&
			    art::ProcessNameSelector(ProcessName)         && 
			    art::ModuleLabelSelector(ModuleLabel)            );
    fEvent->get(selector, shcHandle);
    fEvent->get(selector, chcHandle);
  }
  else {
    art::Selector  selector(art::ProcessNameSelector(ProcessName)         && 
			    art::ModuleLabelSelector(ModuleLabel)            );
    fEvent->get(selector, shcHandle);
    fEvent->get(selector, chcHandle);
  }

  if (shcHandle.isValid()) shc = shcHandle.product();
  else {
    printf(">>> ERROR in %s: Straw Hit Collection by \"%s\" doesn't exist. Bail Out.\n",
	   oname,ModuleLabel);
    return;
  }
//-----------------------------------------------------------------------------
// get straw hit flags (half-hack)
//-----------------------------------------------------------------------------

  if (chcHandle.isValid()) chcol = chcHandle.product();
  else {
    printf(">>> ERROR in %s: ComboHit Collection by \"%s\" doesn't exist. Bail Out.\n",
 	   oname,ModuleLabel);
    return;
  }

  _mcdigis = &fEvent->getByLabel<mu2e::StrawDigiMCCollection>(fStrawDigiMCCollTag.Data());
  if (_mcdigis == nullptr) {
    printf(">>> ERROR in %s: failed to locate StrawDigiMCCollection with tag=%s, BAIL OUT.\n",
	   oname,fStrawDigiMCCollTag.Data());
    return;
  }
 
  int nhits = shc->size();

  const mu2e::StepPointMC*  step;
  const mu2e::StrawHit*     hit;
  int                       flags;
  int                       banner_printed(0);
  for (int i=0; i<nhits; i++) {
    const mu2e::StrawDigiMC* mcdigi = &_mcdigis->at(i);

    if (mcdigi->wireEndTime(mu2e::StrawEnd::cal) < mcdigi->wireEndTime(mu2e::StrawEnd::hv)) {
      step = mcdigi->stepPointMC(mu2e::StrawEnd::cal).get();
    }
    else {
      step = mcdigi->stepPointMC(mu2e::StrawEnd::hv ).get();
    }

    hit   = &shc->at(i);
					// assuming it doesn't move beyond 32 bits
    flags = *((int*) &chcol->at(i).flag());
    if (banner_printed == 0) {
      printStrawHit(hit,step,"banner");
      banner_printed = 1;
    }
    if ((hit->time() >= TMin) && (hit->time() <= TMax)) {
      printStrawHit(hit,step,"data",i,flags);
    }
  }
 
}



//-----------------------------------------------------------------------------
// FlagBgrHitsCollName = 'StrawHits' or 'ComboHits'. Very unfortunate choice!
//-----------------------------------------------------------------------------
void TAnaDump::printComboHitCollection(const char* StrawHitCollTag   , 
				       const char* FlagBgrHitsCollTag,
				       const char* StrawDigiMCCollTag,
				       const char* ProcessName       ,
				       double TMin, double TMax) {

  //  const char* oname = "TAnaDump::printComboHitCollection";
//-----------------------------------------------------------------------------
// get straw hits
//-----------------------------------------------------------------------------
  art::Handle<mu2e::ComboHitCollection> shcH;
  const mu2e::ComboHitCollection* shc(nullptr);
  fEvent->getByLabel(StrawHitCollTag,shcH);

  if (shcH.isValid()) shc = shcH.product();
  else {
    printf("ERROR: cant find StrawHitCollection tag=%s, EXIT\n",StrawHitCollTag);
    return;
  }

  art::Handle<mu2e::StrawHitFlagCollection> shfcH;
  const mu2e::StrawHitFlagCollection*       shfcol(nullptr);
  fEvent->getByLabel<mu2e::StrawHitFlagCollection>(FlagBgrHitsCollTag,shfcH);
  if (shfcH.isValid()) shfcol = shfcH.product();
  else {
    printf("ERROR: cant find StrawHitFlagCollection tag=%s, avalable collections are:\n",FlagBgrHitsCollTag);

    vector<art::Handle<mu2e::StrawHitFlagCollection>> list;
    const art::Provenance*                            prov;
    const art::Handle<mu2e::StrawHitFlagCollection>*  handle;

    art::Selector  selector(art::ProductInstanceNameSelector(""));
    fEvent->getMany(selector,list);

    for (auto it = list.begin(); it != list.end(); it++) {
      handle = it.operator -> ();
      if (handle->isValid()) {
	prov = handle->provenance();
	
	printf("moduleLabel: %-20s, productInstanceName: %-20s, processName:= %-30s\n" ,
	       prov->moduleLabel().data(),
	       prov->productInstanceName().data(),
	       prov->processName().data()
	       );
      }
    }
    return;
  }

  art::Handle<mu2e::StrawDigiMCCollection> mcdH;
  fEvent->getByLabel<mu2e::StrawDigiMCCollection>(StrawDigiMCCollTag,mcdH);
  const mu2e::StrawDigiMCCollection*  mcdigis(nullptr);
  if (mcdH.isValid())   mcdigis = mcdH.product();
  else {
    printf("ERROR: cant find StrawDigiMCCollection tag=%s, EXIT\n",StrawDigiMCCollTag);
    return;
  }

  int nhits = shc->size();

  const mu2e::ComboHit* hit;

  int flags;
  
  const mu2e::ComboHit* hit0 = &shc->at(0);
 
  int banner_printed = 0;
  for (int i=0; i<nhits; i++) {
    hit         = &shc->at(i);
    size_t ish  = hit-hit0;
    vector<StrawDigiIndex> shids;
    shc->fillStrawDigiIndices(*(fEvent),ish,shids);


    const mu2e::StrawDigiMC* mcdigi = &mcdigis->at(shids[0]);
    const mu2e::StepPointMC* Step;

    if (mcdigi->wireEndTime(mu2e::StrawEnd::cal) < mcdigi->wireEndTime(mu2e::StrawEnd::hv)) {
      Step = mcdigi->stepPointMC(mu2e::StrawEnd::cal).get();
    }
    else {
      Step = mcdigi->stepPointMC(mu2e::StrawEnd::hv ).get();
    }
					// assuming it doesn't move beyond 32 bits
    flags = *((int*) &shfcol->at(i));
    if (banner_printed == 0) {
      printComboHit(hit, Step, "banner");
      banner_printed = 1;
    }
    if ((hit->time() >= TMin) && (hit->time() <= TMax)) {
      printComboHit(hit, Step, "data", i, flags);
    }
  }
 
}



//-----------------------------------------------------------------------------------------------------
// if PrintData non-null, it is a pointer to the index ..
//-----------------------------------------------------------------------------------------------------
void TAnaDump::printSimParticle(const mu2e::SimParticle* P, const char* Opt, const void* PrintData) {

    TString opt = Opt;

    if ((opt == "") || (opt == "banner")) {
      printf("-----------------------------------------------------------------------------------------");
      printf("-----------------------------------------");
      printf("---------------------------------------------------------------------------------------\n");
      printf("Index Primary     ID Parent     GenpID        PDG      X0          Y0         Z0         ");
      printf("T0       Px0       Py0      Pz0        E0 ");
      printf("        X1         Y1           Z1        T1         Px1       Py1      Pz1        E1  \n");
      printf("-----------------------------------------------------------------------------------------");
      printf("------------------------------------------");
      printf("---------------------------------------------------------------------------------------\n");
    }
 
    if ((opt == "") || (opt == "data")) {
      int  id        = P->id().asInt();
      
      int  parent_id (-1);

      if (P->parent()) {
	parent_id = P->parent()->id().asInt();
      }
      int  pdg_id    = P->pdgId();
      int  primary   = P->isPrimary();

      int index (-1.);
      if (PrintData) index = *((int*) PrintData);

      printf("%5i %7i %6i %6i %10i %10i",
	     index, primary, id, parent_id, 
	     P->generatorIndex(), pdg_id);

      printf(" %10.3f %10.3f %10.3f %9.3f %9.3f %9.3f %9.3f %9.3f",
	     P->startPosition().x(),
	     P->startPosition().y(),
	     P->startPosition().z(),
	     P->startGlobalTime(),
	     P->startMomentum().x(),
	     P->startMomentum().y(),
	     P->startMomentum().z(),
	     P->startMomentum().e());

      printf(" %10.3f %10.3f %10.3f %10.3f %9.3f %9.3f %9.3f %9.3f\n",
	     P->endPosition().x(),
	     P->endPosition().y(),
	     P->endPosition().z(),
	     P->endGlobalTime(),
	     P->endMomentum().x(),
	     P->endMomentum().y(),
	     P->endMomentum().z(),
	     P->endMomentum().e());
    }
  }



//-----------------------------------------------------------------------------
void TAnaDump::printSimParticleCollection(const char* ModuleLabel, 
					  const char* ProductName, 
					  const char* ProcessName) {

  art::Handle<mu2e::SimParticleCollection> handle;
  const mu2e::SimParticleCollection*       coll(0);
  const mu2e::SimParticle*                 simp(0);

  if (ProductName[0] != 0) {
    art::Selector  selector(art::ProductInstanceNameSelector(ProductName) &&
			    art::ProcessNameSelector        (ProcessName) && 
			    art::ModuleLabelSelector        (ModuleLabel)    );
    fEvent->get(selector, handle);
  }
  else {
    art::Selector  selector(art::ProcessNameSelector(ProcessName)         && 
			    art::ModuleLabelSelector(ModuleLabel)            );
    fEvent->get(selector, handle);
  }

  if (handle.isValid()) coll = handle.product();
  else {
    printf(">>> ERROR in TAnaDump::printSimParticleCollection: failed to locate collection");
    printf(". BAIL OUT. \n");
    return;
  }

  int banner_printed(0);

  int np = coll->size();

  int i = 0;
  for ( mu2e::SimParticleCollection::const_iterator j=coll->begin(); j != coll->end(); ++j) {
    simp = &j->second;

    if (banner_printed == 0) {
      printSimParticle(simp,"banner",&i);
      banner_printed = 1;
    }
    printSimParticle(simp,"data",&i);
    i++;
  }

  if (i != np) {
    printf(" inconsistency in TAnaDump::printSimParticleCollection\n");
  }
}


//-----------------------------------------------------------------------------
void TAnaDump::printStepPointMC(const mu2e::StepPointMC* Step, const char* Detector, const char* Opt) {
  const char* oname = "TAnaDump::printStepPointMC";
    TString opt = Opt;

    if ((opt == "") || (opt.Index("banner") >= 0)) {
      printf("---------------------------------------------------------------------------------------------");
      printf("----------------------------");
      printf("--------------------------------------------------------------------------------------------------------------------\n");
      printf("  Vol          PDG    ID GenIndex PPdg ParentID      X          Y          Z          T      ");
      printf("  X0          Y0         Z0 ");
      printf("  Edep(Tot) Edep(NI)  Edep(I)    Step  EndCode  Energy    EKin     Mom       Pt    doca   Creation       StopProc   \n");
      printf("---------------------------------------------------------------------------------------------");
      printf("----------------------------");
      printf("--------------------------------------------------------------------------------------------------------------------\n");
    }

    mu2e::GeomHandle<mu2e::Tracker> ttHandle;
    const mu2e::Tracker* tracker = ttHandle.get();
  
    art::Ptr<mu2e::SimParticle> const& simptr = Step->simParticle();
    const mu2e::SimParticle* sim  = simptr.operator ->();
    if (sim == NULL) {
      printf(">>> ERROR: %s sim == NULL\n",oname);
    }

    art::Ptr<mu2e::SimParticle> const& parentptr = sim->parent();

    int parent_pdg_id (-1), parent_id(-1);

    const mu2e::SimParticle* par = parentptr.get();
    if (par != NULL) {
      parent_pdg_id = (int) par->pdgId();
      parent_id     = (int) par->id().asInt();
    }

    double doca = -9999.;
    if ((strcmp(Detector,"tracker") == 0) && tracker) {
      const mu2e::Straw* straw = &tracker->getStraw(mu2e::StrawId(Step->volumeId()));

      const Hep3Vector* v1 = &straw->getMidPoint();
      HepPoint p1(v1->x(),v1->y(),v1->z());

      const Hep3Vector* v2 = &Step->position();
      HepPoint    p2(v2->x(),v2->y(),v2->z());

      TrkLineTraj trstraw(p1,straw->getDirection()  ,0.,0.);
      TrkLineTraj trstep (p2,Step->momentum().unit(),0.,0.);

    // 2015-02-16 G. Pezzu and Murat change in the print out to be finished
    // 2015-02-25 P.Murat: fix sign - trajectory is the first !
    //  however, the sign of the disptance of closest approach is invariant
    // wrt the order
      TrkPoca poca(trstep, 0., trstraw, 0.);
    
      doca = poca.doca();
    }
    
    //    art::Ptr<mu2e::GenParticle> const& apgen = sim->genParticle();
    //    mu2e::GenParticle* gen = (mu2e::GenParticle*) &(*sim->genParticle());

    double mass = sim->startMomentum().m();

    double pabs = Step->momentum().mag();
    double energy = sqrt(pabs*pabs+mass*mass);
    double ekin  = energy-mass;
        
    Hep3Vector mom = Step->momentum();
    double pt = sqrt(pabs*pabs - mom.z()*mom.z());

    art::Handle<mu2e::PhysicalVolumeInfoCollection> volumes;
    fEvent->getRun().getByLabel("g4run", volumes);

//2014-26-11 gianipez added the timeoffsets to the steppoints time

    double stepTime(-9999.);
    if (fTimeOffsets) {
      fTimeOffsets->updateMap(*fEvent);
      stepTime = fTimeOffsets->timeWithOffsetsApplied(*Step);
    }
    else {
      stepTime = Step->time();
    }

    //    const mu2e::PhysicalVolumeInfo& pvinfo = volumes->at(sim->startVolumeIndex());
    //    const mu2e::PhysicalVolumeInfo& pvinfo = volumes->at(Step->volumeId()); - sometimes crashes..

    if ((opt == "") || (opt.Index("data") >= 0)) {
      printf("%5i %12i %6i %5i %5i %7i %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %8.2f %8.2f %8.2f %8.3f %4i %10.3f %8.3f %8.3f %8.3f %7.2f %-12s %-s\n",
	     (int) Step->volumeId(),
	     //	     pvinfo.name().data(), // smth is wrong with the name defined by volumeId()....
	     (int) sim->pdgId(),
	     (int) sim->id().asInt(),
	     (int) sim->generatorIndex(),
	     parent_pdg_id,
	     parent_id,
	     Step->position().x(),
	     Step->position().y(),
	     Step->position().z(),
	     stepTime,                             // Step->time(),
	     sim->startPosition().x(),
	     sim->startPosition().y(),
	     sim->startPosition().z(),
	     Step->totalEDep(),
	     Step->nonIonizingEDep(),
	     Step->ionizingEdep(),
	     Step->stepLength(),
	     Step->endProcessCode().id(),
	     energy,
	     ekin,
	     pabs,
	     pt,
	     doca,
	     sim->creationCode().name().data(),
	     Step->endProcessCode().name().data());
    }
}


//-----------------------------------------------------------------------------
void TAnaDump::printStepPointMCCollection(const char* ModuleLabel, 
					  const char* ProductName,
					  const char* ProcessName) {

  art::Handle<mu2e::StepPointMCCollection> handle;
  const mu2e::StepPointMCCollection*       coll(0);

  if (ProductName[0] != 0) {
    art::Selector  selector(art::ProductInstanceNameSelector(ProductName) &&
			    art::ProcessNameSelector        (ProcessName) && 
			    art::ModuleLabelSelector        (ModuleLabel)    );
    fEvent->get(selector, handle);
  }
  else {
    art::Selector  selector(art::ProcessNameSelector(ProcessName)         && 
			    art::ModuleLabelSelector(ModuleLabel)            );
    fEvent->get(selector, handle);
  }

  if (handle.isValid()) coll = handle.product();
  else {
    printf(">>> ERROR in TAnaDump::printStepPointMCCollection: failed to locate collection");
    printf(". BAIL OUT. \n");
    return;
  }

  int nsteps = coll->size();

  const mu2e::StepPointMC* step;


  int banner_printed = 0;
  for (int i=0; i<nsteps; i++) {
    step = &coll->at(i);
    if (banner_printed == 0) {
      printStepPointMC(step,ProductName,"banner");
      banner_printed = 1;
    }
    printStepPointMC(step,ProductName,"data");
  }
 
}

//-----------------------------------------------------------------------------
void TAnaDump::printStrawHitMCTruth(const mu2e::StrawHitMCTruth* Hit, const char* Opt) {
  TString opt = Opt;
  
  if ((opt == "") || (opt == "banner")) {
    printf("--------------------------------------------------------------------\n");
    printf(" Time Distance DistToMid         dt       eDep \n");
    printf("--------------------------------------------------------------------\n");
  }

  if ((opt == "") || (opt == "data")) {
    printf("%12.5f  %12.5f  %12.5f\n",
	   Hit->driftTime(),
	   Hit->driftDistance(),
	   Hit->distanceToMid());
  }
}


//-----------------------------------------------------------------------------
void TAnaDump::printStrawHitMCTruthCollection(const char* ModuleLabel, 
					      const char* ProductName,
					      const char* ProcessName) {

  art::Handle<mu2e::StrawHitMCTruthCollection> shcHandle;
  const mu2e::StrawHitMCTruthCollection*       shc;

  if (ProductName[0] != 0) {
    art::Selector  selector(art::ProductInstanceNameSelector(ProductName) &&
			    art::ProcessNameSelector(ProcessName)         && 
			    art::ModuleLabelSelector(ModuleLabel)            );
    fEvent->get(selector, shcHandle);
  }
  else {
    art::Selector  selector(art::ProcessNameSelector(ProcessName)         && 
			    art::ModuleLabelSelector(ModuleLabel)            );
    fEvent->get(selector, shcHandle);
  }

  shc = shcHandle.product();

  int nhits = shc->size();

  const mu2e::StrawHitMCTruth* hit;


  int banner_printed = 0;
  for (int i=0; i<nhits; i++) {
    hit = &shc->at(i);
    if (banner_printed == 0) {
      printStrawHitMCTruth(hit, "banner");
      banner_printed = 1;
    }
    printStrawHitMCTruth(hit,"data");
  }
 
}

//-----------------------------------------------------------------------------
void TAnaDump::printStrawHitPosition(const mu2e::StrawHitPosition* Hit, const char* Opt) {
  TString opt = Opt;
  
  if ((opt == "") || (opt == "banner")) {
    printf("---------------------------------------------------------------------------\n");
    printf("  STIndex     X        Y        Z        Wdist     Pres     RRes      Flag \n");
    printf("---------------------------------------------------------------------------\n");
  }

  int flag = *((int*) &Hit->flag());

  double wres = Hit->posRes(mu2e::StrawHitPosition::wire);
  if (wres > 999.) wres = 999.;

  double rres = Hit->posRes(mu2e::StrawHitPosition::trans);
  if (rres > 999.) rres = 999.;

  if ((opt == "") || (opt == "data")) {

    printf("   %6i %8.3f %8.3f %9.3f %8.3f  %7.2f  %7.2f  0x%08x\n",
	   -1,//Hit->stereoHitIndex(),
	   Hit->pos().x(),
	   Hit->pos().y(),
	   Hit->pos().z(),
	   Hit->wireDist(),
	   wres,
	   rres,
	   flag);
  }
}


//-----------------------------------------------------------------------------
void TAnaDump::printStrawHitPositionCollection(const char* ModuleLabel, 
					       const char* ProductName,
					       const char* ProcessName) {

  // art::Handle<mu2e::ComboHitCollection>   spcHandle;
  // const mu2e::ComboHitCollection*         spc;

  // if (ProductName[0] != 0) {
  //   art::Selector  selector(art::ProductInstanceNameSelector(ProductName) &&
  // 			    art::ProcessNameSelector(ProcessName)         && 
  // 			    art::ModuleLabelSelector(ModuleLabel)            );
  //   fEvent->get(selector, spcHandle);
  // }
  // else {
  //   art::Selector  selector(art::ProcessNameSelector(ProcessName)         && 
  // 			    art::ModuleLabelSelector(ModuleLabel)            );
  //   fEvent->get(selector, spcHandle);
  // }

  // spc = spcHandle.product();

  // int nhits = spc->size();

  // const mu2e::ComboHit* pos;


  // int banner_printed = 0;
  // for (int i=0; i<nhits; i++) {
  //   pos = &spc->at(i);
  //   if (banner_printed == 0) {
  //     printStrawHitPosition(pos, "banner");
  //     banner_printed = 1;
  //   }
  //   printStrawHitPosition(pos,"data");
  // }
 
}


//-----------------------------------------------------------------------------
void TAnaDump::printTrackClusterMatch(const mu2e::TrackClusterMatch* Tcm, const char* Opt) {

  TString opt = Opt;
  
  if ((opt == "") || (opt == "banner")) {
    printf("--------------------------------------------------------------------------------------\n");
    printf("  Disk         Cluster          Track         chi2     du        dv       dt       E/P\n");
    printf("--------------------------------------------------------------------------------------\n");
  }

  if ((opt == "") || (opt == "data")) {

    const mu2e::CaloCluster*      cl  = Tcm->caloCluster();
    const mu2e::TrkCaloIntersect* tex = Tcm->textrapol  ();

    int disk     = cl->diskId();
    double chi2  = Tcm->chi2();

    printf("%5i %16p  %16p  %8.3f %8.3f %8.3f %8.3f %8.3f\n",
	   disk,  cl,  tex,  chi2,Tcm->du(),Tcm->dv(),Tcm->dt(),Tcm->ep());
  }
}



//-----------------------------------------------------------------------------
void TAnaDump::printTrackClusterMatchCollection(const char* ModuleLabel, 
						const char* ProductName,
						const char* ProcessName) {

  printf(">>>> ModuleLabel = %s\n",ModuleLabel);

  art::Handle<mu2e::TrackClusterMatchCollection> handle;
  const mu2e::TrackClusterMatchCollection*       coll;

  art::Selector  selector(art::ProductInstanceNameSelector(ProductName) &&
			  art::ProcessNameSelector(ProcessName)         && 
			  art::ModuleLabelSelector(ModuleLabel)            );

  fEvent->get(selector,handle);

  if (handle.isValid()) coll = handle.product();
  else {
    printf(">>> ERROR in TAnaDump::printTrackClusterMatchCollection: failed to locate requested collection. Available:");

    vector<art::Handle<mu2e::TrackClusterMatchCollection>> list_of_handles;

    fEvent->getManyByType(list_of_handles);

    for (auto ih=list_of_handles.begin(); ih<list_of_handles.end(); ih++) {
      printf("%s\n", ih->provenance()->moduleLabel().data());
    }

    printf(". BAIL OUT. \n");
    return;
  }

  int nm = coll->size();

  const mu2e::TrackClusterMatch* obj;

  int banner_printed = 0;

  for (int i=0; i<nm; i++) {
    obj = &coll->at(i);
    if (banner_printed == 0) {
      printTrackClusterMatch(obj, "banner");
      banner_printed = 1;
    }
    printTrackClusterMatch(obj,"data");
  }
 
}


//-----------------------------------------------------------------------------
void TAnaDump::refitTrack(void* Trk, double NSig) {
  KalRep* trk = (KalRep*) Trk;

  const TrkHitVector* hits = &trk->hitVector();

  for (auto it=hits->begin(); it!=hits->end(); it++) {

    // TrkStrawHit inherits from TrkHitOnTrk

    TrkHit* hit = (TrkHit*) &(*it);

    const mu2e::TrkStrawHit* straw_hit = (const mu2e::TrkStrawHit*) (*it);

    double res = straw_hit->resid();

    if (fabs(res) > 0.1*NSig) {
      trk->deactivateHit(hit);
    }
  }

  trk->fit();

  printKalRep(trk, "hits");
}


//-----------------------------------------------------------------------------
// emulate calculation of the unbiased residual
//-----------------------------------------------------------------------------
void TAnaDump::Test_000(const KalRep* Krep, mu2e::TrkStrawHit* Hit) {

//  apparently, Hit has (had ?) once to be on the track ?

  // double             s, slen, rdrift, sflt, tflt, doca/*, sig , xdr*/;
  // const mu2e::Straw *straw;
  // int                sign /*, shId, layer*/;
  // HepPoint           spi , tpi , hpos;
 
  // CLHEP::Hep3Vector  spos, sdir;
  // TrkSimpTraj  *ptraj(NULL);

  fTmp[0] = -1;
  fTmp[1] = -1;

//   KalRep* krep = Krep->clone();

//   straw  = &Hit->straw();
//   //  layer  = straw->id().getLayer();
//   rdrift = Hit->driftRadius();
//   //  shId   = straw->index().asInt();
  
//   //  const KalHit* khit = krep->findHotSite(Hit);

//   s      = Hit->fltLen();

//   //  int active = Hit->isActive();

// //   if (active) krep->deactivateHot(Hit);

// //   krep->resetFit();
// //   krep->fit();
// // 					// local track trajectory
// //   ptraj = krep->localTrajectory(s,slen);

//   vector<KalSite*>::const_iterator itt;
//   int found = 0;
//   for (auto /* vector<KalSite*>::const_iterator */ it=krep->siteList().begin();
//        it!= krep->siteList().end(); it++) {
//     const KalHit* kalhit = (*it)->kalHit();
//     if (kalhit && (kalhit->hitOnTrack() == Hit)) {
//       itt   = it;
//       found = 1;
//       break;
//     }
//   }
      
//   if (found == 0) {
//     ptraj = (TrkSimpTraj  *) krep->localTrajectory(s,slen);
//   }
//   else {
//     krep->smoothedTraj(itt,itt,ptraj);
//   }

//   spos = straw->getMidPoint();
//   sdir = straw->getDirection();
// 					// convert Hep3Vector into HepPoint

//   HepPoint    p1(spos.x(),spos.y(),spos.z());

// 					// wire as a trajectory
//   TrkLineTraj st(p1,sdir,0.,0.);

//   TrkPoca poca = TrkPoca(st,0.,*ptraj,0);

//   Hep3Vector        sdi , tdi, u;

//   sflt = poca.flt1();
//   tflt = poca.flt2();

//   st.getInfo(sflt,spi,sdi);
//   ptraj->getInfo(tflt,tpi,tdi);
      
//   u    = sdi.cross(tdi).unit();  // direction towards the center

//   sign     = Hit->ambig();
//   hpos     = spi+u*rdrift*sign;
// 					// hit residal is positive when its residual vector 
// 					// points inside
//   doca     = (tpi-hpos).dot(u);
//   //  sig      = sqrt(rdrift*rdrift +0.1*0.1); // 2.5; // 1.; // hit[ih]->hitRms();
//   //  xdr      = doca/sig;

//   fTmp[0]  = doca;

  //  if (active) krep->activateHot(Hit);

}
