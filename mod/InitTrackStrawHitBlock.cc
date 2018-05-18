///////////////////////////////////////////////////////////////////////////////
// 2014-01-26 P.Murat
///////////////////////////////////////////////////////////////////////////////

#include "Stntuple/mod/InitStntupleDataBlocks.hh"
#include "Stntuple/obj/TTrackStrawHitBlock.hh"

#include "RecoDataProducts/inc/StrawHitCollection.hh"

#include "RecoDataProducts/inc/KalRepPtrCollection.hh"
#include "BTrkData/inc/TrkStrawHit.hh"

#include "MCDataProducts/inc/SimParticle.hh"
#include "MCDataProducts/inc/StepPointMC.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"

//-----------------------------------------------------------------------------
Int_t StntupleInitMu2eTrackStrawHitBlock(TStnDataBlock* Block, AbsEvent* AnEvent, int Mode) 
{
  // initialize COT data block with the `event' data
  // return -1, if bank doesn't exist, 0, if everything is OK

  int ev_number, rn_number, nhits;

  mu2e::KalRepPtrCollection*               list_of_kreps(0);

  static char   strh_module_label[100], strh_description[100];
  static char   stmc_module_label[100], stmc_description[100];
  static char   krep_module_label[100], krep_description[100];

  ev_number = AnEvent->event();
  rn_number = AnEvent->run();

  if (Block->Initialized(ev_number,rn_number)) return 0;

  TTrackStrawHitBlock* data = (TTrackStrawHitBlock*) Block;
  data->Clear();
//-----------------------------------------------------------------------------
//  straw hit information
//-----------------------------------------------------------------------------
  data->GetModuleLabel("mu2e::StrawHitCollection",strh_module_label);
  data->GetDescription("mu2e::StrawHitCollection",strh_description );

  data->GetModuleLabel("mu2e::KalRepCollection"  ,krep_module_label);
  data->GetDescription("mu2e::KalRepCollection"  ,krep_description );

  data->GetModuleLabel("mu2e::PtrStepPointMCVectorCollection",stmc_module_label);
  data->GetDescription("mu2e::PtrStepPointMCVectorCollection",stmc_description );

  art::Handle<mu2e::ComboHitCollection>       strh_handle;
  const mu2e::ComboHitCollection*             list_of_hits(0);

  art::Handle<mu2e::PtrStepPointMCVectorCollection> mcptr_handle;
  const mu2e::PtrStepPointMCVectorCollection* list_of_steps(0);

  art::Handle<mu2e::KalRepPtrCollection> krepsHandle;
  if (krep_module_label[0] != 0) {
    if (krep_description[0] == 0) AnEvent->getByLabel(krep_module_label,krepsHandle);
    else                          AnEvent->getByLabel(krep_module_label,krep_description, krepsHandle);
    if (krepsHandle.isValid())    list_of_kreps = (mu2e::KalRepPtrCollection*) krepsHandle.product();
  }

  if (strh_module_label[0] != 0) {
    if (strh_description[0] == 0) AnEvent->getByLabel(strh_module_label,strh_handle);
    else                          AnEvent->getByLabel(strh_module_label,strh_description,strh_handle);

    if (strh_handle.isValid()) list_of_hits = strh_handle.product();

    AnEvent->getByLabel(stmc_module_label,stmc_description,mcptr_handle);
    if (mcptr_handle.isValid()) list_of_steps = mcptr_handle.product();
  }

  if (list_of_hits == NULL) {
    printf(" >>> ERROR in StntupleInitMu2eTrackStrawHitBlock: no list_of_hits. BAIL OUT\n");
    return -1;
  }
//-----------------------------------------------------------------------------
// 
//-----------------------------------------------------------------------------
  const mu2e::StepPointMC  *step;
  const mu2e::SimParticle  *sim;
  const mu2e::Straw        *straw;
    
  TTrackStrawHitData*      hit; 
  data->fNTracks = list_of_kreps->size();

  const KalRep*            krep;
  const TrkHitVector*      krep_hits; 

  int   pdg_id, mother_pdg_id, sim_id, gen_index;
  float mcdoca, mc_mom;

  data->fNTrackHits->Set(data->fNTracks);
  data->fFirst->Set(data->fNTracks);

  const mu2e::ComboHit      *s_hit0(0), *sh(0); 

  int nhtot = 0;

  if (list_of_hits->size() > 0) {

    s_hit0 = &list_of_hits->at(0);

    for (int i=0; i<data->fNTracks; i++) {
      krep      = list_of_kreps->at(i).get();
      krep_hits = &krep->hitVector();
      nhits     = krep_hits->size();
    
      (*data->fNTrackHits)[i] = nhits;
      (*data->fFirst)     [i] = nhtot;

      nhtot += nhits;

      for (int ih=0; ih<nhits; ih++) {
	const mu2e::TrkStrawHit* tsh = static_cast<const mu2e::TrkStrawHit*> (krep_hits->at(ih));
	sh    = &tsh->comboHit();
	straw = &tsh->straw();
	
	int loc   = sh-s_hit0;
      
	hit  = data->NewHit();

	mu2e::PtrStepPointMCVector const& mcptr (list_of_steps->at(loc));
	step = mcptr[0].get();

	if (step) {
	  art::Ptr<mu2e::SimParticle> const& simptr = step->simParticle(); 
	  art::Ptr<mu2e::SimParticle> mother = simptr;
	
	  while (mother->hasParent()) mother = mother->parent();

	  sim = mother.get();

	  pdg_id        = simptr->pdgId();
	  mother_pdg_id = sim->pdgId();
	
	  if (simptr->fromGenerator()) gen_index = simptr->genParticle()->generatorId().id();
	  else                         gen_index = -1;
	
	  sim_id        = simptr->id().asInt();
	  mc_mom        = step->momentum().mag();


	  const CLHEP::Hep3Vector& v1 = straw->getMidPoint();
	  HepPoint p1(v1.x(),v1.y(),v1.z());
	      
	  const CLHEP::Hep3Vector& v2 = step->position();
	  HepPoint    p2(v2.x(),v2.y(),v2.z());
	      
	  TrkLineTraj trstraw(p1,straw->getDirection()  ,0.,0.);
	  TrkLineTraj trstep (p2,step->momentum().unit(),0.,0.);
	      
	  TrkPoca poca(trstep, 0., trstraw, 0.);
	
	  mcdoca = poca.doca();
	}
	else {
	  pdg_id        = -1;
	  mother_pdg_id = -1;
	  gen_index     = -1;
	  sim_id        = -1;
	  mc_mom        = -1.;
	  mcdoca        = 1.e6;
	}
      
	hit->Set(sh->strawId().asUint16(), sh->time(), 0./*sh->dt()*/, sh->energyDep(),//FIXME!
		 tsh->isActive(),tsh->ambig(),tsh->driftRadius(),
		 pdg_id, mother_pdg_id, gen_index, sim_id, 
		 mcdoca, mc_mom);
      }
    }
  }

  data->fNStrawHits   = nhtot;
  
  data->f_RunNumber   = rn_number;
  data->f_EventNumber = ev_number;
  
  return 0;
}

