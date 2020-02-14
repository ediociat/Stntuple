///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////

#include "Stntuple/mod/InitTriggerBlock.hh"
#include "Stntuple/obj/TStnDBManager.hh"
#include "Stntuple/obj/TStnTriggerTable.hh"

#include "Mu2eUtilities/inc/TriggerResultsNavigator.hh"

//-----------------------------------------------------------------------------
// in this case AbsEvent is just not used
//-----------------------------------------------------------------------------
int StntupleInitTriggerBlock::InitDataBlock(TStnDataBlock* Block, AbsEvent* Event, int Mode) {

  int ev, rn, sr;

  ev = Event->event();
  rn = Event->run();
  sr = Event->subRun();

  if (Block->Initialized(ev,rn,sr)) return 0;
  
  TStnTriggerBlock* block = (TStnTriggerBlock*) Block;

  block->f_EventNumber  = ev;
  block->f_RunNumber    = rn;
  block->f_SubrunNumber = sr;
//-----------------------------------------------------------------------------
// 
//-----------------------------------------------------------------------------
  char   trigProcessName[100];

  block->GetModuleLabel("art::TriggerResults", trigProcessName);

  art::InputTag const tag{Form("TriggerResults::%s", trigProcessName)};  

  auto const trH                = Event->getValidHandle<art::TriggerResults>(tag);
  const art::TriggerResults* tr = trH.product();

  mu2e::TriggerResultsNavigator trn(tr);
  
  TStnDBManager*    dbm  = TStnDBManager::Instance();
  TStnTriggerTable* ttbl = (TStnTriggerTable*) dbm->GetTable("TriggerTable");

  TObjArray list(200);

  block->fNPaths = trn.getTrigPaths().size();

//-----------------------------------------------------------------------------
// number of paths is variable, but paths have predefined bit numbers assigned
// assume bit ID <= 32
//-----------------------------------------------------------------------------
  int nbits = 32;

  block->fPaths.Init(nbits);

  for (int i=0; i<block->fNPaths; i++) {
    const std::string& name = trn.getTrigPath(i);

    if (trn.accept(name)) {
//-----------------------------------------------------------------------------
// the event has passed the logic of the trigger path named 'name'
//-----------------------------------------------------------------------------
      char oname[] = "StntupleInitTriggerBlock::InitDataBlock";

      ttbl->GetListOfTriggers(name.data(),&list);
      int nt = list.GetEntriesFast();
      if (nt == 1) {
	TStnTrigger* tr = (TStnTrigger*) list.At(0);
	int bit = tr->Bit();
	block->fPaths.SetBit(bit);
      }
      else if (nt == 0) {
	printf("%s ERROR: path %s is not in the trigger table for this run. The bit is NOT SET\n",oname,name.data());
      }
      else {
	printf("%s ERROR: path %s defined more than once. The bit is NOT SET\n",oname,name.data());
      }
    }
  }

  return 0;
}


//-----------------------------------------------------------------------------
// keep this function as an example, don't really need it
//-----------------------------------------------------------------------------
int StntupleInitTriggerBlock::ResolveLinks(TStnDataBlock* Block, AbsEvent* Event, int Mode) {

  int ev, rn, sr;
  
  ev = Event->event();
  rn = Event->run();
  sr = Event->subRun();
  
  if (! Block->Initialized(ev,rn,sr)) return -1;

  return 0;
}
