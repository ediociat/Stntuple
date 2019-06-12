///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////

#include "Stntuple/mod/InitCrvDataBlock.hh"

//-----------------------------------------------------------------------------
// in this case AbsEvent is just not used
//-----------------------------------------------------------------------------
int InitCrvDataBlock::InitDataBlock(TStnDataBlock* Block, AbsEvent* Event, int Mode) {

  //  if (Block->Initialized(fScopeEvent->fEventNumber,fScopeEvent->fRunNumber)) return 0;
  
  TCrvDataBlock* block = (TCrvDataBlock*) Block;

  block->f_EventNumber  = -1; // fScopeEvent->fEventNumber;
  block->f_RunNumber    = -1; // fScopeEvent->fRunNumber;
  block->f_SubrunNumber = -1; // fScopeEvent->fSubrunNumber;
  
  return 0;
}
