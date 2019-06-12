///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#ifndef __InitCrvDataBlock__
#define __InitCrvDataBlock__

#include "Stntuple/obj/TStnInitDataBlock.hh"
#include "Stntuple/obj/TCrvDataBlock.hh"

class InitCrvDataBlock : public TStnInitDataBlock {
public:
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
public:
  
  virtual int InitDataBlock(TStnDataBlock* Block, AbsEvent* Evt, int Mode);

};

#endif
