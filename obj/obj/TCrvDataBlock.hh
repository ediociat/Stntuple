#ifndef STNTUPLE_TCrvDataBlock
#define STNTUPLE_TCrvDataBlock

#include "TClonesArray.h"

#include "Stntuple/obj/TStnDataBlock.hh"
#include "Stntuple/obj/TCrvRecoPulse.hh"

#include "Stntuple/mod/InitStntupleDataBlocks.hh"

class TCrvDataBlock: public TStnDataBlock {
  friend class InitCrvDataBlock;

public:
  int            fNPulses;	        // number of reconstructed pulses
  TClonesArray*  fListOfPulses;	        // list of hits
//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
public:
					// ****** constructors and destructor
  TCrvDataBlock();
  virtual ~TCrvDataBlock();
					// ****** accessors

  Int_t          NPulses    () { return fNPulses; }
  TCrvRecoPulse* Pulse (int i) { return (TCrvRecoPulse*) fListOfPulses->UncheckedAt(i); }
  
  TClonesArray* GetListOfPulses () { return fListOfPulses; }
//-----------------------------------------------------------------------------
// modifiers
//-----------------------------------------------------------------------------
                                        //Create hit, increse number of hits

  TCrvRecoPulse* NewPulse() { return new ((*fListOfPulses)[fNPulses++]) TCrvRecoPulse(); } 
//-----------------------------------------------------------------------------
// overloaded methods of TObject
//-----------------------------------------------------------------------------
  void Clear(Option_t* opt="");
  void Print(Option_t* opt="") const;

  ClassDef(TCrvDataBlock,1)	// CRV reco data block
};


#endif


