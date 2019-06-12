#ifndef STNTUPLE_TCrvRecoPulse
#define STNTUPLE_TCrvRecoPulse

#include "TClonesArray.h"


class TCrvRecoPulse: public TObject {
public:
  int       fNpe;
  int       fNPePulseHeight;
  int       fBarIndex;
  int       fSipm;

  float     fTime;
  float     fPulseHeight;
  float     fWidth;
  float     fChi2;
  float     fLeTime;
//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
public:
					// ****** constructors and destructor
  TCrvRecoPulse();
  virtual ~TCrvRecoPulse();
					// ****** accessors

//-----------------------------------------------------------------------------
// modifiers
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// overloaded methods of TObject
//-----------------------------------------------------------------------------
  void Clear(Option_t* opt="");
  void Print(Option_t* opt="") const;

  ClassDef(TCrvRecoPulse,1)	         // CRV reco pulse
};


#endif


