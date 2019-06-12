//
// 2014-01-27 P.Murat
//

#include "Stntuple/obj/TCrvDataBlock.hh"

ClassImp(TCrvDataBlock)

//-----------------------------------------------------------------------------
// R_v is so far unused
//-----------------------------------------------------------------------------
  void TCrvDataBlock::Streamer(TBuffer &R__b) {
  if(R__b.IsReading()) {
    //    Version_t R__v = R__b.ReadVersion();
    R__b.ReadVersion();
    R__b >> fNPulses;
    fListOfPulses->Streamer(R__b);
  }
  else {
    R__b.WriteVersion(TCrvDataBlock::IsA());
    R__b << fNPulses;
    fListOfPulses->Streamer(R__b);
  }
}

//______________________________________________________________________________
TCrvDataBlock::TCrvDataBlock() {
  fListOfPulses = new TClonesArray("TCrvRecoPulse",30240);
  fListOfPulses->BypassStreamer(kFALSE);
  Clear();
}

//______________________________________________________________________________
TCrvDataBlock::~TCrvDataBlock() {
  fListOfPulses->Delete();
  delete fListOfPulses;
}

//______________________________________________________________________________
void TCrvDataBlock::Clear(Option_t* opt) {
  fListOfPulses->Clear();
  fNPulses=0;
}

//______________________________________________________________________________
void TCrvDataBlock::Print(Option_t* opt) const {
  // print all hits in the straw tracker
  printf(" *** print CRV  *** \nNumber of hits: %d\n",fNPulses);
  for(int i=0; i<fNPulses; i++) {
    fListOfPulses->At(i)->Print();
  }
}

