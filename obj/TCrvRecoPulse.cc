///////////////////////////////////////////////////////////////////////////////
//  2014-01-26 P.Murat TCrvRecoPulse
///////////////////////////////////////////////////////////////////////////////
#include "TString.h"
#include "TBuffer.h"

#include "Stntuple/obj/TCrvRecoPulse.hh"

ClassImp(TCrvRecoPulse)

//_____________________________________________________________________________
void TCrvRecoPulse::Streamer(TBuffer &R__b) {

  // int nwi = -1; // ((int*) &fTime) - &fIndex;
  // int nwf = -1; // &fMcMomentum - &fTime +1;
  
  if (R__b.IsReading()) {
    //    Version_t R__v = R__b.ReadVersion();
//-----------------------------------------------------------------------------
// curent version: V1
//-----------------------------------------------------------------------------
    // R__b.ReadFastArray(&fIndex,nwi);
    // R__b.ReadFastArray(&fTime ,nwf);
  }
  else {
    R__b.WriteVersion(TCrvRecoPulse::IsA());
    // R__b.WriteFastArray(&fIndex,nwi);
    // R__b.WriteFastArray(&fTime ,nwf);
  } 
}

//_____________________________________________________________________________
TCrvRecoPulse::TCrvRecoPulse(): TObject() {
  Clear();
}

//_____________________________________________________________________________
TCrvRecoPulse::~TCrvRecoPulse() {
}

// //_____________________________________________________________________________
// void TCrvRecoPulse::Set(int Index, float Time, float Dt, float EnergyDep,
// 			int PdgCode, int MotherPdgCode, int GenCode, 
// 			int SimID, float McMomentum) 
// {
 
// }

//_____________________________________________________________________________
void TCrvRecoPulse::Clear(Option_t* opt) {
}

//_____________________________________________________________________________
void TCrvRecoPulse::Print(Option_t* Option) const {
  // print Straw hit properties
  //  printf("Superlayer: %d, Wire: %d, Cell: %d,\n",fSuperLayer,fWire,fCell);
  
  TString opt = Option;
  opt.ToLower();

  if ((opt == "") || (opt.Index("banner") >= 0)) {
    printf("-----------------------------------------------------------------------------------\n");
    printf("  Index  Time   Dt   Energy  PdgCode  PdgCode(M)  GenCode  SimID  McMom     \n");
    printf("-----------------------------------------------------------------------------------\n");
  }

  if (opt == "banner") return;
  
  printf("not implemented yet\n");
  
}
