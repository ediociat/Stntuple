#ifndef stntuple_ana_scripts_modules_hh
#define stntuple_ana_scripts_modules_hh

// #include "Stntuple/ana/TDFCModule.hh"

class TStnCrvAnaModule;
class TStnDebugModule;
class TStnEventDisplayModule;
class TStnGenAnaModule;
class TStnGeneratorModule;		// event generator base class 
class TStnHelixAnaModule;
class TStnLumiMonModule;
class TStnPhotosAnaModule;
class TStnTrackAnaModule;
class TStnTrackSeedAnaModule;
class TStnTriggerAnaModule;
class TStnValidationModule;

namespace stntuple {
  TStnDebugModule*            m_dbg   = nullptr;
  TDFCModule*                 m_dfc   = nullptr;
  TStnCrvAnaModule*           m_crv   = nullptr;
  TStnLumiMonModule*          m_lum   = nullptr;
  TStnEventDisplayModule*     m_evd   = nullptr;
  TStnGenAnaModule*           m_gen   = nullptr;
  TStnGeneratorModule*        m_stg   = nullptr;
  TStnHelixAnaModule*         m_hel   = nullptr;
  TStnPhotosAnaModule*        m_pho   = nullptr;
  TStnTrackAnaModule*         m_trk   = nullptr;
  TStnTrackSeedAnaModule*     m_trs   = nullptr;
  TStnTriggerAnaModule*       m_trig  = nullptr;
  TStnValidationModule*       m_val   = nullptr;
};

#endif
