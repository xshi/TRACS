// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME srcdITWaveDict

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "/home/jcalvopi/TRACS_Concurrency/include/TWaveform.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_TWaveform(void *p = 0);
   static void *newArray_TWaveform(Long_t size, void *p);
   static void delete_TWaveform(void *p);
   static void deleteArray_TWaveform(void *p);
   static void destruct_TWaveform(void *p);
   static void streamer_TWaveform(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TWaveform*)
   {
      ::TWaveform *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TWaveform >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TWaveform", ::TWaveform::Class_Version(), "include/TWaveform.h", 30,
                  typeid(::TWaveform), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TWaveform::Dictionary, isa_proxy, 16,
                  sizeof(::TWaveform) );
      instance.SetNew(&new_TWaveform);
      instance.SetNewArray(&newArray_TWaveform);
      instance.SetDelete(&delete_TWaveform);
      instance.SetDeleteArray(&deleteArray_TWaveform);
      instance.SetDestructor(&destruct_TWaveform);
      instance.SetStreamerFunc(&streamer_TWaveform);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TWaveform*)
   {
      return GenerateInitInstanceLocal((::TWaveform*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::TWaveform*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr TWaveform::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TWaveform::Class_Name()
{
   return "TWaveform";
}

//______________________________________________________________________________
const char *TWaveform::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TWaveform*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TWaveform::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TWaveform*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TWaveform::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TWaveform*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TWaveform::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TWaveform*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void TWaveform::Streamer(TBuffer &R__b)
{
   // Stream an object of class TWaveform.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      TObject::Streamer(R__b);
      R__b >> LPower;
      R__b >> LNph;
      R__b >> Vmax;
      R__b >> Vmin;
      R__b >> iVmax;
      R__b >> iVmin;
      R__b >> tVmax;
      R__b >> tVmin;
      R__b >> tleft;
      R__b >> tright;
      R__b >> trms;
      R__b >> itleft;
      R__b >> itright;
      R__b >> BlineMean;
      R__b >> BlineRMS;
      R__b >> Q50;
      R__b >> Qtot;
      R__b >> RiseTime;
      R__b.CheckByteCount(R__s, R__c, TWaveform::IsA());
   } else {
      R__c = R__b.WriteVersion(TWaveform::IsA(), kTRUE);
      TObject::Streamer(R__b);
      R__b << LPower;
      R__b << LNph;
      R__b << Vmax;
      R__b << Vmin;
      R__b << iVmax;
      R__b << iVmin;
      R__b << tVmax;
      R__b << tVmin;
      R__b << tleft;
      R__b << tright;
      R__b << trms;
      R__b << itleft;
      R__b << itright;
      R__b << BlineMean;
      R__b << BlineRMS;
      R__b << Q50;
      R__b << Qtot;
      R__b << RiseTime;
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TWaveform(void *p) {
      return  p ? new(p) ::TWaveform : new ::TWaveform;
   }
   static void *newArray_TWaveform(Long_t nElements, void *p) {
      return p ? new(p) ::TWaveform[nElements] : new ::TWaveform[nElements];
   }
   // Wrapper around operator delete
   static void delete_TWaveform(void *p) {
      delete ((::TWaveform*)p);
   }
   static void deleteArray_TWaveform(void *p) {
      delete [] ((::TWaveform*)p);
   }
   static void destruct_TWaveform(void *p) {
      typedef ::TWaveform current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_TWaveform(TBuffer &buf, void *obj) {
      ((::TWaveform*)obj)->::TWaveform::Streamer(buf);
   }
} // end of namespace ROOT for class ::TWaveform

namespace {
  void TriggerDictionaryInitialization_TWaveDict_Impl() {
    static const char* headers[] = {
"include/TWaveform.h",
0
    };
    static const char* includePaths[] = {
"~/TRACS_Concurrency/include/",
"/usr/local/root/include",
"/home/jcalvopi/TRACS_Concurrency/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "TWaveDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate(R"ATTRDUMP(ROOT RTTI)ATTRDUMP"))) __attribute__((annotate(R"ATTRDUMP(ROOT RTTI)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$include/TWaveform.h")))  TWaveform;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "TWaveDict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "include/TWaveform.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"TWaveform", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("TWaveDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_TWaveDict_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_TWaveDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_TWaveDict() {
  TriggerDictionaryInitialization_TWaveDict_Impl();
}
