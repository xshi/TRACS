// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME srcdITMeasDict

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
#include "/home/jcalvopi/TestThreads/include/TMeas.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_TMeas(void *p = 0);
   static void *newArray_TMeas(Long_t size, void *p);
   static void delete_TMeas(void *p);
   static void deleteArray_TMeas(void *p);
   static void destruct_TMeas(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TMeas*)
   {
      ::TMeas *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TMeas >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TMeas", ::TMeas::Class_Version(), "include/TMeas.h", 27,
                  typeid(::TMeas), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TMeas::Dictionary, isa_proxy, 4,
                  sizeof(::TMeas) );
      instance.SetNew(&new_TMeas);
      instance.SetNewArray(&newArray_TMeas);
      instance.SetDelete(&delete_TMeas);
      instance.SetDeleteArray(&deleteArray_TMeas);
      instance.SetDestructor(&destruct_TMeas);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TMeas*)
   {
      return GenerateInitInstanceLocal((::TMeas*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::TMeas*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr TMeas::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TMeas::Class_Name()
{
   return "TMeas";
}

//______________________________________________________________________________
const char *TMeas::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TMeas*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TMeas::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TMeas*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TMeas::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TMeas*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TMeas::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TMeas*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void TMeas::Streamer(TBuffer &R__b)
{
   // Stream an object of class TMeas.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TMeas::Class(),this);
   } else {
      R__b.WriteClassBuffer(TMeas::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TMeas(void *p) {
      return  p ? new(p) ::TMeas : new ::TMeas;
   }
   static void *newArray_TMeas(Long_t nElements, void *p) {
      return p ? new(p) ::TMeas[nElements] : new ::TMeas[nElements];
   }
   // Wrapper around operator delete
   static void delete_TMeas(void *p) {
      delete ((::TMeas*)p);
   }
   static void deleteArray_TMeas(void *p) {
      delete [] ((::TMeas*)p);
   }
   static void destruct_TMeas(void *p) {
      typedef ::TMeas current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TMeas

namespace {
  void TriggerDictionaryInitialization_TMeasDict_Impl() {
    static const char* headers[] = {
"include/TMeas.h",
0
    };
    static const char* includePaths[] = {
"~/TestThreads/include/",
"/usr/local/root/include",
"/home/jcalvopi/TestThreads/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "TMeasDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate(R"ATTRDUMP(Edge-TCT data class)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$include/TMeas.h")))  TMeas;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "TMeasDict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "include/TMeas.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"TMeas", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("TMeasDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_TMeasDict_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_TMeasDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_TMeasDict() {
  TriggerDictionaryInitialization_TMeasDict_Impl();
}
