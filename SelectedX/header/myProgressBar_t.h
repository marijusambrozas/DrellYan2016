#ifndef myProgressBar_t_h
#define myProgressBar_t_h

#include <TROOT.h>
#include <iostream>
#include <math.h>

class myProgressBar_t {
public :
   Double_t nume, deno, progress;
   bool isSet;
   int autoCleanUp;

   myProgressBar_t(UInt_t iMax=0);
   virtual ~myProgressBar_t();

   void SetAutoCleanUp(int yes=1) { autoCleanUp=yes; }
   int GetAutoCleanUp() const { return autoCleanUp; }

   void Draw(UInt_t iEntry=0);
};

#endif

#ifdef myProgressBar_t_cxx

myProgressBar_t::myProgressBar_t(UInt_t iMax) :
    nume(0), deno(0), progress(0), isSet(0),
      autoCleanUp(1)
{
    if (!iMax) return;
    deno = iMax;
    isSet = 1;

}

myProgressBar_t::~myProgressBar_t()
{
   if (!nume && !deno) return;
   if (autoCleanUp) {
     nume = 0;
     deno = 0;
     progress = 0;
     isSet = 0;
   }
}


#undef myProgressBar_t_cxx
#endif // #ifdef myProgressBar_t_cxx
