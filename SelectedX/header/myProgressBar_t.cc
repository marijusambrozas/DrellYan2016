#define myProgressBar_t_cxx

#include "myProgressBar_t.h"
#include <TROOT.h>
#include <iostream>
#include <math.h>


void myProgressBar_t::Draw(UInt_t iEntry)
{
    if(isSet) {
        if (!iEntry && !deno) {cout << "0 events selected.\n"; return;}
        nume = iEntry;
        if (deno<nume || nume<0) {cout << "Wrong iEntry.\n"; return;}
        else {
            progress = 51*nume/deno;
//            cout << progress << endl;
            if(progress-51/deno<round(progress) && progress>=round(progress)) { // Checks progress so it doesnt redraw the same picture over and over
//                cout << progress-51/deno << "   " << round(progress) << "   " << progress;
                cout << "\r|";
                for (int n=1; n<=50; n++){
                    if (n<=progress) { cout << "+";}
                    else cout << "-";
                }
                cout << "|  " << round(progress*1.97) << "% ";
                cout.flush();
            }
//            else cout << "\r";
            if (nume==deno-1) cout << "\r        Finished.                                                         \n";
        }
    }
    else { cout << "Progress bar was not set.\n"; return; }
}
