#define myProgressBar_t_cxx

#include "myProgressBar_t.h"
#include <TROOT.h>
#include <iostream>
#include <math.h>


void myProgressBar_t::Draw(UInt_t iEntry)
{
    if(isSet) {
        if (!iEntry && !deno) {std::cout << "0 events selected.\n"; return;}
        nume = iEntry;
        if (deno<nume || nume<0) {std::cout << "Wrong iEntry.\n"; return;}
        else {
            progress = 51*nume/deno;
//            std::cout << progress << endl;
            if(progress-51/deno<round(progress) && progress>=round(progress)) { // Checks progress so it doesnt redraw the same picture over and over
//                std::cout << progress-51/deno << "   " << round(progress) << "   " << progress;
                std::cout << "\r|";
                for (int n=1; n<=50; n++){
                    if (n<=progress) { std::cout << "+";}
                    else std::cout << "-";
                }
                std::cout << "|  " << round(progress*1.97) << "% ";
                std::cout.flush();
            }
//            else std::cout << "\r";
            if (nume==deno-1) std::cout << "\r        Selection finished.                                                         \n";
        }
    }
    else { std::cout << "Progress bar was not set.\n"; return; }
}
