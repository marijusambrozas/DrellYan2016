#include <TChain.h>
#include <TTree.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLorentzVector.h>
#include <TStopwatch.h>
#include <TTimeStamp.h>
#include <TString.h>
#include <TLegend.h>
#include <THStack.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TAttMarker.h>
#include <TF1.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TApplication.h>
#include <vector>
#include <TMath.h>
#include <THistPainter.h>
#include <TFormula.h>
#include <iostream>
#include <sstream>

// -- Customized Analyzer for Drel-Yan Analysis -- //
#include "./header/DYAnalyzer.h"
#include "./header/SelectedX.h"
#include "./header/myProgressBar_t.h"
#include "./header/FileMgr.h"

void removeNegativeBins(TH1D *h);

void FR_TemplateTest (Int_t bin=0)
{
    if (bin < 0 || bin > nPtBinBarrel_ele) return;

    FileMgr fm;

    TH1D *h_barrel_data_template[nPtBinBarrel_ele],    *h_endcap_data_template[nPtBinEndcap_ele],
         *h_barrel_data_jetTemplate[nPtBinBarrel_ele], *h_endcap_data_jetTemplate[nPtBinEndcap_ele];

    for (Process_t pr=_SinglePhoton_B; pr<=_SinglePhoton_H; pr=next(pr))
    {
        TFile *file;
        file = new TFile("/media/sf_DATA/FR/Electron/FR_Hist_E_"+fm.Procname[pr]+".root", "READ");

        if (pr == _SinglePhoton_B)
        {
            for (Int_t ih=0; ih<nPtBinBarrel_ele; ih++)
            {
                if (ih != bin) continue;
                file->GetObject("h_HoverE_barrel_template_"+TString::Itoa(ih, 10), h_barrel_data_template[ih]);
                file->GetObject("h_HoverE_barrel_jetTemplate_"+TString::Itoa(ih, 10), h_barrel_data_jetTemplate[ih]);

                removeNegativeBins(h_barrel_data_template[ih]);
                removeNegativeBins(h_barrel_data_jetTemplate[ih]);

                h_barrel_data_template[ih]->SetDirectory(0);
                h_barrel_data_jetTemplate[ih]->SetDirectory(0);

                if (ih < nPtBinEndcap_ele)
                {
                    file->GetObject("h_HoverE_endcap_template_"+TString::Itoa(ih, 10), h_endcap_data_template[ih]);
                    file->GetObject("h_HoverE_endcap_jetTemplate_"+TString::Itoa(ih, 10), h_endcap_data_jetTemplate[ih]);

                    removeNegativeBins(h_endcap_data_template[ih]);
                    removeNegativeBins(h_endcap_data_jetTemplate[ih]);

                    h_endcap_data_template[ih]->SetDirectory(0);
                    h_endcap_data_jetTemplate[ih]->SetDirectory(0);
                }
            }
        }
        else
        {
            for (Int_t ih=0; ih<nPtBinBarrel_ele; ih++)
            {
                if (ih != bin) continue;
                TH1D* h_temp[4];
                file->GetObject("h_HoverE_barrel_template_"+TString::Itoa(ih, 10), h_temp[0]);
                file->GetObject("h_HoverE_barrel_jetTemplate_"+TString::Itoa(ih, 10), h_temp[1]);

                removeNegativeBins(h_temp[0]);
                removeNegativeBins(h_temp[1]);

                h_barrel_data_template[ih]->Add(h_temp[0]);
                h_barrel_data_jetTemplate[ih]->Add(h_temp[1]);

                if (ih < nPtBinEndcap_ele)
                {
                    file->GetObject("h_HoverE_endcap_template_"+TString::Itoa(ih, 10), h_temp[2]);
                    file->GetObject("h_HoverE_endcap_jetTemplate_"+TString::Itoa(ih, 10), h_temp[3]);

                    removeNegativeBins(h_temp[2]);
                    removeNegativeBins(h_temp[3]);

                    h_endcap_data_template[ih]->Add(h_temp[2]);
                    h_endcap_data_jetTemplate[ih]->Add(h_temp[3]);
                }
            }
        }
    }

    h_barrel_data_jetTemplate[bin]->Scale(h_barrel_data_template[bin]->Integral(2,4)/h_barrel_data_jetTemplate[bin]->Integral(2,4));

    TCanvas *c_barrel = new TCanvas("c_barrel", "Barrel pT bin "+TString::Itoa(bin, 10), 800, 800);
    h_barrel_data_template[bin]->SetMarkerStyle(kFullDotLarge);
    h_barrel_data_jetTemplate[bin]->SetMarkerStyle(kFullDotLarge);
    h_barrel_data_jetTemplate[bin]->SetLineColor(kRed);
    h_barrel_data_jetTemplate[bin]->SetMarkerColor(kRed);
    h_barrel_data_template[bin]->Draw("hist");
    h_barrel_data_jetTemplate[bin]->Draw("samehist");
    c_barrel->SetLogy();
    c_barrel->Update();

    if (bin < nPtBinEndcap_ele)
    {
        h_endcap_data_jetTemplate[bin]->Scale(h_endcap_data_template[bin]->Integral(2,4)/h_endcap_data_jetTemplate[bin]->Integral(2,4));

        TCanvas *c_endcap = new TCanvas("c_endcap", "Endcap pT bin "+TString::Itoa(bin, 10), 800, 800);
        h_endcap_data_template[bin]->SetMarkerStyle(kFullDotLarge);
        h_endcap_data_jetTemplate[bin]->SetMarkerStyle(kFullDotLarge);
        h_endcap_data_jetTemplate[bin]->SetLineColor(kRed);
        h_endcap_data_jetTemplate[bin]->SetMarkerColor(kRed);
        h_endcap_data_template[bin]->Draw("hist");
        h_endcap_data_jetTemplate[bin]->Draw("samehist");
        c_endcap->SetLogy();
        c_endcap->Update();
    }

} // End of FR_TemplateTEst()


void removeNegativeBins(TH1D *h)
{
    for (int i=0; i<h->GetSize(); i++)
    {
        if (h->GetBinContent(i) < 0)
        {
            h->SetBinContent(i, 0);
            h->SetBinError(i, 0);
        }
    }
}
