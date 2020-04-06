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
#include "./header/myRatioPlot_t.h"
#include "./header/FileMgr.h"

void removeNegativeBins(TH1D *h);

void FR_TemplateTest (Int_t bin=-1, Int_t test=0)
{
    FileMgr fm;
    DYAnalyzer an("");

    TH1D *h_barrel_data_template[nPtBin_ele],    *h_endcap_data_template[nPtBin_ele],
         *h_barrel_data_jetTemplate[nPtBin_ele], *h_endcap_data_jetTemplate[nPtBin_ele],
         *h_barrel_QCD_template[nPtBin_ele],     *h_endcap_QCD_template[nPtBin_ele],
         *h_barrel_QCD_jetTemplate[nPtBin_ele],  *h_endcap_QCD_jetTemplate[nPtBin_ele],
         *h_barrel_MC_template[_SinglePhoton_H][nPtBin_ele],    *h_endcap_MC_template[_SinglePhoton_H][nPtBin_ele],
         *h_barrel_MC_jetTemplate[_SinglePhoton_H][nPtBin_ele], *h_endcap_MC_jetTemplate[_SinglePhoton_H][nPtBin_ele];

    THStack *s_barrel_MC_template[nPtBin_ele],    *s_endcap_MC_template[nPtBin_ele],
            *s_barrel_MC_jetTemplate[nPtBin_ele], *s_endcap_MC_jetTemplate[nPtBin_ele];

    for (Int_t ih=0; ih<nPtBin_ele; ih++)
    {
        s_barrel_MC_template[ih]    = new THStack("s_barrel_MC_template", "");
        s_barrel_MC_jetTemplate[ih] = new THStack("s_barrel_MC_jetTemplate", "");
        s_endcap_MC_template[ih]    = new THStack("s_endcap_MC_template", "");
        s_endcap_MC_jetTemplate[ih] = new THStack("s_endcap_MC_jetTemplate", "");
    }

    for (Process_t pr=_DY_10to50; pr<=_SinglePhoton_H; pr=next(pr))
    {
        TFile *file;
        file = new TFile("/media/sf_DATA/FR/Electron/FR_Hist_E_"+fm.Procname[pr]+".root", "READ");

        if (pr == _SinglePhoton_B)
        {
            for (Int_t ih=0; ih<nPtBin_ele; ih++)
            {
                if (bin >= 0 && ih != bin) continue;
                file->GetObject("h_HoverE_barrel_template_"+TString::Itoa(ih, 10), h_barrel_data_template[ih]);
                file->GetObject("h_HoverE_barrel_jetTemplate_"+TString::Itoa(ih, 10), h_barrel_data_jetTemplate[ih]);
                file->GetObject("h_HoverE_endcap_template_"+TString::Itoa(ih, 10), h_endcap_data_template[ih]);
                file->GetObject("h_HoverE_endcap_jetTemplate_"+TString::Itoa(ih, 10), h_endcap_data_jetTemplate[ih]);

                removeNegativeBins(h_barrel_data_template[ih]);
                removeNegativeBins(h_barrel_data_jetTemplate[ih]);
                removeNegativeBins(h_endcap_data_template[ih]);
                removeNegativeBins(h_endcap_data_jetTemplate[ih]);

                h_barrel_data_template[ih]->SetDirectory(0);
                h_barrel_data_jetTemplate[ih]->SetDirectory(0);
                h_endcap_data_template[ih]->SetDirectory(0);
                h_endcap_data_jetTemplate[ih]->SetDirectory(0);
            }
        }
        else if (pr > _SinglePhoton_B && pr <= _SinglePhoton_H)
        {
            for (Int_t ih=0; ih<nPtBin_ele; ih++)
            {
                if (bin >= 0 && ih != bin) continue;
                TH1D* h_temp[4];
                file->GetObject("h_HoverE_barrel_template_"+TString::Itoa(ih, 10), h_temp[0]);
                file->GetObject("h_HoverE_barrel_jetTemplate_"+TString::Itoa(ih, 10), h_temp[1]);
                file->GetObject("h_HoverE_endcap_template_"+TString::Itoa(ih, 10), h_temp[2]);
                file->GetObject("h_HoverE_endcap_jetTemplate_"+TString::Itoa(ih, 10), h_temp[3]);

                removeNegativeBins(h_temp[0]);
                removeNegativeBins(h_temp[1]);
                removeNegativeBins(h_temp[2]);
                removeNegativeBins(h_temp[3]);

                h_barrel_data_template[ih]->Add(h_temp[0]);
                h_barrel_data_jetTemplate[ih]->Add(h_temp[1]);
                h_endcap_data_template[ih]->Add(h_temp[2]);
                h_endcap_data_jetTemplate[ih]->Add(h_temp[3]);
            }
        }
        else if (pr == _QCDEMEnriched_20to30)
        {
            for (Int_t ih=0; ih<nPtBin_ele; ih++)
            {
                if (bin >= 0 && ih != bin) continue;
                file->GetObject("h_HoverE_barrel_template_"+TString::Itoa(ih, 10),    h_barrel_QCD_template[ih]);
                file->GetObject("h_HoverE_barrel_jetTemplate_"+TString::Itoa(ih, 10), h_barrel_QCD_jetTemplate[ih]);
                file->GetObject("h_HoverE_endcap_template_"+TString::Itoa(ih, 10),    h_endcap_QCD_template[ih]);
                file->GetObject("h_HoverE_endcap_jetTemplate_"+TString::Itoa(ih, 10), h_endcap_QCD_jetTemplate[ih]);

                removeNegativeBins(h_barrel_QCD_template[ih]);
                removeNegativeBins(h_barrel_QCD_jetTemplate[ih]);
                removeNegativeBins(h_endcap_QCD_template[ih]);
                removeNegativeBins(h_endcap_QCD_jetTemplate[ih]);

                h_barrel_QCD_template[ih]->SetDirectory(0);
                h_barrel_QCD_jetTemplate[ih]->SetDirectory(0);
                h_endcap_QCD_template[ih]->SetDirectory(0);
                h_endcap_QCD_jetTemplate[ih]->SetDirectory(0);
            }
        }
        else if (pr > _QCDEMEnriched_20to30 && pr <= _QCDEMEnriched_300toInf)
        {
            for (Int_t ih=0; ih<nPtBin_ele; ih++)
            {
                if (bin >= 0 && ih != bin) continue;
                TH1D* h_temp[4];
                file->GetObject("h_HoverE_barrel_template_"+TString::Itoa(ih, 10), h_temp[0]);
                file->GetObject("h_HoverE_barrel_jetTemplate_"+TString::Itoa(ih, 10), h_temp[1]);
                file->GetObject("h_HoverE_endcap_template_"+TString::Itoa(ih, 10), h_temp[2]);
                file->GetObject("h_HoverE_endcap_jetTemplate_"+TString::Itoa(ih, 10), h_temp[3]);

                removeNegativeBins(h_temp[0]);
                removeNegativeBins(h_temp[1]);
                removeNegativeBins(h_temp[2]);
                removeNegativeBins(h_temp[3]);

                h_barrel_QCD_template[ih]->Add(h_temp[0]);
                h_barrel_QCD_jetTemplate[ih]->Add(h_temp[1]);
                h_endcap_QCD_template[ih]->Add(h_temp[2]);
                h_endcap_QCD_jetTemplate[ih]->Add(h_temp[3]);
            }
        }
        else
        {
            for (Int_t ih=0; ih<nPtBin_ele; ih++)
            {
                if (bin >= 0 && ih != bin) continue;
                file->GetObject("h_HoverE_barrel_template_"+TString::Itoa(ih, 10), h_barrel_MC_template[pr][ih]);
                file->GetObject("h_HoverE_barrel_jetTemplate_"+TString::Itoa(ih, 10), h_barrel_MC_jetTemplate[pr][ih]);
                file->GetObject("h_HoverE_endcap_template_"+TString::Itoa(ih, 10), h_endcap_MC_template[pr][ih]);
                file->GetObject("h_HoverE_endcap_jetTemplate_"+TString::Itoa(ih, 10), h_endcap_MC_jetTemplate[pr][ih]);

                removeNegativeBins(h_barrel_MC_template[pr][ih]);
                removeNegativeBins(h_barrel_MC_jetTemplate[pr][ih]);
                removeNegativeBins(h_endcap_MC_template[pr][ih]);
                removeNegativeBins(h_endcap_MC_jetTemplate[pr][ih]);

                s_barrel_MC_template[ih]->Add(h_barrel_MC_template[pr][ih]);
                s_barrel_MC_jetTemplate[ih]->Add(h_barrel_MC_jetTemplate[pr][ih]);
                s_endcap_MC_template[ih]->Add(h_endcap_MC_template[pr][ih]);
                s_endcap_MC_jetTemplate[ih]->Add(h_endcap_MC_jetTemplate[pr][ih]);
            }
        }

        if (pr == _DY_2000to3000) pr = _EndOf_DYTauTau_Normal; // next -- ttbar
        if (pr == _WJets_ext2v5) pr = _EndOf_QCDMuEnriched_Normal; // next -- QCDEMEnriched
        if (pr == _GJets_2000to5000) pr = _EndOf_SingleElectron_Normal; // next -- SinglePhoton
    } // End of for(process_t)

    TH1D *h_pT_QCD_nume_barrel = new TH1D("h_pT_QCD_nume_barrel", "h_pT_QCD_nume_barrel", nPtBin_ele, an.ptbin_ele);
    TH1D *h_pT_QCD_nume_endcap = new TH1D("h_pT_QCD_nume_endcap", "h_pT_QCD_nume_endcap", nPtBin_ele, an.ptbin_ele);
    TH1D *h_pT_QCD_nume_barrel_original = new TH1D("h_pT_QCD_nume_barrel_original", "h_pT_QCD_nume_barrel_original", nPtBin_ele, an.ptbin_ele);
    TH1D *h_pT_QCD_nume_endcap_original = new TH1D("h_pT_QCD_nume_endcap_original", "h_pT_QCD_nume_endcap_original", nPtBin_ele, an.ptbin_ele);

    for (Int_t ih=0; ih<nPtBin_ele; ih++)
    {
        if (bin >= 0 && ih != bin) continue;
        // Making jet templates
        h_barrel_data_jetTemplate[ih]->Add(((TH1D*)(s_barrel_MC_jetTemplate[ih]->GetStack()->Last())), -1);
        h_endcap_data_jetTemplate[ih]->Add(((TH1D*)(s_endcap_MC_jetTemplate[ih]->GetStack()->Last())), -1);
        h_barrel_data_template[ih]->Add(((TH1D*)(s_barrel_MC_template[ih]->GetStack()->Last())), -1);
        h_endcap_data_template[ih]->Add(((TH1D*)(s_endcap_MC_template[ih]->GetStack()->Last())), -1);

        // Normalizing
        h_barrel_data_jetTemplate[ih]->Scale(h_barrel_data_template[ih]->GetBinContent(5)/h_barrel_data_jetTemplate[ih]->GetBinContent(5));
//        h_endcap_data_jetTemplate[ih]->Scale(h_endcap_data_template[ih]->GetBinContent(5)/h_endcap_data_jetTemplate[ih]->GetBinContent(5));
//        h_barrel_data_jetTemplate[ih]->Scale(h_barrel_data_template[ih]->Integral(4,5)/h_barrel_data_jetTemplate[ih]->Integral(4,5));
        h_endcap_data_jetTemplate[ih]->Scale(h_endcap_data_template[ih]->Integral(4,5)/h_endcap_data_jetTemplate[ih]->Integral(4,5));

        // For test
        h_barrel_QCD_jetTemplate[ih]->Scale(h_barrel_QCD_template[ih]->GetBinContent(5)/h_barrel_QCD_jetTemplate[ih]->GetBinContent(5));
        h_endcap_QCD_jetTemplate[ih]->Scale(h_endcap_QCD_template[ih]->GetBinContent(5)/h_endcap_QCD_jetTemplate[ih]->GetBinContent(5));
//        h_barrel_QCD_jetTemplate[ih]->Scale(h_barrel_QCD_template[ih]->Integral(4,5)/h_barrel_QCD_jetTemplate[ih]->Integral(4,5));
//        h_endcap_QCD_jetTemplate[ih]->Scale(h_endcap_QCD_template[ih]->Integral(3,4)/h_endcap_QCD_jetTemplate[ih]->Integral(3,4));

        // Writing numerator value
        if (bin < 0)
        {
            if (!test)
            {
                h_pT_QCD_nume_barrel->SetBinContent(ih+1, h_barrel_data_jetTemplate[ih]->GetBinContent(1));
                h_pT_QCD_nume_barrel->SetBinError(ih+1,   h_barrel_data_jetTemplate[ih]->GetBinError(1));
                h_pT_QCD_nume_endcap->SetBinContent(ih+1, h_endcap_data_jetTemplate[ih]->GetBinContent(1));
                h_pT_QCD_nume_endcap->SetBinError(ih+1,   h_endcap_data_jetTemplate[ih]->GetBinError(1));

                h_pT_QCD_nume_barrel_original->SetBinContent(ih+1, h_barrel_data_template[ih]->GetBinContent(1));
                h_pT_QCD_nume_barrel_original->SetBinError(ih+1,   h_barrel_data_template[ih]->GetBinError(1));
                h_pT_QCD_nume_endcap_original->SetBinContent(ih+1, h_endcap_data_template[ih]->GetBinContent(1));
                h_pT_QCD_nume_endcap_original->SetBinError(ih+1,   h_endcap_data_template[ih]->GetBinError(1));
            }
            else
            {
                h_pT_QCD_nume_barrel->SetBinContent(ih+1, h_barrel_QCD_jetTemplate[ih]->GetBinContent(1));
                h_pT_QCD_nume_barrel->SetBinError(ih+1,   h_barrel_QCD_jetTemplate[ih]->GetBinError(1));
                h_pT_QCD_nume_endcap->SetBinContent(ih+1, h_endcap_QCD_jetTemplate[ih]->GetBinContent(1));
                h_pT_QCD_nume_endcap->SetBinError(ih+1,   h_endcap_QCD_jetTemplate[ih]->GetBinError(1));

                h_pT_QCD_nume_barrel_original->SetBinContent(ih+1, h_barrel_QCD_template[ih]->GetBinContent(1));
                h_pT_QCD_nume_barrel_original->SetBinError(ih+1,   h_barrel_QCD_template[ih]->GetBinError(1));
                h_pT_QCD_nume_endcap_original->SetBinContent(ih+1, h_endcap_QCD_template[ih]->GetBinContent(1));
                h_pT_QCD_nume_endcap_original->SetBinError(ih+1,   h_endcap_QCD_template[ih]->GetBinError(1));
            }
        }
    }

    if (bin >= 0)
    {
        TCanvas *c_barrel = new TCanvas("c_barrel", "Barrel pT bin "+TString::Itoa(bin, 10), 800, 800);
        h_barrel_data_template[bin]->SetMarkerStyle(kFullDotLarge);
        h_barrel_data_jetTemplate[bin]->SetMarkerStyle(kFullDotLarge);
        h_barrel_data_jetTemplate[bin]->SetLineColor(kRed);
        h_barrel_data_jetTemplate[bin]->SetMarkerColor(kRed);
        h_barrel_QCD_template[bin]->SetMarkerStyle(kFullDotLarge);
        h_barrel_QCD_jetTemplate[bin]->SetMarkerStyle(kFullDotLarge);
        h_barrel_QCD_jetTemplate[bin]->SetLineColor(kRed);
        h_barrel_QCD_jetTemplate[bin]->SetMarkerColor(kRed);

        if (!test)
        {
            h_barrel_data_template[bin]->Draw();
            h_barrel_data_template[bin]->GetYaxis()->SetRangeUser(10, 1e9);
            h_barrel_data_jetTemplate[bin]->Draw("same");
        }
        else
        {
            h_barrel_QCD_template[bin]->Draw();
            h_barrel_QCD_template[bin]->GetYaxis()->SetRangeUser(10, 1e9);
            h_barrel_QCD_jetTemplate[bin]->Draw("same");
        }

        c_barrel->SetLogy();
        c_barrel->Update();

        TCanvas *c_endcap = new TCanvas("c_endcap", "Endcap pT bin "+TString::Itoa(bin, 10), 800, 800);
        h_endcap_data_template[bin]->SetMarkerStyle(kFullDotLarge);
        h_endcap_data_jetTemplate[bin]->SetMarkerStyle(kFullDotLarge);
        h_endcap_data_jetTemplate[bin]->SetLineColor(kRed);
        h_endcap_data_jetTemplate[bin]->SetMarkerColor(kRed);
        h_endcap_QCD_template[bin]->SetMarkerStyle(kFullDotLarge);
        h_endcap_QCD_jetTemplate[bin]->SetMarkerStyle(kFullDotLarge);
        h_endcap_QCD_jetTemplate[bin]->SetLineColor(kRed);
        h_endcap_QCD_jetTemplate[bin]->SetMarkerColor(kRed);

        if (!test)
        {
            h_endcap_data_template[bin]->Draw();
            h_endcap_data_template[bin]->GetYaxis()->SetRangeUser(10, 1e9);
            h_endcap_data_jetTemplate[bin]->Draw("same");
        }
        else
        {
            h_endcap_QCD_template[bin]->Draw();
            h_endcap_QCD_template[bin]->GetYaxis()->SetRangeUser(10, 1e9);
            h_endcap_QCD_jetTemplate[bin]->Draw("same");
        }

        c_endcap->SetLogy();
        c_endcap->Update();
    }
    else
    {
        h_pT_QCD_nume_barrel_original->SetMarkerStyle(kFullDotLarge);
        h_pT_QCD_nume_barrel->SetMarkerStyle(kFullDotLarge);
        h_pT_QCD_nume_barrel->SetLineColor(kWhite);
        h_pT_QCD_nume_barrel->SetMarkerColor(kRed);

        h_pT_QCD_nume_endcap_original->SetMarkerStyle(kFullDotLarge);
        h_pT_QCD_nume_endcap->SetMarkerStyle(kFullDotLarge);
        h_pT_QCD_nume_endcap->SetLineColor(kWhite);
        h_pT_QCD_nume_endcap->SetMarkerColor(kRed);

        TLegend *legend = new TLegend(0.5, 0.75, 0.95, 0.95);
        legend->AddEntry(h_pT_QCD_nume_barrel_original, "Original (A)", "lp");
        legend->AddEntry(h_pT_QCD_nume_barrel, "Normalized template (BC/D)", "lp");

        myRatioPlot_t *RP_barrel = new myRatioPlot_t("RP_barrel", h_pT_QCD_nume_barrel, h_pT_QCD_nume_barrel_original);
        RP_barrel->SetPlots("p_{T} [GeV/c^{2}]", 25, 5000, "AD/BC");
        RP_barrel->ImportLegend(legend);
        RP_barrel->Draw(10, 1e9, 1);
        h_pT_QCD_nume_barrel->SetLineColor(kRed);
        RP_barrel->DrawOnTop(h_pT_QCD_nume_barrel);

        myRatioPlot_t *RP_endcap = new myRatioPlot_t("RP_endcap", h_pT_QCD_nume_endcap, h_pT_QCD_nume_endcap_original);
        RP_endcap->SetPlots("p_{T} [GeV/c^{2}]", 25, 5000, "AD/BC");
        RP_endcap->ImportLegend(legend);
        RP_endcap->Draw(10, 1e9, 1);
        h_pT_QCD_nume_endcap->SetLineColor(kRed);
        RP_endcap->DrawOnTop(h_pT_QCD_nume_endcap);

        if (!test)
        {
            TFile *f_out = new TFile("/media/sf_DATA/FR/Electron/ABCD_hists.root", "RECREATE");
            h_pT_QCD_nume_barrel->Write();
            h_pT_QCD_nume_endcap->Write();
            f_out->Close();
            cout << "Output written to '/media/sf_DATA/FR/Electron/ABCD_hists.root'" << endl;
        }
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
