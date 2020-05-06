#include "header/LocalFileMgr.h"
#include "header/myRatioPlot_t.h"
#include <TChain.h>
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
#include <TVectorT.h>

void ee_Est(Int_t UseFR);
void MuMu_Est(Int_t UseFR);
void removeNegativeBins(TH1D *h);

void EstBkg_EMu(TString whichX = "EE", Int_t UseFR = 1)
{
    TString WhichX = whichX;
    WhichX.ToUpper();

    if (WhichX.Contains("EE"))
        ee_Est(UseFR);
    if (WhichX.Contains("MUMU"))
        MuMu_Est(UseFR);
}


void ee_Est(Int_t UseFR)
{
    LocalFileMgr Mgr;

    Mgr.SetProc(_EE_Bkg_Full);
    cout << "EE Hists location: " << Mgr.HistLocation << endl;

    TFile* f_Bkg_ee = new TFile(Mgr.HistLocation+"Hist_"+Mgr.Procname[_EE_Bkg_Full]+".root", "READ");
    if (f_Bkg_ee->IsOpen()) std::cout << "File Hist_"+Mgr.Procname[_EE_Bkg_Full]+".root opened successfully\n";
    else { std::cout << "File Hist_"+Mgr.Procname[_EE_Bkg_Full]+".root was not opened\n"; return; }

    Mgr.SetProc(_EE_DY_Full);
    TFile* f_DY_ee = new TFile(Mgr.HistLocation+"Hist_"+Mgr.Procname[_EE_DY_Full]+".root", "READ");
    if (f_DY_ee->IsOpen()) std::cout << "File Hist_"+Mgr.Procname[_EE_DY_Full]+".root opened successfully\n";
    else { std::cout << "File Hist_"+Mgr.Procname[_EE_DY_Full]+".root was not opened\n"; return; }

    Mgr.SetProc(_EE_DoubleEG_Full);
    TFile* f_Data_ee = new TFile(Mgr.HistLocation+"Hist_"+Mgr.Procname[_EE_DY_Full]+".root", "READ");
    if (f_Data_ee->IsOpen()) std::cout << "File Hist_"+Mgr.Procname[_EE_DY_Full]+".root opened successfully\n";
    else { std::cout << "File Hist_"+Mgr.Procname[_EE_DY_Full]+".root was not opened\n"; return; }

    Mgr.SetProc(_EMu_Bkg_Full);
    cout << "EMu Hists location: " << Mgr.HistLocation << endl;

    TFile* f_Bkg_EMu = new TFile(Mgr.HistLocation+"Hist_"+Mgr.Procname[_EMu_Bkg_Full]+".root", "READ");
    if (f_Bkg_EMu->IsOpen()) std::cout << "File Hist_"+Mgr.Procname[_EMu_Bkg_Full]+".root opened successfully\n";
    else { std::cout << "File Hist_"+Mgr.Procname[_EMu_Bkg_Full]+".root was not opened\n"; return; }

    Mgr.SetProc(_EMu_SingleMuon_Full);
    TFile* f_Data_EMu = new TFile(Mgr.HistLocation+"Hist_"+Mgr.Procname[_EMu_SingleMuon_Full]+".root", "READ");
    if (f_Data_EMu->IsOpen()) std::cout << "File Hist_"+Mgr.Procname[_EMu_SingleMuon_Full]+".root opened successfully\n";
    else { std::cout << "File Hist_"+Mgr.Procname[_EMu_SingleMuon_Full]+".root was not opened\n"; return; }

//********************************* E MU ***********************************************************

//--------------------------------- Get EMu QCD and WJets from FR -----------------------------------

    TH1D *h_EMu_QCDFR_invm, *h_EMu_SS_QCDFR_invm;
    TH1D *h_EMu_WJetsFR_invm, *h_EMu_SS_WJetsFR_invm;
    if (UseFR)
    {
        TFile *f_QCD_emu = new TFile(Mgr.HistLocation+"EstQCD_EMu.root", "READ");
        if (f_QCD_emu->IsOpen()) std::cout << "File EstQCD_EMu.root opened successfully\n";
        else { std::cout << "File EstQCD_EMu.root was not opened\n"; return; }
        f_QCD_emu->GetObject("h_QCD_est", h_EMu_QCDFR_invm);
        f_QCD_emu->GetObject("h_QCD_est_SS", h_EMu_SS_QCDFR_invm);
        h_EMu_QCDFR_invm->SetFillColor(kRed+3);
        h_EMu_QCDFR_invm->SetLineColor(kRed+3);
        h_EMu_QCDFR_invm->SetDirectory(0);
        h_EMu_SS_QCDFR_invm->SetFillColor(kRed+3);
        h_EMu_SS_QCDFR_invm->SetLineColor(kRed+3);
        h_EMu_SS_QCDFR_invm->SetDirectory(0);
        f_QCD_emu->Close();
        if (!f_QCD_emu->IsOpen()) std::cout << "File EstQCD_EMu.root closed successfully\n";

        TFile *f_WJets_emu = new TFile(Mgr.HistLocation+"EstWJets_EMu.root", "READ");
        if (f_WJets_emu->IsOpen()) std::cout << "File EstWJets_EMu.root opened successfully\n";
        else { std::cout << "File EstWJets_EMu.root was not opened\n"; return; }
        f_WJets_emu->GetObject("h_WJET_est", h_EMu_WJetsFR_invm);
        f_WJets_emu->GetObject("h_WJET_est_SS", h_EMu_SS_WJetsFR_invm);
        h_EMu_WJetsFR_invm->SetFillColor(kRed-2);
        h_EMu_WJetsFR_invm->SetLineColor(kRed-2);
        h_EMu_WJetsFR_invm->SetDirectory(0);
        h_EMu_SS_WJetsFR_invm->SetFillColor(kRed-2);
        h_EMu_SS_WJetsFR_invm->SetLineColor(kRed-2);
        h_EMu_SS_WJetsFR_invm->SetDirectory(0);
        f_WJets_emu->Close();
        if (!f_WJets_emu->IsOpen()) std::cout << "File EstWJets_EMu.root closed successfully\n";
    }

//---------------------------------- MC -----------------------------------------------------------

    TH1D* h_EMu_invm[_EndOf_EMu];
    TH1D* h_EMu_invm2[_EndOf_EMu];
    TH1D* h_EMu_SS_invm[_EndOf_EMu];
    TH1D* h_EMu_SS_invm2[_EndOf_EMu];
    THStack* s_EMu_invm = new THStack("s_emu_invm", "");
    THStack* s_EMu_invm2 = new THStack("s_emu_invm2", "");
    THStack* s_EMu_SS_invm = new THStack("s_emu_SS_invm", "");
    THStack* s_EMu_SS_invm2 = new THStack("s_emu_SS_invm2", "");
    gStyle->SetOptStat(0);

    Int_t isWJ = 1; // Add MC WJets if ==1

    for (SelProc_t pr=_EMu_WJets_Full; pr>_EndOf_EMu_ttbar_Normal; pr=SelProc_t((int)(pr-1)))
    {
        Mgr.SetProc(pr);
        if (pr==_EndOf_EMu_VVnST_Normal) continue;
        if (isWJ == 0 && pr == _EMu_WJets_Full) {pr = _EndOf_EMu_VVnST_Normal; continue;}
        f_Bkg_EMu->GetObject("h_emu_mass_"+Mgr.Procname[pr], h_EMu_invm[pr]);
        f_Bkg_EMu->GetObject("h_emu_mass2_"+Mgr.Procname[pr], h_EMu_invm2[pr]);
        removeNegativeBins(h_EMu_invm[pr]);
        removeNegativeBins(h_EMu_invm2[pr]);

        Color_t color = kBlack;
        if (pr == _EMu_WJets_Full) color = kRed - 2;
        if (pr == _EMu_WW) color = kMagenta - 5;
        if (pr == _EMu_WZ) color = kMagenta - 2;
        if (pr == _EMu_ZZ) color = kMagenta - 6;
        if (pr == _EMu_tbarW) color = kGreen - 2;
        if (pr == _EMu_tW) color = kGreen + 2;
        if (pr == _EMu_ttbar_Full) color = kCyan + 2;
        if (pr == _EMu_DYTauTau_Full) color = kOrange - 5;

        h_EMu_invm[pr]->SetFillColor(color);
        h_EMu_invm[pr]->SetLineColor(color);
        h_EMu_invm[pr]->SetDirectory(0);
        s_EMu_invm->Add(h_EMu_invm[pr]);
        h_EMu_invm2[pr]->SetFillColor(color);
        h_EMu_invm2[pr]->SetLineColor(color);
        h_EMu_invm2[pr]->SetDirectory(0);
        s_EMu_invm2->Add(h_EMu_invm2[pr]);

        f_Bkg_EMu->GetObject("h_emuSS_mass_"+Mgr.Procname[pr], h_EMu_SS_invm[pr]);
        f_Bkg_EMu->GetObject("h_emuSS_mass2_"+Mgr.Procname[pr], h_EMu_SS_invm2[pr]);
        removeNegativeBins(h_EMu_SS_invm[pr]);
        removeNegativeBins(h_EMu_SS_invm2[pr]);

        h_EMu_SS_invm[pr]->SetFillColor(color);
        h_EMu_SS_invm[pr]->SetLineColor(color);
        h_EMu_SS_invm[pr]->SetDirectory(0);
        s_EMu_SS_invm->Add(h_EMu_SS_invm[pr]);
        h_EMu_SS_invm2[pr]->SetFillColor(color);
        h_EMu_SS_invm2[pr]->SetLineColor(color);
        h_EMu_SS_invm2[pr]->SetDirectory(0);
        s_EMu_SS_invm2->Add(h_EMu_SS_invm2[pr]);

        if (pr == _EMu_WJets_Full) pr = _EndOf_EMu_VVnST_Normal; // next - WW
        if (pr == _EMu_tW) pr = _EMu_VVnST; // next -- ttbar
        if (pr == _EMu_DYTauTau_Full) break;
    }

//----------------------------- Data -------------------------------------------------------------

    Mgr.SetProc(_EMu_SingleMuon_Full);
    TH1D *h_EMu_data_invm, *h_EMu_data_invm2;
    f_Data_EMu->GetObject("h_emu_mass_"+Mgr.Procname[_EMu_SingleMuon_Full], h_EMu_data_invm);
    f_Data_EMu->GetObject("h_emu_mass2_"+Mgr.Procname[_EMu_SingleMuon_Full], h_EMu_data_invm2);
    h_EMu_data_invm->SetLineColorAlpha(kBlack, 0);
    h_EMu_data_invm->SetMarkerStyle(kFullDotLarge);
    h_EMu_data_invm->SetMarkerColor(kBlack);
    h_EMu_data_invm->SetDirectory(0);
    h_EMu_data_invm2->SetLineColorAlpha(kBlack, 0);
    h_EMu_data_invm2->SetMarkerStyle(kFullDotLarge);
    h_EMu_data_invm2->SetMarkerColor(kBlack);
    h_EMu_data_invm2->SetDirectory(0);

    TH1D *h_EMu_SS_data_invm, *h_EMu_SS_data_invm2;
    f_Data_EMu->GetObject("h_emuSS_mass_"+Mgr.Procname[_EMu_SingleMuon_Full], h_EMu_SS_data_invm);
    f_Data_EMu->GetObject("h_emuSS_mass2_"+Mgr.Procname[_EMu_SingleMuon_Full], h_EMu_SS_data_invm2);
    h_EMu_SS_data_invm->SetLineColorAlpha(kBlack, 0);
    h_EMu_SS_data_invm->SetMarkerStyle(kFullDotLarge);
    h_EMu_SS_data_invm->SetMarkerColor(kBlack);
    h_EMu_SS_data_invm->SetDirectory(0);
    h_EMu_SS_data_invm2->SetLineColorAlpha(kBlack, 0);
    h_EMu_SS_data_invm2->SetMarkerStyle(kFullDotLarge);
    h_EMu_SS_data_invm2->SetMarkerColor(kBlack);
    h_EMu_SS_data_invm2->SetDirectory(0);

//------------------------------ Est EMu QCD events ---------------------------------------------------

    TH1D* h_EMu_QCD_invm = ((TH1D*)(h_EMu_SS_data_invm->Clone("h_emu_mass_QCD")));
    TH1D* h_EMu_QCD_invm2 = ((TH1D*)(h_EMu_SS_data_invm2->Clone("h_emu_mass_QCD2")));
    h_EMu_QCD_invm->Add(((TH1D*)(s_EMu_SS_invm->GetStack()->Last())), -1);
    h_EMu_QCD_invm2->Add(((TH1D*)(s_EMu_SS_invm2->GetStack()->Last())), -1);
    const double RR = 0.57147108645;
    h_EMu_QCD_invm->Scale(1/RR);
    h_EMu_QCD_invm2->Scale(1/RR);
    removeNegativeBins(h_EMu_QCD_invm);
    removeNegativeBins(h_EMu_QCD_invm2);

    for (Int_t i=1; i<h_EMu_QCD_invm->GetSize()-1; i++)
    {
        if (h_EMu_QCD_invm->GetBinContent(i) < 0)
            h_EMu_QCD_invm->SetBinContent(i, 0);
    }
    for (Int_t i=1; i<h_EMu_QCD_invm2->GetSize()-1; i++)
    {
        if (h_EMu_QCD_invm2->GetBinContent(i) < 0)
            h_EMu_QCD_invm2->SetBinContent(i, 0);
    }

    h_EMu_QCD_invm->SetFillColor(kRed+3);
    h_EMu_QCD_invm->SetLineColor(kRed+3);
    h_EMu_QCD_invm->SetDirectory(0);
    h_EMu_QCD_invm2->SetFillColor(kRed+3);
    h_EMu_QCD_invm2->SetLineColor(kRed+3);
    h_EMu_QCD_invm2->SetDirectory(0);

// ------------------------------------ Draw stack histograms -------------------------------------

    THStack *s_EMu_wQCD_invm = new THStack("s_emu_wQCD_invm", "");
    THStack *s_EMu_wQCD_invm2 = new THStack("s_emu_wQCD_invm2", "");

    if (UseFR)
    {
        s_EMu_wQCD_invm->Add(h_EMu_QCDFR_invm);
        s_EMu_wQCD_invm->Add(h_EMu_WJetsFR_invm);
    }
    else
        s_EMu_wQCD_invm->Add(h_EMu_QCD_invm);
    s_EMu_wQCD_invm2->Add(h_EMu_QCD_invm2);
    for (SelProc_t pr=_EMu_WJets_Full; pr>_EndOf_EMu_ttbar_Normal; pr=SelProc_t((int)(pr-1)))
    {
        if (pr == _EndOf_EMu_VVnST_Normal) continue;
        Mgr.SetProc(pr);
        s_EMu_wQCD_invm2->Add(h_EMu_invm2[pr]);
        if ((isWJ == 0 || UseFR) && pr == _EMu_WJets_Full) {pr = _EndOf_EMu_VVnST_Normal; continue;}
        s_EMu_wQCD_invm->Add(h_EMu_invm[pr]);
        if (pr == _EMu_WJets_Full) pr = _EndOf_EMu_VVnST_Normal; // next -- WW
        if (pr == _EMu_tW) pr = _EMu_VVnST; // next -- ttbar
        if (pr == _EMu_DYTauTau_Full) break;
    }

    myRatioPlot_t* RP_EMu_wQCD_invm = new myRatioPlot_t("EMu_wQCD_mass", s_EMu_wQCD_invm, h_EMu_data_invm);
    myRatioPlot_t* RP_EMu_wQCD_invm2 = new myRatioPlot_t("EMu_wQCD_mass2", s_EMu_wQCD_invm2, h_EMu_data_invm2);

    RP_EMu_wQCD_invm->SetPlots("m_{#lower[-0.2]{#scale[1.2]{#font[12]{e}#mu}}} [GeV/c^{2}]", 15, 3000);
    RP_EMu_wQCD_invm2->SetPlots("m_{#lower[-0.2]{#scale[1.2]{#font[12]{e}#mu}}} [GeV/c^{2}]", 15, 3000);

    TLegend *legend_EMu = new TLegend(0.77, 0.5, 0.95, 0.95);
    legend_EMu->AddEntry(h_EMu_data_invm, "Data", "lp");
    legend_EMu->AddEntry(h_EMu_invm[_EMu_DYTauTau_Full], "DY#rightarrow#tau#tau (MC)","f");
    legend_EMu->AddEntry(h_EMu_invm[_EMu_ttbar_Full], "#kern[0.2]{#font[12]{#scale[1.1]{t#bar{t}}}} (MC)", "f");
    legend_EMu->AddEntry(h_EMu_invm[_EMu_tW], "#kern[0.1]{#font[12]{#scale[1.1]{tW}}} (MC)", "f");
    legend_EMu->AddEntry(h_EMu_invm[_EMu_tbarW], "#kern[0.1]{#font[12]{#scale[1.1]{#bar{t}W}}} (MC)", "f");
    legend_EMu->AddEntry(h_EMu_invm[_EMu_ZZ], "#kern[0.1]{#font[12]{#scale[1.1]{ZZ}}} (MC)", "f");
    legend_EMu->AddEntry(h_EMu_invm[_EMu_WZ], "#font[12]{#scale[1.1]{WZ}} (MC)", "f");
    legend_EMu->AddEntry(h_EMu_invm[_EMu_WW], "#font[12]{#scale[1.1]{WW}} (MC)", "f");
    if (UseFR)
    {
        legend_EMu->AddEntry(h_EMu_WJetsFR_invm, "#font[12]{#scale[1.1]{W}}+Jets (est.)", "f");
        legend_EMu->AddEntry(h_EMu_QCDFR_invm, "#font[12]{#scale[1.1]{QCD}} (est.)", "f");
    }
    else
    {
        if (isWJ) legend_EMu->AddEntry(h_EMu_invm[_EMu_WJets_Full], "#font[12]{#scale[1.1]{W}}+Jets (MC)", "f");
        legend_EMu->AddEntry(h_EMu_QCD_invm, "#font[12]{#scale[1.1]{QCD}} (est.)", "f");
    }

    RP_EMu_wQCD_invm->ImportLegend(legend_EMu);
    RP_EMu_wQCD_invm2->ImportLegend(legend_EMu);

    RP_EMu_wQCD_invm->Draw(4e-1, 7e4, 1);
    RP_EMu_wQCD_invm2->Draw(4e-1, 7e4, 1);

    if (UseFR)
    {
        myRatioPlot_t* RP_EMu_QCD_compare = new myRatioPlot_t("RP_EMu_QCD_compare", h_EMu_QCDFR_invm, h_EMu_QCD_invm);
        RP_EMu_QCD_compare->SetPlots("m_{#lower[-0.2]{#scale[1.2]{#font[12]{e}#mu}}} [GeV/c^{2}]", 15, 3000);
        TLegend *legend_QCD_comp = new TLegend(0.7, 0.7, 0.95, 0.95);
        legend_QCD_comp->AddEntry(h_EMu_QCDFR_invm, "QCD est. (from FR)", "f");
        legend_QCD_comp->AddEntry(h_EMu_QCDFR_invm, "QCD est. (from SS)", "f");
        RP_EMu_QCD_compare->ImportLegend(legend_QCD_comp);
        RP_EMu_QCD_compare->Draw(4e-1, 1e4, 1);
        if (isWJ)
        {
            myRatioPlot_t* RP_EMu_WJets_compare = new myRatioPlot_t("RP_EMu_WJets_compare", h_EMu_WJetsFR_invm, h_EMu_invm[_EMu_WJets_Full]);
            RP_EMu_WJets_compare->SetPlots("m_{#lower[-0.2]{#scale[1.2]{#font[12]{e}#mu}}} [GeV/c^{2}]", 15, 3000);
            TLegend *legend_WJets_comp = new TLegend(0.7, 0.7, 0.95, 0.95);
            legend_WJets_comp->AddEntry(h_EMu_WJetsFR_invm, "W+Jets est. (from FR)", "f");
            legend_WJets_comp->AddEntry(h_EMu_invm[_EMu_WJets_Full], "W+Jets (MC)", "f");
            RP_EMu_WJets_compare->ImportLegend(legend_WJets_comp);
            RP_EMu_WJets_compare->Draw(4e-1, 1e4, 1);
        }
    }

    TH1D *h_EMu_data_invm_old = ((TH1D*)(h_EMu_data_invm->Clone("h_emu_data_invm_old")));
    TH1D *h_EMu_data_invm_old2 = ((TH1D*)(h_EMu_data_invm2->Clone("h_emu_data_invm_old2")));
    if (UseFR)
    {
        h_EMu_data_invm->Add(h_EMu_QCDFR_invm, -1);
        h_EMu_data_invm->Add(h_EMu_WJetsFR_invm, -1);
    }
    else
        h_EMu_data_invm->Add(h_EMu_QCD_invm, -1);
    h_EMu_data_invm2->Add(h_EMu_QCD_invm2, -1);
    removeNegativeBins(h_EMu_QCD_invm);
    removeNegativeBins(h_EMu_QCD_invm2);

    Double_t dataerror_emu, MCerror_emu, dataintegral_emu=2.25081e+07, MCintegral_emu, QCDerror, QCDintegral, QCDFRerror, QCDFRintegral;
    Double_t dataerrorZ_emu, MCerrorZ_emu, dataintegralZ_emu=2.25081e+07, MCintegralZ_emu;
    Double_t dataerror_noZ_emu=0, MCerror_noZ_emu=0, dataintegral_noZ_emu=2.25081e+07, MCintegral_noZ_emu, temp_noZ_emu;
    Double_t WJetsMCerror_emu, WJetsMCintegral_emu, WJetsFRerror_emu, WJetsFRintegral_emu;

    dataintegral_emu = h_EMu_data_invm->IntegralAndError(1, h_EMu_data_invm->GetSize()-2, dataerror_emu);
    MCintegral_emu = ((TH1D*)(s_EMu_invm->GetStack()->Last()))->IntegralAndError(1, h_EMu_data_invm->GetSize()-2, MCerror_emu);

    dataintegralZ_emu = h_EMu_data_invm->IntegralAndError(10, 22, dataerrorZ_emu);
    MCintegralZ_emu = ((TH1D*)(s_EMu_invm->GetStack()->Last()))->IntegralAndError(10, 22, MCerrorZ_emu);

    QCDintegral = h_EMu_QCD_invm->IntegralAndError(1, h_EMu_QCD_invm->GetSize()-2, QCDerror);
    if (UseFR) QCDFRintegral = h_EMu_QCDFR_invm->IntegralAndError(1, h_EMu_QCDFR_invm->GetSize()-2, QCDFRerror);
    WJetsMCintegral_emu = h_EMu_invm[_EMu_WJets_Full]->IntegralAndError(1, h_EMu_invm[_EMu_WJets_Full]->GetSize()-2, WJetsMCerror_emu);
    if (UseFR) WJetsFRintegral_emu = h_EMu_WJetsFR_invm->IntegralAndError(1, h_EMu_WJetsFR_invm->GetSize()-2, WJetsFRerror_emu);

    dataintegral_noZ_emu = h_EMu_data_invm->IntegralAndError(1, 9, temp_noZ_emu);
    dataerror_noZ_emu += temp_noZ_emu * temp_noZ_emu;
    dataintegral_noZ_emu += h_EMu_data_invm->IntegralAndError(23, h_EMu_data_invm->GetSize()-2, temp_noZ_emu);
    dataerror_noZ_emu += temp_noZ_emu * temp_noZ_emu;
    dataerror_noZ_emu = sqrt(dataerror_noZ_emu);

    MCintegral_noZ_emu = ((TH1D*)(s_EMu_invm->GetStack()->Last()))->IntegralAndError(1, 9, temp_noZ_emu);
    MCerror_noZ_emu += temp_noZ_emu * temp_noZ_emu;
    MCintegral_noZ_emu += ((TH1D*)(s_EMu_invm->GetStack()->Last()))->IntegralAndError(23, h_EMu_data_invm->GetSize()-2, temp_noZ_emu);
    MCerror_noZ_emu += temp_noZ_emu * temp_noZ_emu;
    MCerror_noZ_emu = sqrt(MCerror_noZ_emu);

    std::cout << "EMu data events: " << dataintegral_emu << "+-" << dataerror_emu << endl;
    std::cout << "EMu MC events: " << MCintegral_emu << "+-" << MCerror_emu << endl;
    std::cout << "EMu MC/Obs: " << MCintegral_emu / dataintegral_emu << "+-" <<
                 sqrt((dataerror_emu / dataintegral_emu) * (dataerror_emu / dataintegral_emu) +
                       (MCerror_emu / MCintegral_emu) * (MCerror_emu / MCintegral_emu)) << endl;

    std::cout << "EMu data events around Z: " << dataintegralZ_emu << "+-" << dataerrorZ_emu << endl;
    std::cout << "EMu MC events around Z: " << MCintegralZ_emu << "+-" << MCerrorZ_emu << endl;
    std::cout << "EMu data events outside Z: " << dataintegral_noZ_emu << "+-" << dataerror_noZ_emu << endl;
    std::cout << "EMu MC events outside Z: " << MCintegral_noZ_emu << "+-" << MCerror_noZ_emu << endl;
    std::cout << "EMu QCD events (from SS): " << QCDintegral << "+-" << QCDerror << endl;
    if (UseFR) std::cout << "EMu QCD events (FR): " << QCDFRintegral << "+-" << QCDFRerror << endl;
    std::cout << "EMu WJets events (MC): " << WJetsMCintegral_emu << "+-" << WJetsMCerror_emu << endl;
    if (UseFR) std::cout << "EMu WJets events (FR): " << WJetsFRintegral_emu << "+-" << WJetsFRerror_emu << endl;

//************************************ E E **********************************************************

    TH1D *h_ee_invm[_EndOf_EMu], *h_ee_invm2[_EndOf_EMu];
    THStack* s_ee_invm = new THStack("s_ee_invm", "");
    THStack* s_ee_invm2 = new THStack("s_ee_invm2", "");

    for (SelProc_t pr=_EE_WW; pr<=_EE_VVnST; pr=SelProc_t((int)(pr-1)))
    {
        if (pr == _EE_WZ || pr == _EE_ZZ) continue;
        Mgr.SetProc(pr);

        f_Bkg_ee->GetObject("h_mass_"+Mgr.Procname[pr], h_ee_invm[pr]);
        f_Bkg_ee->GetObject("h_mass2_"+Mgr.Procname[pr], h_ee_invm2[pr]);
        removeNegativeBins(h_ee_invm[pr]);
        removeNegativeBins(h_ee_invm2[pr]);

        Color_t color = kBlack;
        if (pr == _EE_WW) color = kMagenta - 5;
        if (pr == _EE_WZ) color = kMagenta - 2;
        if (pr == _EE_ZZ) color = kMagenta - 6;
        if (pr == _EE_tbarW) color = kGreen - 2;
        if (pr == _EE_tW) color = kGreen + 2;
        if (pr == _EE_ttbar_Full) color = kCyan + 2;
        if (pr == _EE_DYTauTau_Full) color = kOrange - 5;

        h_ee_invm[pr]->SetFillColor(color);
        h_ee_invm[pr]->SetLineColor(color);
        h_ee_invm[pr]->SetDirectory(0);
        s_ee_invm->Add(h_ee_invm[pr]);

        h_ee_invm2[pr]->SetFillColor(color);
        h_ee_invm2[pr]->SetLineColor(color);
        h_ee_invm2[pr]->SetDirectory(0);
        s_ee_invm2->Add(h_ee_invm2[pr]);

        if (pr == _EE_tW) pr = _EE_VVnST; // next -- ttbar
        if (pr == _EE_DYTauTau_Full) break;
    }


//####################################### ALL AT ONCE ###############################################

    TH1D* h_ee_Est_invm = ((TH1D*)(h_EMu_data_invm->Clone("h_ee_mass_Est")));
    h_ee_Est_invm->Multiply(((TH1D*)(s_ee_invm->GetStack()->Last())));
    h_ee_Est_invm->Divide(((TH1D*)(s_EMu_invm->GetStack()->Last())));
    removeNegativeBins(h_ee_Est_invm);
    h_ee_Est_invm->SetDirectory(0);
    h_ee_Est_invm->SetLineColor(kBlack);
    h_ee_Est_invm->SetMarkerStyle(kFullTriangleUp);
    h_ee_Est_invm->SetMarkerSize(1.5);
    h_ee_Est_invm->SetMarkerColor(kBlack);

    TH1D* h_ee_Est_invm2 = ((TH1D*)(h_EMu_data_invm2->Clone("h_ee_mass_Est2")));
    h_ee_Est_invm2->Multiply(((TH1D*)(s_ee_invm2->GetStack()->Last())));
    h_ee_Est_invm2->Divide(((TH1D*)(s_EMu_invm2->GetStack()->Last())));
    removeNegativeBins(h_ee_Est_invm2);
    h_ee_Est_invm2->SetDirectory(0);
    h_ee_Est_invm2->SetLineColor(kBlack);
    h_ee_Est_invm2->SetMarkerStyle(kFullTriangleUp);
    h_ee_Est_invm2->SetMarkerSize(1.5);
    h_ee_Est_invm2->SetMarkerColor(kBlack);

/// ############################# Statistical uncertainties (by hand) ##################################### ///

    Double_t UNC_eeMC[43], UNC_emuMC[43], UNC_emuData[43], UNC_eeEst[43], UNC_ee_EstOverMC[43],
             UNC_eeMC2[86], UNC_emuMC2[86], UNC_emuData2[86], UNC_eeEst2[86], UNC_ee_EstOverMC2[86];
    Double_t tmp = 0;
    for (int i=1; i<87; i++)
    {
        if (i < 44)
        {
            // emuMC (N2)
            UNC_emuMC[i-1] = ((TH1D*)(s_EMu_invm->GetStack()->Last()))->GetBinError(i) * ((TH1D*)(s_EMu_invm->GetStack()->Last()))->GetBinError(i);
            tmp = 0;
            for (SelProc_t pr=_EMu_WJets_Full; pr>_EndOf_EMu_ttbar_Normal; pr=SelProc_t((int)(pr-1)))
            {
                Mgr.SetProc(pr);
                if (pr==_EndOf_EMu_VVnST_Normal) continue;
                if (isWJ == 0 && pr == _EMu_WJets_Full) {pr = _EndOf_EMu_VVnST_Normal; continue;}

                tmp += h_EMu_invm[pr]->GetBinError(i) * h_EMu_invm[pr]->GetBinError(i);

                if (pr == _EMu_WJets_Full) // next - WW
                    pr = _EndOf_EMu_VVnST_Normal;
                if (pr == _EMu_tW) pr = _EMu_VVnST; // next -- ttbar
                if (pr == _EMu_DYTauTau_Full) break;
            }
            if ((UNC_emuMC[i-1]-tmp)/tmp > 0.0001)
            {
                cout << "EMU errors from THStack and TH1D sum do not match!!" << endl;
                cout << "Bin " << i << endl;
                cout << "From THStack: " << UNC_emuMC[i-1] << "   From TH1D sum: " << tmp << endl;
                break;
            }

            // eeMC (N1)
            UNC_eeMC[i-1] = ((TH1D*)(s_ee_invm->GetStack()->Last()))->GetBinError(i) * ((TH1D*)(s_ee_invm->GetStack()->Last()))->GetBinError(i);
            tmp = 0;
            for (SelProc_t pr=_EE_WW; pr<=_EE_VVnST; pr=SelProc_t((int)(pr-1)))
            {
                if (pr == _EE_WZ || pr == _EE_ZZ) continue;
                Mgr.SetProc(pr);

                tmp += h_ee_invm[pr]->GetBinError(i) * h_ee_invm[pr]->GetBinError(i);

                if (pr == _EE_tW) pr = _EE_VVnST; // next -- ttbar
                if (pr == _EE_DYTauTau_Full) break;
            }
            if ((UNC_eeMC[i-1]-tmp)/tmp > 0.0001)
            {
                cout << "EE errors from THStack and TH1D sum do not match!!" << endl;
                cout << "Bin " << i << endl;
                cout << "From THStack: " << UNC_eeMC[i-1] << "   From TH1D sum: " << tmp << endl;
                break;
            }

            // emuData (N3)
            UNC_emuData[i-1] = ((TH1D*)(s_EMu_SS_invm->GetStack()->Last()))->GetBinError(i) * ((TH1D*)(s_EMu_SS_invm->GetStack()->Last()))->GetBinError(i);
            tmp = 0;
            for (SelProc_t pr=_EMu_WJets_Full; pr>_EndOf_EMu_ttbar_Normal; pr=SelProc_t((int)(pr-1)))
            {
                Mgr.SetProc(pr);
                if (pr==_EndOf_EMu_VVnST_Normal) continue;
                if (isWJ == 0 && pr == _EMu_WJets_Full) {pr = _EndOf_EMu_VVnST_Normal; continue;}

                tmp += h_EMu_SS_invm[pr]->GetBinError(i) * h_EMu_SS_invm[pr]->GetBinError(i);

                if (pr == _EMu_WJets_Full) // next - WW
                    pr = _EndOf_EMu_VVnST_Normal;
                if (pr == _EMu_tW) pr = _EMu_VVnST; // next -- ttbar
                if (pr == _EMu_DYTauTau_Full) break;
            }
            if ((UNC_emuData[i-1]-tmp)/tmp > 0.0001)
            {
                cout << "EMU SS errors from THStack and TH1D sum do not match!!" << endl;
                cout << "Bin " << i << endl;
                cout << "From THStack: " << UNC_emuData[i-1] << "   From TH1D sum: " << tmp << endl;
                break;
            }
            UNC_emuData[i-1] *= 1/(RR*RR);
            UNC_emuData[i-1] += h_EMu_SS_data_invm->GetBinError(i) * h_EMu_SS_data_invm->GetBinError(i) / (RR*RR);
            UNC_emuData[i-1] += h_EMu_data_invm_old->GetBinError(i) * h_EMu_data_invm_old->GetBinError(i);
            tmp = h_EMu_data_invm->GetBinError(i) * h_EMu_data_invm->GetBinError(i);
            if ((UNC_emuData[i-1]-tmp)/tmp > 0.0001)
            {
                cout << "EMU DATA errors from automatic and manual counting do not match!!" << endl;
                cout << "Bin " << i << endl;
                cout << "From automatic: " << UNC_emuData[i-1] << "   From manual: " << tmp << endl;
                break;
            }

            // FULL UNC
            tmp = UNC_eeMC[i-1] / (((TH1D*)(s_ee_invm->GetStack()->Last()))->GetBinContent(i) * ((TH1D*)(s_ee_invm->GetStack()->Last()))->GetBinContent(i));
            tmp += UNC_emuMC[i-1] / (((TH1D*)(s_EMu_invm->GetStack()->Last()))->GetBinContent(i) * ((TH1D*)(s_EMu_invm->GetStack()->Last()))->GetBinContent(i));
            tmp += UNC_emuData[i-1] / (h_EMu_data_invm->GetBinContent(i) * h_EMu_data_invm->GetBinContent(i));
            UNC_eeEst[i-1] = sqrt(tmp) * fabs(h_ee_Est_invm->GetBinContent(i));
            tmp = h_ee_Est_invm->GetBinError(i) * h_ee_Est_invm->GetBinError(i);
            if ((UNC_eeEst[i-1]-tmp)/tmp > 0.0001)
            {
                cout << "EE EST errors from automatic and manual counting do not match!!" << endl;
                cout << "Bin " << i << endl;
                cout << "From automatic: " << UNC_eeEst[i-1] << "   From manual: " << tmp << endl;
                break;
            }

        } // End of if(i<44)

        // FINER BINNING
        // emuMC (N2)
        UNC_emuMC2[i-1] = ((TH1D*)(s_EMu_invm2->GetStack()->Last()))->GetBinError(i) * ((TH1D*)(s_EMu_invm2->GetStack()->Last()))->GetBinError(i);
        tmp = 0;
        for (SelProc_t pr=_EMu_WJets_Full; pr>_EndOf_EMu_ttbar_Normal; pr=SelProc_t((int)(pr-1)))
        {
            Mgr.SetProc(pr);
            if (pr==_EndOf_EMu_VVnST_Normal) continue;
            if (isWJ == 0 && pr == _EMu_WJets_Full) {pr = _EndOf_EMu_VVnST_Normal; continue;}

            tmp += h_EMu_invm2[pr]->GetBinError(i) * h_EMu_invm2[pr]->GetBinError(i);

            if (pr == _EMu_WJets_Full) // next - WW
                pr = _EndOf_EMu_VVnST_Normal;
            if (pr == _EMu_tW) pr = _EMu_VVnST; // next -- ttbar
            if (pr == _EMu_DYTauTau_Full) break;
        }
        if ((UNC_emuMC2[i-1]-tmp)/tmp > 0.0001)
        {
            cout << "EMU2 errors from THStack and TH1D sum do not match!!" << endl;
            cout << "Bin " << i << endl;
            cout << "From THStack: " << UNC_emuMC2[i-1] << "   From TH1D sum: " << tmp << endl;
            break;
        }

        // eeMC (N1)
        UNC_eeMC2[i-1] = ((TH1D*)(s_ee_invm2->GetStack()->Last()))->GetBinError(i) * ((TH1D*)(s_ee_invm2->GetStack()->Last()))->GetBinError(i);
        tmp = 0;
        for (SelProc_t pr=_EE_WW; pr<=_EE_VVnST; pr=SelProc_t((int)(pr-1)))
        {
            if (pr == _EE_WZ || pr == _EE_ZZ) continue;
            Mgr.SetProc(pr);

            tmp += h_ee_invm2[pr]->GetBinError(i) * h_ee_invm2[pr]->GetBinError(i);

            if (pr == _EE_tW) pr = _EE_VVnST; // next -- ttbar
            if (pr == _EE_DYTauTau_Full) break;
        }
        if ((UNC_eeMC2[i-1]-tmp)/tmp > 0.0001)
        {
            cout << "EE2 errors from THStack and TH1D sum do not match!!" << endl;
            cout << "Bin " << i << endl;
            cout << "From THStack: " << UNC_eeMC2[i-1] << "   From TH1D sum: " << tmp << endl;
            break;
        }

        // emuData (N3)
        UNC_emuData2[i-1] = ((TH1D*)(s_EMu_SS_invm2->GetStack()->Last()))->GetBinError(i) * ((TH1D*)(s_EMu_SS_invm2->GetStack()->Last()))->GetBinError(i);
        tmp = 0;
        for (SelProc_t pr=_EMu_WJets_Full; pr>_EndOf_EMu_ttbar_Normal; pr=SelProc_t((int)(pr-1)))
        {
            Mgr.SetProc(pr);
            if (pr==_EndOf_EMu_VVnST_Normal) continue;
            if (isWJ == 0 && pr == _EMu_WJets_Full) {pr = _EndOf_EMu_VVnST_Normal; continue;}

            tmp += h_EMu_SS_invm2[pr]->GetBinError(i) * h_EMu_SS_invm2[pr]->GetBinError(i);

            if (pr == _EMu_WJets_Full) // next - WW
                pr = _EndOf_EMu_VVnST_Normal;
            if (pr == _EMu_tW) pr = _EMu_VVnST; // next -- ttbar
            if (pr == _EMu_DYTauTau_Full) break;
        }
        if ((UNC_emuData2[i-1]-tmp)/tmp > 0.0001)
        {
            cout << "EMU2 SS errors from THStack and TH1D sum do not match!!" << endl;
            cout << "Bin " << i << endl;
            cout << "From THStack: " << UNC_emuData2[i-1] << "   From TH1D sum: " << tmp << endl;
            break;
        }
        UNC_emuData2[i-1] *= 1/(RR*RR);
        UNC_emuData2[i-1] += h_EMu_SS_data_invm2->GetBinError(i) * h_EMu_SS_data_invm2->GetBinError(i) / (RR*RR);
        UNC_emuData2[i-1] += h_EMu_data_invm_old2->GetBinError(i) * h_EMu_data_invm_old2->GetBinError(i);
        tmp = h_EMu_data_invm2->GetBinError(i) * h_EMu_data_invm2->GetBinError(i);
        if ((UNC_emuData2[i-1]-tmp)/tmp > 0.0001)
        {
            cout << "EMU2 DATA errors from automatic and manual counting do not match!!" << endl;
            cout << "Bin " << i << endl;
            cout << "From automatic: " << UNC_emuData2[i-1] << "   From manual: " << tmp << endl;
            break;
        }

        // FULL UNC
        tmp = UNC_eeMC2[i-1] / (((TH1D*)(s_ee_invm2->GetStack()->Last()))->GetBinContent(i) * ((TH1D*)(s_ee_invm2->GetStack()->Last()))->GetBinContent(i));
        tmp += UNC_emuMC2[i-1] / (((TH1D*)(s_EMu_invm2->GetStack()->Last()))->GetBinContent(i) * ((TH1D*)(s_EMu_invm2->GetStack()->Last()))->GetBinContent(i));
        tmp += UNC_emuData2[i-1] / (h_EMu_data_invm2->GetBinContent(i) * h_EMu_data_invm2->GetBinContent(i));
        UNC_eeEst2[i-1] = sqrt(tmp) * fabs(h_ee_Est_invm2->GetBinContent(i));
        tmp = h_ee_Est_invm2->GetBinError(i) * h_ee_Est_invm2->GetBinError(i);
        if ((UNC_eeEst2[i-1]-tmp)/tmp > 0.0001)
        {
            cout << "EE2 EST errors from automatic and manual counting do not match!!" << endl;
            cout << "Bin " << i << endl;
            cout << "From automatic: " << UNC_eeEst2[i-1] << "   From manual: " << tmp << endl;
            break;
        }
    } // End of for(bins)
//    cout << "Manually calculated errors:" << endl;
//    for (int i=0; i<43; i++)
//    {
//        cout << UNC_eeEst[i] << "  ";
//    }
//    cout << "\nAutomatically calculated errors:" << endl;
//    for (int i=1; i<=43; i++)
//    {
//        cout << h_ee_Est_invm->GetBinError(i) << "  ";
//    }
//    cout << endl;

    // DATA/MC UNC (probably smaller than given automatically)
    for (int i=0; i<86; i++)
    {
        if (i < 43)
        {
            UNC_ee_EstOverMC[i] = ((TH1D*)(s_EMu_invm->GetStack()->Last()))->GetBinError(i+1) / ((TH1D*)(s_EMu_invm->GetStack()->Last()))->GetBinContent(i+1); // N2
            UNC_ee_EstOverMC[i] *= UNC_ee_EstOverMC[i];
            tmp = h_EMu_data_invm->GetBinError(i+1) / h_EMu_data_invm->GetBinContent(i+1); // N3
            tmp *= tmp;
            UNC_ee_EstOverMC[i] = sqrt(UNC_ee_EstOverMC[i] + tmp);
        }
        UNC_ee_EstOverMC2[i] = ((TH1D*)(s_EMu_invm2->GetStack()->Last()))->GetBinError(i+1) / ((TH1D*)(s_EMu_invm2->GetStack()->Last()))->GetBinContent(i+1); // N2
        UNC_ee_EstOverMC2[i] *= UNC_ee_EstOverMC2[i];
        tmp = h_EMu_data_invm2->GetBinError(i+1) / h_EMu_data_invm2->GetBinContent(i+1); // N3
        tmp *= tmp;
        UNC_ee_EstOverMC2[i] = sqrt(UNC_ee_EstOverMC2[i] + tmp);
    }

// ############################# Systematics (DataDriven - MC) ##################################### //

    Double_t systematics[h_ee_Est_invm->GetSize()-2], systematics2[h_ee_Est_invm2->GetSize()-2], systerr=0,
             estSystematics[h_ee_Est_invm->GetSize()-2], estSystematics2[h_ee_Est_invm2->GetSize()-2];
    for (int i=1; i<h_ee_Est_invm2->GetSize()-1; i++)
    {
        if (i < h_ee_Est_invm->GetSize()-1)
        {
            estSystematics[i-1] = fabs(h_ee_Est_invm->GetBinContent(i) - ((TH1D*)(s_ee_invm->GetStack()->Last()))->GetBinContent(i));
            systematics[i-1] = estSystematics[i-1] / fabs(((TH1D*)(s_ee_invm->GetStack()->Last()))->GetBinContent(i));
            systerr += estSystematics[i-1] * estSystematics[i-1];
        }
        estSystematics2[i-1] = fabs(h_ee_Est_invm2->GetBinContent(i) - ((TH1D*)(s_ee_invm2->GetStack()->Last()))->GetBinContent(i));
        systematics2[i-1] = estSystematics2[i-1] / fabs(((TH1D*)(s_ee_invm2->GetStack()->Last()))->GetBinContent(i));
    }

// ###################################### Writing ###################################### //

    TFile* f_NeeEst = new TFile(Mgr.HistLocation+"EstBkg_EE.root", "RECREATE");
    if (f_NeeEst->IsOpen()) std::cout << "File EstBkg_EE.root opened successfully\n";
    else { std::cout << "File EstBkg_EE.root was not opened\n"; return; }
    f_NeeEst->cd();
    h_ee_Est_invm->Write();
    h_ee_Est_invm2->Write();
    TVectorD estSyst(43, estSystematics);
    TVectorD estSyst2(86, estSystematics);

    estSyst.Write("estSystematics");
    estSyst2.Write("estSystematics2");

    Double_t MCerror, esterror, MCintegral, estintegral;
    Double_t MCerrorZ, esterrorZ, MCintegralZ, estintegralZ;
    Double_t MCerror_noZ=0, esterror_noZ=0, MCintegral_noZ, estintegral_noZ, temp_noZ;

    MCintegral = ((TH1D*)(s_ee_invm->GetStack()->Last()))->IntegralAndError(1, h_EMu_data_invm->GetSize()-2, MCerror);
    estintegral = h_ee_Est_invm->IntegralAndError(1, h_ee_Est_invm->GetSize()-2, esterror);

    MCintegralZ = ((TH1D*)(s_ee_invm->GetStack()->Last()))->IntegralAndError(10, 22, MCerrorZ);
    estintegralZ = h_ee_Est_invm->IntegralAndError(10, 22, esterrorZ);

    MCintegral_noZ = ((TH1D*)(s_ee_invm->GetStack()->Last()))->IntegralAndError(1, 9, temp_noZ);
    MCerror_noZ += temp_noZ * temp_noZ;
    MCintegral_noZ += ((TH1D*)(s_ee_invm->GetStack()->Last()))->IntegralAndError(23, h_EMu_data_invm->GetSize()-2, temp_noZ);
    MCerror_noZ += temp_noZ * temp_noZ;
    MCerror_noZ = sqrt(MCerror_noZ);

    estintegral_noZ = h_ee_Est_invm->IntegralAndError(1, 9, temp_noZ);
    esterror_noZ += temp_noZ * temp_noZ;
    estintegral_noZ += h_ee_Est_invm->IntegralAndError(23, h_ee_Est_invm->GetSize()-2, temp_noZ);
    esterror_noZ += temp_noZ * temp_noZ;
    esterror_noZ = sqrt(esterror_noZ);

    systerr = sqrt(systerr);

    std::cout << "ee MC events: " << MCintegral << "+-" << MCerror << endl;
    std::cout << "ee Est. events: " << estintegral << "+-" << esterror << "+-" << systerr << endl << endl;

    std::cout << "ee MC events around Z: " << MCintegralZ << "+-" << MCerrorZ << endl;
    std::cout << "ee Est. events around Z: " << estintegralZ << "+-" << esterrorZ << endl;
    std::cout << "ee MC events outside Z: " << MCintegral_noZ << "+-" << MCerror_noZ << endl;
    std::cout << "ee Est. events outside Z: " << estintegral_noZ << "+-" << esterror_noZ << endl << endl;

//------------------------- Z PEAK ------------------------------------------------------------------------

    Double_t MCZerr;
    Double_t MCZint = ((TH1D*)(s_ee_invm->GetStack()->Last()))->IntegralAndError(10, 23, MCZerr);
    Double_t EMuZerr;
    Double_t EMuZint = h_ee_Est_invm->IntegralAndError(10, 23, EMuZerr);
    std::cout << "MC events around Z: " << MCZint << " +- " << MCZerr << endl;
    std::cout << "EMuEst events around Z: " << EMuZint << " +- " << EMuZerr << endl;

//------------------------------ Drawing ----------------------------------------------------------

    myRatioPlot_t *RP_invm = new myRatioPlot_t("DataDriven_InvariantMass", s_ee_invm, h_ee_Est_invm);
    myRatioPlot_t *RP_invm2 = new myRatioPlot_t("DataDriven_InvariantMass2", s_ee_invm2, h_ee_Est_invm2);
    RP_invm->SetPlots("m_{#lower[-0.2]{#scale[1.2]{#font[12]{ee}}}} [GeV/c^{2}]", 15, 3000, "Est./MC", UNC_ee_EstOverMC);
    RP_invm2->SetPlots("m_{#lower[-0.2]{#scale[1.2]{#font[12]{ee}}}} [GeV/c^{2}]", 15, 3000, "Est./MC", UNC_ee_EstOverMC2);
    RP_invm->SetSystematics(estSystematics, NULL, systematics);
    RP_invm2->SetSystematics(estSystematics2, NULL, systematics2);

    TLegend *legend = new TLegend(0.77, 0.5, 0.95, 0.95);
    legend->AddEntry(h_ee_Est_invm, "Estimation", "lp");
    legend->AddEntry(h_ee_invm[_EE_DYTauTau_Full], "DY#rightarrow#tau#tau","f");
    legend->AddEntry(h_ee_invm[_EE_ttbar_Full], "#kern[0.2]{#font[12]{#scale[1.1]{t#bar{t}}}}", "f");
    legend->AddEntry(h_ee_invm[_EE_tW], "#kern[0.1]{#font[12]{#scale[1.1]{tW}}}", "f");
    legend->AddEntry(h_ee_invm[_EE_tbarW], "#kern[0.1]{#font[12]{#scale[1.1]{#bar{t}W}}}", "f");
    legend->AddEntry(h_ee_invm[_EE_WW], "#font[12]{#scale[1.1]{WW}}", "f");

    RP_invm->ImportLegend(legend);
    RP_invm2->ImportLegend(legend);
    RP_invm->Draw(4e-1, 2e4, 1);
    RP_invm2->Draw(4e-1, 2e4, 1);

//######################################## ONE BY ONE #################################################

    TH1D *h_ee_Est_pr_invm[_EndOf_EMu], *h_ee_Est_pr_invm2[_EndOf_EMu];
    myRatioPlot_t *RP_ee_invm_pr[_EndOf_EMu], *RP_ee_invm_pr2[_EndOf_EMu];
    Double_t emuMCInt[_EndOf_EMu], eeMCInt[_EndOf_EMu], nEEestInt[_EndOf_EMu], emuMCEr[_EndOf_EMu], eeMCEr[_EndOf_EMu], nEEestEr[_EndOf_EMu];

    for (SelProc_t pr=_EE_tW; pr<_EE_VVnST; pr=next(pr))
    {
        if (pr>=_EndOf_EE_VVnST_Normal && pr<_EE_DYTauTau_Full) continue;
        if (pr == _EE_WZ || pr == _EE_ZZ) continue;

        Mgr.SetProc(pr);
        h_ee_Est_pr_invm[pr] = ((TH1D*)(h_ee_invm[pr]->Clone("h_ee_mass_Est_"+Mgr.Procname[pr])));
        removeNegativeBins(h_ee_Est_pr_invm[pr]);
        h_ee_Est_pr_invm[pr]->Multiply(h_EMu_data_invm);
        removeNegativeBins(h_ee_Est_pr_invm[pr]);
        h_ee_Est_pr_invm[pr]->Divide((TH1D*)(s_EMu_invm->GetStack()->Last()));
        removeNegativeBins(h_ee_Est_pr_invm[pr]);

        h_ee_Est_pr_invm[pr]->SetLineColor(kBlack);
        h_ee_Est_pr_invm[pr]->SetDirectory(0);
        h_ee_Est_pr_invm[pr]->Write();

        h_ee_Est_pr_invm2[pr] = ((TH1D*)(h_ee_invm2[pr]->Clone("h_ee_mass_Est2_"+Mgr.Procname[pr])));
        removeNegativeBins(h_ee_Est_pr_invm2[pr]);
        h_ee_Est_pr_invm2[pr]->Multiply(h_EMu_data_invm2);
        removeNegativeBins(h_ee_Est_pr_invm2[pr]);
        h_ee_Est_pr_invm2[pr]->Divide((TH1D*)(s_EMu_invm2->GetStack()->Last()));
        removeNegativeBins(h_ee_Est_pr_invm2[pr]);

        h_ee_Est_pr_invm2[pr]->SetLineColor(kBlack);
        h_ee_Est_pr_invm2[pr]->SetDirectory(0);
        h_ee_Est_pr_invm2[pr]->Write();

        emuMCEr[pr] = 0;
        if (pr <= _EE_WW) // Matching EMu with EE processes
            emuMCInt[pr] = h_EMu_invm[pr+54]->IntegralAndError(1, h_EMu_invm[pr+54]->GetSize()-1, emuMCEr[pr]);
        else
            emuMCInt[pr] = h_EMu_invm[pr+36]->IntegralAndError(1, h_EMu_invm[pr+36]->GetSize()-1, emuMCEr[pr]);
        eeMCEr[pr] = 0;
        eeMCInt[pr] = h_ee_invm[pr]->IntegralAndError(1, h_ee_invm[pr]->GetSize()-1, eeMCEr[pr]);
        nEEestEr[pr] = 0;
        nEEestInt[pr] = h_ee_Est_pr_invm[pr]->IntegralAndError(1, h_ee_Est_pr_invm[pr]->GetSize()-1, nEEestEr[pr]);

        std::cout << Mgr.Procname[pr] << " MC eMu events: " << emuMCInt[pr] << " +- " << emuMCEr[pr] << " (Stat.)\n";
        std::cout << Mgr.Procname[pr] << " MC ee events: " << eeMCInt[pr] << " +- " << eeMCEr[pr] << " (Stat.)\n";
        std::cout << Mgr.Procname[pr] << " ee events estimated from data: " << nEEestInt[pr] << " +- " << nEEestEr[pr] << " (Stat.)\n\n";

//------------------------------ Drawing ----------------------------------------------------------

        RP_ee_invm_pr[pr] = new myRatioPlot_t("h_mass_DataDriven_"+Mgr.Procname[pr], h_ee_invm[pr], h_ee_Est_pr_invm[pr]);
        RP_ee_invm_pr2[pr] = new myRatioPlot_t("h_mass_DataDriven2_"+Mgr.Procname[pr], h_ee_invm2[pr], h_ee_Est_pr_invm2[pr]);
        RP_ee_invm_pr[pr]->SetPlots("m_{#lower[-0.2]{#scale[1.2]{#font[12]{ee}}}} [GeV/c^{2}]", 15, 3000, "Est./MC", UNC_ee_EstOverMC);
        RP_ee_invm_pr2[pr]->SetPlots("m_{#lower[-0.2]{#scale[1.2]{#font[12]{ee}}}} [GeV/c^{2}]", 15, 3000, "Est./MC", UNC_ee_EstOverMC2);
        RP_ee_invm_pr[pr]->SetSystematics(NULL, NULL, systematics);
        RP_ee_invm_pr2[pr]->SetSystematics(NULL, NULL, systematics2);

        TLegend *legend_pr = new TLegend(0.65, 0.75, 0.95, 0.95);

        legend_pr->AddEntry(h_ee_Est_pr_invm[pr], "Estimation", "lp");
//        RP_ee_invm_pr[pr]->AddLegendEntry(h_ee_Est_pr_invm[pr], "Ivertis", "lp");

        if (pr==_EE_tW)
            legend_pr->AddEntry(h_ee_invm[pr], "#kern[0.1]{#font[12]{#scale[1.1]{tW}}} (MC)", "f");
        else if (pr==_EE_tbarW)
            legend_pr->AddEntry(h_ee_invm[pr], "#kern[0.1]{#font[12]{#scale[1.1]{#bar{t}W}}} (MC)", "f");
        else if (pr==_EE_ZZ)
            legend_pr->AddEntry(h_ee_invm[pr], "#kern[0.1]{#font[12]{#scale[1.1]{ZZ}}} (MC)", "f");
        else if (pr==_EE_WZ)
            legend_pr->AddEntry(h_ee_invm[pr], "#font[12]{#scale[1.1]{WZ}} (MC)", "f");
        else if (pr==_EE_WW)
            legend_pr->AddEntry(h_ee_invm[pr], "#font[12]{#scale[1.1]{WW}} (MC)", "f");
        else if (pr==_EE_ttbar_Full)
            legend_pr->AddEntry(h_ee_invm[pr], "#kern[0.2]{#font[12]{#scale[1.1]{t#bar{t}}}} (MC)", "f");
        else if (pr==_EE_DYTauTau_Full)
            legend_pr->AddEntry(h_ee_invm[pr], "DY#rightarrow#tau#tau (MC)", "f");
        else
            legend_pr->AddEntry(h_ee_invm[pr], Mgr.Procname[pr]+" (MC)", "f");

        RP_ee_invm_pr[pr]->ImportLegend(legend_pr);
        RP_ee_invm_pr2[pr]->ImportLegend(legend_pr);
        RP_ee_invm_pr[pr]->Draw(4e-1, 1e4, 1);
        RP_ee_invm_pr2[pr]->Draw(4e-1, 1e4, 1);
    }

//################################## tW + tbarW #######################################################

    THStack* s_tW_invm = new THStack("s_tW", "");
    THStack* s_tW_invm2 = new THStack("s_tW2", "");
    s_tW_invm->Add(h_ee_invm[_EE_tbarW]);
    s_tW_invm->Add(h_ee_invm[_EE_tW]);
    s_tW_invm2->Add(h_ee_invm2[_EE_tbarW]);
    s_tW_invm2->Add(h_ee_invm2[_EE_tW]);

    TH1D* h1_eeTW_invm_Est = ((TH1D*)(h_ee_Est_pr_invm[_EE_tW]->Clone("h_mass_eeTW+TbarW_Est")));
    TH1D* h1_eeTW_invm_Est2 = ((TH1D*)(h_ee_Est_pr_invm2[_EE_tW]->Clone("h_mass_eeTW+TbarW_Est2")));
    h1_eeTW_invm_Est->Add(h_ee_Est_pr_invm[_EE_tbarW]);
    h1_eeTW_invm_Est->SetDirectory(0);
    h1_eeTW_invm_Est->Write();
    h1_eeTW_invm_Est2->Add(h_ee_Est_pr_invm2[_EE_tbarW]);
    h1_eeTW_invm_Est2->SetDirectory(0);
    h1_eeTW_invm_Est2->Write();

    myRatioPlot_t* RP_tW = new myRatioPlot_t("tWest", s_tW_invm, h1_eeTW_invm_Est);
    myRatioPlot_t* RP_tW2 = new myRatioPlot_t("tWest2", s_tW_invm2, h1_eeTW_invm_Est2);
    RP_tW->SetPlots("m_{#lower[-0.2]{#scale[1.2]{#font[12]{ee}}}} [GeV/c^{2}]", 15, 3000, "Est./MC", UNC_ee_EstOverMC);
    RP_tW2->SetPlots("m_{#lower[-0.2]{#scale[1.2]{#font[12]{ee}}}} [GeV/c^{2}]", 15, 3000, "Est./MC", UNC_ee_EstOverMC2);
    RP_tW->SetSystematics(NULL, NULL, systematics);
    RP_tW2->SetSystematics(NULL, NULL, systematics2);

    TLegend *legend_tw = new TLegend(0.65, 0.65, 0.95, 0.95);
    legend_tw->AddEntry(h1_eeTW_invm_Est, "Estimation", "lp");
    legend_tw->AddEntry(h_ee_invm[_EE_tW], "#kern[0.1]{#font[12]{#scale[1.1]{tW}}} (MC)", "f");
    legend_tw->AddEntry(h_ee_invm[_EE_tbarW], "#kern[0.1]{#font[12]{#scale[1.1]{#bar{t}W}}} (MC)", "f");

    RP_tW->ImportLegend(legend_tw);
    RP_tW2->ImportLegend(legend_tw);
    RP_tW->Draw(4e-1, 1e4, 1);
    RP_tW2->Draw(4e-1, 1e4, 1);

//----------------------------------------- Closing up -----------------------------------------------

    f_Bkg_EMu->Close();
    if (!f_Bkg_EMu->IsOpen()) std::cout << "File Hist_"+Mgr.Procname[_EMu_Bkg_Full]+".root closed successfully\n";
    f_Data_EMu->Close();
    if (!f_Data_EMu->IsOpen()) std::cout << "File Hist_"+Mgr.Procname[_EMu_SingleMuon_Full]+".root closed successfully\n";
    f_Bkg_ee->Close();
    if (!f_Bkg_ee->IsOpen()) std::cout << "File Hist_"+Mgr.Procname[_EE_Bkg_Full]+".root closed successfully\n";
    f_DY_ee->Close();
    if(!f_DY_ee->IsOpen()) std::cout << "File Hist_"+Mgr.Procname[_EE_DY_Full]+".root closed successfully.\n";
    f_Data_ee->Close();
    if (!f_Data_ee->IsOpen()) std::cout << "File Hist_"+Mgr.Procname[_EE_DoubleEG_Full]+".root closed successfully\n";
    f_NeeEst->Close();
    if (!f_NeeEst->IsOpen()) std::cout << "File EstBkg_EE.root closed successfully\n\n";

}// End of ee_Est


void MuMu_Est(Int_t UseFR)
{
    LocalFileMgr Mgr;

    Mgr.SetProc(_MuMu_Bkg_Full);
    cout << "MuMu Hists location: " << Mgr.HistLocation << endl;

    TFile* f_Bkg_MuMu = new TFile(Mgr.HistLocation+"Hist_"+Mgr.Procname[_MuMu_Bkg_Full]+".root", "READ");
    if (f_Bkg_MuMu->IsOpen()) std::cout << "File Hist_"+Mgr.Procname[_MuMu_Bkg_Full]+".root opened successfully\n";
    else { std::cout << "File Hist_"+Mgr.Procname[_MuMu_Bkg_Full]+".root was not opened\n"; return; }

    Mgr.SetProc(_MuMu_DY_Full);
    TFile* f_DY_MuMu = new TFile(Mgr.HistLocation+"Hist_"+Mgr.Procname[_MuMu_DY_Full]+".root", "READ");
    if (f_DY_MuMu->IsOpen()) std::cout << "File Hist_"+Mgr.Procname[_MuMu_DY_Full]+".root opened successfully\n";
    else { std::cout << "File Hist_"+Mgr.Procname[_MuMu_DY_Full]+".root was not opened\n"; return; }

    Mgr.SetProc(_EMu_Bkg_Full);
    cout << "EMu Hists location: " << Mgr.HistLocation << endl;

    TFile* f_Bkg_EMu = new TFile(Mgr.HistLocation+"Hist_"+Mgr.Procname[_EMu_Bkg_Full]+".root", "READ");
    if (f_Bkg_EMu->IsOpen()) std::cout << "File Hist_"+Mgr.Procname[_EMu_Bkg_Full]+".root opened successfully\n";
    else { std::cout << "File Hist_"+Mgr.Procname[_EMu_Bkg_Full]+".root was not opened\n"; return; }

    Mgr.SetProc(_EMu_SingleMuon_Full);
    TFile* f_Data_EMu = new TFile(Mgr.HistLocation+"Hist_"+Mgr.Procname[_EMu_SingleMuon_Full]+".root", "READ");
    if (f_Data_EMu->IsOpen()) std::cout << "File Hist_"+Mgr.Procname[_EMu_SingleMuon_Full]+".root opened successfully\n";
    else { std::cout << "File Hist_"+Mgr.Procname[_EMu_SingleMuon_Full]+".root was not opened\n"; return; }

//********************************* E MU ***********************************************************

//--------------------------------- Get EMu QCD and WJets from FR -----------------------------------

    TH1D *h_EMu_QCDFR_invm, *h_EMu_SS_QCDFR_invm;
    TH1D *h_EMu_WJetsFR_invm, *h_EMu_SS_WJetsFR_invm;
    if (UseFR)
    {
        TFile *f_QCD_emu = new TFile(Mgr.HistLocation+"EstQCD_EMu.root", "READ");
        if (f_QCD_emu->IsOpen()) std::cout << "File EstQCD_EMu.root opened successfully\n";
        else { std::cout << "File EstQCD_EMu.root was not opened\n"; return; }
        f_QCD_emu->GetObject("h_QCD_est", h_EMu_QCDFR_invm);
        f_QCD_emu->GetObject("h_QCD_est_SS", h_EMu_SS_QCDFR_invm);
        h_EMu_QCDFR_invm->SetFillColor(kRed+3);
        h_EMu_QCDFR_invm->SetLineColor(kRed+3);
        h_EMu_QCDFR_invm->SetDirectory(0);
        h_EMu_SS_QCDFR_invm->SetFillColor(kRed+3);
        h_EMu_SS_QCDFR_invm->SetLineColor(kRed+3);
        h_EMu_SS_QCDFR_invm->SetDirectory(0);
        f_QCD_emu->Close();
        if (!f_QCD_emu->IsOpen()) std::cout << "File EstQCD_EMu.root closed successfully\n";

        TFile *f_WJets_emu = new TFile(Mgr.HistLocation+"EstWJets_EMu.root", "READ");
        if (f_WJets_emu->IsOpen()) std::cout << "File EstWJets_EMu.root opened successfully\n";
        else { std::cout << "File EstWJets_EMu.root was not opened\n"; return; }
        f_WJets_emu->GetObject("h_WJET_est", h_EMu_WJetsFR_invm);
        f_WJets_emu->GetObject("h_WJET_est_SS", h_EMu_SS_WJetsFR_invm);
        h_EMu_WJetsFR_invm->SetFillColor(kRed-2);
        h_EMu_WJetsFR_invm->SetLineColor(kRed-2);
        h_EMu_WJetsFR_invm->SetDirectory(0);
        h_EMu_SS_WJetsFR_invm->SetFillColor(kRed-2);
        h_EMu_SS_WJetsFR_invm->SetLineColor(kRed-2);
        h_EMu_SS_WJetsFR_invm->SetDirectory(0);
        f_WJets_emu->Close();
        if (!f_WJets_emu->IsOpen()) std::cout << "File EstWJets_EMu.root closed successfully\n";
    }

//---------------------------------- MC -----------------------------------------------------------

    TH1D* h_EMu_invm[_EndOf_EMu];
    TH1D* h_EMu_invm2[_EndOf_EMu];
    TH1D* h_EMu_SS_invm[_EndOf_EMu];
    TH1D* h_EMu_SS_invm2[_EndOf_EMu];
    THStack* s_EMu_invm = new THStack("s_emu_invm", "");
    THStack* s_EMu_invm2 = new THStack("s_emu_invm2", "");
    THStack* s_EMu_SS_invm = new THStack("s_emu_SS_invm", "");
    THStack* s_EMu_SS_invm2 = new THStack("s_emu_SS_invm2", "");
    gStyle->SetOptStat(0);

    Int_t isWJ = 1; // Add MC WJets if ==1
    if (UseFR) isWJ = 0;

    for (SelProc_t pr=_EMu_WJets_Full; pr>_EndOf_EMu_ttbar_Normal; pr=SelProc_t((int)(pr-1)))
    {
        Mgr.SetProc(pr);
        if (pr==_EndOf_EMu_VVnST_Normal) continue;
        if (isWJ == 0 && pr == _EMu_WJets_Full) {pr = _EndOf_EMu_VVnST_Normal; continue;}
        f_Bkg_EMu->GetObject("h_emu_mass_"+Mgr.Procname[pr], h_EMu_invm[pr]);
        f_Bkg_EMu->GetObject("h_emu_mass2_"+Mgr.Procname[pr], h_EMu_invm2[pr]);
        removeNegativeBins(h_EMu_invm[pr]);
        removeNegativeBins(h_EMu_invm2[pr]);

        Color_t color = kBlack;
        if (pr == _EMu_WJets_Full) color = kRed - 2;
        if (pr == _EMu_WW) color = kMagenta - 5;
        if (pr == _EMu_WZ) color = kMagenta - 2;
        if (pr == _EMu_ZZ) color = kMagenta - 6;
        if (pr == _EMu_tbarW) color = kGreen - 2;
        if (pr == _EMu_tW) color = kGreen + 2;
        if (pr == _EMu_ttbar_Full) color = kCyan + 2;
        if (pr == _EMu_DYTauTau_Full) color = kOrange - 5;

        h_EMu_invm[pr]->SetFillColor(color);
        h_EMu_invm[pr]->SetLineColor(color);
        h_EMu_invm[pr]->SetDirectory(0);
        s_EMu_invm->Add(h_EMu_invm[pr]);
        h_EMu_invm2[pr]->SetFillColor(color);
        h_EMu_invm2[pr]->SetLineColor(color);
        h_EMu_invm2[pr]->SetDirectory(0);
        s_EMu_invm2->Add(h_EMu_invm2[pr]);

        f_Bkg_EMu->GetObject("h_emuSS_mass_"+Mgr.Procname[pr], h_EMu_SS_invm[pr]);
        f_Bkg_EMu->GetObject("h_emuSS_mass2_"+Mgr.Procname[pr], h_EMu_SS_invm2[pr]);
        removeNegativeBins(h_EMu_SS_invm[pr]);
        removeNegativeBins(h_EMu_SS_invm2[pr]);

        h_EMu_SS_invm[pr]->SetFillColor(color);
        h_EMu_SS_invm[pr]->SetLineColor(color);
        h_EMu_SS_invm[pr]->SetDirectory(0);
        s_EMu_SS_invm->Add(h_EMu_SS_invm[pr]);
        h_EMu_SS_invm2[pr]->SetFillColor(color);
        h_EMu_SS_invm2[pr]->SetLineColor(color);
        h_EMu_SS_invm2[pr]->SetDirectory(0);
        s_EMu_SS_invm2->Add(h_EMu_SS_invm2[pr]);

        if (pr == _EMu_WJets_Full) pr = _EndOf_EMu_VVnST_Normal; // next -- WW
        if (pr == _EMu_tW) pr = _EMu_VVnST; // next -- ttbar
        if (pr == _EMu_DYTauTau_Full) break;
    }

//----------------------------- Data -------------------------------------------------------------

    Mgr.SetProc(_EMu_SingleMuon_Full);
    TH1D *h_EMu_data_invm, *h_EMu_data_invm2;
    f_Data_EMu->GetObject("h_emu_mass_"+Mgr.Procname[_EMu_SingleMuon_Full], h_EMu_data_invm);
    f_Data_EMu->GetObject("h_emu_mass2_"+Mgr.Procname[_EMu_SingleMuon_Full], h_EMu_data_invm2);
    h_EMu_data_invm->SetLineColorAlpha(kBlack, 0);
    h_EMu_data_invm->SetMarkerStyle(kFullDotLarge);
    h_EMu_data_invm->SetMarkerColor(kBlack);
    h_EMu_data_invm->SetDirectory(0);
    h_EMu_data_invm2->SetLineColorAlpha(kBlack, 0);
    h_EMu_data_invm2->SetMarkerStyle(kFullDotLarge);
    h_EMu_data_invm2->SetMarkerColor(kBlack);
    h_EMu_data_invm2->SetDirectory(0);

    TH1D *h_EMu_SS_data_invm, *h_EMu_SS_data_invm2;
    f_Data_EMu->GetObject("h_emuSS_mass_"+Mgr.Procname[_EMu_SingleMuon_Full], h_EMu_SS_data_invm);
    f_Data_EMu->GetObject("h_emuSS_mass2_"+Mgr.Procname[_EMu_SingleMuon_Full], h_EMu_SS_data_invm2);
    h_EMu_SS_data_invm->SetLineColorAlpha(kBlack, 0);
    h_EMu_SS_data_invm->SetMarkerStyle(kFullDotLarge);
    h_EMu_SS_data_invm->SetMarkerColor(kBlack);
    h_EMu_SS_data_invm->SetDirectory(0);
    h_EMu_SS_data_invm2->SetLineColorAlpha(kBlack, 0);
    h_EMu_SS_data_invm2->SetMarkerStyle(kFullDotLarge);
    h_EMu_SS_data_invm2->SetMarkerColor(kBlack);
    h_EMu_SS_data_invm2->SetDirectory(0);

//------------------------------ Est EMu QCD events ---------------------------------------------------

    TH1D* h_EMu_QCD_invm = ((TH1D*)(h_EMu_SS_data_invm->Clone("h_emu_mass_QCD")));
    TH1D* h_EMu_QCD_invm2 = ((TH1D*)(h_EMu_SS_data_invm2->Clone("h_emu_mass_QCD2")));
    h_EMu_QCD_invm->Add(((TH1D*)(s_EMu_SS_invm->GetStack()->Last())), -1);
    h_EMu_QCD_invm2->Add(((TH1D*)(s_EMu_SS_invm2->GetStack()->Last())), -1);
    const double RR = 0.57147108645;
    h_EMu_QCD_invm->Scale(1/RR);
    h_EMu_QCD_invm2->Scale(1/RR);
    removeNegativeBins(h_EMu_QCD_invm);
    removeNegativeBins(h_EMu_QCD_invm2);

    for (Int_t i=1; i<h_EMu_QCD_invm->GetSize()-1; i++)
    {
        if (h_EMu_QCD_invm->GetBinContent(i) < 0)
            h_EMu_QCD_invm->SetBinContent(i, 0);
    }
    for (Int_t i=1; i<h_EMu_QCD_invm2->GetSize()-1; i++)
    {
        if (h_EMu_QCD_invm2->GetBinContent(i) < 0)
            h_EMu_QCD_invm2->SetBinContent(i, 0);
    }

    h_EMu_QCD_invm->SetFillColor(kRed+3);
    h_EMu_QCD_invm->SetLineColor(kRed+3);
    h_EMu_QCD_invm->SetDirectory(0);
    h_EMu_QCD_invm2->SetFillColor(kRed+3);
    h_EMu_QCD_invm2->SetLineColor(kRed+3);
    h_EMu_QCD_invm2->SetDirectory(0);

// ------------------------------------ Draw stack histograms -------------------------------------

    THStack *s_EMu_wQCD_invm = new THStack("s_emu_wQCD_invm", "");
    THStack *s_EMu_wQCD_invm2 = new THStack("s_emu_wQCD_invm2", "");

    if (UseFR)
    {
        s_EMu_wQCD_invm->Add(h_EMu_QCDFR_invm);
        s_EMu_wQCD_invm->Add(h_EMu_WJetsFR_invm);
    }
    else
        s_EMu_wQCD_invm->Add(h_EMu_QCD_invm);
    s_EMu_wQCD_invm2->Add(h_EMu_QCD_invm2);
    for (SelProc_t pr=_EMu_WJets_Full; pr>_EndOf_EMu_ttbar_Normal; pr=SelProc_t((int)(pr-1)))
    {
        if (pr == _EndOf_EMu_VVnST_Normal) continue;
        if (isWJ == 0 && pr == _EMu_WJets_Full) {pr = _EndOf_EMu_VVnST_Normal; continue;}
        Mgr.SetProc(pr);
        s_EMu_wQCD_invm->Add(h_EMu_invm[pr]);
        s_EMu_wQCD_invm2->Add(h_EMu_invm2[pr]);
        if (pr == _EMu_WJets_Full) pr = _EndOf_EMu_VVnST_Normal; // next -- WW
        if (pr == _EMu_tW) pr = _EMu_VVnST; // next -- ttbar
        if (pr == _EMu_DYTauTau_Full) break;
    }

    myRatioPlot_t* RP_EMu_wQCD_invm = new myRatioPlot_t("EMu_wQCD_mass", s_EMu_wQCD_invm, h_EMu_data_invm);
    myRatioPlot_t* RP_EMu_wQCD_invm2 = new myRatioPlot_t("EMu_wQCD_mass2", s_EMu_wQCD_invm2, h_EMu_data_invm2);
    RP_EMu_wQCD_invm->SetPlots("m_{#lower[-0.2]{#scale[1.2]{#font[12]{e}#mu}}} [GeV/c^{2}]", 15, 3000);
    RP_EMu_wQCD_invm2->SetPlots("m_{#lower[-0.2]{#scale[1.2]{#font[12]{e}#mu}}} [GeV/c^{2}]", 15, 3000);

    TLegend *legend_EMu = new TLegend(0.77, 0.5, 0.95, 0.95);
    legend_EMu->AddEntry(h_EMu_data_invm, "Data", "lp");
    legend_EMu->AddEntry(h_EMu_invm[_EMu_DYTauTau_Full], "DY#rightarrow#tau#tau (MC)","f");
    legend_EMu->AddEntry(h_EMu_invm[_EMu_ttbar_Full], "#kern[0.2]{#font[12]{#scale[1.1]{t#bar{t}}}} (MC)", "f");
    legend_EMu->AddEntry(h_EMu_invm[_EMu_tW], "#kern[0.1]{#font[12]{#scale[1.1]{tW}}} (MC)", "f");
    legend_EMu->AddEntry(h_EMu_invm[_EMu_tbarW], "#kern[0.1]{#font[12]{#scale[1.1]{#bar{t}W}}} (MC)", "f");
    legend_EMu->AddEntry(h_EMu_invm[_EMu_ZZ], "#kern[0.1]{#font[12]{#scale[1.1]{ZZ}}} (MC)", "f");
    legend_EMu->AddEntry(h_EMu_invm[_EMu_WZ], "#font[12]{#scale[1.1]{WZ}} (MC)", "f");
    legend_EMu->AddEntry(h_EMu_invm[_EMu_WW], "#font[12]{#scale[1.1]{WW}} (MC)", "f");
    if (UseFR)
    {
        legend_EMu->AddEntry(h_EMu_WJetsFR_invm, "#font[12]{#scale[1.1]{W}}+Jets (est.)", "f");
        legend_EMu->AddEntry(h_EMu_QCDFR_invm, "#font[12]{#scale[1.1]{QCD}} (est.)", "f");
    }
    else
    {
        if (isWJ) legend_EMu->AddEntry(h_EMu_invm[_EMu_WJets_Full], "#font[12]{#scale[1.1]{W}}+Jets (MC)", "f");
        legend_EMu->AddEntry(h_EMu_QCD_invm, "#font[12]{#scale[1.1]{QCD}} (est.)", "f");
    }

    RP_EMu_wQCD_invm->ImportLegend(legend_EMu);
    RP_EMu_wQCD_invm2->ImportLegend(legend_EMu);
    RP_EMu_wQCD_invm->Draw(4e-1, 2e4, 1);
    RP_EMu_wQCD_invm2->Draw(4e-1, 2e4, 1);

    TH1D *h_EMu_data_invm_old = ((TH1D*)(h_EMu_data_invm->Clone("h_EMu_data_invm_old")));
    TH1D *h_EMu_data_invm_old2 = ((TH1D*)(h_EMu_data_invm2->Clone("h_EMu_data_invm_old2")));
    if (UseFR)
    {
        h_EMu_data_invm->Add(h_EMu_QCDFR_invm, -1);
        h_EMu_data_invm->Add(h_EMu_WJetsFR_invm, -1);
    }
    else
        h_EMu_data_invm->Add(h_EMu_QCD_invm, -1);
    h_EMu_data_invm2->Add(h_EMu_QCD_invm2, -1);
    removeNegativeBins(h_EMu_QCD_invm);
    removeNegativeBins(h_EMu_QCD_invm2);

    Double_t dataerror_emu, MCerror_emu, dataintegral_emu=2.25081e+07, MCintegral_emu, QCDerror, QCDintegral, QCDFRerror, QCDFRintegral;
    Double_t dataerrorZ_emu, MCerrorZ_emu, dataintegralZ_emu=2.25081e+07, MCintegralZ_emu;
    Double_t dataerror_noZ_emu=0, MCerror_noZ_emu=0, dataintegral_noZ_emu=2.25081e+07, MCintegral_noZ_emu, temp_noZ_emu;

    dataintegral_emu = h_EMu_data_invm->IntegralAndError(1, h_EMu_data_invm->GetSize()-2, dataerror_emu);
    MCintegral_emu = ((TH1D*)(s_EMu_invm->GetStack()->Last()))->IntegralAndError(1, h_EMu_data_invm->GetSize()-2, MCerror_emu);

    dataintegralZ_emu = h_EMu_data_invm->IntegralAndError(10, 22, dataerrorZ_emu);
    MCintegralZ_emu = ((TH1D*)(s_EMu_invm->GetStack()->Last()))->IntegralAndError(10, 22, MCerrorZ_emu);

    QCDintegral = h_EMu_QCD_invm->IntegralAndError(1, h_EMu_QCD_invm->GetSize()-2, QCDerror);
    if (UseFR) QCDFRintegral = h_EMu_QCDFR_invm->IntegralAndError(1, h_EMu_QCDFR_invm->GetSize()-2, QCDFRerror);

    dataintegral_noZ_emu = h_EMu_data_invm->IntegralAndError(1, 9, temp_noZ_emu);
    dataerror_noZ_emu += temp_noZ_emu * temp_noZ_emu;
    dataintegral_noZ_emu += h_EMu_data_invm->IntegralAndError(23, h_EMu_data_invm->GetSize()-2, temp_noZ_emu);
    dataerror_noZ_emu += temp_noZ_emu * temp_noZ_emu;
    dataerror_noZ_emu = sqrt(dataerror_noZ_emu);

    MCintegral_noZ_emu = ((TH1D*)(s_EMu_invm->GetStack()->Last()))->IntegralAndError(1, 9, temp_noZ_emu);
    MCerror_noZ_emu += temp_noZ_emu * temp_noZ_emu;
    MCintegral_noZ_emu += ((TH1D*)(s_EMu_invm->GetStack()->Last()))->IntegralAndError(23, h_EMu_data_invm->GetSize()-2, temp_noZ_emu);
    MCerror_noZ_emu += temp_noZ_emu * temp_noZ_emu;
    MCerror_noZ_emu = sqrt(MCerror_noZ_emu);

    std::cout << "EMu data events: " << dataintegral_emu << "+-" << dataerror_emu << endl;
    std::cout << "EMu MC events: " << MCintegral_emu << "+-" << MCerror_emu << endl;
    std::cout << "EMu MC/Obs: " << MCintegral_emu / dataintegral_emu << "+-" <<
                 sqrt((dataerror_emu / dataintegral_emu) * (dataerror_emu / dataintegral_emu) +
                       (MCerror_emu / MCintegral_emu) * (MCerror_emu / MCintegral_emu)) << endl;

    std::cout << "EMu data events around Z: " << dataintegralZ_emu << "+-" << dataerrorZ_emu << endl;
    std::cout << "EMu MC events around Z: " << MCintegralZ_emu << "+-" << MCerrorZ_emu << endl;
    std::cout << "EMu data events outside Z: " << dataintegral_noZ_emu << "+-" << dataerror_noZ_emu << endl;
    std::cout << "EMu MC events outside Z: " << MCintegral_noZ_emu << "+-" << MCerror_noZ_emu << endl << endl;
    std::cout << "EMu QCD events (from SS): " << QCDintegral << "+-" << QCDerror << endl;
    if (UseFR) std::cout << "EMu QCD events (FR): " << QCDFRintegral << "+-" << QCDFRerror << endl;

//************************************ MuMu *********************************************************

    TH1D *h_MuMu_invm[_EndOf_EMu], *h_MuMu_invm2[_EndOf_EMu];
    THStack* s_MuMu_invm = new THStack("s_MuMu_invm", "");
    THStack* s_MuMu_invm2 = new THStack("s_MuMu_invm2", "");

    for (SelProc_t pr=_MuMu_WW; pr<=_MuMu_VVnST; pr=SelProc_t((int)(pr-1)))
    {
        if (pr == _MuMu_WZ || pr == _MuMu_ZZ) continue;
        Mgr.SetProc(pr);

        f_Bkg_MuMu->GetObject("h_mass_"+Mgr.Procname[pr], h_MuMu_invm[pr]);
        f_Bkg_MuMu->GetObject("h_mass2_"+Mgr.Procname[pr], h_MuMu_invm2[pr]);
        removeNegativeBins(h_MuMu_invm[pr]);
        removeNegativeBins(h_MuMu_invm2[pr]);

        Color_t color = kBlack;
        if (pr == _MuMu_WW) color = kMagenta - 5;
        if (pr == _MuMu_WZ) color = kMagenta - 2;
        if (pr == _MuMu_ZZ) color = kMagenta - 6;
        if (pr == _MuMu_tbarW) color = kGreen - 2;
        if (pr == _MuMu_tW) color = kGreen + 2;
        if (pr == _MuMu_ttbar_Full) color = kCyan + 2;
        if (pr == _MuMu_DYTauTau_Full) color = kOrange - 5;

        h_MuMu_invm[pr]->SetFillColor(color);
        h_MuMu_invm[pr]->SetLineColor(color);
        h_MuMu_invm[pr]->SetDirectory(0);
        s_MuMu_invm->Add(h_MuMu_invm[pr]);

        h_MuMu_invm2[pr]->SetFillColor(color);
        h_MuMu_invm2[pr]->SetLineColor(color);
        h_MuMu_invm2[pr]->SetDirectory(0);
        s_MuMu_invm2->Add(h_MuMu_invm2[pr]);

        if (pr == _MuMu_tW) pr = _MuMu_VVnST; // next -- ttbar
        if (pr == _MuMu_DYTauTau_Full) break;
    }

//####################################### ALL AT ONCE ###############################################

    TH1D* h_MuMu_Est_invm = ((TH1D*)(h_EMu_data_invm->Clone("h_MuMu_mass_Est")));
    h_MuMu_Est_invm->Multiply(((TH1D*)(s_MuMu_invm->GetStack()->Last())));
    h_MuMu_Est_invm->Divide(((TH1D*)(s_EMu_invm->GetStack()->Last())));
    removeNegativeBins(h_MuMu_Est_invm);
    h_MuMu_Est_invm->SetDirectory(0);
    h_MuMu_Est_invm->SetLineColor(kBlack);
    h_MuMu_Est_invm->SetMarkerStyle(kFullTriangleUp);
    h_MuMu_Est_invm->SetMarkerSize(1.5);
    h_MuMu_Est_invm->SetMarkerColor(kBlack);

    TH1D* h_MuMu_Est_invm2 = ((TH1D*)(h_EMu_data_invm2->Clone("h_MuMu_mass_Est2")));
    h_MuMu_Est_invm2->Multiply(((TH1D*)(s_MuMu_invm2->GetStack()->Last())));
    h_MuMu_Est_invm2->Divide(((TH1D*)(s_EMu_invm2->GetStack()->Last())));
    removeNegativeBins(h_MuMu_Est_invm2);
    h_MuMu_Est_invm2->SetDirectory(0);
    h_MuMu_Est_invm2->SetLineColor(kBlack);
    h_MuMu_Est_invm2->SetMarkerStyle(kFullTriangleUp);
    h_MuMu_Est_invm2->SetMarkerSize(1.5);
    h_MuMu_Est_invm2->SetMarkerColor(kBlack);

/// ################################# Uncertainties (by hand) ######################################### ///

    Double_t UNC_MuMuMC[43], UNC_emuMC[43], UNC_emuData[43], UNC_MuMuEst[43], UNC_MuMu_EstOverMC[43],
             UNC_MuMuMC2[86], UNC_emuMC2[86], UNC_emuData2[86], UNC_MuMuEst2[86], UNC_MuMu_EstOverMC2[86];
    Double_t tmp = 0;
    for (int i=1; i<87; i++)
    {
        if (i < 44)
        {
            // emuMC (N2)
            UNC_emuMC[i-1] = ((TH1D*)(s_EMu_invm->GetStack()->Last()))->GetBinError(i) * ((TH1D*)(s_EMu_invm->GetStack()->Last()))->GetBinError(i);
            tmp = 0;
            for (SelProc_t pr=_EMu_WJets_Full; pr>_EndOf_EMu_ttbar_Normal; pr=SelProc_t((int)(pr-1)))
            {
                Mgr.SetProc(pr);
                if (pr==_EndOf_EMu_VVnST_Normal) continue;
                if (isWJ == 0 && pr == _EMu_WJets_Full) {pr = _EndOf_EMu_VVnST_Normal; continue;}

                tmp += h_EMu_invm[pr]->GetBinError(i) * h_EMu_invm[pr]->GetBinError(i);

                if (pr == _EMu_WJets_Full) // next - WW
                    pr = _EndOf_EMu_VVnST_Normal;
                if (pr == _EMu_tW) pr = _EMu_VVnST; // next -- ttbar
                if (pr == _EMu_DYTauTau_Full) break;
            }
            if ((UNC_emuMC[i-1]-tmp)/tmp > 0.0001)
            {
                cout << "EMU errors from THStack and TH1D sum do not match!!" << endl;
                cout << "Bin " << i << endl;
                cout << "From THStack: " << UNC_emuMC[i-1] << "   From TH1D sum: " << tmp << endl;
                break;
            }

            // MuMuMC (N1)
            UNC_MuMuMC[i-1] = ((TH1D*)(s_MuMu_invm->GetStack()->Last()))->GetBinError(i) * ((TH1D*)(s_MuMu_invm->GetStack()->Last()))->GetBinError(i);
            tmp = 0;
            for (SelProc_t pr=_MuMu_WW; pr<=_MuMu_VVnST; pr=SelProc_t((int)(pr-1)))
            {
                if (pr == _MuMu_WZ || pr == _MuMu_ZZ) continue;
                Mgr.SetProc(pr);

                tmp += h_MuMu_invm[pr]->GetBinError(i) * h_MuMu_invm[pr]->GetBinError(i);

                if (pr == _MuMu_tW) pr = _MuMu_VVnST; // next -- ttbar
                if (pr == _MuMu_DYTauTau_Full) break;
            }
            if ((UNC_MuMuMC[i-1]-tmp)/tmp > 0.0001)
            {
                cout << "MuMu errors from THStack and TH1D sum do not match!!" << endl;
                cout << "Bin " << i << endl;
                cout << "From THStack: " << UNC_MuMuMC[i-1] << "   From TH1D sum: " << tmp << endl;
                break;
            }

            // emuData (N3)
            UNC_emuData[i-1] = ((TH1D*)(s_EMu_SS_invm->GetStack()->Last()))->GetBinError(i) * ((TH1D*)(s_EMu_SS_invm->GetStack()->Last()))->GetBinError(i);
            tmp = 0;
            for (SelProc_t pr=_EMu_WJets_Full; pr>_EndOf_EMu_ttbar_Normal; pr=SelProc_t((int)(pr-1)))
            {
                Mgr.SetProc(pr);
                if (pr==_EndOf_EMu_VVnST_Normal) continue;
                if (isWJ == 0 && pr == _EMu_WJets_Full) {pr = _EndOf_EMu_VVnST_Normal; continue;}

                tmp += h_EMu_SS_invm[pr]->GetBinError(i) * h_EMu_SS_invm[pr]->GetBinError(i);

                if (pr == _EMu_WJets_Full) // next - WW
                    pr = _EndOf_EMu_VVnST_Normal;
                if (pr == _EMu_tW) pr = _EMu_VVnST; // next -- ttbar
                if (pr == _EMu_DYTauTau_Full) break;
            }
            if ((UNC_emuData[i-1]-tmp)/tmp > 0.0001)
            {
                cout << "EMU SS errors from THStack and TH1D sum do not match!!" << endl;
                cout << "Bin " << i << endl;
                cout << "From THStack: " << UNC_emuData[i-1] << "   From TH1D sum: " << tmp << endl;
                break;
            }
            UNC_emuData[i-1] *= 1/(RR*RR);
            UNC_emuData[i-1] += h_EMu_SS_data_invm->GetBinError(i) * h_EMu_SS_data_invm->GetBinError(i) / (RR*RR);
            UNC_emuData[i-1] += h_EMu_data_invm_old->GetBinError(i) * h_EMu_data_invm_old->GetBinError(i);
            tmp = h_EMu_data_invm->GetBinError(i) * h_EMu_data_invm->GetBinError(i);
            if ((UNC_emuData[i-1]-tmp)/tmp > 0.0001)
            {
                cout << "EMU DATA errors from automatic and manual counting do not match!!" << endl;
                cout << "Bin " << i << endl;
                cout << "From automatic: " << UNC_emuData[i-1] << "   From manual: " << tmp << endl;
                break;
            }

            // FULL UNC
            tmp = UNC_MuMuMC[i-1] / (((TH1D*)(s_MuMu_invm->GetStack()->Last()))->GetBinContent(i) * ((TH1D*)(s_MuMu_invm->GetStack()->Last()))->GetBinContent(i));
            tmp += UNC_emuMC[i-1] / (((TH1D*)(s_EMu_invm->GetStack()->Last()))->GetBinContent(i) * ((TH1D*)(s_EMu_invm->GetStack()->Last()))->GetBinContent(i));
            tmp += UNC_emuData[i-1] / (h_EMu_data_invm->GetBinContent(i) * h_EMu_data_invm->GetBinContent(i));
            UNC_MuMuEst[i-1] = sqrt(tmp) * fabs(h_MuMu_Est_invm->GetBinContent(i));
            tmp = h_MuMu_Est_invm->GetBinError(i) * h_MuMu_Est_invm->GetBinError(i);
            if ((UNC_MuMuEst[i-1]-tmp)/tmp > 0.0001)
            {
                cout << "MuMu EST errors from automatic and manual counting do not match!!" << endl;
                cout << "Bin " << i << endl;
                cout << "From automatic: " << UNC_MuMuEst[i-1] << "   From manual: " << tmp << endl;
                break;
            }

        } // End of if(i<44)

        // emuMC (N2)
        UNC_emuMC2[i-1] = ((TH1D*)(s_EMu_invm2->GetStack()->Last()))->GetBinError(i) * ((TH1D*)(s_EMu_invm2->GetStack()->Last()))->GetBinError(i);
        tmp = 0;
        for (SelProc_t pr=_EMu_WJets_Full; pr>_EndOf_EMu_ttbar_Normal; pr=SelProc_t((int)(pr-1)))
        {
            Mgr.SetProc(pr);
            if (pr==_EndOf_EMu_VVnST_Normal) continue;
            if (isWJ == 0 && pr == _EMu_WJets_Full) {pr = _EndOf_EMu_VVnST_Normal; continue;}

            tmp += h_EMu_invm2[pr]->GetBinError(i) * h_EMu_invm2[pr]->GetBinError(i);

            if (pr == _EMu_WJets_Full) // next - WW
                pr = _EndOf_EMu_VVnST_Normal;
            if (pr == _EMu_tW) pr = _EMu_VVnST; // next -- ttbar
            if (pr == _EMu_DYTauTau_Full) break;
        }
        if ((UNC_emuMC2[i-1]-tmp)/tmp > 0.0001)
        {
            cout << "EMU2 errors from THStack and TH1D sum do not match!!" << endl;
            cout << "Bin " << i << endl;
            cout << "From THStack: " << UNC_emuMC2[i-1] << "   From TH1D sum: " << tmp << endl;
            break;
        }

        // MuMuMC (N1)
        UNC_MuMuMC2[i-1] = ((TH1D*)(s_MuMu_invm2->GetStack()->Last()))->GetBinError(i) * ((TH1D*)(s_MuMu_invm2->GetStack()->Last()))->GetBinError(i);
        tmp = 0;
        for (SelProc_t pr=_MuMu_WW; pr<=_MuMu_VVnST; pr=SelProc_t((int)(pr-1)))
        {
            if (pr == _MuMu_WZ || pr == _MuMu_ZZ) continue;
            Mgr.SetProc(pr);

            tmp += h_MuMu_invm2[pr]->GetBinError(i) * h_MuMu_invm2[pr]->GetBinError(i);

            if (pr == _MuMu_tW) pr = _MuMu_VVnST; // next -- ttbar
            if (pr == _MuMu_DYTauTau_Full) break;
        }
        if ((UNC_MuMuMC2[i-1]-tmp)/tmp > 0.0001)
        {
            cout << "MuMu2 errors from THStack and TH1D sum do not match!!" << endl;
            cout << "Bin " << i << endl;
            cout << "From THStack: " << UNC_MuMuMC2[i-1] << "   From TH1D sum: " << tmp << endl;
            break;
        }

        // emuData (N3)
        UNC_emuData2[i-1] = ((TH1D*)(s_EMu_SS_invm2->GetStack()->Last()))->GetBinError(i) * ((TH1D*)(s_EMu_SS_invm2->GetStack()->Last()))->GetBinError(i);
        tmp = 0;
        for (SelProc_t pr=_EMu_WJets_Full; pr>_EndOf_EMu_ttbar_Normal; pr=SelProc_t((int)(pr-1)))
        {
            Mgr.SetProc(pr);
            if (pr==_EndOf_EMu_VVnST_Normal) continue;
            if (isWJ == 0 && pr == _EMu_WJets_Full) {pr = _EndOf_EMu_VVnST_Normal; continue;}

            tmp += h_EMu_SS_invm2[pr]->GetBinError(i) * h_EMu_SS_invm2[pr]->GetBinError(i);

            if (pr == _EMu_WJets_Full) // next - WW
                pr = _EndOf_EMu_VVnST_Normal;
            if (pr == _EMu_tW) pr = _EMu_VVnST; // next -- ttbar
            if (pr == _EMu_DYTauTau_Full) break;
        }
        if ((UNC_emuData2[i-1]-tmp)/tmp > 0.0001)
        {
            cout << "EMU2 SS errors from THStack and TH1D sum do not match!!" << endl;
            cout << "Bin " << i << endl;
            cout << "From THStack: " << UNC_emuData2[i-1] << "   From TH1D sum: " << tmp << endl;
            break;
        }
        UNC_emuData2[i-1] *= 1/(RR*RR);
        UNC_emuData2[i-1] += h_EMu_SS_data_invm2->GetBinError(i) * h_EMu_SS_data_invm2->GetBinError(i) / (RR*RR);
        UNC_emuData2[i-1] += h_EMu_data_invm_old2->GetBinError(i) * h_EMu_data_invm_old2->GetBinError(i);
        tmp = h_EMu_data_invm2->GetBinError(i) * h_EMu_data_invm2->GetBinError(i);
        if ((UNC_emuData2[i-1]-tmp)/tmp > 0.0001)
        {
            cout << "EMU2 DATA errors from automatic and manual counting do not match!!" << endl;
            cout << "Bin " << i << endl;
            cout << "From automatic: " << UNC_emuData2[i-1] << "   From manual: " << tmp << endl;
            break;
        }

        // FULL UNC
        tmp = UNC_MuMuMC2[i-1] / (((TH1D*)(s_MuMu_invm2->GetStack()->Last()))->GetBinContent(i) * ((TH1D*)(s_MuMu_invm2->GetStack()->Last()))->GetBinContent(i));
        tmp += UNC_emuMC2[i-1] / (((TH1D*)(s_EMu_invm2->GetStack()->Last()))->GetBinContent(i) * ((TH1D*)(s_EMu_invm2->GetStack()->Last()))->GetBinContent(i));
        tmp += UNC_emuData2[i-1] / (h_EMu_data_invm2->GetBinContent(i) * h_EMu_data_invm2->GetBinContent(i));
        UNC_MuMuEst2[i-1] = sqrt(tmp) * fabs(h_MuMu_Est_invm2->GetBinContent(i));
        tmp = h_MuMu_Est_invm2->GetBinError(i) * h_MuMu_Est_invm2->GetBinError(i);
        if ((UNC_MuMuEst2[i-1]-tmp)/tmp > 0.0001)
        {
            cout << "MuMu2 EST errors from automatic and manual counting do not match!!" << endl;
            cout << "Bin " << i << endl;
            cout << "From automatic: " << UNC_MuMuEst2[i-1] << "   From manual: " << tmp << endl;
            break;
        }
    } // End of for(bins)
//    cout << "Manually calculated errors:" << endl;
//    for (int i=0; i<43; i++)
//    {
//        cout << UNC_MuMuEst[i] << "  ";
//    }
//    cout << "\nAutomatically calculated errors:" << endl;
//    for (int i=1; i<=43; i++)
//    {
//        cout << h_MuMu_Est_invm->GetBinError(i) << "  ";
//    }
//    cout << endl;

    // DATA/MC UNC (probably smaller than given automatically)
    for (int i=0; i<86; i++)
    {
        if (i < 43)
        {
            UNC_MuMu_EstOverMC[i] = ((TH1D*)(s_EMu_invm->GetStack()->Last()))->GetBinError(i+1) / ((TH1D*)(s_EMu_invm->GetStack()->Last()))->GetBinContent(i+1); // N2
            UNC_MuMu_EstOverMC[i] *= UNC_MuMu_EstOverMC[i];
            tmp = h_EMu_data_invm->GetBinError(i+1) / h_EMu_data_invm->GetBinContent(i+1); // N3
            tmp *= tmp;
            UNC_MuMu_EstOverMC[i] = sqrt(UNC_MuMu_EstOverMC[i] + tmp);
        }
        UNC_MuMu_EstOverMC2[i] = ((TH1D*)(s_EMu_invm2->GetStack()->Last()))->GetBinError(i+1) / ((TH1D*)(s_EMu_invm2->GetStack()->Last()))->GetBinContent(i+1); // N2
        UNC_MuMu_EstOverMC2[i] *= UNC_MuMu_EstOverMC2[i];
        tmp = h_EMu_data_invm2->GetBinError(i+1) / h_EMu_data_invm2->GetBinContent(i+1); // N3
        tmp *= tmp;
        UNC_MuMu_EstOverMC2[i] = sqrt(UNC_MuMu_EstOverMC2[i] + tmp);
    }

// ############################# Systematics (DataDriven - MC) ##################################### //

    Double_t systematics[h_MuMu_Est_invm->GetSize()-2], systematics2[h_MuMu_Est_invm2->GetSize()-2], systerr=0,
             estSystematics[h_MuMu_Est_invm->GetSize()-2], estSystematics2[h_MuMu_Est_invm2->GetSize()-2];
    for (int i=1; i<h_MuMu_Est_invm2->GetSize()-1; i++)
    {
        if (i < h_MuMu_Est_invm->GetSize()-1)
        {
            estSystematics[i-1] = fabs(h_MuMu_Est_invm->GetBinContent(i) - ((TH1D*)(s_MuMu_invm->GetStack()->Last()))->GetBinContent(i));
            systematics[i-1] = estSystematics[i-1] / fabs(((TH1D*)(s_MuMu_invm->GetStack()->Last()))->GetBinContent(i));
            systerr += systematics[i-1] * systematics[i-1];
        }
        estSystematics2[i-1] = fabs(h_MuMu_Est_invm2->GetBinContent(i) - ((TH1D*)(s_MuMu_invm2->GetStack()->Last()))->GetBinContent(i));
        systematics2[i-1] = estSystematics2[i-1] / fabs(((TH1D*)(s_MuMu_invm2->GetStack()->Last()))->GetBinContent(i));
    }

// ###################################### Writing ######################################### //

    TFile* f_NMuMuEst = new TFile(Mgr.HistLocation+"EstBkg_MuMu.root", "RECREATE");
    if (f_NMuMuEst->IsOpen()) std::cout << "File EstBkg_MuMu.root opened successfully\n";
    else { std::cout << "File EstBkg_MuMu.root was not opened\n"; return; }
    f_NMuMuEst->cd();
    h_MuMu_Est_invm->Write();
    h_MuMu_Est_invm2->Write();
    TVectorD estSyst(43, estSystematics);
    TVectorD estSyst2(86, estSystematics);

    estSyst.Write("estSystematics");
    estSyst2.Write("estSystematics2");

    Double_t MCerror, esterror, MCintegral, estintegral;
    Double_t MCerrorZ, esterrorZ, MCintegralZ, estintegralZ;
    Double_t MCerror_noZ=0, esterror_noZ=0, MCintegral_noZ, estintegral_noZ, temp_noZ;

    MCintegral = ((TH1D*)(s_MuMu_invm->GetStack()->Last()))->IntegralAndError(1, h_EMu_data_invm->GetSize()-2, MCerror);
    estintegral = h_MuMu_Est_invm->IntegralAndError(1, h_MuMu_Est_invm->GetSize()-2, esterror);

    MCintegralZ = ((TH1D*)(s_MuMu_invm->GetStack()->Last()))->IntegralAndError(10, 22, MCerrorZ);
    estintegralZ = h_MuMu_Est_invm->IntegralAndError(10, 22, esterrorZ);

    MCintegral_noZ = ((TH1D*)(s_MuMu_invm->GetStack()->Last()))->IntegralAndError(1, 9, temp_noZ);
    MCerror_noZ += temp_noZ * temp_noZ;
    MCintegral_noZ += ((TH1D*)(s_MuMu_invm->GetStack()->Last()))->IntegralAndError(23, h_EMu_data_invm->GetSize()-2, temp_noZ);
    MCerror_noZ += temp_noZ * temp_noZ;
    MCerror_noZ = sqrt(MCerror_noZ);

    estintegral_noZ = h_MuMu_Est_invm->IntegralAndError(1, 9, temp_noZ);
    esterror_noZ += temp_noZ * temp_noZ;
    estintegral_noZ += h_MuMu_Est_invm->IntegralAndError(23, h_MuMu_Est_invm->GetSize()-2, temp_noZ);
    esterror_noZ += temp_noZ * temp_noZ;
    esterror_noZ = sqrt(esterror_noZ);

    std::cout << "MuMu MC events: " << MCintegral << "+-" << MCerror << endl;
    std::cout << "MuMu Est. events: " << estintegral << "+-" << esterror << "+-" << systerr << endl << endl;

    std::cout << "MuMu MC events around Z: " << MCintegralZ << "+-" << MCerrorZ << endl;
    std::cout << "MuMu Est. events around Z: " << estintegralZ << "+-" << esterrorZ << endl;
    std::cout << "MuMu MC events outside Z: " << MCintegral_noZ << "+-" << MCerror_noZ << endl;
    std::cout << "MuMu Est. events outside Z: " << estintegral_noZ << "+-" << esterror_noZ << endl << endl;

//------------------------- Z PEAK ------------------------------------------------------------------------

    Double_t MCZerr;
    Double_t MCZint = ((TH1D*)(s_MuMu_invm->GetStack()->Last()))->IntegralAndError(10, 23, MCZerr);
    Double_t EMuZerr;
    Double_t EMuZint = h_MuMu_Est_invm->IntegralAndError(10, 23, EMuZerr);
    std::cout << "MC events around Z: " << MCZint << " +- " << MCZerr << endl;
    std::cout << "EMuEst events around Z: " << EMuZint << " +- " << EMuZerr << endl;

//------------------------------ Drawing ----------------------------------------------------------

    myRatioPlot_t *RP_invm = new myRatioPlot_t("DataDriven_InvariantMass", s_MuMu_invm, h_MuMu_Est_invm);
    myRatioPlot_t *RP_invm2 = new myRatioPlot_t("DataDriven_InvariantMass2", s_MuMu_invm2, h_MuMu_Est_invm2);
    RP_invm->SetPlots("m_{#lower[-0.2]{#scale[1.15]{#mu#mu}}} [GeV/c^{2}]", 15, 3000, "Est./MC", UNC_MuMu_EstOverMC);
    RP_invm2->SetPlots("m_{#lower[-0.2]{#scale[1.15]{#mu#mu}}} [GeV/c^{2}]", 15, 3000, "Est./MC", UNC_MuMu_EstOverMC2);
    RP_invm->SetSystematics(estSystematics, NULL, systematics);
    RP_invm2->SetSystematics(estSystematics2, NULL, systematics2);

    TLegend *legend = new TLegend(0.77, 0.5, 0.95, 0.95);
    legend->AddEntry(h_MuMu_Est_invm, "Estimation", "lp");
    legend->AddEntry(h_MuMu_invm[_MuMu_DYTauTau_Full], "DY#rightarrow#tau#tau","f");
    legend->AddEntry(h_MuMu_invm[_MuMu_ttbar_Full], "#kern[0.2]{#font[12]{#scale[1.1]{t#bar{t}}}}", "f");
    legend->AddEntry(h_MuMu_invm[_MuMu_tW], "#kern[0.1]{#font[12]{#scale[1.1]{tW}}}", "f");
    legend->AddEntry(h_MuMu_invm[_MuMu_tbarW], "#kern[0.1]{#font[12]{#scale[1.1]{#bar{t}W}}}", "f");
    legend->AddEntry(h_MuMu_invm[_MuMu_WW], "#font[12]{#scale[1.1]{WW}}", "f");

    RP_invm->ImportLegend(legend);
    RP_invm2->ImportLegend(legend);
    RP_invm->Draw(4e-1, 2e4, 1);
    RP_invm2->Draw(4e-1, 2e4, 1);

//######################################## ONE BY ONE #################################################

    TH1D *h_MuMu_Est_pr_invm[_EndOf_EMu], *h_MuMu_Est_pr_invm2[_EndOf_EMu];
    myRatioPlot_t *RP_MuMu_invm_pr[_EndOf_EMu], *RP_MuMu_invm_pr2[_EndOf_EMu];
    Double_t emuMCInt[_EndOf_EMu], MuMuMCInt[_EndOf_EMu], MuMu_estInt[_EndOf_EMu], emuMCEr[_EndOf_EMu], MuMuMCEr[_EndOf_EMu], MuMu_estEr[_EndOf_EMu];

    for (SelProc_t pr=_MuMu_tW; pr<_MuMu_VVnST; pr=next(pr))
    {
        if (pr>=_EndOf_MuMu_VVnST_Normal && pr<_MuMu_DYTauTau_Full) continue;
        if (pr == _MuMu_WZ || pr == _MuMu_ZZ) continue;

        Mgr.SetProc(pr);
        h_MuMu_Est_pr_invm[pr] = ((TH1D*)(h_MuMu_invm[pr]->Clone("h_MuMu_mass_Est_"+Mgr.Procname[pr])));
        removeNegativeBins(h_MuMu_Est_pr_invm[pr]);
        h_MuMu_Est_pr_invm[pr]->Multiply(h_EMu_data_invm);
        removeNegativeBins(h_MuMu_Est_pr_invm[pr]);
        h_MuMu_Est_pr_invm[pr]->Divide((TH1D*)(s_EMu_invm->GetStack()->Last()));
        removeNegativeBins(h_MuMu_Est_pr_invm[pr]);

        h_MuMu_Est_pr_invm[pr]->SetLineColor(kBlack);
        h_MuMu_Est_pr_invm[pr]->SetDirectory(0);
        h_MuMu_Est_pr_invm[pr]->Write();

        h_MuMu_Est_pr_invm2[pr] = ((TH1D*)(h_MuMu_invm2[pr]->Clone("h_MuMu_mass_Est2_"+Mgr.Procname[pr])));
        removeNegativeBins(h_MuMu_Est_pr_invm2[pr]);
        h_MuMu_Est_pr_invm2[pr]->Multiply(h_EMu_data_invm2);
        removeNegativeBins(h_MuMu_Est_pr_invm2[pr]);
        h_MuMu_Est_pr_invm2[pr]->Divide((TH1D*)(s_EMu_invm2->GetStack()->Last()));
        removeNegativeBins(h_MuMu_Est_pr_invm2[pr]);

        h_MuMu_Est_pr_invm2[pr]->SetLineColor(kBlack);
        h_MuMu_Est_pr_invm2[pr]->SetDirectory(0);
        h_MuMu_Est_pr_invm2[pr]->Write();

        emuMCEr[pr] = 0;
        if (pr <= _MuMu_WW) // Matching EMu with EE processes
            emuMCInt[pr] = h_EMu_invm[pr+115]->IntegralAndError(1, h_EMu_invm[pr+115]->GetSize()-1, emuMCEr[pr]);
        else
            emuMCInt[pr] = h_EMu_invm[pr+101]->IntegralAndError(1, h_EMu_invm[pr+101]->GetSize()-1, emuMCEr[pr]);
        MuMuMCEr[pr] = 0;
        MuMuMCInt[pr] = h_MuMu_invm[pr]->IntegralAndError(1, h_MuMu_invm[pr]->GetSize()-1, MuMuMCEr[pr]);
        MuMu_estEr[pr] = 0;
        MuMu_estInt[pr] = h_MuMu_Est_pr_invm[pr]->IntegralAndError(1, h_MuMu_Est_pr_invm[pr]->GetSize()-1, MuMu_estEr[pr]);

        std::cout << Mgr.Procname[pr] << " MC EMu events: " << emuMCInt[pr] << " +- " << emuMCEr[pr] << " (Stat.)\n";
        std::cout << Mgr.Procname[pr] << " MC MuMu events: " << MuMuMCInt[pr] << " +- " << MuMuMCEr[pr] << " (Stat.)\n";
        std::cout << Mgr.Procname[pr] << " MuMu events estimated from data: " << MuMu_estInt[pr] << " +- " << MuMu_estEr[pr] << " (Stat.)\n\n";

//------------------------------ Drawing ----------------------------------------------------------

        RP_MuMu_invm_pr[pr] = new myRatioPlot_t("h_mass_DataDriven_"+Mgr.Procname[pr], h_MuMu_invm[pr], h_MuMu_Est_pr_invm[pr]);
        RP_MuMu_invm_pr2[pr] = new myRatioPlot_t("h_mass_DataDriven2_"+Mgr.Procname[pr], h_MuMu_invm2[pr], h_MuMu_Est_pr_invm2[pr]);
        RP_MuMu_invm_pr[pr]->SetPlots("m_{#lower[-0.2]{#scale[1.15]{#mu#mu}}} [GeV/c^{2}]", 15, 3000, "Est./MC", UNC_MuMu_EstOverMC);
        RP_MuMu_invm_pr2[pr]->SetPlots("m_{#lower[-0.2]{#scale[1.15]{#mu#mu}}} [GeV/c^{2}]", 15, 3000, "Est./MC", UNC_MuMu_EstOverMC2);
        RP_MuMu_invm_pr[pr]->SetSystematics(NULL, NULL, systematics);
        RP_MuMu_invm_pr2[pr]->SetSystematics(NULL, NULL, systematics2);

        TLegend *legend_pr = new TLegend(0.65, 0.75, 0.95, 0.95);

        legend_pr->AddEntry(h_MuMu_Est_pr_invm[pr], "Estimation", "lp");

        if (pr==_MuMu_tW)
            legend_pr->AddEntry(h_MuMu_invm[pr], "#kern[0.1]{#font[12]{#scale[1.1]{tW}}} (MC)", "f");
        else if (pr==_MuMu_tbarW)
            legend_pr->AddEntry(h_MuMu_invm[pr], "#kern[0.1]{#font[12]{#scale[1.1]{#bar{t}W}}} (MC)", "f");
        else if (pr==_MuMu_ZZ)
            legend_pr->AddEntry(h_MuMu_invm[pr], "#kern[0.1]{#font[12]{#scale[1.1]{ZZ}}} (MC)", "f");
        else if (pr==_MuMu_WZ)
            legend_pr->AddEntry(h_MuMu_invm[pr], "#font[12]{#scale[1.1]{WZ}} (MC)", "f");
        else if (pr==_MuMu_WW)
            legend_pr->AddEntry(h_MuMu_invm[pr], "#font[12]{#scale[1.1]{WW}} (MC)", "f");
        else if (pr==_MuMu_ttbar_Full)
            legend_pr->AddEntry(h_MuMu_invm[pr], "#kern[0.2]{#font[12]{#scale[1.1]{t#bar{t}}}} (MC)", "f");
        else if (pr==_MuMu_DYTauTau_Full)
            legend_pr->AddEntry(h_MuMu_invm[pr], "DY#rightarrow#tau#tau (MC)", "f");
        else
            legend_pr->AddEntry(h_MuMu_invm[pr], Mgr.Procname[pr]+" (MC)", "f");

        RP_MuMu_invm_pr[pr]->ImportLegend(legend_pr);
        RP_MuMu_invm_pr2[pr]->ImportLegend(legend_pr);
        RP_MuMu_invm_pr[pr]->Draw(4e-1, 1e4, 1);
        RP_MuMu_invm_pr2[pr]->Draw(4e-1, 1e4, 1);
    }

//################################## tW + tbarW #######################################################

    THStack* s_tW_invm = new THStack("s_tW", "");
    THStack* s_tW_invm2 = new THStack("s_tW2", "");
    s_tW_invm->Add(h_MuMu_invm[_MuMu_tbarW]);
    s_tW_invm->Add(h_MuMu_invm[_MuMu_tW]);
    s_tW_invm2->Add(h_MuMu_invm2[_MuMu_tbarW]);
    s_tW_invm2->Add(h_MuMu_invm2[_MuMu_tW]);

    TH1D* h1_MuMuTW_invm_Est = ((TH1D*)(h_MuMu_Est_pr_invm[_MuMu_tW]->Clone("h_mass_MuMuTW+TbarW_Est")));
    h1_MuMuTW_invm_Est->Add(h_MuMu_Est_pr_invm[_MuMu_tbarW]);
    h1_MuMuTW_invm_Est->SetDirectory(0);
    h1_MuMuTW_invm_Est->Write();
    TH1D* h1_MuMuTW_invm_Est2 = ((TH1D*)(h_MuMu_Est_pr_invm2[_MuMu_tW]->Clone("h_mass_MuMuTW+TbarW_Est2")));
    h1_MuMuTW_invm_Est2->Add(h_MuMu_Est_pr_invm2[_MuMu_tbarW]);
    h1_MuMuTW_invm_Est2->SetDirectory(0);
    h1_MuMuTW_invm_Est2->Write();

    myRatioPlot_t* RP_tW = new myRatioPlot_t("tWest", s_tW_invm, h1_MuMuTW_invm_Est);
    myRatioPlot_t* RP_tW2 = new myRatioPlot_t("tWest2", s_tW_invm2, h1_MuMuTW_invm_Est2);
    RP_tW->SetPlots("m_{#lower[-0.2]{#scale[1.15]{#mu#mu}}} [GeV/c^{2}]", 15, 3000, "Est./MC");
    RP_tW2->SetPlots("m_{#lower[-0.2]{#scale[1.15]{#mu#mu}}} [GeV/c^{2}]", 15, 3000, "Est./MC");
    RP_tW->SetSystematics(NULL, NULL, systematics);
    RP_tW2->SetSystematics(NULL, NULL, systematics2);

    TLegend *legend_tW = new TLegend(0.65, 0.65, 0.95, 0.95);
    legend_tW->AddEntry(h1_MuMuTW_invm_Est, "Estimation", "lp");
    legend_tW->AddEntry(h_MuMu_invm[_MuMu_tW], "#kern[0.1]{#font[12]{#scale[1.1]{tW}}} (MC)", "f");
    legend_tW->AddEntry(h_MuMu_invm[_MuMu_tbarW], "#kern[0.1]{#font[12]{#scale[1.1]{#bar{t}W}}} (MC)", "f");

    RP_tW->ImportLegend(legend_tW);
    RP_tW2->ImportLegend(legend_tW);
    RP_tW->Draw(4e-1, 1e4, 1);
    RP_tW2->Draw(4e-1, 1e4, 1);

//--------------------------------------- Closing up --------------------------------------------

    f_Bkg_EMu->Close();
    if (!f_Bkg_EMu->IsOpen()) std::cout << "File Hist_"+Mgr.Procname[_EMu_Bkg_Full]+".root closed successfully\n";
    f_Data_EMu->Close();
    if (!f_Data_EMu->IsOpen()) std::cout << "File Hist_"+Mgr.Procname[_EMu_SingleMuon_Full]+".root closed successfully\n";
    f_Bkg_MuMu->Close();
    if (!f_Bkg_MuMu->IsOpen()) std::cout << "File Hist_"+Mgr.Procname[_MuMu_Bkg_Full]+".root closed successfully\n";
    f_DY_MuMu->Close();
    if(!f_DY_MuMu->IsOpen()) std::cout << "File Hist_"+Mgr.Procname[_MuMu_DY_Full]+".root closed successfully.\n";
    f_NMuMuEst->Close();
    if (!f_NMuMuEst->IsOpen()) std::cout << "File EstBkg_MuMu.root closed successfully\n\n";

}// End of MuMu_Est


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
