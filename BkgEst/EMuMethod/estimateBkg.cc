#include <TStyle.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TGraphAsymmErrors.h>
#include <TMultiGraph.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <THStack.h>
#include <TMath.h>
#include <TText.h>
#include <TPad.h>
#include <TPaveText.h>
#include <TLorentzVector.h>
#include <TStopwatch.h>
#include <TColor.h>
#include <TLatex.h>
#include <TEfficiency.h>

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>

const int Nsmooth = 1;
const double AvgWeight = 1.000493;
const double lumi = 35867.0;

const int binnum = 43;
const double bins[44] = {15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 64, 68, 72, 76, 81, 86, 91, 96, 101, 106, 110, 115, 120, 126, 133, 141, 150, 160, 171, 185,  200, 220, 243, 273, 320, 380, 440, 510, 600, 700,  830, 1000, 1500, 3000};

using namespace std;

void fillSystematics( TH1D* data_driven, TH1D* stat, TH1D* systematic, TH1D* total );
void removeNegativeBins( TH1D* hist );
void mkplot( THStack* h_mc, TH1D* h_data, TLegend* leg, TString plotname );

TString inputname = "./root/test/ROOTFile_Muon_Channel_1.root"; // merged root file of all output ones from dimuon and emu event selection
TString output = "TopPtReweighting_on_MuMu";
TString outputname = "./result/estimatedBkg/"+output;


void estimateBkg() {

    int c1 = 3;
    int c2 = 2;
    int c3 = 5;
    int c4 = 13;
    int c5 = 14;
    int c6 = 15;
    int c7 = 16;
    int c8 = 17;

    ////////////////////////////
    // -- Set MC histogram -- //
    ////////////////////////////
    TFile f_input(inputname, "read");
    TH1D *emu_ttbar = (TH1D*)f_input.Get("h_emu_mass_ttbar");
    TH1D *emu_ttbarBackup = (TH1D*)f_input.Get("h_emu_mass_ttbarBackup");
    TH1D *emu_ttbar_M700to1000 = (TH1D*)f_input.Get("h_emu_mass_ttbar_M700to1000");
    TH1D *emu_ttbar_M1000toInf = (TH1D*)f_input.Get("h_emu_mass_ttbar_M1000toInf");
    TH1D *emu_DYtautau = (TH1D*)f_input.Get("h_emu_mass_DYTauTau");
    TH1D *emu_DYtautau_M10to50_v1 = (TH1D*)f_input.Get("h_emu_mass_DYTauTau_M10to50_v1");
    TH1D *emu_DYtautau_M10to50_v2 = (TH1D*)f_input.Get("h_emu_mass_DYTauTau_M10to50_v2");
    TH1D *emu_DYtautau_M10to50_ext1v1 = (TH1D*)f_input.Get("h_emu_mass_DYTauTau_M10to50_ext1v1");
    TH1D *emu_WW = (TH1D*)f_input.Get("h_emu_mass_WW");
    TH1D *emu_WZ = (TH1D*)f_input.Get("h_emu_mass_WZ");
    TH1D *emu_ZZ = (TH1D*)f_input.Get("h_emu_mass_ZZ");
    TH1D *emu_tW = (TH1D*)f_input.Get("h_emu_mass_tW");
    TH1D *emu_antitW = (TH1D*)f_input.Get("h_emu_mass_tbarW");

    TH1D *emuSS_ttbar = (TH1D*)f_input.Get("h_emuSS_mass_ttbar");
    TH1D *emuSS_ttbarBackup = (TH1D*)f_input.Get("h_emuSS_mass_ttbarBackup");
    TH1D *emuSS_ttbar_M700to1000 = (TH1D*)f_input.Get("h_emuSS_mass_ttbar_M700to1000");
    TH1D *emuSS_ttbar_M1000toInf = (TH1D*)f_input.Get("h_emuSS_mass_ttbar_M1000toInf");
    TH1D *emuSS_DYtautau = (TH1D*)f_input.Get("h_emuSS_mass_DYTauTau");
    TH1D *emuSS_DYtautau_M10to50_v1 = (TH1D*)f_input.Get("h_emuSS_mass_DYTauTau_M10to50_v1");
    TH1D *emuSS_DYtautau_M10to50_v2 = (TH1D*)f_input.Get("h_emuSS_mass_DYTauTau_M10to50_v2");
    TH1D *emuSS_DYtautau_M10to50_ext1v1 = (TH1D*)f_input.Get("h_emuSS_mass_DYTauTau_M10to50_ext1v1");
    TH1D *emuSS_WW = (TH1D*)f_input.Get("h_emuSS_mass_WW");
    TH1D *emuSS_WZ = (TH1D*)f_input.Get("h_emuSS_mass_WZ");
    TH1D *emuSS_ZZ = (TH1D*)f_input.Get("h_emuSS_mass_ZZ");
    TH1D *emuSS_tW = (TH1D*)f_input.Get("h_emuSS_mass_tW");
    TH1D *emuSS_antitW = (TH1D*)f_input.Get("h_emuSS_mass_tbarW");

    TH1D *dimu_ttbar = (TH1D*)f_input.Get("h_mass_ttbar");
    TH1D *dimu_ttbarBackup = (TH1D*)f_input.Get("h_mass_ttbarBackup");
    TH1D *dimu_ttbar_M700to1000 = (TH1D*)f_input.Get("h_mass_ttbar_M700to1000");
    TH1D *dimu_ttbar_M1000toInf = (TH1D*)f_input.Get("h_mass_ttbar_M1000toInf");
    TH1D *dimu_DYtautau = (TH1D*)f_input.Get("h_mass_DYTauTau");
    TH1D *dimu_DYtautau_M10to50_v1 = (TH1D*)f_input.Get("h_mass_DYTauTau_M10to50_v1");
    TH1D *dimu_DYtautau_M10to50_v2 = (TH1D*)f_input.Get("h_mass_DYTauTau_M10to50_v2");
    TH1D *dimu_DYtautau_M10to50_ext1v1 = (TH1D*)f_input.Get("h_mass_DYTauTau_M10to50_ext1v1");
    TH1D *dimu_WW = (TH1D*)f_input.Get("h_mass_WW");
    TH1D *dimu_WZ = (TH1D*)f_input.Get("h_mass_WZ");
    TH1D *dimu_ZZ = (TH1D*)f_input.Get("h_mass_ZZ");
    TH1D *dimu_tW = (TH1D*)f_input.Get("h_mass_tW");
    TH1D *dimu_antitW = (TH1D*)f_input.Get("h_mass_tbarW");


    // -- Merge ttbar samples -- //
    emu_ttbar->Add(emu_ttbarBackup);
    emu_ttbar->Add(emu_ttbar_M700to1000);
    emu_ttbar->Add(emu_ttbar_M1000toInf);
    emuSS_ttbar->Add(emuSS_ttbarBackup);
    emuSS_ttbar->Add(emuSS_ttbar_M700to1000);
    emuSS_ttbar->Add(emuSS_ttbar_M1000toInf);
    dimu_ttbar->Add(dimu_ttbarBackup);
    dimu_ttbar->Add(dimu_ttbar_M700to1000);
    dimu_ttbar->Add(dimu_ttbar_M1000toInf);

    // DY M50 + M10to50
    emu_DYtautau->Add(emu_DYtautau_M10to50_v1);
    emu_DYtautau->Add(emu_DYtautau_M10to50_v2);
    emu_DYtautau->Add(emu_DYtautau_M10to50_ext1v1);
    emuSS_DYtautau->Add(emuSS_DYtautau_M10to50_v1);
    emuSS_DYtautau->Add(emuSS_DYtautau_M10to50_v2);
    emuSS_DYtautau->Add(emuSS_DYtautau_M10to50_ext1v1);
    dimu_DYtautau->Add(dimu_DYtautau_M10to50_v1);
    dimu_DYtautau->Add(dimu_DYtautau_M10to50_v2);
    dimu_DYtautau->Add(dimu_DYtautau_M10to50_ext1v1);

    // WW + WZ + ZZ
    TH1D* emu_diboson = emu_WW;
    emu_diboson->Add(emu_WZ);
    emu_diboson->Add(emu_ZZ);
    TH1D* emuSS_diboson = emuSS_WW;
    emuSS_diboson->Add(emuSS_WZ);
    emuSS_diboson->Add(emuSS_ZZ);

    // tW + antitW
    emu_tW->Add(emu_antitW);
    emuSS_tW->Add(emuSS_antitW);
    dimu_tW->Add(dimu_antitW);


    // -- Consider average weight of top pt re-weight -- //
    emu_ttbar->Scale(1/AvgWeight);
    emuSS_ttbar->Scale(1/AvgWeight);
    dimu_ttbar->Scale(1/AvgWeight);


    // -- Set Color -- //
    emu_ttbar->SetFillColor(c1);
    emu_DYtautau->SetFillColor(c2);
    emu_diboson->SetFillColor(c4);
    emu_tW->SetFillColor(c7);

    emuSS_ttbar->SetFillColor(c1);
    emuSS_DYtautau->SetFillColor(c2);
    emuSS_diboson->SetFillColor(c4);
    emuSS_tW->SetFillColor(c7);

    dimu_ttbar->SetFillColor(c1);
    dimu_DYtautau->SetFillColor(c2);
    dimu_WW->SetFillColor(c4);
    dimu_tW->SetFillColor(c7);


    // -- No Stats -- //
    emu_ttbar->SetStats(kFALSE);
    emu_DYtautau->SetStats(kFALSE);
    emu_diboson->SetStats(kFALSE);
    emu_tW->SetStats(kFALSE);

    emuSS_ttbar->SetStats(kFALSE);
    emuSS_DYtautau->SetStats(kFALSE);
    emuSS_diboson->SetStats(kFALSE);
    emuSS_tW->SetStats(kFALSE);

    dimu_ttbar->SetStats(kFALSE);
    dimu_DYtautau->SetStats(kFALSE);
    dimu_WW->SetStats(kFALSE);
    dimu_tW->SetStats(kFALSE);


    // -- Set Data histogram -- //
    TH1D *emu_data = (TH1D*)f_input.Get("h_emu_mass_Data");
    TH1D *emuSS_data = (TH1D*)f_input.Get("h_emuSS_mass_Data");

    emu_data->SetMarkerStyle(20);
    emu_data->SetMarkerSize(0.8);
    emu_data->SetStats(kFALSE);

    emuSS_data->SetMarkerStyle(20);
    emuSS_data->SetMarkerSize(0.8);
    emuSS_data->SetStats(kFALSE);


    // -- Remove negative bins -- //
    removeNegativeBins( emu_DYtautau );
    removeNegativeBins( emuSS_DYtautau );
    removeNegativeBins( dimu_DYtautau );


    // -- Stack histograms -- //
    TH1D* emu_sumBkg = new TH1D("emu_sumBkg","",binnum,bins);
    emu_sumBkg->Add(emu_DYtautau);	
    emu_sumBkg->Add(emu_ttbar);	
    emu_sumBkg->Add(emu_diboson);	
    emu_sumBkg->Add(emu_tW);	

    THStack* emu_stackBkg = new THStack("emu_stackBkg","");
    emu_stackBkg->Add(emu_diboson);
    emu_stackBkg->Add(emu_DYtautau);
    emu_stackBkg->Add(emu_tW);
    emu_stackBkg->Add(emu_ttbar);

    TH1D* emuSS_sumMC = new TH1D("emuSS_sumMC","",binnum,bins);
    emuSS_sumMC->Add(emuSS_DYtautau);  
    emuSS_sumMC->Add(emuSS_ttbar); 
    emuSS_sumMC->Add(emuSS_diboson); 
    emuSS_sumMC->Add(emuSS_tW);  

    THStack* emuSS_stackMC = new THStack("emuSS_stackMC","");
    emuSS_stackMC->Add(emuSS_diboson);
    emuSS_stackMC->Add(emuSS_DYtautau);
    emuSS_stackMC->Add(emuSS_tW);
    emuSS_stackMC->Add(emuSS_ttbar);


    // -- Legend -- //
    TLegend* legend = new TLegend(.75,.75,.95,.89);
    legend->AddEntry(emu_data,"Bkg(data-driven)");
    legend->AddEntry(dimu_ttbar,"ttbar","F");
    legend->AddEntry(dimu_DYtautau,"DY#tau#tau","F");
    legend->AddEntry(dimu_tW,"tW+#bar{t}W","F");
    legend->AddEntry(dimu_WW,"WW","F");
    legend->SetBorderSize(0);  
    legend->SetFillStyle(0);  

    TLegend* legend_ttbar = new TLegend(.75,.75,.95,.89);
    legend_ttbar->AddEntry(emu_data,"Bkg(data-driven)");
    legend_ttbar->AddEntry(dimu_ttbar,"ttbar","F");
    legend_ttbar->SetBorderSize(0);  
    legend_ttbar->SetFillStyle(0);  

    TLegend* legend_DYtautau = new TLegend(.75,.75,.95,.89);
    legend_DYtautau->AddEntry(emu_data,"Bkg(data-driven)");
    legend_DYtautau->AddEntry(dimu_DYtautau,"DY#tau#tau","F");
    legend_DYtautau->SetBorderSize(0);  
    legend_DYtautau->SetFillStyle(0);  

    TLegend* legend_tW = new TLegend(.75,.75,.95,.89);
    legend_tW->AddEntry(emu_data,"Bkg(data-driven)");
    legend_tW->AddEntry(dimu_tW,"tW+#bar{t}W","F");
    legend_tW->SetBorderSize(0);  
    legend_tW->SetFillStyle(0);  

    TLegend* legend_WW = new TLegend(.75,.75,.95,.89);
    legend_WW->AddEntry(emu_data,"Bkg(data-driven)");
    legend_WW->AddEntry(dimu_WW,"WW","F");
    legend_WW->SetBorderSize(0);  
    legend_WW->SetFillStyle(0);  


    // -- Estimate emu QCD -- //
    TH1D* emu_QCD = (TH1D*)emuSS_data->Clone();
    emu_QCD->Add(emuSS_DYtautau,-1.0);
    emu_QCD->Add(emuSS_ttbar,-1.0);
    emu_QCD->Add(emuSS_diboson,-1.0);
    emu_QCD->Add(emuSS_tW,-1.0);
    emu_QCD->SetFillColor(7);

    const double RR = 0.57147108645;
    emu_QCD->Scale(1/RR);

    removeNegativeBins(emu_QCD);

    //emu_QCD->Smooth(Nsmooth,"R"); //Smooting QCD
    //printf("Smoothing QCD : %d times\n", Nsmooth);


    // -- Ratio = emu_Data/emu_MC without QCD -- //
    TH1D* emu_data_nonQCD = (TH1D*)emu_data->Clone();
    emu_data_nonQCD->Add(emu_QCD,-1.0);

    TH1D* emu_ratio = (TH1D*)emu_data_nonQCD->Clone("emu_ratio");
    emu_ratio->Divide(emu_data_nonQCD,emu_sumBkg,1.0,1.0,"B");


    // -- Set data-driven histograms -- //
    removeNegativeBins(emu_ratio);

    TH1D* data_driven_ttbar = (TH1D*)emu_ratio->Clone("data_driven_ttbar");
    TH1D* data_driven_tW = (TH1D*)emu_ratio->Clone("data_driven_tW");
    TH1D* data_driven_WW = (TH1D*)emu_ratio->Clone("data_driven_WW");
    TH1D* data_driven_DYtautau = (TH1D*)emu_ratio->Clone("data_driven_DYtautau");

    data_driven_ttbar->Multiply(dimu_ttbar);
    data_driven_DYtautau->Multiply(dimu_DYtautau);
    data_driven_WW->Multiply(dimu_WW);
    data_driven_tW->Multiply(dimu_tW);


    // -- Replace the highest mass bins with MC (becasue data-driven entries are zero) -- //
    data_driven_ttbar->SetBinContent(binnum,dimu_ttbar->GetBinContent(binnum));
    data_driven_ttbar->SetBinError(binnum,dimu_ttbar->GetBinError(binnum)); 
    data_driven_tW->SetBinContent(binnum,dimu_tW->GetBinContent(binnum));
    data_driven_tW->SetBinError(binnum,dimu_tW->GetBinError(binnum)); 
    data_driven_WW->SetBinContent(binnum,dimu_WW->GetBinContent(binnum));
    data_driven_WW->SetBinError(binnum,dimu_WW->GetBinError(binnum)); 


    // -- Calculate systematic uncertainty -- //
    TH1D* ttbar_total = new TH1D("ttbar_total","",binnum,bins);
    TH1D* ttbar_systematic = new TH1D("ttbar_systematic","",binnum,bins);
    TH1D* ttbar_stat = new TH1D("ttbar_stat","",binnum,bins);

    TH1D* tW_total = new TH1D("tW_total","",binnum,bins);
    TH1D* tW_systematic = new TH1D("tW_systematic","",binnum,bins);
    TH1D* tW_stat = new TH1D("tW_stat","",binnum,bins);

    TH1D* WW_total = new TH1D("WW_total","",binnum,bins);
    TH1D* WW_systematic = new TH1D("WW_systematic","",binnum,bins);
    TH1D* WW_stat = new TH1D("WW_stat","",binnum,bins);

    TH1D* DYtautau_total = new TH1D("DYtautau_total","",binnum,bins);
    TH1D* DYtautau_systematic = new TH1D("DYtautau_systematic","",binnum,bins);
    TH1D* DYtautau_stat = new TH1D("DYtautau_stat","",binnum,bins);

    ttbar_systematic->Add(dimu_ttbar);
    ttbar_systematic->Add(data_driven_ttbar,-1.0);

    DYtautau_systematic->Add(dimu_DYtautau);
    DYtautau_systematic->Add(data_driven_DYtautau,-1.0);

    tW_systematic->Add(dimu_tW);
    tW_systematic->Add(data_driven_tW,-1.0);

    WW_systematic->Add(dimu_WW);
    WW_systematic->Add(data_driven_WW,-1.0);

    fillSystematics( data_driven_ttbar, ttbar_stat, ttbar_systematic, ttbar_total );
    fillSystematics( data_driven_DYtautau, DYtautau_stat, DYtautau_systematic, DYtautau_total );
    fillSystematics( data_driven_tW, tW_stat, tW_systematic, tW_total );
    fillSystematics( data_driven_WW, WW_stat, WW_systematic, WW_total );


    // -- Set Bkg histogram name -- //
    data_driven_ttbar->SetName("ttbar");  
    data_driven_DYtautau->SetName("DYtautau");
    data_driven_tW->SetName("tW");
    data_driven_WW->SetName("WW");

    dimu_ttbar->SetName("ttbar_MC");  
    dimu_DYtautau->SetName("DYtautau_MC");
    dimu_tW->SetName("tW_MC");
    dimu_WW->SetName("WW_MC");


    // -- Output ROOT file -- //
    TFile* g = new TFile("./result/estimatedBkg/emu_"+output+".root","RECREATE");

    data_driven_ttbar->Write();
    data_driven_DYtautau->Write();
    data_driven_tW->Write();
    data_driven_WW->Write();

    dimu_ttbar->Write();
    dimu_DYtautau->Write();
    dimu_tW->Write();
    dimu_WW->Write();

    ttbar_systematic->Write();
    DYtautau_systematic->Write();
    tW_systematic->Write();
    WW_systematic->Write();

    ttbar_stat->Write();
    DYtautau_stat->Write();
    tW_stat->Write();
    WW_stat->Write();

    g->Close();


    // -- Stack MC-Bkg histograms -- //
    THStack* h_ttbar = new THStack("h_ttbar","");
    h_ttbar->Add(dimu_ttbar);

    THStack* h_DYtautau = new THStack("h_DYtautau","");
    h_DYtautau->Add(dimu_DYtautau);

    THStack* h_tW = new THStack("h_tW","");
    h_tW->Add(dimu_tW);

    THStack* h_WW = new THStack("h_WW","");
    h_WW->Add(dimu_WW);

    THStack* h_all = new THStack("h_all","");
    h_all->Add(dimu_WW);
    h_all->Add(dimu_tW);
    h_all->Add(dimu_DYtautau);
    h_all->Add(dimu_ttbar);


    // -- Merge Data -- //
    TList* datalist = new TList;
    datalist->Add(data_driven_ttbar);
    datalist->Add(data_driven_DYtautau);
    datalist->Add(data_driven_tW);
    datalist->Add(data_driven_WW);
    
    TH1D* data_driven_all = new TH1D("data_driven_all", "",binnum, bins);
    data_driven_all->Merge(datalist);


    // -- Make plot -- //
    mkplot( h_ttbar, data_driven_ttbar, legend_ttbar, "ttbar" );
    mkplot( h_DYtautau, data_driven_DYtautau, legend_DYtautau, "DYtautau" );
    mkplot( h_tW, data_driven_tW, legend_tW, "tW" );
    mkplot( h_WW, data_driven_WW, legend_WW, "WW" );
    mkplot( h_all, data_driven_all, legend, "All" );
}


void fillSystematics( TH1D* data_driven, TH1D* stat, TH1D* systematic, TH1D* total ) {

    double binSystematic = 0;
    double binStat = 0;
    double binTotal = 0;

    for(int i=0; i<data_driven->GetNbinsX(); i++) {

        if(data_driven->GetBinContent(i+1)!=0) {  
            binStat = data_driven->GetBinError(i+1);
            binSystematic = fabs(systematic->GetBinContent(i+1));
            binTotal = sqrt(binSystematic*binSystematic + binStat*binStat);
        }
        else{
            binSystematic = 0;
            binStat = 0;
            binTotal = 0;
        }

        systematic->SetBinContent(i+1,binSystematic);
        stat->SetBinContent(i+1,binStat);    
        total->SetBinContent(i+1,binTotal);

        data_driven->SetBinError(i+1,binTotal);

    }

}

void removeNegativeBins( TH1D* hist ) {

    for(int i=0; i<hist->GetNbinsX(); i++) {
        if(hist->GetBinContent(i+1)<0) {
            hist->SetBinContent(i+1,0);
            hist->SetBinError(i+1,0);
        }
    }   

}

void mkplot( THStack* h_mc, TH1D* h_data, TLegend* leg, TString plotname ) {

    double L = 0.0913; // for left margin
    double R = 0.09188888888; // for right margin

    TCanvas *Mass1 = new TCanvas("Mass1", "",600,600);
    Mass1->SetFillColor(0);
    Mass1->SetBorderMode(0);
    Mass1->SetBorderSize(2);
    Mass1->SetFrameBorderMode(0);

    // -- Pad1 -- //
    TPad *InvariantMass = new TPad("InvariantMass", "",0.01,0.01,0.99,0.99);
    InvariantMass->Draw();
    InvariantMass->cd();
    InvariantMass->SetFillColor(0);
    InvariantMass->SetBorderMode(0);
    InvariantMass->SetBorderSize(2);
    InvariantMass->SetLogx();
    InvariantMass->SetLogy();
    InvariantMass->SetGridx();
    InvariantMass->SetGridy();
    InvariantMass->SetTickx();
    InvariantMass->SetTicky();
    InvariantMass->SetBottomMargin(0.3130424);
    InvariantMass->SetFrameBorderMode(0);

    TH1D *h_mass_stack_1;
    h_mass_stack_1 = new TH1D("h_mass_stack_1","",binnum, bins);
    h_mass_stack_1->SetDirectory(0);
    h_mass_stack_1->SetStats(0);
    h_mass_stack_1->GetXaxis()->SetLabelFont(42);
    h_mass_stack_1->GetXaxis()->SetLabelSize(0);
    h_mass_stack_1->GetXaxis()->SetTitleSize(0.035);
    h_mass_stack_1->GetXaxis()->SetTitleFont(42);
    h_mass_stack_1->GetYaxis()->SetLabelFont(42);
    h_mass_stack_1->GetYaxis()->SetLabelSize(0.035);
    h_mass_stack_1->GetYaxis()->SetTitleSize(0.045);
    h_mass_stack_1->GetYaxis()->SetTitleFont(42);
    h_mass_stack_1->GetYaxis()->SetTitle("# of events");
    h_mass_stack_1->GetYaxis()->SetTitleOffset(1.15);
    h_mass_stack_1->GetXaxis()->SetRangeUser(15,3000);

    h_mc->SetHistogram(h_mass_stack_1);
    h_mc->SetTitle("");
    h_mc->SetMaximum(100000);
    h_mc->SetMinimum(1);
    h_mc->Draw("hist");

    h_data->GetXaxis()->SetRangeUser(15,3000);
    h_data->SetMaximum(100000);
    h_data->SetMinimum(1);
    h_data->SetMarkerStyle(20);
    h_data->SetMarkerSize(0.8);
    h_data->Draw("same p e");

    leg->Draw();

    TLatex *    tex = new TLatex(0.74,0.91,"#font[42]{#scale[0.6]{2016, 13TeV}}");
    tex->SetNDC();
    tex->SetLineWidth(2);
    tex->Draw();

    tex = new TLatex(0.1,0.95,"#font[62]{CMS}");
    tex->SetNDC();
    tex->SetLineWidth(2);
    tex->Draw();

    tex = new TLatex(0.1,0.91,"#font[42]{#it{#scale[0.8]{Preliminary}}}");
    tex->SetNDC();
    tex->SetLineWidth(2);
    tex->Draw();

    tex = new TLatex(0.4,0.95,"#font[42]{#scale[0.6]{Run2016}}");
    tex->SetNDC();
    tex->SetLineWidth(2);
    tex->Draw();

    tex = new TLatex(0.4,0.91,"#font[42]{#scale[0.6]{Luminosity = 35.9 fb^{-1}}}");
    tex->SetNDC();
    tex->SetLineWidth(2);
    tex->Draw();

    // -- Bottom pad -- //
    TPad *padc1_2 = new TPad("padc1_2", "padc1_2",0.01,0.01,0.99,0.3);
    padc1_2->Draw();
    padc1_2->cd();
    padc1_2->SetFillColor(0);
    padc1_2->SetBorderMode(0);
    padc1_2->SetBorderSize(2);
    padc1_2->SetGridx();
    padc1_2->SetGridy();
    padc1_2->SetTickx();
    padc1_2->SetTicky();
    padc1_2->SetLogx();
    padc1_2->SetTopMargin(0.07176078);
    padc1_2->SetBottomMargin(0.4233841);
    padc1_2->SetLeftMargin(L);
    padc1_2->SetRightMargin(R);
    padc1_2->SetFrameBorderMode(0);

    // -- ratio plot -- //
    TList *list = (TList*)h_mc->GetHists();

    TH1D *sum_mc;
    sum_mc = new TH1D("sum_mc","",binnum, bins);
    sum_mc->Merge(list);
    sum_mc->Sumw2();

    TH1D *ratio_mass__7 = (TH1D*)h_data->Clone("ratio_mass__7");
    ratio_mass__7->Divide(h_data, sum_mc, 1.0, 1.0, "B");

    ratio_mass__7->SetMarkerStyle(20);
    ratio_mass__7->SetMarkerSize(0.8);
    if( output.Contains("mumu") )
        ratio_mass__7->GetXaxis()->SetTitle("M(#mu^{+}#mu^{-}) [GeV]");
    else if( output.Contains("ee") )
        ratio_mass__7->GetXaxis()->SetTitle("M(e^{+}e^{-}) [GeV]");
    ratio_mass__7->GetXaxis()->SetRangeUser(15,3000);
    ratio_mass__7->GetXaxis()->SetMoreLogLabels();
    ratio_mass__7->GetXaxis()->SetLabelFont(42);
    ratio_mass__7->GetXaxis()->SetLabelOffset(0.007);
    ratio_mass__7->GetXaxis()->SetLabelSize(0.12);
    ratio_mass__7->GetXaxis()->SetTitleSize(0.13);
    ratio_mass__7->GetXaxis()->SetTitleOffset(1.4);
    ratio_mass__7->GetXaxis()->SetTitleFont(42);
    ratio_mass__7->GetYaxis()->SetTitle("Data/MC");
    ratio_mass__7->GetYaxis()->SetRangeUser(0.7,1.3);
    ratio_mass__7->GetYaxis()->SetNdivisions(508);
    ratio_mass__7->GetYaxis()->SetLabelFont(42);
    ratio_mass__7->GetYaxis()->SetLabelSize(0.07);
    ratio_mass__7->GetYaxis()->SetTitleSize(0.13);
    ratio_mass__7->GetYaxis()->SetTitleOffset(0.35);
    ratio_mass__7->GetYaxis()->SetTitleFont(42);
    ratio_mass__7->SetStats(0);
    ratio_mass__7->Draw("p");

    TF1 *f_line1 = new TF1("f_line","1",-10000,10000);
    f_line1->SetLineColor(2);
    f_line1->SetLineWidth(2);
    f_line1->GetXaxis()->SetLabelFont(42);
    f_line1->GetXaxis()->SetLabelSize(0.035);
    f_line1->GetXaxis()->SetTitleSize(0.035);
    f_line1->GetXaxis()->SetTitleFont(42);
    f_line1->GetYaxis()->SetLabelFont(42);
    f_line1->GetYaxis()->SetLabelSize(0.035);
    f_line1->GetYaxis()->SetTitleSize(0.035);
    f_line1->GetYaxis()->SetTitleFont(42);
    f_line1->Draw("SAME");

    padc1_2->Modified();
    InvariantMass->cd();
    InvariantMass->Modified();
    Mass1->cd();
    Mass1->Modified();
    Mass1->cd();
    Mass1->SetSelected(Mass1);

    Mass1->SaveAs(outputname+"_"+plotname+".pdf");
}
