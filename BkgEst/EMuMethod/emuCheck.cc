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

const int rebin = 5;
const int Nsmooth = 1;
const double AvgWeight = 1.000493;
const double lumi = 35867.0;

const int binnum = 43;
const double bins[44] = {15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 64, 68, 72, 76, 81, 86, 91, 96, 101, 106, 110, 115, 120, 126, 133, 141, 150, 160, 171, 185, 200, 220, 243, 273, 320, 380, 440, 510, 600, 700, 830, 1000, 1500, 3000};

using namespace std;

void mkhisto( TString var );
void removeNegativeBins( TH1D* hist );
void massplot( THStack* h_mc, TH1D* h_data, TLegend* leg, TString plotname, Int_t QCDoption );
void commonplot( THStack* h_mc, TH1D* h_data, TLegend* leg, TString plotname, Int_t QCDoption );

TString inputname = "./root/test/ROOTFile_EMuSelection_1.root"; // merged root file of all output ones from emu event selection
TString output = "TopPtReweighting_on";
TString outputname = "./result/"+output;


void emuCheck() {

    cout << "input file : " << inputname << endl;
    cout << endl;

    mkhisto( "emu_mass" );
    mkhisto( "emu_mass_FineBin" );
    mkhisto( "mu_pT"    );
    mkhisto( "mu_eta"   );
    mkhisto( "mu_phi"   );
    mkhisto( "ele_pT"   );
    mkhisto( "ele_eta"  );
    mkhisto( "ele_phi"  );

}

void mkhisto( TString var ) {

    cout << "Running for [" << var << "]..." << endl;
    cout << endl;

    int c1 = 3;
    int c2 = 2;
    int c3 = 5;
    int c4 = 13;
    int c5 = 14;
    int c6 = 15;
    int c7 = 16;
    int c8 = 17;

    
    // -- Set MC histogram -- //
    TFile f_input(inputname, "read");
    TH1D *emu_ttbar = (TH1D*)f_input.Get("h_"+var+"_ttbar");
    TH1D *emu_ttbarBackup = (TH1D*)f_input.Get("h_"+var+"_ttbarBackup");
    TH1D *emu_ttbar_M700to1000 = (TH1D*)f_input.Get("h_"+var+"_ttbar_M700to1000");
    TH1D *emu_ttbar_M1000toInf = (TH1D*)f_input.Get("h_"+var+"_ttbar_M1000toInf");
    TH1D *emu_DYtautau = (TH1D*)f_input.Get("h_"+var+"_DYTauTau");
    TH1D *emu_DYtautau_M10to50_v1 = (TH1D*)f_input.Get("h_"+var+"_DYTauTau_M10to50_v1");
    TH1D *emu_DYtautau_M10to50_v2 = (TH1D*)f_input.Get("h_"+var+"_DYTauTau_M10to50_v2");
    TH1D *emu_DYtautau_M10to50_ext1v1 = (TH1D*)f_input.Get("h_"+var+"_DYTauTau_M10to50_ext1v1");
    TH1D *emu_WW = (TH1D*)f_input.Get("h_"+var+"_WW");
    TH1D *emu_WZ = (TH1D*)f_input.Get("h_"+var+"_WZ");
    TH1D *emu_ZZ = (TH1D*)f_input.Get("h_"+var+"_ZZ");
    TH1D *emu_tW = (TH1D*)f_input.Get("h_"+var+"_tW");
    TH1D *emu_antitW = (TH1D*)f_input.Get("h_"+var+"_tbarW");

    TString varSS;
    if( var == "emu_mass" ) varSS = "emuSS_mass";
    else if( var == "emu_mass_FineBin" ) varSS = "emuSS_mass_FineBin";
    else if( var == "mu_pT" ) varSS = "muSS_pT";
    else if( var == "mu_eta" ) varSS = "muSS_eta";
    else if( var == "mu_phi" ) varSS = "muSS_phi";
    else if( var == "ele_pT" ) varSS = "eleSS_pT";
    else if( var == "ele_eta" ) varSS = "eleSS_eta";
    else if( var == "ele_phi" ) varSS = "eleSS_phi";

    TH1D *emuSS_ttbar = (TH1D*)f_input.Get("h_"+varSS+"_ttbar");
    TH1D *emuSS_ttbarBackup = (TH1D*)f_input.Get("h_"+varSS+"_ttbarBackup");
    TH1D *emuSS_ttbar_M700to1000 = (TH1D*)f_input.Get("h_"+varSS+"_ttbar_M700to1000");
    TH1D *emuSS_ttbar_M1000toInf = (TH1D*)f_input.Get("h_"+varSS+"_ttbar_M1000toInf");
    TH1D *emuSS_DYtautau = (TH1D*)f_input.Get("h_"+varSS+"_DYTauTau");
    TH1D *emuSS_DYtautau_M10to50_v1 = (TH1D*)f_input.Get("h_"+varSS+"_DYTauTau_M10to50_v1");
    TH1D *emuSS_DYtautau_M10to50_v2 = (TH1D*)f_input.Get("h_"+varSS+"_DYTauTau_M10to50_v2");
    TH1D *emuSS_DYtautau_M10to50_ext1v1 = (TH1D*)f_input.Get("h_"+varSS+"_DYTauTau_M10to50_ext1v1");
    TH1D *emuSS_WW = (TH1D*)f_input.Get("h_"+varSS+"_WW");
    TH1D *emuSS_WZ = (TH1D*)f_input.Get("h_"+varSS+"_WZ");
    TH1D *emuSS_ZZ = (TH1D*)f_input.Get("h_"+varSS+"_ZZ");
    TH1D *emuSS_tW = (TH1D*)f_input.Get("h_"+varSS+"_tW");
    TH1D *emuSS_antitW = (TH1D*)f_input.Get("h_"+varSS+"_tbarW");


    // -- Merge ttbar samples -- //
    emu_ttbar->Add(emu_ttbarBackup);
    emu_ttbar->Add(emu_ttbar_M700to1000);
    emu_ttbar->Add(emu_ttbar_M1000toInf);
    emuSS_ttbar->Add(emuSS_ttbarBackup);
    emuSS_ttbar->Add(emuSS_ttbar_M700to1000);
    emuSS_ttbar->Add(emuSS_ttbar_M1000toInf);

    // DY M50 + M10to50
    emu_DYtautau->Add(emu_DYtautau_M10to50_v1);
    emu_DYtautau->Add(emu_DYtautau_M10to50_v2);
    emu_DYtautau->Add(emu_DYtautau_M10to50_ext1v1);
    emuSS_DYtautau->Add(emuSS_DYtautau_M10to50_v1);
    emuSS_DYtautau->Add(emuSS_DYtautau_M10to50_v2);
    emuSS_DYtautau->Add(emuSS_DYtautau_M10to50_ext1v1);

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


    // -- Consider average weight of top pt re-weight -- //
    emu_ttbar->Scale(1/AvgWeight);
    emuSS_ttbar->Scale(1/AvgWeight);


    // -- Set Color -- //
    emu_ttbar->SetFillColor(c1);
    emu_DYtautau->SetFillColor(c2);
    emu_diboson->SetFillColor(c4);
    emu_tW->SetFillColor(c7);

    emuSS_ttbar->SetFillColor(c1);
    emuSS_DYtautau->SetFillColor(c2);
    emuSS_diboson->SetFillColor(c4);
    emuSS_tW->SetFillColor(c7);


    // -- No Stats -- //
    emu_ttbar->SetStats(kFALSE);
    emu_DYtautau->SetStats(kFALSE);
    emu_diboson->SetStats(kFALSE);
    emu_tW->SetStats(kFALSE);

    emuSS_ttbar->SetStats(kFALSE);
    emuSS_DYtautau->SetStats(kFALSE);
    emuSS_diboson->SetStats(kFALSE);
    emuSS_tW->SetStats(kFALSE);


    // -- Set Data histogram -- //
    TH1D *emu_data = (TH1D*)f_input.Get("h_"+var+"_Data");
    TH1D *emuSS_data = (TH1D*)f_input.Get("h_"+varSS+"_Data");

    emu_data->SetMarkerStyle(20);
    emu_data->SetMarkerSize(0.8);
    emu_data->SetStats(kFALSE);

    emuSS_data->SetMarkerStyle(20);
    emuSS_data->SetMarkerSize(0.8);
    emuSS_data->SetStats(kFALSE);


    // -- Remove negative bins -- //
    removeNegativeBins( emu_DYtautau );
    removeNegativeBins( emuSS_DYtautau );


    // -- Rebin -- //
    if( var == "emu_mass_FineBin" ) {
        emu_diboson->Rebin(rebin);
        emu_tW->Rebin(rebin);
        emu_DYtautau->Rebin(rebin);
        emu_ttbar->Rebin(rebin);
        emu_data->Rebin(rebin);

        emuSS_diboson->Rebin(rebin);
        emuSS_tW->Rebin(rebin);
        emuSS_DYtautau->Rebin(rebin);
        emuSS_ttbar->Rebin(rebin);
        emuSS_data->Rebin(rebin);
    }


    // -- Stack histograms -- //
    THStack* emu_stackBkg = new THStack("emu_stackBkg","");
    emu_stackBkg->Add(emu_diboson);
    emu_stackBkg->Add(emu_tW);
    emu_stackBkg->Add(emu_DYtautau);
    emu_stackBkg->Add(emu_ttbar);

    THStack* emuSS_stackMC = new THStack("emuSS_stackMC","");
    emuSS_stackMC->Add(emuSS_diboson);
    emuSS_stackMC->Add(emuSS_tW);
    emuSS_stackMC->Add(emuSS_DYtautau);
    emuSS_stackMC->Add(emuSS_ttbar);


    // -- Legend -- //
    TLegend* legend = new TLegend(.75,.75,.95,.89);
    legend->AddEntry(emu_data,"Data");
    legend->AddEntry(emu_ttbar,"ttbar","F");
    legend->AddEntry(emu_DYtautau,"DY#tau#tau","F");
    legend->AddEntry(emu_tW,"tW+#bar{t}W","F");
    legend->AddEntry(emu_diboson,"VV","F");
    legend->SetBorderSize(0);  
    legend->SetFillStyle(0);  


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
//    emu_QCD->Smooth(Nsmooth,"R"); //Smoothing QCD
//    printf("Smoothing QCD : %d times\n", Nsmooth);


    // -- Make plot without QCD -- //
    if( var == "emu_mass" || var == "emu_mass_FineBin" ) {
        massplot( emu_stackBkg, emu_data, legend, var, 0 );
        massplot( emuSS_stackMC, emuSS_data, legend, varSS, 0 );
    }
    else {
        commonplot( emu_stackBkg, emu_data, legend, var, 0 );
        commonplot( emuSS_stackMC, emuSS_data, legend, varSS, 0 );
    }


    // -- Make plot with QCD -- //	
    THStack* emu_stackBkg_QCD = new THStack("emu_stackBkg_QCD","");
    emu_stackBkg_QCD->Add(emu_QCD);
    emu_stackBkg_QCD->Add(emu_diboson);
    emu_stackBkg_QCD->Add(emu_tW);
    emu_stackBkg_QCD->Add(emu_DYtautau);
    emu_stackBkg_QCD->Add(emu_ttbar);

    legend->AddEntry(emu_QCD,"QCD","F");

    if( var == "emu_mass" || var == "emu_mass_FineBin" )
        massplot( emu_stackBkg_QCD, emu_data, legend, var, 1 );
    else
        commonplot( emu_stackBkg_QCD, emu_data, legend, var, 1 );

/*
    // -- Check emu -- //
    cout<<"[Check emu events...!]"<<endl;
    cout<<endl;

    cout<<"data in emu: "<<emu_data->Integral()<<endl;
    cout<<"ttbar in emu: "<<emu_ttbar->Integral()<<endl;
    cout<<"DYtautau in emu: "<<emu_DYtautau->Integral()<<endl;
    cout<<"tW in emu: "<<emu_tW->Integral()<<endl;
    cout<<"diboson in emu: "<<emu_diboson->Integral()<<endl;
    cout<<"QCD in emu: "<<emu_QCD->Integral()<<endl;
    cout<<endl;

    double data = emu_data->Integral();
    cout<<"ttbar contribution in data: "<<(emu_ttbar->Integral())/data<<endl;
    cout<<"DYtautau contribution in data: "<<(emu_DYtautau->Integral())/data<<endl;
    cout<<"tW contribution in data: "<<(emu_tW->Integral())/data<<endl;
    cout<<"diboson contribution in data: "<<(emu_diboson->Integral())/data<<endl;
    cout<<"QCD contribution in data: "<<(emu_QCD->Integral())/data<<endl;
    cout<<endl;

    // -- Check emuSS -- //
    cout<<"[Check emuSS events...!]"<<endl;
    cout<<endl;

    cout<<"data in emuSS: "<<emuSS_data->Integral()<<endl;
    cout<<"ttbar in emuSS: "<<emuSS_ttbar->Integral()<<endl;
    cout<<"DYtautau in emuSS: "<<emuSS_DYtautau->Integral()<<endl;
    cout<<"tW in emuSS: "<<emuSS_tW->Integral()<<endl;
    cout<<"diboson in emuSS: "<<emuSS_diboson->Integral()<<endl;
    cout<<endl;

    double dataSS = emuSS_data->Integral();
    cout<<"ttbar contribution in data: "<<(emuSS_ttbar->Integral())/dataSS<<endl;
    cout<<"DYtautau contribution in data: "<<(emuSS_DYtautau->Integral())/dataSS<<endl;
    cout<<"tW contribution in data: "<<(emuSS_tW->Integral())/dataSS<<endl;
    cout<<"diboson contribution in data: "<<(emuSS_diboson->Integral())/dataSS<<endl;
    cout<<endl;
*/
}

void removeNegativeBins( TH1D* hist ) {

    for(int i=0; i<hist->GetNbinsX(); i++) {
        if(hist->GetBinContent(i+1)<0) {
            hist->SetBinContent(i+1,0);
            hist->SetBinError(i+1,0);
        }
    }   

}

void massplot( THStack* h_mc, TH1D* h_data, TLegend* leg, TString plotname, Int_t QCDoption ) {

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
    if( plotname == "emu_mass_FineBin" || plotname == "emuSS_mass_FineBin" )
        h_mass_stack_1 = new TH1D("h_mass_stack_1","",10000/rebin,0,10000);
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
    if( plotname == "emu_mass_FineBin" || plotname == "emuSS_mass_FineBin" )
        h_mass_stack_1->GetXaxis()->SetRangeUser(100,3000);
    else
        h_mass_stack_1->GetXaxis()->SetRangeUser(15,3000);

    h_mc->SetHistogram(h_mass_stack_1);
    h_mc->SetTitle("");
    h_mc->SetMaximum(100000);
    h_mc->SetMinimum(1);
    h_mc->Draw("hist");

    if( plotname == "emu_mass_FineBin" || plotname == "emuSS_mass_FineBin" )
        h_data->GetXaxis()->SetRangeUser(100,3000);
    else
        h_data->GetXaxis()->SetRangeUser(15,3000);
    h_data->SetMaximum(100000);
    h_data->SetMinimum(1);
    h_data->SetMarkerStyle(20);
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

    // -- Bottom pad - //
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
    if( plotname == "emu_mass_FineBin" || plotname == "emuSS_mass_FineBin" )
        sum_mc = new TH1D("sum_mc","",10000/rebin,0,10000);
    sum_mc->Merge(list);
    sum_mc->Sumw2();

    TH1D *ratio_mass__7 = (TH1D*)h_data->Clone("ratio_mass__7");
    ratio_mass__7->Divide(h_data, sum_mc, 1.0, 1.0, "B");

    ratio_mass__7->SetMarkerStyle(20);
    ratio_mass__7->SetMarkerSize(0.8);
    if( plotname == "emu_mass" || plotname == "emu_mass_FineBin" )
        ratio_mass__7->GetXaxis()->SetTitle("e#mu mass [GeV]");
    else if( plotname == "emuSS_mass" || plotname == "emuSS_mass_FineBin" )
        ratio_mass__7->GetXaxis()->SetTitle("e#muSS mass [GeV]");
    if( plotname == "emu_mass_FineBin" || plotname == "emuSS_mass_FineBin" )
        ratio_mass__7->GetXaxis()->SetRangeUser(100,3000);
    else
        ratio_mass__7->GetXaxis()->SetRangeUser(15,3000);
    ratio_mass__7->GetXaxis()->SetMoreLogLabels();
    ratio_mass__7->GetXaxis()->SetLabelFont(42);
    ratio_mass__7->GetXaxis()->SetLabelOffset(0.007);
    ratio_mass__7->GetXaxis()->SetLabelSize(0.12);
    ratio_mass__7->GetXaxis()->SetTitleSize(0.13);
    ratio_mass__7->GetXaxis()->SetTitleOffset(1.4);
    ratio_mass__7->GetXaxis()->SetTitleFont(42);
    ratio_mass__7->GetYaxis()->SetTitle("Data/MC");
    if( plotname == "emu_mass_FineBin" || plotname == "emuSS_mass_FineBin" )
        ratio_mass__7->GetYaxis()->SetRangeUser(0.5,1.5);
    else
        ratio_mass__7->GetYaxis()->SetRangeUser(0.7,1.3);
        //ratio_mass__7->GetYaxis()->SetRangeUser(0.9,1.1);
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

    if( QCDoption == 0 )
        Mass1->SaveAs(outputname+"_"+plotname+".pdf");
    else if( QCDoption == 1 )
        plotname = plotname+"_withQCD";
        Mass1->SaveAs(outputname+"_"+plotname+".pdf");
}

void commonplot( THStack* h_mc, TH1D* h_data, TLegend* leg, TString plotname, Int_t QCDoption ) {

    double nbins; double xmin; double xmax;

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
    //InvariantMass->SetLogx();
    InvariantMass->SetLogy();
    InvariantMass->SetGridx();
    InvariantMass->SetGridy();
    InvariantMass->SetTickx();
    InvariantMass->SetTicky();
    InvariantMass->SetBottomMargin(0.3130424);
    InvariantMass->SetFrameBorderMode(0);

    if( plotname == "mu_pT" ) {nbins = 300; xmin = 0; xmax = 600;}
    if( plotname == "mu_eta" ) {nbins = 56; xmin = -2.8; xmax = 2.8;}
    if( plotname == "mu_phi" ) {nbins = 50; xmin = -3.5; xmax = 3.5;}
    if( plotname == "muSS_pT" ) {nbins = 300; xmin = 0; xmax = 600;}
    if( plotname == "muSS_eta" ) {nbins = 56; xmin = -2.8; xmax = 2.8;}
    if( plotname == "muSS_phi" ) {nbins = 50; xmin = -3.5; xmax = 3.5;}
    if( plotname == "ele_pT" ) {nbins = 300; xmin = 0; xmax = 600;}
    if( plotname == "ele_eta" ) {nbins = 56; xmin = -2.8; xmax = 2.8;}
    if( plotname == "ele_phi" ) {nbins = 50; xmin = -3.5; xmax = 3.5;}
    if( plotname == "eleSS_pT" ) {nbins = 300; xmin = 0; xmax = 600;}
    if( plotname == "eleSS_eta" ) {nbins = 56; xmin = -2.8; xmax = 2.8;}
    if( plotname == "eleSS_phi" ) {nbins = 50; xmin = -3.5; xmax = 3.5;}

    TH1D *h_mass_stack_1 = new TH1D("h_mass_stack_1","", nbins, xmin, xmax);
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

    h_mc->SetHistogram(h_mass_stack_1);
    h_mc->SetTitle("");
    h_mc->SetMaximum(100000);
    h_mc->SetMinimum(1);
    h_mc->Draw("hist");

    h_data->SetMaximum(100000);
    h_data->SetMinimum(1);
    h_data->SetMarkerStyle(20);
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
    //padc1_2->SetLogx();
    padc1_2->SetTopMargin(0.07176078);
    padc1_2->SetBottomMargin(0.4233841);
    padc1_2->SetLeftMargin(L);
    padc1_2->SetRightMargin(R);
    padc1_2->SetFrameBorderMode(0);

    // -- ratio plot -- //
    TList *list = (TList*)h_mc->GetHists();

    if( plotname == "mu_eta" ) {nbins = 100; xmin = -5; xmax = 5;}
    if( plotname == "mu_phi" ) {nbins = 100; xmin = -5; xmax = 5;}
    if( plotname == "muSS_eta" ) {nbins = 100; xmin = -5; xmax = 5;}
    if( plotname == "muSS_phi" ) {nbins = 100; xmin = -5; xmax = 5;}
    if( plotname == "ele_eta" ) {nbins = 100; xmin = -5; xmax = 5;}
    if( plotname == "ele_phi" ) {nbins = 100; xmin = -5; xmax = 5;}
    if( plotname == "eleSS_eta" ) {nbins = 100; xmin = -5; xmax = 5;}
    if( plotname == "eleSS_phi" ) {nbins = 100; xmin = -5; xmax = 5;}

    TH1D *sum_mc = new TH1D("sum_mc", "", nbins, xmin, xmax);
    sum_mc->Merge(list);
    sum_mc->Sumw2();

    TH1D *ratio_mass__7 = (TH1D*)h_data->Clone("ratio_mass__7");
    ratio_mass__7->Divide(h_data, sum_mc, 1.0, 1.0, "B");

    ratio_mass__7->SetMarkerStyle(20);
    ratio_mass__7->SetMarkerSize(0.8);
    ratio_mass__7->GetXaxis()->SetTitle(plotname);
    if( plotname == "mu_eta" || plotname == "muSS_eta" || plotname == "ele_eta" || plotname == "eleSS_eta" )
        ratio_mass__7->GetXaxis()->SetRangeUser(-2.8,2.8);
    if( plotname == "mu_phi" || plotname == "muSS_phi" || plotname == "ele_phi" || plotname == "eleSS_phi" )
        ratio_mass__7->GetXaxis()->SetRangeUser(-3.5,3.5);
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

    if( QCDoption == 0 )
        Mass1->SaveAs(outputname+"_"+plotname+".pdf");
    else if( QCDoption == 1 )
        plotname = plotname+"_withQCD";
        Mass1->SaveAs(outputname+"_"+plotname+".pdf");
}

