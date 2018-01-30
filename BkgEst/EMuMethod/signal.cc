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

//const double ScaleFactor = 1.0338;
const double ScaleFactor = 1.0;
const double lumi = 35867.0;

const int binnum = 43;
const double bins[44] = {15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 64, 68, 72, 76, 81, 86, 91, 96, 101, 106, 110, 115, 120, 126, 133, 141, 150, 160, 171, 185, 200, 220, 243, 273, 320, 380, 440, 510, 600, 700, 830, 1000, 1500, 3000};

using namespace std;


void massplot( THStack* h_mc, TH1D* h_data, TLegend* leg, TString plotname );


// -- Channel : select "MuMu" or "EE" -- //
TString channel = "MuMu";
//TString channel = "EE";

TString inputname1 = "./root/test/ROOTFile_Muon_Channel_1.root"; // merged root file of all output ones from dimuon and emu event selection
TString inputname2 = "./result/estimatedBkg/emu_TopPtReweighting_on_"+channel+".root"; // output root file of "estimateBkg.cc"
TString outputname = "./result/estimatedBkg/signal_"+channel;


void signal() {

    cout << "==============================================================================" << endl;
    cout << "input file 1 : " << inputname1 << endl;
    cout << "input file 2 : " << inputname2 << endl;
    cout << "output file  : " << outputname << endl;
    cout << "==============================================================================" << endl;
    cout << endl;

    TString var = "mass";
    
    // -- Set MC histogram -- //
    TFile f_input(inputname1, "read");
    TH1D *h_DY_M10to50_v1 = (TH1D*)f_input.Get("h_"+var+"_DY"+channel+"_M10to50_v1");
    TH1D *h_DY_M10to50_v2 = (TH1D*)f_input.Get("h_"+var+"_DY"+channel+"_M10to50_v2");
    TH1D *h_DY_M10to50_ext1v1 = (TH1D*)f_input.Get("h_"+var+"_DY"+channel+"_M10to50_ext1v1");
    TH1D *h_DY_M50to100 = (TH1D*)f_input.Get("h_"+var+"_DY"+channel+"_M50to100");
    TH1D *h_DY_M100to200 = (TH1D*)f_input.Get("h_"+var+"_DY"+channel+"_M100to200");
    TH1D *h_DY_M100to200_ext = (TH1D*)f_input.Get("h_"+var+"_DY"+channel+"_M100to200_ext");
    TH1D *h_DY_M200to400 = (TH1D*)f_input.Get("h_"+var+"_DY"+channel+"_M200to400");
    TH1D *h_DY_M400to500 = (TH1D*)f_input.Get("h_"+var+"_DY"+channel+"_M400to500");
    TH1D *h_DY_M500to700 = (TH1D*)f_input.Get("h_"+var+"_DY"+channel+"_M500to700");
    TH1D *h_DY_M700to800 = (TH1D*)f_input.Get("h_"+var+"_DY"+channel+"_M700to800");
    TH1D *h_DY_M800to1000 = (TH1D*)f_input.Get("h_"+var+"_DY"+channel+"_M800to1000");
    TH1D *h_DY_M1000to1500 = (TH1D*)f_input.Get("h_"+var+"_DY"+channel+"_M1000to1500");
    TH1D *h_DY_M1500to2000 = (TH1D*)f_input.Get("h_"+var+"_DY"+channel+"_M1500to2000");
    TH1D *h_DY_M2000to3000 = (TH1D*)f_input.Get("h_"+var+"_DY"+channel+"_M2000to3000");
    TH1D *h_WZ = (TH1D*)f_input.Get("h_"+var+"_WZ");
    TH1D *h_ZZ = (TH1D*)f_input.Get("h_"+var+"_ZZ");
    TH1D *h_WJets = (TH1D*)f_input.Get("h_"+var+"_WJetsToLNu");
    TH1D *h_WJets_ext = (TH1D*)f_input.Get("h_"+var+"_WJetsToLNu_ext");

    // Data-driven Bkg
    TFile f_input2(inputname2, "read");
    TH1D *h_ttbar = (TH1D*)f_input2.Get("ttbar");
    TH1D *h_DYtautau = (TH1D*)f_input2.Get("DYtautau");
    TH1D *h_WW = (TH1D*)f_input2.Get("WW");
    TH1D *h_tW = (TH1D*)f_input2.Get("tW");

    // -- Single Scale Factor -- //
    h_DY_M100to200->Scale(ScaleFactor);
    h_DY_M100to200_ext->Scale(ScaleFactor);
    h_DY_M200to400->Scale(ScaleFactor);
    h_DY_M400to500->Scale(ScaleFactor);
    h_DY_M500to700->Scale(ScaleFactor);
    h_DY_M700to800->Scale(ScaleFactor);
    h_DY_M800to1000->Scale(ScaleFactor);
    h_DY_M1000to1500->Scale(ScaleFactor);
    h_DY_M1500to2000->Scale(ScaleFactor);
    h_DY_M2000to3000->Scale(ScaleFactor);

    // Merge All DY sig
    TH1D *h_DY = (TH1D*)h_DY_M10to50_v1->Clone();
    h_DY->Add(h_DY_M10to50_v2);
    h_DY->Add(h_DY_M10to50_ext1v1);
    h_DY->Add(h_DY_M50to100);
    h_DY->Add(h_DY_M100to200);
    h_DY->Add(h_DY_M100to200_ext);
    h_DY->Add(h_DY_M200to400);
    h_DY->Add(h_DY_M400to500);
    h_DY->Add(h_DY_M500to700);
    h_DY->Add(h_DY_M700to800);
    h_DY->Add(h_DY_M800to1000);
    h_DY->Add(h_DY_M1000to1500);
    h_DY->Add(h_DY_M1500to2000);
    h_DY->Add(h_DY_M2000to3000);

    // WW + WZ + ZZ
    TH1D* h_diboson = h_WW;
    h_diboson->Add(h_WZ);
    h_diboson->Add(h_ZZ);

    // WJets + WJets_ext
    h_WJets->Add(h_WJets_ext);


    // -- Normalization -- //
    //h_WJets->Scale(1/86731698.0);


    // -- Merge processes: DY, Top, EW, QCD -- //
    TH1D* h_Top = h_ttbar;
    h_Top->Add(h_tW);
    TH1D* h_EW = h_DYtautau;
    h_EW->Add(h_diboson);
    TH1D* h_Fakes = h_WJets;


    // -- Set Fill Color -- //
    h_DY->SetFillColor(kOrange-2);
    h_Top->SetFillColor(kRed+2);
    h_EW->SetFillColor(kOrange+10);
    h_Fakes->SetFillColor(kViolet-5);


    // -- Set Line Color -- //
    h_DY->SetLineColor(kOrange+3);
    h_Top->SetLineColor(kRed+4);
    h_EW->SetLineColor(kOrange+3);
    h_Fakes->SetLineColor(kOrange+3);


    // -- No Stats -- //
    h_DY->SetStats(kFALSE);
    h_Top->SetStats(kFALSE);
    h_EW->SetStats(kFALSE);
    h_Fakes->SetStats(kFALSE);


    // -- Set Data histogram -- //
    TH1D *h_data = (TH1D*)f_input.Get("h_"+var+"_Data");
    h_data->SetMarkerStyle(20);
    h_data->SetMarkerSize(0.8);
    h_data->SetStats(kFALSE);


    // -- Stack histograms -- //
    THStack* h_stack = new THStack("h_stack","");
    h_stack->Add(h_Fakes);
    h_stack->Add(h_EW);
    h_stack->Add(h_Top);
    h_stack->Add(h_DY);


    // -- Legend -- //
    TLegend* legend = new TLegend(.65,.75,.85,.89);
    legend->AddEntry(h_data,"Data");
    legend->AddEntry(h_DY,"DY"+channel,"F");
    legend->AddEntry(h_Top,"t#bar{t}+tW+#bar{t}W","F");
    legend->AddEntry(h_EW,"EW","F");
    legend->AddEntry(h_Fakes,"Fakes(WJets)","F");
    legend->SetBorderSize(0);  
    legend->SetFillStyle(0);  


    // -- Set signal -- //
    TH1D *h_sig = (TH1D*)h_data->Clone();
    h_sig->Add(h_Top,-1.0); 
    h_sig->Add(h_EW,-1.0); 
    h_sig->Add(h_Fakes,-1.0); 

    THStack* h_dy = new THStack("h_dy", "");
    h_dy->Add(h_DY);

    // output file
    TFile* g = new TFile(outputname+".root","RECREATE");
    h_sig->SetName("h_sig");
    h_sig->Write();
    g->Close();

    // legend
    TLegend* leg = new TLegend(.65,.75,.85,.89);
    leg->AddEntry(h_sig,"Signal(Data)");
    leg->AddEntry(h_DY,"MC(DY)","F");
    leg->SetBorderSize(0);  
    leg->SetFillStyle(0);  


    // -- Make plot -- //
    massplot( h_stack, h_data, legend, "All" );
    massplot( h_dy, h_sig, leg, "signal" );
}

void massplot( THStack* h_mc, TH1D* h_data, TLegend* leg, TString plotname ) {

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
    if( plotname.Contains("Zpeak") == kFALSE )
        InvariantMass->SetLogx();
    InvariantMass->SetLogy();
    InvariantMass->SetGridx();
    InvariantMass->SetGridy();
    InvariantMass->SetTickx();
    InvariantMass->SetTicky();
    InvariantMass->SetBottomMargin(0.3130424);
    InvariantMass->SetFrameBorderMode(0);

    TH1D *h_mass_stack_1;
    if( plotname.Contains("Zpeak") )
        h_mass_stack_1 = new TH1D("h_mass_stack_1","",60,60,120);
    else
        h_mass_stack_1 = new TH1D("h_mass_stack_1","",binnum,bins);
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
    if( plotname.Contains("Zpeak") )
        h_mass_stack_1->GetXaxis()->SetRangeUser(60,120);
    else
        h_mass_stack_1->GetXaxis()->SetRangeUser(15,3000);

    h_mc->SetHistogram(h_mass_stack_1);
    h_mc->SetTitle("");
    h_mc->SetMaximum(50000000);
    h_mc->SetMinimum(1);
    h_mc->Draw("hist");

    if( plotname.Contains("Zpeak") )
        h_data->GetXaxis()->SetRangeUser(60,120);
    else
        h_data->GetXaxis()->SetRangeUser(15,3000);
    h_data->SetMaximum(50000000);
    h_data->SetMinimum(1);
    h_data->SetMarkerStyle(20);
    h_data->Draw("same p e");

    TList *list = (TList*)h_mc->GetHists();

    TH1D *sum_mc;
    if( plotname.Contains("Zpeak") )
        sum_mc = new TH1D("sum_mc","",60,60,120);
    else
        sum_mc = new TH1D("sum_mc","",binnum,bins);
    sum_mc->Merge(list);
    sum_mc->Sumw2();

    // -- Grid -- //
    TLine grid_;
    grid_.SetLineColor(kGray+2);
    grid_.SetLineStyle(kSolid);
    for( size_t ii=0; ii<44; ii++ )
        grid_.DrawLine(bins[ii],0.2,bins[ii],sum_mc->GetBinContent(ii+1));

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
    if( plotname.Contains("Zpeak") == kFALSE )
        padc1_2->SetLogx();
    padc1_2->SetTopMargin(0.07176078);
    padc1_2->SetBottomMargin(0.4233841);
    padc1_2->SetLeftMargin(L);
    padc1_2->SetRightMargin(R);
    padc1_2->SetFrameBorderMode(0);

    // -- ratio plot -- //
    TH1D *ratio_mass__7 = (TH1D*)h_data->Clone("ratio_mass__7");
    ratio_mass__7->Divide(h_data, sum_mc, 1.0, 1.0, "B");

    ratio_mass__7->SetMarkerStyle(20);
    ratio_mass__7->SetMarkerSize(0.8);
    if( channel == "MuMu" )
        ratio_mass__7->GetXaxis()->SetTitle("M(#mu#mu) [GeV]");
    else if( channel == "EE" )
        ratio_mass__7->GetXaxis()->SetTitle("M(ee) [GeV]");
    if( plotname.Contains("Zpeak") )
        ratio_mass__7->GetXaxis()->SetRangeUser(60,120);
    else
        ratio_mass__7->GetXaxis()->SetRangeUser(15,3000);
    ratio_mass__7->GetXaxis()->SetMoreLogLabels();
    ratio_mass__7->GetXaxis()->SetLabelFont(42);
    ratio_mass__7->GetXaxis()->SetLabelOffset(0.007);
    ratio_mass__7->GetXaxis()->SetLabelSize(0.12);
    ratio_mass__7->GetXaxis()->SetTitleSize(0.13);
    ratio_mass__7->GetXaxis()->SetTitleOffset(1.4);
    ratio_mass__7->GetXaxis()->SetTitleFont(42);
    if( plotname == "All" )
        ratio_mass__7->GetYaxis()->SetTitle("Data/(Sig+Bkg)");
    else if( plotname == "signal" )
        ratio_mass__7->GetYaxis()->SetTitle("Sig/MC");
    //ratio_mass__7->GetYaxis()->SetRangeUser(0.7,1.3);
    ratio_mass__7->GetYaxis()->SetRangeUser(0.9,1.1);
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

