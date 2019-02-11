#define myRatioPlot_t_cxx

#include "myRatioPlot_t.h"
#include <TROOT.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TString.h>
#include <TGraph.h>
#include <TLegend.h>
#include <THStack.h>
#include <TLine.h>
#include <TF1.h>
#include <TLatex.h>

void myRatioPlot_t::SetPlots(TString xAxisName, Double_t xmin, Double_t xmax)
{
    x1=xmin; x2=xmax;
    // Calculating data/MC
    h1_dataovermc = ((TH1D*)h1_data->Clone("h1_dataovermc"));
    TH1D* h1_deno_invm=((TH1D*)(s_stackedProcesses->GetStack()->Last()));
    h1_dataovermc->Divide(h1_deno_invm);

    //Setting hist options
    h1_data->SetMarkerStyle(kFullDotLarge);
    h1_data->SetMarkerColor(kBlack);

    h1_dataovermc->SetLineColorAlpha(kBlack, 0);
    h1_dataovermc->SetMarkerStyle(kFullDotLarge);
    h1_dataovermc->SetMarkerColor(kBlack);
    h1_dataovermc->SetTitle(" ");
    h1_dataovermc->GetXaxis()->SetTitle(xAxisName);
    h1_dataovermc->GetYaxis()->SetTitle("Data/MC");
//    h1_dataovermc->GetYaxis()->SetTitle("Eksp./Mod.");
//    h1_dataovermc->GetYaxis()->SetTitle("Eksp./Ivert.");
//    h1_dataovermc->GetYaxis()->SetTitle("Ivert./Mod.");
    h1_dataovermc->GetXaxis()->SetNoExponent(1);
    h1_dataovermc->GetXaxis()->SetMoreLogLabels(1);    
    if (h1_dataovermc->GetBinContent(h1_dataovermc->GetMinimumBin())>0.3 && h1_dataovermc->GetBinContent(h1_dataovermc->GetMaximumBin())<2) {
        h1_dataovermc->GetYaxis()->SetRangeUser(h1_dataovermc->GetBinContent(h1_dataovermc->GetMinimumBin())-0.05,
                                                h1_dataovermc->GetBinContent(h1_dataovermc->GetMaximumBin())+0.05);
    }
    else{
        Double_t minv=1000, maxv=-1000;
        for(Int_t n=1; n<h1_dataovermc->GetSize()-1; n++) {
            if (h1_dataovermc->GetBinContent(n)>maxv && h1_dataovermc->GetBinContent(n)<2.5) maxv=h1_dataovermc->GetBinContent(n);
            else if (h1_dataovermc->GetBinContent(n)<minv && h1_dataovermc->GetBinContent(n)>0.3) minv=h1_dataovermc->GetBinContent(n);
        }
        h1_dataovermc->GetYaxis()->SetRangeUser(minv-0.05, maxv+0.05);
    }
    // ######### DELETE THIS WHEN NOT NEEDED ########
    h1_dataovermc->GetYaxis()->SetRangeUser(0.75, 1.25);
    // ##############################################

    h1_dataovermc->GetXaxis()->SetRangeUser(xmin, xmax);
    h1_dataovermc->GetXaxis()->SetLabelSize(0.13);
    h1_dataovermc->GetXaxis()->SetTitleSize(0.17);
    h1_dataovermc->GetXaxis()->SetTitleOffset(0.75);
    h1_dataovermc->GetYaxis()->SetLabelSize(0.11);
    h1_dataovermc->GetYaxis()->SetTitleSize(0.11);
    h1_dataovermc->GetYaxis()->SetTitleOffset(0.3);
    h1_dataovermc->GetYaxis()->CenterTitle();
    h1_dataovermc->GetYaxis()->SetNdivisions(8);
    h1_dataovermc->GetXaxis()->SetTickLength(0.01);
    h1_dataovermc->GetYaxis()->SetTickLength(0.02);

    PlotsSet++;
}

void myRatioPlot_t::Draw(Double_t ymin, Double_t ymax, UInt_t logX)
{
    if(PlotsSet==1) {
//        canvas = new TCanvas(CanvasName, CanvasName, 1000, 1000);
        canvas = new TCanvas(CanvasName, CanvasName, 850, 1000);
//        canvas = new TCanvas(CanvasName, CanvasName, 1500, 1000);
        TPad* pad1 = new TPad("pad1", "pad1", 0, 0.255, 1, 1);
        pad1->SetBottomMargin(0.001);
        pad1->Draw();
        pad1->cd();
        s_stackedProcesses->Draw("HIST");
        s_stackedProcesses->GetYaxis()->SetTitle("Number of events");
//        s_stackedProcesses->GetYaxis()->SetTitle("Ivykiu skai#check{c}ius");
        s_stackedProcesses->GetYaxis()->SetTitleOffset(1);
        s_stackedProcesses->GetYaxis()->SetTitleSize(0.04);
        s_stackedProcesses->GetXaxis()->SetNoExponent(1);
        s_stackedProcesses->SetMinimum(ymin);
        s_stackedProcesses->SetMaximum(ymax);
        if(logX) s_stackedProcesses->GetXaxis()->SetMoreLogLabels(1);
        s_stackedProcesses->GetXaxis()->SetRangeUser(x1, x2);
        s_stackedProcesses->GetYaxis()->SetLabelSize(0.04);
        s_stackedProcesses->Draw("sameaxis");
        h1_data->Draw("same E");
        h1_data->Draw("sameaxis");
        h1_data->SetDirectory(0);

        if(legendSet==1) {
            if(legendEntries) legend->Draw();
            else std::cout << "myRatioPlot: Legend has " << legendEntries << " entries and cannot be drawn.\n";
        }
        else std::cout << "myRatioPlot: Legend was not set properly.\n";

        if(logX) pad1->SetLogx();
        pad1->SetLogy();
        pad1->SetTickx(1);
        pad1->SetTicky(1);
        pad1->SetGridx(1);
        pad1->SetGridy(1);
//        text = new TText (.1, .91, "CMS Preliminary");
//        text->SetTextAlign(11);
//        text->SetTextSize(0.05);
//        text->SetNDC(true);
//        text->Draw();
        pad1->Update();
        canvas->cd();

        pad2 = new TPad("pad2", "pad2", 0, 0, 1, 0.24);
        pad2->SetTopMargin(0);
        pad2->Draw();
        pad2->cd();
        h1_dataovermc->Draw("E");
        h1_dataovermc->SetDirectory(0);        
        if(logX) l_one = new TLine(x1, 1, x2+5, 1);
        else l_one = new TLine (x1, 1, x2, 1);
        l_one->SetLineColor(kRed);
        l_one->SetLineWidth(2);
        l_one->Draw("same");
        l_one->Draw("sameaxis");
        h1_dataovermc->Draw("sameaxis");
        if(logX) pad2->SetLogx();
        pad2->SetTickx(1);
        pad2->SetTicky(1);
        pad2->SetGridx(1);
        pad2->SetGridy(1);
        pad2->SetBottomMargin(0.3);
        pad2->Update();
        canvas->cd();

        canvas->Update();
    }
    else { std::cout << "myRatioPlot: Plots were not set properly.\n"; return; }
}

void myRatioPlot_t::SetLegend(Double_t xb, Double_t yb, Double_t xe, Double_t ye) {
    if(xb>=0 && yb>=0 && xe<=1 && ye<=1 && xb<xe && yb<ye) {
        legend = new TLegend(xb, yb, xe, ye);
        legendSet++;
    }
}

void myRatioPlot_t::AddLegendEntry(TH1D *h, TString name, TString option) {
    if (legendSet==1) {
        legend->AddEntry(h, name, option);
        legendEntries++;
    }
    else return;
}


