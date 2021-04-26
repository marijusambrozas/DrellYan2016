// Class for drawing histograms with ratio plots. Created by Marijus Ambrozas

#ifndef myRatioPlot_t_h
#define myRatioPlot_t_h

#include <TROOT.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TString.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TLine.h>
#include <THStack.h>
#include <TF1.h>
#include <TLatex.h>

class myRatioPlot_t {
public :
   TH1D *h1_data, *h1_dataovermc, *h1_fullunc, *h1_one;
   THStack *s_stackedProcesses;
   TLegend *legend;
   TCanvas *canvas;
   TPad *pad1, *pad2;
   TText *text;
   TString CanvasName;
   Double_t x1, x2;
   Int_t legendSet, legendEntries, PlotsSet, systSet;
   int autoCleanUp;

   myRatioPlot_t(TString cName=" ", const THStack *stack=NULL, const TH1D *hist=NULL);
   myRatioPlot_t(TString cName, const TH1D *stack=NULL, const TH1D *hist=NULL);
   virtual ~myRatioPlot_t();

   void SetAutoCleanUp(int yes=1) { autoCleanUp=yes; }
   int GetAutoCleanUp() const { return autoCleanUp; }

   void Draw(Double_t ymin, Double_t ymax, UInt_t logX=0, TString option="HIST", TString yAxisName="Number of events", UInt_t logY=1);
   void DrawWithExpectationBar(TH1D *h_errs, Double_t ymin, Double_t ymax, UInt_t logX=0, TString option="HIST", TString yAxisName="Number of events", UInt_t logY=1); // NOT WORKING PROPERLY!!!
   void SetPlots(TString xAxisName, Double_t xmin, Double_t xmax, TString DataMCname="Data/MC", Double_t *statunc=NULL);
   void SetLegend(Double_t xb=0.7, Double_t yb=0.5, Double_t xe=0.9, Double_t ye=0.9);
   void AddLegendEntry(TH1D *h, TString name, TString option);
   void ImportLegend(TLegend *newLegend, Int_t overwrite=0);
   void SetSystematics(Double_t *dataSystematics=NULL, Double_t *MCsystematics=NULL, Double_t *dataOverMCsystematics=NULL);
   void DrawOnTop(TH1D* h1_onTop, TString option="");
};

#endif

#ifdef myRatioPlot_t_cxx

myRatioPlot_t::myRatioPlot_t(TString cName, const THStack *stack, const TH1D *hist) :
    h1_data(0), h1_dataovermc(0), h1_fullunc(0), h1_one(0),
    s_stackedProcesses(0), legend(0),
    canvas(0), pad1(0), pad2(0),
      x1(0.), x2(0.), legendSet(0), legendEntries(0), PlotsSet(0), systSet(0),
      autoCleanUp(1)
{
    if (!cName || !stack || !hist) return;
    h1_data = ((TH1D*)hist->Clone("h1_data"));
    s_stackedProcesses = new THStack(*stack);
    CanvasName = cName;
    gStyle->SetOptStat(0);
}

myRatioPlot_t::myRatioPlot_t(TString cName, const TH1D *stack, const TH1D *hist) :
    h1_data(0), h1_dataovermc(0), h1_fullunc(0), h1_one(0),
    s_stackedProcesses(0), legend(0),
    canvas(0), pad1(0), pad2(0),
      x1(0.), x2(0.), legendSet(0), legendEntries(0), PlotsSet(0), systSet(0),
      autoCleanUp(1)
{
    if (!cName || !stack || !hist) return;
    h1_data = ((TH1D*)hist->Clone("h1_data"));
    s_stackedProcesses = new THStack("stack", "");
    TH1D* hclone = ((TH1D*)stack->Clone(stack->GetName()+TString("_clone")));
    hclone->SetDirectory(0);
    s_stackedProcesses->Add(hclone);
    CanvasName = cName;
    gStyle->SetOptStat(0);
}

myRatioPlot_t::~myRatioPlot_t()
{
   if (!h1_data && !s_stackedProcesses) return;
   if (autoCleanUp) {
     delete h1_data;
     delete h1_dataovermc;
     delete h1_fullunc;
     delete h1_one;
     delete s_stackedProcesses;
     delete legend;
     delete canvas;
     delete pad1;
     delete pad2;
   }
}


#undef myRatioPlot_t_cxx
#endif // #ifdef myRatioPlot_t_cxx
