#include <Include/PlotTools.h>

TH2D* CalculateFractionPerBin(TH2D *h_nEvents, TH1D* h_Truth);
TH1D* MakeDiagonalHist( TH2D* h2D);

void DrawPlot_2015BinMigration()
{
	TString fileName = "/Users/KyeongPil_Lee/Physics/DYAnalysis_76X_LumiUpdate/Include/Results_ROOTFiles_76X/ROOTFile_Histogram_ResponseM_1D_aMCNLO_IsoMu20_OR_IsoTkMu20.root";

	TH2D* h2D_nEvent = Get_Hist_2D( fileName, "h_RespM_RooUnfold");
	TH1D* h_Truth = Get_Hist( fileName, "h_Truth_RooUnfold");

	TH2D* h_RespM = CalculateFractionPerBin( h2D_nEvent, h_Truth );
	h_RespM->Draw("COLZ");
	TH1D* h_diag = MakeDiagonalHist( h_RespM );

	// -- drawing -- //
	SampleInfo* sampleInfo = new SampleInfo();
	sampleInfo->SetColor( kBlack );

	HistInfo* histInfo = new HistInfo("m [GeV]", "Diagonal term");
	histInfo->SetYRange( 0, 1.2 );

	TH1Ext *hext_diag = new TH1Ext(sampleInfo, histInfo, h_diag );

	TCanvas *c;
	TString CanvasName = "Local/c_diag_2015_MuonChannel";
	SetCanvas_Square( c, CanvasName, 1, 0 );

	c->cd();

	hext_diag->DrawAndSet( "HISTLPSAME" );

	TLatex latex;
	Latex_Simulation( latex );

	TString channelInfo = "#mu channel (2015 analysis)";
	latex.DrawLatexNDC(0.16, 0.91, "#font[62]{#scale[0.6]{"+channelInfo+"}}");

	// c->SaveAs(".png");
	c->SaveAs(".pdf");
}

TH2D* CalculateFractionPerBin(TH2D *h_nEvents, TH1D* h_Truth)
{
	TH2D* h_Response = (TH2D*)h_nEvents->Clone();

	Int_t nBins = h_Truth->GetNbinsX();
	for(Int_t i_genbin=0; i_genbin <= nBins+1; i_genbin++) // -- Include under/overflow -- //
	{
		Double_t SumN_truth = h_Truth->GetBinContent(i_genbin);

		for(Int_t i_recobin=0; i_recobin <= nBins+1; i_recobin++) // -- Include under/overflow -- //
		{
			// Double_t nEvent = h_nEvents->GetBinContent(i_genbin, i_recobin);
			Double_t nEvent = h_nEvents->GetBinContent(i_recobin, i_genbin);

			Double_t fraction = 0;
			if( SumN_truth == 0 )
				fraction = 0;
			else
				fraction = nEvent / SumN_truth;
			
			if( fraction < 0 && fabs(fraction) < 1e-3 ) fraction = 0;

			h_Response->SetBinContent( i_genbin, i_recobin, fraction );
		}
	}

	return h_Response;
}

TH1D* MakeDiagonalHist( TH2D* h2D)
{
	const Int_t nMassBin = 43;
	Double_t MassBinEdges[nMassBin+1] = {15, 20, 25, 30, 35, 40, 45, 50, 55, 60,
										 64, 68, 72, 76, 81, 86, 91, 96, 101, 106,
										 110, 115, 120, 126, 133, 141, 150, 160, 171, 185,
										 200, 220, 243, 273, 320, 380, 440, 510, 600, 700,
										 830, 1000, 1500, 3000};

	TH1D* h_diag = new TH1D("h_diag", "", nMassBin, MassBinEdges);

	Int_t nBin = h2D->GetNbinsX();
	for(Int_t i=0; i<nBin; i++)
	{
		Int_t i_bin = i+1;
		Double_t diag = h2D->GetBinContent(i_bin, i_bin);
		Double_t err_diag = h2D->GetBinError(i_bin, i_bin);
		// printf("[%2d bin] (diag, err_diag) = (%lf, %lf)\n", i_bin, diag, err_diag);
		h_diag->SetBinContent(i_bin, diag);
		h_diag->SetBinError(i_bin, err_diag);
	}

	return h_diag;
}