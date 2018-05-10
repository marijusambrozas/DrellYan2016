#include <Include/PlotTools.h>

class HistUnit
{
public:
	Int_t binSize;
	TString fileName;
	TString sampleTag;

	TH1D* h_mass_preFSR;
	TH1D* h_mass_dressed;
	TH1D* h_mass_postFSR;
	TH1D* h_mass_reco;

	TH1D* h_RelUnc_Stat;

	TH1D* h_diag_preFSR_reco;
	TH1D* h_diag_dressed_reco;
	TH1D* h_diag_postFSR_reco;

	TH2D* h2D_preFSR_reco;
	TH2D* h2D_dressed_reco;
	TH2D* h2D_postFSR_reco;

	HistUnit()
	{

	}

	HistUnit( Int_t _binSize, TString _sampleTag, TString _fileName )
	{
		this->binSize = _binSize;
		this->sampleTag = _sampleTag;
		this->fileName = _fileName;
		this->Init();
	}

	void Save()
	{

	}

	void Print2DHist_All()
	{
		this->Print2DHist( this->h2D_preFSR_reco );
		this->Print2DHist( this->h2D_dressed_reco );
		this->Print2DHist( this->h2D_postFSR_reco );
	}

private:
	void Print2DHist( TH2D* h )
	{
		gStyle->SetPalette(1);

		SampleInfo *sampleInfo = new SampleInfo(0, this->sampleTag, this->sampleTag);

		TString histName = h->GetName();
		TString XTitle = "";
		if( histName.Contains("postFSR") ) XTitle = "m (post-FSR) [GeV]";
		if( histName.Contains("dressed") ) XTitle = "m (dressed) [GeV]";
		if( histName.Contains("preFSR") ) XTitle = "m (pre-FSR) [GeV]";
		HistInfo* histInfo = new HistInfo( h->GetName(), XTitle, "m (reco) [GeV]");
		TH2Ext* hext = new TH2Ext( sampleInfo, histInfo, h );

		TCanvas *c;
		TString CanvasName = histInfo->name;
		CanvasName.ReplaceAll( "h2D_", "Local/c2D_"+sampleInfo->name+"_" );

		SetCanvas_Square2D( c, CanvasName, 0, 0 );

		hext->DrawAndSet( "COLZ" );
		hext->h->SetMinimum(0.001);
		hext->h->SetMaximum(1.0);
		hext->h->GetXaxis()->SetNoExponent();
		hext->h->GetXaxis()->SetMoreLogLabels();
		hext->h->GetYaxis()->SetNoExponent();
		hext->h->GetYaxis()->SetMoreLogLabels();

		TLatex latex;
		Latex_Simulation( latex );

		c->SaveAs(".png");

		hext->h->GetXaxis()->SetRangeUser(10, 310);
		hext->h->GetYaxis()->SetRangeUser(10, 310);
		c->SaveAs( TString::Format("%s_Zoomin.png", CanvasName.Data()) );
	}

	void Init()
	{
		TString	BinType = TString::Format("%dGeV", this->binSize);
		this->h_mass_preFSR = Get_Hist(this->fileName, sampleTag+"/h_mass_preFSR_"+BinType);
		this->h_mass_dressed = Get_Hist(this->fileName, sampleTag+"/h_mass_dressed_"+BinType);
		this->h_mass_postFSR = Get_Hist(this->fileName, sampleTag+"/h_mass_postFSR_"+BinType);
		this->h_mass_reco = Get_Hist(this->fileName, sampleTag+"/h_mass_reco_"+BinType);

		this->h2D_preFSR_reco = Get_Hist_2D(this->fileName, sampleTag+"/h2D_preFSR_reco_"+BinType);
		this->h2D_dressed_reco = Get_Hist_2D(this->fileName, sampleTag+"/h2D_dressed_reco_"+BinType);
		this->h2D_postFSR_reco = Get_Hist_2D(this->fileName, sampleTag+"/h2D_postFSR_reco_"+BinType);

		this->h2D_preFSR_reco = this->CalculateFractionPerBin( this->h2D_preFSR_reco, this->h_mass_preFSR );
		this->h2D_dressed_reco = this->CalculateFractionPerBin( this->h2D_dressed_reco, this->h_mass_dressed );
		this->h2D_postFSR_reco = this->CalculateFractionPerBin( this->h2D_postFSR_reco, this->h_mass_postFSR );

		this->h_RelUnc_Stat = this->MakeStatUncHist();

		this->h_diag_preFSR_reco = this->MakeDiagonalHist( this->h2D_preFSR_reco );
		this->h_diag_dressed_reco = this->MakeDiagonalHist( this->h2D_dressed_reco );
		this->h_diag_postFSR_reco = this->MakeDiagonalHist( this->h2D_postFSR_reco );
	}

	TH1D* MakeStatUncHist()
	{
		TH1D* h_Unc = (TH1D*)this->h_mass_reco->Clone();
		Int_t nBin = this->h_mass_reco->GetNbinsX();
		for(Int_t i=0; i<nBin; i++)
		{
			Int_t i_bin = i+1;
			Double_t nEvent = this->h_mass_reco->GetBinContent(i_bin);
			if( nEvent < 0 ) nEvent = 0;

			Double_t RelUnc_Stat = 0;
			if( nEvent > 0 ) RelUnc_Stat = 1.0 / sqrt(nEvent);

			h_Unc->SetBinContent(i_bin, RelUnc_Stat);
			h_Unc->SetBinError(i_bin, 0);
		}

		return h_Unc;
	}

	TH1D* MakeDiagonalHist( TH2D* h2D)
	{
		TH1D* h_diag = (TH1D*)this->h_mass_reco->Clone();

		Int_t nBin = h2D->GetNbinsX();
		for(Int_t i=0; i<nBin; i++)
		{
			Int_t i_bin = i+1;
			Double_t diag = h2D->GetBinContent(i_bin, i_bin);
			Double_t err_diag = h2D->GetBinError(i_bin, i_bin);
			h_diag->SetBinContent(i_bin, diag);
			h_diag->SetBinError(i_bin, err_diag);
		}

		return h_diag;
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
				Double_t nEvent = h_nEvents->GetBinContent(i_genbin, i_recobin);

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
};

class DrawingTool
{
public:
	Int_t binSize;
	TString fileName_MuMu;
	TString fileName_EE;

	HistUnit *histUnit_MuMu;
	HistUnit *histUnit_EE;

	DrawingTool( Int_t _binSize, TString _fileName_MuMu, TString _fileName_EE )
	{
		this->binSize = _binSize;
		this->fileName_MuMu = _fileName_MuMu;
		this->fileName_EE = _fileName_EE;

		this->Init();
	}

	void Draw()
	{
		SampleInfo* sampleInfo_MuMu = new SampleInfo( 0, "DYMuMu", "Z/#gamma* #rightarrow #mu#mu (post-FSR vs. reco)");
		sampleInfo_MuMu->SetColor( kBlue );

		SampleInfo* sampleInfo_EE2 = new SampleInfo( 0, "DYEE", "Z/#gamma* #rightarrow ee (post-FSR vs. reco)");
		sampleInfo_EE2->SetColor( kViolet+2 );

		SampleInfo* sampleInfo_EE = new SampleInfo( 0, "DYEE", "Z/#gamma* #rightarrow ee (dressed vs. reco)");
		sampleInfo_EE->SetColor( kGreen+2 );

		HistInfo* histInfo = new HistInfo( "", "m [GeV]", "diagonal terms");
		histInfo->SetYRange(0.7, 1.1);

		TH1D* h_MuMu = (TH1D*)this->histUnit_MuMu->h_diag_postFSR_reco->Clone();
		TH1D* h_EE2 = (TH1D*)this->histUnit_EE->h_diag_postFSR_reco->Clone();
		TH1D* h_EE = (TH1D*)this->histUnit_EE->h_diag_dressed_reco->Clone();
		Double_t maxX_MuMu = this->CalcMaxRangeX( h_MuMu );
		Double_t maxX_EE = this->CalcMaxRangeX( h_EE2 );
		Double_t maxX = maxX_MuMu > maxX_EE ? maxX_MuMu : maxX_EE;
		histInfo->SetXRange( 15, maxX );

		TH1Ext* hext_MuMu = new TH1Ext( sampleInfo_MuMu, histInfo, h_MuMu );
		TH1Ext* hext_EE2 = new TH1Ext( sampleInfo_EE2, histInfo, h_EE2 );
		TH1Ext* hext_EE = new TH1Ext( sampleInfo_EE, histInfo, h_EE );

		TCanvas *c;
		TString canvasName = "Local/c_diag_MuMu_EE_"+TString::Format("%02dGeV", this->binSize);
		SetCanvas_Square( c, canvasName );
		c->cd();

		hext_MuMu->DrawAndSet("HISTLPSAME");
		hext_EE2->DrawAndSet("HISTLPSAME");
		hext_EE->DrawAndSet("HISTLPSAME");

		TLegend *legend;
		SetLegend( legend, 0.50, 0.80, 0.97, 0.95 );
		hext_MuMu->AddToLegend( legend );
		hext_EE2->AddToLegend( legend );
		hext_EE->AddToLegend( legend );
		legend->Draw();

		TLatex latex;
		Latex_Simulation( latex );

		TString binSizeInfo = TString::Format("Bin size = %d GeV", this->binSize);
		latex.DrawLatexNDC(0.16, 0.91, "#scale[0.7]{#font[42]{"+binSizeInfo+"}}");

		TF1 *f_line;
		f_line = new TF1("f_line", "0.8", -10000, 10000);
		f_line->SetLineColor(kRed);
		f_line->SetLineWidth(1);
		f_line->Draw("PSAME");

		c->SaveAs(".pdf");
	}

	Double_t CalcMaxRangeX( TH1D* h )
	{
		Double_t maxX = 0;
		Int_t nBin = h->GetNbinsX();
		for(Int_t i=0; i<nBin; i++)
		{
			Int_t i_bin = i+1;
			Double_t value = h->GetBinContent(i_bin);
			if( value < 0.7 )
			{
				maxX = h->GetBinLowEdge( i_bin+1 ); // -- upper edge -- //
				printf("value = %lf -> maxX = %lf\n", value, maxX );
				break;
			}
		}

		return maxX;
	}

private:
	void Init()
	{
		this->histUnit_MuMu = new HistUnit( this->binSize, "DYMuMu", this->fileName_MuMu );
		this->histUnit_EE = new HistUnit( this->binSize, "DYEE", this->fileName_EE );
	}
};


void DrawPlot_DiagTermsForEachBinSize()
{
	TString fileName_MuMu = "./Local/ROOTFile_MakeHist_VariousBinning_MuMu_BinSize_21to40.root";
	TString fileName_EE = "./Local/ROOTFile_MakeHist_VariousBinning_EE_BinSize_21to40.root";

	for(Int_t i_binSize=21; i_binSize<41; i_binSize++)
	{
		DrawingTool *tool = new DrawingTool( i_binSize, fileName_MuMu, fileName_EE );
		tool->Draw();
	}
}