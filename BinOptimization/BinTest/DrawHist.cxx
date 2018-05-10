#include <Include/PlotTools.h>

class HistContainer
{
public:
	TString subDirName;
	// vector< TH1D* > vec_Hist;
	TH1D* h_mass_preFSR;
	TH1D* h_mass_dressed;
	TH1D* h_mass_postFSR;
	TH1D* h_mass_reco;

	// vector< TH2D* > vec_2DHist;
	TH2D* h2D_preFSR_reco;
	TH2D* h2D_dressed_reco;
	TH2D* h2D_postFSR_reco;

	TH1D* h_diag_preFSR_reco;
	TH1D* h_diag_dressed_reco;
	TH1D* h_diag_postFSR_reco;

	TH1D* h_RelUnc_Stat;

	HistContainer( TString fileName )
	{
		this->subDirName = "";
		if( fileName.Contains("MuMu") )
			this->subDirName = "DYMuMu";
		else if( fileName.Contains("EE") )
			this->subDirName = "DYEE";

		this->h_mass_preFSR = Get_Hist( fileName, subDirName+"/h_mass_preFSR" );
		this->h_mass_dressed = Get_Hist( fileName, subDirName+"/h_mass_dressed" );
		this->h_mass_postFSR = Get_Hist( fileName, subDirName+"/h_mass_postFSR" );
		this->h_mass_reco = Get_Hist( fileName, subDirName+"/h_mass_reco" );

		this->h2D_preFSR_reco = Get_Hist_2D( fileName, subDirName+"/h2D_preFSR_reco" );
		this->h2D_dressed_reco = Get_Hist_2D( fileName, subDirName+"/h2D_dressed_reco" );
		this->h2D_postFSR_reco = Get_Hist_2D( fileName, subDirName+"/h2D_postFSR_reco" );

		this->h2D_preFSR_reco = this->CalculateFractionPerBin( this->h2D_preFSR_reco, this->h_mass_preFSR );
		this->h2D_dressed_reco = this->CalculateFractionPerBin( this->h2D_dressed_reco, this->h_mass_dressed );
		this->h2D_postFSR_reco = this->CalculateFractionPerBin( this->h2D_postFSR_reco, this->h_mass_postFSR );

		this->h_diag_preFSR_reco = MakeDiagonalHist( h2D_preFSR_reco );
		this->h_diag_dressed_reco = MakeDiagonalHist( h2D_dressed_reco );
		this->h_diag_postFSR_reco = MakeDiagonalHist( h2D_postFSR_reco );

		this->h_RelUnc_Stat = this->MakeStatUncHist();
	}

	void PrintHist_All()
	{
		this->Print2DHist( this->h2D_preFSR_reco );
		this->Print2DHist( this->h2D_dressed_reco );
		this->Print2DHist( this->h2D_postFSR_reco );

		this->Draw_RelUnc_Stat();
		this->Draw_DiagonalTerms( this->h_diag_preFSR_reco, "preFSR" );
		this->Draw_DiagonalTerms( this->h_diag_dressed_reco, "dressed" );
		this->Draw_DiagonalTerms( this->h_diag_postFSR_reco, "postFSR" );
	}

	TH1D* MassHist( TString type )
	{
		if( type == "reco" ) return this->h_mass_reco;
		if( type == "postFSR" ) return this->h_mass_postFSR;
		if( type == "dressed" ) return this->h_mass_dressed;
		if( type == "preFSR" ) return this->h_mass_preFSR;

		return NULL;
	}

private:
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

	void Print2DHist( TH2D* h )
	{
		gStyle->SetPalette(1);

		TString HistName = h->GetName();

		TCanvas *c;
		TString CanvasName = HistName;
		CanvasName.ReplaceAll( "h2D_", "Local/c2D_"+this->subDirName+"_" );

		SetCanvas_Square( c, CanvasName, 1, 1 );
		c->cd();
		c->SetRightMargin( 0.12 );
		c->SetLeftMargin( 0.15 );
		h->Draw("COLZ");

		h->SetStats(kFALSE);

		TString XTitle = "";
		if( HistName.Contains("postFSR") ) XTitle = "m (post-FSR) [GeV]";
		if( HistName.Contains("dressed") ) XTitle = "m (dressed) [GeV]";
		if( HistName.Contains("preFSR") ) XTitle = "m (pre-FSR) [GeV]";
		SetAxis_SinglePad( h->GetXaxis(), h->GetYaxis(), XTitle, "m (reco) [GeV]");
		h->GetYaxis()->SetLabelSize( 0.04 );
		h->GetYaxis()->SetTitleOffset( 1.45 );
		h->GetZaxis()->SetRangeUser(0.001, 1.0);

		c->SaveAs(".png");

		// h->GetXaxis()->SetRangeUser(10, 310);
		// h->GetYaxis()->SetRangeUser(10, 310);
		// c->SaveAs( TString::Format("%s_Zoomin.png", CanvasName.Data()) );
	}

	void Draw_RelUnc_Stat()
	{
		SampleInfo* sampleInfo = new SampleInfo(); // -- dummy -- //
		sampleInfo->SetColor( kBlack );
		HistInfo *histInfo = new HistInfo( "", "m [GeV]", "Stat. unc. (%)");

		TH1Ext *hext = new TH1Ext( sampleInfo, histInfo, this->h_RelUnc_Stat );
		hext->h->Scale( 100 );

		TCanvas *c;
		TString CanvasName = "Local/c_RelUnc_Stat_"+this->subDirName;
		SetCanvas_Square( c, CanvasName, 1, 1 );
		c->cd();

		hext->DrawAndSet( "HISTLPSAME" );

		TLatex latex;
		Latex_Simulation( latex );
		TString channelInfo = "";
		if( this->subDirName == "DYEE" ) channelInfo = "e channel";
		if( this->subDirName == "DYMuMu" ) channelInfo = "#mu channel";
		latex.DrawLatexNDC( 0.16, 0.91, "#font[42]{#scale[0.7]{"+channelInfo+"}}");

		c->SaveAs(".pdf");
	}
	
	void Draw_DiagonalTerms( TH1D* h_diag, TString type )
	{
		SampleInfo* sampleInfo = new SampleInfo(); // -- dummy -- //
		sampleInfo->SetColor( kBlack );

		HistInfo *histInfo = new HistInfo( "", "m [GeV]", "Diagonal term");
		histInfo->SetYRange(0, 1.1);

		TH1Ext *hext = new TH1Ext( sampleInfo, histInfo, h_diag );
		// hext->h->Scale( 100 );

		TCanvas *c;
		TString CanvasName = "Local/c_Diag_"+this->subDirName+"_"+type;;
		SetCanvas_Square( c, CanvasName, 1, 0 );
		c->cd();

		hext->DrawAndSet( "HISTLPSAME" );

		TLatex latex;
		Latex_Simulation( latex );
		TString channelInfo = "";
		if( this->subDirName == "DYEE" ) channelInfo = "e channel";
		if( this->subDirName == "DYMuMu" ) channelInfo = "#mu channel";
		latex.DrawLatexNDC( 0.16, 0.91, "#font[42]{#scale[0.7]{"+channelInfo+"}}");

		c->SaveAs(".pdf");
	}
};

class DrawingTool
{
public:
	HistContainer* hist_MuMu;
	HistContainer* hist_EE;

	DrawingTool()
	{
		this->hist_MuMu = new HistContainer( "./Local/ROOTFile_MakeHist_FixBinnning_MuMu.root" );
		hist_MuMu->PrintHist_All();

		this->hist_EE = new HistContainer( "./Local/ROOTFile_MakeHist_FixBinnning_EE.root" );
		hist_EE->PrintHist_All();
	}

	void DrawDiagonalTerms()
	{
		SampleInfo *sampleInfo_MuMu = new SampleInfo( 0, "DYMuMu_postFSR", "Z/#gamma* #rightarrow #mu#mu (postFSR vs. reco)");
		sampleInfo_MuMu->SetColor( kBlue );

		SampleInfo *sampleInfo_EE = new SampleInfo( 0, "DYEE_postFSR", "Z/#gamma* #rightarrow ee (postFSR vs. reco)");
		sampleInfo_EE->SetColor( kViolet );

		SampleInfo *sampleInfo_EEDressed = new SampleInfo( 0, "DYEE_dressed", "Z/#gamma* #rightarrow ee (dressed vs. reco)");
		sampleInfo_EEDressed->SetColor( kGreen+2 );

		HistInfo *histInfo = new HistInfo("", "m [GeV]", "Diagonal term");
		histInfo->SetYRange(0.4, 1.1);

		TH1Ext *hext_MuMu = new TH1Ext( sampleInfo_MuMu, histInfo, this->hist_MuMu->h_diag_postFSR_reco );
		TH1Ext *hext_EE = new TH1Ext( sampleInfo_EE, histInfo, this->hist_EE->h_diag_postFSR_reco );
		TH1Ext *hext_EEDressed = new TH1Ext( sampleInfo_EEDressed, histInfo, this->hist_EE->h_diag_dressed_reco );

		TCanvas *c;
		TString CanvasName = "Local/c_Diag_Comparison";
		SetCanvas_Square( c, CanvasName, 1, 0 );
		c->cd();

		hext_MuMu->DrawAndSet( "HISTLPSAME" );
		hext_EE->DrawAndSet( "HISTLPSAME" );
		hext_EEDressed->DrawAndSet( "HISTLPSAME" );

		TLegend *legend;
		SetLegend( legend, 0.50, 0.75, 0.97, 0.95 );
		hext_MuMu->AddToLegend( legend );
		hext_EE->AddToLegend( legend );
		hext_EEDressed->AddToLegend( legend );
		legend->Draw();

		TF1 *f_line;
		f_line = new TF1("f_line", "0.8", -10000, 10000);
		f_line->SetLineColor(kRed);
		f_line->SetLineWidth(1);
		f_line->Draw("PSAME");

		TLatex latex;
		Latex_Simulation( latex );

		c->SaveAs(".pdf");
	}

	void DrawMassDistribution( TString type )
	{
		SampleInfo *sampleInfo_MuMu = new SampleInfo( 0, "DYMuMu", "Z/#gamma* #rightarrow #mu#mu");
		sampleInfo_MuMu->SetColor( kBlue );

		SampleInfo *sampleInfo_EE = new SampleInfo( 0, "DYEE", "Z/#gamma* #rightarrow ee");
		sampleInfo_EE->SetColor( kGreen+2 );

		HistInfo *histInfo = new HistInfo( "", "m [GeV]", "Entries per bin");
		
		TH1Ext *hext_MuMu = new TH1Ext( sampleInfo_MuMu, histInfo, this->hist_MuMu->MassHist(type) ); 
		TH1Ext *hext_EE = new TH1Ext( sampleInfo_EE, histInfo, this->hist_EE->MassHist(type) );

		TCanvas *c;
		TString CanvasName = "Local/c_MassHist_Comparison_"+type;
		SetCanvas_Square( c, CanvasName, 1, 1 );
		c->cd();

		hext_MuMu->DrawAndSet( "HISTLPSAME" );
		hext_EE->DrawAndSet( "HISTLPSAME" );

		TLegend *legend;
		SetLegend( legend, 0.70, 0.80, 0.97, 0.95 );
		hext_MuMu->AddToLegend( legend );
		hext_EE->AddToLegend( legend );
		legend->Draw();

		TLatex latex;
		Latex_Simulation( latex );
		TString typeInfo = "";
		if( type == "reco" ) typeInfo = "Reconstruction level";
		if( type == "postFSR" ) typeInfo = "post-FSR level";
		if( type == "dressed" ) typeInfo = "dressed level";
		if( type == "preFSR" ) typeInfo = "pre-FSR level";
		latex.DrawLatexNDC( 0.16, 0.91, "#scale[0.7]{#font[42]{"+typeInfo+"}}");

		c->SaveAs(".pdf");
	}

	void DrawStatUnc()
	{
		SampleInfo *sampleInfo_MuMu = new SampleInfo( 0, "DYMuMu", "Z/#gamma* #rightarrow #mu#mu");
		sampleInfo_MuMu->SetColor( kBlue );

		SampleInfo *sampleInfo_EE = new SampleInfo( 0, "DYEE", "Z/#gamma* #rightarrow ee");
		sampleInfo_EE->SetColor( kGreen+2 );

		HistInfo *histInfo = new HistInfo( "", "m [GeV]", "Stat. unc. (%)");
		histInfo->SetYRange( 1e-2, 10 );
		
		TH1Ext *hext_MuMu = new TH1Ext( sampleInfo_MuMu, histInfo, this->hist_MuMu->h_RelUnc_Stat ); hext_MuMu->h->Scale( 100 );
		TH1Ext *hext_EE = new TH1Ext( sampleInfo_EE, histInfo, this->hist_EE->h_RelUnc_Stat ); hext_EE->h->Scale( 100 );

		TCanvas *c;
		TString CanvasName = "Local/c_RelUncStat_Comparison";
		SetCanvas_Square( c, CanvasName, 1, 1 );
		c->cd();

		hext_MuMu->DrawAndSet( "HISTLPSAME" );
		hext_EE->DrawAndSet( "HISTLPSAME" );

		TLegend *legend;
		SetLegend( legend, 0.70, 0.80, 0.97, 0.95 );
		hext_MuMu->AddToLegend( legend );
		hext_EE->AddToLegend( legend );
		legend->Draw();

		TLatex latex;
		Latex_Simulation( latex );

		c->SaveAs(".pdf");
	}

};


void DrawHist()
{
	DrawingTool *tool = new DrawingTool();
	tool->DrawDiagonalTerms();
	tool->DrawMassDistribution( "reco" );
	tool->DrawMassDistribution( "postFSR" );
	tool->DrawMassDistribution( "dressed" );
	tool->DrawMassDistribution( "preFSR" );
	tool->DrawStatUnc();
}