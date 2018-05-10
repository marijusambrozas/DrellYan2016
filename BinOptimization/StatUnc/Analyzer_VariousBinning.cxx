#include <Include/PlotTools_Old.h>

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

		TString HistName = h->GetName();

		TCanvas *c;
		TString CanvasName = HistName;
		CanvasName.ReplaceAll( "h2D_", "Local/c2D_"+this->sampleTag+"_" );

		SetCanvas_Square( c, CanvasName );
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
		h->GetZaxis()->SetRangeUser(0, 1.0);

		c->SaveAs(".png");

		h->GetXaxis()->SetRangeUser(10, 310);
		h->GetYaxis()->SetRangeUser(10, 310);
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

class Analyzer
{
public:
	TString channelType;

	vector<Int_t> vec_BinSize;
	vector<Int_t> vec_Color;
	vector<HistUnit*> vec_HistUnit;


	Analyzer()
	{

	}

	Analyzer( TString _channelType )
	{
		this->SetChannelType( _channelType );
		this->Init();
	}

	void SetChannelType( TString _channelType )
	{
		this->channelType = _channelType;
	}

	void Print2DHist()
	{
		for(const auto& histUnit : vec_HistUnit )
			histUnit->Print2DHist_All();
	}

	void Draw_RelUnc_Stat()
	{
		vector<HistInfo*> vec_HistInfo;

		Int_t nHist = (Int_t)vec_HistUnit.size();
		for(Int_t i=0; i<nHist; i++)
		{
			TString Legend = TString::Format("Bin size = %d GeV", this->vec_BinSize[i]);
			TH1D* h_temp = this->vec_HistUnit[i]->h_RelUnc_Stat;
			h_temp->Scale( 100 ); // -- convert to % -- //

			HistInfo *histInfo = new HistInfo( this->vec_Color[i], Legend, h_temp );
			vec_HistInfo.push_back( histInfo );
		}

		TCanvas *c;
		TString CanvasName = "Local/c_RelUnc_Stat_"+this->channelType;
		SetCanvas_Square( c, CanvasName, 1, 1 );
		c->cd();

		TLegend *legend;
		SetLegend( legend, 0.40, 0.80, 0.95, 0.95 );
		legend->SetNColumns(2);

		for(Int_t i=0; i<nHist; i++)
		{
			vec_HistInfo[i]->Draw("HISTLPSAME");
			vec_HistInfo[i]->AddToLegend( legend );
		}
		legend->Draw();
		SetHistFormat_SinglePad( vec_HistInfo[0]->h, "m [GeV]", "Stat. Unc. (%)" );

		TLatex latex;
		Latex_Simulation( latex );
		TString channelInfo = "";
		if( this->channelType == "EE" ) channelInfo = "e channel";
		if( this->channelType == "MuMu" ) channelInfo = "#mu channel";
		latex.DrawLatexNDC( 0.16, 0.91, "#font[42]{#scale[0.7]{"+channelInfo+"}}");
		vec_HistInfo[0]->h->GetYaxis()->SetRangeUser(0.01, 200);

		c->SaveAs(".pdf");

		for(Int_t i=0; i<nHist; i++)
		{
			c->SetLogx(kFALSE);
			c->SetLogy(kFALSE);
			vec_HistInfo[i]->h->GetXaxis()->SetRangeUser(15, 310);
			vec_HistInfo[i]->h->GetYaxis()->SetRangeUser(0, 10);
		}
		c->SaveAs(CanvasName+"_Zoomin.pdf");
	}

	void Draw_DiagonalHist( TString GenType = "postFSR" )
	{
		vector<HistInfo*> vec_HistInfo;

		Int_t nHist = (Int_t)vec_HistUnit.size();
		for(Int_t i=0; i<nHist; i++)
		{
			TString Legend = TString::Format("Bin size = %d GeV", this->vec_BinSize[i]);
			TH1D* h_temp;
			if( GenType == "postFSR" ) h_temp = this->vec_HistUnit[i]->h_diag_postFSR_reco;
			if( GenType == "dressed" ) h_temp = this->vec_HistUnit[i]->h_diag_dressed_reco;
			if( GenType == "preFSR" ) h_temp = this->vec_HistUnit[i]->h_diag_preFSR_reco;

			HistInfo *histInfo = new HistInfo( this->vec_Color[i], Legend, h_temp );
			vec_HistInfo.push_back( histInfo );
		}
		this->Print_ReferencePoint( 0.9, vec_HistInfo );
		this->Print_ReferencePoint( 0.8, vec_HistInfo );

		TCanvas *c;
		TString CanvasName = "Local/c_diag_"+GenType+"_reco_"+this->channelType;
		SetCanvas_Square( c, CanvasName );
		c->cd();

		TLegend *legend;
		SetLegend( legend, 0.40, 0.80, 0.95, 0.95 );
		legend->SetNColumns(2);

		for(Int_t i=0; i<nHist; i++)
		{
			vec_HistInfo[i]->Draw("HISTLPSAME");
			vec_HistInfo[i]->AddToLegend( legend );
		}
		legend->Draw();
		SetHistFormat_SinglePad( vec_HistInfo[0]->h, "m [GeV]", "Diagonal term" );

		TLatex latex;
		Latex_Simulation( latex );
		TString channelInfo = "";
		if( this->channelType == "EE" ) channelInfo = "e channel";
		if( this->channelType == "MuMu" ) channelInfo = "#mu channel";
		latex.DrawLatexNDC( 0.16, 0.91, "#font[42]{#scale[0.7]{"+channelInfo+"}}");
		vec_HistInfo[0]->h->GetYaxis()->SetRangeUser(0, 1.2);

		c->SaveAs(".pdf");

		for(Int_t i=0; i<nHist; i++)
		{
			vec_HistInfo[i]->h->GetXaxis()->SetRangeUser(15, 310);
		}
		c->SaveAs(CanvasName+"_Zoomin.pdf");
	}

private:
	void Init()
	{
		// this->vec_BinSize.push_back( 1 );
		// this->vec_BinSize.push_back( 2 );
		// this->vec_BinSize.push_back( 3 );
		// this->vec_BinSize.push_back( 4 );
		// this->vec_BinSize.push_back( 5 );

		// this->vec_BinSize.push_back( 6 );
		// this->vec_BinSize.push_back( 7 );
		// this->vec_BinSize.push_back( 8 );
		// this->vec_BinSize.push_back( 9 );
		// this->vec_BinSize.push_back( 10 );

		// this->vec_BinSize.push_back( 11 );
		// this->vec_BinSize.push_back( 12 );
		// this->vec_BinSize.push_back( 13 );
		// this->vec_BinSize.push_back( 14 );
		// this->vec_BinSize.push_back( 15 );

		this->vec_BinSize.push_back( 16 );
		this->vec_BinSize.push_back( 17 );
		this->vec_BinSize.push_back( 18 );
		this->vec_BinSize.push_back( 19 );
		this->vec_BinSize.push_back( 20 );

		this->vec_Color.push_back( kBlack );
		this->vec_Color.push_back( kRed );
		this->vec_Color.push_back( kBlue );
		this->vec_Color.push_back( kGreen+2 );
		this->vec_Color.push_back( kViolet );

		TString fileName = TString::Format("Local/ROOTFile_MakeHist_VariousBinning_%s_BinSize_1to20.root", this->channelType.Data());
		TString sampleTag = "DY"+this->channelType;
		for( const auto& binSize : this->vec_BinSize )
		{
			HistUnit *histUnit = new HistUnit( binSize, sampleTag, fileName );
			this->vec_HistUnit.push_back( histUnit );
		}
	}

	void Print_ReferencePoint( Double_t Ref, vector<HistInfo*>& vec_HistInfo )
	{
		TString HistName = vec_HistInfo[0]->h->GetName();
		cout << "[" << this->channelType << "] HistName: " << HistName << ", Ref value = " << Ref << endl;
		Int_t nHist = (Int_t)vec_HistInfo.size();
		for(Int_t i_hist=0; i_hist<nHist; i_hist++)
		{
			TH1D* h_temp = vec_HistInfo[i_hist]->h;
			Int_t nBin = (Int_t)h_temp->GetNbinsX();
			for(Int_t i=0; i<nBin; i++)
			{
				Int_t i_bin = i+1;
				Double_t LowerBinEdge = h_temp->GetBinLowEdge(i_bin);
				if( LowerBinEdge < 15 ) continue;

				Double_t value = h_temp->GetBinContent(i_bin);
				// cout << "value = " << value << endl;
				if( value < Ref )
				{

					printf("%.0lf GeV\t", LowerBinEdge);
					break;
				}
			} // -- iteration over bins -- //
		}
		cout << endl;
		cout << endl;
	}
};

void Analyzer_VariousBinning()
{
	Analyzer *analyzer_MuMu = new Analyzer("MuMu");
	analyzer_MuMu->Print2DHist();
	analyzer_MuMu->Draw_DiagonalHist();
	analyzer_MuMu->Draw_RelUnc_Stat();

	Analyzer *analyzer_EE = new Analyzer("EE");
	analyzer_EE->Print2DHist();
	analyzer_EE->Draw_DiagonalHist();
	analyzer_EE->Draw_DiagonalHist("dressed");
	analyzer_EE->Draw_RelUnc_Stat();
}