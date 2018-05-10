#include <Include/PlotTools_Old.h>

class PrintingTool
{
public:
	vector< TH1D* > vec_Hist_EE;
	vector< TH1D* > vec_Hist_MuMu;
	TString fileName;

	PrintingTool(TString _fileName = "ROOTFile_GetResolution.root")
	{
		this->fileName = _fileName;

		this->Init( "EE", this->vec_Hist_EE );
		this->Init( "MuMu", this->vec_Hist_MuMu );
	}

	void PrintAll()
	{
		Int_t nHist = (Int_t)this->vec_Hist_EE.size();
		for( Int_t i_hist=0; i_hist<nHist; i_hist++)
		{
			this->PrintHistogram_EachChannel( this->vec_Hist_EE[i_hist] );
			this->PrintHistogram_EachChannel( this->vec_Hist_MuMu[i_hist] );
			this->PrintHistogram_Superimpose( this->vec_Hist_EE[i_hist], this->vec_Hist_MuMu[i_hist] );
		}
	}

private:
	void Init( TString channelType, vector<TH1D*>& vec_Hist )
	{
		TH1D* h_Mean = Get_Hist(this->fileName, "h_Mean_"+channelType);
		vec_Hist.push_back( h_Mean );

		TH1D* h_Resolution = Get_Hist(this->fileName, "h_Resolution_"+channelType);
		vec_Hist.push_back( h_Resolution );

		TH1D* h_RMS = Get_Hist(this->fileName, "h_RMS_"+channelType);
		vec_Hist.push_back( h_RMS );

		TH1D* h_NormChi2 = Get_Hist(this->fileName, "h_NormChi2_"+channelType);
		vec_Hist.push_back( h_NormChi2 );
	}

	void PrintHistogram_EachChannel( TH1D* h )
	{
		TString HistName = h->GetName();
		TString channelType = "";
		if( HistName.Contains("EE") ) channelType = "EE";
		if( HistName.Contains("MuMu") ) channelType = "MuMu";

		TCanvas *c;
		TString CanvasName = HistName;
		
		CanvasName.ReplaceAll("h_", "Local/c_"+channelType+"_");

		HistInfo *Hist = new HistInfo( kBlack, "temp", h );

		SetCanvas_Square( c, CanvasName, 1, 0 );
		c->cd();
		Hist->Draw("HISTLPSAME");

		TString titleY = this->GetTitle_Y( HistName );
		SetHistFormat_SinglePad(Hist->h, "m [GeV]", titleY);
		if( HistName.Contains("Resolution") )
			Hist->h->GetYaxis()->SetRangeUser(0, 0.05);

		TLatex latex;
		Latex_Simulation( latex );
		TString channelInfo = "";
		if( channelType == "EE" ) channelInfo = "e channel";
		if( channelType == "MuMu" ) channelInfo = "#mu channel";
		latex.DrawLatexNDC(0.16, 0.91, "#font[42]{#scale[0.6]{"+channelInfo+"}}");

		c->SaveAs(".pdf");
	}

	void PrintHistogram_Superimpose( TH1D* h_EE, TH1D* h_MuMu )
	{
		TString HistName = h_EE->GetName();

		TCanvas *c;
		TString CanvasName = HistName;
		
		CanvasName.ReplaceAll("h_", "Local/c_");
		CanvasName.ReplaceAll("_EE", "");

		HistInfo *Hist_EE = new HistInfo( kGreen+2, "e channel", h_EE );
		HistInfo *Hist_MuMu = new HistInfo( kBlue, "#mu channel", h_MuMu );

		SetCanvas_Square( c, CanvasName, 1, 0 );
		c->cd();
		Hist_EE->Draw("HISTLPSAME");
		Hist_MuMu->Draw("HISTLPSAME");

		TLegend *legend;
		SetLegend( legend, 0.75, 0.80, 0.95, 0.95 );
		Hist_EE->AddToLegend( legend );
		Hist_MuMu->AddToLegend( legend );
		legend->Draw();

		TString titleY = this->GetTitle_Y( HistName );
		SetHistFormat_SinglePad(Hist_EE->h, "m [GeV]", titleY);
		if( HistName.Contains("Resolution") )
			Hist_EE->h->GetYaxis()->SetRangeUser(0, 0.08);

		TLatex latex;
		Latex_Simulation( latex );

		c->SaveAs(".pdf");
	}

	TString GetTitle_Y( TString HistName )
	{
		TString titleY = "";
		if( HistName.Contains("NormChi2") )
			titleY = "Normalized #chi^{2}";
		else if( HistName.Contains("RMS") )
			titleY = "RMS";
		else if( HistName.Contains("Mean") )
			titleY = "Mean";
		else if( HistName.Contains("Resolution") )
			titleY = "Resolution";

		return titleY;
	}

};

void PrintHistogram()
{
	PrintingTool* tool = new PrintingTool("ROOTFile_GetResolution.root");
	tool->PrintAll();
}