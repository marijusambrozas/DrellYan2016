// -- for nBin, MassBinEdges -- //
#include "MakeHist_MassResolution.h"
#include <Include/PlotTools_Old.h>

// -- sigmaL and R have same value! -- //
Double_t Cruijff2( Double_t *x, Double_t* par )
{
	// -- par[0] = A, par[1] = m, par[2] = sigma, par[3] = alphaL, par[4] = alphaR -- //
	Double_t sigma = par[2];

	Double_t alpha = 0.0;
	Double_t dx = (x[0] - par[1]);
	if(dx<0){
	  // sigma = par[2];
	  alpha = par[3];
	} else {
	  // sigma = par[3];
	  alpha = par[4];
	}

	Double_t f = 2*sigma*sigma + alpha*dx*dx ;
	return par[0]*exp(-dx*dx/f);
}

class FittingTool
{
public:
	TString ChannelType;
	TString FileName;
	vector<TH1D*> vec_Hist;

	Double_t MassBinEdges_Double[nBin+1];
	TH1D* h_NormChi2;
	TH1D* h_RMS;
	TH1D* h_Mean;
	TH1D* h_Resolution;

	FittingTool()
	{
		for(Int_t i=0; i<nBin+1; i++)
			MassBinEdges_Double[i] = (Double_t)MassBinEdges[i];
	}

	void Set_ChannelType( TString value )
	{
		this->ChannelType = value;
		this->Set_Histogram();	
	}

	void GetResolution()
	{
		UInt_t nHist = this->vec_Hist.size();
		for(UInt_t i_hist=0; i_hist<nHist; i_hist++)
		{
			this->Fit_GetResolution( i_hist, this->vec_Hist[i_hist] );
		}
	}

	void Save( TFile *f_output )
	{
		f_output->cd();
		this->h_NormChi2->Write();
		this->h_RMS->Write();
		this->h_Mean->Write();
		this->h_Resolution->Write();
	}

private:
	void Set_Histogram()
	{
		// -- make histograms -- //
		this->h_NormChi2 = new TH1D("h_NormChi2_"+this->ChannelType, "", nBin, MassBinEdges_Double);
		this->h_RMS = new TH1D("h_RMS_"+this->ChannelType, "", nBin, MassBinEdges_Double);
		this->h_Mean = new TH1D("h_Mean_"+this->ChannelType, "", nBin, MassBinEdges_Double);
		this->h_Resolution = new TH1D("h_Resolution_"+this->ChannelType, "", nBin, MassBinEdges_Double);

		// -- get histograms -- //
		this->FileName = TString::Format("Local/ROOTFile_MakeHist_MassResolution_%s.root", this->ChannelType.Data() );

		for(Int_t i=0; i<nBin; i++)
		{
			TString TStr_MassRange = TString::Format("M%dto%d", MassBinEdges[i], MassBinEdges[i+1]);
			TString HistName = TString::Format("DY%s/h_RelDiff_%s", this->ChannelType.Data(), TStr_MassRange.Data());

			TH1D* h_temp = Get_Hist( this->FileName, HistName );
			// h_temp->Rebin(10);
			this->vec_Hist.push_back( h_temp );
		}
	}

	void Fit_GetResolution( Int_t i_hist, TH1D* h )
	{
		TH1D* h_Clone = (TH1D*)h->Clone();

		TString fitFuncName = "Cruijff2";
		Int_t nParam = 5;

		Double_t rangeMax = 2.0*h->GetRMS();

		TF1* fitFunc_Gauss = new TF1("gaus", "gaus", (-1)*rangeMax, rangeMax);
		fitFunc_Gauss->SetParameters(h->Integral(), h->GetMean(), h->GetRMS());
		h->Fit( fitFunc_Gauss, "R", "", (-1)*rangeMax, rangeMax );

		Double_t constInit = fitFunc_Gauss->GetParameter(0);
		Double_t meanInit = fitFunc_Gauss->GetParameter(1);
		Double_t sigmaInit = fitFunc_Gauss->GetParameter(2);

		TF1 *fitFunc = NULL;
		if( fitFuncName == "Cruijff2")
		{
			fitFunc = new TF1(fitFuncName, Cruijff2, (-1)*rangeMax, rangeMax, nParam);
			fitFunc->SetParameter(0, constInit); fitFunc->SetParName(0, "const");
			fitFunc->SetParameter(1, meanInit); fitFunc->SetParName(1, "mean");
			fitFunc->SetParameter(2, sigmaInit); fitFunc->SetParName(2, "sigma");
			fitFunc->SetParameter(3, 0.1); fitFunc->SetParName(3, "alphaL");
			fitFunc->SetParameter(4, 0.1); fitFunc->SetParName(4, "alphaR");
		}

		// -- fit -- //
		h->Fit( fitFuncName, "R", "", (-1)*rangeMax, rangeMax );

		Double_t chi2 = fitFunc->GetChisquare();
		Double_t nDOF = fitFunc->GetNDF();
		Double_t NormChi2 = chi2/nDOF;

		Double_t RMS = h->GetRMS();
		Double_t RMS_err = h->GetRMSError();

		Double_t mean = fitFunc->GetParameter(1);
		Double_t mean_err = fitFunc->GetParError(1);

		Double_t sigma = fitFunc->GetParameter(2);
		Double_t sigma_err = fitFunc->GetParError(2);

		Int_t i_bin = i_hist+1;
		this->h_NormChi2->SetBinContent(i_bin, NormChi2);
		this->h_NormChi2->SetBinError(i_bin, 0);

		this->h_RMS->SetBinContent(i_bin, RMS);
		this->h_RMS->SetBinError(i_bin, RMS_err);

		this->h_Mean->SetBinContent(i_bin, mean);
		this->h_Mean->SetBinError(i_bin, mean_err);

		this->h_Resolution->SetBinContent(i_bin, sigma);
		this->h_Resolution->SetBinError(i_bin, sigma_err);

		TCanvas *c;
		TString CanvasName = h->GetName();
		CanvasName.ReplaceAll("h_", "Local/c_"+this->ChannelType+"_");

		HistInfo* Hist = new HistInfo( kBlack, "temp", h_Clone );
		fitFunc_Gauss->SetLineColor(kRed);
		fitFunc_Gauss->SetMarkerColorAlpha(kWhite, 0);
		fitFunc_Gauss->SetFillColorAlpha(kWhite, 0);
		fitFunc->SetLineColor(kBlue);
		fitFunc->SetMarkerColorAlpha(kWhite, 0);
		fitFunc->SetFillColorAlpha(kWhite, 0);

		SetCanvas_Square( c, CanvasName );
		c->cd();
		Hist->Draw("EPSAME");
		fitFunc_Gauss->Draw("SAME");
		fitFunc->Draw("SAME");

		TLegend *legend;
		SetLegend(legend, 0.80, 0.75, 0.95, 0.95 );
		legend->AddEntry( fitFunc_Gauss, "Gaussian" );
		legend->AddEntry( fitFunc, fitFuncName );
		legend->Draw();

		SetHistFormat_SinglePad(Hist->h, "(M^{Reco}-M^{Gen})/M^{Gen}", "Entries");
		Hist->h->GetXaxis()->SetRangeUser( (-1)*rangeMax, rangeMax );

		TLatex latex;
		Latex_Simulation( latex );
		TString MassRange = TString::Format("%d < M^{Gen} < %d GeV", MassBinEdges[i_hist], MassBinEdges[i_hist+1]);
		latex.DrawLatexNDC(0.16, 0.91, "#font[42]{#scale[0.7]{"+MassRange+"}}");
		TString normChi2Info = TString::Format("#chi2/ndof = %.2lf/%.0lf = %.2lf", chi2, nDOF, NormChi2);
		latex.DrawLatexNDC(0.16, 0.88, "#font[42]{#scale[0.5]{"+normChi2Info+"}}");

		if( fitFuncName == "gaus" )
		{
			TString meanInfo = TString::Format("mean = %.5lf", mean);
			latex.DrawLatexNDC(0.16, 0.84, "#font[42]{#scale[0.5]{"+meanInfo+"}}");
			TString sigmaInfo = TString::Format("resolution(#sigma) = %.3lf", sigma);
			latex.DrawLatexNDC(0.16, 0.80, "#font[42]{#scale[0.5]{"+sigmaInfo+"}}");
		}

		if( fitFuncName == "Cruijff2" )
		{
			TString meanInfo = TString::Format("mean = %.5lf", mean);
			latex.DrawLatexNDC(0.16, 0.84, "#font[42]{#scale[0.5]{"+meanInfo+"}}");
			TString sigmaInfo = TString::Format("resolution(#sigma) = %.3lf", sigma);
			latex.DrawLatexNDC(0.16, 0.80, "#font[42]{#scale[0.5]{"+sigmaInfo+"}}");

			Double_t alphaL = fitFunc->GetParameter(3);
			Double_t alphaR = fitFunc->GetParameter(4);

			TString alphaLInfo = TString::Format("#alpha_{L} = %.3lf", alphaL);
			latex.DrawLatexNDC(0.16, 0.76, "#font[42]{#scale[0.5]{"+alphaLInfo+"}}");
			TString alphaRInfo = TString::Format("#alpha_{R} = %.3lf", alphaR);
			latex.DrawLatexNDC(0.16, 0.72, "#font[42]{#scale[0.5]{"+alphaRInfo+"}}");
		}

		c->SaveAs(".pdf");
	}

};

void GetResolution()
{
	FittingTool *tool = new FittingTool();
	tool->Set_ChannelType("MuMu");
	tool->GetResolution();

	TFile *f_output = TFile::Open("ROOTFile_GetResolution.root", "RECREATE");
	tool->Save( f_output );
	f_output->Close();
}