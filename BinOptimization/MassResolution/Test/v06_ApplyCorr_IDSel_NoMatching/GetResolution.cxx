// -- for nBin, MassBinEdges -- //
#include "MakeHist_MassResolution.h"
#include <Include/PlotTools_Old.h>


Double_t Cruijff( Double_t *x, Double_t* par )
{
	// -- par[0] = A, par[1] = m, par[2] = sigmaL, par[3] = sigmaR, par[4] = alphaL, par[5] = alphaR -- //
	Double_t sigma = 0.0;
	Double_t alpha = 0.0;
	Double_t dx = (x[0] - par[1]);
	if(dx<0){
	  sigma = par[2];
	  alpha = par[4];
	} else {
	  sigma = par[3];
	  alpha = par[5];
	}

	Double_t f = 2*sigma*sigma + alpha*dx*dx ;
	return par[0]*exp(-dx*dx/f);
}

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

Double_t CruijffGauss( Double_t *x, Double_t* par )
{
	// -- par[0] = A, par[1] = m, par[2] = sigmaL, par[3] = sigmaR, par[4] = alphaL, par[5] = alphaR, par[6] = mRes, par[7] = sigmaRes -- //
	Double_t sigma = 0.0;
	Double_t alpha = 0.0;
	Double_t dx = (x[0] - par[1]);
	if(dx<0){
	  sigma = par[2];
	  alpha = par[4];
	} else {
	  sigma = par[3];
	  alpha = par[5];
	}

	Double_t f = 2*sigma*sigma + alpha*dx*dx ;
	Double_t cruijff = par[0]*exp(-dx*dx/f);

	Double_t dxRes = (x[0] - par[6]);
	Double_t sigmaRes = par[7];
	Double_t fRes = 2*sigmaRes*sigmaRes;
	Double_t gaussRes = exp( -dxRes*dxRes/fRes );

	return cruijff*gaussRes;
}

Double_t Cruijff2Gauss( Double_t *x, Double_t* par )
{
	// -- par[0] = A, par[1] = m, par[2] = sigma, par[3] = alphaL, par[4] = alphaR, par[5] = mRes, par[6] = sigmaRes -- //
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
	Double_t cruijff = par[0]*exp(-dx*dx/f);

	Double_t dxRes = (x[0] - par[5]);
	Double_t sigmaRes = par[6];
	Double_t fRes = 2*sigmaRes*sigmaRes;
	Double_t gaussRes = exp( -dxRes*dxRes/fRes );

	return cruijff*gaussRes;
}

class FittingTool
{
public:
	TString channelType;
	TString genType;
	TString FileName;
	vector<TH1D*> vec_Hist;

	Double_t MassBinEdges_Double[nBin+1];
	TH1D* h_NormChi2;
	TH1D* h_RMS;
	TH1D* h_Mean;
	TH1D* h_Resolution;

	FittingTool()
	{
		TH1::AddDirectory(kFALSE);
		for(Int_t i=0; i<nBin+1; i++)
			MassBinEdges_Double[i] = (Double_t)MassBinEdges[i];

	}

	void Set_channelType( TString value )
	{
		this->channelType = value;
	}

	void Set_genType( TString value )
	{
		this->genType = value;
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
		this->h_NormChi2 = new TH1D("h_NormChi2_"+this->channelType, "", nBin, MassBinEdges_Double);
		this->h_RMS = new TH1D("h_RMS_"+this->channelType, "", nBin, MassBinEdges_Double);
		this->h_Mean = new TH1D("h_Mean_"+this->channelType, "", nBin, MassBinEdges_Double);
		this->h_Resolution = new TH1D("h_Resolution_"+this->channelType, "", nBin, MassBinEdges_Double);

		// -- get histograms -- //
		this->FileName = TString::Format("Local/ROOTFile_MakeHist_MassResolution_%s.root", this->channelType.Data() );

		for(Int_t i=0; i<nBin; i++)
		{
			TString TStr_MassRange = TString::Format("M%dto%d", MassBinEdges[i], MassBinEdges[i+1]);
			TString HistName = TString::Format("DY%s/h_RelDiff_%s_%s", this->channelType.Data(), this->genType.Data(), TStr_MassRange.Data());

			TH1D* h_temp = Get_Hist( this->FileName, HistName );
			// h_temp->Rebin(5);
			this->vec_Hist.push_back( h_temp );
		}
	}

	void Fit_GetResolution( Int_t i_hist, TH1D* h )
	{
		TH1D* h_Clone = (TH1D*)h->Clone();

		TString fitFuncName = "Cruijff";
		if( this->channelType == "EE" && this->genType == "dressed" ) fitFuncName = "Cruijff2Gauss";
		if( this->channelType == "MuMu" && this->genType == "postFSR" ) fitFuncName = "Cruijff2";

		Int_t nParam = -999;
		if( fitFuncName == "Cruijff" ) nParam = 6;
		if( fitFuncName == "CruijffGauss" ) nParam = 8;
		if( fitFuncName == "Cruijff2" ) nParam = 5;
		if( fitFuncName == "Cruijff2Gauss" ) nParam = 7;

		// Double_t rangeMax = 0.1;
		Double_t rangeMax = 2.0*h->GetRMS();

		// -- guassian fit to get the initial values for cruijff fit -- //
		TF1* fitFunc_Gauss = new TF1("gaus", "gaus", (-1)*rangeMax, rangeMax);
		fitFunc_Gauss->SetParameters(h->Integral(), h->GetMean(), h->GetRMS());
		h->Fit( fitFunc_Gauss, "R", "", (-1)*rangeMax, rangeMax );

		Double_t constInit = fitFunc_Gauss->GetParameter(0);
		Double_t meanInit = fitFunc_Gauss->GetParameter(1);
		Double_t sigmaInit = fitFunc_Gauss->GetParameter(2);

		TF1* func = NULL;
		if( fitFuncName == "Cruijff" )
		{
			func = new TF1(fitFuncName, Cruijff, (-1)*rangeMax, rangeMax, nParam);
			func->SetParameter(0, constInit); func->SetParName(0, "const");
			func->SetParameter(1, meanInit); func->SetParName(1, "mean");
			func->SetParameter(2, sigmaInit); func->SetParName(2, "sigmaL");
			func->SetParameter(3, sigmaInit); func->SetParName(3, "sigmaR");
			func->SetParameter(4, 0); func->SetParName(4, "alphaL");
			func->SetParameter(5, 0); func->SetParName(5, "alphaR");
		}
		else if( fitFuncName == "Cruijff2")
		{
			func = new TF1(fitFuncName, Cruijff2, (-1)*rangeMax, rangeMax, nParam);
			func->SetParameter(0, h->Integral()); func->SetParName(0, "const");
			func->SetParameter(1, meanInit); func->SetParName(1, "mean");
			func->SetParameter(2, sigmaInit); func->SetParName(2, "sigma");
			func->SetParameter(3, 0.1); func->SetParName(3, "alphaL");
			func->SetParameter(4, 0.1); func->SetParName(4, "alphaR");
		}
		else if( fitFuncName == "CruijffGauss" )
		{
			func = new TF1(fitFuncName, CruijffGauss, (-1)*rangeMax, rangeMax, nParam);
			func->SetParameter(0, constInit); func->SetParName(0, "const");
			func->SetParameter(1, meanInit); func->SetParName(1, "mean");
			func->SetParameter(2, sigmaInit); func->SetParName(2, "sigmaL");
			func->SetParameter(3, sigmaInit); func->SetParName(3, "sigmaR");
			func->SetParameter(4, 0); func->SetParName(4, "alphaL");
			func->SetParameter(5, 0); func->SetParName(5, "alphaR");
			func->SetParameter(6, 0); func->SetParName(6, "mRes"); // -- initial mean value = 0: describe resolution effect at each point -- //
			func->SetParameter(7, 0.01); func->SetParName(7, "sigmaRes");
		}
		else if( fitFuncName == "Cruijff2Gauss" )
		{
			func = new TF1(fitFuncName, Cruijff2Gauss, (-1)*rangeMax, rangeMax, nParam);
			func->SetParameter(0, constInit); func->SetParName(0, "const");
			func->SetParameter(1, meanInit); func->SetParName(1, "mean");
			func->SetParameter(2, sigmaInit); func->SetParName(2, "sigma");
			func->SetParameter(3, 0); func->SetParName(3, "alphaL");
			func->SetParameter(4, 0); func->SetParName(4, "alphaR");
			func->SetParameter(5, 0); func->SetParName(5, "mRes"); // -- initial mean value = 0: describe resolution effect at each point -- //
			func->SetParameter(6, sigmaInit); func->SetParName(6, "sigmaRes");
		}

		h->Fit( fitFuncName, "R", "", (-1)*rangeMax, rangeMax );

		TF1 *FitFunc = h->GetFunction( fitFuncName );
		Double_t chi2 = FitFunc->GetChisquare();
		Double_t nDOF = FitFunc->GetNDF();
		Double_t NormChi2 = chi2/nDOF;

		Double_t RMS = h->GetRMS();
		Double_t RMS_err = h->GetRMSError();

		Double_t mean = FitFunc->GetParameter(1);
		Double_t mean_err = FitFunc->GetParError(1);

		Double_t sigma = FitFunc->GetParameter(2);
		Double_t sigma_err = FitFunc->GetParError(2);

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
		TString canvasName = h->GetName();
		canvasName.ReplaceAll("h_", "Local/c_"+this->channelType+"_");

		HistInfo* Hist = new HistInfo( kBlack, "temp", h_Clone );
		fitFunc_Gauss->SetLineColor(kRed);
		fitFunc_Gauss->SetMarkerColorAlpha(kWhite, 0);
		fitFunc_Gauss->SetFillColorAlpha(kWhite, 0);
		FitFunc->SetLineColor(kBlue);
		FitFunc->SetMarkerColorAlpha(kWhite, 0);
		FitFunc->SetFillColorAlpha(kWhite, 0);

		SetCanvas_Square( c, canvasName );
		c->cd();
		Hist->h->Draw("EPSAME");
		fitFunc_Gauss->Draw("SAME");
		FitFunc->Draw("SAME");

		TLegend *legend;
		SetLegend(legend, 0.80, 0.75, 0.95, 0.95 );
		legend->AddEntry( fitFunc_Gauss, "Gaussian" );
		legend->AddEntry( FitFunc, fitFuncName );
		legend->Draw();

		SetHistFormat_SinglePad(Hist->h, "(M^{Reco}-M^{Gen})/M^{Gen}", "Entries");
		Hist->h->GetXaxis()->SetRangeUser( (-1)*rangeMax, rangeMax );

		TLatex latex;
		Latex_Simulation( latex );
		TString MassRange = TString::Format("%d < M^{Gen} < %d GeV", MassBinEdges[i_hist], MassBinEdges[i_hist+1]);
		latex.DrawLatexNDC(0.16, 0.91, "#font[42]{#scale[0.7]{"+MassRange+"}}");
		TString normChi2Info = TString::Format("#chi2/ndof = %.2lf/%.0lf = %.2lf", chi2, nDOF, NormChi2);
		latex.DrawLatexNDC(0.16, 0.88, "#font[42]{#scale[0.5]{"+normChi2Info+"}}");

		if( fitFuncName == "Cruijff2" )
		{
			TString meanInfo = TString::Format("mean = %.5lf", mean);
			latex.DrawLatexNDC(0.16, 0.84, "#font[42]{#scale[0.5]{"+meanInfo+"}}");
			TString sigmaInfo = TString::Format("resolution(#sigma) = %.3lf", sigma);
			latex.DrawLatexNDC(0.16, 0.80, "#font[42]{#scale[0.5]{"+sigmaInfo+"}}");

			Double_t alphaL = FitFunc->GetParameter(3);
			Double_t alphaR = FitFunc->GetParameter(4);

			TString alphaLInfo = TString::Format("#alpha_{L} = %.3lf", alphaL);
			latex.DrawLatexNDC(0.16, 0.76, "#font[42]{#scale[0.5]{"+alphaLInfo+"}}");
			TString alphaRInfo = TString::Format("#alpha_{R} = %.3lf", alphaR);
			latex.DrawLatexNDC(0.16, 0.72, "#font[42]{#scale[0.5]{"+alphaRInfo+"}}");
		}

		if( fitFuncName == "Cruijff" )
		{
			TString meanInfo = TString::Format("mean = %.5lf", mean);
			latex.DrawLatexNDC(0.16, 0.84, "#font[42]{#scale[0.5]{"+meanInfo+"}}");

			Double_t sigmaL = FitFunc->GetParameter(2);
			TString sigmaLInfo = TString::Format("resolution(#sigma_{L}) = %.3lf", sigmaL);
			latex.DrawLatexNDC(0.16, 0.80, "#font[42]{#scale[0.5]{"+sigmaLInfo+"}}");

			Double_t sigmaR = FitFunc->GetParameter(3);
			TString sigmaRInfo = TString::Format("resolution(#sigma_{R}) = %.3lf", sigmaR);
			latex.DrawLatexNDC(0.16, 0.76, "#font[42]{#scale[0.5]{"+sigmaRInfo+"}}");

			Double_t alphaL = FitFunc->GetParameter(4);
			Double_t alphaR = FitFunc->GetParameter(5);

			TString alphaLInfo = TString::Format("#alpha_{L} = %.3lf", alphaL);
			latex.DrawLatexNDC(0.16, 0.72, "#font[42]{#scale[0.5]{"+alphaLInfo+"}}");
			TString alphaRInfo = TString::Format("#alpha_{R} = %.3lf", alphaR);
			latex.DrawLatexNDC(0.16, 0.68, "#font[42]{#scale[0.5]{"+alphaRInfo+"}}");
		}

		if( fitFuncName == "CruijffGauss" )
		{
			TString meanInfo = TString::Format("mean = %.5lf", mean);
			latex.DrawLatexNDC(0.16, 0.84, "#font[42]{#scale[0.5]{"+meanInfo+"}}");

			Double_t sigmaL = FitFunc->GetParameter(2);
			TString sigmaLInfo = TString::Format("resolution(#sigma_{L}) = %.3lf", sigmaL);
			latex.DrawLatexNDC(0.16, 0.80, "#font[42]{#scale[0.5]{"+sigmaLInfo+"}}");

			Double_t sigmaR = FitFunc->GetParameter(3);
			TString sigmaRInfo = TString::Format("resolution(#sigma_{R}) = %.3lf", sigmaR);
			latex.DrawLatexNDC(0.16, 0.76, "#font[42]{#scale[0.5]{"+sigmaRInfo+"}}");

			Double_t alphaL = FitFunc->GetParameter(4);
			Double_t alphaR = FitFunc->GetParameter(5);

			TString alphaLInfo = TString::Format("#alpha_{L} = %.3lf", alphaL);
			latex.DrawLatexNDC(0.16, 0.72, "#font[42]{#scale[0.5]{"+alphaLInfo+"}}");
			TString alphaRInfo = TString::Format("#alpha_{R} = %.3lf", alphaR);
			latex.DrawLatexNDC(0.16, 0.68, "#font[42]{#scale[0.5]{"+alphaRInfo+"}}");

			Double_t mRes = FitFunc->GetParameter(6);
			TString mResInfo = TString::Format("m_{Res} = %.3lf", mRes);
			latex.DrawLatexNDC(0.16, 0.64, "#font[42]{#scale[0.5]{"+mResInfo+"}}");

			Double_t sigmaRes = FitFunc->GetParameter(7);
			TString sigmaResInfo = TString::Format("#sigma_{Res} = %.3lf", sigmaRes);
			latex.DrawLatexNDC(0.16, 0.60, "#font[42]{#scale[0.5]{"+sigmaResInfo+"}}");
		}

		if( fitFuncName == "Cruijff2Gauss" )
		{
			TString meanInfo = TString::Format("mean = %.5lf", mean);
			latex.DrawLatexNDC(0.16, 0.84, "#font[42]{#scale[0.5]{"+meanInfo+"}}");

			Double_t sigma = FitFunc->GetParameter(2);
			TString sigmaInfo = TString::Format("resolution(#sigma) = %.3lf", sigma);
			latex.DrawLatexNDC(0.16, 0.80, "#font[42]{#scale[0.5]{"+sigmaInfo+"}}");

			Double_t alphaL = FitFunc->GetParameter(3);
			Double_t alphaR = FitFunc->GetParameter(4);

			TString alphaLInfo = TString::Format("#alpha_{L} = %.3lf", alphaL);
			latex.DrawLatexNDC(0.16, 0.72, "#font[42]{#scale[0.5]{"+alphaLInfo+"}}");
			TString alphaRInfo = TString::Format("#alpha_{R} = %.3lf", alphaR);
			latex.DrawLatexNDC(0.16, 0.68, "#font[42]{#scale[0.5]{"+alphaRInfo+"}}");

			Double_t mRes = FitFunc->GetParameter(5);
			TString mResInfo = TString::Format("m_{Res} = %.3lf", mRes);
			latex.DrawLatexNDC(0.16, 0.64, "#font[42]{#scale[0.5]{"+mResInfo+"}}");

			Double_t sigmaRes = FitFunc->GetParameter(6);
			TString sigmaResInfo = TString::Format("#sigma_{Res} = %.3lf", sigmaRes);
			latex.DrawLatexNDC(0.16, 0.60, "#font[42]{#scale[0.5]{"+sigmaResInfo+"}}");
		}

		c->SaveAs(".pdf");
	}

};

void GetResolution()
{
	FittingTool *tool_EE = new FittingTool();
	tool_EE->Set_channelType("EE");
	tool_EE->Set_genType( "dressed" );
	tool_EE->GetResolution();

	FittingTool *tool_MuMu = new FittingTool();
	tool_MuMu->Set_channelType("MuMu");
	tool_MuMu->Set_genType( "postFSR" );
	tool_MuMu->GetResolution();

	TFile *f_output = TFile::Open("ROOTFile_GetResolution.root", "RECREATE");
	tool_EE->Save( f_output );
	tool_MuMu->Save( f_output );
	f_output->Close();
}