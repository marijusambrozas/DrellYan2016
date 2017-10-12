#include <BinOptimization/MassResolution/Test/v02_preFSR_vs_Reco/MakeHist_MassResolution.h>

void MakeHist_MassResolution_EE( TString FileName_ROOTFileList, TString Tag, Int_t IsMC, Double_t NormFactor )
{
	HistProducer *producer = new HistProducer( FileName_ROOTFileList, Tag, IsMC );
	producer->Set_NormFactor( NormFactor );
	producer->Set_ChannelType( "EE" );
	producer->Produce();

	TFile *f_output = TFile::Open("ROOTFile_MakeHist_MassResolution_EE.root", "RECREATE");
	producer->Save( f_output );
	f_output->Close();
}