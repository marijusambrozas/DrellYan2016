#include <GenLevelDistribution/MakeHist_GenLevel.h>

void MakeHist_GenLevel_EE( TString FileName_ROOTFileList, TString Tag, Int_t IsMC, Double_t NormFactor )
{
	TString channelType = "MuMu";
	HistProducer *producer = new HistProducer( FileName_ROOTFileList, Tag, IsMC );
	producer->Set_NormFactor( NormFactor );
	producer->Set_ChannelType( channelType );
	producer->Produce();

	TFile *f_output = TFile::Open("ROOTFile_MakeHist_GenLevel_"+channelType+".root", "RECREATE");
	producer->Save( f_output );
	f_output->Close();
}