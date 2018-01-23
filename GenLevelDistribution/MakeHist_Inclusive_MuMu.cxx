#include <GenLevelDistribution/MakeHist_Inclusive.h>

void MakeHist_Inclusive_MuMu( TString FileName_ROOTFileList, TString Tag, Int_t IsMC, Double_t NormFactor )
{
	TString channelType = "MuMu";
	HistProducer *producer = new HistProducer( FileName_ROOTFileList, Tag, IsMC );
	producer->Set_NormFactor( NormFactor );
	producer->Set_ChannelType( channelType );
	producer->Produce();

	TFile *f_output = TFile::Open("ROOTFile_MakeHist_Inclusive_"+channelType+".root", "RECREATE");
	producer->Save( f_output );
	f_output->Close();
}