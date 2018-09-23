#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLorentzVector.h>
#include <TStopwatch.h>
#include <TTimeStamp.h>
#include <TString.h>
#include <TLegend.h>
#include <THStack.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TAttMarker.h>
#include <TF1.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TApplication.h>
#include <vector>
#include <TMath.h>
#include <THistPainter.h>
#include <TFormula.h>
#include <iostream>
#include <sstream>

// -- Customized Analyzer for Drel-Yan Analysis -- //
#include "./header/DYAnalyzer.h"
#include "./header/SelectedX.h"
#include "./header/myProgressBar_t.h"
#include "./header/LocalFileMgr.h"
#include "./etc/RoccoR/RoccoR.cc"

void EE_HistMaker ( TString type, TString HLTname );
void MuMu_HistMaker ( TString type, Bool_t SwitchROCCORR, TString HLTname );
void EMu_HistMaker ( TString type, Bool_t SwitchROCCORR, TString HLTname );
void MuMu_merge();

// -- Drell-Yan mass bins -- //
const Double_t massbins[44] = {15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 64, 68, 72, 76, 81, 86, 91, 96, 101, 106, 110, 115, 120, 126,
                               133, 141, 150, 160, 171, 185, 200, 220, 243, 273, 320, 380, 440, 510, 600, 700, 830, 1000, 1500, 3000};


void HistMaker ( TString whichX = "", TString Type = "", Bool_t SwitchROCCORR = kFALSE, TString HLTname = "DEFAULT" )
{
    TString HLT;
    Int_t Xselected = 0;
    TString type = whichX+"_"+Type;
    if ( whichX.Contains("EE") || whichX.Contains("ee") )
    {
        Xselected++;
        if ( HLTname == "DEFAULT" ) HLT = "Ele23Ele12";
        else HLT = HLTname;
        cout << "\n*******      EE_HistMaker ( " << type << " )      *******" << endl;
        EE_HistMaker( type, HLT );
    }
    if ( whichX.Contains("MuMu") || whichX.Contains("mumu") || whichX.Contains("MUMU") )
    {
        Xselected++;
        if ( HLTname == "DEFAULT" ) HLT = "IsoMu24_OR_IsoTkMu24";
        else HLT = HLTname;
        cout << "\n*****  MuMu_HistMaker ( " << type << " )  *****" << endl;
        MuMu_HistMaker( type, SwitchROCCORR, HLT );
    }
    if ( whichX.Contains("EMu") || whichX.Contains("emu") || whichX.Contains("Emu") || whichX.Contains("eMu") || whichX.Contains("EMU") )
    {
        Xselected++;
        if ( HLTname == "DEFAULT" ) HLT = "IsoMu24_OR_IsoTkMu24";
        else HLT = HLTname;
        cout << "\n*****   EMu_HistMaker ( " << type << " )  *****" << endl;
        EMu_HistMaker( type, SwitchROCCORR, HLT );
    }
    if ( whichX.Contains("MuMu_merge") || whichX.Contains("MuMu_Merge") )
    {
        Xselected++;
        cout << "\n*********   MuMu_merge  *********" << endl;
        MuMu_merge();
    }
    if ( Xselected == 0 ) cout << "Wrong arument! \nType in: >> .x HistMaker.C+(\"whichX\", \"whichProcess\", SwitchROCCORR)" << endl;

} // End of HistMaker()


/// ----------------------------- Electron Channel ------------------------------ ///
void EE_HistMaker ( TString type, TString HLTname )
{
    if ( !type.Length() )
    {
        cout << "Error: no type specified!" << endl;
        return;
    }

    // -- Run2016 luminosity [/pb] -- //
    Double_t L_B2F = 19721.0, L_G2H = 16146.0, L_B2H = 35867.0, L = 0;
    L = L_B2H;

    TTimeStamp ts_start;
    cout << "[Start Time(local time): " << ts_start.AsString("l") << "]" << endl;
    TStopwatch totaltime;
    totaltime.Start();

    LocalFileMgr Mgr;
    vector<SelProc_t> Processes = Mgr.FindProc( type );
    TFile *f;
    TString OutputDir;
    Bool_t isBkgFull = kFALSE;  // To tell if this is the _EE_Bkg_Full process (it is handled differently)

    if ( !Processes.size() )
    {
        cout << "Error: no processes!" << endl;
        return;
    }

    if ( Processes[0] == _EE_Bkg_Full )
    {
        isBkgFull = kTRUE;
        Mgr.GetProc( _EE_Bkg_Full );
        // -- Output ROOTFile -- //
        OutputDir = Mgr.HistLocation;
        f = new TFile( OutputDir+"Hist_"+Mgr.Procname[_EE_Bkg_Full]+".root", "RECREATE" );
        Processes.clear();
        Processes.push_back( _EE_DYTauTau_Full );
        Processes.push_back( _EE_ttbar_Full );
        Processes.push_back( _EE_VVnST );
        Processes.push_back( _EE_WJets );
        Processes.push_back( _EE_QCDEMEnriched_Full );
    }

    for ( Int_t i_proc=0; i_proc<((int)(Processes.size())); i_proc++ )
    {
        if ( Processes[i_proc] <= _EndOf_MuMu || Processes[i_proc] > _EndOf_EE )
        {
            cout << "Error: process " << Mgr.Procname[Processes[i_proc]] << " is not EE!" << endl;
            continue;
        }

        Mgr.GetProc( Processes[i_proc] );

        // -- Output ROOTFile -- //
        if ( isBkgFull == kFALSE )
        {
            OutputDir = Mgr.HistLocation;
            f = new TFile( OutputDir+"Hist_"+Mgr.Procname[Mgr.CurrentProc]+".root", "RECREATE" );
        }

        cout << "Process: " << Mgr.Procname[Mgr.CurrentProc] << endl;
        cout << "Type: " << Mgr.Type << endl;
        cout << "DATA location: " << Mgr.BaseLocation << endl;
        cout << "Output directory: " << OutputDir << endl;

        TStopwatch totaltime;
        totaltime.Start();

        DYAnalyzer *analyzer = new DYAnalyzer( HLTname );

        // -- For PU re-weighting -- //
        analyzer->SetupPileUpReWeighting_80X( Mgr.isMC, "ROOTFile_PUReWeight_80X_v20170817_64mb.root" );

        // -- For efficiency SF -- //
        analyzer->SetupEfficiencyScaleFactor_electron();        

        // -- Creating Histograms -- //
        TH1D *h_mass_fine_before_PUCorr = new TH1D("h_mass_fine_before_PUCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 10000, 0, 10000);
        TH1D *h_mass_fine_before_EffCorr = new TH1D("h_mass_fine_before_EffCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 10000, 0, 10000);
        TH1D *h_mass_fine = new TH1D("h_mass_fine_"+Mgr.Procname[Mgr.CurrentProc], "", 10000, 0, 10000);
        TH1D *h_mass = new TH1D("h_mass_"+Mgr.Procname[Mgr.CurrentProc], "", 43, massbins);
        TH1D *h_Pt = new TH1D("h_Pt_"+Mgr.Procname[Mgr.CurrentProc], "", 300, 0, 600);
        TH1D *h_rapi = new TH1D("h_rapi_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -5, 5);

        TH1D *h_nPU_beforePUCorr = new TH1D( "h_nPU_before_PUCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 50, 0, 50 );
        TH1D *h_nPU_beforeEffCorr = new TH1D( "h_nPU_before_EffCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 50, 0, 50 );
        TH1D *h_nPU = new TH1D( "h_nPU_"+Mgr.Procname[Mgr.CurrentProc], "", 50, 0, 50 );

        TH1D *h_pT = new TH1D("h_pT_"+Mgr.Procname[Mgr.CurrentProc], "", 300, 0, 600);
        TH1D *h_eta = new TH1D("h_eta_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -5, 5);
        TH1D *h_phi = new TH1D("h_phi_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -5, 5);

        //Loop for all samples
        const Int_t Ntup = Mgr.FileLocation.size();
        for(Int_t i_tup = 0; i_tup<Ntup; i_tup++)
        {
            TStopwatch looptime;
            looptime.Start();

            cout << "\t<" << Mgr.Tag[i_tup] << ">" << endl;

            TChain *chain = new TChain( Mgr.TreeName[i_tup] );
            chain->Add( Mgr.FullLocation[i_tup] );

            SelectedEE_t *EE = new SelectedEE_t();
            EE->CreateFromChain( chain );

            Int_t NEvents = chain->GetEntries();
            if ( NEvents != Mgr.nEvents[i_tup] ) cout << "\tEvent numbers do not match!!!" << endl;
            cout << "\t[Total events: " << Mgr.Wsum[i_tup] << "]" << endl;
            cout << "\t[Selected Events: " << NEvents << "]" << endl;

            myProgressBar_t bar( NEvents );

            for ( Int_t i=0; i<NEvents; i++ )
            {
                EE->GetEvent(i);

                // -- Bring the weights -- //
                Double_t GenWeight = EE->GENEvt_weight;

                // -- Pileup-Reweighting -- //
                Double_t PUWeight = 1;
                if( Mgr.isMC == kTRUE ) PUWeight = analyzer->PileUpWeightValue_80X( EE->nPileUp );

                // -- efficiency weights -- //
                Double_t effweight = 1;

                // -- Normalization -- //
                Double_t TotWeight = GenWeight;
                if( Mgr.isMC == kTRUE ) TotWeight = ( L * Mgr.Xsec[i_tup] / Mgr.Wsum[i_tup] ) * GenWeight;

                if( EE->isSelPassed == kTRUE )
                {
                    TLorentzVector ele1, ele2;
                    ele1.SetPtEtaPhiE( EE->Electron_pT->at(0), EE->Electron_eta->at(0), EE->Electron_phi->at(0), EE->Electron_Energy->at(0) );
                    ele1.SetPtEtaPhiE( EE->Electron_pT->at(1), EE->Electron_eta->at(1), EE->Electron_phi->at(1), EE->Electron_Energy->at(1) );
                    Double_t reco_Pt = ( ele1 + ele2 ).Pt();
                    Double_t reco_rapi = ( ele1 + ele2 ).Rapidity();

                    // -- Apply efficiency correcion -- //
                    if( Mgr.isMC == kTRUE )
                            effweight = analyzer->EfficiencySF_EventWeight_electron( EE );

                    h_mass_fine_before_PUCorr->Fill( EE->Electron_InvM, TotWeight );
                    h_mass_fine_before_EffCorr->Fill( EE->Electron_InvM, TotWeight * PUWeight );
                    h_mass_fine->Fill( EE->Electron_InvM, TotWeight * PUWeight * effweight );
                    h_mass->Fill( EE->Electron_InvM, TotWeight * PUWeight * effweight );
                    h_Pt->Fill( reco_Pt, TotWeight * PUWeight * effweight );
                    h_rapi->Fill( reco_rapi, TotWeight * PUWeight * effweight );

                    h_nPU_beforePUCorr->Fill( EE->nPileUp, TotWeight );
                    h_nPU_beforeEffCorr->Fill( EE->nPileUp, TotWeight * PUWeight );
                    h_nPU->Fill( EE->nPileUp, TotWeight * PUWeight * effweight );

                    h_pT->Fill( EE->Electron_pT->at(0), TotWeight * PUWeight * effweight );
                    h_pT->Fill( EE->Electron_pT->at(1), TotWeight * PUWeight * effweight );
                    h_eta->Fill( EE->Electron_eta->at(0), TotWeight * PUWeight * effweight );
                    h_eta->Fill( EE->Electron_eta->at(1), TotWeight * PUWeight * effweight );
                    h_phi->Fill( EE->Electron_phi->at(0), TotWeight * PUWeight * effweight );
                    h_phi->Fill( EE->Electron_phi->at(1), TotWeight * PUWeight * effweight );

                } // End of event selection
                bar.Draw(i);

            } // End of event iteration

            if( Mgr.isMC == kTRUE ) printf( "\tNormalization factor: %.8f\n", L*Mgr.Xsec[i_tup]/Mgr.Wsum[i_tup] );

            Double_t LoopRunTime = looptime.CpuTime();
            cout << "\tLoop RunTime(" << Mgr.Tag[i_tup] << "): " << LoopRunTime << " seconds\n" << endl;

        }// End of i_tup iteration

        f->cd();
        cout << "\tWriting into file...";

        h_mass_fine_before_PUCorr->Write();
        h_mass_fine_before_EffCorr->Write();
        h_mass_fine->Write();
        h_mass->Write();
        h_Pt->Write();
        h_rapi->Write();

        h_nPU_beforePUCorr->Write();
        h_nPU_beforeEffCorr->Write();
        h_nPU->Write();

        h_pT->Write();
        h_eta->Write();
        h_phi->Write();

        cout << " Finished." << endl;

    }// End of i_proc iteration

    f->Close();
    if ( isBkgFull == kTRUE )
    {
        if ( !f->IsOpen() ) cout << "File Hist_" << Mgr.Procname[_EE_Bkg_Full] << ".root has been closed successfully.\n" << endl;
        else cout << "FILE Hist_" << Mgr.Procname[_EE_Bkg_Full] << ".root COULD NOT BE CLOSED!\n" << endl;
    }
    else
    {
        if ( !f->IsOpen() ) cout << "File Hist_" << Mgr.Procname[Mgr.CurrentProc] << ".root has been closed successfully.\n" << endl;
        else cout << "FILE Hist_" << Mgr.Procname[Mgr.CurrentProc] << ".root COULD NOT BE CLOSED!\n" << endl;
    }

    Double_t TotalRunTime = totaltime.CpuTime();
    cout << "Total RunTime: " << TotalRunTime << " seconds" << endl;

    TTimeStamp ts_end;
    cout << "[End Time(local time): " << ts_end.AsString("l") << "]" << endl;

} // End of EE_HistMaker()


/// -------------------------------- Muon Channel ------------------------------------ ///
void MuMu_HistMaker ( TString type, Bool_t SwitchROCCORR, TString HLTname )
{
    if ( !type.Length() )
    {
        cout << "Error: no type specified!" << endl;
        return;
    }

    // -- Run2016 luminosity [/pb] -- //
    Double_t L_B2F = 19721.0, L_G2H = 16146.0, L_B2H = 35867.0, L = 0;
    L = L_B2H;

    TTimeStamp ts_start;
    cout << "[Start Time(local time): " << ts_start.AsString("l") << "]" << endl;
    TStopwatch totaltime;
    totaltime.Start();

    LocalFileMgr Mgr;
    vector<SelProc_t> Processes = Mgr.FindProc( type );
    if ( SwitchROCCORR == kTRUE ) Mgr.SwitchROCCORR();  // If kTRUE, this will go through events that were selected without Rochester correction applied
    TFile *f;
    TString OutputDir;
    Bool_t isBkgFull = kFALSE;  // To tell if this is the _EE_Bkg_Full process (it is handled differently)

    if ( !Processes.size() )
    {
        cout << "Error: no processes!" << endl;
        return;
    }

    if ( Processes[0] == _MuMu_Bkg_Full )
    {
        isBkgFull = kTRUE;
        Mgr.GetProc( _MuMu_Bkg_Full );
        // -- Output ROOTFile -- //
        OutputDir = Mgr.HistLocation;
        f = new TFile( OutputDir+"Hist_"+Mgr.Procname[_MuMu_Bkg_Full]+".root", "RECREATE" );
        Processes.clear();
        Processes.push_back( _MuMu_DYTauTau_Full );
        Processes.push_back( _MuMu_ttbar_Full );
        Processes.push_back( _MuMu_VVnST );
        Processes.push_back( _MuMu_WJets );
        Processes.push_back( _MuMu_QCDMuEnriched_Full );
    }

    for ( Int_t i_proc; i_proc<((int)(Processes.size())); i_proc++ )
    {
        if ( Processes[i_proc] > _EndOf_MuMu )
        {
            cout << "Error: process " << Mgr.Procname[Processes[i_proc]] << " is not MuMu!" << endl;
            continue;
        }

        Mgr.GetProc( Processes[i_proc] );

        // -- Output ROOTFile -- //
        if ( isBkgFull == kFALSE )
        {
            OutputDir = Mgr.HistLocation;
            f = new TFile( OutputDir+"Hist_"+Mgr.Procname[Mgr.CurrentProc]+".root", "RECREATE" );
        }

        cout << "Process: " << Mgr.Procname[Mgr.CurrentProc] << endl;
        cout << "Type: " << Mgr.Type << endl;
        cout << "DATA location: " << Mgr.BaseLocation << endl;
        cout << "Output directory: " << OutputDir << endl;

        TStopwatch totaltime;
        totaltime.Start();

        DYAnalyzer *analyzer = new DYAnalyzer( HLTname );

        // -- For PU re-weighting -- //
        analyzer->SetupPileUpReWeighting_80X( Mgr.isMC, "ROOTFile_PUReWeight_80X_v20170817_64mb.root" );

        // -- For Rochester correction -- //
        TRandom3 *r1 = new TRandom3(0);

        // -- For efficiency SF -- //
        analyzer->SetupEfficiencyScaleFactor_BtoF();
        analyzer->SetupEfficiencyScaleFactor_GtoH();

        // -- Creating Histograms -- //
        TH1D *h_mass_fine_before_PUCorr = new TH1D( "h_mass_fine_before_PUCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 10000, 0, 10000 );
        TH1D *h_mass_fine_before_RoccoR = new TH1D( "h_mass_fine_before_RoccoR_"+Mgr.Procname[Mgr.CurrentProc], "", 10000, 0, 10000 );
        TH1D *h_mass_fine_before_EffCorr = new TH1D( "h_mass_fine_before_EffCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 10000, 0, 10000 );
        TH1D *h_mass_fine = new TH1D( "h_mass_fine_"+Mgr.Procname[Mgr.CurrentProc], "", 10000, 0, 10000 );
        TH1D *h_mass = new TH1D( "h_mass_"+Mgr.Procname[Mgr.CurrentProc], "", 43, massbins );
        TH1D *h_Pt = new TH1D( "h_Pt_"+Mgr.Procname[Mgr.CurrentProc], "", 300, 0, 600 );
        TH1D *h_rapi = new TH1D( "h_rapi_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -5, 5 );

        TH1D *h_nPU_beforePUCorr = new TH1D( "h_nPU_before_PUCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 50, 0, 50 );
        TH1D *h_nPU_beforeEffCorr = new TH1D( "h_nPU_before_EffCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 50, 0, 50 );
        TH1D *h_nPU = new TH1D( "h_nPU_"+Mgr.Procname[Mgr.CurrentProc], "", 50, 0, 50 );

        TH1D *h_pT = new TH1D( "h_pT_"+Mgr.Procname[Mgr.CurrentProc], "", 300, 0, 600 );
        TH1D *h_eta = new TH1D( "h_eta_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -5, 5 );
        TH1D *h_phi = new TH1D( "h_phi_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -5, 5 );

        //Loop for all samples
        const Int_t Ntup = Mgr.FileLocation.size();
        for(Int_t i_tup = 0; i_tup<Ntup; i_tup++)
        {
            TStopwatch looptime;
            looptime.Start();

            cout << "\t<" << Mgr.Tag[i_tup] << ">" << endl;

            TChain *chain = new TChain( Mgr.TreeName[i_tup] );
            chain->Add( Mgr.FullLocation[i_tup] );

            SelectedMuMu_t *MuMu = new SelectedMuMu_t();
            MuMu->CreateFromChain( chain );

            RoccoR rc("./etc/RoccoR/rcdata.2016.v3");

            Int_t NEvents = chain->GetEntries();
            if ( NEvents != Mgr.nEvents[i_tup] ) cout << "\tEvent numbers do not match!!!" << endl;
            cout << "\t[Total events: " << Mgr.Wsum[i_tup] << "]" << endl;
            cout << "\t[Selected Events: " << NEvents << "]" << endl;

            myProgressBar_t bar( NEvents );

            for(Int_t i=0; i<NEvents; i++)
            {
                MuMu->GetEvent(i);

                Double_t GenWeight = MuMu->GENEvt_weight;

                // -- Pileup-Reweighting -- //
                Double_t PUWeight = 1;
                if( Mgr.isMC == kTRUE ) PUWeight = analyzer->PileUpWeightValue_80X( MuMu->nPileUp );

                // -- efficiency weights -- //
                Double_t weight1 = 0, weight2 = 0, effweight = 1;

                // -- Normalization -- //
                Double_t TotWeight = GenWeight;
                if( Mgr.isMC == kTRUE ) TotWeight = ( L * Mgr.Xsec[i_tup] / Mgr.Wsum[i_tup] ) * GenWeight;

                // -- Rochester correction -- //
                for ( Int_t iter = 0; i<=1; i++ )
                {
                    Double_t rndm[2], SF=0; r1->RndmArray(2, rndm);
                    Int_t s, m;

                    if( Mgr.Type == "DATA" )
                            SF = rc.kScaleDT(MuMu->Muon_charge->at(iter), MuMu->Muon_TuneP_pT->at(iter), MuMu->Muon_TuneP_eta->at(iter),
                                             MuMu->Muon_TuneP_phi->at(iter), s=0, m=0);
                    else
                            SF = rc.kScaleAndSmearMC(MuMu->Muon_charge->at(iter), MuMu->Muon_TuneP_pT->at(iter),
                                                     MuMu->Muon_TuneP_eta->at(iter), MuMu->Muon_TuneP_phi->at(iter),
                                                     MuMu->Muon_trackerLayers->at(iter), rndm[0], rndm[1], s=0, m=0);

                    MuMu->Muon_TuneP_pT->at(iter) = SF*MuMu->Muon_TuneP_pT->at(iter);
                }

                if( MuMu->isSelPassed == kTRUE )
                {
                    TLorentzVector mu1, mu2;
                    mu1.SetPtEtaPhiE( MuMu->Muon_TuneP_pT->at(0), MuMu->Muon_TuneP_eta->at(0), MuMu->Muon_TuneP_phi->at(0), MuMu->Muon_Energy->at(0) );
                    mu2.SetPtEtaPhiE( MuMu->Muon_TuneP_pT->at(1), MuMu->Muon_TuneP_eta->at(1), MuMu->Muon_TuneP_phi->at(1), MuMu->Muon_Energy->at(1) );
                    // -- Apply efficiency scale factor -- //
                    if( Mgr.isMC == kTRUE )
                    {
                        weight1 = analyzer->EfficiencySF_EventWeight_HLT_BtoF( MuMu );
                        weight2 = analyzer->EfficiencySF_EventWeight_HLT_GtoH( MuMu );
                        effweight = (L_B2F*weight1 + L_G2H*weight2)/L_B2H;
                    }

                    Double_t reco_Pt = ( mu1 + mu2 ).Pt();
                    Double_t reco_rapi = ( mu1 + mu2 ).Rapidity();

                    h_mass_fine_before_PUCorr->Fill( MuMu->Muon_InvM, TotWeight );
                    h_mass_fine_before_EffCorr->Fill( MuMu->Muon_InvM, TotWeight * PUWeight );
                    h_mass_fine->Fill( MuMu->Muon_InvM, TotWeight * PUWeight * effweight );
                    h_mass->Fill( MuMu->Muon_InvM, TotWeight * PUWeight * effweight );
                    h_Pt->Fill( reco_Pt, TotWeight * PUWeight * effweight );
                    h_rapi->Fill( reco_rapi, TotWeight * PUWeight * effweight );

                    h_nPU_beforePUCorr->Fill( MuMu->nPileUp, TotWeight );
                    h_nPU_beforeEffCorr->Fill( MuMu->nPileUp, TotWeight * PUWeight );
                    h_nPU->Fill( MuMu->nPileUp, TotWeight * PUWeight * effweight );

                    h_pT->Fill( MuMu->Muon_TuneP_pT->at(0), TotWeight * PUWeight * effweight );
                    h_pT->Fill( MuMu->Muon_TuneP_pT->at(1), TotWeight * PUWeight * effweight );
                    h_eta->Fill( MuMu->Muon_TuneP_eta->at(0), TotWeight * PUWeight * effweight );
                    h_eta->Fill( MuMu->Muon_TuneP_eta->at(1), TotWeight * PUWeight * effweight );
                    h_phi->Fill( MuMu->Muon_TuneP_phi->at(0), TotWeight * PUWeight * effweight );
                    h_phi->Fill( MuMu->Muon_TuneP_phi->at(1), TotWeight * PUWeight * effweight );

                }// End of event selection
                bar.Draw(i);

            }// End of event iteration

            if( Mgr.isMC == kTRUE ) printf( "\tNormalization factor: %.8f\n", L*Mgr.Xsec[i_tup]/Mgr.Wsum[i_tup] );

            Double_t LoopRunTime = looptime.CpuTime();
            cout << "\tLoop RunTime(" << Mgr.Tag[i_tup] << "): " << LoopRunTime << " seconds\n" << endl;

        }// End of i_tup iteration

        f->cd();
        cout << "\tWriting into file...";

        h_mass_fine_before_PUCorr->Write();
        h_mass_fine_before_RoccoR->Write();
        h_mass_fine_before_EffCorr->Write();
        h_mass_fine->Write();
        h_mass->Write();
        h_Pt->Write();
        h_rapi->Write();

        h_nPU_beforePUCorr->Write();
        h_nPU_beforeEffCorr->Write();
        h_nPU->Write();

        h_pT->Write();
        h_eta->Write();
        h_phi->Write();

        cout << " Finished." << endl;

    } // End of i_proc iteration

    f->Close();
    if ( isBkgFull == kTRUE )
    {
        if ( !f->IsOpen() ) cout << "File Hist_" << Mgr.Procname[_MuMu_Bkg_Full] << ".root has been closed successfully.\n" << endl;
        else cout << "FILE Hist_" << Mgr.Procname[_MuMu_Bkg_Full] << ".root COULD NOT BE CLOSED!\n" << endl;
    }
    else
    {
        if ( !f->IsOpen() ) cout << "File Hist_" << Mgr.Procname[Mgr.CurrentProc] << ".root has been closed successfully.\n" << endl;
        else cout << "FILE Hist_" << Mgr.Procname[Mgr.CurrentProc] << ".root COULD NOT BE CLOSED!\n" << endl;
    }

    Double_t TotalRunTime = totaltime.CpuTime();
    cout << "Total RunTime: " << TotalRunTime << " seconds" << endl;

    TTimeStamp ts_end;
    cout << "[End Time(local time): " << ts_end.AsString("l") << "]" << endl;

} // End of MuMu_HistMaker()


/// --------------------------------- EMu events --------------------------------- ///
void EMu_HistMaker ( TString type, Bool_t SwitchROCCORR, TString HLTname )
{
    if ( !type.Length() )
    {
        cout << "Error: no type specified!" << endl;
        return;
    }

    // -- Run2016 luminosity [/pb] -- //
    Double_t L_B2F = 19721.0, L_G2H = 16146.0, L_B2H = 35867.0, L = 0;
    L = L_B2H;

    TTimeStamp ts_start;
    cout << "[Start Time(local time): " << ts_start.AsString("l") << "]" << endl;
    TStopwatch totaltime;
    totaltime.Start();

    LocalFileMgr Mgr;
    vector<SelProc_t> Processes = Mgr.FindProc( type );
    if ( SwitchROCCORR == kTRUE ) Mgr.SwitchROCCORR(); // If kTRUE, this will go through events that were selected without Rochester correction applied
    TFile *f;
    TString OutputDir;
    Bool_t isBkgFull = kFALSE;  // To tell if this is the _EE_Bkg_Full process (it is handled differently)

    if ( !Processes.size() )
    {
        cout << "Error: no processes!" << endl;
        return;
    }

    if ( Processes[0] == _EMu_Bkg_Full )
    {
        isBkgFull = kTRUE;
        Mgr.GetProc( _EMu_Bkg_Full );
        // -- Output ROOTFile -- //
        OutputDir = Mgr.HistLocation;
        f = new TFile( OutputDir+"Hist_"+Mgr.Procname[_EMu_Bkg_Full]+".root", "RECREATE" );
        Processes.clear();
        Processes.push_back( _EMu_DYTauTau_Full );
        Processes.push_back( _EMu_ttbar_Full );
        Processes.push_back( _EMu_VVnST );
        Processes.push_back( _EMu_WJets );
//        Processes.push_back( _EMu_QCDMuEnriched_Full );
    }

    for ( Int_t i_proc; i_proc<((int)(Processes.size())); i_proc++ )
    {
        if ( Processes[i_proc] <= _EndOf_EE )
        {
            cout << "Error: process " << Mgr.Procname[Processes[i_proc]] << " is not EMu!" << endl;
            continue;
        }

        Mgr.GetProc( Processes[i_proc] );

        // -- Output ROOTFile -- //
        if ( isBkgFull == kFALSE )
        {
            OutputDir = Mgr.HistLocation;
            f = new TFile( OutputDir+"Hist_"+Mgr.Procname[Mgr.CurrentProc]+".root", "RECREATE" );
        }

        cout << "Process: " << Mgr.Procname[Mgr.CurrentProc] << endl;
        cout << "Type: " << Mgr.Type << endl;
        cout << "DATA location: " << Mgr.BaseLocation << endl;
        cout << "Output directory: " << OutputDir << endl;

        TStopwatch totaltime;
        totaltime.Start();

        DYAnalyzer *analyzer = new DYAnalyzer( HLTname );

        // -- For PU re-weighting -- //
        analyzer->SetupPileUpReWeighting_80X( Mgr.isMC, "ROOTFile_PUReWeight_80X_v20170817_64mb.root" );

        // -- For Rochester correction -- //
        TRandom3 *r1 = new TRandom3(0);

        // -- For efficiency SF -- //
        analyzer->SetupEfficiencyScaleFactor_BtoF();
        analyzer->SetupEfficiencyScaleFactor_GtoH();

        // -- Creating Histograms -- //
        TH1D *h_emu_mass_fine_before_PUCorr = new TH1D( "h_emu_mass_fine_before_PUCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 10000, 0, 10000 );
        TH1D *h_emu_mass_fine_before_RoccoR = new TH1D( "h_emu_mass_fine_before_RoccoR_"+Mgr.Procname[Mgr.CurrentProc], "", 10000, 0, 10000 );
        TH1D *h_emu_mass_fine_before_EffCorr = new TH1D( "h_emu_mass_fine_before_EffCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 10000, 0, 10000 );
        TH1D *h_emu_mass_fine = new TH1D( "h_emu_mass_fine_"+Mgr.Procname[Mgr.CurrentProc], "", 10000, 0, 10000 );
        TH1D *h_emu_mass = new TH1D( "h_emu_mass_"+Mgr.Procname[Mgr.CurrentProc], "", 43, massbins );
        TH1D *h_emuSS_mass_fine_before_PUCorr = new TH1D( "h_emuSS_mass_fine_before_PUCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 10000, 0, 10000 );
        TH1D *h_emuSS_mass_fine_before_RoccoR = new TH1D( "h_emuSS_mass_fine_before_RoccoR_"+Mgr.Procname[Mgr.CurrentProc], "", 10000, 0, 10000 );
        TH1D *h_emuSS_mass_fine_before_EffCorr = new TH1D( "h_emuSS_mass_fine_before_EffCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 10000, 0, 10000 );
        TH1D *h_emuSS_mass_fine = new TH1D( "h_emuSS_mass_fine_"+Mgr.Procname[Mgr.CurrentProc], "", 10000, 0, 10000 );
        TH1D *h_emuSS_mass = new TH1D( "h_emuSS_mass_"+Mgr.Procname[Mgr.CurrentProc], "", 43, massbins );

        TH1D *h_nPU_beforePUCorr = new TH1D( "h_nPU_before_PUCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 50, 0, 50 );
        TH1D *h_nPU_beforeEffCorr = new TH1D( "h_nPU_before_EffCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 50, 0, 50 );
        TH1D *h_nPU = new TH1D( "h_nPU_"+Mgr.Procname[Mgr.CurrentProc], "", 50, 0, 50 );

        TH1D *h_ele_pT = new TH1D( "h_ele_pT_"+Mgr.Procname[Mgr.CurrentProc], "", 300, 0, 600 );
        TH1D *h_ele_eta = new TH1D( "h_ele_eta_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -5, 5 );
        TH1D *h_ele_phi = new TH1D( "h_ele_phi_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -5, 5 );
        TH1D *h_eleSS_pT = new TH1D( "h_eleSS_pT_"+Mgr.Procname[Mgr.CurrentProc], "", 300, 0, 600 );
        TH1D *h_eleSS_eta = new TH1D( "h_eleSS_eta_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -5, 5 );
        TH1D *h_eleSS_phi = new TH1D( "h_eleSS_phi_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -5, 5 );

        TH1D *h_mu_pT = new TH1D( "h_mu_pT_"+Mgr.Procname[Mgr.CurrentProc], "", 300, 0, 600 );
        TH1D *h_mu_eta = new TH1D( "h_mu_eta_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -5, 5 );
        TH1D *h_mu_phi = new TH1D( "h_mu_phi_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -5, 5 );
        TH1D *h_muSS_pT = new TH1D( "h_muSS_pT_"+Mgr.Procname[Mgr.CurrentProc], "", 300, 0, 600 );
        TH1D *h_muSS_eta = new TH1D( "h_muSS_eta_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -5, 5 );
        TH1D *h_muSS_phi = new TH1D( "h_muSS_phi_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -5, 5 );

        //Loop for all samples
        const Int_t Ntup = Mgr.FileLocation.size();
        for(Int_t i_tup = 0; i_tup<Ntup; i_tup++)
        {
            TStopwatch looptime;
            looptime.Start();

            cout << "\t<" << Mgr.Tag[i_tup] << ">" << endl;

            TChain *chain = new TChain( Mgr.TreeName[i_tup] );
            chain->Add( Mgr.FullLocation[i_tup] );

            SelectedEMu_t *EMu = new SelectedEMu_t();
            EMu->CreateFromChain( chain );

            RoccoR rc("./etc/RoccoR/rcdata.2016.v3");           

            Int_t NEvents = chain->GetEntries();
            if ( NEvents != Mgr.nEvents[i_tup] ) cout << "\tEvent numbers do not match!!!" << endl;
            cout << "\t[Total events: " << Mgr.Wsum[i_tup] << "]" << endl;
            cout << "\t[Selected Events: " << NEvents << "]" << endl;

            myProgressBar_t bar( NEvents );

            for(Int_t i=0; i<NEvents; i++)
            {
                EMu->GetEvent(i);

                Double_t GenWeight = EMu->GENEvt_weight;

                // -- Pileup-Reweighting -- //
                Double_t PUWeight = 1;
                if( Mgr.isMC == kTRUE ) PUWeight = analyzer->PileUpWeightValue_80X( EMu->nPileUp );

                // -- efficiency weights -- //
                Double_t weight1 = 0, weight2 = 0, effweight = 1;

                // -- Normalization -- //
                Double_t TotWeight = GenWeight;
                if( Mgr.isMC == kTRUE ) TotWeight = ( L * Mgr.Xsec[i_tup] / Mgr.Wsum[i_tup] ) * GenWeight;

                // -- Rochester correction -- //
                Double_t rndm[2], SF=0; r1->RndmArray(2, rndm);
                Int_t s, m;

                if( Mgr.Type == "DATA" )
                        SF = rc.kScaleDT(EMu->Muon_charge, EMu->Muon_TuneP_pT, EMu->Muon_TuneP_eta,
                                         EMu->Muon_TuneP_phi, s=0, m=0);
                else
                        SF = rc.kScaleAndSmearMC(EMu->Muon_charge, EMu->Muon_TuneP_pT,
                                                 EMu->Muon_TuneP_eta, EMu->Muon_TuneP_phi,
                                                 EMu->Muon_trackerLayers, rndm[0], rndm[1], s=0, m=0);

                EMu->Muon_TuneP_pT = SF*EMu->Muon_TuneP_pT;

                if( EMu->isSelPassed == kTRUE )
                {
                    TLorentzVector mu, ele;
                    mu.SetPtEtaPhiE( EMu->Muon_TuneP_pT, EMu->Muon_TuneP_eta, EMu->Muon_TuneP_phi, EMu->Muon_Energy );
                    ele.SetPtEtaPhiE( EMu->Electron_pT, EMu->Electron_eta, EMu->Electron_phi, EMu->Electron_Energy );
                    // -- Apply efficiency scale factor -- //
                    if( Mgr.isMC == kTRUE )
                    {
                        weight1 = analyzer->EfficiencySF_EventWeight_HLT_BtoF( EMu );
                        weight2 = analyzer->EfficiencySF_EventWeight_HLT_GtoH( EMu );
                        effweight = (L_B2F*weight1 + L_G2H*weight2)/L_B2H;
                    }

                    Double_t reco_Pt = ( mu + ele ).Pt();
                    Double_t reco_rapi = ( mu + ele ).Rapidity();

                    if ( EMu->Electron_charge != EMu->Muon_charge )
                    {
                        h_emu_mass_fine_before_PUCorr->Fill( EMu->EMu_InvM, TotWeight );
                        h_emu_mass_fine_before_EffCorr->Fill( EMu->EMu_InvM, TotWeight * PUWeight );
                        h_emu_mass_fine->Fill( EMu->EMu_InvM, TotWeight * PUWeight * effweight );
                        h_emu_mass->Fill( EMu->EMu_InvM, TotWeight * PUWeight * effweight );

                        h_ele_pT->Fill( EMu->Electron_pT, TotWeight * PUWeight * effweight );
                        h_ele_eta->Fill( EMu->Electron_eta, TotWeight * PUWeight * effweight );
                        h_ele_phi->Fill( EMu->Electron_phi, TotWeight * PUWeight * effweight );

                        h_mu_pT->Fill( EMu->Muon_TuneP_pT, TotWeight * PUWeight * effweight );
                        h_mu_eta->Fill( EMu->Muon_TuneP_eta, TotWeight * PUWeight * effweight );
                        h_mu_phi->Fill( EMu->Muon_TuneP_phi, TotWeight * PUWeight * effweight );
                    }
                    else
                    {
                        h_emuSS_mass_fine_before_PUCorr->Fill( EMu->EMu_InvM, TotWeight );
                        h_emuSS_mass_fine_before_EffCorr->Fill( EMu->EMu_InvM, TotWeight * PUWeight );
                        h_emuSS_mass_fine->Fill( EMu->EMu_InvM, TotWeight * PUWeight * effweight );
                        h_emuSS_mass->Fill( EMu->EMu_InvM, TotWeight * PUWeight * effweight );

                        h_eleSS_pT->Fill( EMu->Electron_pT, TotWeight * PUWeight * effweight );
                        h_eleSS_eta->Fill( EMu->Electron_eta, TotWeight * PUWeight * effweight );
                        h_eleSS_phi->Fill( EMu->Electron_phi, TotWeight * PUWeight * effweight );

                        h_muSS_pT->Fill( EMu->Muon_TuneP_pT, TotWeight * PUWeight * effweight );
                        h_muSS_eta->Fill( EMu->Muon_TuneP_eta, TotWeight * PUWeight * effweight );
                        h_muSS_phi->Fill( EMu->Muon_TuneP_phi, TotWeight * PUWeight * effweight );
                    }
                    h_nPU_beforePUCorr->Fill( EMu->nPileUp, TotWeight );
                    h_nPU_beforeEffCorr->Fill( EMu->nPileUp, TotWeight * PUWeight );
                    h_nPU->Fill( EMu->nPileUp, TotWeight * PUWeight * effweight );

                }// End of event selection
                bar.Draw(i);

            }// End of event iteration

            if( Mgr.isMC == kTRUE ) printf( "\tNormalization factor: %.8f\n", L*Mgr.Xsec[i_tup]/Mgr.Wsum[i_tup] );

            Double_t LoopRunTime = looptime.CpuTime();
            cout << "\tLoop RunTime(" << Mgr.Tag[i_tup] << "): " << LoopRunTime << " seconds\n" << endl;

        }// End of i_tup iteration

        f->cd();
        cout << "\tWriting into file...";

        h_emu_mass_fine_before_PUCorr->Write();
        h_emu_mass_fine_before_RoccoR->Write();
        h_emu_mass_fine_before_EffCorr->Write();
        h_emu_mass_fine->Write();
        h_emu_mass->Write();
        h_emuSS_mass_fine_before_PUCorr->Write();
        h_emuSS_mass_fine_before_RoccoR->Write();
        h_emuSS_mass_fine_before_EffCorr->Write();
        h_emuSS_mass_fine->Write();
        h_emuSS_mass->Write();

        h_nPU_beforePUCorr->Write();
        h_nPU_beforeEffCorr->Write();
        h_nPU->Write();

        h_ele_pT->Write();
        h_ele_eta->Write();
        h_ele_phi->Write();
        h_eleSS_pT->Write();
        h_eleSS_eta->Write();
        h_eleSS_phi->Write();

        h_mu_pT->Write();
        h_mu_eta->Write();
        h_mu_phi->Write();
        h_muSS_pT->Write();
        h_muSS_eta->Write();
        h_muSS_phi->Write();

        cout << " Finished." << endl;       

    }// End of i_proc iteration

    f->Close();
    if ( isBkgFull == kTRUE )
    {
        if ( !f->IsOpen() ) cout << "File Hist_" << Mgr.Procname[_EMu_Bkg_Full] << ".root has been closed successfully.\n" << endl;
        else cout << "FILE Hist_" << Mgr.Procname[_EMu_Bkg_Full] << ".root COULD NOT BE CLOSED!\n" << endl;
    }
    else
    {
        if ( !f->IsOpen() ) cout << "File Hist_" << Mgr.Procname[Mgr.CurrentProc] << ".root has been closed successfully.\n" << endl;
        else cout << "FILE Hist_" << Mgr.Procname[Mgr.CurrentProc] << ".root COULD NOT BE CLOSED!\n" << endl;
    }

    Double_t TotalRunTime = totaltime.CpuTime();
    cout << "Total RunTime: " << TotalRunTime << " seconds" << endl;

    TTimeStamp ts_end;
    cout << "[End Time(local time): " << ts_end.AsString("l") << "]" << endl;

} // End of EMu_HistMaker()


/// ----------------------- MuMu_merge ----------------------------- ///
void MuMu_merge()
{
    LocalFileMgr Mgr;

    TH1D *h_mass_fine_before_PUCorr[7], *h_mass_fine_before_RoccoR[7], *h_mass_fine_before_EffCorr[7], *h_mass_fine[7], *h_mass[7], *h_Pt[7], *h_rapi[7],
         *h_nPU_before_PUCorr[7], *h_nPU_before_EffCorr[7], *h_nPU[7], *h_pT[7], *h_eta[7], *h_phi[7];

    TFile* files[7];
    Mgr.GetProc(_MuMu_SingleMuon_Full);

    TH1D *h_mass_fine_before_PUCorr_full  = new TH1D( "h_mass_fine_before_PUCorr_"+Mgr.Procname[Mgr.CurrentProc],
                                                      "h_mass_fine_before_PUCorr_"+Mgr.Procname[Mgr.CurrentProc], 10000, 0, 10000 );
    TH1D *h_mass_fine_before_RoccoR_full  = new TH1D( "h_mass_fine_before_RoccoR_"+Mgr.Procname[Mgr.CurrentProc],
                                                      "h_mass_fine_before_RoccoR_"+Mgr.Procname[Mgr.CurrentProc], 10000, 0, 10000 );
    TH1D *h_mass_fine_before_EffCorr_full = new TH1D( "h_mass_fine_before_EffCorr_"+Mgr.Procname[Mgr.CurrentProc],
                                                      "h_mass_fine_before_EffCorr_"+Mgr.Procname[Mgr.CurrentProc], 10000, 0, 10000 );
    TH1D *h_mass_fine_full                = new TH1D( "h_mass_fine_"+Mgr.Procname[Mgr.CurrentProc],
                                                      "h_mass_fine_"+Mgr.Procname[Mgr.CurrentProc], 10000, 0, 10000 );
    TH1D *h_mass_full                     = new TH1D( "h_mass_"+Mgr.Procname[Mgr.CurrentProc],
                                                      "h_mass_"+Mgr.Procname[Mgr.CurrentProc], 43, massbins );
    TH1D *h_Pt_full                       = new TH1D( "h_Pt_"+Mgr.Procname[Mgr.CurrentProc],
                                                      "h_Pt_"+Mgr.Procname[Mgr.CurrentProc], 300, 0, 600 );
    TH1D *h_rapi_full                     = new TH1D( "h_rapi_"+Mgr.Procname[Mgr.CurrentProc],
                                                      "h_rapi_"+Mgr.Procname[Mgr.CurrentProc], 100, -5, 5 );
    TH1D *h_nPU_before_PUCorr_full         = new TH1D( "h_nPU_before_PUCorr_"+Mgr.Procname[Mgr.CurrentProc],
                                                      "h_nPU_before_PUCorr_"+Mgr.Procname[Mgr.CurrentProc], 50, 0, 50 );
    TH1D *h_nPU_before_EffCorr_full        = new TH1D( "h_nPU_before_EffCorr_"+Mgr.Procname[Mgr.CurrentProc],
                                                      "h_nPU_before_EffCorr_"+Mgr.Procname[Mgr.CurrentProc], 50, 0, 50 );
    TH1D *h_nPU_full                      = new TH1D( "h_nPU_"+Mgr.Procname[Mgr.CurrentProc],
                                                      "h_nPU_"+Mgr.Procname[Mgr.CurrentProc], 50, 0, 50 );
    TH1D *h_pT_full                       = new TH1D( "h_pT_"+Mgr.Procname[Mgr.CurrentProc],
                                                      "h_pT_"+Mgr.Procname[Mgr.CurrentProc], 300, 0, 600 );
    TH1D *h_eta_full                      = new TH1D( "h_eta_"+Mgr.Procname[Mgr.CurrentProc],
                                                      "h_eta_"+Mgr.Procname[Mgr.CurrentProc], 100, -5, 5 );
    TH1D *h_phi_full                      = new TH1D( "h_phi_"+Mgr.Procname[Mgr.CurrentProc],
                                                      "h_phi_"+Mgr.Procname[Mgr.CurrentProc], 100, -5, 5 );

    TFile* newFile = new TFile( Mgr.HistLocation+"Hist_"+Mgr.Procname[Mgr.CurrentProc]+".root", "RECREATE" );
    Int_t iter = 0;
    for ( SelProc_t pr=_MuMu_SingleMuon_B; pr<=_MuMu_SingleMuon_H; pr=next(pr) )
    {
        Mgr.GetProc(pr);
        files[iter] = new TFile( Mgr.HistLocation+"Hist_"+Mgr.Procname[Mgr.CurrentProc]+".root", "READ" );

        files[iter]->GetObject( "h_mass_fine_before_PUCorr_"+Mgr.Procname[Mgr.CurrentProc],  h_mass_fine_before_PUCorr[iter] );
        files[iter]->GetObject( "h_mass_fine_before_RoccoR_"+Mgr.Procname[Mgr.CurrentProc],  h_mass_fine_before_RoccoR[iter] );
        files[iter]->GetObject( "h_mass_fine_before_EffCorr_"+Mgr.Procname[Mgr.CurrentProc], h_mass_fine_before_EffCorr[iter] );
        files[iter]->GetObject( "h_mass_fine_"+Mgr.Procname[Mgr.CurrentProc],                h_mass_fine[iter] );
        files[iter]->GetObject( "h_mass_"+Mgr.Procname[Mgr.CurrentProc],                     h_mass[iter] );
        files[iter]->GetObject( "h_Pt_"+Mgr.Procname[Mgr.CurrentProc],                       h_Pt[iter] );
        files[iter]->GetObject( "h_rapi_"+Mgr.Procname[Mgr.CurrentProc],                     h_rapi[iter] );
        files[iter]->GetObject( "h_nPU_before_PUCorr_"+Mgr.Procname[Mgr.CurrentProc],        h_nPU_before_PUCorr[iter] );
        files[iter]->GetObject( "h_nPU_before_EffCorr_"+Mgr.Procname[Mgr.CurrentProc],       h_nPU_before_EffCorr[iter] );
        files[iter]->GetObject( "h_nPU_"+Mgr.Procname[Mgr.CurrentProc],                      h_nPU[iter] );
        files[iter]->GetObject( "h_pT_"+Mgr.Procname[Mgr.CurrentProc],                       h_pT[iter] );
        files[iter]->GetObject( "h_eta_"+Mgr.Procname[Mgr.CurrentProc],                      h_eta[iter] );
        files[iter]->GetObject( "h_phi_"+Mgr.Procname[Mgr.CurrentProc],                      h_phi[iter] );

        h_mass_fine_before_PUCorr[iter] ->SetDirectory(0);
        h_mass_fine_before_RoccoR[iter] ->SetDirectory(0);
        h_mass_fine_before_EffCorr[iter]->SetDirectory(0);
        h_mass_fine[iter]               ->SetDirectory(0);
        h_mass[iter]                    ->SetDirectory(0);
        h_Pt[iter]                      ->SetDirectory(0);
        h_rapi[iter]                    ->SetDirectory(0);
        h_nPU_before_PUCorr[iter]       ->SetDirectory(0);
        h_nPU_before_EffCorr[iter]      ->SetDirectory(0);
        h_nPU[iter]                     ->SetDirectory(0);
        h_pT[iter]                      ->SetDirectory(0);
        h_eta[iter]                     ->SetDirectory(0);
        h_phi[iter]                     ->SetDirectory(0);

        h_mass_fine_before_PUCorr_full ->Add( h_mass_fine_before_PUCorr[iter] );
        h_mass_fine_before_RoccoR_full ->Add( h_mass_fine_before_RoccoR[iter] );
        h_mass_fine_before_EffCorr_full->Add( h_mass_fine_before_EffCorr[iter] );
        h_mass_fine_full               ->Add( h_mass_fine[iter] );
        h_mass_full                    ->Add( h_mass[iter] );
        h_Pt_full                      ->Add( h_Pt[iter] );
        h_rapi_full                    ->Add( h_rapi[iter] );
        h_nPU_before_PUCorr_full       ->Add( h_nPU_before_PUCorr[iter] );
        h_nPU_before_EffCorr_full      ->Add( h_nPU_before_EffCorr[iter] );
        h_nPU_full                     ->Add( h_nPU[iter] );
        h_pT_full                      ->Add( h_pT[iter] );
        h_eta_full                     ->Add( h_eta[iter] );
        h_phi_full                     ->Add( h_phi[iter] );

        iter++;
    }// End of for(SingleMuon_BtoH)

    newFile->cd();
    h_mass_fine_before_PUCorr_full  ->Write();
    h_mass_fine_before_RoccoR_full  ->Write();
    h_mass_fine_before_EffCorr_full ->Write();
    h_mass_fine_full                ->Write();
    h_mass_full                     ->Write();
    h_Pt_full                       ->Write();
    h_rapi_full                     ->Write();
    h_nPU_before_PUCorr_full        ->Write();
    h_nPU_before_EffCorr_full       ->Write();
    h_nPU_full                      ->Write();
    h_pT_full                       ->Write();
    h_eta_full                      ->Write();
    h_phi_full                      ->Write();

}// End of MuMu_Merge
