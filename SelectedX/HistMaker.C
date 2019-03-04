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

void EE_HistMaker ( TString type, TString HLTname, Bool_t DEBUG );
void MuMu_HistMaker ( TString type, Bool_t SwitchROCCORR, TString HLTname, Bool_t DEBUG );
void EMu_HistMaker ( TString type, Bool_t SwitchROCCORR, TString HLTname, Bool_t DEBUG );
void MuMu_merge();

// -- Drell-Yan mass bins -- //
const Double_t massbins[44] = {15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 64, 68, 72, 76, 81, 86, 91, 96, 101, 106, 110, 115, 120, 126,
                               133, 141, 150, 160, 171, 185, 200, 220, 243, 273, 320, 380, 440, 510, 600, 700, 830, 1000, 1500, 3000};


void HistMaker ( TString WhichX = "", TString Type = "", Bool_t SwitchROCCORR = kFALSE, TString HLTname = "DEFAULT" )
{
    TString whichX = WhichX;
    whichX.ToUpper();
    TString HLT;
    Int_t Xselected = 0;
    Bool_t DEBUG = kFALSE;
    if ( whichX.Contains("DEBUG") )
    {
        DEBUG = kTRUE;
        cout << "**** DEBUG MODE: Running with 10 events only. ****" << endl;
    }
    TString type = whichX+"_"+Type;
    if ( whichX.Contains("EE") )
    {
        Xselected++;
        if ( HLTname == "DEFAULT" ) HLT = "Ele23Ele12";
        else HLT = HLTname;
        cout << "\n*******      EE_HistMaker ( " << type << " )      *******" << endl;
        EE_HistMaker( type, HLT, DEBUG );
    }
    if ( whichX.Contains("MUMU") )
    {
        Xselected++;
        if ( HLTname == "DEFAULT" ) HLT = "IsoMu24_OR_IsoTkMu24";
        else HLT = HLTname;
        cout << "\n*****  MuMu_HistMaker ( " << type << " )  *****" << endl;
        MuMu_HistMaker( type, SwitchROCCORR, HLT, DEBUG );
    }
    if ( whichX.Contains("EMU") )
    {
        Xselected++;
        if ( HLTname == "DEFAULT" ) HLT = "IsoMu24_OR_IsoTkMu24";
        else HLT = HLTname;
        cout << "\n*****   EMu_HistMaker ( " << type << " )  *****" << endl;
        EMu_HistMaker( type, SwitchROCCORR, HLT, DEBUG );
    }
    if ( whichX.Contains("MUMU_MERGE") )
    {
        Xselected++;
        cout << "\n*********   MuMu_merge  *********" << endl;
        MuMu_merge();
    }
    if ( Xselected == 0 ) cout << "Wrong arument! \nType in: >> .x HistMaker.C+(\"whichX\", \"whichProcess\", SwitchROCCORR)" << endl;

} // End of HistMaker()


/// ----------------------------- Electron Channel ------------------------------ ///
void EE_HistMaker (TString type, TString HLTname , Bool_t DEBUG)
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
    TString debug = "";
    if ( DEBUG == kTRUE ) debug = "_DEBUG";
    Bool_t isBkgFull = kFALSE;  // To tell if this is the _EE_Bkg_Full process (it is handled differently)

    if ( !Processes.size() )
    {
        cout << "Error: no processes!" << endl;
        return;
    }

    if ( Processes[0] == _EE_Bkg_Full )
    {
        isBkgFull = kTRUE;
        Mgr.SetProc( _EE_Bkg_Full );
        // -- Output ROOTFile -- //
        OutputDir = Mgr.HistLocation;
        f = new TFile( OutputDir+"Hist_"+Mgr.Procname[_EE_Bkg_Full]+debug+".root", "RECREATE" );
        Processes.clear();
        Processes.push_back( _EE_DYTauTau_Full );
        Processes.push_back( _EE_ttbar_Full );
        Processes.push_back( _EE_tW );
        Processes.push_back( _EE_tbarW );
        Processes.push_back( _EE_WW );
        Processes.push_back( _EE_WZ );
        Processes.push_back( _EE_ZZ );
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

        Mgr.SetProc( Processes[i_proc] );

        // -- Output ROOTFile -- //
        if ( isBkgFull == kFALSE )
        {
            OutputDir = Mgr.HistLocation;
            f = new TFile( OutputDir+"Hist_"+Mgr.Procname[Mgr.CurrentProc]+debug+".root", "RECREATE" );
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

        // -- For PVz reweighting -- //
//        analyzer->SetupPVzWeights(Mgr.isMC, "ee", "./etc/PVzWeights.root");
        analyzer->SetupPVzWeights(Mgr.isMC, "mumu", "./etc/PVzWeights.root");


        // -- Creating Histograms -- //
        TH1D *h_mass_before_PUCorr = new TH1D("h_mass_before_PUCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 43, massbins);
        TH1D *h_mass_before_EffCorr = new TH1D("h_mass_before_EffCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 43, massbins);
        TH1D *h_mass_before_PVzCorr = new TH1D("h_mass_before_PVzCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 43, massbins);
        TH1D *h_mass_before_L1Corr = new TH1D("h_mass_before_L1Corr_"+Mgr.Procname[Mgr.CurrentProc], "", 43, massbins);
        TH1D *h_mass_fine = new TH1D("h_mass_fine_"+Mgr.Procname[Mgr.CurrentProc], "", 10000, 0, 10000);
        TH1D *h_mass = new TH1D("h_mass_"+Mgr.Procname[Mgr.CurrentProc], "", 43, massbins);
        TH1D *h_pT = new TH1D("h_pT_"+Mgr.Procname[Mgr.CurrentProc], "", 100, 0, 1000);
        TH1D *h_rapi = new TH1D("h_rapi_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -5, 5);

        TH1D *h_nVTX_before_PUCorr = new TH1D( "h_nVTX_before_PUCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 75, 0, 75 );
        TH1D *h_nVTX_before_EffCorr = new TH1D( "h_nVTX_before_EffCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 75, 0, 75 );
        TH1D *h_nVTX = new TH1D( "h_nVTX_"+Mgr.Procname[Mgr.CurrentProc], "", 75, 0, 75 );

        TH1D *h_pT_lead_before_PUCorr = new TH1D("h_pT_lead_before_PUCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, 0, 1000);
        TH1D *h_pT_sublead_before_PUCorr = new TH1D("h_pT_sublead_before_PUCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, 0, 1000);
        TH1D *h_eta_lead_before_PUCorr = new TH1D("h_eta_lead_before_PUCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4);
        TH1D *h_eta_sublead_before_PUCorr = new TH1D("h_eta_sublead_before_PUCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4);
        TH1D *h_phi_lead_before_PUCorr = new TH1D("h_phi_lead_before_PUCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4);
        TH1D *h_phi_sublead_before_PUCorr = new TH1D("h_phi_sublead_before_PUCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4);

        TH1D *h_pT_lead_before_EffCorr = new TH1D("h_pT_lead_before_EffCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, 0, 1000);
        TH1D *h_pT_sublead_before_EffCorr = new TH1D("h_pT_sublead_before_EffCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, 0, 1000);
        TH1D *h_eta_lead_before_EffCorr = new TH1D("h_eta_lead_before_EffCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4);
        TH1D *h_eta_sublead_before_EffCorr = new TH1D("h_eta_sublead_before_EffCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4);
        TH1D *h_phi_lead_before_EffCorr = new TH1D("h_phi_lead_before_EffCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4);
        TH1D *h_phi_sublead_before_EffCorr = new TH1D("h_phi_sublead_before_EffCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4);

        TH1D *h_pT_lead_before_PVzCorr = new TH1D("h_pT_lead_before_PVzCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, 0, 1000);
        TH1D *h_pT_sublead_before_PVzCorr = new TH1D("h_pT_sublead_before_PVzCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, 0, 1000);
        TH1D *h_eta_lead_before_PVzCorr = new TH1D("h_eta_lead_before_PVzCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4);
        TH1D *h_eta_sublead_before_PVzCorr = new TH1D("h_eta_sublead_before_PVzCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4);
        TH1D *h_phi_lead_before_PVzCorr = new TH1D("h_phi_lead_before_PVzCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4);
        TH1D *h_phi_sublead_before_PVzCorr = new TH1D("h_phi_sublead_before_PVzCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4);

        TH1D *h_pT_lead_before_L1Corr = new TH1D("h_pT_lead_before_L1Corr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, 0, 1000);
        TH1D *h_pT_sublead_before_L1Corr = new TH1D("h_pT_sublead_before_L1Corr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, 0, 1000);
        TH1D *h_eta_lead_before_L1Corr = new TH1D("h_eta_lead_before_L1Corr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4);
        TH1D *h_eta_sublead_before_L1Corr = new TH1D("h_eta_sublead_before_L1Corr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4);
        TH1D *h_phi_lead_before_L1Corr = new TH1D("h_phi_lead_before_L1Corr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4);
        TH1D *h_phi_sublead_before_L1Corr = new TH1D("h_phi_sublead_before_L1Corr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4);

        TH1D *h_pT_lead = new TH1D("h_pT_lead_"+Mgr.Procname[Mgr.CurrentProc], "", 100, 0, 1000);
        TH1D *h_pT_sublead = new TH1D("h_pT_sublead_"+Mgr.Procname[Mgr.CurrentProc], "", 100, 0, 1000);
        TH1D *h_eta_lead = new TH1D("h_eta_lead_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4);
        TH1D *h_eta_sublead = new TH1D("h_eta_sublead_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4);
        TH1D *h_phi_lead = new TH1D("h_phi_lead_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4);
        TH1D *h_phi_sublead = new TH1D("h_phi_sublead_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4);

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
            cout << "\t[Sum of weights:: " << Mgr.Wsum[i_tup] << "]" << endl;
            cout << "\t[Selected Events: " << NEvents << "]" << endl;

            if ( DEBUG == kTRUE ) NEvents = 10;
            myProgressBar_t bar( NEvents );

            for ( Int_t i=0; i<NEvents; i++ )
            {
                EE->GetEvent(i);

                // -- Bring the weights -- //
                Double_t GenWeight = EE->GENEvt_weight;
                if ( GenWeight > 1 ) GenWeight = 1;
                else if ( GenWeight < -1 ) GenWeight = -1;


                // -- Pileup-Reweighting -- //
                Double_t PUWeight = 1;
                if( Mgr.isMC == kTRUE ) PUWeight = analyzer->PileUpWeightValue_80X( EE->nPileUp );

                // -- efficiency weights -- //
                Double_t effweight = 1;

                // -- PVz weights -- //
                Double_t PVzWeight = 1;
                if( Mgr.isMC == kTRUE ) PVzWeight = analyzer->PVzWeightValue( EE->PVz );

                // -- L1 prefiring weights -- //
                Double_t L1weight = 1;
                if ( Mgr.isMC == kTRUE ) L1weight = EE->_prefiringweight;
//                if ( Mgr.isMC == kTRUE ) L1weight = EE->_prefiringweightup;
//                if ( Mgr.isMC == kTRUE ) L1weight = EE->_prefiringweightdown;

                // -- Normalization -- //
                Double_t TotWeight = GenWeight;
                if ( Mgr.isMC == kTRUE ) TotWeight = ( L * Mgr.Xsec[i_tup] / Mgr.Wsum[i_tup] ) * GenWeight;

                if( EE->isSelPassed == kTRUE )
                {
                    TLorentzVector ele1, ele2;
                    ele1.SetPtEtaPhiE( EE->Electron_pT->at(0), EE->Electron_eta->at(0), EE->Electron_phi->at(0), EE->Electron_Energy->at(0) );
                    ele2.SetPtEtaPhiE( EE->Electron_pT->at(1), EE->Electron_eta->at(1), EE->Electron_phi->at(1), EE->Electron_Energy->at(1) );
                    Double_t reco_Pt = ( ele1 + ele2 ).Pt();
                    Double_t reco_rapi = ( ele1 + ele2 ).Rapidity();

                    // -- Apply efficiency correcion -- //
                    if( Mgr.isMC == kTRUE )
                            effweight = analyzer->EfficiencySF_EventWeight_electron( EE );
//                    cout << effweight << endl;

                    h_mass_before_PUCorr->Fill( EE->Electron_InvM, TotWeight );
                    h_mass_before_EffCorr->Fill( EE->Electron_InvM, TotWeight * PUWeight );
                    h_mass_before_PVzCorr->Fill( EE->Electron_InvM, TotWeight * PUWeight * effweight );
                    h_mass_before_L1Corr->Fill( EE->Electron_InvM, TotWeight * PUWeight * effweight * PVzWeight );
                    h_mass_fine->Fill( EE->Electron_InvM, TotWeight * PUWeight * effweight * PVzWeight * L1weight );
                    h_mass->Fill( EE->Electron_InvM, TotWeight * PUWeight * effweight * PVzWeight * L1weight );
                    h_pT->Fill( reco_Pt, TotWeight * PUWeight * effweight * PVzWeight * L1weight );
                    h_rapi->Fill( reco_rapi, TotWeight * PUWeight * effweight * PVzWeight * L1weight );

                    h_nVTX_before_PUCorr->Fill( EE->nVertices, TotWeight );
                    h_nVTX_before_EffCorr->Fill( EE->nVertices, TotWeight * PUWeight );
                    h_nVTX->Fill( EE->nVertices, TotWeight * PUWeight * effweight * PVzWeight * L1weight );

                    int lead=0, sublead=1;
                    if ( EE->Electron_pT->at(0) < EE->Electron_pT->at(1) )
                    {
                        lead = 1;
                        sublead = 0;
                    }

                    h_pT_lead_before_PUCorr->Fill( EE->Electron_pT->at(lead), TotWeight );
                    h_pT_sublead_before_PUCorr->Fill( EE->Electron_pT->at(sublead), TotWeight );
//                    h_eta_lead_before_PUCorr->Fill( EE->Electron_eta->at(lead), TotWeight );
//                    h_eta_sublead_before_PUCorr->Fill( EE->Electron_eta->at(sublead), TotWeight );
                    h_eta_lead_before_PUCorr->Fill( EE->Electron_etaSC->at(lead), TotWeight );
                    h_eta_sublead_before_PUCorr->Fill( EE->Electron_etaSC->at(sublead), TotWeight );
                    h_phi_lead_before_PUCorr->Fill( EE->Electron_phi->at(lead), TotWeight );
                    h_phi_sublead_before_PUCorr->Fill( EE->Electron_phi->at(sublead), TotWeight );

                    h_pT_lead_before_EffCorr->Fill( EE->Electron_pT->at(lead), TotWeight * PUWeight );
                    h_pT_sublead_before_EffCorr->Fill( EE->Electron_pT->at(sublead), TotWeight * PUWeight );
//                    h_eta_lead_before_EffCorr->Fill( EE->Electron_eta->at(lead), TotWeight * PUWeight );
//                    h_eta_sublead_before_EffCorr->Fill( EE->Electron_eta->at(sublead), TotWeight * PUWeight );
                    h_eta_lead_before_EffCorr->Fill( EE->Electron_etaSC->at(lead), TotWeight * PUWeight );
                    h_eta_sublead_before_EffCorr->Fill( EE->Electron_etaSC->at(sublead), TotWeight * PUWeight );
                    h_phi_lead_before_EffCorr->Fill( EE->Electron_phi->at(lead), TotWeight * PUWeight );
                    h_phi_sublead_before_EffCorr->Fill( EE->Electron_phi->at(sublead), TotWeight * PUWeight );

                    h_pT_lead_before_PVzCorr->Fill( EE->Electron_pT->at(lead), TotWeight * PUWeight * effweight );
                    h_pT_sublead_before_PVzCorr->Fill( EE->Electron_pT->at(sublead), TotWeight * PUWeight * effweight );
//                    h_eta_lead_before_PVzCorr->Fill( EE->Electron_eta->at(lead), TotWeight * PUWeight * effweight );
//                    h_eta_sublead_before_PVzCorr->Fill( EE->Electron_eta->at(sublead), TotWeight * PUWeight * effweight );
                    h_eta_lead_before_PVzCorr->Fill( EE->Electron_etaSC->at(lead), TotWeight * PUWeight * effweight );
                    h_eta_sublead_before_PVzCorr->Fill( EE->Electron_etaSC->at(sublead), TotWeight * PUWeight * effweight );
                    h_phi_lead_before_PVzCorr->Fill( EE->Electron_phi->at(lead), TotWeight * PUWeight * effweight );
                    h_phi_sublead_before_PVzCorr->Fill( EE->Electron_phi->at(sublead), TotWeight * PUWeight * effweight );

                    h_pT_lead_before_L1Corr->Fill( EE->Electron_pT->at(lead), TotWeight * PUWeight * effweight * PVzWeight );
                    h_pT_sublead_before_L1Corr->Fill( EE->Electron_pT->at(sublead), TotWeight * PUWeight * effweight * PVzWeight );
//                    h_eta_lead_before_L1Corr->Fill( EE->Electron_eta->at(lead), TotWeight * PUWeight * effweight * PVzWeight );
//                    h_eta_sublead_before_L1Corr->Fill( EE->Electron_eta->at(sublead), TotWeight * PUWeight * effweight * PVzWeight );
                    h_eta_lead_before_L1Corr->Fill( EE->Electron_etaSC->at(lead), TotWeight * PUWeight * effweight * PVzWeight );
                    h_eta_sublead_before_L1Corr->Fill( EE->Electron_etaSC->at(sublead), TotWeight * PUWeight * effweight * PVzWeight );
                    h_phi_lead_before_L1Corr->Fill( EE->Electron_phi->at(lead), TotWeight * PUWeight * effweight * PVzWeight );
                    h_phi_sublead_before_L1Corr->Fill( EE->Electron_phi->at(sublead), TotWeight * PUWeight * effweight * PVzWeight );

                    h_pT_lead->Fill( EE->Electron_pT->at(lead), TotWeight * PUWeight * effweight * PVzWeight * L1weight );
                    h_pT_sublead->Fill( EE->Electron_pT->at(sublead), TotWeight * PUWeight * effweight * PVzWeight * L1weight );
//                    h_eta_lead->Fill( EE->Electron_eta->at(lead), TotWeight * PUWeight * effweight * PVzWeight * L1weight );
//                    h_eta_sublead->Fill( EE->Electron_eta->at(sublead), TotWeight * PUWeight * effweight * PVzWeight * L1weight );
                    h_eta_lead->Fill( EE->Electron_etaSC->at(lead), TotWeight * PUWeight * effweight * PVzWeight * L1weight );
                    h_eta_sublead->Fill( EE->Electron_etaSC->at(sublead), TotWeight * PUWeight * effweight * PVzWeight * L1weight );
                    h_phi_lead->Fill( EE->Electron_phi->at(lead), TotWeight * PUWeight * effweight * PVzWeight * L1weight );
                    h_phi_sublead->Fill( EE->Electron_phi->at(sublead), TotWeight * PUWeight * effweight * PVzWeight * L1weight );

                } // End of event selection
                bar.Draw(i);

            } // End of event iteration

            if( Mgr.isMC == kTRUE ) printf( "\tNormalization factor: %.8f\n", L*Mgr.Xsec[i_tup]/Mgr.Wsum[i_tup] );

            Double_t LoopRunTime = looptime.CpuTime();
            cout << "\tLoop RunTime(" << Mgr.Tag[i_tup] << "): " << LoopRunTime << " seconds\n" << endl;

        }// End of i_tup iteration

        f->cd();
        cout << "\tWriting into file...";

        h_mass_before_PUCorr->Write();
        h_mass_before_EffCorr->Write();
        h_mass_before_PVzCorr->Write();
        h_mass_before_L1Corr->Write();
        h_mass_fine->Write();
        h_mass->Write();
        h_pT->Write();
        h_rapi->Write();

        h_nVTX_before_PUCorr->Write();
        h_nVTX_before_EffCorr->Write();
        h_nVTX->Write();

        h_pT_lead_before_PUCorr->Write();
        h_pT_sublead_before_PUCorr->Write();
        h_eta_lead_before_PUCorr->Write();
        h_eta_sublead_before_PUCorr->Write();
        h_phi_lead_before_PUCorr->Write();
        h_phi_sublead_before_PUCorr->Write();

        h_pT_lead_before_EffCorr->Write();
        h_pT_sublead_before_EffCorr->Write();
        h_eta_lead_before_EffCorr->Write();
        h_eta_sublead_before_EffCorr->Write();
        h_phi_lead_before_EffCorr->Write();
        h_phi_sublead_before_EffCorr->Write();

        h_pT_lead_before_PVzCorr->Write();
        h_pT_sublead_before_PVzCorr->Write();
        h_eta_lead_before_PVzCorr->Write();
        h_eta_sublead_before_PVzCorr->Write();
        h_phi_lead_before_PVzCorr->Write();
        h_phi_sublead_before_PVzCorr->Write();

        h_pT_lead_before_L1Corr->Write();
        h_pT_sublead_before_L1Corr->Write();
        h_eta_lead_before_L1Corr->Write();
        h_eta_sublead_before_L1Corr->Write();
        h_phi_lead_before_L1Corr->Write();
        h_phi_sublead_before_L1Corr->Write();

        h_pT_lead->Write();
        h_pT_sublead->Write();
        h_eta_lead->Write();
        h_eta_sublead->Write();
        h_phi_lead->Write();
        h_phi_sublead->Write();

        cout << " Finished." << endl;

    }// End of i_proc iteration

    f->Close();
    if ( isBkgFull == kTRUE )
    {
        if ( !f->IsOpen() ) cout << "File Hist_" << Mgr.Procname[_EE_Bkg_Full]+debug << ".root has been closed successfully.\n" << endl;
        else cout << "FILE Hist_" << Mgr.Procname[_EE_Bkg_Full]+debug << ".root COULD NOT BE CLOSED!\n" << endl;
    }
    else
    {
        if ( !f->IsOpen() ) cout << "File Hist_" << Mgr.Procname[Mgr.CurrentProc]+debug << ".root has been closed successfully.\n" << endl;
        else cout << "FILE Hist_" << Mgr.Procname[Mgr.CurrentProc]+debug << ".root COULD NOT BE CLOSED!\n" << endl;
    }

    Double_t TotalRunTime = totaltime.CpuTime();
    cout << "Total RunTime: " << TotalRunTime << " seconds" << endl;

    TTimeStamp ts_end;
    cout << "[End Time(local time): " << ts_end.AsString("l") << "]" << endl;

} // End of EE_HistMaker()


/// -------------------------------- Muon Channel ------------------------------------ ///
void MuMu_HistMaker (TString type, Bool_t SwitchROCCORR, TString HLTname , Bool_t DEBUG)
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
    TString RocCor = "_roccor";
    if ( SwitchROCCORR == kTRUE )
    {
        Mgr.SwitchROCCORR();  // If kTRUE, this will go through events that were selected without Rochester correction applied
        RocCor = "";
    }
    TFile *f;
    TString OutputDir;
    Bool_t isBkgFull = kFALSE;  // To tell if this is the _EE_Bkg_Full process (it is handled differently)
    TString debug = "";
    if ( DEBUG == kTRUE ) debug = "_DEBUG";

    if ( !Processes.size() )
    {
        cout << "Error: no processes!" << endl;
        return;
    }

    if ( Processes[0] == _MuMu_Bkg_Full )
    {
        isBkgFull = kTRUE;
        Mgr.SetProc( _MuMu_Bkg_Full );
        // -- Output ROOTFile -- //
        OutputDir = Mgr.HistLocation;
        f = new TFile( OutputDir+"Hist_"+Mgr.Procname[_MuMu_Bkg_Full]+RocCor+debug+".root", "RECREATE" );
        Processes.clear();
        Processes.push_back( _MuMu_DYTauTau_Full );
        Processes.push_back( _MuMu_ttbar_Full );
        Processes.push_back( _MuMu_tW );
        Processes.push_back( _MuMu_tbarW );
        Processes.push_back( _MuMu_WW );
        Processes.push_back( _MuMu_WZ );
        Processes.push_back( _MuMu_ZZ );
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

        Mgr.SetProc( Processes[i_proc] );

        // -- Output ROOTFile -- //
        if ( isBkgFull == kFALSE )
        {
            OutputDir = Mgr.HistLocation;
            f = new TFile( OutputDir+"Hist_"+Mgr.Procname[Mgr.CurrentProc]+RocCor+debug+".root", "RECREATE" );
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
        analyzer->SetupEfficiencyScaleFactor_BtoF();
        analyzer->SetupEfficiencyScaleFactor_GtoH();

        // -- For PVz reweighting -- //
        analyzer->SetupPVzWeights( Mgr.isMC, "mumu", "./etc/PVzWeights.root");

        // -- Creating Histograms -- //
        TH1D *h_mass_before_PUCorr = new TH1D( "h_mass_before_PUCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 43, massbins );
        TH1D *h_mass_before_EffCorr = new TH1D( "h_mass_before_EffCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 43, massbins );
        TH1D *h_mass_before_PVzCorr = new TH1D( "h_mass_before_PVzCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 43, massbins );
        TH1D *h_mass_before_L1Corr = new TH1D( "h_mass_before_L1Corr_"+Mgr.Procname[Mgr.CurrentProc], "", 43, massbins );
        TH1D *h_mass_fine = new TH1D( "h_mass_fine_"+Mgr.Procname[Mgr.CurrentProc], "", 10000, 0, 10000 );
        TH1D *h_mass = new TH1D( "h_mass_"+Mgr.Procname[Mgr.CurrentProc], "", 43, massbins );
        TH1D *h_pT = new TH1D( "h_pT_"+Mgr.Procname[Mgr.CurrentProc], "", 100, 0, 1000 );
        TH1D *h_rapi = new TH1D( "h_rapi_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -5, 5 );

        TH1D *h_nVTX_before_PUCorr = new TH1D( "h_nVTX_before_PUCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 75, 0, 75 );
        TH1D *h_nVTX_before_EffCorr = new TH1D( "h_nVTX_before_EffCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 75, 0, 75 );
        TH1D *h_nVTX = new TH1D( "h_nVTX_"+Mgr.Procname[Mgr.CurrentProc], "", 75, 0, 75 );

        TH1D *h_pT_lead_before_PUCorr = new TH1D( "h_pT_lead_before_PUCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, 0, 1000 );
        TH1D *h_pT_sublead_before_PUCorr = new TH1D( "h_pT_sublead_before_PUCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, 0, 1000 );
        TH1D *h_eta_lead_before_PUCorr = new TH1D( "h_eta_lead_before_PUCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4 );
        TH1D *h_eta_sublead_before_PUCorr = new TH1D( "h_eta_sublead_before_PUCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4 );
        TH1D *h_phi_lead_before_PUCorr = new TH1D( "h_phi_lead_before_PUCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4 );
        TH1D *h_phi_sublead_before_PUCorr = new TH1D( "h_phi_sublead_before_PUCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4 );

        TH1D *h_pT_lead_before_EffCorr = new TH1D( "h_pT_lead_before_EffCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, 0, 1000 );
        TH1D *h_pT_sublead_before_EffCorr = new TH1D( "h_pT_sublead_before_EffCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, 0, 1000 );
        TH1D *h_eta_lead_before_EffCorr = new TH1D( "h_eta_lead_before_EffCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4 );
        TH1D *h_eta_sublead_before_EffCorr = new TH1D( "h_eta_sublead_before_EffCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4 );
        TH1D *h_phi_lead_before_EffCorr = new TH1D( "h_phi_lead_before_EffCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4 );
        TH1D *h_phi_sublead_before_EffCorr = new TH1D( "h_phi_sublead_before_EffCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4 );

        TH1D *h_pT_lead_before_PVzCorr = new TH1D( "h_pT_lead_before_PVzCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, 0, 1000 );
        TH1D *h_pT_sublead_before_PVzCorr = new TH1D( "h_pT_sublead_before_PVzCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, 0, 1000 );
        TH1D *h_eta_lead_before_PVzCorr = new TH1D( "h_eta_lead_before_PVzCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4 );
        TH1D *h_eta_sublead_before_PVzCorr = new TH1D( "h_eta_sublead_before_PVzCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4 );
        TH1D *h_phi_lead_before_PVzCorr = new TH1D( "h_phi_lead_before_PVzCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4 );
        TH1D *h_phi_sublead_before_PVzCorr = new TH1D( "h_phi_sublead_before_PVzCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4 );

        TH1D *h_pT_lead_before_L1Corr = new TH1D( "h_pT_lead_before_L1Corr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, 0, 1000 );
        TH1D *h_pT_sublead_before_L1Corr = new TH1D( "h_pT_sublead_before_L1Corr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, 0, 1000 );
        TH1D *h_eta_lead_before_L1Corr = new TH1D( "h_eta_lead_before_L1Corr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4 );
        TH1D *h_eta_sublead_before_L1Corr = new TH1D( "h_eta_sublead_before_L1Corr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4 );
        TH1D *h_phi_lead_before_L1Corr = new TH1D( "h_phi_lead_before_L1Corr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4 );
        TH1D *h_phi_sublead_before_L1Corr = new TH1D( "h_phi_sublead_before_L1Corr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4 );

        TH1D *h_pT_lead = new TH1D( "h_pT_lead_"+Mgr.Procname[Mgr.CurrentProc], "", 100, 0, 1000 );
        TH1D *h_pT_sublead = new TH1D( "h_pT_sublead_"+Mgr.Procname[Mgr.CurrentProc], "", 100, 0, 1000 );
        TH1D *h_eta_lead = new TH1D( "h_eta_lead_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4 );
        TH1D *h_eta_sublead = new TH1D( "h_eta_sublead_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4 );
        TH1D *h_phi_lead = new TH1D( "h_phi_lead_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4 );
        TH1D *h_phi_sublead = new TH1D( "h_phi_sublead_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4 );

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

            Int_t NEvents = chain->GetEntries();
            if ( NEvents != Mgr.nEvents[i_tup] ) cout << "\tEvent numbers do not match!!!" << endl;
            cout << "\t[Sum of weights: " << Mgr.Wsum[i_tup] << "]" << endl;
            cout << "\t[Selected Events: " << NEvents << "]" << endl;

            if ( DEBUG == kTRUE ) NEvents = 10;

            myProgressBar_t bar( NEvents );

            for(Int_t i=0; i<NEvents; i++)
            {              
                MuMu->GetEvent(i);

                if ( MuMu->Muon_eta->size() != 2 )
                {
                    cout << "eta vector has " << MuMu->Muon_eta->size() << " elements!!!" << endl;
                    break;
                }
                else if ( fabs(MuMu->Muon_eta->at(0)) > 2.4 || fabs(MuMu->Muon_eta->at(1)) > 2.4 )
                {
                    cout << "Muon etas:  [0] " << MuMu->Muon_eta->at(0) << "    [1] " << MuMu->Muon_eta->at(1) << endl;
                    continue;
                }
                if ( MuMu->Muon_charge->size() != 2 )
                {
                    cout << "Charge vector has " << MuMu->Muon_charge->size() << " elements!!!" << endl;
                    break;
                }
                if ( MuMu->Muon_pT->size() != 2 )
                {
                    cout << "pT vector has " << MuMu->Muon_pT->size() << " elements!!!" << endl;
                    break;
                }
                if ( MuMu->Muon_phi->size() != 2 )
                {
                    cout << "phi vector has " << MuMu->Muon_phi->size() << " elements!!!" << endl;
                    break;
                }
                if ( MuMu->Muon_Energy->size() != 2 )
                {
                    cout << "Energy vector has " << MuMu->Muon_Energy->size() << " elements!!!" << endl;
                    break;
                }
                if ( MuMu->Muon_TuneP_pT->size() != 2 )
                {
                    cout << "TuneP_pT vector has " << MuMu->Muon_TuneP_pT->size() << " elements!!!" << endl;
                    break;
                }
                if ( MuMu->Muon_TuneP_eta->size() != 2 )
                {
                    cout << "TuneP_eta vector has " << MuMu->Muon_TuneP_eta->size() << " elements!!!" << endl;
                    break;
                }
                if ( MuMu->Muon_TuneP_phi->size() != 2 )
                {
                    cout << "TuneP_phi vector has " << MuMu->Muon_TuneP_phi->size() << " elements!!!" << endl;
                    break;
                }


                Double_t GenWeight = MuMu->GENEvt_weight;            

                // -- Pileup-Reweighting -- //
                Double_t PUWeight = 1;
                if( Mgr.isMC == kTRUE ) PUWeight = analyzer->PileUpWeightValue_80X( MuMu->nPileUp );

                // -- efficiency weights -- //
                Double_t weight1 = 0, weight2 = 0, effweight = 1;

                // -- PVz weights -- //
                Double_t PVzWeight = 1;
                if ( Mgr.isMC == kTRUE ) PVzWeight = analyzer->PVzWeightValue( MuMu->PVz );

                // -- L1 prefiring weights -- //
                Double_t L1weight = 1;
                if ( Mgr.isMC == kTRUE ) L1weight = MuMu->_prefiringweight;
//                if ( Mgr.isMC == kTRUE ) L1weight = MuMu->_prefiringweightup;
//                if ( Mgr.isMC == kTRUE ) L1weight = MuMu->_prefiringweightdown;

                // -- Normalization -- //
                Double_t TotWeight = GenWeight;
                if( Mgr.isMC == kTRUE ) TotWeight = ( L * Mgr.Xsec[i_tup] / Mgr.Wsum[i_tup] ) * GenWeight;

                if( MuMu->isSelPassed == kTRUE )
                {
                    TLorentzVector mu1, mu2;
                    // Temporarily using TuneP_phi instead of Phi, because of the bug (wrote eta into phi)
                    mu1.SetPtEtaPhiE( MuMu->Muon_pT->at(0), MuMu->Muon_eta->at(0), MuMu->Muon_TuneP_phi->at(0), MuMu->Muon_Energy->at(0) );
                    mu2.SetPtEtaPhiE( MuMu->Muon_pT->at(1), MuMu->Muon_eta->at(1), MuMu->Muon_TuneP_phi->at(1), MuMu->Muon_Energy->at(1) );

                    // -- Apply efficiency scale factor -- //
                    if( Mgr.isMC == kTRUE )
                    {
                        weight1 = analyzer->EfficiencySF_EventWeight_HLT_BtoF( MuMu );
                        weight2 = analyzer->EfficiencySF_EventWeight_HLT_GtoH( MuMu );
                        effweight = ( L_B2F * weight1 + L_G2H * weight2 ) / L_B2H;
                    }

                    Double_t reco_Pt = ( mu1 + mu2 ).Pt();
                    Double_t reco_rapi = ( mu1 + mu2 ).Rapidity();
                    Double_t reco_mass = MuMu->Muon_InvM;

                    h_mass_before_PUCorr->Fill( reco_mass, TotWeight );
                    h_mass_before_EffCorr->Fill( reco_mass, TotWeight * PUWeight );
                    h_mass_before_PVzCorr->Fill( reco_mass, TotWeight * PUWeight * effweight );
                    h_mass_before_L1Corr->Fill( reco_mass, TotWeight * PUWeight * effweight * PVzWeight );
                    h_mass_fine->Fill( reco_mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight );
                    h_mass->Fill( reco_mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight );
                    h_pT->Fill( reco_Pt, TotWeight * PUWeight * effweight * PVzWeight * L1weight );
                    h_rapi->Fill( reco_rapi, TotWeight * PUWeight * effweight * PVzWeight * L1weight );

                    h_nVTX_before_PUCorr->Fill( MuMu->nVertices, TotWeight );
                    h_nVTX_before_EffCorr->Fill( MuMu->nVertices, TotWeight * PUWeight );
                    h_nVTX->Fill( MuMu->nVertices, TotWeight * PUWeight * effweight * PVzWeight * L1weight );

                    int lead=0, sublead=1;
                    if ( MuMu->Muon_pT->at(0) < MuMu->Muon_pT->at(1) )
                    {
                        lead = 1;
                        sublead = 0;
                    }

                    h_pT_lead_before_PUCorr->Fill( MuMu->Muon_pT->at(lead), TotWeight );
                    h_pT_sublead_before_PUCorr->Fill( MuMu->Muon_pT->at(sublead), TotWeight );
                    h_eta_lead_before_PUCorr->Fill( MuMu->Muon_eta->at(lead), TotWeight );
                    h_eta_sublead_before_PUCorr->Fill( MuMu->Muon_eta->at(sublead), TotWeight );
                    h_phi_lead_before_PUCorr->Fill( MuMu->Muon_TuneP_phi->at(lead), TotWeight );
                    h_phi_sublead_before_PUCorr->Fill( MuMu->Muon_TuneP_phi->at(sublead), TotWeight );

                    h_pT_lead_before_EffCorr->Fill( MuMu->Muon_pT->at(lead), TotWeight * PUWeight );
                    h_pT_sublead_before_EffCorr->Fill( MuMu->Muon_pT->at(sublead), TotWeight * PUWeight );
                    h_eta_lead_before_EffCorr->Fill( MuMu->Muon_eta->at(lead), TotWeight * PUWeight );
                    h_eta_sublead_before_EffCorr->Fill( MuMu->Muon_eta->at(sublead), TotWeight * PUWeight );
                    h_phi_lead_before_EffCorr->Fill( MuMu->Muon_TuneP_phi->at(lead), TotWeight * PUWeight );
                    h_phi_sublead_before_EffCorr->Fill( MuMu->Muon_TuneP_phi->at(sublead), TotWeight * PUWeight );

                    h_pT_lead_before_PVzCorr->Fill( MuMu->Muon_pT->at(lead), TotWeight * PUWeight * effweight );
                    h_pT_sublead_before_PVzCorr->Fill( MuMu->Muon_pT->at(sublead), TotWeight * PUWeight * effweight );
                    h_eta_lead_before_PVzCorr->Fill( MuMu->Muon_eta->at(lead), TotWeight * PUWeight * effweight );
                    h_eta_sublead_before_PVzCorr->Fill( MuMu->Muon_eta->at(sublead), TotWeight * PUWeight * effweight );
                    h_phi_lead_before_PVzCorr->Fill( MuMu->Muon_TuneP_phi->at(lead), TotWeight * PUWeight * effweight );
                    h_phi_sublead_before_PVzCorr->Fill( MuMu->Muon_TuneP_phi->at(sublead), TotWeight * PUWeight * effweight );

                    h_pT_lead_before_L1Corr->Fill( MuMu->Muon_pT->at(lead), TotWeight * PUWeight * effweight * PVzWeight );
                    h_pT_sublead_before_L1Corr->Fill( MuMu->Muon_pT->at(sublead), TotWeight * PUWeight * effweight * PVzWeight );
                    h_eta_lead_before_L1Corr->Fill( MuMu->Muon_eta->at(lead), TotWeight * PUWeight * effweight * PVzWeight );
                    h_eta_sublead_before_L1Corr->Fill( MuMu->Muon_eta->at(sublead), TotWeight * PUWeight * effweight * PVzWeight );
                    h_phi_lead_before_L1Corr->Fill( MuMu->Muon_TuneP_phi->at(lead), TotWeight * PUWeight * effweight * PVzWeight );
                    h_phi_sublead_before_L1Corr->Fill( MuMu->Muon_TuneP_phi->at(sublead), TotWeight * PUWeight * effweight * PVzWeight );

                    h_pT_lead->Fill( MuMu->Muon_pT->at(lead), TotWeight * PUWeight * effweight * PVzWeight * L1weight );
                    h_pT_sublead->Fill( MuMu->Muon_pT->at(sublead), TotWeight * PUWeight * effweight * PVzWeight * L1weight );
                    h_eta_lead->Fill( MuMu->Muon_eta->at(lead), TotWeight * PUWeight * effweight * PVzWeight * L1weight );
                    h_eta_sublead->Fill( MuMu->Muon_eta->at(sublead), TotWeight * PUWeight * effweight * PVzWeight * L1weight );
                    h_phi_lead->Fill( MuMu->Muon_TuneP_phi->at(lead), TotWeight * PUWeight * effweight * PVzWeight * L1weight );
                    h_phi_sublead->Fill( MuMu->Muon_TuneP_phi->at(sublead), TotWeight * PUWeight * effweight * PVzWeight * L1weight );

                }// End of event selection
                bar.Draw(i);

            }// End of event iteration

            if( Mgr.isMC == kTRUE ) printf( "\tNormalization factor: %.8f\n", L*Mgr.Xsec[i_tup]/Mgr.Wsum[i_tup] );

            Double_t LoopRunTime = looptime.CpuTime();
            cout << "\tLoop RunTime(" << Mgr.Tag[i_tup] << "): " << LoopRunTime << " seconds\n" << endl;

        }// End of i_tup iteration

        f->cd();
        cout << "\tWriting into file...";

        h_mass_before_PUCorr->Write();
        h_mass_before_EffCorr->Write();
        h_mass_before_PVzCorr->Write();
        h_mass_before_L1Corr->Write();
        h_mass_fine->Write();
        h_mass->Write();
        h_pT->Write();
        h_rapi->Write();

        h_nVTX_before_PUCorr->Write();
        h_nVTX_before_EffCorr->Write();
        h_nVTX->Write();

        h_pT_lead_before_PUCorr->Write();
        h_pT_sublead_before_PUCorr->Write();
        h_eta_lead_before_PUCorr->Write();
        h_eta_sublead_before_PUCorr->Write();
        h_phi_lead_before_PUCorr->Write();
        h_phi_sublead_before_PUCorr->Write();

        h_pT_lead_before_EffCorr->Write();
        h_pT_sublead_before_EffCorr->Write();
        h_eta_lead_before_EffCorr->Write();
        h_eta_sublead_before_EffCorr->Write();
        h_phi_lead_before_EffCorr->Write();
        h_phi_sublead_before_EffCorr->Write();

        h_pT_lead_before_PVzCorr->Write();
        h_pT_sublead_before_PVzCorr->Write();
        h_eta_lead_before_PVzCorr->Write();
        h_eta_sublead_before_PVzCorr->Write();
        h_phi_lead_before_PVzCorr->Write();
        h_phi_sublead_before_PVzCorr->Write();

        h_pT_lead_before_L1Corr->Write();
        h_pT_sublead_before_L1Corr->Write();
        h_eta_lead_before_L1Corr->Write();
        h_eta_sublead_before_L1Corr->Write();
        h_phi_lead_before_L1Corr->Write();
        h_phi_sublead_before_L1Corr->Write();

        h_pT_lead->Write();
        h_pT_sublead->Write();
        h_eta_lead->Write();
        h_eta_sublead->Write();
        h_phi_lead->Write();
        h_phi_sublead->Write();

        cout << " Finished." << endl;

    } // End of i_proc iteration

    f->Close();
    if ( isBkgFull == kTRUE )
    {
        if ( !f->IsOpen() ) cout << "File Hist_" << Mgr.Procname[_MuMu_Bkg_Full]+RocCor+debug << ".root has been closed successfully.\n" << endl;
        else cout << "FILE Hist_" << Mgr.Procname[_MuMu_Bkg_Full]+RocCor+debug << ".root COULD NOT BE CLOSED!\n" << endl;
    }
    else
    {
        if ( !f->IsOpen() ) cout << "File Hist_" << Mgr.Procname[Mgr.CurrentProc]+RocCor+debug << ".root has been closed successfully.\n" << endl;
        else cout << "FILE Hist_" << Mgr.Procname[Mgr.CurrentProc]+RocCor+debug << ".root COULD NOT BE CLOSED!\n" << endl;
    }

    Double_t TotalRunTime = totaltime.CpuTime();
    cout << "Total RunTime: " << TotalRunTime << " seconds" << endl;

    TTimeStamp ts_end;
    cout << "[End Time(local time): " << ts_end.AsString("l") << "]" << endl;

} // End of MuMu_HistMaker()


/// --------------------------------- EMu events --------------------------------- ///
void EMu_HistMaker (TString type, Bool_t SwitchROCCORR, TString HLTname , Bool_t DEBUG)
{
    if ( !type.Length() )
    {
        cout << "Error: no type specified!" << endl;
        return;
    }

    // -- Run2016 luminosity [/pb] -- //
    Double_t L_B2F = 19721.0, L_G2H = 16146.0, L_B2H = 35867.0, L = 0;
    L = L_B2H;
    Int_t nWJetsSS=0, nWJetsSS_weighted=0;

    TTimeStamp ts_start;
    cout << "[Start Time(local time): " << ts_start.AsString("l") << "]" << endl;
    TStopwatch totaltime;
    totaltime.Start();

    LocalFileMgr Mgr;
    vector<SelProc_t> Processes = Mgr.FindProc( type );
    TString RocCor = "_roccor";
    if ( SwitchROCCORR == kTRUE )
    {
        Mgr.SwitchROCCORR();  // If kTRUE, this will go through events that were selected without Rochester correction applied
        RocCor = "";
    }
    TFile *f;
    TString OutputDir;
    Bool_t isBkgFull = kFALSE;  // To tell if this is the _EE_Bkg_Full process (it is handled differently)
    TString debug = "";
    if ( DEBUG == kTRUE ) debug = "_DEBUG";

    if ( !Processes.size() )
    {
        cout << "Error: no processes!" << endl;
        return;
    }

    if ( Processes[0] == _EMu_Bkg_Full )
    {
        isBkgFull = kTRUE;
        Mgr.SetProc( _EMu_Bkg_Full );
        // -- Output ROOTFile -- //
        OutputDir = Mgr.HistLocation;
        f = new TFile( OutputDir+"Hist_"+Mgr.Procname[_EMu_Bkg_Full]+RocCor+debug+".root", "RECREATE" );
        Processes.clear();
        Processes.push_back( _EMu_DYTauTau_Full );
        Processes.push_back( _EMu_ttbar_Full );
        Processes.push_back( _EMu_tW );
        Processes.push_back( _EMu_tbarW );
        Processes.push_back( _EMu_WW );
        Processes.push_back( _EMu_WZ );
        Processes.push_back( _EMu_ZZ );
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

        Mgr.SetProc( Processes[i_proc] );

        // -- Output ROOTFile -- //
        if ( isBkgFull == kFALSE )
        {
            OutputDir = Mgr.HistLocation;
            f = new TFile( OutputDir+"Hist_"+Mgr.Procname[Mgr.CurrentProc]+RocCor+debug+".root", "RECREATE" );
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
        analyzer->SetupEfficiencyScaleFactor_BtoF();
        analyzer->SetupEfficiencyScaleFactor_GtoH();
        analyzer->SetupEfficiencyScaleFactor_electron();

        // -- For PVz reweighting -- //
        analyzer->SetupPVzWeights( Mgr.isMC, "mumu", "./etc/PVzWeights.root");

        // -- Creating Histograms -- //
        TH1D *h_emu_mass_before_PUCorr = new TH1D( "h_emu_mass_before_PUCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 43, massbins );
        TH1D *h_emu_mass_before_EffCorr = new TH1D( "h_emu_mass_before_EffCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 43, massbins );
        TH1D *h_emu_mass_before_PVzCorr = new TH1D( "h_emu_mass_before_PVzCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 43, massbins );
        TH1D *h_emu_mass_before_L1Corr = new TH1D( "h_emu_mass_before_L1Corr_"+Mgr.Procname[Mgr.CurrentProc], "", 43, massbins );
        TH1D *h_emu_mass_fine = new TH1D( "h_emu_mass_fine_"+Mgr.Procname[Mgr.CurrentProc], "", 10000, 0, 10000 );
        TH1D *h_emu_mass = new TH1D( "h_emu_mass_"+Mgr.Procname[Mgr.CurrentProc], "", 43, massbins );
        TH1D *h_emuSS_mass_before_PUCorr = new TH1D( "h_emuSS_mass_before_PUCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 43, massbins );
        TH1D *h_emuSS_mass_before_EffCorr = new TH1D( "h_emuSS_mass_before_EffCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 43, massbins );
        TH1D *h_emuSS_mass_before_PVzCorr = new TH1D( "h_emuSS_mass_before_PVzCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 43, massbins );
        TH1D *h_emuSS_mass_before_L1Corr = new TH1D( "h_emuSS_mass_before_L1Corr_"+Mgr.Procname[Mgr.CurrentProc], "", 43, massbins );
        TH1D *h_emuSS_mass_fine = new TH1D( "h_emuSS_mass_fine_"+Mgr.Procname[Mgr.CurrentProc], "", 10000, 0, 10000 );
        TH1D *h_emuSS_mass = new TH1D( "h_emuSS_mass_"+Mgr.Procname[Mgr.CurrentProc], "", 43, massbins );

        TH1D *h_nVTX_before_PUCorr = new TH1D( "h_nVTX_before_PUCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 75, 0, 75 );
        TH1D *h_nVTX_before_EffCorr = new TH1D( "h_nVTX_before_EffCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 75, 0, 75 );
        TH1D *h_nVTX = new TH1D( "h_nVTX_"+Mgr.Procname[Mgr.CurrentProc], "", 75, 0, 75 );

        TH1D *h_ele_pT_before_PUCorr = new TH1D( "h_ele_pT_before_PUCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, 0, 1000 );
        TH1D *h_ele_eta_before_PUCorr = new TH1D( "h_ele_eta_before_PUCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4 );
        TH1D *h_ele_phi_before_PUCorr = new TH1D( "h_ele_phi_before_PUCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4 );
        TH1D *h_eleSS_pT_before_PUCorr = new TH1D( "h_eleSS_pT_before_PUCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, 0, 1000 );
        TH1D *h_eleSS_eta_before_PUCorr = new TH1D( "h_eleSS_eta_before_PUCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4 );
        TH1D *h_eleSS_phi_before_PUCorr = new TH1D( "h_eleSS_phi_before_PUCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4 );

        TH1D *h_ele_pT_before_EffCorr = new TH1D( "h_ele_pT_before_EffCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, 0, 1000 );
        TH1D *h_ele_eta_before_EffCorr = new TH1D( "h_ele_eta_before_EffCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4 );
        TH1D *h_ele_phi_before_EffCorr = new TH1D( "h_ele_phi_before_EffCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4 );
        TH1D *h_eleSS_pT_before_EffCorr = new TH1D( "h_eleSS_pT_before_EffCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, 0, 1000 );
        TH1D *h_eleSS_eta_before_EffCorr = new TH1D( "h_eleSS_eta_before_EffCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4 );
        TH1D *h_eleSS_phi_before_EffCorr = new TH1D( "h_eleSS_phi_before_EffCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4 );

        TH1D *h_ele_pT_before_PVzCorr = new TH1D( "h_ele_pT_before_PVzCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, 0, 1000 );
        TH1D *h_ele_eta_before_PVzCorr = new TH1D( "h_ele_eta_before_PVzCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4 );
        TH1D *h_ele_phi_before_PVzCorr = new TH1D( "h_ele_phi_before_PVzCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4 );
        TH1D *h_eleSS_pT_before_PVzCorr = new TH1D( "h_eleSS_pT_before_PVzCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, 0, 1000 );
        TH1D *h_eleSS_eta_before_PVzCorr = new TH1D( "h_eleSS_eta_before_PVzCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4 );
        TH1D *h_eleSS_phi_before_PVzCorr = new TH1D( "h_eleSS_phi_before_PVzCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4 );

        TH1D *h_ele_pT_before_L1Corr = new TH1D( "h_ele_pT_before_L1Corr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, 0, 1000 );
        TH1D *h_ele_eta_before_L1Corr = new TH1D( "h_ele_eta_before_L1Corr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4 );
        TH1D *h_ele_phi_before_L1Corr = new TH1D( "h_ele_phi_before_L1Corr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4 );
        TH1D *h_eleSS_pT_before_L1Corr = new TH1D( "h_eleSS_pT_before_L1Corr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, 0, 1000 );
        TH1D *h_eleSS_eta_before_L1Corr = new TH1D( "h_eleSS_eta_before_L1Corr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4 );
        TH1D *h_eleSS_phi_before_L1Corr = new TH1D( "h_eleSS_phi_before_L1Corr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4 );

        TH1D *h_ele_pT = new TH1D( "h_ele_pT_"+Mgr.Procname[Mgr.CurrentProc], "", 100, 0, 1000 );
        TH1D *h_ele_eta = new TH1D( "h_ele_eta_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4 );
        TH1D *h_ele_phi = new TH1D( "h_ele_phi_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4 );
        TH1D *h_eleSS_pT = new TH1D( "h_eleSS_pT_"+Mgr.Procname[Mgr.CurrentProc], "", 100, 0, 1000 );
        TH1D *h_eleSS_eta = new TH1D( "h_eleSS_eta_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4 );
        TH1D *h_eleSS_phi = new TH1D( "h_eleSS_phi_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4 );

        TH1D *h_mu_pT_before_PUCorr = new TH1D( "h_mu_pT_before_PUCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, 0, 1000 );
        TH1D *h_mu_eta_before_PUCorr = new TH1D( "h_mu_eta_before_PUCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4 );
        TH1D *h_mu_phi_before_PUCorr = new TH1D( "h_mu_phi_before_PUCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4 );
        TH1D *h_muSS_pT_before_PUCorr = new TH1D( "h_muSS_pT_before_PUCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, 0, 1000 );
        TH1D *h_muSS_eta_before_PUCorr = new TH1D( "h_muSS_eta_before_PUCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4 );
        TH1D *h_muSS_phi_before_PUCorr = new TH1D( "h_muSS_phi_before_PUCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4 );

        TH1D *h_mu_pT_before_EffCorr = new TH1D( "h_mu_pT_before_EffCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, 0, 1000 );
        TH1D *h_mu_eta_before_EffCorr = new TH1D( "h_mu_eta_before_EffCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4 );
        TH1D *h_mu_phi_before_EffCorr = new TH1D( "h_mu_phi_before_EffCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4 );
        TH1D *h_muSS_pT_before_EffCorr = new TH1D( "h_muSS_pT_before_EffCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, 0, 1000 );
        TH1D *h_muSS_eta_before_EffCorr = new TH1D( "h_muSS_eta_before_EffCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4 );
        TH1D *h_muSS_phi_before_EffCorr = new TH1D( "h_muSS_phi_before_EffCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4 );

        TH1D *h_mu_pT_before_PVzCorr = new TH1D( "h_mu_pT_before_PVzCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, 0, 1000 );
        TH1D *h_mu_eta_before_PVzCorr = new TH1D( "h_mu_eta_before_PVzCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4 );
        TH1D *h_mu_phi_before_PVzCorr = new TH1D( "h_mu_phi_before_PVzCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4 );
        TH1D *h_muSS_pT_before_PVzCorr = new TH1D( "h_muSS_pT_before_PVzCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, 0, 1000 );
        TH1D *h_muSS_eta_before_PVzCorr = new TH1D( "h_muSS_eta_before_PVzCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4 );
        TH1D *h_muSS_phi_before_PVzCorr = new TH1D( "h_muSS_phi_before_PVzCorr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4 );

        TH1D *h_mu_pT_before_L1Corr = new TH1D( "h_mu_pT_before_L1Corr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, 0, 1000 );
        TH1D *h_mu_eta_before_L1Corr = new TH1D( "h_mu_eta_before_L1Corr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4 );
        TH1D *h_mu_phi_before_L1Corr = new TH1D( "h_mu_phi_before_L1Corr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4 );
        TH1D *h_muSS_pT_before_L1Corr = new TH1D( "h_muSS_pT_before_L1Corr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, 0, 1000 );
        TH1D *h_muSS_eta_before_L1Corr = new TH1D( "h_muSS_eta_before_L1Corr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4 );
        TH1D *h_muSS_phi_before_L1Corr = new TH1D( "h_muSS_phi_before_L1Corr_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4 );

        TH1D *h_mu_pT = new TH1D( "h_mu_pT_"+Mgr.Procname[Mgr.CurrentProc], "", 100, 0, 1000 );
        TH1D *h_mu_eta = new TH1D( "h_mu_eta_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4 );
        TH1D *h_mu_phi = new TH1D( "h_mu_phi_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4 );
        TH1D *h_muSS_pT = new TH1D( "h_muSS_pT_"+Mgr.Procname[Mgr.CurrentProc], "", 100, 0, 1000 );
        TH1D *h_muSS_eta = new TH1D( "h_muSS_eta_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4 );
        TH1D *h_muSS_phi = new TH1D( "h_muSS_phi_"+Mgr.Procname[Mgr.CurrentProc], "", 100, -4, 4 );

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

            Int_t NEvents = chain->GetEntries();
            if ( NEvents != Mgr.nEvents[i_tup] ) cout << "\tEvent numbers do not match!!!" << endl;
            cout << "\t[Sum of weights: " << Mgr.Wsum[i_tup] << "]" << endl;
            cout << "\t[Selected Events: " << NEvents << "]" << endl;

            if ( DEBUG == kTRUE ) NEvents = 10;

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

                // -- PVz weights -- //
                Double_t PVzWeight = 1;
                if( Mgr.isMC == kTRUE ) PVzWeight = analyzer->PVzWeightValue( EMu->PVz );

                // -- L1 weights --//
                Double_t L1weight = 1;
                if( Mgr.isMC == kTRUE ) L1weight = EMu->_prefiringweight;
//                if( Mgr.isMC == kTRUE ) L1weight = EMu->_prefiringweightup;
//                if( Mgr.isMC == kTRUE ) L1weight = EMu->_prefiringweightdown;

                // -- Normalization -- //
                Double_t TotWeight = GenWeight;
                if( Mgr.isMC == kTRUE ) TotWeight = ( L * Mgr.Xsec[i_tup] / Mgr.Wsum[i_tup] ) * GenWeight;

                if( EMu->isSelPassed == kTRUE )
                {
                    TLorentzVector mu, ele;
                    mu.SetPtEtaPhiE( EMu->Muon_pT, EMu->Muon_eta, EMu->Muon_phi, EMu->Muon_Energy );
                    ele.SetPtEtaPhiE( EMu->Electron_pT, EMu->Electron_eta, EMu->Electron_phi, EMu->Electron_Energy );

                    // -- Apply efficiency scale factor -- //
                    if( Mgr.isMC == kTRUE )
                    {
                        weight1 = analyzer->EfficiencySF_EventWeight_emu_BtoF( EMu );
                        weight2 = analyzer->EfficiencySF_EventWeight_emu_GtoH( EMu );
                        effweight = (L_B2F*weight1 + L_G2H*weight2)/L_B2H;
                    }

                    Double_t reco_mass = ( mu + ele ).M();

                    if ( EMu->Electron_charge != EMu->Muon_charge )
                    {
                        h_emu_mass_before_PUCorr->Fill( reco_mass, TotWeight );
                        h_emu_mass_before_EffCorr->Fill( reco_mass, TotWeight * PUWeight );
                        h_emu_mass_before_PVzCorr->Fill( reco_mass, TotWeight * PUWeight * effweight );
                        h_emu_mass_before_L1Corr->Fill( reco_mass, TotWeight * PUWeight * effweight * PVzWeight );
                        h_emu_mass_fine->Fill( reco_mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight );
                        h_emu_mass->Fill( reco_mass, TotWeight * PUWeight * effweight* PVzWeight * L1weight );

                        h_ele_pT_before_PUCorr->Fill( EMu->Electron_pT, TotWeight );
                        h_ele_eta_before_PUCorr->Fill( EMu->Electron_etaSC, TotWeight );
                        h_ele_phi_before_PUCorr->Fill( EMu->Electron_phi, TotWeight );

                        h_ele_pT_before_EffCorr->Fill( EMu->Electron_pT, TotWeight * PUWeight );
                        h_ele_eta_before_EffCorr->Fill( EMu->Electron_etaSC, TotWeight * PUWeight );
                        h_ele_phi_before_EffCorr->Fill( EMu->Electron_phi, TotWeight * PUWeight );

                        h_ele_pT_before_PVzCorr->Fill( EMu->Electron_pT, TotWeight * PUWeight * effweight );
                        h_ele_eta_before_PVzCorr->Fill( EMu->Electron_etaSC, TotWeight * PUWeight * effweight );
                        h_ele_phi_before_PVzCorr->Fill( EMu->Electron_phi, TotWeight * PUWeight * effweight );

                        h_ele_pT_before_L1Corr->Fill( EMu->Electron_pT, TotWeight * PUWeight * effweight * PVzWeight );
                        h_ele_eta_before_L1Corr->Fill( EMu->Electron_etaSC, TotWeight * PUWeight * effweight * PVzWeight );
                        h_ele_phi_before_L1Corr->Fill( EMu->Electron_phi, TotWeight * PUWeight * effweight * PVzWeight );

                        h_ele_pT->Fill( EMu->Electron_pT, TotWeight * PUWeight * effweight * PVzWeight * L1weight );
                        h_ele_eta->Fill( EMu->Electron_etaSC, TotWeight * PUWeight * effweight * PVzWeight * L1weight );
                        h_ele_phi->Fill( EMu->Electron_phi, TotWeight * PUWeight * effweight * PVzWeight * L1weight );

                        h_mu_pT_before_PUCorr->Fill( EMu->Muon_pT, TotWeight );
                        h_mu_eta_before_PUCorr->Fill( EMu->Muon_eta, TotWeight );
                        h_mu_phi_before_PUCorr->Fill( EMu->Muon_phi, TotWeight );

                        h_mu_pT_before_EffCorr->Fill( EMu->Muon_pT, TotWeight * PUWeight );
                        h_mu_eta_before_EffCorr->Fill( EMu->Muon_eta, TotWeight * PUWeight );
                        h_mu_phi_before_EffCorr->Fill( EMu->Muon_phi, TotWeight * PUWeight );

                        h_mu_pT_before_PVzCorr->Fill( EMu->Muon_pT, TotWeight * PUWeight * effweight );
                        h_mu_eta_before_PVzCorr->Fill( EMu->Muon_eta, TotWeight * PUWeight * effweight );
                        h_mu_phi_before_PVzCorr->Fill( EMu->Muon_phi, TotWeight * PUWeight * effweight );

                        h_mu_pT_before_L1Corr->Fill( EMu->Muon_pT, TotWeight * PUWeight * effweight * PVzWeight );
                        h_mu_eta_before_L1Corr->Fill( EMu->Muon_eta, TotWeight * PUWeight * effweight * PVzWeight );
                        h_mu_phi_before_L1Corr->Fill( EMu->Muon_phi, TotWeight * PUWeight * effweight * PVzWeight );

                        h_mu_pT->Fill( EMu->Muon_pT, TotWeight * PUWeight * effweight * PVzWeight * L1weight );
                        h_mu_eta->Fill( EMu->Muon_eta, TotWeight * PUWeight * effweight * PVzWeight * L1weight );
                        h_mu_phi->Fill( EMu->Muon_phi, TotWeight * PUWeight * effweight * PVzWeight * L1weight );
                    }
                    else
                    {
                        h_emuSS_mass_before_PUCorr->Fill( reco_mass, TotWeight );
                        h_emuSS_mass_before_EffCorr->Fill( reco_mass, TotWeight * PUWeight );
                        h_emuSS_mass_before_PVzCorr->Fill( reco_mass, TotWeight * PUWeight * effweight );
                        h_emuSS_mass_before_L1Corr->Fill( reco_mass, TotWeight * PUWeight * effweight * PVzWeight );
                        h_emuSS_mass_fine->Fill( reco_mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight );
                        h_emuSS_mass->Fill( reco_mass, TotWeight * PUWeight * effweight * PVzWeight * L1weight );

                        h_eleSS_pT_before_PUCorr->Fill( EMu->Electron_pT, TotWeight );
                        h_eleSS_eta_before_PUCorr->Fill( EMu->Electron_etaSC, TotWeight );
                        h_eleSS_phi_before_PUCorr->Fill( EMu->Electron_phi, TotWeight );

                        h_eleSS_pT_before_EffCorr->Fill( EMu->Electron_pT, TotWeight * PUWeight );
                        h_eleSS_eta_before_EffCorr->Fill( EMu->Electron_etaSC, TotWeight * PUWeight );
                        h_eleSS_phi_before_EffCorr->Fill( EMu->Electron_phi, TotWeight * PUWeight );

                        h_eleSS_pT_before_PVzCorr->Fill( EMu->Electron_pT, TotWeight * PUWeight * effweight );
                        h_eleSS_eta_before_PVzCorr->Fill( EMu->Electron_etaSC, TotWeight * PUWeight * effweight );
                        h_eleSS_phi_before_PVzCorr->Fill( EMu->Electron_phi, TotWeight * PUWeight * effweight );

                        h_eleSS_pT_before_L1Corr->Fill( EMu->Electron_pT, TotWeight * PUWeight * effweight * PVzWeight );
                        h_eleSS_eta_before_L1Corr->Fill( EMu->Electron_etaSC, TotWeight * PUWeight * effweight * PVzWeight );
                        h_eleSS_phi_before_L1Corr->Fill( EMu->Electron_phi, TotWeight * PUWeight * effweight * PVzWeight );

                        h_eleSS_pT->Fill( EMu->Electron_pT, TotWeight * PUWeight * effweight * PVzWeight * L1weight );
                        h_eleSS_eta->Fill( EMu->Electron_etaSC, TotWeight * PUWeight * effweight * PVzWeight * L1weight );
                        h_eleSS_phi->Fill( EMu->Electron_phi, TotWeight * PUWeight * effweight * PVzWeight * L1weight );

                        h_muSS_pT_before_PUCorr->Fill( EMu->Muon_TuneP_pT, TotWeight );
                        h_muSS_eta_before_PUCorr->Fill( EMu->Muon_TuneP_eta, TotWeight );
                        h_muSS_phi_before_PUCorr->Fill( EMu->Muon_TuneP_phi, TotWeight );

                        h_muSS_pT_before_EffCorr->Fill( EMu->Muon_TuneP_pT, TotWeight * PUWeight );
                        h_muSS_eta_before_EffCorr->Fill( EMu->Muon_TuneP_eta, TotWeight * PUWeight );
                        h_muSS_phi_before_EffCorr->Fill( EMu->Muon_TuneP_phi, TotWeight * PUWeight );

                        h_muSS_pT_before_PVzCorr->Fill( EMu->Muon_TuneP_pT, TotWeight * PUWeight * effweight );
                        h_muSS_eta_before_PVzCorr->Fill( EMu->Muon_TuneP_eta, TotWeight * PUWeight * effweight );
                        h_muSS_phi_before_PVzCorr->Fill( EMu->Muon_TuneP_phi, TotWeight * PUWeight * effweight );

                        h_muSS_pT_before_L1Corr->Fill( EMu->Muon_TuneP_pT, TotWeight * PUWeight * effweight * PVzWeight );
                        h_muSS_eta_before_L1Corr->Fill( EMu->Muon_TuneP_eta, TotWeight * PUWeight * effweight * PVzWeight );
                        h_muSS_phi_before_L1Corr->Fill( EMu->Muon_TuneP_phi, TotWeight * PUWeight * effweight * PVzWeight );

                        h_muSS_pT->Fill( EMu->Muon_TuneP_pT, TotWeight * PUWeight * effweight * PVzWeight * L1weight );
                        h_muSS_eta->Fill( EMu->Muon_TuneP_eta, TotWeight * PUWeight * effweight * PVzWeight * L1weight );
                        h_muSS_phi->Fill( EMu->Muon_TuneP_phi, TotWeight * PUWeight * effweight * PVzWeight * L1weight );

                        if ( Mgr.CurrentProc == _EMu_WJets ) { nWJetsSS++; nWJetsSS_weighted += TotWeight * PUWeight * effweight * PVzWeight * L1weight; }
                    }
                    h_nVTX_before_PUCorr->Fill( EMu->nVertices, TotWeight );
                    h_nVTX_before_EffCorr->Fill( EMu->nVertices, TotWeight * PUWeight );
                    h_nVTX->Fill( EMu->nVertices, TotWeight * PUWeight * effweight * PVzWeight * L1weight );

                }// End of event selection
                bar.Draw(i);

            }// End of event iteration

            if( Mgr.isMC == kTRUE ) printf( "\tNormalization factor: %.8f\n", L*Mgr.Xsec[i_tup]/Mgr.Wsum[i_tup] );

            Double_t LoopRunTime = looptime.CpuTime();
            cout << "\tLoop RunTime(" << Mgr.Tag[i_tup] << "): " << LoopRunTime << " seconds\n" << endl;

        }// End of i_tup iteration

        cout << "Physical SS WJets events: " << nWJetsSS << endl;
        cout << "Weighted SS WJets events: " << nWJetsSS_weighted << endl;

        f->cd();
        cout << "\tWriting into file...";

        h_emu_mass_before_PUCorr->Write();
        h_emu_mass_before_EffCorr->Write();
        h_emu_mass_before_PVzCorr->Write();
        h_emu_mass_before_L1Corr->Write();
        h_emu_mass_fine->Write();
        h_emu_mass->Write();
        h_emuSS_mass_before_PUCorr->Write();
        h_emuSS_mass_before_EffCorr->Write();
        h_emuSS_mass_before_PVzCorr->Write();
        h_emuSS_mass_before_L1Corr->Write();
        h_emuSS_mass_fine->Write();
        h_emuSS_mass->Write();

        h_nVTX_before_PUCorr->Write();
        h_nVTX_before_EffCorr->Write();
        h_nVTX->Write();

        h_ele_pT_before_PUCorr->Write();
        h_ele_eta_before_PUCorr->Write();
        h_ele_phi_before_PUCorr->Write();
        h_eleSS_pT_before_PUCorr->Write();
        h_eleSS_eta_before_PUCorr->Write();
        h_eleSS_phi_before_PUCorr->Write();

        h_ele_pT_before_EffCorr->Write();
        h_ele_eta_before_EffCorr->Write();
        h_ele_phi_before_EffCorr->Write();
        h_eleSS_pT_before_EffCorr->Write();
        h_eleSS_eta_before_EffCorr->Write();
        h_eleSS_phi_before_EffCorr->Write();

        h_ele_pT_before_PVzCorr->Write();
        h_ele_eta_before_PVzCorr->Write();
        h_ele_phi_before_PVzCorr->Write();
        h_eleSS_pT_before_PVzCorr->Write();
        h_eleSS_eta_before_PVzCorr->Write();
        h_eleSS_phi_before_PVzCorr->Write();

        h_ele_pT_before_L1Corr->Write();
        h_ele_eta_before_L1Corr->Write();
        h_ele_phi_before_L1Corr->Write();
        h_eleSS_pT_before_L1Corr->Write();
        h_eleSS_eta_before_L1Corr->Write();
        h_eleSS_phi_before_L1Corr->Write();

        h_ele_pT->Write();
        h_ele_eta->Write();
        h_ele_phi->Write();
        h_eleSS_pT->Write();
        h_eleSS_eta->Write();
        h_eleSS_phi->Write();

        h_mu_pT_before_PUCorr->Write();
        h_mu_eta_before_PUCorr->Write();
        h_mu_phi_before_PUCorr->Write();
        h_muSS_pT_before_PUCorr->Write();
        h_muSS_eta_before_PUCorr->Write();
        h_muSS_phi_before_PUCorr->Write();

        h_mu_pT_before_EffCorr->Write();
        h_mu_eta_before_EffCorr->Write();
        h_mu_phi_before_EffCorr->Write();
        h_muSS_pT_before_EffCorr->Write();
        h_muSS_eta_before_EffCorr->Write();
        h_muSS_phi_before_EffCorr->Write();

        h_mu_pT_before_PVzCorr->Write();
        h_mu_eta_before_PVzCorr->Write();
        h_mu_phi_before_PVzCorr->Write();
        h_muSS_pT_before_PVzCorr->Write();
        h_muSS_eta_before_PVzCorr->Write();
        h_muSS_phi_before_PVzCorr->Write();

        h_mu_pT_before_L1Corr->Write();
        h_mu_eta_before_L1Corr->Write();
        h_mu_phi_before_L1Corr->Write();
        h_muSS_pT_before_L1Corr->Write();
        h_muSS_eta_before_L1Corr->Write();
        h_muSS_phi_before_L1Corr->Write();

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
        if ( !f->IsOpen() ) cout << "File Hist_" << Mgr.Procname[_EMu_Bkg_Full]+RocCor+debug << ".root has been closed successfully.\n" << endl;
        else cout << "FILE Hist_" << Mgr.Procname[_EMu_Bkg_Full]+RocCor+debug << ".root COULD NOT BE CLOSED!\n" << endl;
    }
    else
    {
        if ( !f->IsOpen() ) cout << "File Hist_" << Mgr.Procname[Mgr.CurrentProc]+RocCor+debug << ".root has been closed successfully.\n" << endl;
        else cout << "FILE Hist_" << Mgr.Procname[Mgr.CurrentProc]+RocCor+debug << ".root COULD NOT BE CLOSED!\n" << endl;
    }

    Double_t TotalRunTime = totaltime.CpuTime();
    cout << "Total RunTime: " << TotalRunTime << " seconds" << endl;

    TTimeStamp ts_end;
    cout << "[End Time(local time): " << ts_end.AsString("l") << "]" << endl;

} // End of EMu_HistMaker()


// ######################## NOT NEEDED ############################## //
/// ----------------------- MuMu_merge ----------------------------- ///
void MuMu_merge()
{
    LocalFileMgr Mgr;

    TH1D *h_mass_before_PUCorr[7], *h_mass_before_RoccoR[7], *h_mass_before_EffCorr[7], *h_mass_fine[7], *h_mass[7], *h_Pt[7], *h_rapi[7],
         *h_nPU_before_PUCorr[7], *h_nPU_before_EffCorr[7], *h_nPU[7], *h_nVTX_before_PUCorr[7], *h_nVTX_before_EffCorr[7], *h_nVTX[7], *h_pT[7], *h_eta[7], *h_phi[7];

    TFile* files[7];
    Mgr.SetProc(_MuMu_SingleMuon_Full);

    TH1D *h_mass_before_PUCorr_full  = new TH1D( "h_mass_before_PUCorr_"+Mgr.Procname[Mgr.CurrentProc],
                                                      "h_mass_before_PUCorr_"+Mgr.Procname[Mgr.CurrentProc], 10000, 0, 10000 );
    TH1D *h_mass_before_RoccoR_full  = new TH1D( "h_mass_before_RoccoR_"+Mgr.Procname[Mgr.CurrentProc],
                                                      "h_mass_before_RoccoR_"+Mgr.Procname[Mgr.CurrentProc], 10000, 0, 10000 );
    TH1D *h_mass_before_EffCorr_full = new TH1D( "h_mass_before_EffCorr_"+Mgr.Procname[Mgr.CurrentProc],
                                                      "h_mass_before_EffCorr_"+Mgr.Procname[Mgr.CurrentProc], 10000, 0, 10000 );
    TH1D *h_mass_fine_full                = new TH1D( "h_mass_fine_"+Mgr.Procname[Mgr.CurrentProc],
                                                      "h_mass_fine_"+Mgr.Procname[Mgr.CurrentProc], 10000, 0, 10000 );
    TH1D *h_mass_full                     = new TH1D( "h_mass_"+Mgr.Procname[Mgr.CurrentProc],
                                                      "h_mass_"+Mgr.Procname[Mgr.CurrentProc], 43, massbins );
    TH1D *h_Pt_full                       = new TH1D( "h_Pt_"+Mgr.Procname[Mgr.CurrentProc],
                                                      "h_Pt_"+Mgr.Procname[Mgr.CurrentProc], 300, 0, 600 );
    TH1D *h_rapi_full                     = new TH1D( "h_rapi_"+Mgr.Procname[Mgr.CurrentProc],
                                                      "h_rapi_"+Mgr.Procname[Mgr.CurrentProc], 100, -5, 5 );
    TH1D *h_nPU_before_PUCorr_full        = new TH1D( "h_nPU_before_PUCorr_"+Mgr.Procname[Mgr.CurrentProc],
                                                      "h_nPU_before_PUCorr_"+Mgr.Procname[Mgr.CurrentProc], 75, 0, 75 );
    TH1D *h_nPU_before_EffCorr_full       = new TH1D( "h_nPU_before_EffCorr_"+Mgr.Procname[Mgr.CurrentProc],
                                                      "h_nPU_before_EffCorr_"+Mgr.Procname[Mgr.CurrentProc], 75, 0, 75 );
    TH1D *h_nPU_full                      = new TH1D( "h_nPU_"+Mgr.Procname[Mgr.CurrentProc],
                                                      "h_nPU_"+Mgr.Procname[Mgr.CurrentProc], 75, 0, 75 );
    TH1D *h_nVTX_before_PUCorr_full       = new TH1D( "h_nVTX_before_PUCorr_"+Mgr.Procname[Mgr.CurrentProc],
                                                      "h_nVTX_before_PUCorr_"+Mgr.Procname[Mgr.CurrentProc], 75, 0, 75 );
    TH1D *h_nVTX_before_EffCorr_full      = new TH1D( "h_nVTX_before_EffCorr_"+Mgr.Procname[Mgr.CurrentProc],
                                                      "h_nVTX_before_EffCorr_"+Mgr.Procname[Mgr.CurrentProc], 75, 0, 75 );
    TH1D *h_nVTX_full                     = new TH1D( "h_nVTX_"+Mgr.Procname[Mgr.CurrentProc],
                                                      "h_nVTX_"+Mgr.Procname[Mgr.CurrentProc], 75, 0, 75 );
    TH1D *h_pT_full                       = new TH1D( "h_pT_"+Mgr.Procname[Mgr.CurrentProc],
                                                      "h_pT_"+Mgr.Procname[Mgr.CurrentProc], 300, 0, 600 );
    TH1D *h_eta_full                      = new TH1D( "h_eta_"+Mgr.Procname[Mgr.CurrentProc],
                                                      "h_eta_"+Mgr.Procname[Mgr.CurrentProc], 100, -5, 5 );
    TH1D *h_phi_full                      = new TH1D( "h_phi_"+Mgr.Procname[Mgr.CurrentProc],
                                                      "h_phi_"+Mgr.Procname[Mgr.CurrentProc], 100, -5, 5 );

    cout << "Merging histograms...";

    TFile* newFile = new TFile( Mgr.HistLocation+"Hist_"+Mgr.Procname[Mgr.CurrentProc]+".root", "RECREATE" );
    Int_t iter = 0;
    for ( SelProc_t pr=_MuMu_SingleMuon_B; pr<=_MuMu_SingleMuon_H; pr=next(pr) )
    {
        Mgr.SetProc(pr);
        files[iter] = new TFile( Mgr.HistLocation+"Hist_"+Mgr.Procname[Mgr.CurrentProc]+".root", "READ" );

        files[iter]->GetObject( "h_mass_before_PUCorr_"+Mgr.Procname[Mgr.CurrentProc],  h_mass_before_PUCorr[iter] );
        files[iter]->GetObject( "h_mass_before_RoccoR_"+Mgr.Procname[Mgr.CurrentProc],  h_mass_before_RoccoR[iter] );
        files[iter]->GetObject( "h_mass_before_EffCorr_"+Mgr.Procname[Mgr.CurrentProc], h_mass_before_EffCorr[iter] );
        files[iter]->GetObject( "h_mass_fine_"+Mgr.Procname[Mgr.CurrentProc],                h_mass_fine[iter] );
        files[iter]->GetObject( "h_mass_"+Mgr.Procname[Mgr.CurrentProc],                     h_mass[iter] );
        files[iter]->GetObject( "h_Pt_"+Mgr.Procname[Mgr.CurrentProc],                       h_Pt[iter] );
        files[iter]->GetObject( "h_rapi_"+Mgr.Procname[Mgr.CurrentProc],                     h_rapi[iter] );
        files[iter]->GetObject( "h_nPU_before_PUCorr_"+Mgr.Procname[Mgr.CurrentProc],        h_nPU_before_PUCorr[iter] );
        files[iter]->GetObject( "h_nPU_before_EffCorr_"+Mgr.Procname[Mgr.CurrentProc],       h_nPU_before_EffCorr[iter] );
        files[iter]->GetObject( "h_nPU_"+Mgr.Procname[Mgr.CurrentProc],                      h_nPU[iter] );
        files[iter]->GetObject( "h_nVTX_before_PUCorr_"+Mgr.Procname[Mgr.CurrentProc],       h_nVTX_before_PUCorr[iter] );
        files[iter]->GetObject( "h_nVTX_before_EffCorr_"+Mgr.Procname[Mgr.CurrentProc],      h_nVTX_before_EffCorr[iter] );
        files[iter]->GetObject( "h_nVTX_"+Mgr.Procname[Mgr.CurrentProc],                     h_nVTX[iter] );
        files[iter]->GetObject( "h_pT_"+Mgr.Procname[Mgr.CurrentProc],                       h_pT[iter] );
        files[iter]->GetObject( "h_eta_"+Mgr.Procname[Mgr.CurrentProc],                      h_eta[iter] );
        files[iter]->GetObject( "h_phi_"+Mgr.Procname[Mgr.CurrentProc],                      h_phi[iter] );

        h_mass_before_PUCorr[iter] ->SetDirectory(0);
        h_mass_before_RoccoR[iter] ->SetDirectory(0);
        h_mass_before_EffCorr[iter]->SetDirectory(0);
        h_mass_fine[iter]               ->SetDirectory(0);
        h_mass[iter]                    ->SetDirectory(0);
        h_Pt[iter]                      ->SetDirectory(0);
        h_rapi[iter]                    ->SetDirectory(0);
        h_nPU_before_PUCorr[iter]       ->SetDirectory(0);
        h_nPU_before_EffCorr[iter]      ->SetDirectory(0);
        h_nPU[iter]                     ->SetDirectory(0);
        h_nVTX_before_PUCorr[iter]      ->SetDirectory(0);
        h_nVTX_before_EffCorr[iter]     ->SetDirectory(0);
        h_nVTX[iter]                    ->SetDirectory(0);
        h_pT[iter]                      ->SetDirectory(0);
        h_eta[iter]                     ->SetDirectory(0);
        h_phi[iter]                     ->SetDirectory(0);

        h_mass_before_PUCorr_full ->Add( h_mass_before_PUCorr[iter] );
        h_mass_before_RoccoR_full ->Add( h_mass_before_RoccoR[iter] );
        h_mass_before_EffCorr_full->Add( h_mass_before_EffCorr[iter] );
        h_mass_fine_full               ->Add( h_mass_fine[iter] );
        h_mass_full                    ->Add( h_mass[iter] );
        h_Pt_full                      ->Add( h_Pt[iter] );
        h_rapi_full                    ->Add( h_rapi[iter] );
        h_nPU_before_PUCorr_full       ->Add( h_nPU_before_PUCorr[iter] );
        h_nPU_before_EffCorr_full      ->Add( h_nPU_before_EffCorr[iter] );
        h_nPU_full                     ->Add( h_nPU[iter] );
        h_nVTX_before_PUCorr_full      ->Add( h_nVTX_before_PUCorr[iter] );
        h_nVTX_before_EffCorr_full     ->Add( h_nVTX_before_EffCorr[iter] );
        h_nVTX_full                    ->Add( h_nVTX[iter] );
        h_pT_full                      ->Add( h_pT[iter] );
        h_eta_full                     ->Add( h_eta[iter] );
        h_phi_full                     ->Add( h_phi[iter] );

        iter++;
    }// End of for(SingleMuon_BtoH)

    newFile->cd();
    h_mass_before_PUCorr_full  ->Write();
    h_mass_before_RoccoR_full  ->Write();
    h_mass_before_EffCorr_full ->Write();
    h_mass_fine_full                ->Write();
    h_mass_full                     ->Write();
    h_Pt_full                       ->Write();
    h_rapi_full                     ->Write();
    h_nPU_before_PUCorr_full        ->Write();
    h_nPU_before_EffCorr_full       ->Write();
    h_nPU_full                      ->Write();
    h_nVTX_before_PUCorr_full       ->Write();
    h_nVTX_before_EffCorr_full      ->Write();
    h_nVTX_full                     ->Write();
    h_pT_full                       ->Write();
    h_eta_full                      ->Write();
    h_phi_full                      ->Write();
    newFile->Close();

    cout << "Finished." << endl;

}// End of MuMu_Merge
