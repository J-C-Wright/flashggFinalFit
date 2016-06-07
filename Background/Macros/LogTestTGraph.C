#include <iostream>
#include <vector>
#include <iomanip>
#include <TMath>

using namespace std;

void LogTestTGraph() {

    TFile *_file0 = TFile::Open("/afs/cern.ch/user/m/musella/public/workspace/exo/full_analysis_spring15_7415v2_sync_v5_data_ecorr_cic2_final_ws.root");
    t_EBEB=_file0->Get("tree_data_cic2_EBEB");
    TTree *tree = (TTree*)_file0->Get("tree_data_cic2_EBEB");
    gStyle->SetOptStat(11111);
    gStyle->SetOptFit(11111);

    t_EBEB->Draw("mgg>>h_EBEB(34,230,910)");
    TH1F *h_EBEB = (TH1F*)gPad->GetPrimitive("h_EBEB");
    h_EBEB->SetMarkerColor(kBlack);
    h_EBEB->SetStats(kFALSE);

    double xMin(230.0);
    double xMax(910.0);
    const Int_t nBins(34);

    //logData
    //Bin limit test
    Float_t logbins[35];
    for (unsigned i=0;i<=nBins;i++){
        logbins[i] = log(xMin + i*(xMax-xMin)/(double)nBins);
    }
    TH1F* loglogHist = new TH1F("logMassHist","logMassHist",34,logbins);
    for (unsigned i=0;i<tree->GetEntries();i++) {
        tree->GetEntry(i);
        loglogHist->Fill(log(tree->GetBranch("mgg")->GetLeaf("mgg")->GetValue()));
    }

    //Polynomial in loglog
    TF1* poly1 = new TF1("Poly1","[0] + [1]*x",log(xMin),log(xMax));
    poly1->SetLineColor(kBlue);
    TF1* poly1Alt = new TF1("Poly1Alt","[0] + [1]*x",log(xMin),log(xMax));
    poly1Alt->SetLineColor(kRed);

    //Normal dijet
    TF1* dijet1 = new TF1("Dijet1","TMath::Exp([0] + [1]*TMath::Log(x))",xMin,xMax);
    dijet1->SetLineColor(kRed);
    TF1* dijet1Alt = new TF1("Dijet1Alt","TMath::Exp([0] + [1]*TMath::Log(x))",xMin,xMax);
    dijet1Alt->SetLineColor(kBlue);

    Float_t binContent[34];
    Float_t binCentre[34];
    Float_t binXErrTop[34];
    Float_t binXErrBottom[34];
    Float_t binYErrTop[34];
    Float_t binYErrBottom[34];

    for(unsigned i=1;i<=loglogHist->GetNbinsX();i++){

        if (loglogHist->GetBinContent(i) != 0.0) {
            binContent[i-1] = log( loglogHist->GetBinContent(i) );
            binCentre[i-1] = loglogHist->GetBinCenter(i);
            binYErrTop[i-1] = log( loglogHist->GetBinContent(i) + loglogHist->GetBinError(i) ) - binContent[i-1];
            if ( loglogHist->GetBinContent(i) != loglogHist->GetBinError(i)) {
                binYErrBottom[i-1] = binContent[i-1] - log( loglogHist->GetBinContent(i) - loglogHist->GetBinError(i) );
            }else{
                binYErrBottom[i-1] = binContent[i-1];
            }
        }else{
            binContent[i-1] = 0;
            binCentre[i-1] = loglogHist->GetBinCenter(i); 
            binYErrTop[i-1] = 0;
            binYErrBottom[i-1] = 0;
        }
        binXErrTop[i-1] = fabs(binCentre[i-1] - logbins[i]);
        binXErrBottom[i-1] = fabs(binCentre[i-1] - logbins[i-1]);

        bool setXErrsZero = true;
        if (setXErrsZero){
            binXErrTop[i-1] = 0;
            binXErrBottom[i-1] = 0;
        }
    }

    TGraphAsymmErrors *logGraph = new TGraphAsymmErrors(34,binCentre,binContent,binXErrBottom,binXErrTop,binYErrBottom,binYErrTop);
    TCanvas c1("c1");
    logGraph->Draw("AP");
    c1.Print("Plots/TGraphErrsTest.pdf");

    logGraph->Fit("Poly1","RM");
    h_EBEB->Fit("Dijet1","RM");

    poly1Alt->SetParameter(0,dijet1->GetParameter(0));
    poly1Alt->SetParameter(1,dijet1->GetParameter(1));
    dijet1Alt->SetParameter(0,poly1->GetParameter(0));
    dijet1Alt->SetParameter(1,poly1->GetParameter(1));

    c1.Clear();
    logGraph->Draw("AP");
    poly1->Draw("same");
    poly1Alt->Draw("same");
    c1.Print("Plots/TGraphErrsTestPoly.pdf");

    c1.Clear();
    h_EBEB->Draw("EP");
    dijet1->Draw("same");
    dijet1Alt->Draw("same");
    c1.Print("Plots/TGraphErrsTestDijet.pdf");

}

