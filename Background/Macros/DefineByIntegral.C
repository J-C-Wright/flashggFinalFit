#include <vector>
#include<iostream>

using namespace std;

void DefineByIntegral() {

    TFile *output = new TFile("FitPlots/EBEB_BackgroundFits.root","RECREATE");

    TFile *_file0 = TFile::Open("/afs/cern.ch/user/m/musella/public/workspace/exo/full_analysis_spring15_7415v2_sync_v5_data_ecorr_cic2_final_ws.root");
    t_EBEB=_file0->Get("tree_data_cic2_EBEB");
    gStyle->SetOptStat(11111);
    gStyle->SetOptFit(11111);

    t_EBEB->Draw("mgg>>h_EBEB(34,230,910)");
    TH1F *h_EBEB = (TH1F*)gPad->GetPrimitive("h_EBEB");
    h_EBEB->SetMarkerColor(kBlack);

//-------------------
//Functions and fits
//-------------------
    TCanvas c1("c1");
    c1.SetLogy();
    c1.SetLogx();

    h_EBEB->SetStats(kFALSE);

    unsigned order(5);
    vector<TF1*> dijets(order);
    vector<TF1*> dijetAlts(order);

    for (unsigned i=1;i<order;i++) {

        cout << "-------" << endl;
        cout << DijetAltString(i) << endl;
        cout << DijetString(i) << endl;
        cout << "-------" << endl;

        dijets[i] = new TF1(Form("Dijet%d",i),DijetString(i),230,910);
        dijets[i]->SetLineColor(i);
        dijetAlts[i] = new TF1(Form("DijetAlt%d",i),DijetAltString(i),230,910);
        dijetAlts[i]->SetLineColor(i);
        dijetAlts[i]->SetLineStyle(2);

        h_EBEB->Fit(Form("Dijet%d",i),"QRMI0");

        cout << dijets[i]->GetNpar() << endl;
        for (unsigned j=0;j<dijets[i]->GetNpar();j++){
            cout << setw(16) << dijets[i]->GetParameter(j);
            //dijetAlts[i]->SetParameter(j,dijets[i]->GetParameter(j));
        }
        cout << endl;

        h_EBEB->Fit(Form("DijetAlt%d",i),"RMI0");
    }

    c1.Clear();
    h_EBEB->Draw("EP");
    h_EBEB->GetXaxis()->SetTitle("m_{#gamma#gamma} (GeV)");
    for (unsigned k(1);k<order;k++){
        dijets[k]->Draw("same");
        dijetAlts[k]->Draw("same");
    }
    c1.Print("DijetAllOrders.pdf");

}


TString DijetAltString(unsigned k) {

    TString form;
    form += "(";
    for (unsigned i=1;i<=k;i++) {
        if (i==1){
            form += "[1]/x";
        }else if (i == 2) {
            form += " + 2*[2]*TMath::Log(x)/x";
        }else if (i > 2) {
            form += Form(" + %d*[%d]*pow(TMath::Log(x),%d)/x",i,i,i-1);
        }
    }
    form += ")*TMath::Exp(";
    for (unsigned i=0;i<=k;i++) {
        if (i == 0) {
            form += "[0]";
        }else if (i == 1) {
            form += " + [1]*TMath::Log(x)";
        }else if (i > 1) {
            form += Form(" + [%d]*pow(TMath::Log(x),%d)",i,i);
        }
    }
    form += ")";
    return form;
}

TString DijetString(unsigned k) {
    std::ostringstream dijet;
    dijet << "TMath::Exp(";
    if (k == 0) {
        dijet << " [0]";
    }else{
        dijet << " [0]";
        for (unsigned i(1);i<k+1;i++) {
            dijet << " + [" << i << "]*pow(TMath::Log(x)," << i << ")";
        }
    }
    dijet << " )";
    TString output = dijet.str();
    return output;
}



