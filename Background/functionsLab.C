#include <vector>
void functionsLab() {

    TFile *_file0 = TFile::Open("/afs/cern.ch/user/m/musella/public/workspace/exo/full_analysis_spring15_7415v2_sync_v5_data_ecorr_cic2_final_ws.root");
    t_EBEE=_file0->Get("tree_data_cic2_EBEE");
    t_EBEB=_file0->Get("tree_data_cic2_EBEB");
    //TCanvas *t = new TCanvas ("functionLab","functionLab",800,400);
    //t->Divide(2,1);
    gStyle->SetOptStat(11111);
    gStyle->SetOptFit(11111);

    //t->cd(1);
    t_EBEE->Draw("mgg>>h_EBEE(100,320,1600)");
    TH1F *h_EBEE = (TH1F*)gPad->GetPrimitive("h_EBEE");
    h_EBEE->SetMarkerColor(kBlack);
    h_EBEE->SetMarkerSize(4);
    h_EBEE->Draw("EP");

    //t->cd(2);
    t_EBEB->Draw("mgg>>h_EBEB(100,230,1600)");
    TH1F *h_EBEB = (TH1F*)gPad->GetPrimitive("h_EBEB");
    h_EBEB->SetMarkerColor(kBlack);
    h_EBEB->Draw("EP");




//-------------------
//Functions and fits
//-------------------
    TCanvas c1("c1");
    h_EBEB->SetStats(kFALSE);

//EBEB Region

//Polynomial fits

    unsigned maxPolyDegree(10);

    //Power basis
    std::vector<TF1*> powerBasis(maxPolyDegree);
    for (unsigned i(0);i<maxPolyDegree;i++) {
        powerBasis[i] = new TF1(Form("polyPB%d",i+1),Form("pol%d",i+1),0,1600);
        powerBasis[i]->SetLineColor(i+1);

        h_EBEB->Fit(Form("polyPB%d",i+1),"RM");
    }

    h_EBEB->Draw("EP");
    for (unsigned i(0);i<maxPolyDegree;i++) {
        powerBasis[i]->Draw("SAME");
    }
    c1.Print("FitPlots/Poly_PowerBasis.pdf");

    //Bernstein basis
    std::vector<TF1*> bernsteinBasis(maxPolyDegree);
    h_EBEB->Draw("EP");
    for (unsigned i(0);i<maxPolyDegree;i++) {

        std::cout << "\nMaking and fitting Bernstein basis polynomial of degree " << i+1 << " function string is:" << std::endl;
        TString functionString = BernsteinString(i+1);
        std::cout << functionString << std::endl;

        bernsteinBasis[i] = new TF1(Form("Bern_D%d",i+1),functionString,0,1600);
        bernsteinBasis[i]->SetLineColor(i+1);

        h_EBEB->Fit(Form("Bern_D%d",i+1),"RM");
    }

    h_EBEB->Draw("EP");
    for (unsigned i(0);i<maxPolyDegree;i++) {
        bernsteinBasis[i]->Draw("SAME");
    }
    c1.Print("FitPlots/Poly_BernsteinBasis.pdf");

//Power Laws

//Exponentials

//Laurent


    
}

TString BernsteinString(unsigned n) {

    if (n == 0) {return TString("[0]");}

    std::ostringstream poly;

    for (unsigned i(0);i<n+1;i++) {

        unsigned factorialPart = (unsigned)TMath::Factorial(n)/( TMath::Factorial(i)*TMath::Factorial(n-i) );

        poly << "[" << i << "]*";
        if (factorialPart > 1) {poly << factorialPart << ".0*";}
        if (i > 1) {poly << "pow(x," << i << ")";} else if (i == 1) {poly << "x";}
        if (n != i && i != 0) {poly << "*";}
        if (n-i > 1) {poly << "pow(1-x," << n-i << ")";} else if (n-i == 1) {poly << "(1-x)";}
        if (i < n) {poly << " + ";}
    }
    TString output = poly.str();

    return output;
};




