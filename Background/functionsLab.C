#include <vector>
void functionsLab() {

    TFile *output = new TFile("FitPlots/BackgroundFits.root","RECREATE");

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

    output->cd();
    h_EBEE->Write();
    h_EBEB->Write();
    output->Close();



//-------------------
//Functions and fits
//-------------------
    TCanvas c1("c1");
    h_EBEB->SetStats(kFALSE);

//EBEB Region

//Polynomial fits

    

    //Power basis
    std::vector<TF1*> powerBasis(2);
    for (unsigned i(0);i<3;i++) {
        powerBasis[i] = new TF1(Form("polyPB%d",i+1),Form("pol%d",i+1),230,910);
        powerBasis[i]->SetLineColor(i+1);

        h_EBEB->Fit(Form("polyPB%d",i+1),"RM");
    }

    h_EBEB->Draw("EP");
    for (unsigned i(0);i<3;i++) {
        powerBasis[i]->Draw("SAME");
    }
    h_EBEB->Draw("SAME");
    c1.Print("FitPlots/Polys_PowerBasis.pdf");

    //Bernstein basis
    std::vector<TF1*> bernsteinBasis(4);
    h_EBEB->Draw("EP");
    for (unsigned i(0);i<5;i++) {

        std::cout << "\nMaking and fitting Bernstein basis polynomial of degree " << i+1 << " function string is:" << std::endl;
        TString functionString = SubbasisBernsteinString(i+1,230,910);
        std::cout << functionString << std::endl;

        bernsteinBasis[i] = new TF1(Form("Bern_D%d",i+1),functionString,230,910);
        bernsteinBasis[i]->SetLineColor(i+1);

        h_EBEB->Fit(Form("Bern_D%d",i+1),"RM");
    }

    h_EBEB->Draw("EP");
    for (unsigned i(0);i<5;i++) {
        bernsteinBasis[i]->Draw("SAME");
    }
    h_EBEB->Draw("SAME");
    c1.Print("FitPlots/Polys_BernsteinBasis.pdf");

    //(1-t)^n functions
    std::vector<TF1*> mdBasis(3);
    h_EBEB->Draw("EP");
    for (unsigned i(0);i<3;i++) {
        std::cout << MonoDecreaseString(i,230,910) << std::endl << std::endl;
        mdBasis[i] = new TF1(Form("mdBasis%d",i),MonoDecreaseString(i,230,910),230,910);
        mdBasis[i]->SetLineColor(i+1);
        h_EBEB->Fit(Form("mdBasis%d",i),"RM");
        mdBasis[i]->Draw("SAME");
    }
    h_EBEB->Draw("SAME");
    c1.Print("FitPlots/Polys_Convex.pdf");

//Power Laws

    TF1 *powerlaw = new TF1("power","[0]*pow(x,[1])",230,910);
    h_EBEB->Fit("power","RM");
    h_EBEB->Draw("EP");
    powerlaw->Draw("SAME");
    h_EBEB->Draw("SAME");
    c1.Print("FitPlots/PowerLaws.pdf");

//Exponentials

    TF1 *exponent = new TF1("exponent","exp([0] + x*[1])",230,910);
    h_EBEB->Fit("exponent","RM");
    h_EBEB->Draw("EP");
    exponent->Draw("SAME");
    h_EBEB->Draw("SAME");
    c1.Print("FitPlots/Exponentials.pdf");

//Laurent



//From outside
    TF1 *other = new TF1("other","x*( [0] + [1]*TMath::Log(x) )",230,910);
    h_EBEB->Fit("other","RM");
    h_EBEB->Draw("EP");
    other->Draw("SAME");
    c1.Print("FitPlots/Other.pdf");
    
}

TString BernsteinString(unsigned n, float xMin, float xMax) {

    if (n == 0) {return TString("[0]");}

    std::ostringstream poly;

    for (unsigned i(0);i<n+1;i++) {

        unsigned factorialPart = (unsigned)TMath::Factorial(n)/( TMath::Factorial(i)*TMath::Factorial(n-i) );

        poly << "[" << i << "]*";
        if (factorialPart > 1) {poly << factorialPart << ".0*";}
        if (i > 1) {poly << "pow((x - " << xMin << ")/" << xMax << "," << i << ")";} else if (i == 1) {poly << "((x-" << xMin << ")/" << xMax << ")";}
        if (n != i && i != 0) {poly << "*";}
        if (n-i > 1) {poly << "pow(1-(x-" << xMin << ")/" << xMax << "," << n-i << ")";} else if (n-i == 1) {poly << "(1-(x-" << xMin << ")/" << xMax  << ")";}
        if (i < n) {poly << " + ";}
    }
    TString output = poly.str();

    return output;
};

TString SubbasisBernsteinString(unsigned n, float xMin, float xMax) {

    if (n == 0) {return TString("[0]");}
    
    std::ostringstream poly;

    unsigned iLimit = n/2;

    for (unsigned i(0);i<iLimit+1;i++) {

        unsigned factorialPart = (unsigned)TMath::Factorial(n)/( TMath::Factorial(i)*TMath::Factorial(n-i) );

        poly << "[" << i << "]*";
        if (factorialPart > 1) {poly << factorialPart << ".0*";}
        if (i > 1) {poly << "pow((x - " << xMin << ")/" << xMax << "," << i << ")";} else if (i == 1) {poly << "((x-" << xMin << ")/" << xMax << ")";}
        if (n != i && i != 0) {poly << "*";}
        if (n-i > 1) {poly << "pow(1-(x-" << xMin << ")/" << xMax << "," << n-i << ")";} else if (n-i == 1) {poly << "(1-(x-" << xMin << ")/" << xMax  << ")";}
        if (i < iLimit) {poly << " + ";}
    }
    TString output = poly.str();

    return output;
};

TString MonoDecreaseString(unsigned n, float xMin, float xMax) {

    std::ostringstream poly;
    for (unsigned i(0);i<n+1;i++) {
        poly << "[" << i << "]*pow(1-(x-" << xMin << ")/" << xMax << "," << i << ")";
        if (i != n) poly << " + ";
    }
    TString output = poly.str();

    return output;
};
    










