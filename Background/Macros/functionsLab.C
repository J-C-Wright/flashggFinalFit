{
TFile *_file0 = TFile::Open("/afs/cern.ch/work/j/jwright/public/Louie/full_analysis_spring15_7415v2_sync_v5_data_ecorr_cic2_final_ws.root");
t_EBEE=_file0->Get("tree_data_cic2_EBEE");
t_EBEB=_file0->Get("tree_data_cic2_EBEB");
//TCanvas *t = new TCanvas ("functionLab","functionLab",800,400);
//t->Divide(2,1);
gStyle->SetOptStat(11111);
gStyle->SetOptFit(11111);

//t->cd(1);
t_EBEE->Draw("mgg>>h_EBEE(28,320,910)");
TH1F *h_EBEE = (TH1F*)gPad->GetPrimitive("h_EBEE");
h_EBEE->SetMarkerColor(kBlack);
h_EBEE->SetMarkerSize(4);
h_EBEE->Draw("EP");

//t->cd(2);
t_EBEB->Draw("mgg>>h_EBEB(32,230,910)");
TH1F *h_EBEB = (TH1F*)gPad->GetPrimitive("h_EBEB");
h_EBEB->SetMarkerColor(kBlack);
h_EBEB->Draw("EP");

std::cout << " --------------------------------------------------------------------------- " << std::endl;
std::cout << " Usage : The two categories considered in this analysis are EBEE and EBEB.   " << std::endl;
std::cout << " The mgg dists are loaded with correct ranges as TH1Fs : h_EBEB and h_EBEE   " << std::endl;
std::cout << " You can try to fit them using CINT the usual format for ROOT fitting:       " << std::endl;
std::cout << " h_EBEE->Fit(\"expo\")                                                       " << std::endl;
std::cout << " or, define arbitrary functions, using any number of params [0],[1],[2]...:  " << std::endl;
std::cout << " TF1 *f1 = new TF1(\"f1\",\"[0]*exp(-x*[1])\", 0, 1600);                     " << std::endl;
std::cout << " h_EBEE->Fit(\"f1\")                                                         " << std::endl;
std::cout << " --------------------------------------------------------------------------- " << std::endl;

        TCanvas c1("c1");

        TF1 *newFunc = new TF1("newFunc","[0]*pow((x/13000.0),[1]*pow(1-(x/13000.0),[2]))",230,910);
        TF1 *dijet2  = new TF1("dijet2","TMath::Exp([0] + [1]*TMath::Log(x) + [2]*pow(TMath::Log(x),2))",230,910);
        TF1 *newFunc_EE = new TF1("newFunc_EE","[0]*pow((x/13000.0),[1]*pow(1-(x/13000.0),[2]))",320,910);
        TF1 *dijet2_EE  = new TF1("dijet2_EE","TMath::Exp([0] + [1]*TMath::Log(x) + [2]*pow(TMath::Log(x),2))",320,910);

        newFunc->SetLineColor(kBlue);
        dijet2->SetLineColor(kRed);
        newFunc_EE->SetLineColor(kBlue);
        dijet2_EE->SetLineColor(kRed);

        h_EBEB->Fit("newFunc","RM0");
        h_EBEB->Fit("dijet2","RM0");

        h_EBEB->Draw();
        newFunc->Draw("same");
        dijet2->Draw("same");
        c1.SetLogy();
        c1.Print("~/www/NewFunction/Comparison_EBEB.pdf");
        c1.Print("~/www/NewFunction/Comparison_EBEB.png");

        std::cout << "ChiSqPerDOF New    " << newFunc->GetChisquare()/newFunc->GetNDF() << std::endl;
        std::cout << "ChiSqPerDOF Dijet2 " << dijet2->GetChisquare()/dijet2->GetNDF() << std::endl;

        h_EBEE->Fit("newFunc_EE","RM0");
        h_EBEE->Fit("dijet2_EE","RM0");

        h_EBEE->Draw();
        newFunc_EE->Draw("same");
        dijet2_EE->Draw("same");
        c1.SetLogy();
        c1.Print("~/www/NewFunction/Comparison_EBEE.pdf");
        c1.Print("~/www/NewFunction/Comparison_EBEE.png");

        std::cout << "ChiSqPerDOF New    " << newFunc_EE->GetChisquare()/newFunc_EE->GetNDF() << std::endl;
        std::cout << "ChiSqPerDOF Dijet2 " << dijet2_EE->GetChisquare()/dijet2_EE->GetNDF() << std::endl;





}

