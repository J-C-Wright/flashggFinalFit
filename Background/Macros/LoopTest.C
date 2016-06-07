#include <iostream>
#include <vector>
#include <iomanip>
#include <TMath>

using namespace std;


void LoopTest() {

    TFile *_file0 = TFile::Open("/afs/cern.ch/user/m/musella/public/workspace/exo/full_analysis_spring15_7415v2_sync_v5_data_ecorr_cic2_final_ws.root");
    t_EBEB=_file0->Get("tree_data_cic2_EBEB");
    gStyle->SetOptStat(11111);
    gStyle->SetOptFit(11111);

    t_EBEB->Draw("mgg>>h_EBEB(34,230,910)");
    TH1F *h_EBEB = (TH1F*)gPad->GetPrimitive("h_EBEB");
    h_EBEB->SetMarkerColor(kBlack);

    TCanvas c1("c1");
    c1.SetLogy();
    c1.SetLogx();

    h_EBEB->SetStats(kFALSE);

    double xMin(230.0);
    double xMax(910.0);

    unsigned order(5);
    vector<TF1*> dijets(order);
    vector<vector<double>> p(order);

    h_EBEB->Draw();
    for (unsigned i=0;i<order;i++) {

        dijets[i] = new TF1(Form("Dijet%d",i+1),ScaledDijetString(i+1,230,910),xMin,xMax);
        cout << ScaledDijetString(i+1,230,910) << endl;
        dijets[i]->SetLineColor(i+1);
        h_EBEB->Fit(Form("Dijet%d",i+1),"RMI0");

        for (unsigned j=0;j<dijets[i]->GetNpar();j++){
            p[i].push_back(dijets[i]->GetParameter(j));
        }
        dijets[i]->Draw("same");
    }
    c1.Print("AllScaledDijets.pdf");

    /*
    unsigned test(1);
    double minI = dijetIntegral(p[test],rescale(xMin,xMin,xMax));
    double maxI = dijetIntegral(p[test],rescale(xMax,xMin,xMax));
    cout << setw(12) << minI << setw(12) << maxI << setw(24) << maxI - minI << endl;
    cout << setw(12) << dijet(p[test],rescale(xMin,xMin,xMax)) << setw(12) << dijet(p[test],rescale(xMax,xMin,xMax)) << endl;
    cout << setw(12) << numericalIntegral(p[test],rescale(xMin,xMin,xMax),rescale(xMax,xMin,xMax),1e2) << endl;
    */
}








double rescale(double x, double xMin, double xMax) {
    return (x-xMin)/(xMax-xMin)*1.7182818284 + 1; //Maps x to between 1 and e
}

double fact(unsigned &n) {
    double out(1);
    for (unsigned i(0);i<n;i++) out *= n-i;
    return out;
}

double I_N_iteration(vector<double> &p, vector<unsigned> &it, double &x) {

    bool debug(false);

    double part1 = pow(2.7182818284,p[0])*x;
    for (unsigned i(0);i<it.size();i++) {
        part1 *= pow(p[i+1],it[i])/fact(it[i]);
    }
        
    unsigned alpha(0);
    for (unsigned i(0);i<it.size();i++) {
        alpha += (i+1)*it[i];
    }

    double part2(0.0);
    for (unsigned i(0);i<=alpha;i++) {
        part2 += pow(-1,alpha-i)*fact(alpha)*pow(log(x),i)/fact(i);
        cout << "---- ---- " << part2 << setw(16) << pow(-1,alpha-i)*fact(alpha)*pow(log(x),i)/fact(i) << endl;
    }

    double out = part1*part2;

    if (debug) cout << "----" << setw(16) << part1 << setw(6) << alpha << setw(16) << part2 << setw(16) << out << endl;

    if (TMath::Finite(out)) {
        return out;
    }else{
        return 0.0;
    }

}

double dijetIntegral(vector<double> &p, double &x) {

    bool skip(false);
    bool debug(false);
    unsigned itMax(5);
    double precision(1e-7);

    vector<unsigned> it(p.size()-1,0);

    double out(0.0);
    while (it[0] < itMax) {

        double step = I_N_iteration(p,it,x);
        out += step;

        if (debug) {
            for (unsigned i(0);i<it.size();i++) {
                cout << setw(6) << it[i];
            }
            cout << setw(16) << step << setw(12) << out << endl;
        }

        unsigned zeroCount(0);
        if (fabs(step) < precision && it.back() == 0 && skip) {
            for (unsigned i(it.size()-1);i>=0;i--) {
                if (it[i] == 0) {
                    zeroCount++;
                }else{
                    break;
                }
            }
            if (zeroCount > 0) it[it.size()-1-zeroCount] = itMax;
        }

        it.back()++;
        if (fabs(step) < precision && skip) it.back() = itMax;
        if (it.back() == itMax) {
            for (unsigned i(0);i<it.size()-1;i++) {
                if (it[it.size()-1-i] >= itMax) {
                    it[it.size()-1-i] = 0;
                    it[it.size()-2-i]++;
                }
            }
        }
    }

    return out;
}

double dijet(vector<double> &p, double &x) {

    bool debug(false);
    double out(0.0);
    for (unsigned i(0);i<p.size();i++) {

        if (i == 0) {
            out += p[0];
        }else if (i == 1) {
            out += p[1]*log(x);
        }else{
            out += p[i]*pow(log(x),i);
        }
        if (debug) {
            cout << setw(16) << out;
        }

    }
    if (debug) cout << " --> " << exp(out) << endl;
    return exp(out);
}

double numericalIntegral(vector<double> &p, double &xMin, double &xMax, unsigned slices) {
    
    double step = (xMax - xMin)/(double)slices;
    double out(0.0);
    for (unsigned i(0);i<slices;i++) {
        out += dijet(p,xMin + i*step) * step;
    }
    return out;

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

TString ScaledDijetString(unsigned k, float xMin, float xMax) {
    ostringstream dijet;
    dijet << "TMath::Exp(";
    for (unsigned i=0;i<=k;i++){
        if (i==0){
            dijet << "[0]";
        }else if (i==1){
            dijet << Form("+[1]*TMath::Log((x-%f)*1.7182818284/%f+1)",xMin,xMax-xMin);
        }else if (i>1){
            dijet << Form("+[%d]*pow(TMath::Log((x-%f)*1.7182818284/%f+1), %d)",i,xMin,xMax-xMin,i);
        }
    }
    dijet << ")";
    return dijet.str();
}

