#include <vector>
#include<iostream>

    void functionsLab() {

//        TFile *output = new TFile("FitPlots/BackgroundFits.root","RECREATE");

        TFile *_file0 = TFile::Open("/afs/cern.ch/user/m/musella/public/workspace/exo/full_analysis_spring15_7415v2_sync_v5_data_ecorr_cic2_final_ws.root");
        t_EBEE=_file0->Get("tree_data_cic2_EBEE");
        t_EBEB=_file0->Get("tree_data_cic2_EBEB");
        //TCanvas *t = new TCanvas ("functionLab","functionLab",800,400);
        //t->Divide(2,1);
        gStyle->SetOptStat(11111);
        gStyle->SetOptFit(11111);

        //t->cd(1);
        t_EBEE->Draw("mgg>>h_EBEE(34,320,1000)");
        TH1F *h_EBEE = (TH1F*)gPad->GetPrimitive("h_EBEE");
        h_EBEE->SetMarkerColor(kBlack);
        h_EBEE->SetMarkerSize(4);
        h_EBEE->Draw("EP");

        //t->cd(2);
        t_EBEB->Draw("mgg>>h_EBEB(34,230,910)");
        TH1F *h_EBEB = (TH1F*)gPad->GetPrimitive("h_EBEB");
        h_EBEB->SetMarkerColor(kBlack);
        h_EBEB->Draw("EP");


    //-------------------
    //Functions and fits
    //-------------------
        TCanvas c1("c1");
        c1.SetLogy();

        h_EBEB->SetStats(kFALSE);

    //EBEB Region
        std::ofstream EBEB_FitsLog;
        EBEB_FitsLog.open("FitPlots/EBEB_Fits.txt");
        EBEB_FitsLog << setw(24) << "Name" << setw(12) << "ChiSq" << setw(12) << "NDOF" << setw(12) << "ChiSq/NDOF" << std::endl;

    //Polynomial fits
        //Power basis
        EBEB_FitsLog << "Polynomials - Power basis" << std::endl;
        std::vector<TF1*> powerBasis(2);
        for (unsigned i(0);i<3;i++) {
            powerBasis[i] = new TF1(Form("polyPB%d",i+1),Form("pol%d",i+1),230,910);
            powerBasis[i]->SetLineColor(i+1);

            h_EBEB->Fit(Form("polyPB%d",i+1),"QRM");
            EBEB_FitsLog << setw(24) << powerBasis[i]->GetName() << setw(12) << powerBasis[i]->GetChisquare() << setw(12) << powerBasis[i]->GetNDF();
            EBEB_FitsLog << setw(12) << powerBasis[i]->GetChisquare()/powerBasis[i]->GetNDF() << std::endl;
        }

        h_EBEB->Draw("EP");
        for (unsigned i(0);i<3;i++) {
            powerBasis[i]->Draw("SAME");
        }
        c1.Print("FitPlots/Polys_PowerBasis.pdf");

        //Bernstein basis
        EBEB_FitsLog << "Polynomials - Bernstein sub-basis (i<n/2)" << std::endl;
        std::vector<TF1*> bernsteinBasis(4);
        h_EBEB->Draw("EP");
        for (unsigned i(0);i<5;i++) {

            TString functionString = SubbasisBernsteinString(i+1,230,910);
            bernsteinBasis[i] = new TF1(Form("Bern_D%d",i+1),functionString,230,910);
            bernsteinBasis[i]->SetLineColor(i+1);

            h_EBEB->Fit(Form("Bern_D%d",i+1),"QRM");
            EBEB_FitsLog << setw(24) << bernsteinBasis[i]->GetName() << setw(12) << bernsteinBasis[i]->GetChisquare() << setw(12) << bernsteinBasis[i]->GetNDF();
            EBEB_FitsLog << setw(12) << bernsteinBasis[i]->GetChisquare()/bernsteinBasis[i]->GetNDF() << std::endl;
        }

        h_EBEB->Draw("EP");
        for (unsigned i(0);i<5;i++) {
            bernsteinBasis[i]->Draw("SAME");
        }
        c1.Print("FitPlots/Polys_BernsteinBasis.pdf");

        //Convex basis pairs
        EBEB_FitsLog << "Polynomials - Convex sub-basis" << std::endl;
        std::vector<TF1*> convexPairs(45); 
        unsigned count(0);
        for (unsigned i(0);i<10;i++) {
            for (unsigned j(0);j<10;j++) {
                if ( j >= i ) continue;
                convexPairs[count] = new TF1(Form("Convex_%d_%d",i,j),ConvexPolyPairString(i,j,230,910),230,910);
                count++;
            }
        }
        for (unsigned i(0);i<convexPairs.size();i++) {
            h_EBEB->Fit( convexPairs[i]->GetName(),"QRM");
            EBEB_FitsLog << setw(24) << convexPairs[i]->GetName() << setw(12) << convexPairs[i]->GetChisquare() << setw(12) << convexPairs[i]->GetNDF();
            EBEB_FitsLog << setw(12) << convexPairs[i]->GetChisquare()/convexPairs[i]->GetNDF() << std::endl;
        }
        h_EBEB->Draw("EP");
        for (unsigned i(0);i<convexPairs.size();i++) {
            convexPairs[i]->Draw("SAME");
        }
        h_EBEB->Draw("EPSAME");
        c1.Print("FitPlots/Polys_ConvexPairs.pdf");

    //Power Laws
        //Single power law
        EBEB_FitsLog << "Power Laws" << std::endl;
        TF1 *powerlaw = new TF1("power","[0]*pow(x,[1])",230,910);
        h_EBEB->Fit("power","QRM");
        EBEB_FitsLog << setw(24) << powerlaw->GetName() << setw(12) << powerlaw->GetChisquare() << setw(12) << powerlaw->GetNDF();
        EBEB_FitsLog << setw(12) << powerlaw->GetChisquare()/powerlaw->GetNDF() << std::endl;
        h_EBEB->Draw("EP");
        powerlaw->Draw("SAME");
        c1.Print("FitPlots/PowerLaws.pdf");

    //Exponentials
        //Single exponential
        EBEB_FitsLog << "Exponentials" << std::endl;
        TF1 *exponent = new TF1("exponent","exp([0] + x*[1])",230,910);
        h_EBEB->Fit("exponent","QRM");
        EBEB_FitsLog << setw(24) << exponent->GetName() << setw(12) << exponent->GetChisquare() << setw(12) << exponent->GetNDF();
        EBEB_FitsLog << setw(12) << exponent->GetChisquare()/exponent->GetNDF() << std::endl;
        h_EBEB->Draw("EP");
        exponent->Draw("SAME");
        c1.Print("FitPlots/Exponentials.pdf");

    //Laurent




    //ATLAS Function
        //N parameters
        EBEB_FitsLog << "ATLAS functions" << std::endl;
        std::vector<TF1*> atlasFunctions(3);
        for (unsigned k(0);k<atlasFunctions.size();k++){
            std::cout << ATLASString(k) << std::endl;
            atlasFunctions[k] = new TF1(Form("Atlas%d",k),ATLASString(k),230,910);
            h_EBEB->Fit(Form("Atlas%d",k),"RM");
            EBEB_FitsLog << setw(24) << atlasFunctions[k]->GetName() << setw(12) << atlasFunctions[k]->GetChisquare() << setw(12) << atlasFunctions[k]->GetNDF();
            EBEB_FitsLog << setw(12) << atlasFunctions[k]->GetChisquare()/atlasFunctions[k]->GetNDF() << std::endl;
        }
        h_EBEB->Draw("EP");
        for (unsigned k(0);k<atlasFunctions.size();k++) {
            atlasFunctions[k]->SetLineColor(k+1);
            atlasFunctions[k]->Draw("SAME");
        }
        c1.Print("FitPlots/ATLASFunctions.pdf");

        //ATLAS Pair search
        //

    //Dijet Function
        //N parameters
        EBEB_FitsLog << "Dijet functions" << std::endl;
        std::vector<TF1*> dijetFunctions(4);
        for (unsigned k(0);k<4;k++) {
            std::cout << DijetString(k+1) << std::endl;
            dijetFunctions[k] = new TF1(Form("Dijet%d",k+1),DijetString(k+1),230,910);
            h_EBEB->Fit(Form("Dijet%d",k+1),"RM");
            EBEB_FitsLog << setw(24) << dijetFunctions[k]->GetName() << setw(12) << dijetFunctions[k]->GetChisquare() << setw(12) << dijetFunctions[k]->GetNDF();
            EBEB_FitsLog << setw(12) << dijetFunctions[k]->GetChisquare()/dijetFunctions[k]->GetNDF() << std::endl;
        }
        h_EBEB->Draw("EP");
        for (unsigned k(0);k<dijetFunctions.size();k++) {
            dijetFunctions[k]->SetLineColor(k+1);
            dijetFunctions[k]->Draw("SAME");
        }
        c1.Print("FitPlots/DijetFunctions.pdf");

        //Dijet Pair search
        //
       





        std::cout << "testing out exponent sum" << std::endl;
        for (unsigned i(0);i<4;i++) {
            std::cout << exponentialSumString(i) << std::endl;
        }
        std::cout << "testing out exponentiated polynomial" << std::endl;
        for (unsigned i(0);i<4;i++) {
            std::cout << exponentialPoly(i) << std::endl;
        }
        std::cout << "testing out exponentiated power laws" << std::endl;
        for (unsigned i(0);i<4;i++) {
            std::cout << exponentialPowerLaw(i) << std::endl;
        }
        std::cout << "Testing out power law sum" << std::endl;
        for (unsigned i(0);i<4;i++) {
            std::cout << powerLawSum(i) << std::endl;
        }
        std::cout << "Testing out power law of polynomial" << std::endl;
        for (unsigned i(2);i<6;i++) {
            std::cout << powerLawPoly(i) << std::endl;
        }
        std::cout << "Testing out power law of polynomial in the power part" << std::endl;
        for (unsigned i(0);i<6;i++) {
            std::cout << polyPowerLaw(i) << std::endl;
        }
        std::cout << "Testing out power law of convex basis polynomial" << std::endl;
        std::cout << powerLawConvexBasis(9,1,230,910) << std::endl;
        std::cout << powerLawConvexBasis(3,6,230,910) << std::endl;
        std::cout << "Testing out Laurent series function" << std::endl;
        for (unsigned i(0);i<6;i++) {
            std::cout << laurentString(i) << std::endl;
        }
        std::cout << "Testing out Laurent series pairs" << std::endl;
        std::cout << laurentPairString(4,5) << std::endl;
        std::cout << laurentPairString(3,6) << std::endl;
        std::cout << "Testing out ATLAS function pairs" << std::endl;
        std::cout << ATLASPair(4,5) << std::endl;
        std::cout << ATLASPair(1,2) << std::endl;
        std::cout << ATLASPair(0,1) << std::endl;
        std::cout << "Testing out Dijet function pairs" << std::endl;
        std::cout << DijetPair(4,5) << std::endl;
        std::cout << DijetPair(1,2) << std::endl;
        std::cout << DijetPair(0,1) << std::endl;


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

    TString ConvexPolyString(unsigned start, unsigned end, float xMin, float xMax) {
        std::ostringstream poly;
        for (unsigned i(start);i<end+1;i++) {
            poly << "[" << i-start << "]*pow(1-(x-" << xMin << ")/" << xMax << "," << i << ")";
            if (i != end) poly << " + ";
        }
        TString output = poly.str();
        return output;
    };
        
    TString ConvexPolyPairString(unsigned n, unsigned m, float xMin, float xMax) {
        std::ostringstream poly;
        poly << "[0]*pow(1-(x-" << xMin << ")/" << xMax << ", " << n << ")";
        poly << " + [1]*pow(1-(x-" << xMin << ")/" << xMax << ", " << m << ")";
        TString output = poly.str();
        return output;
    };

    TString ATLASString(unsigned k) {
        std::ostringstream atlas;
        atlas << "pow( (1-pow(x/13000.0,1.0/3.0)), [0])";
        if (k==0){
            atlas << "*pow(x/13000.0,[1])";
        }else{
            atlas << "*pow(x/13000.0,[1]";
            for (unsigned i(1);i<k+1;i++){
                atlas << " + [" << i+1 << "]*pow(TMath::Log(x/13000.0)," << i << ")";
            }
            atlas << ")";
        }
        TString output = atlas.str();
        return output;
    };

    TString ATLASPair(unsigned n, unsigned m) {
        std::ostringstream atlas;
        atlas << "pow((1-pow(x/13000.0,1.0/3.0)),[0])";
        atlas << "*pow(x/13000.0,";
        atlas << "[1]*pow(TMath::Log(x/13000.0)," << n << ")";
        atlas << "[2]*pow(TMath::Log(x/13000.0)," << m << ")";
        atlas << ")";
        TString output = atlas.str();
        return output;
    };

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
    };

    TString DijetPair(unsigned n, unsigned m) {
        std::ostringstream dijet;
        dijet << "TMath::Exp(";
        dijet << "[0]*pow(TMath::Log(x)," << n << ")";
        dijet << " + [1]*pow(TMath::Log(x)," << m << ")";
        dijet << ")";
        TString output = dijet.str();
        return output;
    };












    TString exponentialSumString(unsigned n) {
        std::ostringstream expo;
        expo << "exp( [0] + [1]*x )";
        if (n>0) {
            for (unsigned i(1);i<n+1;i++) {
                expo << " + exp( [" << i*2 << "] + [" << i*2+1 << "]*x )";
            }
        }
        TString output = expo.str();
        return output;
    }
        
    TString exponentialPoly(unsigned n) {
        std::ostringstream expo;
        expo << "exp( [0] ";
        if (n > 0) {
            for (unsigned i(1);i<n+1;i++) {
                expo << " + [" << i << "]*pow(x," << i << ")";
            }
        }
        expo << " )";
        TString output = expo.str();
        return output;
    }

    TString exponentialPowerLaw(unsigned n) {
        std::ostringstream expo;
        expo << "exp( [0]*pow(x,[1])";
        if (n >0) {
            for (unsigned i(1);i<n+1;i++){
                expo << " + [" << i*2 << "]*pow(x,[" << i*2 + 1 << "])";
            }
        }
        expo << " )";
        TString output = expo.str();
        return output;
    }

    TString powerLawSum(unsigned n) {
        std::ostringstream power;
        power << "[0]*pow(x,[1])";
        if (n>0) {
            for (unsigned i(0);i<n+1;i++) {
                power << " + [" << i*2 << "]*pow(x,[" << i*2+1 << "])";
            }
        }
        TString output = power.str();
        return output;
    }

    TString powerLawPoly(unsigned n) {
        std::ostringstream power;
        power << "pow( ";
        for (unsigned i(1);i<n;i++) {
            if (i-1 != 0) {power << " + ";}
            power << "[" << i-1 << "]*pow(x," << i << ")";
        }
        power << ", [" << n-1 << "])";
        TString output = power.str();
        return output;
    }
    
    TString polyPowerLaw(unsigned n) {
        std::ostringstream power;
        power << "[0]";
        if (n > 0) {
            power << "*pow(x,-1.0*([1]";
            for (unsigned i(2);i<n+1;i++) {
                power << " + [" << i << "]*pow(x," << i-1 << ")";
            }
            power << "))";
        }
        TString output = power.str();
        return output;
    }

    TString powerLawConvexBasis(unsigned n, unsigned m, float xMin, float xMax) {
        if (m > n) {std::swap(n,m);}
        std::ostringstream power;
        power << "[2]*pow( x, -1.0*( ";
        power << "[0]*pow(1-(x-" << xMin << ")/" << xMax << ", " << n << ")";
        power << " + [1]*pow(1-(x-" << xMin << ")/" << xMax << ", " << m << ")";
        power << " ) )";
        TString output = power.str();
        return output;
    }

    TString laurentString(unsigned n) {
        std::ostringstream laurent;
        laurent << "[0]";
        if (n > 0) {
            for (unsigned i(1);i<n+1;i++) {
                laurent << " + [" << i << "]*pow(x,-" << i << ")";
            }
        }
        TString output = laurent.str();
        return output;
    }

    TString laurentPairString(unsigned n, unsigned m) {
        std::ostringstream laurent;
        laurent << "[0]*pow(x,-" << n << ")";
        laurent << " + [1]*pow(x,-" << m << ")";
        TString output = laurent.str();
        return output;
    }














