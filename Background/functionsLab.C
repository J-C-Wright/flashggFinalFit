#include <vector>
#include<iostream>

    void functionsLab() {

        TFile *output = new TFile("FitPlots/EBEB_BackgroundFits.root","RECREATE");

        TFile *_file0 = TFile::Open("/afs/cern.ch/user/m/musella/public/workspace/exo/full_analysis_spring15_7415v2_sync_v5_data_ecorr_cic2_final_ws.root");
        t_EBEB=_file0->Get("tree_data_cic2_EBEB");
        gStyle->SetOptStat(11111);
        gStyle->SetOptFit(11111);

        t_EBEB->Draw("mgg>>h_EBEB(34,230,910)");
        TH1F *h_EBEB = (TH1F*)gPad->GetPrimitive("h_EBEB");
        h_EBEB->SetMarkerColor(kBlack);
        h_EBEB->Draw("EP");


    //-------------------
    //Functions and fits
    //-------------------
        TCanvas c1("c1");
        c1.SetLogy();
        //c1.SetLogx();

        h_EBEB->SetStats(kFALSE);

        std::vector<TF1*> functions;

    //EBEB Region
        std::ofstream EBEB_FitsLog;
        EBEB_FitsLog.open("FitPlots/EBEB_Fits.txt");
        EBEB_FitsLog << "\n----------------\n";
        EBEB_FitsLog << "Fits to EBEB set\n";
        EBEB_FitsLog << "----------------\n\n";
        EBEB_FitsLog << setw(12) << "Name" << setw(12) << "ChiSq" << setw(12) << "NDOF";
        EBEB_FitsLog << setw(12) << "ChiSq/NDOF"; 
        EBEB_FitsLog << setw(12) << "Num Param" << "    Formula String" << std::endl;

    //Polynomial fits
        //Power basis
        EBEB_FitsLog << "\nPolynomials - Power basis" << std::endl;
        std::vector<TF1*> powerBasis(2);
        for (unsigned i(0);i<3;i++) {
            powerBasis[i] = new TF1(Form("polyPB%d",i+1),Form("pol%d",i+1),230,910);

            h_EBEB->Fit(Form("polyPB%d",i+1),"QRM");

            EBEB_FitsLog << setw(12) << powerBasis[i]->GetName() << setw(12) << powerBasis[i]->GetChisquare() << setw(12) << powerBasis[i]->GetNDF();
            EBEB_FitsLog << setw(12) << powerBasis[i]->GetChisquare()/powerBasis[i]->GetNDF();
            EBEB_FitsLog << setw(12) << powerBasis[i]->GetNumberFreeParameters() << "    " << Form("pol%d",i+1) << std::endl;
        }

        h_EBEB->Draw("EP");
        for (unsigned i(0);i<3;i++) {
            powerBasis[i]->Draw("SAME");
            functions.push_back(powerBasis[i]);
        }
        c1.Print("FitPlots/Polys_PowerBasis.pdf");

        //Bernstein basis
        EBEB_FitsLog << "\nPolynomials - Bernstein sub-basis (i<n/2)" << std::endl;
        std::vector<TF1*> bernsteinBasis(4);
        h_EBEB->Draw("EP");
        for (unsigned i(0);i<5;i++) {

            TString functionString = SubbasisBernsteinString(i+1,230,910);
            bernsteinBasis[i] = new TF1(Form("Bern_D%d",i+1),functionString,230,910);

            h_EBEB->Fit(Form("Bern_D%d",i+1),"QRM");

            EBEB_FitsLog << setw(12) << bernsteinBasis[i]->GetName() << setw(12) << bernsteinBasis[i]->GetChisquare() << setw(12) << bernsteinBasis[i]->GetNDF();
            EBEB_FitsLog << setw(12) << bernsteinBasis[i]->GetChisquare()/bernsteinBasis[i]->GetNDF();
            EBEB_FitsLog << setw(12) << bernsteinBasis[i]->GetNumberFreeParameters() << "    " << functionString << std::endl;
        }

        h_EBEB->Draw("EP");
        for (unsigned i(0);i<5;i++) {
            bernsteinBasis[i]->Draw("SAME");
            functions.push_back(bernsteinBasis[i]);
        }
        c1.Print("FitPlots/Polys_BernsteinBasis.pdf");






    //Power Laws
        //Power Law Sum (includes single power law
        EBEB_FitsLog << "\nPower Law Sum" << std::endl;
        std::vector<TF1*> powerLawSums(4);
        for (unsigned i(0);i<powerLawSums.size();i++) {
            TString funcString =  powerLawSum(i);
            powerLawSums[i] = new TF1(Form("plSums%d",i),funcString,230,910);

            h_EBEB->Fit(Form("plSums%d",i),"QRM");

            EBEB_FitsLog << setw(12) << powerLawSums[i]->GetName() << setw(12) << powerLawSums[i]->GetChisquare() << setw(12) << powerLawSums[i]->GetNDF();
            EBEB_FitsLog << setw(12) << powerLawSums[i]->GetChisquare()/powerLawSums[i]->GetNDF();
            EBEB_FitsLog << setw(12) << powerLawSums[i]->GetNumberFreeParameters() << "    " << funcString << std::endl;
        }
        h_EBEB->Draw("EP");
        for (unsigned i(0);i<powerLawSums.size();i++) {
            powerLawSums[i]->SetLineColor(i+1);
            powerLawSums[i]->Draw("SAME");
            functions.push_back(powerLawSums[i]);
        }
        c1.Print("FitPlots/PowerLawSums.pdf");

        //Polynomial Power Law
        EBEB_FitsLog << "\nPolynomial power law (Poly in place of x)" << std::endl;
        std::vector<TF1*> polyPowerLaws(4);
        for (unsigned i(0);i<polyPowerLaws.size();i++) {
            TString funcString = powerLawPoly(i+2);
            polyPowerLaws[i] = new TF1(Form("polyPL%d",i),funcString,230,910);

            h_EBEB->Fit(Form("polyPL%d",i),"QRM");

            EBEB_FitsLog << setw(12) << polyPowerLaws[i]->GetName() << setw(12) << polyPowerLaws[i]->GetChisquare() << setw(12) << polyPowerLaws[i]->GetNDF();
            EBEB_FitsLog << setw(12) << polyPowerLaws[i]->GetChisquare()/polyPowerLaws[i]->GetNDF();
            EBEB_FitsLog << setw(12) << polyPowerLaws[i]->GetNumberFreeParameters() << "    " << funcString << std::endl;
        }
        h_EBEB->Draw("EP");
        for (unsigned i(0);i<polyPowerLaws.size();i++) {
            polyPowerLaws[i]->SetLineColor(i+1);
            polyPowerLaws[i]->Draw("SAME");
            functions.push_back(polyPowerLaws[i]);
        }
        c1.Print("FitPlots/PolynomialPowerLaw.pdf");

        //Poly-power Law
        EBEB_FitsLog << "\nPolynomial-power law (Poly in the power)" << std::endl;
        std::vector<TF1*> polyPowerLawsAlt(4);
        for (unsigned i(0);i<polyPowerLawsAlt.size();i++) {
            TString funcString = polyPowerLaw(i+1);
            polyPowerLawsAlt[i] = new TF1(Form("polyPLAlt%d",i+1),funcString,230,910);

            h_EBEB->Fit(Form("polyPLAlt%d",i+1),"QRM");

            EBEB_FitsLog << setw(12) << polyPowerLawsAlt[i]->GetName() << setw(12) << polyPowerLawsAlt[i]->GetChisquare() << setw(12) << polyPowerLawsAlt[i]->GetNDF();
            EBEB_FitsLog << setw(12) << polyPowerLawsAlt[i]->GetChisquare()/polyPowerLawsAlt[i]->GetNDF();
            EBEB_FitsLog << setw(12) << polyPowerLawsAlt[i]->GetNumberFreeParameters() << "    " << funcString << std::endl;
        }
        h_EBEB->Draw("EP");
        for (unsigned i(0);i<polyPowerLawsAlt.size();i++) {
            polyPowerLawsAlt[i]->SetLineColor(i+1);
            polyPowerLawsAlt[i]->Draw("SAME");
            functions.push_back(polyPowerLawsAlt[i]);
        }
        c1.Print("FitPlots/PolynomialPowerLawAlt.pdf");

    //Exponentials
        //Exponential Sum
        EBEB_FitsLog << "\nExponential Sums" << std::endl;
        std::vector<TF1*> expoSums(4);
        for (unsigned i(0);i<expoSums.size();i++) {
            TString funcString =  exponentialSumString(i);
            expoSums[i] = new TF1(Form("ExpoSum%d",i),funcString,230,910);

            h_EBEB->Fit(Form("ExpoSum%d",i),"QRM");

            EBEB_FitsLog << setw(12) << expoSums[i]->GetName() << setw(12) << expoSums[i]->GetChisquare() << setw(12) << expoSums[i]->GetNDF();
            EBEB_FitsLog << setw(12) << expoSums[i]->GetChisquare()/expoSums[i]->GetNDF();
            EBEB_FitsLog << setw(12) << expoSums[i]->GetNumberFreeParameters() << "    " << funcString << std::endl;
        }
        h_EBEB->Draw("EP");
        for (unsigned i(0);i<expoSums.size();i++) {
            expoSums[i]->SetLineColor(i+1);
            expoSums[i]->Draw("SAME");
            functions.push_back(expoSums[i]); 
        }
        c1.Print("FitPlots/ExpoSums.pdf");

        //Exponential Polynomial
        EBEB_FitsLog << "\nExponentiated Polynomial" << std::endl;
        std::vector<TF1*> expoPolys(4);
        for (unsigned i(0);i<expoPolys.size();i++) {
            TString funcString = exponentialPoly(i+1);
            expoPolys[i] = new TF1(Form("ExpoPoly%d",i+1),funcString,230,910);

            h_EBEB->Fit(Form("ExpoPoly%d",i+1),"QRM");

            EBEB_FitsLog << setw(12) << expoPolys[i]->GetName() << setw(12) << expoPolys[i]->GetChisquare() << setw(12) << expoPolys[i]->GetNDF();
            EBEB_FitsLog << setw(12) << expoPolys[i]->GetChisquare()/expoPolys[i]->GetNDF();
            EBEB_FitsLog << setw(12) << expoPolys[i]->GetNumberFreeParameters() << "    " << funcString << std::endl;
        }
        h_EBEB->Draw("EP");
        for (unsigned i(0);i<expoPolys.size();i++) {
            expoPolys[i]->SetLineColor(i+1);
            expoPolys[i]->Draw("SAME");
            functions.push_back(expoPolys[i]);
        }
        c1.Print("FitPlots/ExpoPolys.pdf");

        //Exponential Power Law Sum
        EBEB_FitsLog << "\nExponentiated Power Law" << std::endl;
        std::vector<TF1*> expoPowers(4);
        for (unsigned i(0);i<expoPowers.size();i++) {
            TString funcString = exponentialPowerLaw(i);
            expoPowers[i] = new TF1(Form("ExpoPow%d",i),funcString,230,910);

            h_EBEB->Fit(Form("ExpoPow%d",i),"QRM");

            EBEB_FitsLog << setw(12) << expoPowers[i]->GetName() << setw(12) << expoPowers[i]->GetChisquare() << setw(12) << expoPowers[i]->GetNDF();
            EBEB_FitsLog << setw(12) << expoPowers[i]->GetChisquare()/expoPowers[i]->GetNDF();
            EBEB_FitsLog << setw(12) << expoPowers[i]->GetNumberFreeParameters() << "    " << funcString << std::endl;
        }
        h_EBEB->Draw("EP");
        for (unsigned i(0);i<expoPowers.size();i++) {
            expoPowers[i]->SetLineColor(i+1);
            expoPowers[i]->Draw("SAME");
            functions.push_back(expoPowers[i]); 
        }
        c1.Print("FitPlots/ExpoPowers.pdf");

        //Exponential times power law
        TString expTimesPowString = "exp([0]+[1])*pow(x,[2])";
        TF1 *expTimesPow = new TF1("ExpTimesPow",expTimesPowString,230,910);

        h_EBEB->Fit("ExpTimesPow","QRM");

        h_EBEB->Draw("EP");
        expTimesPow->Draw("SAME");
        c1.Print("FitPlots/ExpTimesPowerLaw.pdf");

        EBEB_FitsLog << "\nExponential times power law" << std::endl;
        EBEB_FitsLog << setw(12) << expTimesPow->GetName() << setw(12) << expTimesPow->GetChisquare();
        EBEB_FitsLog << setw(12) << expTimesPow->GetNDF();
        EBEB_FitsLog << setw(12) << expTimesPow->GetChisquare()/expTimesPow->GetNDF();
        EBEB_FitsLog << setw(12) << expTimesPow->GetNumberFreeParameters();
        EBEB_FitsLog << "    " << expTimesPowString << std::endl;
        functions.push_back(expTimesPow);

    //Laurent
        //Laurent Series Function
        EBEB_FitsLog << "\nLaurent series function" << std::endl;
        std::vector<TF1*> laurents(6);
        for (unsigned i(0);i<laurents.size();i++) {
            TString funcString = laurentString(i+1);
            laurents[i] = new TF1(Form("Laurent%d",i+1),funcString,230,910);
            h_EBEB->Fit(Form("Laurent%d",i+1),"QRM");
            EBEB_FitsLog << setw(12) << laurents[i]->GetName() << setw(12) << laurents[i]->GetChisquare() << setw(12) << laurents[i]->GetNDF();
            EBEB_FitsLog << setw(12) << laurents[i]->GetChisquare()/laurents[i]->GetNDF();
            EBEB_FitsLog << setw(12) << laurents[i]->GetNumberFreeParameters() << "    " << funcString << std::endl;
        }
        h_EBEB->Draw("EP");
        for (unsigned i(0);i<laurents.size();i++) {
            laurents[i]->SetLineColor(i+1);
            laurents[i]->Draw("SAME");
            functions.push_back(laurents[i]); 
        }
        c1.Print("FitPlots/Laurents.pdf");

    //ATLAS Function
        //N parameters
        EBEB_FitsLog << "\nATLAS functions" << std::endl;
        std::vector<TF1*> atlasFunctions(3);
        for (unsigned k(0);k<atlasFunctions.size();k++){
            atlasFunctions[k] = new TF1(Form("Atlas%d",k),ATLASString(k),230,910);

            h_EBEB->Fit(Form("Atlas%d",k),"QRM");

            EBEB_FitsLog << setw(12) << atlasFunctions[k]->GetName() << setw(12) << atlasFunctions[k]->GetChisquare() << setw(12) << atlasFunctions[k]->GetNDF();
            EBEB_FitsLog << setw(12) << atlasFunctions[k]->GetChisquare()/atlasFunctions[k]->GetNDF();
            EBEB_FitsLog << setw(12) << atlasFunctions[k]->GetNumberFreeParameters() << "    " << ATLASString(k) << std::endl;
        }
        h_EBEB->Draw("EP");
        for (unsigned k(0);k<atlasFunctions.size();k++) {
            atlasFunctions[k]->SetLineColor(k+1);
            atlasFunctions[k]->Draw("SAME");
            functions.push_back(atlasFunctions[k]); 
        }
        c1.Print("FitPlots/ATLASFunctions.pdf");

    //Dijet Function
        //N parameters
        EBEB_FitsLog << "\nDijet functions" << std::endl;
        std::vector<TF1*> dijetFunctions(4);
        for (unsigned k(0);k<4;k++) {
            dijetFunctions[k] = new TF1(Form("Dijet%d",k+1),DijetString(k+1),230,910);
            h_EBEB->Fit(Form("Dijet%d",k+1),"QRM");
            EBEB_FitsLog << setw(12) << dijetFunctions[k]->GetName() << setw(12) << dijetFunctions[k]->GetChisquare() << setw(12) << dijetFunctions[k]->GetNDF();
            EBEB_FitsLog << setw(12) << dijetFunctions[k]->GetChisquare()/dijetFunctions[k]->GetNDF();
            EBEB_FitsLog << setw(12) << dijetFunctions[k]->GetNumberFreeParameters() << "    " << DijetString(k+1) << std::endl;
        }
        h_EBEB->Draw("EP");
        for (unsigned k(0);k<dijetFunctions.size();k++) {
            dijetFunctions[k]->SetLineColor(k+1);
            dijetFunctions[k]->Draw("SAME");
            functions.push_back(dijetFunctions[k]); 
        }
        c1.Print("FitPlots/DijetFunctions.pdf");

    //From analysis note
        //Benchmark 
        TString benchmarkString = "pow( x, [0] + [1]*log(x) )";
        TF1 *benchmark = new TF1("Bench",benchmarkString,230,910);
        h_EBEB->Fit("Bench","QRM");
        h_EBEB->Draw("EP");
        benchmark->Draw("same");
        c1.Print("FitPlots/Benchmark.pdf");
        EBEB_FitsLog << "\nBenchmark function from paper" << std::endl;
        EBEB_FitsLog << setw(12) << benchmark->GetName() << setw(12) << benchmark->GetChisquare();
        EBEB_FitsLog << setw(12) << benchmark->GetNDF();
        EBEB_FitsLog << setw(12) << benchmark->GetChisquare()/benchmark->GetNDF() << setw(12) << benchmark->GetNumberFreeParameters();
        EBEB_FitsLog << "    " << benchmarkString << std::endl;
        functions.push_back(benchmark);

    //From Monte Carlo studies
        //Laurent
        TString lPairString = laurentPairString(5,6);
        TF1* laurentMCPair = new TF1("LPairMC",lPairString,230,910);
        h_EBEB->Fit("LPairMC","RM");

        EBEB_FitsLog << "\nLaurent Pair chosen from MC" << std::endl;
        EBEB_FitsLog << setw(12) << laurentMCPair->GetName() << setw(12) << laurentMCPair->GetChisquare();
        EBEB_FitsLog << setw(12) << laurentMCPair->GetNDF();
        EBEB_FitsLog << setw(12) << laurentMCPair->GetChisquare()/laurentMCPair->GetNDF() << setw(12) << laurentMCPair->GetNumberFreeParameters();
        EBEB_FitsLog << "    " << laurentPairString(5,6) << std::endl;
        functions.push_back(laurentMCPair);

        //ATLAS
        TString atlasPairString = ATLASPair(6,7);
        TF1* ATLASMCPair = new TF1("AtlasPairMC",atlasPairString,230,910);
        h_EBEB->Fit("AtlasPairMC","RM");

        EBEB_FitsLog << "\nLaurent Pair chosen from MC" << std::endl;
        EBEB_FitsLog << setw(12) << ATLASMCPair->GetName() << setw(12) << ATLASMCPair->GetChisquare();
        EBEB_FitsLog << setw(12) << ATLASMCPair->GetNDF();
        EBEB_FitsLog << setw(12) << ATLASMCPair->GetChisquare()/ATLASMCPair->GetNDF() << setw(12) << ATLASMCPair->GetNumberFreeParameters();
        EBEB_FitsLog << "    " << atlasPairString << std::endl;
        functions.push_back(ATLASMCPair);
        
        //Dijet
        TString dijetPairString = DijetPair(2,1);
        TF1* DijetMCPair = new TF1("DijetPairMC",dijetPairString,230,910);
        h_EBEB->Fit("DijetPairMC","RM");

        EBEB_FitsLog << "\nLaurent Pair chosen from MC" << std::endl;
        EBEB_FitsLog << setw(12) << DijetMCPair->GetName() << setw(12) << DijetMCPair->GetChisquare();
        EBEB_FitsLog << setw(12) << DijetMCPair->GetNDF();
        EBEB_FitsLog << setw(12) << DijetMCPair->GetChisquare()/DijetMCPair->GetNDF() << setw(12) << DijetMCPair->GetNumberFreeParameters();
        EBEB_FitsLog << "    " << dijetPairString << std::endl;
        functions.push_back(DijetMCPair);
 
        

        //Expo poly










        
//Ranking and analysis
    //Ranking by ChiSq/NDOF and separating by number of parameters
        std::vector<TF1*> goodFits;
        for (unsigned i(0);i<functions.size();i++) {
            if (functions[i]->GetChisquare()/functions[i]->GetNDF() >= 0.8) continue;
            float funcChiSqPerDof = functions[i]->GetChisquare()/functions[i]->GetNDF();
            unsigned insertionIndex(0);
            for (unsigned j(0);j<goodFits.size();j++) {
                float goodChiSqPerDof = goodFits[j]->GetChisquare()/goodFits[j]->GetNDF();
                if (funcChiSqPerDof >= goodChiSqPerDof) insertionIndex = j+1;
            }
            goodFits.insert( goodFits.begin() + insertionIndex, functions[i] );
        }
        output->cd();
        h_EBEB->Write();
        for (unsigned i(0);i<goodFits.size();i++) {
            goodFits[i]->Write();
        }
        output->Close();

        EBEB_FitsLog << "\n\nSummary of 2/3 parameter fits that perform better than threshold" << std::endl;
        for (unsigned j(2);j<5;j++) {
            EBEB_FitsLog << setw(12) << j << "-parameter functions:" << std::endl;
            EBEB_FitsLog << setw(18) << "Name" << setw(12) << "ChiSq/DOF";
            EBEB_FitsLog << setw(12) << "Num Param" << "  Formula" << std::endl;
            count = 0;
            for (unsigned i(0);i<goodFits.size();i++) {
                if (goodFits[i]->GetNumberFreeParameters() != j) continue;
                count++;
                EBEB_FitsLog << setw(18) << goodFits[i]->GetName();
                EBEB_FitsLog << setw(12) << goodFits[i]->GetChisquare()/goodFits[i]->GetNDF();
                EBEB_FitsLog << setw(12) << goodFits[i]->GetNumberFreeParameters();
                EBEB_FitsLog << "  " << goodFits[i]->GetExpFormula("CLING") << std::endl;
            }
            EBEB_FitsLog << "There are " << count << " " << j << "-parameter functions with ChiSq/DOF < 0.8" << std::endl;
        }


        //Comparisons of benchmark to set of 2param and 3param with better ChiSq/DOF
        h_EBEB->Draw("EP");
        std::vector<TF1*> bestTwoParams;
        for (unsigned i(0);i<goodFits.size();i++) {
            if (goodFits[i]->GetNumberFreeParameters() != 2) continue;
            goodFits[i]->SetLineColor(kBlue);
            goodFits[i]->Draw("SAME");
            bestTwoParams.push_back(goodFits[i]);
        }
        benchmark->SetLineColor(kRed);
        benchmark->Draw("SAME");
        c1.Print("FitPlots/TwoParameterFits.pdf");
        
 
        //Sim fits of best 3 parameter functions
        h_EBEB->Draw("EP");
        std::vector<TF1*> bestThreeParams;
        for (unsigned i(0);i<goodFits.size();i++) {
            if (goodFits[i]->GetNumberFreeParameters() != 3) continue;
            goodFits[i]->SetLineColor(kBlue);
            goodFits[i]->Draw("SAME");
            bestThreeParams.push_back(goodFits[i]);
        }
        benchmark->SetLineColor(kRed);
        benchmark->Draw("SAME");
        c1.Print("FitPlots/ThreeParameterFits.pdf");      

    //Sim fit
        //Sim fits of best 2 parameter functions
        std::vector<TF1*> bestTwoParamSimFits;
        EBEB_FitsLog << "\nSimultaneous fits of the best 2-parameter functions" << std::endl;
        EBEB_FitsLog << setw(18) << "Name" << setw(12) << "ChiSq/DOF";
        EBEB_FitsLog << setw(12) << "Num Param" << "    Formula" << std::endl;
        for (unsigned i(0);i<bestTwoParams.size();i++) {

            //Convex ones actually suck
            if (bestTwoParams[i]->GetName() == TString("Convex_9_1")) continue;
            if (bestTwoParams[i]->GetName() == TString("Convex_8_1")) continue;

            float sqrtS(0.0);
            if (bestTwoParams[i]->GetName() == TString("Atlas0")) sqrtS = 13000.0;

            TF1* simFit = simultaneousFit(bestTwoParams[i],0,0,sqrtS);
            bestTwoParamSimFits.push_back(simFit);
            h_EBEB->Fit(bestTwoParams[i]->GetName() + TString("_SF"),"IRM");
            EBEB_FitsLog << setw(18) << simFit->GetName(); 
            EBEB_FitsLog << setw(12) << simFit->GetChisquare()/simFit->GetNDF();
            EBEB_FitsLog << setw(12) << simFit->GetNumberFreeParameters() << "    " << simFit->GetExpFormula("CLING") << std::endl;
        }

        h_EBEB->Draw("EP");
        leg = new TLegend(0.6,0.6,0.9,0.9);
        for (unsigned i(0);i<bestTwoParamSimFits.size();i++) {
            bestTwoParamSimFits[i]->Draw("SAME");
            bestTwoParamSimFits[i]->SetLineColor(i+1);
            leg->AddEntry(bestTwoParamSimFits[i]->GetName(),bestTwoParamSimFits[i]->GetName(),"l");
        }
        leg->Draw();
        c1.Print("FitPlots/SimFitBestTwoParam.pdf");

 
        EBEB_FitsLog << "\nSimultaneous fits of the best 3-parameter functions" << std::endl;
        EBEB_FitsLog << setw(18) << "Name" << setw(12) << "ChiSq/DOF";
        EBEB_FitsLog << setw(12) << "Num Param" << "    Formula" << std::endl;
        std::vector<TF1*> bestThreeParamSimFits;
        for (unsigned i(0);i<bestThreeParams.size();i++) {
            TString name =  bestThreeParams[i]->GetName();
            std::cout << name << std::endl;
            float sqrtS(0.0);
            if (name.Contains("Atlas")) {
                sqrtS = 13000.0;
            }
            TF1* simFit = simultaneousFit(bestThreeParams[i],0,0,sqrtS);

            bestThreeParamSimFits.push_back(simFit);
            h_EBEB->Fit(name + TString("_SF"),"IRM");
            EBEB_FitsLog << setw(18) << simFit->GetName(); 
            EBEB_FitsLog << setw(12) << simFit->GetChisquare()/simFit->GetNDF();
            EBEB_FitsLog << setw(12) << simFit->GetNumberFreeParameters() << "    " << simFit->GetExpFormula("CLING") << std::endl;
        }

        h_EBEB->Draw("EP");
        leg2 = new TLegend(0.6,0.6,0.9,0.9);
        for (unsigned i(0);i<bestThreeParamSimFits.size();i++) {
            bestThreeParamSimFits[i]->Draw("SAME");
            bestThreeParamSimFits[i]->SetLineColor(i+1);
            leg2->AddEntry(bestThreeParamSimFits[i]->GetName(),bestThreeParamSimFits[i]->GetName(),"l");
        }
        leg2->Draw();
        c1.Print("FitPlots/SimFitBestThreeParam.pdf");






    }//End of main








//---------
//Functions
//---------

    TString BernsteinString(unsigned n, float xMin, float xMax) {
        if (n == 0) {return TString("[0]");}
        std::ostringstream poly;
        for (unsigned i(0);i<n+1;i++) {
            unsigned factorialPart = (unsigned)TMath::Factorial(n)/( TMath::Factorial(i)*TMath::Factorial(n-i) );
            poly << "[" << i << "]*";
            if (factorialPart > 1) {poly << factorialPart << ".0*";}
            if (i > 1) {poly << "pow((x - " << xMin << ")/" << xMax-xMin << "," << i << ")";} else if (i == 1) {poly << "((x-" << xMin << ")/" << xMax-xMin << ")";}
            if (n != i && i != 0) {poly << "*";}
            if (n-i > 1) {poly << "pow(1-(x-" << xMin << ")/" << xMax-xMin << "," << n-i << ")";} else if (n-i == 1) {poly << "(1-(x-" << xMin << ")/" << xMax-xMin  << ")";}
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
            if (i > 1) {poly << "pow((x - " << xMin << ")/" << xMax-xMin << "," << i << ")";} else if (i == 1) {poly << "((x-" << xMin << ")/" << xMax-xMin << ")";}
            if (n != i && i != 0) {poly << "*";}
            if (n-i > 1) {poly << "pow(1-(x-" << xMin << ")/" << xMax-xMin << "," << n-i << ")";} else if (n-i == 1) {poly << "(1-(x-" << xMin << ")/" << xMax-xMin << ")";}
            if (i < iLimit) {poly << " + ";}
        }
        TString output = poly.str();
        return output;
    };

    TString ConvexPolyString(unsigned start, unsigned end, float xMin, float xMax) {
        std::ostringstream poly;
        for (unsigned i(start);i<end+1;i++) {
            poly << "fabs([" << i-start << "])*pow(1-(x-" << xMin << ")/" << xMax-xMin << "," << i << ")";
            if (i != end) poly << " + ";
        }
        TString output = poly.str();
        return output;
    };
        
    TString ConvexPolyPairString(unsigned n, unsigned m, float xMin, float xMax) {
        std::ostringstream poly;
        poly << "fabs([0])*pow(1-(x-" << xMin << ")/" << xMax-xMin << ", " << n << ")";
        poly << " + fabs([1])*pow(1-(x-" << xMin << ")/" << xMax-xMin << ", " << m << ")";
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
        atlas << "+[2]*pow(TMath::Log(x/13000.0)," << m << ")";
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

    TString exponentialPolyPair(unsigned n, unsigned m){
        std::ostringstream expo;
        expo << "exp( [0] ";
        expo << " + [1]*pow(x," << n << ")";
        expo << " + [2]*pow(x," << m << ")";
        expo << ")";
        TString output = expo.str();
        return output;
    }

    TString exponentialPowerLaw(unsigned n) {
        std::ostringstream expo;
        expo << "exp( [0] + [1]*pow(x,[2])";
        if (n >0) {
            for (unsigned i(1);i<n+1;i++){
                expo << " + [" << i*2 + 1 << "]*pow(x,[" << i*2 + 2 << "])";
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
        power << "[0]*pow( ";
        for (unsigned i(1);i<n;i++) {
            if (i-1 != 0) {power << " + [" << i-1 << "]*";}
            power << "pow(x," << i << ")";
        }
        power << ", [" << n-1 << "])";
        TString output = power.str();
        return output;
    }

    TString powerLawPolyPair(unsigned n, unsigned m) {
        std::ostringstream power;
        power << "[0]*pow( ";
        power << "[1]*pow(x," << n << ")";
        power << " + [2]*pow(x," << m << ")";
        power << ",[3])";
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

    TString polyPowerLawPair(unsigned n, unsigned m) {
        std::ostringstream power;
        power << "[0]*pow(x,-1.0*(";
        power << "[1]*pow(x," << n << ")";
        power << " + [2]*pow(x," << m << ")";
        power << "))";
        TString output = power.str();
        return output;
    }

    TString powerLawConvexBasis(unsigned n, unsigned m, float xMin, float xMax) {
        if (m > n) {std::swap(n,m);}
        std::ostringstream power;
        power << "[2]*pow( x, -1.0*( ";
        power << "fabs([0])*pow(1-(x-" << xMin << ")/" << xMax-xMax << ", " << n << ")";
        power << " + fabs([1])*pow(1-(x-" << xMin << ")/" << xMax-xMax << ", " << m << ")";
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

    TF1* simultaneousFit(TF1* input, float xMin, float xMax, float sqrtS) {

        unsigned numParams = input->GetNumberFreeParameters();

        float width = 20.5/2.35;

        TString formulaString = input->GetExpFormula("CLING");
        std::ostringstream simFit;
        simFit << formulaString;
        if (xMin == 0.0 && xMax == 0.0 && sqrtS == 0.0) {
            //This is the case where the variable is not scaled
            simFit << " + [" << numParams << "]*exp( -1.0*pow( x - 750, 2)/(2*pow(" << width << ",2)))";
        }else if (xMin != 0.0 && xMax != 0.0 && sqrtS == 0.0) {
            //This is the case where the m range is scaled to (0,1)
            simFit << " + [" << numParams << "]*exp( -1.0*pow( (x-" << xMin << ")/" << xMax-xMin << " - (750-" << xMin << ")/" << xMax-xMin << ", 2 )";
            simFit << "/( 2*pow(" << width/(xMax-xMin) << ",2)))";
        }else if (xMin == 0.0 && xMax == 0.0 && sqrtS != 0.0) {
            //This is the case when m is scaled by 1/sqrt(s)
            simFit << " + [" << numParams << "]*exp( -1.0*pow( x/" << sqrtS << " - 750/" << sqrtS << ",2)";
            simFit << "/( 2*pow(" << width/sqrtS << ",2)))";
        }
        TString simFitString = simFit.str(); 

        //Make new TF1
        TF1* simFitFunction = new TF1(input->GetName() + TString("_SF"), simFitString, 230, 910);
        //Set parameters
        for (unsigned parameter(0);parameter<numParams;parameter++) {
            simFitFunction->SetParameter(parameter, input->GetParameter(parameter));
        }
        //Gaussian parameters
        if (xMin == 0.0 && xMax == 0.0 && sqrtS == 0.0) {
            simFitFunction->SetParameter( numParams, input(750)*3.0 );
        }else if (xMin != 0.0 && xMax != 0.0 && sqrtS == 0.0) {
            simFitFunction->SetParameter( numParams, input((750-xMin)/(xMax-xMin))*3.0 );
        }else if (xMin == 0.0 && xMax == 0.0 && sqrtS != 0.0) {
            simFitFunction->SetParameter( numParams, input(750/sqrtS)*3.0 );
        }

        return simFitFunction;
    }











