#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TString.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TMath.h>
#include <TLatex.h>
#include <TMatrixD.h>
#include <TVectorD.h>
#include <THStack.h>

#include <iostream>
#include <vector>
#include <sstream>
#include <string>

#include "Util.h"
#include "helper_functions.cxx"
#include "constants.h"

using namespace std;
using namespace constants;

//----------------------------------------//

void data_mc_overlay(bool PlotGENIE = true, bool PlotGen = false, bool plot_closure = false) {

	//----------------------------------------//

	TH1D::SetDefaultSumw2();
	gStyle->SetEndErrorSize(6);		

	TString PathToFiles = xsec_path;

	//----------------------------------------//

	TString Extra = "";
	if (PlotGENIE) { Extra = "Genie"; }
	if (PlotGen) { Extra = "Gene"; }
	if (plot_closure) { Extra = "Closure"; }

	//----------------------------------------//

	vector<TString> PlotNames;
	PlotNames.push_back("SingleBinPlot");
	//PlotNames.push_back("Pi0MomentumPlot");	

	const int N1DPlots = PlotNames.size();
	cout << "Number of 1D Plots = " << N1DPlots << endl;

	//----------------------------------------//

	int NRuns = (int)(xsec_Runs.size());
	cout << "Number of Runs = " << NRuns << endl;

	//----------------------------------------//

	for (int WhichRun = 0; WhichRun < NRuns; WhichRun++) {

		vector<vector<TH1D*> > PlotsStatReco; PlotsStatReco.clear();
		vector<vector<TH1D*> > PlotsNormOnly; PlotsNormOnly.clear();		
		vector<vector<TH1D*> > PlotsReco; PlotsReco.clear();
		vector<vector<TH1D*> > PlotsTrue; PlotsTrue.clear();
		vector<vector<TH1D*> > QEPlotsTrue; QEPlotsTrue.clear();
		vector<vector<TH1D*> > MECPlotsTrue; MECPlotsTrue.clear();
		vector<vector<TH1D*> > RESPlotsTrue; RESPlotsTrue.clear();
		vector<vector<TH1D*> > DISPlotsTrue; DISPlotsTrue.clear();	
		vector<vector<TH1D*> > COHPlotsTrue; COHPlotsTrue.clear();									

		gStyle->SetPalette(55); const Int_t NCont = 999; gStyle->SetNumberContours(NCont); gStyle->SetTitleSize(0.07,"t");

		vector<TString> NameOfSamples; NameOfSamples.clear();
		vector<int> Colors; Colors.clear();		
		vector<TString> Labels; Labels.clear();
		vector<int> LineStyle; LineStyle.clear();
		vector<TString> weighted; weighted.clear();

		// CV

		NameOfSamples.push_back("mcc9_10_Overlay9"); Colors.push_back(OverlayColor); Labels.push_back("G18T "); LineStyle.push_back(kSolid); weighted.push_back("");

		//----------------------------------------//	

		if (PlotGENIE) {

			NameOfSamples.push_back("GENIE_v3_6_0_AR23");  Colors.push_back(kOrange+7); Labels.push_back("AR23 "); weighted.push_back("");
			NameOfSamples.push_back("GENIE_v3_0_6"); Colors.push_back(kGreen+2); Labels.push_back("G18 "); LineStyle.push_back(kSolid); weighted.push_back("");

		}

		//----------------------------------------//		

		if (PlotGen) {

			NameOfSamples.push_back("GENIE_v3_6_0_RS");  Colors.push_back(kGreen+1); Labels.push_back("G18 RS "); weighted.push_back("");
			NameOfSamples.push_back("NuWro_25_03_1"); Colors.push_back(kOrange+7); Labels.push_back("NuWro "); LineStyle.push_back(kSolid); weighted.push_back("");
			NameOfSamples.push_back("GiBUU_2025"); Colors.push_back(kMagenta+1); Labels.push_back("GiBUU "); LineStyle.push_back(kSolid); weighted.push_back(""); 
			NameOfSamples.push_back("NEUT_5_6_0"); Colors.push_back(kYellow-6); Labels.push_back("NEUT "); LineStyle.push_back(kSolid); weighted.push_back("");
			//NameOfSamples.push_back("NuWro_25_03_1_delta"); Colors.push_back(kBlue); Labels.push_back("NuWro Delta"); LineStyle.push_back(kSolid); weighted.push_back("");			
			//NameOfSamples.push_back("NuWro_25_03_1_tweaked_coh"); Colors.push_back(kBlue); Labels.push_back("NuWro Tweaked"); LineStyle.push_back(kSolid); weighted.push_back("");
			NameOfSamples.push_back("achilles"); Colors.push_back(kBlue); Labels.push_back("ACHILLES"); LineStyle.push_back(kSolid); weighted.push_back("");

		}			

		//----------------------------------------//		

		if (plot_closure) {

			NameOfSamples.push_back("GENIE_v3_0_6"); Colors.push_back(kGreen+2); Labels.push_back("G18 "); LineStyle.push_back(kSolid); weighted.push_back("");

		}

		//----------------------------------------//

		const int NSamples = NameOfSamples.size();
		vector<TFile*> FileSample; FileSample.clear();

		//----------------------------------------//

		// Open the files and grap the relevant plots

		for (int WhichSample = 0; WhichSample < NSamples; WhichSample ++) {

			vector<TH1D*> CurrentPlotsStatReco; CurrentPlotsStatReco.clear();
			vector<TH1D*> CurrentPlotsNormOnly; CurrentPlotsNormOnly.clear();			
			vector<TH1D*> CurrentPlotsReco; CurrentPlotsReco.clear();
			vector<TH1D*> CurrentPlotsTrue; CurrentPlotsTrue.clear();
			vector<TH1D*> QECurrentPlotsTrue; QECurrentPlotsTrue.clear();
			vector<TH1D*> MECCurrentPlotsTrue; MECCurrentPlotsTrue.clear();	
			vector<TH1D*> RESCurrentPlotsTrue; RESCurrentPlotsTrue.clear();
			vector<TH1D*> DISCurrentPlotsTrue; DISCurrentPlotsTrue.clear();	
			vector<TH1D*> COHCurrentPlotsTrue; COHCurrentPlotsTrue.clear();													

			// CV With Statistical Uncertainties

			if (NameOfSamples[WhichSample] == "mcc9_10_Overlay9") { // CV with statistical uncertainties only for now

				TString FileSampleName = PathToFiles+"/WienerSVD_ExtractedXSec_"+NameOfSamples[WhichSample]+"_"+xsec_Runs[WhichRun]+".root";
				FileSample.push_back(TFile::Open(FileSampleName,"readonly")); 

				for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++) {

					TH1D* histTotalReco = (TH1D*)(FileSample[WhichSample]->Get("StatReco"+PlotNames[WhichPlot]));
					CurrentPlotsStatReco.push_back(histTotalReco);								

					TH1D* histNormOnly = (TH1D*)(FileSample[WhichSample]->Get("NormOnlyReco"+PlotNames[WhichPlot]));
					CurrentPlotsNormOnly.push_back(histNormOnly);					

					TH1D* histReco = (TH1D*)(FileSample[WhichSample]->Get("RecoFullUnc"+PlotNames[WhichPlot]));
					CurrentPlotsReco.push_back(histReco);

					TString TrueString = "True";

					TH1D* histTrue = (TH1D*)(FileSample[WhichSample]->Get(TrueString+PlotNames[WhichPlot]));
					CurrentPlotsTrue.push_back(histTrue);

					TH1D* QEhistTrue = (TH1D*)(FileSample[WhichSample]->Get("QE"+TrueString+PlotNames[WhichPlot]));
					QECurrentPlotsTrue.push_back(QEhistTrue);
					TH1D* MEChistTrue = (TH1D*)(FileSample[WhichSample]->Get("MEC"+TrueString+PlotNames[WhichPlot]));
					MECCurrentPlotsTrue.push_back(MEChistTrue);
					TH1D* REShistTrue = (TH1D*)(FileSample[WhichSample]->Get("RES"+TrueString+PlotNames[WhichPlot]));
					RESCurrentPlotsTrue.push_back(REShistTrue);
					TH1D* DIShistTrue = (TH1D*)(FileSample[WhichSample]->Get("DIS"+TrueString+PlotNames[WhichPlot]));
					DISCurrentPlotsTrue.push_back(DIShistTrue);
					TH1D* COHhistTrue = (TH1D*)(FileSample[WhichSample]->Get("COH"+TrueString+PlotNames[WhichPlot]));
					COHCurrentPlotsTrue.push_back(COHhistTrue);																								     
		
				}

			}

			else {
		
			  FileSample.push_back(TFile::Open("output_files/" + weighted[WhichSample] + "analyzer_ncpi0_"+NameOfSamples[WhichSample]+".root")); 

				for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++) {

					TH1D* histTotalReco = nullptr;
					CurrentPlotsStatReco.push_back(histTotalReco);									

					TH1D* histNormOnly = nullptr;
					CurrentPlotsNormOnly.push_back(histNormOnly);					

					TH1D* histReco = nullptr;
					CurrentPlotsReco.push_back(histReco);

					TString TrueString = "True";
					TH1D* histTrue = (TH1D*)(FileSample[WhichSample]->Get(TrueString+PlotNames[WhichPlot]));
					CurrentPlotsTrue.push_back(histTrue);

					TH1D* QEhistTrue = (TH1D*)(FileSample[WhichSample]->Get("QE"+TrueString+PlotNames[WhichPlot]));
					QECurrentPlotsTrue.push_back(QEhistTrue);
					TH1D* MEChistTrue = (TH1D*)(FileSample[WhichSample]->Get("MEC"+TrueString+PlotNames[WhichPlot]));
					MECCurrentPlotsTrue.push_back(MEChistTrue);
					TH1D* REShistTrue = (TH1D*)(FileSample[WhichSample]->Get("RES"+TrueString+PlotNames[WhichPlot]));
					RESCurrentPlotsTrue.push_back(REShistTrue);
					TH1D* DIShistTrue = (TH1D*)(FileSample[WhichSample]->Get("DIS"+TrueString+PlotNames[WhichPlot]));
					DISCurrentPlotsTrue.push_back(DIShistTrue);
					TH1D* COHhistTrue = (TH1D*)(FileSample[WhichSample]->Get("COH"+TrueString+PlotNames[WhichPlot]));
					COHCurrentPlotsTrue.push_back(COHhistTrue);					
		
				}

			}

			PlotsStatReco.push_back(CurrentPlotsStatReco);
			PlotsNormOnly.push_back(CurrentPlotsNormOnly);					
			PlotsReco.push_back(CurrentPlotsReco);
			PlotsTrue.push_back(CurrentPlotsTrue);

			QEPlotsTrue.push_back(QECurrentPlotsTrue);			
			MECPlotsTrue.push_back(MECCurrentPlotsTrue);
			RESPlotsTrue.push_back(RESCurrentPlotsTrue);
			DISPlotsTrue.push_back(DISCurrentPlotsTrue);
			COHPlotsTrue.push_back(COHCurrentPlotsTrue);		

		}

		//----------------------------------------//

		// Loop over the plots

		for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++) {	

			//----------------------------------------//

			TH2D* Ac = (TH2D*)FileSample[0]->Get("Ac"+PlotNames[WhichPlot]);

			TString CovString = "UnfCov"+PlotNames[WhichPlot];
			TH2D* Cov = (TH2D*)FileSample[0]->Get(CovString);	
			TH2D* NormCov = (TH2D*)FileSample[0]->Get("NormUnfCov"+PlotNames[WhichPlot]);	
			TH2D* ShapeCov = (TH2D*)FileSample[0]->Get("ShapeUnfCov"+PlotNames[WhichPlot]);	

			//----------------------------------------//

			// The covariance matrix needs to be scaled by the 2D bin width

			TH2D* CovClone = (TH2D*)(Cov->Clone());
			TH2D* NormCovClone = (TH2D*)(NormCov->Clone());	
			TH2D* ShapeCovClone = (TH2D*)(ShapeCov->Clone());					 

			int n = Cov->GetXaxis()->GetNbins();

			for (int ix = 1; ix <= n; ix++) {

				for (int iy = 1; iy <= n; iy++) {

					double WidthX = Cov->GetXaxis()->GetBinWidth(ix);
					double WidthY = Cov->GetYaxis()->GetBinWidth(iy);

					double TwoDWidth = WidthX * WidthY;

					double BinContent = Cov->GetBinContent(ix,iy);
					double NewBinContent = BinContent/TwoDWidth;

					double NormBinContent = NormCov->GetBinContent(ix,iy);
					double NormNewBinContent = NormBinContent/TwoDWidth;

					double ShapeBinContent = ShapeCov->GetBinContent(ix,iy);
					double ShapeNewBinContent = ShapeBinContent/TwoDWidth;													

					// Only for the diagonal elements
					// Add the unfolding uncertainty
					// On top of everything else
					// That is done both for the final xsec result and for the unfolded covariance
					if (ix == iy) { 

						// unfolded covariance matrix
//						double UnfUncBin = UncHist->GetBinContent(ix);
						double UnfUncBin = 0.;

						NewBinContent = NewBinContent + TMath::Power(UnfUncBin,2.) ;
						ShapeNewBinContent = ShapeNewBinContent + TMath::Power(UnfUncBin,2.) ;						 

						// xsec uncertainty
						double CurrentUnc = PlotsReco[0][WhichPlot]->GetBinError(ix);
						double NewError = TMath::Sqrt( TMath::Power(CurrentUnc,2.) + TMath::Power(UnfUncBin,2.) ) ;
						PlotsReco[0][WhichPlot]->SetBinError(ix,NewError);					
						
					}

					CovClone->SetBinContent(ix,iy,NewBinContent);
					ShapeCovClone->SetBinContent(ix,iy,ShapeNewBinContent);
					NormCovClone->SetBinContent(ix,iy,NormNewBinContent);										

				}					

			}	

			//----------------------------------------//

			TCanvas* PlotCanvas = new TCanvas(PlotNames[WhichPlot]+"_"+xsec_Runs[WhichRun],PlotNames[WhichPlot]+"_"+xsec_Runs[WhichRun],205,34,1024,768);
			PlotCanvas->cd();

			TPad *midPad = new TPad("midPad", "",0.005, 0., 0.995, 0.995);
			midPad->SetBottomMargin(0.16);
			midPad->SetTopMargin(0.12);
			midPad->SetLeftMargin(0.2);
			midPad->SetRightMargin(0.02);			
			midPad->Draw();

			TLegend* leg = new TLegend(0.26,0.65,0.59,0.82);
			TLegend* legMC = new TLegend(0.6,0.64,0.7,0.82);

			leg->SetBorderSize(0);
			leg->SetTextSize(0.05);
			leg->SetTextFont(FontStyle);
			leg->SetNColumns(1);
			leg->SetMargin(0.15);
			leg->SetFillStyle(0);

			legMC->SetBorderSize(0);
			legMC->SetTextSize(0.05);
			legMC->SetTextFont(FontStyle);
			legMC->SetNColumns(1);
			legMC->SetMargin(0.3);	
			legMC->SetFillStyle(0);					

			//----------------------------------------//

			// BeamOn Total Uncertainty

			double MaxValue = PlotsReco[0][WhichPlot]->GetMaximum();
			int MaxValueBin = locate_bin_with_value(PlotsReco[0][WhichPlot],MaxValue);
			double MaxValueError = PlotsReco[0][WhichPlot]->GetBinError(MaxValueBin);

			double MinValue = PlotsReco[0][WhichPlot]->GetMinimum();
													
			PlotsReco[0][WhichPlot]->GetYaxis()->SetRangeUser(XSecRange[PlotNames[WhichPlot]].first,XSecRange[PlotNames[WhichPlot]].second);

			PlotsReco[0][WhichPlot]->SetLineColor(BeamOnColor);
			PlotsReco[0][WhichPlot]->SetMarkerColor(BeamOnColor);
			PlotsReco[0][WhichPlot]->SetMarkerSize(1.);
			PlotsReco[0][WhichPlot]->SetMarkerStyle(20);
			PlotsReco[0][WhichPlot]->SetLineWidth(1);	
			PlotsReco[0][WhichPlot]->GetYaxis()->SetTitle(VarLabel[PlotNames[WhichPlot]]);					
		
			midPad->cd();	

			//------------------------------//

			PlotsReco[0][WhichPlot]->GetXaxis()->SetTitleFont(FontStyle);
			PlotsReco[0][WhichPlot]->GetXaxis()->SetTitleSize(TextSize);
			PlotsReco[0][WhichPlot]->GetXaxis()->SetLabelFont(FontStyle);
			PlotsReco[0][WhichPlot]->GetXaxis()->SetLabelSize(TextSize);						
			PlotsReco[0][WhichPlot]->GetXaxis()->CenterTitle();	
			PlotsReco[0][WhichPlot]->GetXaxis()->SetNdivisions(8);					

			PlotsReco[0][WhichPlot]->GetYaxis()->SetTitleFont(FontStyle);	
			PlotsReco[0][WhichPlot]->GetYaxis()->SetTitleSize(TextSize);	
			PlotsReco[0][WhichPlot]->GetYaxis()->SetLabelFont(FontStyle);
			PlotsReco[0][WhichPlot]->GetYaxis()->SetLabelSize(TextSize);							
			PlotsReco[0][WhichPlot]->GetYaxis()->CenterTitle();
			PlotsReco[0][WhichPlot]->GetYaxis()->SetNdivisions(8);			

			PlotsReco[0][WhichPlot]->Draw("e1x0 same"); // Total Unc

			PlotsStatReco[0][WhichPlot]->SetLineColor(BeamOnColor);
			PlotsStatReco[0][WhichPlot]->SetMarkerColor(BeamOnColor);
			PlotsStatReco[0][WhichPlot]->SetLineWidth(1);			
			PlotsStatReco[0][WhichPlot]->Draw("e1x0 same"); // Stat Only	
					
			// PlotsNormOnly[0][WhichPlot]->SetFillColorAlpha(kGray+1, 0.75);	
			// PlotsNormOnly[0][WhichPlot]->SetLineColor(kGray+1);
			// PlotsNormOnly[0][WhichPlot]->SetMarkerColor(kGray+1);
			// if (PlotNames[WhichPlot] != "SingleBinPlot") { PlotsNormOnly[0][WhichPlot]->Draw("e2 hist same"); } // Norm unc Only					

			//----------------------------------------//

			// Overlay GENIE v3 + Tune

			PlotsTrue[0][WhichPlot]->SetLineColor(Colors[0]);
			PlotsTrue[0][WhichPlot]->SetMarkerColor(Colors[0]);
			PlotsTrue[0][WhichPlot]->SetLineStyle(kSolid);

			//----------------------------------------//

			// arrays for NSamples

			double Chi2[NSamples];			
			int Ndof[NSamples];
			double pval[NSamples];
			double sigma[NSamples];

			// Clones for the NSamples-1 model predictions
			// index 0 corresponds to nominal overlay / CV

			TH1D* Clone[NSamples-1];			

			for (int WhichSample = 1; WhichSample < NSamples; WhichSample++) {

				// Apply the additional smearing matrix Ac
				Clone[WhichSample-1] = Multiply(PlotsTrue[WhichSample][WhichPlot],Ac);
				QEPlotsTrue[WhichSample][WhichPlot] = Multiply(QEPlotsTrue[WhichSample][WhichPlot],Ac);		
				MECPlotsTrue[WhichSample][WhichPlot] = Multiply(MECPlotsTrue[WhichSample][WhichPlot],Ac);	
				RESPlotsTrue[WhichSample][WhichPlot] = Multiply(RESPlotsTrue[WhichSample][WhichPlot],Ac);	
				DISPlotsTrue[WhichSample][WhichPlot] = Multiply(DISPlotsTrue[WhichSample][WhichPlot],Ac);
				COHPlotsTrue[WhichSample][WhichPlot] = Multiply(COHPlotsTrue[WhichSample][WhichPlot],Ac);												
				// Divide by the bin width
				rm_bin_width(Clone[WhichSample-1],1.);
				rm_bin_width(QEPlotsTrue[WhichSample][WhichPlot],1.);
				rm_bin_width(MECPlotsTrue[WhichSample][WhichPlot],1.);
				rm_bin_width(RESPlotsTrue[WhichSample][WhichPlot],1.);
				rm_bin_width(DISPlotsTrue[WhichSample][WhichPlot],1.);
				rm_bin_width(COHPlotsTrue[WhichSample][WhichPlot],1.);							
			
				Clone[WhichSample-1]->SetLineColor(Colors[WhichSample]);
				Clone[WhichSample-1]->SetLineStyle(kSolid);
				Clone[WhichSample-1]->SetMarkerColor(Colors[WhichSample]);

				Clone[WhichSample-1]->SetLineWidth(3);		
				Clone[WhichSample-1]->Draw("hist same");		

				calc_chi2(Clone[WhichSample-1],PlotsReco[0][WhichPlot],CovClone,Chi2[WhichSample],Ndof[WhichSample],pval[WhichSample],sigma[WhichSample]);
				TString Chi2NdofAlt = "(" + to_string_with_precision(Chi2[WhichSample],1) + "/" + TString(std::to_string(Ndof[WhichSample])) +")";

				TLegendEntry* lGenie = legMC->AddEntry(Clone[WhichSample-1],Labels[WhichSample] + Chi2NdofAlt,"l");
				lGenie->SetTextColor(Colors[WhichSample]); 										

			}

			//----------------------------------------//

			// Legend & Run / POT

			double tor860_wcut = PeLEE_ReturnBeamOnRunPOT(xsec_Runs[WhichRun]);
			TString Label = ToString(tor860_wcut).ReplaceAll("e"," #times 10").ReplaceAll("+","^{")+"} POT";	
	
			//----------------------------------------//

			// chi2
			calc_chi2(PlotsTrue[0][WhichPlot],PlotsReco[0][WhichPlot],CovClone,Chi2[0],Ndof[0],pval[0],sigma[0]);
			TString Chi2NdofNom = "(" + to_string_with_precision(Chi2[0],1) + "/" + TString(std::to_string(Ndof[0])) +")";

			TLegendEntry* lGenie_GenieOverlay = legMC->AddEntry(PlotsTrue[0][WhichPlot],Labels[0]+Chi2NdofNom,"l");
			PlotsTrue[0][WhichPlot]->SetLineWidth(3); 
			PlotsTrue[0][WhichPlot]->Draw("hist same"); 
			lGenie_GenieOverlay->SetTextColor(Colors[0]); 

			//----------------------------------------//
			//----------------------------------------//

			PlotsStatReco[0][WhichPlot]->Draw("e1x0 same"); // Stat Only
			PlotsReco[0][WhichPlot]->Draw("e1x0 same"); // BeamOn Stat Total

			leg->AddEntry(PlotsReco[0][WhichPlot],"MicroBooNE Data","");
			leg->AddEntry(PlotsReco[0][WhichPlot],Label,"");
			leg->AddEntry(PlotsReco[0][WhichPlot],"Total unc","ep");
			//leg->AddEntry(PlotsReco[0][WhichPlot],"Stat #oplus Shape","ep");			
			//leg->AddEntry(PlotsNormOnly[0][WhichPlot],"Norm (#chi^{2}/ndf)","f"); 
			leg->Draw();			

			legMC->Draw();			

			TLatex *textSlice = new TLatex();
			textSlice->SetTextFont(FontStyle);
			textSlice->SetTextSize(0.06);
			TString PlotNameDuplicate = PlotNames[WhichPlot];
			TString ReducedPlotName = PlotNameDuplicate.ReplaceAll("Reco","") ;
			TString ReducedLatexLabel = LatexLabel[ ReducedPlotName ];
			textSlice->DrawLatexNDC(0.2, 0.92, ReducedLatexLabel.ReplaceAll("All events","") );

			TLatex *textPanel = new TLatex();
			textPanel->SetTextFont(FontStyle);
			textPanel->SetTextSize(TextSize);

			//----------------------------------------//

			// Saving the canvas with the data (total uncertainties) vs overlay & generator predictions

			PlotCanvas->SaveAs("/exp/uboone/data/users/"+UserID+"/ncpi0/myPlots/flat_ttree/"+Extra+"XSections_"+PlotNames[WhichPlot]+"_"+xsec_Runs[WhichRun]+".pdf");
			delete PlotCanvas;

			//----------------------------------------//
			//----------------------------------------//

			// interaction breakdown canvas			
			// start of the loop over the mc samples

			for (int WhichSample = 0; WhichSample < NSamples; WhichSample++) {
			
				// Canvas with interaction breakdown for each generator
			
				TCanvas* inte_can = new TCanvas(NameOfSamples[WhichSample]+"_" + PlotNames[WhichPlot]+"_"+xsec_Runs[WhichRun],NameOfSamples[WhichSample]+"_" + PlotNames[WhichPlot]+"_"+xsec_Runs[WhichRun],205,34,1024,768);
				inte_can->cd();
				inte_can->SetBottomMargin(0.14);
				inte_can->SetLeftMargin(0.16);

				//----------------------------------------//

				TLegend* ilegmc = new TLegend(0.3,0.59,0.73,0.89);
				ilegmc->SetBorderSize(0);
				ilegmc->SetTextSize(0.05);
				ilegmc->SetTextFont(FontStyle);
				ilegmc->SetNColumns(2);
				ilegmc->SetMargin(0.18);
				ilegmc->SetFillStyle(0);

				//----------------------------------------//
				
				TString THStackName = NameOfSamples[WhichSample] + "THStack_" + PlotNames[WhichPlot]+"_"+xsec_Runs[WhichRun];
				THStack* thstack = new THStack(THStackName,"");	

				//----------------------------------------//

				// plot data for the first time
				
				PlotsReco[0][WhichPlot]->GetYaxis()->SetTitleOffset(1.1);
				PlotsReco[0][WhichPlot]->Draw("e1x0 same"); // BeamOn Stat Total

				//----------------------------------------//

				// QE

				QEPlotsTrue[WhichSample][WhichPlot]->SetLineColor(OverlayColor);
				QEPlotsTrue[WhichSample][WhichPlot]->SetFillColor(OverlayColor);
				thstack->Add(QEPlotsTrue[WhichSample][WhichPlot],"hist");
				thstack->Draw("same");

				//----------------------------------------//

				// MEC

				MECPlotsTrue[WhichSample][WhichPlot]->SetLineColor(kOrange-3);
				MECPlotsTrue[WhichSample][WhichPlot]->SetFillColor(kOrange-3);
				thstack->Add(MECPlotsTrue[WhichSample][WhichPlot],"hist");
				thstack->Draw("same");

				//----------------------------------------//

				// RES

				RESPlotsTrue[WhichSample][WhichPlot]->SetLineColor(kGreen+1);
				RESPlotsTrue[WhichSample][WhichPlot]->SetFillColor(kGreen+1);
				thstack->Add(RESPlotsTrue[WhichSample][WhichPlot],"hist");
				thstack->Draw("same");

				//----------------------------------------//

				// DIS

				DISPlotsTrue[WhichSample][WhichPlot]->SetLineColor(kRed+1);
				DISPlotsTrue[WhichSample][WhichPlot]->SetFillColor(kRed+1);
				thstack->Add(DISPlotsTrue[WhichSample][WhichPlot],"hist");
				thstack->Draw("same");

				//----------------------------------------//

				// COH

				COHPlotsTrue[WhichSample][WhichPlot]->SetLineColor(kMagenta);
				COHPlotsTrue[WhichSample][WhichPlot]->SetFillColor(kMagenta);
				thstack->Add(COHPlotsTrue[WhichSample][WhichPlot],"hist");
				thstack->Draw("same");				

				//----------------------------------------//
				
				// plot the data points again so that they can be on top

				PlotsReco[0][WhichPlot]->Draw("e1x0 same"); // BeamOn Stat Total
				PlotsStatReco[0][WhichPlot]->Draw("e1x0 same"); // Stat Only
				//PlotsNormOnly[0][WhichPlot]->Draw("e2 hist same"); // norm only	

				//----------------------------------------//
			
				TH1D* hstack = (TH1D*)(thstack->GetStack()->Last());
	
				ilegmc->AddEntry(PlotsReco[0][WhichPlot],"MicroBooNE Data","");	
				
				ilegmc->AddEntry(PlotsReco[0][WhichPlot],"","");	
				
				ilegmc->AddEntry(PlotsReco[0][WhichPlot],Label,"");
				
				ilegmc->AddEntry(PlotsReco[0][WhichPlot],"","");	

//				ilegmc->AddEntry(PlotsReco[0][WhichPlot],"Stat#oplusShape","ep");
				ilegmc->AddEntry(PlotsReco[0][WhichPlot],"Total unc","ep");			
				//ilegmc->AddEntry(PlotsNormOnly[0][WhichPlot],"Norm","f");
	
				ilegmc->AddEntry(PlotsReco[0][WhichPlot],"#chi^{2}/ndf = " + to_string_with_precision(Chi2[WhichSample],1.) + "/"+ to_string_with_precision(Ndof[WhichSample],0) ,"");
				
				ilegmc->AddEntry(PlotsReco[0][WhichPlot],"p = " + to_string_with_precision(pval[WhichSample],2.),"");
				
				double qe_frac = QEPlotsTrue[WhichSample][WhichPlot]->Integral() / hstack->Integral() * 100.;
				TLegendEntry* lqe = ilegmc->AddEntry(QEPlotsTrue[WhichSample][WhichPlot],"QE (" + to_string_with_precision(qe_frac,1.) + "%)","f");
				lqe->SetTextColor(OverlayColor);	

				double mec_frac = MECPlotsTrue[WhichSample][WhichPlot]->Integral() / hstack->Integral() * 100.;
				TLegendEntry* lmec = ilegmc->AddEntry(MECPlotsTrue[WhichSample][WhichPlot],"MEC (" + to_string_with_precision(mec_frac,1.) + "%)","f");
				lmec->SetTextColor(kOrange-3);	

				double res_frac = RESPlotsTrue[WhichSample][WhichPlot]->Integral() / hstack->Integral() * 100.;
				TLegendEntry* lres = ilegmc->AddEntry(RESPlotsTrue[WhichSample][WhichPlot],"RES (" + to_string_with_precision(res_frac,1.) + "%)","f");
				lres->SetTextColor(kGreen+1);	

				double dis_frac = DISPlotsTrue[WhichSample][WhichPlot]->Integral() / hstack->Integral() * 100.;
				TLegendEntry* ldis = ilegmc->AddEntry(DISPlotsTrue[WhichSample][WhichPlot],"DIS (" + to_string_with_precision(dis_frac,1.) + "%)","f");
				ldis->SetTextColor(kRed+1);	

				double coh_frac = COHPlotsTrue[WhichSample][WhichPlot]->Integral() / hstack->Integral() * 100.;
				TLegendEntry* lcoh = ilegmc->AddEntry(COHPlotsTrue[WhichSample][WhichPlot],"COH (" + to_string_with_precision(coh_frac,1.) + "%)","f");
				lcoh->SetTextColor(kMagenta);					

				textSlice->DrawLatexNDC(0.17, 0.92, LatexLabel[ ReducedPlotName ] );
				ilegmc->Draw();

				//----------------------------------------//

				// Save canvas with interaction breakdown
	
				gPad->RedrawAxis();
				inte_can->SaveAs("/exp/uboone/data/users/"+UserID+"/ncpi0/myPlots/flat_ttree/intebreak_" + NameOfSamples[WhichSample] + "_XSections_"+PlotNames[WhichPlot]+"_"+xsec_Runs[WhichRun]+".pdf");
				delete inte_can;

			} // end of the loop over the mc samples

			//----------------------------------------//

		} // End of the loop over the plots

		//----------------------------------------//					

	} // End of the loop over the runs	

} // End of the program