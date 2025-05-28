#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TLatex.h>

#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <fstream>
#include <stdlib.h>

#include "constants.h"

using namespace std;
using namespace constants;

void generator_overlay(TString Tag = "") {

	//------------------------------//

	TH1D::SetDefaultSumw2();
	TH2D::SetDefaultSumw2();
	gStyle->SetOptStat(0);		

	TString OutFilePath = "output_files/";

	//------------------------------//

	// Event generators

	std::vector<TString> Names; std::vector<TString> Labels; 

	vector<int> Colors;  vector<int> LineStyle;

	if (Tag == "") {

	  Names.push_back(OutFilePath+"analyzer_ncpi0_GENIE_v3_0_6.root"); Labels.push_back("G18"); Colors.push_back(kBlack); LineStyle.push_back(kSolid);
	  //Names.push_back(OutFilePath+"analyzer_ncpi0_GENIE_v3_4_0_AR23_20i_00_000.root"); Labels.push_back("AR23"); Colors.push_back(kBlue); LineStyle.push_back(Gv2LineStyle);
	  Names.push_back(OutFilePath+"analyzer_ncpi0_NEUT_5_4_0_1.root"); Labels.push_back("NEUT"); Colors.push_back(kMagenta); LineStyle.push_back(kSolid); // kMagenta - 9
	  //Names.push_back(OutFilePath+"analyzer_ncpi0_NuWro_21_09_02.root"); Labels.push_back("NuWro"); Colors.push_back(NEUTColor); LineStyle.push_back(NuWroLineStyle);
	  //Names.push_back(OutFilePath+"analyzer_ncpi0_GiBUU_2023.root"); Labels.push_back("GiBUU"); Colors.push_back(GiBUUColor); LineStyle.push_back(GiBUULineStyle);

	}
	
	const int NSamples = Names.size();
	const int NColors = Colors.size();

	// Sanity check
	if (NColors < NSamples) { cout << "Give me some more colors!" << endl; return; }

	std::vector<TFile*> Files; Files.resize(NSamples);

	//------------------------------//

	// Plots to overlay

	std::vector<TString> PlotNames;
	std::vector<TString> YAxisLabel;

	//------------------------------//

	// Post FSI

	// 1D

	PlotNames.push_back("TruePi0CosThetaPlot"); YAxisLabel.push_back("#frac{d#sigma}{dcos#theta_{#pi^{0}}} #left[10^{-38} #frac{cm^{2}}{Ar}#right]");
	PlotNames.push_back("TruePi0MomentumPlot"); YAxisLabel.push_back("#frac{d#sigma}{dp_{#pi^{0}}}  #left[10^{-38} #frac{cm^{2}}{(GeV/c) Ar}#right]");
	PlotNames.push_back("TrueSingleBinPlot"); YAxisLabel.push_back("sigma #left[10^{-38} #frac{cm^{2}}{Ar}#right]");

	//------------------------------//

	const int NPlots = PlotNames.size();
	const int NLabels = YAxisLabel.size();

	// sanity check
	if ( NPlots != NLabels) { cout << "Inconsistent number of plots and labels! Aborting !" << endl; return; }

	//------------------------------//	

	// Loop over the samples to open the files and the TTree

	for (int iSample = 0; iSample < NSamples; iSample++) {

		Files[iSample] = new TFile(Names[iSample],"readonly");

	} // End of the loop over the samples

	//------------------------------//

	// Loop over the plots to be compared

	for (int iPlot = 0; iPlot < NPlots; iPlot++) {

		TString CanvasName = "ThreeDKI_GeneratorOverlay_" + PlotNames[iPlot];
		TCanvas* PlotCanvas = new TCanvas(CanvasName,CanvasName,205,34,1024,768);
		PlotCanvas->cd();
		PlotCanvas->SetTopMargin(0.15);
		PlotCanvas->SetLeftMargin(0.17);
		PlotCanvas->SetRightMargin(0.05);
		PlotCanvas->SetBottomMargin(0.16);		
		PlotCanvas->Draw();	

		TLegend* leg = new TLegend(0.17,0.855,0.99,0.985);
		if (Tag == "Honda" && PlotNames[iPlot] == "TrueFineBinEvPlot") { leg = new TLegend(0.75,0.7,0.9,0.8);  }
		leg->SetBorderSize(0);
		leg->SetNColumns(6);
		if (Tag == "Honda" && PlotNames[iPlot] == "TrueFineBinEvPlot") { leg->SetNColumns(1); }
		leg->SetTextSize(TextSize);	
		leg->SetTextFont(FontStyle);						
		leg->SetMargin(0.3);						

		// Loop over the samples to open the files and to get the corresponding plot

		std::vector<TH1D*> Histos; Histos.resize(NSamples);

		for (int iSample = 0; iSample < NSamples; iSample++) {	

		        Histos[iSample] = (TH1D*)(Files[iSample]->Get(PlotNames[iPlot]));

			//Histos[iSample]->Scale(Histos[0]->Integral() / Histos[iSample]->Integral());

			Histos[iSample]->SetLineWidth(3);
			Histos[iSample]->SetLineColor( Colors.at(iSample) );	
			Histos[iSample]->SetLineStyle( LineStyle.at(iSample) );	

			Histos[iSample]->GetXaxis()->SetTitleFont(FontStyle);
			Histos[iSample]->GetXaxis()->SetLabelFont(FontStyle);
			Histos[iSample]->GetXaxis()->SetNdivisions(8);
			Histos[iSample]->GetXaxis()->SetLabelSize(TextSize);
			Histos[iSample]->GetXaxis()->SetTitleSize(TextSize);	
			Histos[iSample]->GetXaxis()->SetTitleOffset(1.1);					
			Histos[iSample]->GetXaxis()->SetLabelOffset(0.01);					
			Histos[iSample]->GetXaxis()->CenterTitle();						

			Histos[iSample]->GetYaxis()->SetTitleFont(FontStyle);
			Histos[iSample]->GetYaxis()->SetLabelFont(FontStyle);
			Histos[iSample]->GetYaxis()->SetNdivisions(6);
			Histos[iSample]->GetYaxis()->SetLabelSize(TextSize);
			Histos[iSample]->GetYaxis()->SetTitle(YAxisLabel.at(iPlot));
			Histos[iSample]->GetYaxis()->SetTitleSize(TextSize);
			Histos[iSample]->GetYaxis()->SetTitleOffset(1.2);
			//Histos[iSample]->GetYaxis()->SetTickSize(0);
			Histos[iSample]->GetYaxis()->CenterTitle();	

			double imax = TMath::Max(Histos[iSample]->GetMaximum(),Histos[0]->GetMaximum());
			double MaxScale = 1.05;
			Histos[iSample]->GetYaxis()->SetRangeUser(0.,MaxScale*imax);
			Histos[0]->GetYaxis()->SetRangeUser(0.,MaxScale*imax);			

			PlotCanvas->cd();
			Histos[iSample]->Draw("hist same");
			//Histos[0]->Draw("hist same");	

			TLegendEntry* legColor = leg->AddEntry(Histos[iSample],Labels[iSample],"l");
			legColor->SetTextColor( Colors.at(iSample) ); 

			TLatex *textSlice = new TLatex();
			textSlice->SetTextFont(FontStyle);
			textSlice->SetTextSize(TextSize);
			TString PlotNameDuplicate = PlotNames[iPlot];
			TString ReducedPlotName = PlotNameDuplicate.ReplaceAll("TrueFineBin","") ;
			if (Tag == "Honda" && PlotNames[iPlot] == "TrueFineBinEvPlot") {

				textSlice->DrawLatexNDC(0.3, 0.9, "Flux-averaged cross section");
			
			} else {

				textSlice->DrawLatexNDC(0.2, 0.8, LatexLabel[ReducedPlotName].ReplaceAll("All events",""));

			}

			//----------------------------------------//					

		} // End of the loop over the samples grabing the plots	

		PlotCanvas->cd();
		leg->Draw();

		PlotCanvas->SaveAs("/exp/uboone/data/users/"+UserID+"/PeLEETuples_Atmospherics/FlatTTreePlots/"+Tag+CanvasName+".pdf");
		delete PlotCanvas;

	} // End of the loop over the plots

	//------------------------------//

} // End of the program
