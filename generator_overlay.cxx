#include <TFile.h>
#include <TString.h>
#include <TH1D.h>
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
	gStyle->SetOptStat(0);		

	TString OutFilePath = "output_files/";

	//------------------------------//

	// Event generators

	std::vector<TString> Names; std::vector<TString> Labels; 

	vector<int> Colors;  vector<int> LineStyle;

	if (Tag == "") {

		Names.push_back(OutFilePath+"analyzer_ncpi0_GENIE_v3_0_6.root"); Labels.push_back("G18"); Colors.push_back(kBlack); LineStyle.push_back(kSolid);
		Names.push_back(OutFilePath+"analyzer_ncpi0_GENIE_v3_6_0_AR23.root"); Labels.push_back("AR23"); Colors.push_back(kRed+1); LineStyle.push_back(kSolid);
		Names.push_back(OutFilePath+"analyzer_ncpi0_GENIE_v3_6_0_RS.root"); Labels.push_back("G18_RS"); Colors.push_back(kMagenta+1); LineStyle.push_back(kSolid);		
		Names.push_back(OutFilePath+"analyzer_ncpi0_NEUT_5_6_0.root"); Labels.push_back("NEUT"); Colors.push_back(kGreen+1); LineStyle.push_back(kSolid); // kMagenta - 9
		Names.push_back(OutFilePath+"analyzer_ncpi0_NuWro_25_03_1.root"); Labels.push_back("NuWro"); Colors.push_back(kOrange+7); LineStyle.push_back(kSolid);
		Names.push_back(OutFilePath+"analyzer_ncpi0_GiBUU_2025.root"); Labels.push_back("GiBUU"); Colors.push_back(kAzure+7); LineStyle.push_back(kSolid);

	}
	
	const int NSamples = Names.size();
	const int NColors = Colors.size();

	std::vector<TFile*> Files; Files.resize(NSamples);

	//------------------------------//

	// Plots to overlay

	std::vector<TString> PlotNames;

	// 1D

	PlotNames.push_back("TrueSingleBinPlot"); 
	//PlotNames.push_back("TruePi0CosThetaPlot");
	PlotNames.push_back("TruePi0MomentumPlot"); 
	
	const int NPlots = PlotNames.size();

	//------------------------------//	

	// Loop over the samples to open the files

	for (int iSample = 0; iSample < NSamples; iSample++) {

		Files[iSample] = new TFile(Names[iSample],"readonly");

	} // End of the loop over the samples

	//------------------------------//

	// Loop over the plots to be compared

	for (int iPlot = 0; iPlot < NPlots; iPlot++) {

		TString CanvasName = "ncpi0_generator_overlay_" + PlotNames[iPlot];
		TCanvas* PlotCanvas = new TCanvas(CanvasName,CanvasName,205,34,1024,768);
		PlotCanvas->cd();
		PlotCanvas->SetTopMargin(0.15);
		PlotCanvas->SetLeftMargin(0.19);
		PlotCanvas->SetRightMargin(0.04);
		PlotCanvas->SetBottomMargin(0.16);		
		PlotCanvas->Draw();	

		TLegend* leg = new TLegend(0.17,0.865,0.99,0.985);
		leg->SetBorderSize(0);
		leg->SetNColumns(3);
		leg->SetTextSize(TextSize);	
		leg->SetTextFont(FontStyle);						
		leg->SetMargin(0.3);						

		// Loop over the samples to open the files and to get the corresponding plot

		std::vector<TH1D*> Histos; Histos.resize(NSamples);

		for (int iSample = 0; iSample < NSamples; iSample++) {	

		    Histos[iSample] = (TH1D*)(Files[iSample]->Get(PlotNames[iPlot]));

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
			TString reduced_plot_name = PlotNames[iPlot];
			reduced_plot_name.ReplaceAll("True","");
			Histos[iSample]->GetYaxis()->SetTitle( VarLabel[reduced_plot_name] );
			Histos[iSample]->GetYaxis()->SetTitleSize(TextSize);
			Histos[iSample]->GetYaxis()->SetTitleOffset(1.2);
			Histos[iSample]->GetYaxis()->CenterTitle();	

			double imax = TMath::Max(Histos[iSample]->GetMaximum(),Histos[0]->GetMaximum());
			double MaxScale = 1.05;
			Histos[iSample]->GetYaxis()->SetRangeUser(0.,MaxScale*imax);
			Histos[0]->GetYaxis()->SetRangeUser(0.,MaxScale*imax);			

			PlotCanvas->cd();
			Histos[iSample]->Draw("hist same");

			TLegendEntry* legColor = leg->AddEntry(Histos[iSample],Labels[iSample],"l");
			legColor->SetTextColor( Colors.at(iSample) ); 

			//----------------------------------------//					

		} // End of the loop over the samples grabing the plots	

		PlotCanvas->cd();
		leg->Draw();

		PlotCanvas->SaveAs("/exp/uboone/data/users/"+UserID+"/ncpi0/myPlots/"+Tag+CanvasName+".pdf");
		delete PlotCanvas;

	} // End of the loop over the plots

	//------------------------------//

} // End of the program