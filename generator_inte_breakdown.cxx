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
#include "helper_functions.cxx"

using namespace std;
using namespace constants;

void generator_inte_breakdown() {

	//------------------------------//

	TH1D::SetDefaultSumw2();
	gStyle->SetOptStat(0);			

	TString OutFilePath = "output_files/";

	//------------------------------//

	// Event generators

	std::vector<TString> Names; std::vector<TString> Labels; 

	std::vector<int> Colors{kBlack,kBlue,kRed+1,kOrange+7,kGreen+1, kMagenta+1};
	std::vector<TString> Process{"","QE","MEC","RES","DIS","COH"};

	Names.push_back(OutFilePath+"analyzer_ncpi0_GENIE_v3_0_6.root"); Labels.push_back("G18");
	Names.push_back(OutFilePath+"analyzer_ncpi0_GENIE_v3_6_0_RS.root"); Labels.push_back("G18_RS");	
	Names.push_back(OutFilePath+"analyzer_ncpi0_NEUT_5_6_0.root"); Labels.push_back("NEUT");
	Names.push_back(OutFilePath+"analyzer_ncpi0_NuWro_25_03_1.root"); Labels.push_back("NuWro");
	Names.push_back(OutFilePath+"analyzer_ncpi0_GiBUU_2025.root"); Labels.push_back("GiBUU");			
	
	const int NSamples = Names.size();
	const int NColors = Colors.size();
	const int NProcesses = Process.size();

	std::vector<TFile*> Files; Files.resize(NSamples);

	//------------------------------//

	// Plots to overlay

	std::vector<TString> PlotNames;

	// 1D

	PlotNames.push_back("TrueSingleBinPlot");
	//PlotNames.push_back("TruePi0CosThetaPlot");
	PlotNames.push_back("TruePi0MomentumPlot");

	//------------------------------//

	const int NPlots = PlotNames.size();

	//------------------------------//	

	// Loop over the samples to open the files

	for (int iSample = 0; iSample < NSamples; iSample++) {

		Files[iSample] = new TFile(Names[iSample],"readonly");

	} // End of the loop over the samples

	//------------------------------//

	// Loop over the plots to be compared

	for (int iPlot = 0; iPlot < NPlots; iPlot++) {

		for (int iSample = 0; iSample < NSamples; iSample++) {			

			TString LabelCopy = Labels[iSample];
		  	TString CanvasName = "ncpi0_"+LabelCopy.ReplaceAll(" ","_")+"_InteBreakDown_" + PlotNames[iPlot];
		  	TCanvas* PlotCanvas = new TCanvas(CanvasName,CanvasName,205,34,1024,768);
		  	PlotCanvas->cd();
		  	PlotCanvas->SetTopMargin(0.14);
		  	PlotCanvas->SetLeftMargin(0.19);
		  	PlotCanvas->SetRightMargin(0.03);
		  	PlotCanvas->SetBottomMargin(0.16);		
		  	PlotCanvas->Draw();	

			TLegend* leg = new TLegend(0.15,0.87,0.99,0.99);
			leg->SetBorderSize(0);
			leg->SetNColumns(3);
			leg->SetTextSize(TextSize);	
			leg->SetTextFont(FontStyle);						
			leg->SetMargin(0.1);				
			leg->SetFillColor(0);				

			// Loop over the interaction processes

		 	std::vector<TH1D*> Histos; Histos.resize(NProcesses);

		  	for (int iprocess = 0; iprocess < NProcesses; iprocess++) {	

				Histos[iprocess] = (TH1D*)(Files[iSample]->Get(Process[iprocess]+PlotNames[iPlot]));

				Histos[iprocess]->SetLineWidth(4);
				Histos[iprocess]->SetLineColor( Colors.at(iprocess) );	

				Histos[iprocess]->GetXaxis()->SetTitleFont(FontStyle);
				Histos[iprocess]->GetXaxis()->SetLabelFont(FontStyle);
				Histos[iprocess]->GetXaxis()->SetNdivisions(8);
				Histos[iprocess]->GetXaxis()->SetLabelSize(TextSize);
				Histos[iprocess]->GetXaxis()->SetTitleSize(TextSize);	
				Histos[iprocess]->GetXaxis()->SetTitleOffset(1.1);					
				Histos[iprocess]->GetXaxis()->CenterTitle();						

				Histos[iprocess]->GetYaxis()->SetTitleFont(FontStyle);
				Histos[iprocess]->GetYaxis()->SetLabelFont(FontStyle);
				Histos[iprocess]->GetYaxis()->SetNdivisions(6);
				Histos[iprocess]->GetYaxis()->SetLabelSize(TextSize);
				TString reduced_plot_name = PlotNames[iPlot];
				reduced_plot_name.ReplaceAll("True","");
				Histos[iprocess]->GetYaxis()->SetTitle( VarLabel[reduced_plot_name] );
				Histos[iprocess]->GetYaxis()->SetTitleSize(TextSize);
				Histos[iprocess]->GetYaxis()->SetTitleOffset(1.2);
				Histos[iprocess]->GetYaxis()->CenterTitle();	
				Histos[iprocess]->GetYaxis()->SetRangeUser(0.,1.15*Histos[0]->GetMaximum());

				Histos[iprocess]->Draw("hist same");
				Histos[0]->Draw("hist same");	

				double frac = Histos[iprocess]->Integral("width")/Histos[0]->Integral("width") * 100.;
				TString LegLabel = Process[iprocess] + " (" + to_string_with_precision(frac,1) + "%)";
				if (iprocess == 0) { LegLabel = "Total (" + to_string_with_precision(frac,1) + "%)"; }
				TLegendEntry* legColor = leg->AddEntry(Histos[iprocess],LegLabel,"l");
				legColor->SetTextColor( Colors.at(iprocess) ); 

			} // End of the loop over the processes

			PlotCanvas->cd();
			leg->Draw();

			TLatex *textSlice = new TLatex();
			textSlice->SetTextFont(FontStyle);
			textSlice->SetTextSize(TextSize);
			textSlice->DrawLatexNDC(0.23, 0.81, Labels[iSample] + "      " + LatexLabel[PlotNames[iPlot]]);

			gPad->RedrawAxis();

			PlotCanvas->SaveAs("/exp/uboone/data/users/"+UserID+"/ncpi0/myPlots/flat_ttree/"+CanvasName+".pdf");
			delete PlotCanvas;

		} // End of the loop over the samples grabing the plots	

	} // End of the loop over the plots

} // End of the program