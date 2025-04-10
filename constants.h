#ifndef CONSTANTS_H
#define CONSTANTS_H

#include "TString.h"
#include "TMath.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <map>

using namespace std;

namespace constants {

	//----------------------------------------//

	// Kerberos user name
  
	TString UserID = getenv("USER");

	//----------------------------------------//

	// constants
	// Argon 

	static const double A = 40.;
	static const double Z = 18.;

	const int FontStyle = 132;
	const double TextSize = 0.07;
	const int NCont = 999; 
	
	//----------------------------------------//

	// Labels / Ranges & Label  map
	// max values

	static std::map<TString,std::pair<double,double> > XSecRange =
	{

		{ "SingleBinPlot",  std::make_pair(0, 16.9) },
		{ "Pi0CosThetaPlot",  std::make_pair(0, 109.) },		
		{ "Pi0MomentumPlot",  std::make_pair(0, 34.) },		

	};	
	
	//----------------------------------------//

	static std::map<TString,TString> VarLabel =
	{

		{ "SingleBinPlot",  "#frac{d#sigma}{dcos#theta_{#mu}} #left[10^{-38} #frac{cm^{2}}{Ar}#right]" },
		{ "Pi0CosThetaPlot",  "#frac{d#sigma}{dcos#theta_{vis}} #left[10^{-38} #frac{cm^{2}}{Ar}#right]" },
		{ "Pi0MomentumPlot",  "#frac{d#sigma}{dp_{miss}} #left[10^{-38} #frac{cm^{2}}{(GeV/c) Ar}#right]" },

	};

	static std::map<TString,TString> LatexLabel =
	{

		{ "MuonCosThetaSingleBinPlot",  "All events" },
		{ "Pi0CosThetaPlot", "All events" },	
		{ "Pi0MomentumPlot", "All events" },	
	
	};	
	
	//----------------------------------------//
	
	static std::map<TString,TString> MapUncorCor =
	{

		{ "SingleBinPlot", "SingleBinPlot" },
		{ "Pi0CosThetaPlot", "Pi0CosThetaPlot" },
		{ "Pi0MomentumPlot", "Pi0MomentumPlot" },
	
	};
					
	//----------------------------------------//

	// Global Constants

	static const double Units = 1E38; // so that the extracted cross-section is in 10^{-38} cm^{2}
	
	static const int NuMuPdg = 14, MuonPdg = 13, ProtonPdg = 2212, AbsChargedPionPdg = 211, NeutralPionPdg = 111;
	static const int ElectronPdg = 11, nue_pdg = 12, PhotonPdg = 22, NeutronPdg = 2112, KaonPdg = 321, NeutralKaonPdg = 311;
	static const int SigmaPlusPdg = 3222, SigmaMinusPdg = 3112, NeutralSigmaPdg = 3212, LambdaPdg = 3122;
	static const int DeuteriumPdg = 1000010020, HeliumPdg = 1000020040, ArgonPdg = 1000180400;
	static const int rho_pdg = 113, charged_rho_pdg = 213, eta_pdg = 221;	
	static const int hydrogen_cluster_pdg = 2000000101, nucleon_pair = 2000000201;
	static const int proton_pair = 2000000202, neutron_pair = 2000000200;

	static const double MuonMass = 106, ProtonMass = 938.272, NeutronMass = 939.565; // MeV
	static const double MuonMass_GeV = 0.106, ProtonMass_GeV = 0.938272, NeutronMass_GeV = 0.939565; // GeV
	static const double DeltaM2 = TMath::Power(NeutronMass_GeV,2.) - TMath::Power(ProtonMass_GeV,2.); // GeV^2	

	//----------------------------------------//

	// Plots to be included for xsec extraction purposes

	vector<TString> PlotNames{
				 "SingleBinPlot"
				 ,"Pi0CosThetaPlot"
				 ,"Pi0MomentumPlot"
			 
	};
	
	//----------------------------------------//
	
	vector<TString> OneDimXSec = {
				 "SingleBinPlot"
				 ,"Pi0CosThetaPlot"
				 ,"Pi0MomentumPlot"
				 				 
	};	
	
	//----------------------------------------//

	// Binning

	static const int NBinsSingleBin = 1; static const double ArrayNBinsSingleBin[NBinsSingleBin+1] = {0.,1.};
	
	static const int NBinsPi0CosTheta = 18;
	static const double ArrayNBinsPi0CosTheta[NBinsPi0CosTheta+1] = { -1.,-0.85,-0.7,-0.57,-0.45,-0.32,-0.2,-0.1,0.,0.1,0.2,0.3,0.4,0.5,0.6,0.72,0.84,0.95,1.};

	static const int NBinsPi0Momentum = 13; static const double ArrayNBinsPi0Momentum[NBinsPi0Momentum+1] = {0.,0.1,0.2,0.3,0.38,0.45,0.5,0.55,0.625,0.7,0.75,0.8,0.87,1.};

	//----------------------------------------//

	// Labels for 1D plots
	
	static TString LabelXAxisSingleBin = ";single bin"; static TString LabelXAxisTrueSingleBin = ";true single bin";
	static TString LabelXAxisPi0CosTheta = ";cos#theta_{#pi^{0}}"; static TString LabelXAxisTruePi0CosTheta = ";true cos#theta_{#pi^{0}}";
	static TString LabelXAxisPi0Momentum = ";#p^{0} momentum [GeV/c]"; static TString LabelXAxisTruePi0Momentum = ";true #pi^{0} momentum [GeV/c]";
	
	//----------------------------------------//

	// Labels for 2D Plots

	static TString LabelXAxisSingleBin2D = LabelXAxisTrueSingleBin+";reco single bin";
	static TString LabelXAxisPi0CosTheta2D = LabelXAxisTruePi0CosTheta+";reco cos#theta_{#pi^{0}}";
	static TString LabelXAxisPi0Momentum2D = LabelXAxisTruePi0Momentum+";reco #pi^{0} momentum [GeV/c]";

	//----------------------------------------//
	
	// Interaction labels
	
	const std::vector<int> InteBreakColors{kBlack,kAzure-4,kOrange-3,kGreen+1,kRed+1,kBlue};		
	std::vector<TString> InteractionLabels = {"","QE","MEC","RES","DIS","COH"};
	const int NInte = InteractionLabels.size();
	
	//----------------------------------------//		


}
#endif
