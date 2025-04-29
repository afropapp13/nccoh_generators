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

	// export paths

	TString preselection_file_path = "/exp/uboone/data/users/"+UserID+"/ncpi0";

	//----------------------------------------//

	// constants
	// Argon 

	static const double A = 40.;
	static const double Z = 18.;

	const int FontStyle = 132;
	const double TextSize = 0.07;
	const int NCont = 999; 

	double CosmicPID = -99.;
	
	//----------------------------------------//

	double neutron_ke_thres = 0.02; //GeV
	double proton_ke_thres = 0.02; //GeV
	double pi0_costheta_thres = 0.5; //GeV	

	// ------------------------------------ //

	// Run 1 

	static const double tor860_wcut_Run1 = 4.566e+19;
	static const double E1DCNT_wcut_Run1 = 10127594.0;
	static const double EXT_Run1 = 32514217.0;

	static const double Fulltor860_wcut_Run1 = 1.67e20;
	static const double FullE1DCNT_wcut_Run1 = 37273255.0;
	static const double FullEXT_Run1 = 65744587.0;
	
	// ------------------------------------ //

	// Run 1A Open Trigger 

	static const double tor860_wcut_Run1A_open_trigger = 5.885e+19;
	static const double E1DCNT_wcut_Run1A_open_trigger = 15639842.0;
	static const double EXT_Run1A_open_trigger = 65744587.0;

	static const double Fulltor860_wcut_Run1A_open_trigger = 5.885e+19;
	static const double FullE1DCNT_wcut_Run1A_open_trigger = 15639842.0;
	static const double FullEXT_Run1A_open_trigger = 65744587.0;
	
	// ------------------------------------ //

	// Run 1B Open Trigger 

	static const double tor860_wcut_Run1B_open_trigger = 8.932e+19;
	static const double E1DCNT_wcut_Run1B_open_trigger = 19789737.0;
	static const double EXT_Run1B_open_trigger = 65744587.0;

	static const double Fulltor860_wcut_Run1B_open_trigger = 8.932e+19;
	static const double FullE1DCNT_wcut_Run1B_open_trigger = 19789737.0;
	static const double FullEXT_Run1B_open_trigger = 65744587.0;
	
	// ------------------------------------ //	
	
	// Run 2 

	static const double tor860_wcut_Run2 = 0.;
	static const double E1DCNT_wcut_Run2 = 0.;
	static const double EXT_Run2 = 1.;

	static const double Fulltor860_wcut_Run2 = 2.61e20;
	static const double FullE1DCNT_wcut_Run2 = 61882791.0;
	static const double FullEXT_Run2 = 153905891.0;
	
	// ------------------------------------ //	
	
	// Run 3 

	static const double tor860_wcut_Run3 = 9.513e+18;
	static const double E1DCNT_wcut_Run3 = 2299517.0;
	static const double EXT_Run3 = 131356320.0;

	static const double Fulltor860_wcut_Run3 = 2.57e20;
	static const double FullE1DCNT_wcut_Run3 = 61326657.0;
	static const double FullEXT_Run3 = 207224641.0;
	
	// ------------------------------------ //	
	
	// Run 4a 

	static const double tor860_wcut_Run4a = 0.;
	static const double E1DCNT_wcut_Run4a = 0.;
	static const double EXT_Run4a = 1.;

	static const double Fulltor860_wcut_Run4a = 4.5e19;
	static const double FullE1DCNT_wcut_Run4a = 9897624.;
	static const double FullEXT_Run4a = 27596585.;	

	// ------------------------------------ //	
	
	// Run 4b 

	static const double tor860_wcut_Run4b = 0.;
	static const double E1DCNT_wcut_Run4b = 0.;
	static const double EXT_Run4b = 1.;

	static const double Fulltor860_wcut_Run4b = 1.36e20;
	static const double FullE1DCNT_wcut_Run4b = 32305463.0;
	static const double FullEXT_Run4b = 89244940.0;	

	// ------------------------------------ //	
	
	// mcc9_10 Run 4b standalone

	static const double tor860_wcut_mcc9_10_Run4b_standalone = 0.;
	static const double E1DCNT_wcut_mcc9_10_Run4b_standalone = 0.;
	static const double EXT_mcc9_10_Run4b_standalone = 1.;

	// good run list applied
	// produced by production team
	//static const double Fulltor860_wcut_mcc9_10_Run4b_standalone = 3.92e+19;
	//static const double FullE1DCNT_wcut_mcc9_10_Run4b_standalone = 9515547.0;
	//static const double FullEXT_mcc9_10_Run4b_standalone = 27868945.0;	

	//// good run list applied and rse matched
	//// Afro filtered it so proceed with caution
	static const double Fulltor860_wcut_mcc9_10_Run4b_standalone = 3.776e+19;
	static const double FullE1DCNT_wcut_mcc9_10_Run4b_standalone = 9176822.0;
	static const double FullEXT_mcc9_10_Run4b_standalone = 27868945.0;		

	// ------------------------------------ //	
	
	// mcc9_10 Run 4b unified

	static const double tor860_wcut_mcc9_10_Run4b_unified = 0.;
	static const double E1DCNT_wcut_mcc9_10_Run4b_unified = 0.;
	static const double EXT_mcc9_10_Run4b_unified = 1.;

	// //Erin filtered them and added her vars and good run list applied
	//static const double Fulltor860_wcut_mcc9_10_Run4b_unified = 4.28e19;
	//static const double FullE1DCNT_wcut_mcc9_10_Run4b_unified = 10398793.;
	//static const double FullEXT_mcc9_10_Run4b_unified = 28964045.;	
	
	//// Afro rse's the WC processed files from Erin
	static const double Fulltor860_wcut_mcc9_10_Run4b_unified = 3.776e+19;
	static const double FullE1DCNT_wcut_mcc9_10_Run4b_unified = 9176822.0;
	static const double FullEXT_mcc9_10_Run4b_unified = 28964045.;	


	// ------------------------------------ //
		
	// Run 4c

	static const double tor860_wcut_Run4c = 0.;
	static const double E1DCNT_wcut_Run4c = 0.;
	static const double EXT_Run4c = 1.;

	static const double Fulltor860_wcut_Run4c = 8.95e19;
	static const double FullE1DCNT_wcut_Run4c = 20273291.0;
	static const double FullEXT_Run4c = 47178301.0;	

	// ------------------------------------ //

	// Run 4d 

	static const double tor860_wcut_Run4d = 0.;
	static const double E1DCNT_wcut_Run4d = 0.;
	static const double EXT_Run4d = 1.;

	static const double Fulltor860_wcut_Run4d = 4.93e19;
	static const double FullE1DCNT_wcut_Run4d = 11192660.0;
	static const double FullEXT_Run4d = 74409530.0;	

	// ------------------------------------ //
		
	// Run 5 

	static const double tor860_wcut_Run5 = 0.;
	static const double E1DCNT_wcut_Run5 = 0.;
	static const double EXT_Run5 = 1.;

	static const double Fulltor860_wcut_Run5 = 1.48e20;
	static const double FullE1DCNT_wcut_Run5 = 35265730.0;
	static const double FullEXT_Run5 = 107466402.0;	

	// ------------------------------------ //	
	
	// Combined POT

	static const double Fulltor860_wcut_Run1all = Fulltor860_wcut_Run1 + Fulltor860_wcut_Run1A_open_trigger + Fulltor860_wcut_Run1B_open_trigger;	
	static const double Fulltor860_wcut_Run4 = Fulltor860_wcut_Run4a + Fulltor860_wcut_Run4b + Fulltor860_wcut_Run4c + Fulltor860_wcut_Run4d;	
	static const double Fulltor860_wcut_Combined = Fulltor860_wcut_Run1all + Fulltor860_wcut_Run2 + Fulltor860_wcut_Run3 + Fulltor860_wcut_Run4 + Fulltor860_wcut_Run5;	
	
	static const double FullEXT_Run1all = FullEXT_Run1 + FullEXT_Run1A_open_trigger + FullEXT_Run1B_open_trigger;
	static const double FullEXT_Run4 = FullEXT_Run4a + FullEXT_Run4b + FullEXT_Run4c + FullEXT_Run4d;
	static const double FullEXT_Combined = FullEXT_Run1all + FullEXT_Run2 + FullEXT_Run3 + FullEXT_Run4 + FullEXT_Run5;	

	static const double FullE1DCNT_wcut_Run1all = FullE1DCNT_wcut_Run1 + FullE1DCNT_wcut_Run1A_open_trigger + FullE1DCNT_wcut_Run1B_open_trigger;	
	static const double FullE1DCNT_wcut_Run4 =  FullE1DCNT_wcut_Run4a + FullE1DCNT_wcut_Run4b + FullE1DCNT_wcut_Run4c + FullE1DCNT_wcut_Run4d;	
	static const double FullE1DCNT_wcut_Combined = FullE1DCNT_wcut_Run1all + FullE1DCNT_wcut_Run2 + FullE1DCNT_wcut_Run3 + FullE1DCNT_wcut_Run4 + FullE1DCNT_wcut_Run5;
	
	//----------------------------------------//	

	// Labels / Ranges & Label  map
	// max values

	static std::map<TString,std::pair<double,double> > XSecRange =
	{

		{ "SingleBinPlot",  std::make_pair(0, 16.9) },
		{ "Pi0CosThetaPlot",  std::make_pair(0, 109.) },		
		{ "Pi0MomentumPlot",  std::make_pair(0, 34.) },		
		{ "Pi0InvMassPlot",  std::make_pair(0, 34.) },		

	};	
	
	//----------------------------------------//

	static std::map<TString,TString> VarLabel =
	{

		{ "SingleBinPlot",  "#sigma #left[10^{-38} #frac{cm^{2}}{Ar}#right]" },
		{ "Pi0CosThetaPlot",  "#frac{d#sigma}{dcos#theta_{#pi^{ 0}}} #left[10^{-38} #frac{cm^{2}}{Ar}#right]" },
		{ "Pi0MomentumPlot",  "#frac{d#sigma}{dp_{#pi^{ 0}}} #left[10^{-38} #frac{cm^{2}}{(GeV/c) Ar}#right]" },
		{ "Pi0InvMassPlot",  "#frac{d#sigma}{dp_{#pi^{ 0}}} #left[10^{-38} #frac{cm^{2}}{(GeV/c^{2}) Ar}#right]" },

	};

	static std::map<TString,TString> LatexLabel =
	{

		{ "MuonCosThetaSingleBinPlot",  "all events" },
		{ "Pi0CosThetaPlot", "all events" },	
		{ "Pi0MomentumPlot", "all events" },	
		{ "Pi0InvMassPlot", "all events" },	
	
	};	
	
	//----------------------------------------//
	
	static std::map<TString,TString> MapUncorCor =
	{

		{ "SingleBinPlot", "SingleBinPlot" },
		{ "Pi0CosThetaPlot", "Pi0CosThetaPlot" },
		{ "Pi0MomentumPlot", "Pi0MomentumPlot" },
		{ "Pi0InvMassPlot", "Pi0InvMassPlot" },
	
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
				 ,"Pi0InvMassPlot"
			 
	};
	
	//----------------------------------------//
	
	vector<TString> OneDimXSec = {
				 "SingleBinPlot"
				 ,"Pi0CosThetaPlot"
				 ,"Pi0MomentumPlot"
				 ,"Pi0InvMassPlot"
				 				 
	};	
	
	//----------------------------------------//

	// Binning

	static const int NBinsSingleBin = 1; static const double ArrayNBinsSingleBin[NBinsSingleBin+1] = {0.,1.};
	
	static const int NBinsPi0CosTheta = 18;
	static const double ArrayNBinsPi0CosTheta[NBinsPi0CosTheta+1] = { -1.,-0.85,-0.7,-0.57,-0.45,-0.32,-0.2,-0.1,0.,0.1,0.2,0.3,0.4,0.5,0.6,0.72,0.84,0.95,1.};

	static const int NBinsPi0Momentum = 15; static const double ArrayNBinsPi0Momentum[NBinsPi0Momentum+1] = {0.,0.1,0.2,0.3,0.38,0.45,0.5,0.55,0.625,0.7,0.75,0.8,0.87,1.,1.1,1.2};

	static const int NBinsPi0InvMass = 15; static const double ArrayNBinsPi0InvMass[NBinsPi0InvMass+1] = {0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.,1.1,1.2,1.3,1.4,1.5};

	//----------------------------------------//

	// Labels for 1D plots
	
	static TString LabelXAxisSingleBin = ";single bin"; static TString LabelXAxisTrueSingleBin = ";true single bin";
	static TString LabelXAxisPi0CosTheta = ";cos#theta_{#pi^{ 0}}"; static TString LabelXAxisTruePi0CosTheta = ";true cos#theta_{#pi^{ 0}}";
	static TString LabelXAxisPi0Momentum = ";#pi^{ 0} momentum [GeV/c]"; static TString LabelXAxisTruePi0Momentum = ";true #pi^{0} momentum [GeV/c]";
	static TString LabelXAxisPi0InvMass= ";#pi^{ 0} invariant mass [GeV/c^{2}]"; static TString LabelXAxisTruePi0InvMass = ";true #pi^{0} invariant mass [GeV/c^{2}]";
	
	//----------------------------------------//

	// Labels for 2D Plots

	static TString LabelXAxisSingleBin2D = LabelXAxisTrueSingleBin+";reco single bin";
	static TString LabelXAxisPi0CosTheta2D = LabelXAxisTruePi0CosTheta+";reco cos#theta_{#pi^{ 0}}";
	static TString LabelXAxisPi0Momentum2D = LabelXAxisTruePi0Momentum+";reco #pi^{ 0} momentum [GeV/c]";
	static TString LabelXAxisPi0InvMass2D = LabelXAxisTruePi0InvMass+";reco #pi^{ 0} invariant mass [GeV/c^{2}]";

	//----------------------------------------//
	
	// Interaction labels
	
	const std::vector<int> InteBreakColors{kBlack,kAzure-4,kOrange-3,kGreen+1,kRed+1,kBlue};		
	std::vector<TString> InteractionLabels = {"","QE","MEC","RES","DIS","COH"};
	const int NInte = InteractionLabels.size();
	
	//----------------------------------------//		


}

#endif
