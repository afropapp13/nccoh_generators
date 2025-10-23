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

	TString preselection_file_path = "/exp/uboone/data/users/"+UserID+"/ncpi0/";
	TString event_selection_file_path = "/exp/uboone/data/users/"+UserID+"/ncpi0/output_files/";
	TString plot_path = "/exp/uboone/data/users/"+UserID+"/ncpi0/myPlots/"; 
	TString efficiency_path = "/exp/uboone/data/users/"+UserID+"/ncpi0/efficiencies/";
	TString migration_matrices_path = "/exp/uboone/data/users/"+UserID+"/ncpi0/migration_matrices/";	
	TString xsec_path = "/exp/uboone/data/users/"+UserID+"/ncpi0/xsec/";	

	//----------------------------------------//
	
	// Colors

	int BeamOnColor = kBlack;
	int OverlayColor = kAzure+7;

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

	// Zarko Pavlovic, Jun 22 2020
	// /pnfs/uboone/persistent/uboonebeam/bnb_gsimple/bnb_gsimple_fluxes_01.09.2019_463_hist/readme.txt

	double Nominal_UB_XY_Surface = 256.35*233.; // cm2
	double SoftFidSurface = 236. * 210.;  // cm2
	double POTPerSpill = 4997.*5e8;

	//----------------------------------------//

	double neutron_ke_thres = 0.02; //GeV
	double proton_ke_thres = 0.02; //GeV
	double pi0_costheta_thres = 0.85; //GeV	
	double gamma1_costheta_thres = 0.6; //GeV	
	double gamma2_costheta_thres = 0.4; //GeV		
	double TRACK_SCORE_CUT = 0.5;

	// ------------------------------------ //

	vector<TString> xsec_Runs = {"Run4b_unified"};

	TString CutExtension = "_nocuts";

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
	
	// Run 4a unified

	static const double tor860_wcut_Run4a_unified = 0.;
	static const double E1DCNT_wcut_Run4a_unified = 0.;
	static const double EXT_Run4a_unified = 1.;

	static const double Fulltor860_wcut_Run4a_unified = 4.5e19;
	static const double FullE1DCNT_wcut_Run4a_unified = 9897624.;
	static const double FullEXT_Run4a_unified = 27596585.;		

	// ------------------------------------ //	
	
	// Run 4b 

	static const double tor860_wcut_Run4b = 0.;
	static const double E1DCNT_wcut_Run4b = 0.;
	static const double EXT_Run4b = 1.;

	static const double Fulltor860_wcut_Run4b = 1.332e20;
	static const double FullE1DCNT_wcut_Run4b = 31582916;
	static const double FullEXT_Run4b = 88445969;	

	// ------------------------------------ //	
	
	// mcc9_10 Run 4b standalone

	static const double tor860_wcut_mcc9_10_Run4b_standalone = 0.;
	static const double E1DCNT_wcut_mcc9_10_Run4b_standalone = 0.;
	static const double EXT_mcc9_10_Run4b_standalone = 1.;	

	// Aug 12 
	static const double Fulltor860_wcut_mcc9_10_Run4b_standalone = 1.2e20;
	static const double FullE1DCNT_wcut_mcc9_10_Run4b_standalone = 28396891;
	static const double FullEXT_mcc9_10_Run4b_standalone = 46393391;	

	// ------------------------------------ //	
	
	// mcc9_10 Run 4b unified

	static const double tor860_wcut_mcc9_10_Run4b_unified = 0.;
	static const double E1DCNT_wcut_mcc9_10_Run4b_unified = 0.;
	static const double EXT_mcc9_10_Run4b_unified = 1.;

	// Aug 12
	static const double Fulltor860_wcut_mcc9_10_Run4b_unified = 1.332e+20;
	static const double FullE1DCNT_wcut_mcc9_10_Run4b_unified = 31582916.0;
	static const double FullEXT_mcc9_10_Run4b_unified = 88445969.0;

	// ------------------------------------ //
		
	// Run 4c

	static const double tor860_wcut_Run4c = 0.;
	static const double E1DCNT_wcut_Run4c = 0.;
	static const double EXT_Run4c = 1.;

	static const double Fulltor860_wcut_Run4c = 8.95e19;
	static const double FullE1DCNT_wcut_Run4c = 20273291.0;
	static const double FullEXT_Run4c = 47178301.0;	

	// ------------------------------------ //
		
	// Run 4c unified

	static const double tor860_wcut_Run4c_unified = 0.;
	static const double E1DCNT_wcut_Run4c_unified = 0.;
	static const double EXT_Run4c_unified = 1.;

	static const double Fulltor860_wcut_Run4c_unified = 8.95e19;
	static const double FullE1DCNT_wcut_Run4c_unified = 20273291.0;
	static const double FullEXT_Run4c_unified = 47178301.0;		

	// ------------------------------------ //

	// Run 4d 

	static const double tor860_wcut_Run4d = 0.;
	static const double E1DCNT_wcut_Run4d = 0.;
	static const double EXT_Run4d = 1.;

	static const double Fulltor860_wcut_Run4d = 4.93e19;
	static const double FullE1DCNT_wcut_Run4d = 11192660.0;
	static const double FullEXT_Run4d = 74409530.0;	

	// ------------------------------------ //

	// Run 4d unified

	static const double tor860_wcut_Run4d_unified = 0.;
	static const double E1DCNT_wcut_Run4d_unified = 0.;
	static const double EXT_Run4d_unified = 1.;

	static const double Fulltor860_wcut_Run4d_unified = 4.93e19;
	static const double FullE1DCNT_wcut_Run4d_unified = 11192660.0;
	static const double FullEXT_Run4d_unified = 74409530.0;		

	// ------------------------------------ //
		
	// Run 5 

	static const double tor860_wcut_Run5 = 0.;
	static const double E1DCNT_wcut_Run5 = 0.;
	static const double EXT_Run5 = 1.;

	static const double Fulltor860_wcut_Run5 = 1.48e20;
	static const double FullE1DCNT_wcut_Run5 = 35265730.0;
	static const double FullEXT_Run5 = 107466402.0;	

	// ------------------------------------ //
		
	// Run 5 unified

	static const double tor860_wcut_Run5_unified = 0.;
	static const double E1DCNT_wcut_Run5_unified = 0.;
	static const double EXT_Run5_unified = 1.;

	static const double Fulltor860_wcut_Run5_unified = 1.48e20;
	static const double FullE1DCNT_wcut_Run5_unified = 35265730.0;
	static const double FullEXT_Run5_unified = 107466402.0;		

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

		{ "SingleBinPlot",  std::make_pair(0, 0.099) },
		//{ "SingleBinPlot",  std::make_pair(0, 0.159) },	
		{ "Pi0MomentumPlot",  std::make_pair(0, 0.34) },			

	};	
	
	//----------------------------------------//

	static std::map<TString,TString> VarLabel =
	{

		{ "SingleBinPlot",  "#sigma #left[10^{-38} #frac{cm^{2}}{Ar}#right]" },
		{ "Pi0MomentumPlot",  "#frac{d#sigma}{dp_{#pi^{ 0}}} #left[10^{-38} #frac{cm^{2}}{(GeV/c) Ar}#right]" },

	};

	static std::map<TString,TString> LatexLabel =
	{

		{ "MuonCosThetaSingleBinPlot",  "all events" },
		{ "Pi0MomentumPlot", "all events" },	
	
	};	
	
	//----------------------------------------//
	
	static std::map<TString,TString> MapUncorCor =
	{

		{ "SingleBinPlot", "SingleBinPlot" },
		{ "Pi0MomentumPlot", "Pi0MomentumPlot" },
	
	};
					
	//----------------------------------------//

	// Global Constants

	static const double Units = 1E38; // so that the extracted cross-section is in 10^{-38} cm^{2}
	static const double NTargets = 1.05E30; // Argon nuclei, not nucleons	

	static double POTUncertainty = 0.02; // 2% POT Uncertainty		

	static double NTargetUncertainty = 0.01; // 1% NTarget Uncertainty	
	
	static const int NuMuPdg = 14, MuonPdg = 13, ProtonPdg = 2212, AbsChargedPionPdg = 211, NeutralPionPdg = 111;
	static const int ElectronPdg = 11, nue_pdg = 12, PhotonPdg = 22, NeutronPdg = 2112;
	static const int KaonPdg = 321, NeutralKaonPdg = 311, NeutralKaonLongPdg = 130, NeutralKaonShortPdg = 310;
	static const int xi_pdg = 3312, xi0_pdg = 3322;
	static const int d0_pdg = 38, dp_pdg = 39, dm_pdg = 40;
	static const int SigmaPlusPdg = 3222, SigmaMinusPdg = 3112, NeutralSigmaPdg = 3212, LambdaPdg = 3122, omega_pdg = 223;
	static const int DeuteriumPdg = 1000010020, HeliumPdg = 1000020040, ArgonPdg = 1000180400;
	static const int rho_pdg = 113, charged_rho_pdg = 213, eta_pdg = 221;	
	static const int hydrogen_cluster_pdg = 2000000101, nucleon_pair = 2000000201;
	static const int proton_pair = 2000000202, neutron_pair = 2000000200;

	static const double MuonMass = 106, ProtonMass = 938.272, NeutronMass = 939.565; // MeV
	static const double MuonMass_GeV = 0.106, ProtonMass_GeV = 0.938272, NeutronMass_GeV = 0.939565, pi0_mass_gev = 0.135; // GeV
	static const double DeltaM2 = TMath::Power(NeutronMass_GeV,2.) - TMath::Power(ProtonMass_GeV,2.); // GeV^2	

	//----------------------------------------//

	// Plots to be included for xsec extraction purposes

	vector<TString> PlotNames{
				 "SingleBinPlot"
				 ,"Pi0MomentumPlot"
			 
	};
	
	//----------------------------------------//
	
	vector<TString> OneDimXSec = {
				 "SingleBinPlot"
				 ,"Pi0MomentumPlot"
				 				 
	};	
	
	//----------------------------------------//

	// Binning

	static const int NBinsSingleBin = 1; static const double ArrayNBinsSingleBin[NBinsSingleBin+1] = {0.,1.};
	
	static const int NBinsPi0CosTheta = 18;
	static const double ArrayNBinsPi0CosTheta[NBinsPi0CosTheta+1] = { -1.,-0.85,-0.7,-0.57,-0.45,-0.32,-0.2,-0.1,0.,0.1,0.2,0.3,0.4,0.5,0.6,0.72,0.84,0.95,1.};

	static const int NBinsPi0Momentum = 7; static const double ArrayNBinsPi0Momentum[NBinsPi0Momentum+1] = {0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8};

	static const int NBinsPi0InvMass = 15; static const double ArrayNBinsPi0InvMass[NBinsPi0InvMass+1] = {0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.,1.1,1.2,1.3,1.4,1.5};

	static const int NBinsg1CosTheta = 18;
	static const double ArrayNBinsg1CosTheta[NBinsg1CosTheta+1] = { -1.,-0.85,-0.7,-0.57,-0.45,-0.32,-0.2,-0.1,0.,0.1,0.2,0.3,0.4,0.5,0.6,0.72,0.84,0.95,1.};

	static const int NBinsg2CosTheta = 18;
	static const double ArrayNBinsg2CosTheta[NBinsg2CosTheta+1] = { -1.,-0.85,-0.7,-0.57,-0.45,-0.32,-0.2,-0.1,0.,0.1,0.2,0.3,0.4,0.5,0.6,0.72,0.84,0.95,1.};

	static const int NBinsg1Momentum = 16; static const double ArrayNBinsg1Momentum[NBinsg1Momentum+1] = {0.,0.1,0.2,0.3,0.38,0.45,0.5,0.55,0.625,0.7,0.75,0.8,0.87,1.,1.2,1.5,1.8};

	static const int NBinsg2Momentum = 13; static const double ArrayNBinsg2Momentum[NBinsg2Momentum+1] = {0.,0.1,0.2,0.3,0.38,0.45,0.5,0.55,0.625,0.7,0.75,0.8,0.87,1.};

	static const int NBinstwo_shower_angle = 18; static const double ArrayNBinstwo_shower_angle[NBinstwo_shower_angle+1] = {0.,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180};

	static const int NBinstwo_shower_start_dist = 18; static const double ArrayNBinstwo_shower_start_dist[NBinstwo_shower_start_dist+1] = {0.,20,40,60,80,100,120,140,160,180,200,220,240,260,280,300,320,340,360};
	
	static const int NBinsnshowers = 8; static const double ArrayNBinsnshowers[NBinsnshowers+1] = {-0.5,0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5};
	static const int NBinspd_nshowers = 8; static const double ArrayNBinspd_nshowers[NBinspd_nshowers+1] = {-0.5,0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5};		

	//----------------------------------------//

	// Labels for 1D plots
	
	static TString LabelXAxisSingleBin = ";"; static TString LabelXAxisTrueSingleBin = ";true";
	static TString LabelXAxisPi0CosTheta = ";cos#theta_{#pi^{ 0}}"; static TString LabelXAxisTruePi0CosTheta = ";true cos#theta_{#pi^{ 0}}";
	static TString LabelXAxisPi0Momentum = ";#pi^{ 0} momentum [GeV/c]"; static TString LabelXAxisTruePi0Momentum = ";true #pi^{0} momentum [GeV/c]";
	static TString LabelXAxisPi0InvMass= ";#pi^{ 0} invariant mass [GeV/c^{2}]"; static TString LabelXAxisTruePi0InvMass = ";true #pi^{0} invariant mass [GeV/c^{2}]";
	static TString LabelXAxisg1CosTheta = ";cos#theta_{#gamma_{1}}"; static TString LabelXAxisTrueg1CosTheta = ";true cos#theta_{#gamma^{1}}";
	static TString LabelXAxisg2CosTheta = ";cos#theta_{#gamma_{2}}"; static TString LabelXAxisTrueg2CosTheta = ";true cos#theta_{#gamma^{2}}";
	static TString LabelXAxisg1Momentum = ";#gamma_{1} momentum [GeV/c]"; static TString LabelXAxisTrueg1Momentum = ";true #gamma_{1} momentum [GeV/c]";
	static TString LabelXAxisg2Momentum = ";#gamma_{2} momentum [GeV/c]"; static TString LabelXAxisTrueg2Momentum = ";true #gamma_{2} momentum [GeV/c]";
	static TString LabelXAxistwo_shower_angle = ";two-shower angle [deg]"; static TString LabelXAxisTruetwo_shower_angle = ";true two-shower angle [deg]";
	static TString LabelXAxistwo_shower_start_dist = ";two-shower start distance [cm]"; static TString LabelXAxisTruetwo_shower_start_dist = ";true two-shower start distance [cm]";	
	static TString LabelXAxisnshowers = ";# wc showers"; static TString LabelXAxisTruenshowers = ";true # wc showers";
	static TString LabelXAxispd_nshowers = ";# pd showers"; static TString LabelXAxisTruepd_nshowers = ";true # pd showers";		

	//----------------------------------------//

	// Labels for 2D Plots

	static TString LabelXAxisSingleBin2D = LabelXAxisTrueSingleBin+";reco single bin";
	static TString LabelXAxisPi0CosTheta2D = LabelXAxisTruePi0CosTheta+";reco cos#theta_{#pi^{ 0}}";
	static TString LabelXAxisPi0Momentum2D = LabelXAxisTruePi0Momentum+";reco #pi^{ 0} momentum [GeV/c]";
	static TString LabelXAxisPi0InvMass2D = LabelXAxisTruePi0InvMass+";reco #pi^{ 0} invariant mass [GeV/c^{2}]";
	static TString LabelXAxisg1CosTheta2D = LabelXAxisTrueg1CosTheta+";reco cos#theta_{#gamma_{1}}";
	static TString LabelXAxisg2CosTheta2D = LabelXAxisTrueg2CosTheta+";reco cos#theta_{#gamma_{2}}";
	static TString LabelXAxisg1Momentum2D = LabelXAxisTrueg1Momentum+";reco #gamma_{1} momentum [GeV/c]";
	static TString LabelXAxisg2Momentum2D = LabelXAxisTrueg2Momentum+";reco #gamma_{2} momentum [GeV/c]";
	static TString LabelXAxistwo_shower_angle2D = LabelXAxisTruetwo_shower_angle+";reco two-shower angle [deg]";
	static TString LabelXAxistwo_shower_start_dist2D = LabelXAxisTruetwo_shower_start_dist+";reco two-shower start distance [cm]";	
	static TString LabelXAxisnshowers2D = LabelXAxisTruenshowers+";reco # wc showers";
	static TString LabelXAxispd_nshowers2D = LabelXAxisTruepd_nshowers+";reco # pd showers";	

	//----------------------------------------//
	
	// Interaction labels
	
	const std::vector<int> InteBreakColors{kBlack,kAzure-4,kOrange-3,kGreen+1,kRed+1,kBlue};		
	std::vector<TString> InteractionLabels = {"","QE","MEC","RES","DIS","COH","OTHER"};
	const int NInte = InteractionLabels.size();
	
	//----------------------------------------//
	
	// FV

	double FVx = 256., FVy = 232., FVz = 1037.;
	double borderx = 10., bordery = 10., borderz = 10.;

	//----------------------------------------//

	// neutron counters

	static TString LabelXAxisneutron_counter = ";# off-vertex objects"; 
	static TString LabelXAxisTrueneutron_counter = ";true off-vertex objects";	
	static const int NBinsneutron_counter = 6;
	static TString LabelXAxisneutron_counter2D = LabelXAxisTrueneutron_counter+";reco off-vertex objects";
	double min_neutron_counter = -0.5, max_neutron_counter = 5.5;

	//----------------------------------------//

	// gamma start / end points

	// start g1

	static TString LabelXAxisg1_start_x = ";#gamma_{1} start x [cm]"; 
	static TString LabelXAxisTrueg1_start_x = ";True #gamma_{1} x [cm]";	
	static const int NBinsg1_start_x = 20;
	static TString LabelXAxisg1_start_x2D = LabelXAxisTrueg1_start_x+";Reco #gamma_{1} x [cm]";
	double min_g1_start_x = 0., max_g1_start_x = FVx;

	static TString LabelXAxisg1_start_y = ";#gamma_{1} start y [cm]"; 
	static TString LabelXAxisTrueg1_start_y = ";True #gamma_{1} start y [cm]";	
	static const int NBinsg1_start_y = 20;
	static TString LabelXAxisg1_start_y2D = LabelXAxisTrueg1_start_y+";Reco #gamma_{1} start [cm]";
	double min_g1_start_y = -FVy/2., max_g1_start_y = FVy/2.;
	
	static TString LabelXAxisg1_start_z = ";#gamma_{1} start z [cm]"; 
	static TString LabelXAxisTrueg1_start_z = ";True #gamma_{1} z start [cm]";	
	static const int NBinsg1_start_z = 20;
	static TString LabelXAxisg1_start_z2D = LabelXAxisTrueg1_start_z+";Reco #gamma_{1} start z [cm]";
	double min_g1_start_z = 0., max_g1_start_z = FVz;	

	// end g1

	static TString LabelXAxisg1_end_x = ";#gamma_{1} end x [cm]"; 
	static TString LabelXAxisTrueg1_end_x = ";True #gamma_{1} end x [cm]";	
	static const int NBinsg1_end_x = 20;
	static TString LabelXAxisg1_end_x2D = LabelXAxisTrueg1_end_x+";Reco #gamma_{1} end x [cm]";
	double min_g1_end_x = 0., max_g1_end_x = FVx;

	static TString LabelXAxisg1_end_y = ";#gamma_{1} end y [cm]"; 
	static TString LabelXAxisTrueg1_end_y = ";True #gamma_{1} end y [cm]";	
	static const int NBinsg1_end_y = 20;
	static TString LabelXAxisg1_end_y2D = LabelXAxisTrueg1_end_y+";Reco #gamma_{1} end y [cm]";
	double min_g1_end_y = -FVy/2., max_g1_end_y = FVy/2.;
	
	static TString LabelXAxisg1_end_z = ";#gamma_{1} end z [cm]"; 
	static TString LabelXAxisTrueg1_end_z = ";True #gamma_{1} end z [cm]";	
	static const int NBinsg1_end_z = 20;
	static TString LabelXAxisg1_end_z2D = LabelXAxisTrueg1_end_z+";Reco #gamma_{1} end z [cm]";
	double min_g1_end_z = 0., max_g1_end_z = FVz;		

	// start g2

	static TString LabelXAxisg2_start_x = ";#gamma_{2} start x [cm]"; 
	static TString LabelXAxisTrueg2_start_x = ";True #gamma_{2} start x [cm]";	
	static const int NBinsg2_start_x = 20;
	static TString LabelXAxisg2_start_x2D = LabelXAxisTrueg2_start_x+";Reco #gamma_{2} start x [cm]";
	double min_g2_start_x = 0., max_g2_start_x = FVx;

	static TString LabelXAxisg2_start_y = ";#gamma_{2} start y [cm]"; 
	static TString LabelXAxisTrueg2_start_y = ";True #gamma_{2} start y [cm]";	
	static const int NBinsg2_start_y = 20;
	static TString LabelXAxisg2_start_y2D = LabelXAxisTrueg2_start_y+";Reco #gamma_{2} start y [cm]";
	double min_g2_start_y = -FVy/2., max_g2_start_y = FVy/2.;
	
	static TString LabelXAxisg2_start_z = ";#gamma_{2} start z [cm]"; 
	static TString LabelXAxisTrueg2_start_z = ";True #gamma_{2} start z [cm]";	
	static const int NBinsg2_start_z = 20;
	static TString LabelXAxisg2_start_z2D = LabelXAxisTrueg2_start_z+";Reco #gamma_{2} start z [cm]";
	double min_g2_start_z = 0., max_g2_start_z = FVz;	
	
	// end g2

	static TString LabelXAxisg2_end_x = ";#gamma_{2} end x [cm]"; 
	static TString LabelXAxisTrueg2_end_x = ";True #gamma_{2} end x [cm]";	
	static const int NBinsg2_end_x = 20;
	static TString LabelXAxisg2_end_x2D = LabelXAxisTrueg2_end_x+";Reco #gamma_{2} end x [cm]";
	double min_g2_end_x = 0., max_g2_end_x = FVx;

	static TString LabelXAxisg2_end_y = ";#gamma_{2} end y [cm]"; 
	static TString LabelXAxisTrueg2_end_y = ";True #gamma_{2} end y [cm]";	
	static const int NBinsg2_end_y = 20;
	static TString LabelXAxisg2_end_y2D = LabelXAxisTrueg2_end_y+";Reco #gamma_{2} end y [cm]";
	double min_g2_end_y = -FVy/2., max_g2_end_y = FVy/2.;
	
	static TString LabelXAxisg2_end_z = ";#gamma_{2} end z [cm]"; 
	static TString LabelXAxisTrueg2_end_z = ";True #gamma_{2} end z [cm]";	
	static const int NBinsg2_end_z = 20;
	static TString LabelXAxisg2_end_z2D = LabelXAxisTrueg2_end_z+";Reco #gamma_{2} end z [cm]";
	double min_g2_end_z = 0., max_g2_end_z = FVz;		

	//----------------------------------------//		

	static TString LabelXAxisvertex_x = ";vertex x [cm]"; 
	static TString LabelXAxisTruevertex_x = ";True vertex x [cm]";	
	static const int NBinsvertex_x = 50;
	static TString LabelXAxisvertex_x2D = LabelXAxisTruevertex_x+";Reco vertex x [cm]";
	double min_vertex_x = 0., max_vertex_x = FVx;	

	static TString LabelXAxisvertex_y = ";vertex y [cm]"; 
	static TString LabelXAxisTruevertex_y = ";True vertex y [cm]";	
	static const int NBinsvertex_y = 50;
	static TString LabelXAxisvertex_y2D = LabelXAxisTruevertex_y+";Reco vertex y [cm]";
	double min_vertex_y = -FVy/2., max_vertex_y = FVy/2.;

	static TString LabelXAxisvertex_z = ";vertex z [cm]"; 
	static TString LabelXAxisTruevertex_z = ";True vertex z [cm]";	
	static const int NBinsvertex_z = 50;
	static TString LabelXAxisvertex_z2D = LabelXAxisTruevertex_z+";Reco vertex z [cm]";
	double min_vertex_z = 0., max_vertex_z = FVz;	

	static TString LabelXAxiskine_pio_vtx_dis = ";pio-vertex distance [cm]"; 
	static TString LabelXAxisTruekine_pio_vtx_dis = ";True pio-vertex distance [cm]";	
	static const int NBinskine_pio_vtx_dis = 20;
	static TString LabelXAxiskine_pio_vtx_dis2D = LabelXAxisTruekine_pio_vtx_dis+";Reco pio-vertex distance [cm]";
	double min_kine_pio_vtx_dis = 0, max_kine_pio_vtx_dis = 20;

	static TString LabelXAxissingle_photon_numu_score = ";single-photon numu score"; 
	static TString LabelXAxisTruesingle_photon_numu_score = ";True single-photon numu score";	
	static const int NBinssingle_photon_numu_score = 20;
	static TString LabelXAxissingle_photon_numu_score2D = LabelXAxisTruesingle_photon_numu_score+";Reco single-photon numu score";
	double min_single_photon_numu_score = -4, max_single_photon_numu_score = 4;

	static TString LabelXAxissingle_photon_other_score = ";single-photon other score"; 
	static TString LabelXAxisTruesingle_photon_other_score = ";True single-photon other score";	
	static const int NBinssingle_photon_other_score = 20;
	static TString LabelXAxissingle_photon_other_score2D = LabelXAxisTruesingle_photon_other_score+";Reco single-photon other score";
	double min_single_photon_other_score = -5, max_single_photon_other_score = 4;

	static TString LabelXAxissingle_photon_ncpi0_score = ";single-photon pio score"; 
	static TString LabelXAxisTruesingle_photon_ncpi0_score = ";True single-photon pio score";	
	static const int NBinssingle_photon_ncpi0_score = 20;
	static TString LabelXAxissingle_photon_ncpi0_score2D = LabelXAxisTruesingle_photon_ncpi0_score+";Reco single-photon pio score";
	double min_single_photon_ncpi0_score = -3, max_single_photon_ncpi0_score = 3;

	static TString LabelXAxissingle_photon_nue_score = ";single-photon nue score"; 
	static TString LabelXAxisTruesingle_photon_nue_score = ";True single-photon nue score";	
	static const int NBinssingle_photon_nue_score = 20;
	static TString LabelXAxissingle_photon_nue_score2D = LabelXAxisTruesingle_photon_nue_score+";Reco single-photon nue score";
	double min_single_photon_nue_score = -5, max_single_photon_nue_score = 2;

	static TString LabelXAxisnc_pio_score = ";nc pio score"; 
	static TString LabelXAxisTruenc_pio_score = ";True nc pio score";	
	static const int NBinsnc_pio_score = 20;
	static TString LabelXAxisnc_pio_score2D = LabelXAxisTruenc_pio_score+";Reco nc pio score";
	double min_nc_pio_score = -10, max_nc_pio_score = 6;

	static TString LabelXAxisnumu_score = ";numu score"; 
	static TString LabelXAxisTruenumu_score = ";True numu score";	
	static const int NBinsnumu_score = 50;
	static TString LabelXAxisnumu_score2D = LabelXAxisTruenumu_score+";Reco numu score";
	double min_numu_score = -3., max_numu_score = 3.;

	static TString LabelXAxiskine_pio_flag = ";kine pio flag"; 
	static TString LabelXAxisTruekine_pio_flag = ";True kine pio flag";	
	static const int NBinskine_pio_flag = 3;
	static TString LabelXAxiskine_pio_flag2D = LabelXAxisTruekine_pio_flag+";Reco kine pio flag";
	double min_kine_pio_flag = -.5, max_kine_pio_flag = 2.5;	

	// Blips

	static TString LabelXAxisBlip_x = ";x-blip [cm]"; 
	static TString LabelXAxisTrueBlip_x = ";True x-blip [cm]";	
	static const int NBinsBlip_x = 20;
	static TString LabelXAxisBlip_x2D = LabelXAxisTrueBlip_x+";Reco x-blip [cm]";

	static TString LabelXAxisBlip_y = ";y-blip [cm]"; 
	static TString LabelXAxisTrueBlip_y = ";True y-blip [cm]";	
	static const int NBinsBlip_y = 20;
	static TString LabelXAxisBlip_y2D = LabelXAxisTrueBlip_y+";Reco y-blip [cm]";	

	static TString LabelXAxisBlip_z = ";z-blip [cm]"; 
	static TString LabelXAxisTrueBlip_z = ";True z-blip [cm]";	
	static const int NBinsBlip_z = 20;
	static TString LabelXAxisBlip_z2D = LabelXAxisTrueBlip_z+";Reco z-blip [cm]";		

	TString LabelXAxisnBlips_saved = ";# saved blips", LabelXAxisnBlips_saved2D = ";true # saved blips"; 
	int NBinsnBlips_saved = 16;
	double nBlips_saved_min = -0.5, nBlips_saved_max = 15.5;

	int radius = 70; // cm, since neutron interaction length is 70 cm
	TString LabelXAxisnBlips_radius = ";# blips within " + TString( to_string(radius) ) + "cm from vertex", LabelXAxisnBlips_radius2D = ";true # radius blips within " + TString( to_string(radius) ) + "cm from vertex"; 
	int NBinsnBlips_radius = 8;
	double nBlips_radius_min = -0.5, nBlips_radius_max = 7.5;
	
	TString LabelXAxisBlip_size = ";blip size [cm]", LabelXAxisBlip_size2D = ";true blip size [cm]"; 
	int NBinsBlip_size = 15;
	double Blip_size_min = 0., Blip_size_max = 3.;	

	TString LabelXAxisBlip_energy = ";blip energy [MeV]", LabelXAxisBlip_energy2D = ";true blip energy [MeV]"; 
	int NBinsBlip_energy = 20;
	double Blip_energy_min = 0., Blip_energy_max = 5.;	

	TString LabelXAxisBlip_proxtrkdist = ";blip-proxtrack distance [cm]", LabelXAxisBlip_proxtrkdist2D = ";true blip-proxtrack distance [cm]"; 
	int NBinsBlip_proxtrkdist = 50;
	double Blip_proxtrkdist_min = 0., Blip_proxtrkdist_max = 250.;
	
	TString LabelXAxisBlip_pairdist = ";blip-pair track min distance [cm]", LabelXAxisBlip_pairdist2D = ";true blip-pair track min distance [cm]"; 
	int NBinsBlip_pairdist = 50;
	double Blip_pairdist_min = 0., Blip_pairdist_max = 250.;	

	TString LabelXAxisblip_vrt = ";blip-vertex distance [cm]", LabelXAxisblip_vrt2D = ";true blip-vertex distance [cm]"; 
	int NBinsblip_vrt = 50;
	double blip_vrt_min = 0., blip_vrt_max = 250.;		

	TString LabelXAxisblip_cos_alphapi0 = ";blip cos#alpha_{#pi^{0}} within " + TString( to_string(radius) ) + "cm from vertex", LabelXAxisblip_cos_alphapi02D = ";true blip cos#alpha_{#pi^{0}} within " + TString( to_string(radius) ) + "cm from vertex"; 
	int NBinsblip_cos_alphapi0 = 50;
	double blip_cos_alphapi0_min = -1., blip_cos_alphapi0_max = 1.;

	TString LabelXAxisblip_cos_alphag1 = ";blip cos#alpha_{#gamma_{1}} within " + TString( to_string(radius) ) + "cm from vertex", LabelXAxisblip_cos_alphag12D = ";true blip cos#alpha_{#gamma_{1}} within " + TString( to_string(radius) ) + "cm from vertex"; 
	int NBinsblip_cos_alphag1 = 50;
	double blip_cos_alphag1_min = -1., blip_cos_alphag1_max = 1.;
	
	TString LabelXAxisblip_cos_alphag2 = ";blip cos#alpha_{#gamma_{2}} within " + TString( to_string(radius) ) + "cm from vertex", LabelXAxisblip_cos_alphag22D = ";true blip cos#alpha_{#gamma_{2}} within " + TString( to_string(radius) ) + "cm from vertex"; 
	int NBinsblip_cos_alphag2 = 50;
	double blip_cos_alphag2_min = -1., blip_cos_alphag2_max = 1.;

	//----------------------------------------//	


}

#endif