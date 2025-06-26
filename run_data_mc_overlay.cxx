{

	gROOT->ProcessLine(".L Util.C");
	gROOT->ProcessLine(".L data_mc_overlay.cxx++");
	
	//GENIE versions
	gROOT->ProcessLine("data_mc_overlay()");
	//AltGen
	gROOT->ProcessLine("data_mc_overlay(false,true)");
	//Closure test
	gROOT->ProcessLine("data_mc_overlay(false,false,true)");

}
