{

	gROOT->ProcessLine(".L generator_overlay.cxx+");
	gROOT->ProcessLine("generator_overlay()");
	//gROOT->ProcessLine("generator_overlay(\"NuclearModel\")");

};
