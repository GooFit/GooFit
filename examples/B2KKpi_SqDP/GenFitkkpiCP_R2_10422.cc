#include "TSystem.h"

#include <cstdlib>
#include <iostream>
#include <vector>

#include "TFile.h"
#include "TH2.h"
#include "TString.h"
#include "TTree.h"

#include "LauCPFitModel.hh"
#include "LauBkgndDPModel.hh"
#include "LauDaughters.hh"
#include "LauEffModel.hh"
#include "LauIsobarDynamics.hh"
#include "LauCartesianCPCoeffSet.hh"
#include "LauVetoes.hh"
#include "LauResonanceMaker.hh"
#include "LauRescatteringRes.hh"
#include "LauMagPhaseCPCoeffSet.hh"
#include "LauAbsCoeffSet.hh"

void usage( std::ostream& out, const TString& progName )
{
	out<<"Usage:\n";
	out<<progName<<" gen [nExpt = 1] [firstExpt = 0]\n";
	out<<"or\n";
	out<<progName<<" fit <iFit> [nExpt = 1] [firstExpt = 0]"<<std::endl;
}

int main( int argc, char** argv )
{
	// Process command-line arguments
	// Usage:
	// ./GenFit3pi gen [nExpt = 1] [firstExpt = 0]
	// or
	// ./GenFit3pi fit <iFit> [nExpt = 1] [firstExpt = 0]
	if ( argc < 2 ) {
		usage( std::cerr, argv[0] );
		return EXIT_FAILURE;
	}

	TString command = argv[1];
	command.ToLower();
	Int_t iFit(0);
	Int_t nExpt(1);
	Int_t firstExpt(0);
	if ( command == "gen" ) {
		if ( argc > 2 ) {
			nExpt = atoi( argv[2] );
			if ( argc > 3 ) {
				firstExpt = atoi( argv[3] );
			}
		}
	} else if ( command == "fit" ) {
		if ( argc < 3 ) {
			usage( std::cerr, argv[0] );
			return EXIT_FAILURE;
		}
		iFit = atoi( argv[2] );
		if ( argc > 3 ) {
			nExpt = atoi( argv[3] );
			if ( argc > 4 ) {
				firstExpt = atoi( argv[4] );
			}
		}
	} else {
		usage( std::cerr, argv[0] );
		return EXIT_FAILURE;
	}

    /************************************************************************/
    /*** Inputs globais (inicio)                                         ****/

        TString AccFileBpn   = "../../inputs/ACC_kkpi.root";
		//RUN1 -> TString BkgFile      = "../../inputs/CombBkg_v5.root";
        TString BkgFile_plus      = "../../inputs/histos_kkpi_SqDP_Adaptive_2x2bins_iter3_GaussFree_Plus.root";
        TString BkgFile_minus      = "../../inputs/histos_kkpi_SqDP_Adaptive_2x2bins_iter3_GaussFree_Minus.root";		
        TString BkgKpipiFile = "../../inputs/PeakingBkg_v5.root";
		TString BkgPhiFile = "../../inputs/Histo_phiBkg_withNewACC.root";
        TString DataFile     = "../../inputs/dataLHCbForLaura-KKPi-TOS_or_TIS-New_veto_RUN2.root";
        //
        // Create the folder where the results will be saved
        TDatime a;
	TString resDir = "tmp/";
	resDir+=command;
	if (command == "fit") {resDir += iFit;}	
        TString cont  = Form( "_kkpinotav3_%d-%d--%dh%dm%ds", a.GetDay(), a.GetMonth(), a.GetHour(), a.GetMinute(), a.GetSecond() );
        resDir += cont;
	if (command == "fit") {resDir += iFit;}	
        resDir += "/";
        gSystem->Exec("mkdir "+resDir);

    /*** Fim do "Inputs globais"                                         ****/
    /************************************************************************/
 

	// If you want to use square DP histograms for efficiency,
	// backgrounds or you just want the square DP co-ordinates
	// stored in the toy MC ntuple then set this to kTRUE
	Bool_t squareDP = kTRUE;
	
	// Set this to kFALSE if you want to remove the DP from the fit
	Bool_t doDP = kTRUE;

	// Set this to kTRUE if you want to fix the CPV parameters
	//Bool_t fixCP = kTRUE;
	Bool_t doTwoStageFit = kFALSE;
		

	// This defines the DP => decay is B+ -> K+ pi+ K-
	// Particle 1 = K+
	// Particle 2 = pi+ //hhh group convention, For Laura++ we will use the described below.
	// Particle 3 = K+
	// The DP is defined in terms of m12Sq and m12Sq**
	LauDaughters* posdaughters = new LauDaughters("B+", "K+", "pi+", "K-", squareDP); 
	LauDaughters* negdaughters = new LauDaughters("B-", "K-", "pi-", "K+", squareDP);



	// Signal and continuum yields - Obtidos do Fit de massa 1D B+ + B-  na região de sinal de +- 17 MeV
	Double_t nSigEvents = 36280 *0.76;
//	Double_t nSigSigma = 100; //erro do numero de eventos.
	Bool_t fixNSigEvents = kTRUE; //Fixed or not
	Double_t nBgEvents =   36280 *0.24;
	Bool_t fixNBgEvents = kTRUE; //Fixed or not
	Double_t nBgKpipiEvents    = 0; //run1 -> 130
	Bool_t   fixNBg3KpipiEvents = kTRUE;//Fixed or not.
	Double_t nBgPhiEvents    = 0; // run1 -> 30
	Bool_t   fixNBgPhiEvents = kTRUE;//Fixed or not.

	// Signal and continuum asymmetries
	// NB the signal asymmetry value here is only used if the DP is NOT
	// in the fit, otherwise the asymmetry is included within the
	// isobar model.
	
	Double_t sigAsym = -0.123;
	Bool_t fixSigAsym =  kFALSE; // -->Fixar ou nao a asimetria.
	Double_t bgAsym =  -0.009401;
	Bool_t fixBgAsym = kTRUE;
	Double_t bgAsymKpipi = 0.025;
	Bool_t fixBgAsymKpipi = kTRUE;
	Double_t bgAsymPhi = 0.00;
	Bool_t fixBgAsymPhi = kTRUE;



	//Bachelor
	double kpi_bachelor = 1;
	double kk_bachelor = 2;

	// Optionally apply some vetoes to the DP
	LauVetoes* vetoes = new LauVetoes();
	Double_t DMin23 = 1.83;//1.834; (1.834 era o antigo, 1.83 é o novo sugerido pela Melissa)
	Double_t DMax23 = 2.00;//1.894;

	Double_t DMin31 = 1.83;//1.834;
	Double_t DMax31 = 2.00;//2.000;

	vetoes->addMassVeto(kk_bachelor, DMin31, DMax31); // D0 veto, mKK
	vetoes->addMassVeto(kpi_bachelor, DMin23, DMax23); // D0 veto, mpiK


	
	Bool_t useInterpolation = kTRUE;
	Bool_t fluctuateBins = kFALSE;
	Bool_t useUpperHalfOnly = kFALSE;




	printf("*********************************ACEPTANCIA-Se comentado a parte onde pega acc, entao eff = 1 em todo DP************************************************************* \n");



	// Define the efficiency model (defaults to unity everywhere)------> If I don't give an Acceptance histogram, the program use eff=1 for everyone.
	// Can optionally provide a histogram to model variation over DP
	// (example syntax given in commented-out section)
	// General histogram booleans
	
	LauEffModel* effModelpos = new LauEffModel(posdaughters, vetoes);
	LauEffModel* effModelneg = new LauEffModel(negdaughters, vetoes);


	TString effFileNameBpn(AccFileBpn);
	TFile *effHistFileBpn = TFile::Open(effFileNameBpn, "read");
	std::cout<<" using EFFILE:" << effFileNameBpn << std::endl;
	TH2* effHistBp = dynamic_cast<TH2*>(effHistFileBpn->Get("Hacc_trg_kkpi_151617_Plus_NoDAsyCorr"));
	TH2* effHistBn = dynamic_cast<TH2*>(effHistFileBpn->Get("Hacc_trg_kkpi_151617_Minus_NoDAsyCorr"));
	
	effModelneg->setEffHisto(effHistBn, useInterpolation, fluctuateBins, 0.0, 0.0, useUpperHalfOnly, squareDP);
	effModelpos->setEffHisto(effHistBp, useInterpolation, fluctuateBins, 0.0, 0.0, useUpperHalfOnly, squareDP);


	//***********************************************************************************************************//


	// Create the isobar model

	// Set the values of the Blatt-Weisskopf barrier radii and whether they are fixed or floating
	LauResonanceMaker& resMaker = LauResonanceMaker::get();
	resMaker.setDefaultBWRadius( LauBlattWeisskopfFactor::Parent,     4.0 ); //Meson radius parameter used in the Blatt-weisskopf barrier factors
	resMaker.setDefaultBWRadius( LauBlattWeisskopfFactor::Light,      4.0 ); // values can be found on PDG.
	resMaker.fixBWRadius( LauBlattWeisskopfFactor::Parent,  kTRUE );
	resMaker.fixBWRadius( LauBlattWeisskopfFactor::Light,  kTRUE );


	printf("*********************************Signal Model Positive ************************************************************* \n");

	LauIsobarDynamics* sigModelpos = new LauIsobarDynamics(posdaughters, effModelpos);

	sigModelpos->setASqMaxValue(0.5);  


	LauAbsResonance* reson(0);
 	reson = sigModelpos->addResonance("K*0(892)",    kpi_bachelor, LauAbsResonance::RelBW); // bachelor da melissa era 3 o meu é o 1 (K+)
	reson = sigModelpos->addResonance("K*0_0(1430)", kpi_bachelor, LauAbsResonance::RelBW);
	reson = sigModelpos->addResonance("PolarFFNR",     kpi_bachelor, LauAbsResonance::PolarFFNR);
	reson->setResonanceParameter("lambda", 1.18);
 	reson = sigModelpos->addResonance("rho0(1700)", kk_bachelor, LauAbsResonance::RelBW);
 	reson = sigModelpos->addResonance("rho0(1450)", kk_bachelor, LauAbsResonance::RelBW);
	reson = sigModelpos->addResonance("f'_2(1525)", kk_bachelor, LauAbsResonance::RelBW);
	reson = sigModelpos->addResonance("f_2(1270)", kk_bachelor, LauAbsResonance::RelBW);
 	reson = sigModelpos->addResonance("Rescattering",    kk_bachelor, LauAbsResonance::Rescattering2);//Rescattering mais moderno, pelaez 2018
	reson = sigModelpos->addResonance("phi(1020)",   kk_bachelor, LauAbsResonance::RelBW);
	reson = sigModelpos->addResonance("chi_c0",   kk_bachelor, LauAbsResonance::RelBW);



	printf("*********************************Signal Model Negative modelo negativo************************************************************* \n");

	LauIsobarDynamics* sigModelneg = new LauIsobarDynamics(negdaughters, effModelneg);
	sigModelneg->setASqMaxValue(0.5);  


	reson = sigModelneg->addResonance("K*0(892)",    kpi_bachelor, LauAbsResonance::RelBW); // bachelor da melissa era 3 o meu é o 1 (K+)
	reson = sigModelneg->addResonance("K*0_0(1430)", kpi_bachelor, LauAbsResonance::RelBW); //Same settings for positive model.
	reson = sigModelneg->addResonance("PolarFFNR",     kpi_bachelor, LauAbsResonance::PolarFFNR);
	reson = sigModelneg->addResonance("rho0(1700)", kk_bachelor, LauAbsResonance::RelBW);
	reson = sigModelneg->addResonance("rho0(1450)", kk_bachelor, LauAbsResonance::RelBW);
	reson = sigModelneg->addResonance("f'_2(1525)", kk_bachelor, LauAbsResonance::RelBW);
	reson = sigModelneg->addResonance("f_2(1270)", kk_bachelor, LauAbsResonance::RelBW);
    reson = sigModelneg->addResonance("Rescattering",    kk_bachelor, LauAbsResonance::Rescattering2);
 	reson = sigModelneg->addResonance("phi(1020)",   kk_bachelor, LauAbsResonance::RelBW);
 	reson = sigModelneg->addResonance("chi_c0",   kk_bachelor, LauAbsResonance::RelBW);


 
	// Create the fit model
	
	 LauCPFitModel* fitModel = new LauCPFitModel(sigModelneg,sigModelpos);
	 fitModel->useDP(doDP);
	 


	// Create the complex coefficients for the isobar model
	std::vector<LauAbsCoeffSet*> coeffset;
	printf("*********************************Magnitude and Phase definition for the fitl Model ************************************************************* \n");

        coeffset.push_back(new LauCartesianCPCoeffSet("K*0(892)",     1.000,  0.00, -2.67326e-01, 0.000,  kTRUE,  kTRUE, kFALSE,  kTRUE, doTwoStageFit, doTwoStageFit));
        coeffset.push_back(new LauCartesianCPCoeffSet("K*0_0(1430)",  -1.82039e+00,  -6.31532e-01, -1.23741e-01,-4.59960e-01,  kFALSE,  kFALSE, kFALSE,  kFALSE, doTwoStageFit, doTwoStageFit));
        coeffset.push_back(new LauCartesianCPCoeffSet("PolarFFNR",     2.59302e+00,  -1.42689e+00, -1.15660e-01,-4.43618e-02,  kFALSE,  kFALSE, kFALSE,  kFALSE, doTwoStageFit, doTwoStageFit));
        coeffset.push_back(new LauCartesianCPCoeffSet("rho0(1700)", -1.92876e+00,  -2.06752e+00, 1.94301e-02, -3.63772e-01,  kFALSE,  kFALSE, kFALSE,  kFALSE, doTwoStageFit, doTwoStageFit));
        coeffset.push_back(new LauCartesianCPCoeffSet("rho0(1450)", -1.84713e+00,  -1.21421e-01, -4.69240e-01, 4.99873e-01,  kFALSE,  kFALSE, kFALSE,  kFALSE, doTwoStageFit, doTwoStageFit));
        coeffset.push_back(new LauCartesianCPCoeffSet("f'_2(1525)", -2.98060e-01,  -2.48489e+00, 6.95060e-02, -1.06819e-01, kFALSE,  kFALSE, kFALSE,  kFALSE, doTwoStageFit, doTwoStageFit));
        coeffset.push_back(new LauCartesianCPCoeffSet("f_2(1270)", -2.12933e+00,  -2.62605e-01, -4.81535e-01, 1.19802e-01, kFALSE,  kFALSE, kFALSE,  kFALSE, doTwoStageFit, doTwoStageFit));
        coeffset.push_back(new LauCartesianCPCoeffSet("Rescattering",     -2.59438e+00,  2.20930e+00, 1.34448e-01, -3.65647e-01,  kFALSE,  kFALSE, kFALSE,  kFALSE, doTwoStageFit, doTwoStageFit));
        coeffset.push_back(new LauCartesianCPCoeffSet("phi(1020)", 2.12347e+00,  9.53155e-01, -9.08623e-02, 3.84786e-01, kFALSE,  kFALSE, kFALSE,  kFALSE, doTwoStageFit, doTwoStageFit));////
        coeffset.push_back(new LauCartesianCPCoeffSet("chi_c0",-2.60479e+00,  -1.63962e+00, -4.09123e-01, -9.91773e-02, kFALSE,  kFALSE, kFALSE,  kFALSE, doTwoStageFit, doTwoStageFit));////


	for (std::vector<LauAbsCoeffSet*>::iterator iter=coeffset.begin(); iter!=coeffset.end(); ++iter) {
		fitModel->setAmpCoeffSet(*iter);
	}


	// Set the number of signal events and the number of experiments
	LauParameter * signalEvents = new    
	LauParameter("signalEvents",nSigEvents,-2.0*nSigEvents,2.0*nSigEvents,fixNSigEvents);
	LauParameter * signalAsym = new LauParameter("signalAsym",sigAsym,-1.0,1.0,fixSigAsym);
	fitModel->setNSigEvents(signalEvents, signalAsym);
	fitModel->setNExpts(nExpt, firstExpt);


	// Optionally load in continuum background DP model histogram
	// (example syntax given in commented-out section)
	std::vector<TString> bkgndNames(3); //changing 1 -> 2
	bkgndNames[0] = "qqbar";
	bkgndNames[1] = "peakingkpipi";
	bkgndNames[2] = "phiBkg";

	fitModel->setBkgndClassNames( bkgndNames );



	//bkg histos
/*	
	TString qqFileName(BkgFile);
	std::cout <<" using BKGFILE :" << qqFileName << std::endl;
	TFile* qqFile = TFile::Open(qqFileName, "read");
	TH2* negqqDP = dynamic_cast<TH2*>(qqFile->Get("h2_BkgDistrSplot_GaussFree_NoNegCut_Norm")); 
	TH2* posqqDP = dynamic_cast<TH2*>(qqFile->Get("h2_BkgDistrSplot_GaussFree_NoNegCut_Norm")); 
	LauBkgndDPModel* qqbarModelneg = new LauBkgndDPModel(negdaughters, vetoes);
	LauBkgndDPModel* qqbarModelpos = new LauBkgndDPModel(posdaughters, vetoes);
	qqbarModelneg->setBkgndHisto( negqqDP, useInterpolation, fluctuateBins, useUpperHalfOnly, squareDP );
	qqbarModelpos->setBkgndHisto( posqqDP, useInterpolation, fluctuateBins, useUpperHalfOnly, squareDP );

*/

	TString qqFileName_plus(BkgFile_plus);
	std::cout <<" using BKGFILE plus :" << qqFileName_plus << std::endl;
	TFile* qqFile_plus = TFile::Open(qqFileName_plus, "read");
	TH2* posqqDP = dynamic_cast<TH2*>(qqFile_plus->Get("h2_BkgDistrSplot_GaussFree_NoNegCut_Norm")); 
	LauBkgndDPModel* qqbarModelpos = new LauBkgndDPModel(posdaughters, vetoes);
	qqbarModelpos->setBkgndHisto( posqqDP, useInterpolation, fluctuateBins, useUpperHalfOnly, squareDP );



	TString qqFileName_minus(BkgFile_minus);
	std::cout <<" using BKGFILE minus :" << qqFileName_minus << std::endl;
	TFile* qqFile_minus = TFile::Open(qqFileName_minus, "read");
	TH2* negqqDP = dynamic_cast<TH2*>(qqFile_minus->Get("h2_BkgDistrSplot_GaussFree_NoNegCut_Norm")); 
	LauBkgndDPModel* qqbarModelneg = new LauBkgndDPModel(negdaughters, vetoes);
	qqbarModelneg->setBkgndHisto( negqqDP, useInterpolation, fluctuateBins, useUpperHalfOnly, squareDP );





	
	printf("Modelo de Background Combinatorial aceito \n");


//***************************************************************Load Kpipi DP histograms********************************//


	TString kpipiFileName(BkgKpipiFile);
	std::cout <<" using BKGFILE peaking kpipi :" << kpipiFileName << std::endl;
	TFile* kpipiFile = TFile::Open(kpipiFileName, "read");
	TH2* kpipiDPBpn = dynamic_cast<TH2*>(kpipiFile->Get("h1")); // m', theta'
	LauBkgndDPModel* KpipiModelpos = new LauBkgndDPModel(posdaughters, vetoes);
	LauBkgndDPModel* KpipiModelneg = new LauBkgndDPModel(negdaughters, vetoes);
	KpipiModelneg->setBkgndHisto(kpipiDPBpn, useInterpolation, fluctuateBins, useUpperHalfOnly, squareDP);
	KpipiModelpos->setBkgndHisto(kpipiDPBpn, useInterpolation, fluctuateBins, useUpperHalfOnly, squareDP);

	
	printf("Modelo de Peaking Background aceito \n");



	//************************************************************************************************************************************//	
	

	TString phiFileName(BkgPhiFile);
	std::cout <<" using BKGFILE phi :" << phiFileName << std::endl;
	TFile* phiFile = TFile::Open(phiFileName, "read");
	TH2* phiDPBpn = dynamic_cast<TH2*>(phiFile->Get("h1")); 
	LauBkgndDPModel* PhiBkgModelpos = new LauBkgndDPModel(posdaughters, vetoes);
	LauBkgndDPModel* PhiBkgModelneg = new LauBkgndDPModel(negdaughters, vetoes);
	PhiBkgModelneg->setBkgndHisto(phiDPBpn, useInterpolation, fluctuateBins, useUpperHalfOnly, squareDP);
	PhiBkgModelpos->setBkgndHisto(phiDPBpn, useInterpolation, fluctuateBins, useUpperHalfOnly, squareDP);


	printf("nBgEvents %f \n",nBgEvents );

  //************************************************************************************************************************************//	


	printf("Modelo de PhiBkg Background aceito \n");

	LauParameter* nBkgndEvents = new LauParameter("qqbar",nBgEvents,-2.0*nBgEvents,2.0*nBgEvents,fixNBgEvents);// Modificar aqui 2*sigma, 974.0,1198.0 em este momento está fixo o numero de eventos
	LauParameter* nKpipiBkgndEvents = new LauParameter("peakingkpipi",nBgKpipiEvents, 0.0, 260.0,fixNBg3KpipiEvents);//  -2*nBgKpipiEvents, 2*nBgKpipiEvents
	LauParameter* nPhiBkgndEvents  = new LauParameter("phiBkg",nBgPhiEvents,0.0,3*nBgPhiEvents,fixNBgPhiEvents);

	LauParameter* backgroundAsym = new LauParameter("backgroundAsym",bgAsym,-1.0, 1.0,fixBgAsym); //-0.036727,0.000383
	LauParameter* backgroundKpipiAsym = new LauParameter("backgroundAsymKpipi",bgAsymKpipi,-1.,1.,fixBgAsymKpipi);//0.016,0.034
	LauParameter* backgroundPhiAsym = new LauParameter("backgroundAsymPhi",bgAsymPhi,-1.,1.,fixBgAsymPhi);//0.016,0.034


	fitModel->setNBkgndEvents( nBkgndEvents, backgroundAsym  );
	fitModel->setBkgndDPModels( "qqbar", qqbarModelneg, qqbarModelpos);

	fitModel->setNBkgndEvents( nKpipiBkgndEvents, backgroundKpipiAsym ); //< ===============================
	fitModel->setBkgndDPModels( bkgndNames[1], KpipiModelneg, KpipiModelpos );


 	fitModel->setNBkgndEvents( nPhiBkgndEvents, backgroundPhiAsym ); //< ===============================
	fitModel->setBkgndDPModels( bkgndNames[2], PhiBkgModelneg, PhiBkgModelpos );


	// Switch on/off calculation of asymmetric errors.
	fitModel->useAsymmFitErrors(kFALSE);

	// Randomise initial fit values for the signal mode
	//fitModel->useRandomInitFitPars(kFALSE);

	// Switch on/off Poissonian smearing of total number of events
	fitModel->doPoissonSmearing(kTRUE);

	// Switch on/off Extended ML Fit option
	Bool_t emlFit = ( fitModel->nBkgndClasses() > 0 );
	fitModel->doEMLFit(emlFit);



    /************************************************************************/
    /*** Outputs globais (inicio)                                        ****/

	// Set the names of the files to read/write
	TString dataFile=("");
	TString treeName("fitTree");
	TString rootFileName("");
	TString tableFileName("");
	TString fitToyFileName("");
	TString fitToyTableName("");
	TString splotFileName("");
	if (command == "fit") {
	        dataFile=DataFile;
		rootFileName = resDir+"fit"; 
		rootFileName += "_expt_"; rootFileName += firstExpt;
		rootFileName += "-"; rootFileName += (firstExpt+nExpt-1);
		rootFileName += ".root";
		tableFileName = resDir+"fitResults"; 
                //tableFileName += iFit;
                fitToyTableName += resDir+"ToyMCTable_"; 
		fitToyFileName += resDir+"ToyMC"; 
		//fitToyFileName += iFit;
		fitToyFileName += ".root";
		splotFileName += resDir+"splot_";
		//splotFileName += iFit;
		splotFileName += ".root";
	} else {
                dataFile=resDir+"GenKKpi.root";
		rootFileName = resDir+"gen_";
                //rootFileName += iFit;
                rootFileName += ".root";
                //
		tableFileName = resDir+"genResults_";
                tableFileName += iFit;
	}

    /*** Fim do "Outputs globais"                                        ****/
    /************************************************************************/

        fitModel->writeLatexTable(true);

	// Generate toy from the fitted parameters
	fitModel->compareFitData(5, fitToyFileName, fitToyTableName);

	// Write out per-event likelihoods and sWeights
	//fitModel->writeSPlotData(splotFileName, "splot", kFALSE);

	// Execute the generation/fit
	fitModel->run( command, dataFile, treeName, rootFileName, tableFileName );
        gSystem->Exec("cp   Integ_FF.dat "+resDir);
        gSystem->Exec("cp   Integ_FFNorm.dat "+resDir);

	return EXIT_SUCCESS;
}
