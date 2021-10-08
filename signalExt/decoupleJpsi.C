#include "../common/headers.h"
#include "../common/function.C"
#include "../common/funUtil.h"

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooFormulaVar.h"
#include "RooHistPdf.h"
#include "RooGenericPdf.h"
#include "RooAddPdf.h"
#include "RooPolynomial.h"
#include "RooChi2Var.h"
#include "RooMinimizer.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"
#include "RooPlot.h"
#include "RooFitResult.h"

using namespace RooFit;

const Bool_t  mStorePDF = kFALSE;

const Double_t mTinyNum = 1.e-6;
const Double_t mOffSet  = 0.1;

Int_t    mTextFont  = 42;
Double_t mTextSize  = 0.045;
Int_t    mTextColor = 1;

Double_t mMarkerStyle = 20;
Double_t mMarkerSize  = 0.8;

Double_t mTitleSize    = 0.06;
Double_t mXTitleOffset = 0.95;
Double_t mYTitleOffset = 0.95;
Double_t mLabelSize    = 0.05;
Double_t mTickLength   = 0.02;
Int_t    mXNdivisions  = 210;
Int_t    mYNdivisions  = 208;

Int_t    mLineWidth = 1;
Int_t    cohJpsiColor   = kBlue,        cohJpsiStyle = 1;
Int_t    incohJpsiColor = kViolet-1,    incohJpsiStyle = 1;
Int_t    dissoJpsiColor = kMagenta-4,   dissoJpsiStyle = 1;
Int_t    feeddownJpsiColor = kAzure+10, feeddownJpsiStyle = 1;
Int_t    qedColor = kGreen+3, qedStyle = 1;

const Double_t mJpsi = 3.0969, mPsi = 3.686097;

void decoupleJpsi(Bool_t incHadron = kFALSE, TString hfVetoType="Default", Double_t massLow=2.5, Double_t massHi=3.5, Double_t ptLow=0, Double_t ptHi=3.6)
{
	gStyle->SetOptFit(1111);

	if(!hfVetoType.EqualTo("Default") && !hfVetoType.EqualTo("Tight") && !hfVetoType.EqualTo("Loose") && !hfVetoType.EqualTo("removeHF")){
		cout<<"Please input the correct hfVetoType string: 'Default' OR 'Tight' OR 'Loose' OR 'removeHF'"<<endl;
		return;
	}

	TString histoDir = "../anaData/jpsiHistos";
	TString plotDir  = "jpsiPlots";

	TString fileName = Form("%s/rawSig", histoDir.Data());
	TString dirName  = Form("%s/rawSig", plotDir.Data());

	if(incHadron)
	{
		fileName += ".incHadron";
		dirName  += "_incHadron";
	}

	if(hfVetoType.EqualTo("Tight"))
	{
		fileName += ".tightHF";
		dirName  += "_tightHF";
	}
	else if(hfVetoType.EqualTo("Loose"))
	{
		fileName += ".looseHF";
		dirName  += "_looseHF";
	}
	else if(hfVetoType.EqualTo("removeHF"))
	{
		fileName += ".removeHF";
		dirName  += "_removeHF";
	}

	cout<<"fileName: "<<fileName<<endl;
	cout<<"dirName: " <<dirName<<endl;

	TFile *f = TFile::Open(Form("%s.root", fileName.Data()));

	system(Form("mkdir -p %s", dirName.Data()));
	system(Form("rm -rf %s/*", dirName.Data()));

	TCanvas* c1 = new TCanvas("c1", "c1", 0, 0, 800, 600);
	setPad(0.12, 0.08, 0.07, 0.13);

	TH3D *hMvsPtvsRap = (TH3D *)f->Get("hMvsPtvsRap");

	TH1D *hMass = nullptr, *hPt = nullptr, *hPt_SideBand = nullptr;
	TH1D *hMass_Rap[nRapBins], *hPt_Rap[nRapBins], *hPt_SideBand_Rap[nRapBins], *hPt_LowMassBand_Rap[nRapBins], *hPt_HiMassBand_Rap[nRapBins];
	
	for(Int_t irap=0; irap<nRapBins; irap++)
	{
		Int_t ptBinLow   = 1;
		Int_t ptBinHi    = hMvsPtvsRap->GetNbinsY();
		Int_t massBinLow = hMvsPtvsRap->GetZaxis()->FindBin(mJpsiMassLow  + mTinyNum);
		Int_t massBinHi  = hMvsPtvsRap->GetZaxis()->FindBin(mJpsiMassHi   - mTinyNum);
		Int_t rapBinLow  = hMvsPtvsRap->GetXaxis()->FindBin(mRapLow[irap] + mTinyNum);
		Int_t rapBinHi   = hMvsPtvsRap->GetXaxis()->FindBin(mRapHi[irap]  - mTinyNum);
		
		hPt_Rap[irap]    = (TH1D *)hMvsPtvsRap->ProjectionY(Form("hPt_RapBin%d",   irap), rapBinLow, rapBinHi, massBinLow, massBinHi);
		hMass_Rap[irap]  = (TH1D *)hMvsPtvsRap->ProjectionZ(Form("hMass_RapBin%d", irap), rapBinLow, rapBinHi, ptBinLow,   ptBinHi  );

		massBinLow = hMvsPtvsRap->GetZaxis()->FindBin(mLowMassBandLow + mTinyNum);
		massBinHi  = hMvsPtvsRap->GetZaxis()->FindBin(mLowMassBandHi - mTinyNum);
		hPt_LowMassBand_Rap[irap]   = (TH1D *)hMvsPtvsRap->ProjectionY(Form("hPt_LowMassBand_RapBin%d", irap), rapBinLow, rapBinHi, massBinLow, massBinHi);

		massBinLow = hMvsPtvsRap->GetZaxis()->FindBin(mHiMassBandLow + mTinyNum);
		massBinHi  = hMvsPtvsRap->GetZaxis()->FindBin(mHiMassBandHi - mTinyNum);
		hPt_HiMassBand_Rap[irap]   = (TH1D *)hMvsPtvsRap->ProjectionY(Form("hPt_HiMassBand_RapBin%d", irap), rapBinLow, rapBinHi, massBinLow, massBinHi);

		hPt_SideBand_Rap[irap] = (TH1D *)hPt_LowMassBand_Rap[irap]->Clone(Form("hPt_SideBand_RapBin%d", irap));
		hPt_SideBand_Rap[irap]->Add(hPt_HiMassBand_Rap[irap]);

		if(irap==0)
		{
			hPt          = (TH1D *)hPt_Rap[irap]->Clone("hPt");
			hPt          ->SetTitle(Form("%1.1f < |y| < %1.1f", mRapLow[nRapBins/2], mRapHi[nRapBins-1]));

			hPt_SideBand = (TH1D*)hPt_SideBand_Rap[irap]->Clone("hPt_SideBand");
			hPt_SideBand ->SetTitle(Form("%1.1f < |y| < %1.1f", mRapLow[nRapBins/2], mRapHi[nRapBins-1]));

			hMass        = (TH1D *)hMass_Rap[irap]->Clone("hMass");
			hMass        ->SetTitle(Form("%1.1f < |y| < %1.1f", mRapLow[nRapBins/2], mRapHi[nRapBins-1]));
		}
		else
		{
			hPt          ->Add(hPt_Rap[irap]);
			hPt_SideBand ->Add(hPt_SideBand_Rap[irap]);
			hMass        ->Add(hMass_Rap[irap]);
		}
	}

	const int mRebMass = 1, mRebMass_Rap = 4;
	const int mRebPt   = 1, mRebPt_Rap   = 5;
	hMass       ->Rebin(mRebMass);
	hPt         ->Rebin(mRebPt);
	hPt_SideBand->Rebin(mRebPt);

	for(Int_t irap=0; irap<nRapBins; irap++)
	{
		hMass_Rap[irap]        ->Rebin(mRebMass_Rap);
		hPt_Rap[irap]          ->Rebin(mRebPt_Rap);
		hPt_SideBand_Rap[irap] ->Rebin(mRebPt_Rap);
	}

	TFile *fTemps = TFile::Open("../simulation/effAndTemp/MassPtTemp_AllSpecs.root");

	TF1  *fCohJpsiTemp  = (TF1  *) fTemps->Get("fCohJpsiTemp");
	TH1D *hJpsiMassHist = (TH1D *) fTemps->Get("hCohJpsiMass");
	TH1D *hQEDMassHist  = (TH1D *) fTemps->Get("hLowMassGammaGammaMass");

	RooRealVar mMass("mMass", "m_{#mu#mu} (GeV)", massLow, massHi);
	//RooRealVar cbAlpha("cbAlpha", "cbAlpha", fCohJpsiTemp->GetParameter(1), 0, 10);
	//RooRealVar cbN("cbN", "cbN", fCohJpsiTemp->GetParameter(2), 0, 10);
	//RooRealVar sigmaRatio("sigmaRatio", "sigmaRatio", fCohJpsiTemp->GetParameter(3), 1, 10);
	RooConstVar cbAlpha(    "cbAlpha",     "cbAlpha",    fCohJpsiTemp->GetParameter(1));
	RooConstVar cbN(        "cbN",         "cbN",        fCohJpsiTemp->GetParameter(2));
	RooConstVar sigmaRatio( "sigmaRatio",  "sigmaRatio", fCohJpsiTemp->GetParameter(3));
	RooRealVar  jpsiMu(     "jpsiMu",      "jpsiMu",     3.096, 3.0, 3.2);
	RooRealVar  gausN(      "gausN",       "gausN",      4.6,    -5, 20 );
	RooRealVar  jpsiSigma(  "jpsiSigma",   "jpsiSigma",  0.045, 0, 0.1  );
	RooConstVar massRatio(  "massRatio",   "massRatio",  mPsi/mJpsi     );

	//RooDataHist hJpsiMassRooHist("hJpsiMassRooHist", "hJpsiMassRooHist", mMass, hJpsiMassHist);
	//RooHistPdf  jpsiPdf("jpsiPdf", "jpsiPdf", mMass, hJpsiMassRooHist, 2); // RebinX and interpolation order to make the Jpsi pdf smooth

	RooGenericPdf *jpsiPdf = new RooGenericPdf("jpsiPdf", "jpsiPdf", "ROOT::Math::crystalball_function(mMass,cbAlpha,cbN,jpsiSigma*sigmaRatio,jpsiMu) + gausN*TMath::Gaus(mMass, jpsiMu, jpsiSigma)", RooArgSet(mMass, cbAlpha, cbN, jpsiSigma, sigmaRatio, jpsiMu, gausN));

	//RooGenericPdf *psiPdf = new RooGenericPdf("psiPdf", "psiPdf", "ROOT::Math::crystalball_function(mMass,cbAlpha,cbN,jpsiSigma*massRatio*sigmaRatio,jpsiMu*massRatio) + gausN*TMath::Gaus(mMass, jpsiMu*massRatio, jpsiSigma*massRatio)", RooArgSet(mMass, cbAlpha, cbN, jpsiSigma, sigmaRatio, jpsiMu, massRatio, gausN)); // psiMu = jpsiMu * massRatio; psiSigma = jpsiSigma * massRatio

	TF1 *fQED = new TF1("fQED", "[0] + [1]*x + [2]*x*x +[3]*x*x*x", 0, 5);
	hQEDMassHist->Fit(fQED, "R", "", massLow, massHi); // Using the parameters extracted from simulation to initialize qedPdf

	//RooRealVar mP0("mP0", "mP0", -3.e4, -1.e5, 1.e5);
	//RooRealVar mP1("mP1", "mP1",   3e4, -1.e5, 1.e5);
	//RooRealVar mP2("mP2", "mP2", -7000, -5.e4, 5.e4);
	//RooRealVar mP3("mP3", "mP3",   700, -5.e3, 5.e3);
	RooRealVar mP0(  "mP0", "mP0",  fQED->GetParameter(0), -1, 1);
	RooRealVar mP1(  "mP1", "mP1",  fQED->GetParameter(1), -1, 1);
	RooRealVar mP2(  "mP2", "mP2",  fQED->GetParameter(2), -1, 1);
	RooRealVar mP3(  "mP3", "mP3",  fQED->GetParameter(3), -1, 1);
	RooRealVar nQED( "nQED","nQED", 3.e4,                0, 8.e4);
	RooGenericPdf *qedPdf = new RooGenericPdf("qedPdf", "qedPdf", "mP0 + mP1*mMass + mP2*mMass*mMass + mP3*mMass*mMass*mMass", RooArgSet(mP0, mP1, mP2, mP3, mMass));

	//// directly use QED template from simulation
	//Int_t massBinLow = hQEDMassHist->GetXaxis()->FindBin(massLow + mTinyNum);
	//Int_t massBinHi  = hQEDMassHist->GetXaxis()->FindBin(massHi - mTinyNum);
	//Int_t jpsiMassBinLow = hQEDMassHist->GetXaxis()->FindBin(mJpsiMassLow + mTinyNum);
	//Int_t jpsiMassBinHi  = hQEDMassHist->GetXaxis()->FindBin(mJpsiMassHi - mTinyNum);
	//Double_t mQEDFrac = hQEDMassHist->Integral(jpsiMassBinLow, jpsiMassBinHi)*1./hQEDMassHist->Integral(massBinLow, massBinHi);

	//hQEDMassHist->RebinX(5);
	//RooDataHist hQEDMassRooHist("hQEDMassRooHist", "hQEDMassRooHist", mMass, hQEDMassHist);
	//RooHistPdf  qedPdf("qedPdf", "qedPdf", mMass, hQEDMassRooHist, 2); // RebinX and interpolation order to make the QED pdf smooth 

	RooRealVar nJpsi("nJpsi", "nJpsi", 6.e4, 0, 1.e5);
	//RooRealVar nPsi("nPsi", "nPsi", 2.e3, 0, 1.e4);
	RooAddPdf totMassPdf("totMassPdf", "totMassPdf", RooArgList(*jpsiPdf, *qedPdf), RooArgList(nJpsi, nQED)); 
	//RooAddPdf totMassPdf("totMassPdf", "totMassPdf", RooArgList(jpsiPdf, *qedPdf), RooArgList(nJpsi, nQED)); 
	//RooAddPdf totMassPdf("totMassPdf", "totMassPdf", RooArgList(*jpsiPdf, *psiPdf, qedPdf), RooArgList(nJpsi, nPsi, nQED)); 

	RooDataHist dataMass("dataMass", "dataMass", mMass, hMass); 
	totMassPdf.fitTo(dataMass,Extended(kTRUE),SumW2Error(kTRUE),Hesse(kTRUE),Minos(kFALSE),Save());
	//totMassPdf.fitTo(dataMass,Range(massLow, massHi),Extended(kTRUE),SumW2Error(kTRUE),Hesse(kTRUE),Minos(kFALSE),Save());

	Int_t nFrameMBins = (massHi - massLow) / hMass->GetBinWidth(1);

	c1->cd();
	RooPlot *frameMass = mMass.frame(Range(massLow, massHi), Title(""), Bins(nFrameMBins));
	dataMass.plotOn(frameMass, MarkerStyle(20), MarkerSize(1), MarkerColor(1), LineColor(1), LineWidth(2), DrawOption("pz"));
	totMassPdf.plotOn(frameMass, LineColor(2), LineStyle(1), LineWidth(2));
	//totMassPdf.plotOn(frameMass, Components(RooArgSet(jpsiPdf)), LineColor(kBlue), LineStyle(5), LineWidth(2));
	totMassPdf.plotOn(frameMass, Components(RooArgSet(*jpsiPdf)), LineColor(kBlue), LineStyle(5), LineWidth(2));
	//totMassPdf.plotOn(frameMass, Components(RooArgSet(*psiPdf)), LineColor(kOrange+1), LineStyle(6), LineWidth(2));
	totMassPdf.plotOn(frameMass, Components(RooArgSet(*qedPdf)), LineColor(kMagenta), LineStyle(2), LineWidth(2));

	//cout<<endl;
	//cout<<"******** Print frame ********"<<endl;
	//frameMass->Print();
	//cout<<"******** End ********"<<endl;
	//cout<<endl;
	
	Double_t chi2ndf = frameMass->chiSquare("totMassPdf_Norm[mMass]", "h_dataMass", 9);
	//Double_t chi2ndf = frameMass->chiSquare("","h_dataMass", 9);

	frameMass->Draw() ;
	drawLatex(0.18, 0.84, Form("%1.1f < |y| < %1.1f", mRapLow[nRapBins/2], mRapHi[nRapBins-1]), mTextFont, 0.05, mTextColor);
	drawLatex(0.18, 0.72, Form("#chi^{2}/ndf = %1.1f", chi2ndf),                                mTextFont, mTextSize, mTextColor);
	drawLatex(0.18, 0.66, Form("N_{J/#psi} = %d #pm %d", TMath::Nint(nJpsi.getVal()), TMath::Nint(nJpsi.getError())), mTextFont, mTextSize, mTextColor);
	//drawLatex(0.18, 0.60, Form("N_{#psi(2S)} = %d #pm %d", TMath::Nint(nPsi.getVal()), TMath::Nint(nPsi.getError())), mTextFont, mTextSize, mTextColor);
	drawLatex(0.18, 0.60, Form("N_{QED} = %d #pm %d", TMath::Nint(nQED.getVal()), TMath::Nint(nQED.getError())), mTextFont, mTextSize, mTextColor);
	c1->SaveAs(Form("%s/massSpec.pdf", dirName.Data()));

	//----------------------------------------------------------------------------------------------------------------------------------------------------
	//----------------------------------------------------------------------------------------------------------------------------------------------------
	// *** fit to pt distribution ***
	//----------------------------------------------------------------------------------------------------------------------------------------------------
	//----------------------------------------------------------------------------------------------------------------------------------------------------
	
	TCanvas* c2 = new TCanvas("c2", "c2", 900, 900);
		
	TH1D *hCohJpsiPtHist      = (TH1D *)fTemps->Get("hCohJpsiPt"          );
	TH1D *hInCohJpsiPtHist    = (TH1D *)fTemps->Get("hInCohJpsiPt"        );
	TH1D *hFeeddownJpsiPtHist = (TH1D *)fTemps->Get("hCohPsi2SFeeddownPt" );
	TH1D *hQEDPtHist          = (TH1D *)fTemps->Get("hLowMassGammaGammaPt");

	RooRealVar mPt("mPt", "p_{T} (GeV)", ptLow, ptHi);

	hCohJpsiPtHist->RebinX(mRebPt);
	RooDataHist hCohJpsiPtRooHist(     "hCohJpsiPtRooHist",      "hCohJpsiPtRooHist",      mPt, hCohJpsiPtHist);
	RooHistPdf  cohJpsiPdf(            "cohJpsiPdf",             "cohJpsiPdf",             mPt, hCohJpsiPtRooHist,      0);

	hInCohJpsiPtHist->RebinX(mRebPt);
	RooDataHist hInCohJpsiPtRooHist(   "hInCohJpsiPtRooHist",    "hInCohJpsiPtRooHist",    mPt, hInCohJpsiPtHist);
	RooHistPdf  incohJpsiPdf(          "incohJpsiPdf",           "incohJpsiPdf",           mPt, hInCohJpsiPtRooHist,    0);

	hFeeddownJpsiPtHist->RebinX(mRebPt);
	RooDataHist hFeeddownJpsiPtRooHist("hFeeddownJpsiPtRooHist", "hFeeddownJpsiPtRooHist", mPt, hFeeddownJpsiPtHist);
	RooHistPdf  feeddownJpsiPdf(       "feeddownJpsiPdf",        "feeddownJpsiPdf",        mPt, hFeeddownJpsiPtRooHist, 0);

	//hQEDPtHist->RebinX(mRebPt);
	//RooDataHist hQEDPtRooHist("hQEDPtRooHist", "hQEDPtRooHist", mPt, hQEDPtHist);
	//RooHistPdf  qedPtPdf("qedPtPdf", "qedPtPdf", mPt, hQEDPtRooHist, 0);
	RooDataHist hQEDPtRooHist("hQEDPtRooHist", "hQEDPtRooHist", mPt, hPt_SideBand);
	RooHistPdf  qedPtPdf(     "qedPtPdf",      "qedPtPdf",      mPt, hQEDPtRooHist, 0);

	//RooConstVar bpd("bpd", "bpd", 1.79);
	//RooConstVar npd("npd", "npd", 3.58);
	RooRealVar bpd("bpd", "bpd", 1.79, 0, 10);
	RooRealVar npd("npd", "npd", 3.58, 0, 10);
	RooGenericPdf *dissoJpsiPdf = new RooGenericPdf("dissoJpsiPdf", "dissoJpsiPdf", "mPt*TMath::Power(1+bpd/npd*mPt*mPt, -npd)", RooArgSet(mPt, bpd, npd));

	fQED->SetParameters(mP0.getVal(), mP1.getVal(), mP2.getVal(), mP3.getVal());
	Double_t mQEDFrac = fQED->Integral(mJpsiMassLow, mJpsiMassHi) / fQED->Integral(massLow, massHi);

	RooRealVar  nCohJpsi(     "nCohJpsi",       "nCohJpsi",      3.6e4, 0, 1.e5);
	RooRealVar  nInCohJpsi(   "nInCohJpsi",     "nInCohJpsi",    1.2e4, 0, 5.e4);
	RooRealVar  nFeeddownJpsi("nFeeddownJpsi",  "nFeeddownJpsi", 6.e2,  0, 1.e4);
	RooRealVar  nDissoJpsi(   "nDissoJpsi",     "nDissoJpsi",    1.e4,  0, 5.e4);
	RooRealVar  nQEDBg(       "nQEDBg",         "nQEDBg",        nQED.getVal()*mQEDFrac, 0, 5.e4);
	//RooConstVar nQEDBg(       "nQEDBg",         "nQEDBg",        nQED.getVal()*mQEDFrac);

	RooAddPdf totPtPdf("totPtPdf", "totPtPdf", 
			RooArgList(cohJpsiPdf, incohJpsiPdf, feeddownJpsiPdf, *dissoJpsiPdf, qedPtPdf),
			RooArgList(nCohJpsi,   nInCohJpsi,   nFeeddownJpsi,    nDissoJpsi,   nQEDBg));

	RooDataHist dataPt("dataPt", "dataPt", mPt, hPt); 

	//cout<<endl;
	//cout<<"******** Print QED yield ********"<<endl;
	//cout<<nQED.getVal()<<"        "<<mQEDFrac<<endl;
	//cout<<"nQED(SideBand): "<<hPt_SideBand->Integral()<<endl;
	//cout<<"******** Print end ********"<<endl;
	//cout<<endl;

	totPtPdf.fitTo(dataPt,Extended(kTRUE),SumW2Error(kTRUE),Hesse(kTRUE),Minos(kFALSE),Save());

	Int_t nFramePtBins = (ptHi - ptLow)/hPt->GetBinWidth(1);
	cout<<nFramePtBins<<endl;

	RooPlot *framePt = mPt.frame(Range(ptLow, ptHi), Title(""), Bins(nFramePtBins));
	dataPt  .plotOn(framePt, MarkerStyle(20), MarkerSize(0.6), MarkerColor(1), LineColor(1), LineWidth(mLineWidth), DrawOption("pz"));
	totPtPdf.plotOn(framePt, LineColor(1), LineStyle(1), LineWidth(1));
	totPtPdf.plotOn(framePt, Components(RooArgSet(cohJpsiPdf)),      LineColor(cohJpsiColor),      LineStyle(cohJpsiStyle),      LineWidth(mLineWidth));
	totPtPdf.plotOn(framePt, Components(RooArgSet(incohJpsiPdf)),    LineColor(incohJpsiColor),    LineStyle(incohJpsiStyle),    LineWidth(mLineWidth));
	totPtPdf.plotOn(framePt, Components(RooArgSet(*dissoJpsiPdf)),   LineColor(dissoJpsiColor),    LineStyle(dissoJpsiStyle),    LineWidth(mLineWidth));
	totPtPdf.plotOn(framePt, Components(RooArgSet(feeddownJpsiPdf)), LineColor(feeddownJpsiColor), LineStyle(feeddownJpsiStyle), LineWidth(mLineWidth));
	totPtPdf.plotOn(framePt, Components(RooArgSet(qedPtPdf)),        LineColor(qedColor),          LineStyle(qedStyle),          LineWidth(mLineWidth));

	cout<<endl;
	cout<<"******** Print frame ********"<<endl;
	framePt->Print();
	cout<<"******** End ********"<<endl;
	cout<<endl;
	
	//chi2ndf = framePt->chiSquare("totPtPdf_Norm[mPt]", "h_dataPt", 6);
	chi2ndf = framePt->chiSquare("totPtPdf_Norm[mPt]", "h_dataPt", 7);

	c2 -> Divide(1,2);
	c2 -> cd(1);
	gPad->SetPad(0.0,0.25,1.0,0.96);
	gPad->SetBottomMargin(0);
	gPad->SetRightMargin(0.05);
	gPad->SetTopMargin(0);
	//setPad(0.12, 0.08, 0.07, 0.13);

	gPad->SetLogy(1);
	framePt->GetYaxis()->SetRangeUser(0.5, hPt->GetMaximum()*5);
	framePt->Draw() ;
	drawLatex(0.18, 0.82, Form("%1.1f < |y| < %1.1f",  mRapLow[nRapBins/2], mRapHi[nRapBins-1]), mTextFont, 0.05,      mTextColor);
	drawLatex(0.15, 0.92, "CMS Pb-Pb #sqrt{s_{NN}} = 5.02 TeV UPC",                              42,       0.06,      mTextColor);
	//drawLatex(0.24, 0.78, Form("#chi^{2}/ndf = %1.1f", chi2ndf),                                 mTextFont, mTextSize, mTextColor);

	//	TF1  *cohJpsi      = new TF1("cohJpsi",      "1", 0, 1);
	//	TF1  *incohJpsi    = new TF1("incohJpsi",    "1", 0, 1);
	//	TF1  *dissoJpsi    = new TF1("dissoJpsi",    "1", 0, 1);
	//	TF1  *feeddownJpsi = new TF1("feeddownJpsi", "1", 0, 1);
	//	TF1  *qed          = new TF1("qed",          "1", 0, 1);
	//
	//	setFun(cohJpsi,       cohJpsiColor,      mLineWidth, cohJpsiStyle     );
	//	setFun(incohJpsi,     incohJpsiColor,    mLineWidth, incohJpsiStyle   );
	//	setFun(dissoJpsi,     dissoJpsiColor,    mLineWidth, dissoJpsiStyle   );
	//	setFun(feeddownJpsi,  feeddownJpsiColor, mLineWidth, feeddownJpsiStyle);
	//	setFun(qed, qedColor, mLineWidth, qedStyle);

	TLegend  *leg =  new TLegend(0.54, 0.45, 0.84, 0.80);
	leg->SetFillStyle(0);
	leg->SetFillColor(0);
	leg->SetTextFont(mTextFont);
	leg->SetTextSize(0.042);
	leg->AddEntry(framePt->findObject("h_dataPt"),            "Data",        "p");
	leg->AddEntry(framePt->findObject("totPtPdf_Norm[mPt]"),   Form("Total fit: #chi^{2}/ndf = %1.1f", chi2ndf),   "l");
	
	const TString curveName[5]  = {"totPtPdf_Norm[mPt]_Comp[cohJpsiPdf]","totPtPdf_Norm[mPt]_Comp[incohJpsiPdf]", "totPtPdf_Norm[mPt]_Comp[dissoJpsiPdf]", "totPtPdf_Norm[mPt]_Comp[feeddownJpsiPdf]","totPtPdf_Norm[mPt]_Comp[qedPtPdf]"};
	const TString curveTitle[5] = {"Coherent J/#psi", "Incoherent J/#psi", "Incoherent J/#psi with disso.", "Coherent #psi' #rightarrow J/#psi+X", "#gamma#gamma #rightarrow #mu#mu"};
	for(int icv=0; icv<5; icv++) 
	{
		leg->AddEntry(framePt->findObject( curveName[icv]), curveTitle[icv], "l");
	}

	//leg->AddEntry(cohJpsi,      "Coherent J/#psi",                     "l");
	//leg->AddEntry(incohJpsi,    "Incoherent J/#psi",                   "l");
	//leg->AddEntry(dissoJpsi,    "Incoherent J/#psi with",              "l");
	//leg->AddEntry((TObject*)0,  "nucleon dissociation",                "" );
	//leg->AddEntry(feeddownJpsi, "Coherent #psi' #rightarrow J/#psi+X", "l");
	//leg->AddEntry(qed,          "#gamma#gamma #rightarrow #mu#mu",     "l");
	leg->Draw("same");

	c2	->cd(2);
	gPad->SetPad(0.0,0.0,1.0,0.25);
	gPad->SetTopMargin(0);
	gPad->SetRightMargin(0.05);
	gPad->SetBottomMargin(0.45);

	//calculate the fiting pull: (Data-FitCurve)/#Sigma_{Data}
	RooHist *hpull_pt = framePt->pullHist("h_dataPt", "totPtPdf_Norm[mPt]");
	//TH1D* histPull_pt = (TH1D*) hpull_pt->GetHistogram();

	RooPlot *framePtPull = mPt.frame(Range(ptLow, ptHi), Title(""), Bins(nFramePtBins));
	
	hpull_pt ->SetMarkerStyle(24);
	hpull_pt ->SetMarkerSize(0.6);
	hpull_pt ->SetMarkerColor(1);
	hpull_pt ->SetLineColor(1);
	hpull_pt ->SetLineWidth(1);

	framePtPull ->addPlotable(hpull_pt, "pz");
	framePtPull ->GetYaxis()->SetRangeUser(-14.0, 14.0);
	framePtPull ->SetYTitle("#frac{Data-Fit}{#sigma_{Data}}");
	framePtPull ->SetXTitle("p_{T} (GeV/c)");
	framePtPull ->GetYaxis()->CenterTitle();
	framePtPull ->GetYaxis()->SetNdivisions(6);
	framePtPull ->GetYaxis()->SetTitleSize(0.15);
	framePtPull ->GetYaxis()->SetTitleOffset(0.33);
	framePtPull ->GetYaxis()->SetLabelSize(0.10);
	//framePtPull ->GetYaxis()->SetLabelFont(20);
	framePtPull ->GetXaxis()->SetTitleSize(0.20);
	framePtPull ->GetXaxis()->SetTitleOffset(0.95);
	framePtPull ->GetXaxis()->SetLabelSize(0.16);
	framePtPull ->GetXaxis()->SetLabelFont(40);
	framePtPull ->SetTickLength(0.08);

	framePtPull ->Draw() ;

	//	
	//	TH2D* htem2d_4pt = new TH2D("htem2d_4pt", "", nFramePtBins, ptLow, ptHi, 10, -0.001, 0.001);
	//	htem2d_4pt->SetYTitle("Data/Fit");
	//	htem2d_4pt->SetXTitle("p_{T} (GeV/c)");
	//	htem2d_4pt->GetYaxis()->SetNdivisions(4);
	//	htem2d_4pt->GetYaxis()->SetTitleSize(0.22);
	//	htem2d_4pt->GetYaxis()->SetTitleOffset(0.25);
	//	htem2d_4pt->GetYaxis()->SetLabelSize(0.15);
	//	htem2d_4pt->GetYaxis()->SetLabelFont(40);
	//	htem2d_4pt->GetXaxis()->SetTitleSize(0.20);
	//	htem2d_4pt->GetXaxis()->SetTitleOffset(0.95);
	//	htem2d_4pt->GetXaxis()->SetLabelSize(0.16);
	//	htem2d_4pt->GetXaxis()->SetLabelFont(40);
	//	htem2d_4pt->SetTickLength(0.08);
	//	htem2d_4pt->Draw() ;
	//
	//	histPull_pt->Draw("lsame");
	//	histPull_pt->Draw("pesame");
	//
	c2->SaveAs(Form("%s/ptSpec.pdf", dirName.Data()));

	//	//RooHist is TGraphError
	//	for(int i=0; i<hpull_pt->GetN(); i++)
	//	{
	//		double x,y,xer,yer;
	//		hpull_pt->GetPoint(i,x,y);
	//		xer = hpull_pt->GetErrorX(i);
	//		yer = hpull_pt->GetErrorY(i);
	//
	//		cout<<"i: "<<i<<" x:"<<x<<"+/-"<<xer<<" y:"<<y<<"+/-"<<yer<<endl;
	//	}

	cout << "End of program !" << endl;
}
