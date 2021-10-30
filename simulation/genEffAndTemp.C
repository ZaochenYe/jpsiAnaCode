#include "../common/headers.h"
#include "../common/function.C"
#include "../common/funUtil.h"

const Double_t mTinyNum = 1.e-6;

Int_t    mGenMarker = 20;
Int_t    mGenColor = 1;
Int_t    mGenWidth = 2;
Double_t mGenSize = 0.8;

//------------------------------------------------------ //bremsstrahlung tail
Double_t DoubleNGaus(Double_t *x, Double_t *par);
Double_t CrystalBall(Double_t *x, Double_t *par);
Double_t CrystalBallAndGaus(Double_t *x, Double_t *par);

//const Double_t mFitMassLow[nRapBins] = {2.3, 2.3, 2.3, 2.4, 2.4, 2.6, 2.6, 2.4, 2.4, 2.3, 2.3, 2.3};

void genEffAndTemp()
{
	gStyle->SetOptFit(1111);
	Int_t mFont = 42;

	const Int_t   nSpecs            = 9;
	const TString specName[nSpecs]  = {"CohJpsi", "CohJpsi_0n0n", "CohJpsi_0nXn", "CohJpsi_XnXn", "InCohJpsi", "CohPsi2SFeeddown", "CohPsi2S", "InCohPsi2S", "LowMassGammaGamma"};
	const TString specTitle[nSpecs] = {"Coherent J/#psi", "Coherent J/#psi (0n0n)", "Coherent J/#psi (0nXn)", "Coherent J/#psi (XnXn)", "Incoherent J/#psi", "Coherent #psi(2S) #rightarrow J/#psi + X", "Coherent #psi(2S)", "Incoherent #psi(2S)", "#gamma#gamma#rightarrow#mu#mu"};
	const Double_t mMass[nSpecs-1]  = {3.096, 3.096, 3.096, 3.096, 3.096, 3.096, 3.686, 3.686};

	//const Int_t nSpecs = 2;
	//TString specName[nSpecs]  = {"CohJpsi", "LowMassGammaGamma"};
	//TString specTitle[nSpecs] = {"Coherent J/#psi", "#gamma#gamma#rightarrow#mu#mu"};
	//const Double_t    mMass[nSpecs-1] = {3.096};

	TFile* f[nSpecs];
	TH3D *hMvsPtvsRap_Gen[nSpecs];
	TH3D *hMvsPtvsRap_woEvtSel[nSpecs];
	TH3D *hMvsPtvsRap[nSpecs];
	TH1D *hRap_Gen[nSpecs];
	TH1D *hRap_woEvtSel[nSpecs];
	TH1D *hRap[nSpecs];

	TH1D *hEffvsRap_woEvtSel[nSpecs];
	TH1D *hEffvsRap[nSpecs];
	TH1D *hEvtSelEffvsRap[nSpecs];

	TH1D *hMass[nSpecs];
	TH1D *hPt[nSpecs];
	TH1D *hMass_Rap[nSpecs][nRapBins];
	TH1D *hPt_Rap[nSpecs][nRapBins];

	for(Int_t i=0; i<nSpecs; i++)
	{
		f[i] = TFile::Open(Form("mcHistos/dimuonHistos.%s.root", specName[i].Data()));
		hMvsPtvsRap_Gen[i]      = (TH3D *)f[i]->Get("hMvsPtvsRap_Gen");      //no smearing
		hMvsPtvsRap_woEvtSel[i] = (TH3D *)f[i]->Get("hMvsPtvsRap_woEvtSel"); //with smearing, no event level cut
		hMvsPtvsRap[i]          = (TH3D *)f[i]->Get("hMvsPtvsRap");          //with smearing, with event level cut

		Int_t massBinLow, massBinHi;
		if(i<nSpecs-1)
		{
			massBinLow = 1;
			massBinHi  = hMvsPtvsRap_Gen[i]->GetNbinsZ(); 
		}
		else
		{
			massBinLow = hMvsPtvsRap_Gen[i]->GetZaxis()->FindBin( mJpsiMassLow + mTinyNum );
			massBinHi  = hMvsPtvsRap_Gen[i]->GetZaxis()->FindBin( mJpsiMassHi  - mTinyNum );
		}

		hRap_Gen[i]           = (TH1D *)hMvsPtvsRap_Gen[i]      ->ProjectionX( Form("hRap_Gen_%s",      specName[i].Data()), 0, -1, massBinLow, massBinHi );
		hRap_woEvtSel[i]      = (TH1D *)hMvsPtvsRap_woEvtSel[i] ->ProjectionX( Form("hRap_woEvtSel_%s", specName[i].Data()), 0, -1, massBinLow, massBinHi );
		hRap[i]               = (TH1D *)hMvsPtvsRap[i]          ->ProjectionX( Form("hRap_%s",          specName[i].Data()), 0, -1, massBinLow, massBinHi );
		
		//calculate efficiencies
		hEffvsRap_woEvtSel[i] = (TH1D *)hRap_woEvtSel[i]->Clone(Form("hEffvsRap_woEvtSel_%s", specName[i].Data()));
		hEffvsRap_woEvtSel[i] ->Divide( hRap_woEvtSel[i], hRap_Gen[i], 1, 1, "B");
		hEffvsRap_woEvtSel[i] ->SetTitle( specTitle[i].Data() );
		hEffvsRap_woEvtSel[i] ->GetYaxis()->SetTitle("Efficiency");

		hEffvsRap[i]          = (TH1D *)hRap[i]->Clone(Form("hEffvsRap_%s", specName[i].Data()));
		hEffvsRap[i]          ->Divide( hRap[i], hRap_Gen[i], 1, 1, "B");
		hEffvsRap[i]          ->SetTitle( specTitle[i].Data());
		hEffvsRap[i]          ->GetYaxis()->SetTitle("Efficiency");

		hEvtSelEffvsRap[i]    = (TH1D *)hRap[i]->Clone(Form("hEvtSelEffvsRap_%s", specName[i].Data())); 
		hEvtSelEffvsRap[i]    ->Divide(hRap[i], hRap_woEvtSel[i], 1, 1, "B");
		hEvtSelEffvsRap[i]    ->SetTitle(specTitle[i].Data());
		hEvtSelEffvsRap[i]    ->GetYaxis()->SetTitle("Efficiency");

		for(Int_t irap=0; irap<nRapBins; irap++)
		{
			Int_t rapBinLow  = hMvsPtvsRap[i]->GetXaxis()->FindBin(mRapLow[irap] + mTinyNum);
			Int_t rapBinHi   = hMvsPtvsRap[i]->GetXaxis()->FindBin(mRapHi[irap]  - mTinyNum);
			hPt_Rap[i][irap] = (TH1D *)hMvsPtvsRap[i]->ProjectionY(Form("h%sPt_RapBin%d", specName[i].Data(), irap), rapBinLow, rapBinHi, massBinLow, massBinHi);
			
			if(i<nSpecs-1)
				hPt_Rap[i][irap]->SetTitle(Form("%1.1f < y < %1.1f", mRapLow[irap], mRapHi[irap]));
			else
				hPt_Rap[i][irap]->SetTitle(Form("%1.1f < y < %1.1f, %1.2f < mass < %1.2f", mRapLow[irap], mRapHi[irap], mJpsiMassLow, mJpsiMassHi));

			hMass_Rap[i][irap] = (TH1D *)hMvsPtvsRap[i]->ProjectionZ(Form("h%sMass_RapBin%d", specName[i].Data(), irap), rapBinLow, rapBinHi, 0, -1);
			hMass_Rap[i][irap]->SetTitle(Form("%1.1f < y < %1.1f", mRapLow[irap], mRapHi[irap]));

			if(irap==0)
			{
				hPt[i] = (TH1D *)hPt_Rap[i][irap]->Clone(Form("h%sPt", specName[i].Data()));
				
				if(i<nSpecs-1)
					hPt[i]->SetTitle(Form("%1.1f < |y| < %1.1f", mRapLow[nRapBins/2], mRapHi[nRapBins-1]));
				else
					hPt[i]->SetTitle(Form("%1.1f < |y| < %1.1f, %1.2f < mass < %1.2f", mRapLow[nRapBins/2], mRapHi[nRapBins-1], mJpsiMassLow, mJpsiMassHi));

				hMass[i] = (TH1D *)hMass_Rap[i][irap]->Clone(Form("h%sMass", specName[i].Data()));
				hMass[i] ->SetTitle(Form("%1.1f < |y| < %1.1f",  mRapLow[nRapBins/2], mRapHi[nRapBins-1]));
			}
			else
			{
				hPt[i]   ->Add( hPt_Rap[i][irap]   );
				hMass[i] ->Add( hMass_Rap[i][irap] );
			}
		}
	}

	//----------------------------------------------------------------------------------------------------------------------------------------------------------
	// write down efficiencies
	//----------------------------------------------------------------------------------------------------------------------------------------------------------
	TString dir = "effAndTemp";
	system(Form("mkdir -p %s", dir.Data()));
	system(Form("rm -rf %s/*", dir.Data()));

	TFile *fOut = new TFile(Form("%s/Efficiency_AllSpecs.root", dir.Data()), "recreate");
	fOut->cd();
	
	for(Int_t i=0; i<nSpecs; i++)
	{
		hEffvsRap_woEvtSel[i]->Write();
		hEffvsRap[i]->Write();
		hEvtSelEffvsRap[i]->Write();
	}
	
	fOut->Close();
	//----------------------------------------------------------------------------------------------------------------------------------------------------------

	//----------------------------------------------------------------------------------------------------------------------------------------------------------
	//draw efficiencies and templates
	//----------------------------------------------------------------------------------------------------------------------------------------------------------
	TCanvas* c1 = new TCanvas("c1", "c1", 0, 0, 800, 600);

	TCanvas* c2 = new TCanvas("c2", "c2", 0, 0, 1200, 450);
	c2->Divide(2, 1);

	TCanvas* c3 = new TCanvas("c3", "c3", 0, 0, 800, 600);
	c3->Divide(2, 2);

	Double_t xPos = 0.15;
	Double_t yPos = 0.84;

	TLegend* leg1 = new TLegend(0.36, 0.65, 0.6, 0.86);
	leg1->SetFillStyle(0);
	leg1->SetTextSize(0.06);

	TLegend* leg2 = new TLegend(0.36, 0.72, 0.6, 0.86);
	leg1->SetFillStyle(0);
	leg1->SetTextSize(0.06);

	TF1 *fTemp[nSpecs];
	TF1 *fTemp_Rap[nSpecs][nRapBins];
	for(Int_t i=0; i<nSpecs; i++)
	{
		c2->cd(1);
		gPad->SetLogy(1);
		setHisto(hRap_Gen[i], 20, 1, 1, 1, 2);
		setHisto(hRap_woEvtSel[i], 24, 1, 2, 2, 2);
		setHisto(hRap[i], 25, 1, 4, 4, 2);
		hRap_Gen[i]      ->GetYaxis()->SetTitle("Entries");
		hRap_Gen[i]      ->SetMinimum(1);
		hRap_Gen[i]      ->Draw("p");
		hRap_woEvtSel[i] ->Draw("psame");
		hRap[i]          ->Draw("psame");
		
		if(i==0)
		{
			leg1->AddEntry( hRap_Gen[i],      "GEN",                                                                  "pl");
			leg1->AddEntry( hRap_woEvtSel[i], "#varepsilon_{reco}#times#varepsilon_{trig}",                           "pl");
			leg1->AddEntry( hRap[i],          "#varepsilon_{reco}#times#varepsilon_{trig}#times#varepsilon_{evtSel}", "pl");
		}
		leg1->Draw("same");
		
		if(     specName[i].EqualTo("CohPsi2SFeeddown"))  drawLatex(0.27, 0.95, Form("%s", specTitle[i].Data()), mFont, 0.06, 1);
		else if(specName[i].EqualTo("LowMassGammaGamma")) drawLatex(0.43, 0.95, Form("%s", specTitle[i].Data()), mFont, 0.06, 1);
		else                                              drawLatex(0.38, 0.95, Form("%s", specTitle[i].Data()), mFont, 0.06, 1);

		c2->cd(2);
		gPad->SetLogy(1);
		setHisto(hEffvsRap_woEvtSel[i], 24, 1, 2, 2, 2);
		setHisto(hEffvsRap[i],          25, 1, 4, 4, 2);

		hEffvsRap_woEvtSel[i] ->SetMaximum(1);
		hEffvsRap_woEvtSel[i] ->Draw("p");
		hEffvsRap[i]          ->Draw("psame");

		if(i==0)
		{
			leg2->AddEntry(hEffvsRap_woEvtSel[i], "#varepsilon_{reco}#times#varepsilon_{trig}",                           "pl");
			leg2->AddEntry(hEffvsRap[i],          "#varepsilon_{reco}#times#varepsilon_{trig}#times#varepsilon_{evtSel}", "pl");
		}
		leg2->Draw("same");
		
		if(     specName[i].EqualTo("CohPsi2SFeeddown") ) drawLatex(0.27, 0.95, Form("%s", specTitle[i].Data()), mFont, 0.06, 1);
		else if(specName[i].EqualTo("LowMassGammaGamma")) drawLatex(0.43, 0.95, Form("%s", specTitle[i].Data()), mFont, 0.06, 1);
		else                                              drawLatex(0.38, 0.95, Form("%s", specTitle[i].Data()), mFont, 0.06, 1);

		c2->SaveAs(Form("%s/EffvsRap_%s.pdf", dir.Data(), specName[i].Data()));
		c2->SaveAs(Form("%s/EffvsRap_%s.png", dir.Data(), specName[i].Data()));
		//----------------------------------------------------------------------------------------------------------------------------------------------------------
		//----------------------------------------------------------------------------------------------------------------------------------------------------------


		//----------------------------------------------------------------------------------------------------------------------------------------------------------
		// working for templates
		//----------------------------------------------------------------------------------------------------------------------------------------------------------
		TString tempDir = Form("%s/%sTemp", dir.Data(), specName[i].Data());
		system(Form("mkdir -p %s", tempDir.Data()));

		if(i<nSpecs-1)
		{
			//fTemp[i] = new TF1(Form("f%sTemp",specName[i].Data()), "([0]*TMath::Gaus(x,[1],[2],1) + [3]*TMath::Gaus(x,[1],[2]*[4],1))*[5]", 0, 5);
			//fTemp[i]->SetParNames("N1","#mu","#sigma1","N2","#sigma2/#sigma1","binWidth");
			//fTemp[i]->SetParameters(0.8, mMass[i], 0.04, 0.2, 1.5, hMass[i]->GetBinWidth(1));
			//fTemp[i]->FixParameter(5, hMass[i]->GetBinWidth(1));

			//fTemp[i] = new TF1(Form("f%sTemp",specName[i].Data()), "[0]*ROOT::Math::crystalball_function(x,[1],[2],[3],[4])", 0, 5);
			//fTemp[i]->SetParNames("N","#alpha","n","#sigma","#mu");
			//fTemp[i]->SetParameters(0.1, 2, 5, 0.05, mMass[i]);

			fTemp[i] = new TF1(Form("f%sTemp", specName[i].Data()), "[0]*(ROOT::Math::crystalball_function(x,[1],[2],[3]*[6],[4]) + [5]*TMath::Gaus(x, [4], [6], 0))", 0, 5);
			fTemp[i]->SetParNames("N","#alpha","n","#sigma_{cb}/#sigma_{gaus}","#mu","N_{gaus}","#sigma_{gaus}");
			fTemp[i]->SetParameters(0.02, 3, 6, 1.5, mMass[i], 4, 0.04);
			//fTemp[i]->FixParameter(5, 1);
		}
		else
		{
			// try to parameterize gg->mumu mass shape
			fTemp[i] = new TF1(Form("f%sTemp", specName[i].Data()), "[0]*(x-[1])/(TMath::Exp([2]*(x-[1])+[3]/(x-[1]))+[4])", 0, 5);
			fTemp[i]->SetParameters(0.1, 2.5, 1, -1, 1);
		}

		fTemp[i]->SetTitle(Form("%1.1f < |y| < %1.1f", mRapLow[nRapBins/2], mRapHi[nRapBins-1]));
		fTemp[i]->SetNpx(1000);

		c2->cd(1);
		setHisto(hMass[i], 20, 0.6, 1, 1, 2);
		hMass[i]->Scale(1./hMass[i]->GetEntries());

		if(specName[i].Contains("Jpsi"))
		{
			hMass[i]->GetXaxis()->SetRangeUser(2.5, 4);
		}
		else if(specName[i].Contains("Psi2S"))
		{
			if(specName[i].Contains("Feeddown")) hMass[i]->GetXaxis()->SetRangeUser(2.5, 4);
			else                                 hMass[i]->GetXaxis()->SetRangeUser(3, 4.5);
		}
		
		hMass[i]->GetXaxis()->SetTitleOffset(0.9);
		hMass[i]->GetYaxis()->SetTitle("a.u.");
		if(i<nSpecs-1)
		{
			gPad->SetLogy(1);
			hMass[i]->Fit(fTemp[i], "RQ0", "", mMass[i]-0.4, mMass[i]+0.4);
			hMass[i]->Draw("p");
			fTemp[i]->Draw("same");
		}
		else
		{
			gPad->SetLogy(0);
			//hMass[i]->Fit(fTemp[i], "RQ0", "", 2.3, 5);
			hMass[i]->Draw("p");
			//fTemp[i]->Draw("same");
		}
		drawLatex(0.39, 0.95, Form("%1.1f < |y| < %1.1f", mRapLow[nRapBins/2], mRapHi[nRapBins-1]), mFont, 0.06, 1);

		c2->cd(2);
		gPad->SetLogy(1);
		hPt[i]->Scale(1./hPt[i]->GetEntries());
		setHisto(hPt[i], 20, 0.6, 1, 1, 2);
		if(specName[i].Contains("InCoh"))
		{
			hPt[i]->GetXaxis()->SetRangeUser(0, 2);
		}
		else if(specName[i].Contains("Feeddown"))
		{
			hPt[i]->GetXaxis()->SetRangeUser(0, 1);
		}
		else
		{
			hPt[i]->GetXaxis()->SetRangeUser(0, 0.5);
		}
		hPt[i]->GetYaxis()->SetRangeUser(0.5, hPt[i]->GetMaximum()*2);
		hPt[i]->GetXaxis()->SetTitleOffset(0.9);
		hPt[i]->GetYaxis()->SetTitle("a.u.");
		hPt[i]->Draw("p");
		drawLatex(0.39, 0.95, Form("%1.1f < |y| < %1.1f", mRapLow[nRapBins/2], mRapHi[nRapBins-1]), mFont, 0.06, 1);

		c2->SaveAs(Form("%s/%sTemp.pdf", tempDir.Data(), specName[i].Data()));
		c2->SaveAs(Form("%s/%sTemp.png", tempDir.Data(), specName[i].Data()));

		for(Int_t irap=0; irap<nRapBins; irap++)
		{
			if(i<nSpecs-1)
			{
				//fTemp_Rap[i][irap] = new TF1(Form("f%sTemp_RapBin%d",specName[i].Data(),irap), "([0]*TMath::Gaus(x,[1],[2],1) + [3]*TMath::Gaus(x,[1],[2]*[4],1))*[5]", 0, 5);
				//fTemp_Rap[i][irap]->SetParNames("N1","#mu","#sigma1","N2","#sigma2/#sigma1","binWidth");
				//fTemp_Rap[i][irap]->SetParameters(0.8, mMass[i], 0.04, 0.2, 1.5, hMass_Rap[i][irap]->GetBinWidth(1));
				//fTemp_Rap[i][irap]->FixParameter(5, hMass_Rap[i][irap]->GetBinWidth(1));

				//fTemp_Rap[i][irap] = new TF1(Form("f%sTemp_RapBin%d",specName[i].Data(),irap), "[0]*ROOT::Math::crystalball_function(x,[1],[2],[3],[4])", 0, 5);
				//fTemp_Rap[i][irap]->SetTitle(Form("%1.1f < y < %1.1f", mRapLow[irap], mRapHi[irap]));
				//fTemp_Rap[i][irap]->SetParNames("N","#alpha","n","#sigma","#mu");
				//fTemp_Rap[i][irap]->SetParameters(0.1, 2, 5, 0.05, mMass[i]);
				//fTemp_Rap[i][irap]->FixParameter(2, fTemp[i]->GetParameter(2));

				fTemp_Rap[i][irap] = new TF1(Form("f%sTemp_RapBin%d", specName[i].Data(),irap), "[0]*(ROOT::Math::crystalball_function(x,[1],[2],[3]*[6],[4]) + [5]*TMath::Gaus(x, [4], [6], 0))", 0, 5);
				fTemp_Rap[i][irap]->SetTitle(Form("%1.1f < y < %1.1f", mRapLow[irap], mRapHi[irap]));
				fTemp_Rap[i][irap]->SetParNames("N","#alpha","n","#sigma_{cb}/#sigma_{gaus}","#mu","N_{gaus}","#sigma_{Gaus}");
				fTemp_Rap[i][irap]->SetParameters(0.02, 3, 6, 1.5, mMass[i], 4, 0.04);
				fTemp_Rap[i][irap]->FixParameter(1, fTemp[i]->GetParameter(1));
				fTemp_Rap[i][irap]->FixParameter(2, fTemp[i]->GetParameter(2));
				//fTemp_Rap[i][irap]->FixParameter(5, 1);
			}
			else
			{
				// try to parameterize gg->mumu mass shape
				fTemp_Rap[i][irap] = new TF1(Form("f%sTemp_RapBin%d", specName[i].Data(),irap), "pol4", 0, 5);
			}
			fTemp_Rap[i][irap]->SetTitle(Form("%1.1f < y < %1.1f", mRapLow[irap], mRapHi[irap]));
			fTemp_Rap[i][irap]->SetNpx(1000);

			c2->cd(1);
			gPad->SetLogy(1);
			setHisto(hMass_Rap[i][irap], 20, 0.6, 1, 1, 2);
			hMass_Rap[i][irap]->Scale(1./hMass_Rap[i][irap]->GetEntries());
			if(specName[i].Contains("Jpsi"))
			{
				hMass_Rap[i][irap]->GetXaxis()->SetRangeUser(2.5, 4);
			}
			else if(specName[i].Contains("Psi2S"))
			{
				if(specName[i].Contains("Feeddown")) hMass_Rap[i][irap]->GetXaxis()->SetRangeUser(2.5, 4.0);
				else                                 hMass_Rap[i][irap]->GetXaxis()->SetRangeUser(3.0, 4.5);
			}
			
			//hMass_Rap[i][irap]->GetYaxis()->SetRangeUser(0.5, hMass_Rap[i][irap]->GetMaximum()*2);
			hMass_Rap[i][irap]->GetXaxis()->SetTitleOffset(0.9);
			hMass_Rap[i][irap]->GetYaxis()->SetTitle("a.u.");
			if(i<nSpecs-1)
			{
				hMass_Rap[i][irap]->Fit(fTemp_Rap[i][irap], "RQ0", "", mMass[i]-0.4, mMass[i]+0.4);
				hMass_Rap[i][irap]->Draw("p");
				fTemp_Rap[i][irap]->Draw("same");
			}
			else
			{
				gPad->SetLogy(0);
				//hMass_Rap[i][irap]->Fit(fTemp_Rap[i][irap], "RQ0", "", mFitMassLow[irap], 5);
				hMass_Rap[i][irap]->Draw("p");
				//fTemp_Rap[i][irap]->Draw("same");
			}
			drawLatex(0.39, 0.95, Form("%1.1f < y < %1.1f", mRapLow[irap], mRapHi[irap]), mFont, 0.06, 1);

			c2->cd(2);
			gPad->SetLogy(1);
			hPt_Rap[i][irap]->Scale(1./hPt_Rap[i][irap]->GetEntries());
			setHisto(hPt_Rap[i][irap], 20, 0.6, 1, 1, 2);
			if(specName[i].Contains("InCoh"))
			{
				hPt_Rap[i][irap]->GetXaxis()->SetRangeUser(0, 2);
			}
			else if(specName[i].Contains("Feeddown"))
			{
				hPt_Rap[i][irap]->GetXaxis()->SetRangeUser(0, 1);
			}
			else
			{
				hPt_Rap[i][irap]->GetXaxis()->SetRangeUser(0, 0.5);
			}
			//hPt_Rap[i][irap]->GetYaxis()->SetRangeUser(0.5, hPt_Rap[i][irap]->GetMaximum()*2);
			hPt_Rap[i][irap]->GetXaxis()->SetTitleOffset(0.9);
			hPt_Rap[i][irap]->GetYaxis()->SetTitle("a.u.");
			hPt_Rap[i][irap]->Draw("p");
			drawLatex(0.39, 0.95, Form("%1.1f < y < %1.1f", mRapLow[irap], mRapHi[irap]), mFont, 0.06, 1);

			c2->SaveAs(Form("%s/%sTemp_RapBin%d.pdf", tempDir.Data(), specName[i].Data(), irap));
			c2->SaveAs(Form("%s/%sTemp_RapBin%d.png", tempDir.Data(), specName[i].Data(), irap));
		}//iy
	}//ispec

	TFile *fOutTemp = new TFile(Form("%s/MassPtTemp_AllSpecs.root", dir.Data()), "recreate");
	fOutTemp->cd();
	for(Int_t i=0; i<nSpecs; i++)
	{
		hMass[i]->Write();
		
		if(i<nSpecs-1) fTemp[i]->Write();

		hPt[i]->GetXaxis()->UnZoom();
		hPt[i]->Write();
		
		for(Int_t irap=0; irap<nRapBins; irap++)
		{
			hMass_Rap[i][irap]->Write();
			if(i<nSpecs-1) fTemp_Rap[i][irap]->Write();
			hPt_Rap[i][irap]->GetXaxis()->UnZoom();
			hPt_Rap[i][irap]->Write();
		}
	}
	fOutTemp->Close();

	cout << "End of program !" << endl;
}
//------------------------------------------------------ //bremsstrahlung tail
Double_t DoubleNGaus(Double_t *x, Double_t *par)
{
	Double_t N1       = par[0];
	Double_t mu       = par[1];
	Double_t s1       = par[2];
	Double_t N2       = par[3];
	Double_t s2       = par[4]*s1;
	Double_t binWidth = par[5];

	Double_t norm1    = (x[0]-mu)/s1;
	Double_t norm2    = (x[0]-mu)/s2;

	return ( N1/TMath::Sqrt(2*PI)/s1*TMath::Exp(-0.5*norm1*norm1) + N2/TMath::Sqrt(2*PI)/s2*TMath::Exp(-0.5*norm2*norm2) ) * binWidth;
}
//------------------------------------------------------ //bremsstrahlung tail
Double_t CrystalBall(Double_t *x, Double_t *par)
{
	Double_t N     = par[0];
	Double_t mu    = par[1];
	Double_t s     = par[2];
	Double_t alpha = par[3];
	Double_t n     = par[4];

	Double_t A     = TMath::Power(n/fabs(alpha), n) * TMath::Exp(-alpha*alpha/2.);
	Double_t B     = n/fabs(alpha) - fabs(alpha);
	Double_t norm  = (x[0]-mu)/s;

	if(norm > -alpha)
	{
		return N * TMath::Exp(-0.5*norm*norm);
	}
	else
	{
		return N * A * TMath::Power(B-norm, -n);
	}
}
//------------------------------------------------------ //bremsstrahlung tail
Double_t CrystalBallAndGaus(Double_t *x, Double_t *par)
{
	Double_t N     = par[0];
	Double_t mu    = par[1];
	Double_t s     = par[2];
	Double_t alpha = par[3];
	Double_t n     = par[4];

	Double_t A     = TMath::Power(n/fabs(alpha), n) * TMath::Exp(-alpha*alpha/2.);
	Double_t B     = n/fabs(alpha) - fabs(alpha);
	Double_t norm  = (x[0]-mu)/s;

	Double_t gausNorm = (x[0]-mu)/par[6];
	Double_t gausCom  = par[5] * TMath::Exp(-0.5*gausNorm*gausNorm);

	if(norm > -alpha)
	{
		return N * TMath::Exp(-0.5*norm*norm) + gausCom;
	}
	else
	{
		return N * A * TMath::Power(B-norm, -n) + gausCom;
	}
}
