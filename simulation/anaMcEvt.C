#include "../common/headers.h"
#include "../common/funUtil.h"
#include "../common/VertexCompositeTree.h"

TH1D *hnEvts;
TH3D *hVzvsVyvsVx;
TH2D *hHFMinusvsHFPlus;
TH3D *hNCanvsNtrkvsCen;
TH2D *hNtrkofflinevsNtrkHP;

TH3D *hPhivsEtavsPt_Gen;
TH3D *hForward_PhivsEtavsPt_Gen;
TH2D *hNegPtvsPosPt_Gen;
TH2D *hNegEtavsPosEta_Gen;
TH2D *hNegPhivsPosPhi_Gen;
TH1D *hDeltaR;
TH1D *hDeltaPt;
TH2D *hPtResvsGenPt;
TH2D *hEtaResvsGenEta;
TH2D *hPhiResvsGenPhi;

TH2D *hNRcMuvsNGenMu;
TH3D *hPhivsEtavsPt;

TH2D *hRcPairPtvsGenPairPt;
TH2D *hRcPairEtavsGenPairEta;
TH2D *hRcPairPhivsGenPairPhi;
TH2D *hMassResvsGenMass;
TH2D *hPairPtResvsGenPairPt;
TH2D *hRapResvsGenRap;

TH3D *hMvsPtvsRap_Gen;
TH3D *hMvsAsyPhivsRap_Gen;
TH2D *hDeltaPhivsM_Gen;

TH3D *hMvsPtvsRap_woEvtSel_woSmear;
TH3D *hMvsPtvsRap_woSmear;

TH3D *hMvsPtvsRap_woEvtSel;
TH3D *hMvsPtvsRap;

TH3D *hMvsAsyPhivsRap;
TH2D *hDeltaPhivsM;

// To calculate 3D efficiency
TH3D *hPosMuPhivsEtavsPt_Gen;
TH3D *hMthPosMuPhivsEtavsPt_Gen;
TH3D *hMthPosMuPhivsEtavsPt;
TH3D *hTrigPosMuPhivsEtavsPt;

TH3D *hNegMuPhivsEtavsPt_Gen;
TH3D *hMthNegMuPhivsEtavsPt_Gen;
TH3D *hMthNegMuPhivsEtavsPt;
TH3D *hTrigNegMuPhivsEtavsPt;

TH3D *hMthPosMuPhivsEtavsPtInpair;
TH3D *hTrigPosMuPhivsEtavsPtInpair;
TH3D *hMthNegMuPhivsEtavsPtInpair;
TH3D *hTrigNegMuPhivsEtavsPtInpair;

Bool_t   goodMuPair(const TLorentzVector posFourMom, const Bool_t isPosTrig, const TLorentzVector negFourMom, const Bool_t isNegTrig);
Bool_t   goodMuPair(VertexCompositeTree& evtTree, const int icand);
Double_t shiftDeltaPhi(Double_t dPhi);

void bookHistos();
void writeHistos(TString fileName = "test");

void anaMcEvt(TString fileName = "CohJpsi")
{
	cout<<"running with: "<<fileName<<endl;


	TH1::SetDefaultSumw2(kTRUE);

	std::string inputFile;
	if(fileName.EqualTo("LowMassGammaGamma"))     inputFile = "../rootfiles/VertexCompositeTree_STARLIGHT_LowMassGammaGammaToMuMu_GenFilter_DiMuMC_20200906.root";
	else if(fileName.EqualTo("CohJpsi"))          inputFile = "../rootfiles/VertexCompositeTree_STARLIGHT_CohJpsiToMuMu_GenFilter_DiMuMC_20200906.root";
	else if(fileName.EqualTo("CohJpsi_0n0n"))     inputFile = "../rootfiles/VertexCompositeTree_STARLIGHT_CohJpsiToMuMu_0n0n_GenFilter_DiMuMC_20210131.root";
	else if(fileName.EqualTo("CohJpsi_0nXn"))     inputFile = "../rootfiles/VertexCompositeTree_STARLIGHT_CohJpsiToMuMu_0nXn_GenFilter_DiMuMC_20210131.root";
	else if(fileName.EqualTo("CohJpsi_XnXn"))     inputFile = "../rootfiles/VertexCompositeTree_STARLIGHT_CohJpsiToMuMu_XnXn_GenFilter_DiMuMC_20210131.root";
	else if(fileName.EqualTo("InCohJpsi"))        inputFile = "../rootfiles/VertexCompositeTree_STARLIGHT_InCohJpsiToMuMu_GenFilter_DiMuMC_20200906.root";
	else if(fileName.EqualTo("CohPsi2S"))         inputFile = "../rootfiles/VertexCompositeTree_STARLIGHT_CohPsi2SToMuMu_GenFilter_DiMuMC_20200906.root";
	else if(fileName.EqualTo("CohPsi2SFeeddown")) inputFile = "../rootfiles/VertexCompositeTree_STARLIGHT_CohPsi2SFeeddownToMuMu_GenFilter_DiMuMC_20201002.root";
	else if(fileName.EqualTo("InCohPsi2S"))       inputFile = "../rootfiles/VertexCompositeTree_STARLIGHT_InCohPsi2SToMuMu_GenFilter_DiMuMC_20200906.root";
	else
	{
		cout<<"The filename string is wrong! Should be 'CohJpsi', 'CohJpsi_0n0n', 'CohJpsi_0nXn', 'CohJpsi_XnXn', 'InCohJpsi', 'CohPsi2S', 'CohPsi2SFeeddown', 'InCohPsi2S', OR 'LowMassGammaGamma' "<<endl;
		return;
	}

	const auto& csTreeDir = "dimucontana_mc";

	// Extract the tree
	VertexCompositeTree csTree;
	if (!csTree.GetTree(inputFile, csTreeDir)) 
	{
		cout << "Invalid Correct-Sign tree!" << endl;
		return;
	}

	if(!init()) 
	{
		cout<<"Initialization failed !"<<endl;
		return;
	}

	bookHistos();


	for (Long64_t jentry = 1; jentry < csTree.GetEntries(); jentry++) 
	{
		if (jentry % (csTree.GetEntries() / 10) == 0)
			cout << "begin " << jentry << "th entry...." << endl;

		// Get the entry
		if (csTree.GetEntry(jentry) < 0) 
		{
			cout << "Invalid correct-sign entry!" << endl;
			return;
		}

		Int_t   cen          = csTree.centrality();
		Float_t vtxX         = csTree.bestvtxX();
		Float_t vtxY         = csTree.bestvtxY();
		Float_t vtxZ         = csTree.bestvtxZ();
		Float_t hfsumETPlus  = csTree.HFsumETPlus();
		Float_t hfsumETMinus = csTree.HFsumETMinus();
		Int_t   nTrkoffline  = csTree.Ntrkoffline();
		Int_t   nTrkHP       = csTree.NtrkHP();

		hnEvts               ->Fill(0.5);
		hVzvsVyvsVx          ->Fill(vtxX, vtxY, vtxZ);
		hHFMinusvsHFPlus     ->Fill(hfsumETPlus, hfsumETMinus);
		hNCanvsNtrkvsCen     ->Fill(cen, nTrkoffline, csTree.candSize());
		hNtrkofflinevsNtrkHP ->Fill(nTrkHP, nTrkoffline);

		// Select Event - require this event has a valid vertex and is not beam-halo event
		//if(!csTree.evtSel()[2] || !csTree.evtSel()[3]) continue;
		Bool_t goodVtx = (csTree.evtSel()[2] && csTree.evtSel()[3]);
		if(goodVtx) hnEvts->Fill(1.5);

		//evtSel[4-15]
		//[4]=0: HFPlusMaxTower < 3 GeV;  [4]=1: HFPlusMaxTower > 3 GeV
		//[5]=0: HFMinusMaxTower < 3 GeV;  [5]=1: HFMinusMaxTower > 3 GeV
		//[6] is for Plus & [7] is for Minus; Threshold = 4 GeV
		//[8] is for Plus & [9] is for Minus; Threshold = 5 GeV
		//[10] is for Plus & [11] is for Minus; Threshold = 6 GeV
		//[12] is for Plus & [13] is for Minus; Threshold = 7 GeV
		//[14] is for Plus & [15] is for Minus; Threshold = 8 GeV
		//[16] is for Plus (Th = 7.6 GeV) & [17] is for Minus (Th = 7.3 GeV); 

		//if(csTree.evtSel()[16] || csTree.evtSel()[17]) continue;
		Bool_t goodHFVeto = (!csTree.evtSel()[16] && !csTree.evtSel()[17]);
		if(goodVtx & goodHFVeto) hnEvts->Fill(2.5);

		//if(nTrkHP != 2) continue;

		Bool_t passEvtSel = goodVtx && goodHFVeto && (nTrkHP==2);

		if(passEvtSel) hnEvts->Fill(3.5);

		// Loop over the correct-sign candidates
		Int_t nSoftMuon = 0;
		for (UInt_t icand = 0; icand < csTree.candSize_gen(); icand++) 
		{
			Double_t posPt_gen  = csTree.chargeD1_gen()[icand] > 0 ? csTree.pTD1_gen()[icand] : csTree.pTD2_gen()[icand]; 
			Double_t posEta_gen = csTree.chargeD1_gen()[icand] > 0 ? csTree.EtaD1_gen()[icand] : csTree.EtaD2_gen()[icand]; 
			Double_t posPhi_gen = csTree.chargeD1_gen()[icand] > 0 ? csTree.PhiD1_gen()[icand] : csTree.PhiD2_gen()[icand]; 
			Double_t negPt_gen  = csTree.chargeD1_gen()[icand] < 0 ? csTree.pTD1_gen()[icand] : csTree.pTD2_gen()[icand]; 
			Double_t negEta_gen = csTree.chargeD1_gen()[icand] < 0 ? csTree.EtaD1_gen()[icand] : csTree.EtaD2_gen()[icand]; 
			Double_t negPhi_gen = csTree.chargeD1_gen()[icand] < 0 ? csTree.PhiD1_gen()[icand] : csTree.PhiD2_gen()[icand]; 

			hPhivsEtavsPt_Gen      ->Fill(posPt_gen, posEta_gen, posPhi_gen);
			hPhivsEtavsPt_Gen      ->Fill(negPt_gen, negEta_gen, negPhi_gen);

			hPosMuPhivsEtavsPt_Gen ->Fill(posPt_gen, posEta_gen, posPhi_gen);
			hNegMuPhivsEtavsPt_Gen ->Fill(negPt_gen, negEta_gen, negPhi_gen);

			hNegPtvsPosPt_Gen      ->Fill(posPt_gen, negPt_gen);
			hNegEtavsPosEta_Gen    ->Fill(posEta_gen, negEta_gen);
			hNegPhivsPosPhi_Gen    ->Fill(posPhi_gen, negPhi_gen);

			TLorentzVector posFourMom_gen, negFourMom_gen, pairFourMom_gen;
			posFourMom_gen.SetPtEtaPhiM(posPt_gen, posEta_gen, posPhi_gen, Mmuon);
			negFourMom_gen.SetPtEtaPhiM(negPt_gen, negEta_gen, negPhi_gen, Mmuon);
			pairFourMom_gen = posFourMom_gen + negFourMom_gen;

			Double_t pt_gen   = pairFourMom_gen.Pt();
			Double_t eta_gen  = pairFourMom_gen.Eta();
			Double_t phi_gen  = pairFourMom_gen.Phi();
			Double_t mass_gen = pairFourMom_gen.M();
			Double_t y_gen    = pairFourMom_gen.Rapidity();

			if(y_gen>1.5)
			{
				hForward_PhivsEtavsPt_Gen->Fill(posPt_gen, posEta_gen, posPhi_gen);
				hForward_PhivsEtavsPt_Gen->Fill(negPt_gen, negEta_gen, negPhi_gen);
			}

			Double_t asyPhi_gen = 1 - TMath::Abs(shiftDeltaPhi(posFourMom_gen.DeltaPhi(negFourMom_gen))) / PI; //acoplanarity

			TVector3 muMomDiff_gen = posFourMom_gen.Vect() - negFourMom_gen.Vect();
			TVector3 pairMom_gen   = pairFourMom_gen.Vect();
			Double_t phiDiff_gen   = shiftDeltaPhi(pairMom_gen.DeltaPhi(muMomDiff_gen));

			hMvsPtvsRap_Gen        ->Fill(y_gen, pt_gen, mass_gen);
			hMvsAsyPhivsRap_Gen    ->Fill(y_gen, asyPhi_gen, mass_gen);
			hDeltaPhivsM_Gen       ->Fill(mass_gen, phiDiff_gen);

			Double_t posMthDeltaR = 99999999.;
			Double_t negMthDeltaR = 99999999.;
			Int_t    posRecoIdx = -1;
			Int_t    negRecoIdx = -1;
			
			for(UInt_t imu=0; imu<csTree.candSize_mu(); imu++)
			{
				if(!csTree.softMuon_mu()[imu]) continue;

				if(icand==0) nSoftMuon++;

				Double_t muPt  = csTree.pT_mu()[imu];
				Double_t muEta = csTree.eta_mu()[imu];
				Double_t muPhi = csTree.phi_mu()[imu];

				TVector3 recoMom; recoMom.SetPtEtaPhi(muPt, muEta, muPhi);
				Double_t posDeltaR = posFourMom_gen.Vect().DeltaR(recoMom);
				Double_t negDeltaR = negFourMom_gen.Vect().DeltaR(recoMom);

				// DeltaPt/genPt < 0.15 ensures to keep 99.9842% events
				// DeltaR < 0.05 ensures to keep 99.9763% events
				if(
						TMath::Abs((muPt-posPt_gen)/posPt_gen) < 0.15
						&& posDeltaR < 0.05 
						&& posDeltaR < posMthDeltaR
				  )
				{
					posMthDeltaR = posDeltaR;
					posRecoIdx = imu;
				}

				if(
						TMath::Abs((muPt-negPt_gen)/negPt_gen) < 0.15
						&& negDeltaR < 0.05
						&& negDeltaR < negMthDeltaR
				  )
				{
					negMthDeltaR = negDeltaR;
					negRecoIdx = imu;
				}
			}

			if(posRecoIdx>=0 && negRecoIdx>=0 && posRecoIdx == negRecoIdx){
				cout<<"One reco-track is matched to multiple gen-tracks !"<<endl;
			}

			if(posRecoIdx>=0){
				Double_t muPt  = csTree.pT_mu()[posRecoIdx];
				Double_t muEta = csTree.eta_mu()[posRecoIdx];
				Double_t muPhi = csTree.phi_mu()[posRecoIdx];
				Bool_t   isTrigMu = csTree.trigMuon_mu()[trigIdx][posRecoIdx];

				hDeltaR->Fill(posMthDeltaR);
				hDeltaPt->Fill((muPt-posPt_gen)/posPt_gen);

				hPtResvsGenPt->Fill(posPt_gen, (muPt-posPt_gen)/posPt_gen);
				hEtaResvsGenEta->Fill(posEta_gen, muEta-posEta_gen);
				hPhiResvsGenPhi->Fill(posPhi_gen, muPhi-posPhi_gen);

				hMthPosMuPhivsEtavsPt_Gen->Fill(posPt_gen, posEta_gen, posPhi_gen);
				hMthPosMuPhivsEtavsPt->Fill(muPt, muEta, muPhi);

				if(isTrigMu){
					hTrigPosMuPhivsEtavsPt->Fill(muPt, muEta, muPhi);
				}
			}

			if(negRecoIdx>=0){
				Double_t muPt  = csTree.pT_mu()[negRecoIdx];
				Double_t muEta = csTree.eta_mu()[negRecoIdx];
				Double_t muPhi = csTree.phi_mu()[negRecoIdx];
				Bool_t   isTrigMu = csTree.trigMuon_mu()[trigIdx][negRecoIdx];

				hDeltaR->Fill(negMthDeltaR);
				hDeltaPt->Fill((muPt-negPt_gen)/negPt_gen);

				hPtResvsGenPt->Fill(negPt_gen, (muPt-negPt_gen)/negPt_gen);
				hEtaResvsGenEta->Fill(negEta_gen, muEta-negEta_gen);
				hPhiResvsGenPhi->Fill(negPhi_gen, muPhi-negPhi_gen);

				hMthNegMuPhivsEtavsPt_Gen->Fill(negPt_gen, negEta_gen, negPhi_gen);
				hMthNegMuPhivsEtavsPt->Fill(muPt, muEta, muPhi);

				if(isTrigMu){
					hTrigNegMuPhivsEtavsPt->Fill(muPt, muEta, muPhi);
				}
			}

			Int_t recoIdx = csTree.RecIdx_gen()[icand];

			if(recoIdx<0 || csTree.candSize()<=0) continue;

			Double_t posPt     = csTree.chargeD1()[recoIdx] > 0 ? csTree.pTD1()[recoIdx] : csTree.pTD2()[recoIdx]; 
			Double_t posEta    = csTree.chargeD1()[recoIdx] > 0 ? csTree.EtaD1()[recoIdx] : csTree.EtaD2()[recoIdx]; 
			Double_t posPhi    = csTree.chargeD1()[recoIdx] > 0 ? csTree.PhiD1()[recoIdx] : csTree.PhiD2()[recoIdx]; 
			Double_t isPosTrig = csTree.chargeD1()[recoIdx] > 0 ? csTree.trigMuon1()[trigIdx][recoIdx] : csTree.trigMuon2()[trigIdx][recoIdx]; 
			Double_t negPt     = csTree.chargeD1()[recoIdx] < 0 ? csTree.pTD1()[recoIdx] : csTree.pTD2()[recoIdx]; 
			Double_t negEta    = csTree.chargeD1()[recoIdx] < 0 ? csTree.EtaD1()[recoIdx] : csTree.EtaD2()[recoIdx]; 
			Double_t negPhi    = csTree.chargeD1()[recoIdx] < 0 ? csTree.PhiD1()[recoIdx] : csTree.PhiD2()[recoIdx]; 
			Double_t isNegTrig = csTree.chargeD1()[recoIdx] < 0 ? csTree.trigMuon1()[trigIdx][recoIdx] : csTree.trigMuon2()[trigIdx][recoIdx]; 

			if(!csTree.softCand(recoIdx)) continue;

			// shift mean
			Double_t shiftPosPt = posPt - funPtMeanShift->Eval(posPt_gen)*posPt_gen;
			Double_t shiftNegPt = negPt - funPtMeanShift->Eval(negPt_gen)*negPt_gen;

			// smear width
			Double_t posPtNew = posPt_gen + (shiftPosPt - posPt_gen)/funRawPtRes->Eval(posPt_gen)*funTunedPtRes->Eval(posPt_gen);
			Double_t negPtNew = negPt_gen + (shiftNegPt - negPt_gen)/funRawPtRes->Eval(negPt_gen)*funTunedPtRes->Eval(negPt_gen);

			hMthPosMuPhivsEtavsPtInpair->Fill(posPtNew, posEta, posPhi);
			if(isPosTrig) hTrigPosMuPhivsEtavsPtInpair->Fill(posPtNew, posEta, posPhi);

			hMthNegMuPhivsEtavsPtInpair->Fill(negPtNew, negEta, negPhi);
			if(isNegTrig) hTrigNegMuPhivsEtavsPtInpair->Fill(negPtNew, negEta, negPhi);

			TLorentzVector posFourMom, negFourMom, pairFourMom;
			posFourMom.SetPtEtaPhiM(posPtNew, posEta, posPhi, Mmuon);
			negFourMom.SetPtEtaPhiM(negPtNew, negEta, negPhi, Mmuon);
			pairFourMom = posFourMom + negFourMom;
			//cout<<"posPhi: "<<posPhi<<"   negPhi: "<<negPhi<<"    phiAngel:"<<posFourMom.DeltaPhi(negFourMom)<<endl;

			// without pt smear
			Double_t pt_woSmear   = csTree.pT()[recoIdx];
			Double_t eta_woSmear  = csTree.eta()[recoIdx];
			Double_t phi_woSmear  = csTree.phi()[recoIdx];
			Double_t mass_woSmear = csTree.mass()[recoIdx];
			Double_t y_woSmear    = csTree.y()[recoIdx];

			// with pt smear
			Double_t pt   = pairFourMom.Pt();
			Double_t eta  = pairFourMom.Eta();
			Double_t phi  = pairFourMom.Phi();
			Double_t mass = pairFourMom.M();
			Double_t y    = pairFourMom.Rapidity();
			//cout<<"pt: "<<pt<<"   eta: "<<eta<<"   phi: "<<phi<<"   M: "<<mass<<"   y: "<<y<<endl;

			hRcPairPtvsGenPairPt     ->Fill(pt_gen,  pt );
			hRcPairEtavsGenPairEta   ->Fill(phi_gen, phi);
			hRcPairPhivsGenPairPhi   ->Fill(eta_gen, eta);
			hMassResvsGenMass        ->Fill(mass_gen, (mass-mass_gen)/mass_gen);
			hPairPtResvsGenPairPt    ->Fill(pt_gen,  pt - pt_gen);
			hRapResvsGenRap          ->Fill(y_gen,   y-y_gen    );

			Double_t asyPhi = 1 - TMath::Abs(shiftDeltaPhi(posFourMom.DeltaPhi(negFourMom))) / PI; //acoplanarity

			// already required soft muon pair before
			
			//----------------------------------------------------------------------------------------
			// without pt smear
			if(goodMuPair(csTree, recoIdx))
			{
				hMvsPtvsRap_woEvtSel_woSmear->Fill(y_woSmear, pt_woSmear, mass_woSmear);

				if( passEvtSel )
				{
					hMvsPtvsRap_woSmear->Fill(y_woSmear, pt_woSmear, mass_woSmear);
				}
			}
			//----------------------------------------------------------------------------------------

			//----------------------------------------------------------------------------------------
			// with pt smear
			if(!goodMuPair(posFourMom, isPosTrig, negFourMom, isNegTrig)) continue;
			hMvsPtvsRap_woEvtSel->Fill(y, pt, mass);
			//----------------------------------------------------------------------------------------

			if(!passEvtSel) continue;
			
			hMvsPtvsRap      ->Fill(y, pt,     mass);
			hMvsAsyPhivsRap  ->Fill(y, asyPhi, mass);

			TVector3 muMomDiff = posFourMom.Vect() - negFourMom.Vect();
			TVector3 pairMom   = pairFourMom.Vect();

			Double_t phiDiff   = shiftDeltaPhi(pairMom.DeltaPhi(muMomDiff));
			hDeltaPhivsM->Fill(mass, phiDiff);
		}

		hNRcMuvsNGenMu->Fill( csTree.candSize_gen()*2, nSoftMuon);
	}

	system("mkdir -p mcHistos");
	writeHistos(Form("mcHistos/dimuonHistos.%s", fileName.Data()));
}

void bookHistos()
{
	const Int_t    mHistRapBins     = 60;
	const Double_t mHistRapLow      = -3;
	const Double_t mHistRapHi       = 3;

	const Int_t    mHistPtBins      = 800;
	const Double_t mHistPtLow       = 0;
	const Double_t mHistPtHi        = 4;
	
	const Int_t    mHistMassBins    = 300;
	const Double_t mHistMassLow     = 2;
	const Double_t mHistMassHi      = 5;
	
	const Int_t    mHistAsyPhiBins  = 300;
	const Double_t mHistAsyPhiLow   = 0;
	const Double_t mHistAsyPhiHi    = 0.3;
	
	const Int_t    mHistAsyPtBins   = 500;
	const Double_t mHistAsyPtLow    = 0;
	const Double_t mHistAsyPtHi     = 0.5;

	// event level
	hnEvts = new TH1D("hnEvts", "hnEvts;", 5, 0, 5);
	hnEvts->GetXaxis()->SetBinLabel(1, "trigEvt");
	hnEvts->GetXaxis()->SetBinLabel(2, "validVtx & !Beam-halo");
	hnEvts->GetXaxis()->SetBinLabel(3, "HFMaxE <= 7.6(7.3) GeV");
	hnEvts->GetXaxis()->SetBinLabel(4, "N_{trk}^{HP} == 2");
	hnEvts->GetXaxis()->LabelsOption("d");
	hnEvts->GetXaxis()->SetLabelSize(0.06);

	hVzvsVyvsVx = new TH3D("hVzvsVyvsVx", "hVzvsVyvsVx; V_{x} (cm); V_{y} (cm); V_{z} (cm)", 100, -1, 1, 100, -1, 1, 300, -30, 30);
	hHFMinusvsHFPlus = new TH2D("hHFMinusvsHFPlus", "hHFMinusvsHFPlus; HFsumETPlus; HFsumETMinus", 200, 0, 20, 200, 0, 20);
	hNCanvsNtrkvsCen = new TH3D("hNCanvsNtrkvsCen", "hNCanvsNtrkvsCen; Centrality; N_{trk}^{offline}; # CS pair", 50, 150, 200, 10, 0, 10, 5, 0, 5);
	hNtrkofflinevsNtrkHP = new TH2D("hNtrkofflinevsNtrkHP", "hNtrkofflinevsNtrkHP; N_{trk}^{HP}; N_{trk}^{offline}", 10, 0, 10, 10, 0, 10);

	hPhivsEtavsPt_Gen   = new TH3D("hPhivsEtavsPt_Gen", "hPhivsEtavsPt_Gen; p_{T} (GeV/c); #eta; #phi", 500, 0, 5, 300, -3, 3, 180, -PI, PI);
	hForward_PhivsEtavsPt_Gen   = new TH3D("hForward_PhivsEtavsPt_Gen", "hForward_PhivsEtavsPt_Gen; p_{T} (GeV/c); #eta; #phi", 500, 0, 5, 300, -3, 3, 180, -PI, PI);
	hNegPtvsPosPt_Gen    = new TH2D("hNegPtvsPosPt_Gen", "hNegPtvsPosPt_Gen; #mu^{+} p_{T} (GeV/c); #mu^{-} p_{T} (GeV/c);", 500, 0, 5, 500, 0, 5);
	hNegEtavsPosEta_Gen  = new TH2D("hNegEtavsPosEta_Gen", "hNegEtavsPosEta_Gen; #mu^{+} #eta; #mu^{-} #eta;", 300, -3, 3, 300, -3, 3);
	hNegPhivsPosPhi_Gen  = new TH2D("hNegPhivsPosPhi_Gen", "hNegPhivsPosPhi_Gen; #mu^{+} #phi; #mu^{-} #phi;", 180, -PI, PI, 180, -PI, PI);
	hDeltaR              = new TH1D("hDeltaR", "hDeltaR; #DeltaR", 2000, 0, 0.2);
	hDeltaPt             = new TH1D("hDeltaPt", "hDeltaPt; (p_{T}^{Rc}-p_{T}^{Gen})/p_{T}^{Gen}", 5000, -0.5, 0.5);
	hPtResvsGenPt        = new TH2D("hPtResvsGenPt", "hPtResvsGenPt; p_{T}^{Gen} (GeV/c); (p_{T}^{Rc}-p_{T}^{Gen})/p_{T}^{Gen};", 500, 0, 5, 300, -0.15, 0.15);
	hEtaResvsGenEta      = new TH2D("hEtaResvsGenEta", "hEtaResvsGenEta; #eta^{Gen}; #eta^{Rc}-#eta^{Gen};", 300, -3, 3, 200, -0.01, 0.01);
	hPhiResvsGenPhi      = new TH2D("hPhiResvsGenPhi", "hPhiResvsGenPhi; #phi^{Gen}; #phi^{Rc}-#phi^{Gen};", 180, -PI, PI, 200, -0.01, 0.01);

	hNRcMuvsNGenMu          = new TH2D("hNRcMuvsNGenMu", "hNRcMuvsNGenMu; # of muon (GEN); # of muon (Soft)", 10, 0, 10, 10, 0, 10);
	hPhivsEtavsPt           = new TH3D("hPhivsEtavsPt", "hPhivsEtavsPt; p_{T} (GeV/c); #eta; #phi", 500, 0, 5, 300, -3, 3, 180, -PI, PI);

	hRcPairPtvsGenPairPt   = new TH2D("hRcPairPtvsGenPairPt", "hRcPairPtvsGenPairPt; pair p_{T}^{Gen} (GeV/c); pair p_{T}^{Rc} (GeV/c);", mHistPtBins, mHistPtLow, mHistPtHi, mHistPtBins, mHistPtLow, mHistPtHi);
	hRcPairEtavsGenPairEta = new TH2D("hRcPairEtavsGenPairEta", "hRcPairEtavsGenPairEta; pair #eta^{Gen}; pair #eta^{Rc};", 160, -4, 4, 160, -4, 4);
	hRcPairPhivsGenPairPhi = new TH2D("hRcPairPhivsGenPairPhi", "hRcPairPhivsGenPairPhi; pair #phi^{Gen}; pair #phi^{Rc};", 120, -PI, PI, 120, -PI, PI);
	hMassResvsGenMass      = new TH2D("hMassResvsGenMass", "hMassResvsGenMass; M_{#mu#mu}^{Gen} (GeV/c^{2}); (M_{#mu#mu}^{Rc} - M_{#mu#mu}^{Gen})/M_{#mu#mu}^{Gen};", mHistMassBins, mHistMassLow, mHistMassHi, 300, -0.15, 0.15);
	hPairPtResvsGenPairPt  = new TH2D("hPairPtResvsGenPairPt", "hPairPtResvsGenPairPt; J/#psi p_{T}^{Gen} (GeV/c); J/#psi p_{T}^{Rc} - p_{T}^{Gen};", 300, 0, 0.3, 300, -0.3, 0.3);
	hRapResvsGenRap        = new TH2D("hRapResvsGenRap", "hRapResvsGenRap; pair y^{Gen}; y^{Rc} - y^{Gen};", 600, -3, 3, 1000, -0.1, 0.1);

	hMvsPtvsRap_Gen     = new TH3D("hMvsPtvsRap_Gen", "hMvsPtvsRap_Gen; Rapidity; p_{T} (GeV/c); M_{#mu#mu} (GeV/c^{2})", mHistRapBins, mHistRapLow, mHistRapHi, mHistPtBins, mHistPtLow, mHistPtHi, mHistMassBins, mHistMassLow, mHistMassHi);
	hMvsAsyPhivsRap_Gen = new TH3D("hMvsAsyPhivsRap_Gen", "hMvsAsyPhivsRap_Gen; Rapidity; #alpha; M_{#mu#mu} (GeV/c^{2})", mHistRapBins, mHistRapLow, mHistRapHi, mHistAsyPhiBins, mHistAsyPhiLow, mHistAsyPhiHi, mHistMassBins, mHistMassLow, mHistMassHi);
	hDeltaPhivsM_Gen   = new TH2D("hDeltaPhivsM_Gen", "hDeltaPhivsM_Gen; M_{#mu#mu} (GeV/c^{2}); #phi_{#mu^{+}+#mu^{-}} - #phi_{#mu^{+}-#mu^{-}}", mHistMassBins, mHistMassLow, mHistMassHi, 120, -PI, PI);

	hMvsPtvsRap_woEvtSel_woSmear = new TH3D("hMvsPtvsRap_woEvtSel_woSmear", "hMvsPtvsRap_woEvtSel_woSmear; Rapidity; p_{T} (GeV/c); M_{#mu#mu} (GeV/c^{2})", mHistRapBins, mHistRapLow, mHistRapHi, mHistPtBins, mHistPtLow, mHistPtHi, mHistMassBins, mHistMassLow, mHistMassHi);
	hMvsPtvsRap_woSmear = new TH3D("hMvsPtvsRap_woSmear", "hMvsPtvsRap_woSmear; Rapidity; p_{T} (GeV/c); M_{#mu#mu} (GeV/c^{2})", mHistRapBins, mHistRapLow, mHistRapHi, mHistPtBins, mHistPtLow, mHistPtHi, mHistMassBins, mHistMassLow, mHistMassHi);

	hMvsPtvsRap_woEvtSel = new TH3D("hMvsPtvsRap_woEvtSel", "hMvsPtvsRap_woEvtSel; Rapidity; p_{T} (GeV/c); M_{#mu#mu} (GeV/c^{2})", mHistRapBins, mHistRapLow, mHistRapHi, mHistPtBins, mHistPtLow, mHistPtHi, mHistMassBins, mHistMassLow, mHistMassHi);
	hMvsPtvsRap = new TH3D("hMvsPtvsRap", "hMvsPtvsRap; Rapidity; p_{T} (GeV/c); M_{#mu#mu} (GeV/c^{2})", mHistRapBins, mHistRapLow, mHistRapHi, mHistPtBins, mHistPtLow, mHistPtHi, mHistMassBins, mHistMassLow, mHistMassHi);
	hMvsAsyPhivsRap = new TH3D("hMvsAsyPhivsRap", "hMvsAsyPhivsRap; Rapidity; #alpha; M_{#mu#mu} (GeV/c^{2})", mHistRapBins, mHistRapLow, mHistRapHi, mHistAsyPhiBins, mHistAsyPhiLow, mHistAsyPhiHi, mHistMassBins, mHistMassLow, mHistMassHi);
	hDeltaPhivsM = new TH2D("hDeltaPhivsM", "hDeltaPhivsM; M_{#mu#mu} (GeV/c^{2}); #phi_{#mu^{+}+#mu^{-}} - #phi_{#mu^{+}-#mu^{-}}", mHistMassBins, mHistMassLow, mHistMassHi, 120, -PI, PI);

	const Int_t nPtBins = 25;
	Double_t    Pt[nPtBins+1] = {0, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 3.0, 3.2, 3.6, 4.0};

	const Int_t nEtaBins = 100;
	Double_t Eta[nEtaBins+1];
	Double_t mEtaLow = -2.5, mEtaHi = 2.5;
	Double_t mEtaStep = (mEtaHi - mEtaLow)/nEtaBins;
	
	for(Int_t ieta=0; ieta<=nEtaBins; ieta++)
	{
		Eta[ieta] = mEtaLow + ieta*mEtaStep;
	}

	const Int_t nPhiBins = 60;
	Double_t Phi[nPhiBins+1];
	Double_t mPhiLow = -PI, mPhiHi = PI;
	Double_t mPhiStep = (mPhiHi - mPhiLow)/nPhiBins;
	for(Int_t iphi=0; iphi<=nPhiBins; iphi++){
		Phi[iphi] = mPhiLow + iphi*mPhiStep;
	}

	cout<<endl;
	cout<<"Pt Bin Boundary: "<<endl;
	for(Int_t ibin=0; ibin<nPtBins; ibin++){
		cout<<Pt[ibin]<<", ";
	}
	cout<<Pt[nPtBins]<<endl;
	cout<<endl;

	cout<<"Eta Bin Boundary: "<<endl;
	for(Int_t ibin=0; ibin<nEtaBins; ibin++){
		cout<<Eta[ibin]<<", ";
	}
	cout<<Eta[nEtaBins]<<endl;
	cout<<endl;

	cout<<"Phi Bin Boundary: "<<endl;
	for(Int_t ibin=0; ibin<nPhiBins; ibin++){
		cout<<Phi[ibin]<<", ";
	}
	cout<<Phi[nPhiBins]<<endl;
	cout<<endl;

	hPosMuPhivsEtavsPt_Gen    = new TH3D("hPosMuPhivsEtavsPt_Gen", "hPosMuPhivsEtavsPt_Gen; p_{T} (GeV/c); #eta; #phi", nPtBins, Pt, nEtaBins, Eta, nPhiBins, Phi);
	hMthPosMuPhivsEtavsPt_Gen = new TH3D("hMthPosMuPhivsEtavsPt_Gen", "hMthPosMuPhivsEtavsPt_Gen; p_{T} (GeV/c); #eta; #phi", nPtBins, Pt, nEtaBins, Eta, nPhiBins, Phi);
	hMthPosMuPhivsEtavsPt     = new TH3D("hMthPosMuPhivsEtavsPt", "hMthPosMuPhivsEtavsPt; p_{T} (GeV/c); #eta; #phi", nPtBins, Pt, nEtaBins, Eta, nPhiBins, Phi);
	hTrigPosMuPhivsEtavsPt    = new TH3D("hTrigPosMuPhivsEtavsPt", "hTrigPosMuPhivsEtavsPt; p_{T} (GeV/c); #eta; #phi", nPtBins, Pt, nEtaBins, Eta, nPhiBins, Phi);
	hNegMuPhivsEtavsPt_Gen    = new TH3D("hNegMuPhivsEtavsPt_Gen", "hNegMuPhivsEtavsPt_Gen; p_{T} (GeV/c); #eta; #phi", nPtBins, Pt, nEtaBins, Eta, nPhiBins, Phi);
	hMthNegMuPhivsEtavsPt_Gen = new TH3D("hMthNegMuPhivsEtavsPt_Gen", "hMthNegMuPhivsEtavsPt_Gen; p_{T} (GeV/c); #eta; #phi", nPtBins, Pt, nEtaBins, Eta, nPhiBins, Phi);
	hMthNegMuPhivsEtavsPt     = new TH3D("hMthNegMuPhivsEtavsPt", "hMthNegMuPhivsEtavsPt; p_{T} (GeV/c); #eta; #phi", nPtBins, Pt, nEtaBins, Eta, nPhiBins, Phi);
	hTrigNegMuPhivsEtavsPt    = new TH3D("hTrigNegMuPhivsEtavsPt", "hTrigNegMuPhivsEtavsPt; p_{T} (GeV/c); #eta; #phi", nPtBins, Pt, nEtaBins, Eta, nPhiBins, Phi);

	hMthPosMuPhivsEtavsPtInpair  = new TH3D("hMthPosMuPhivsEtavsPtInpair", "hMthPosMuPhivsEtavsPtInpair; p_{T} (GeV/c); #eta; #phi", nPtBins, Pt, nEtaBins, Eta, nPhiBins, Phi);
	hTrigPosMuPhivsEtavsPtInpair = new TH3D("hTrigPosMuPhivsEtavsPtInpair", "hTrigPosMuPhivsEtavsPtInpair; p_{T} (GeV/c); #eta; #phi", nPtBins, Pt, nEtaBins, Eta, nPhiBins, Phi);
	hMthNegMuPhivsEtavsPtInpair  = new TH3D("hMthNegMuPhivsEtavsPtInpair", "hMthNegMuPhivsEtavsPtInpair; p_{T} (GeV/c); #eta; #phi", nPtBins, Pt, nEtaBins, Eta, nPhiBins, Phi);
	hTrigNegMuPhivsEtavsPtInpair = new TH3D("hTrigNegMuPhivsEtavsPtInpair", "hTrigNegMuPhivsEtavsPtInpair; p_{T} (GeV/c); #eta; #phi", nPtBins, Pt, nEtaBins, Eta, nPhiBins, Phi);
}

void writeHistos(TString fileName)
{
	TFile* fOut = new TFile(Form("%s.root", fileName.Data()), "recreate");

	fOut->cd();

	hnEvts->Write();
	hVzvsVyvsVx->Write();
	hHFMinusvsHFPlus->Write();
	hNCanvsNtrkvsCen->Write();
	hNtrkofflinevsNtrkHP->Write();

	hPhivsEtavsPt_Gen->Write();
	hForward_PhivsEtavsPt_Gen->Write();
	hNegPtvsPosPt_Gen->Write();
	hNegEtavsPosEta_Gen->Write();
	hNegPhivsPosPhi_Gen->Write();
	hDeltaR->Write();
	hDeltaPt->Write();
	hPtResvsGenPt->Write();
	hEtaResvsGenEta->Write();
	hPhiResvsGenPhi->Write();

	hNRcMuvsNGenMu->Write();
	hPhivsEtavsPt->Write();

	hRcPairPtvsGenPairPt->Write();
	hRcPairEtavsGenPairEta->Write();
	hRcPairPhivsGenPairPhi->Write();
	hMassResvsGenMass->Write();
	hPairPtResvsGenPairPt->Write();
	hRapResvsGenRap->Write();

	hMvsPtvsRap_Gen->Write();
	hMvsAsyPhivsRap_Gen->Write();
	hDeltaPhivsM_Gen->Write();

	hMvsPtvsRap_woEvtSel_woSmear->Write();
	hMvsPtvsRap_woSmear->Write();

	hMvsPtvsRap_woEvtSel->Write();
	hMvsPtvsRap->Write();
	hMvsAsyPhivsRap->Write();
	hDeltaPhivsM->Write();

	hPosMuPhivsEtavsPt_Gen->Write();
	hMthPosMuPhivsEtavsPt_Gen->Write();
	hMthPosMuPhivsEtavsPt->Write();
	hTrigPosMuPhivsEtavsPt->Write();

	hNegMuPhivsEtavsPt_Gen->Write();
	hMthNegMuPhivsEtavsPt_Gen->Write();
	hMthNegMuPhivsEtavsPt->Write();
	hTrigNegMuPhivsEtavsPt->Write();

	hMthPosMuPhivsEtavsPtInpair->Write();
	hTrigPosMuPhivsEtavsPtInpair->Write();
	hMthNegMuPhivsEtavsPtInpair->Write();
	hTrigNegMuPhivsEtavsPtInpair->Write();

	fOut->Close();
}

Bool_t goodMuPair(VertexCompositeTree& evtTree, const int icand)
{
	Double_t mTrkPtTh1  = fTrkAcc  ->Eval(evtTree.EtaD1()[icand]);
	Double_t mTrkPtTh2  = fTrkAcc  ->Eval(evtTree.EtaD2()[icand]);
	Double_t mTrigPtTh1 = fTrigAcc ->Eval(evtTree.EtaD1()[icand]);
	Double_t mTrigPtTh2 = fTrigAcc ->Eval(evtTree.EtaD2()[icand]);

	if(evtTree.pTD1()[icand] < mTrkPtTh1 || evtTree.pTD2()[icand] < mTrkPtTh2) return kFALSE;

	Bool_t isTrigAcc1 = kFALSE, isTrigAcc2 = kFALSE;
	if(evtTree.trigMuon1()[trigIdx][icand] && evtTree.pTD1()[icand] >= mTrigPtTh1) isTrigAcc1 = kTRUE;
	if(evtTree.trigMuon2()[trigIdx][icand] && evtTree.pTD2()[icand] >= mTrigPtTh2) isTrigAcc2 = kTRUE;

	if(!isTrigAcc1 && !isTrigAcc2) return kFALSE;

	return kTRUE;
}

Bool_t goodMuPair(const TLorentzVector posFourMom, const Bool_t isPosTrig, const TLorentzVector negFourMom, const Bool_t isNegTrig)
{
	Double_t mPosTrkPtTh  = fTrkAcc->Eval(posFourMom.Eta());
	Double_t mNegTrkPtTh  = fTrkAcc->Eval(negFourMom.Eta());
	Double_t mPosTrigPtTh = fTrigAcc->Eval(posFourMom.Eta());
	Double_t mNegTrigPtTh = fTrigAcc->Eval(negFourMom.Eta());

	if(posFourMom.Pt() < mPosTrkPtTh || negFourMom.Pt() < mNegTrkPtTh) return kFALSE;

	Bool_t isPosTrigAcc = kFALSE, isNegTrigAcc = kFALSE;
	if(isPosTrig && posFourMom.Pt() >= mPosTrigPtTh) isPosTrigAcc = kTRUE;
	if(isNegTrig && negFourMom.Pt() >= mNegTrigPtTh) isNegTrigAcc = kTRUE;

	if(!isPosTrigAcc && !isNegTrigAcc) return kFALSE;

	return kTRUE;
}

Double_t shiftDeltaPhi(Double_t dPhi)
{
	if (isnan(dPhi))
		return -999;

	while (dPhi < -PI)
		dPhi += 2 * PI;
	while (dPhi >= PI)
		dPhi -= 2 * PI;

	return dPhi;
}
