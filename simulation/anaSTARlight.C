#include "TParticle.h" 
#include "../common/headers.h" 

const Double_t mRapLow = 1.45, mRapHi = 2.45;
const Double_t mMuPtTh = 0.7, mMuEtaLow = 1.28, mMuEtaHi = 2.42;
const Double_t PI = TMath::Pi();

void anaSTARlight(TString parSpec = "CohJpsi", const Int_t nAnaEvts=1000000)
{
	string inFileName = "";
	if(parSpec.EqualTo("CohJpsi"))           inFileName = "/Users/syang/work/run2/STARlight/genSTARlightLHE_v313/CohJpsi/slight.CohJpsi.0_30M.out";
	else if(parSpec.EqualTo("CohJpsi_0n0n")) inFileName = "/Users/syang/work/run2/STARlight/genSTARlightLHE_v313/CohJpsi_0n0n/cohJpsi_0n0n_0_60M.out";
	else if(parSpec.EqualTo("CohJpsi_0nXn")) inFileName = "/Users/syang/work/run2/STARlight/genSTARlightLHE_v313/CohJpsi_0nXn/cohJpsi_0nXn_0_50M.out";
	else if(parSpec.EqualTo("CohJpsi_XnXn")) inFileName = "/Users/syang/work/run2/STARlight/genSTARlightLHE_v313/CohJpsi_XnXn/cohJpsi_XnXn_0_50M.out";
	else if(parSpec.EqualTo("InCohJpsi"))    inFileName = "/Users/syang/work/run2/STARlight/genSTARlightLHE_v313/InCohJpsi/slight.InCohJpsi.0_30M.out";
	else if(parSpec.EqualTo("CohPsi2S"))     inFileName = "/Users/syang/work/run2/STARlight/genSTARlightLHE_v313/CohPsi2S/slight.CohPsi2S.0_25M.out";
	else if(parSpec.EqualTo("InCohPsi2S"))   inFileName = "/Users/syang/work/run2/STARlight/genSTARlightLHE_v313/InCohPsi2S/slight.InCohPsi2S.0_25M.out";
	else 
	{
		cout<<"The particle species string is wrong! Should be 'CohJpsi', 'CohJpsi_0n0n', 'CohJpsi_0nXn', 'CohJpsi_XnXn', 'InCohJpsi', 'CohPsi2S' OR 'InCohPsi2S'"<<endl;
		return;
	}

	TH1D *hnEvts                = new TH1D("hnEvts", "hnEvts", 5, -0.5, 4.5);
	TH3D *hMvsPtvsRap_STARlight = new TH3D("hMvsPtvsRap_STARlight", "hMvsPtvsRap_STARlight; y; p_{T} (GeV); M_{#mu#mu} (GeV); Entries", 100, -2.5, 2.5, 2000, 0, 2, 300, 2, 5);
	TH3D *hMvsPtvsRap_muFilter  = new TH3D("hMvsPtvsRap_muFilter", "hMvsPtvsRap_muFilter; y; p_{T} (GeV); M_{#mu#mu} (GeV); Entries", 100, -2.5, 2.5, 2000, 0, 2, 300, 2, 5);
	TH2D *hNegPtvsPosPt         = new TH2D("hNegPtvsPosPt", "hNegPtvsPosPt; #mu^{+} p_{T} (GeV); #mu^{-} p_{T} (GeV)", 200, 0, 10, 200, 0, 10);
	TH2D *hNegEtavsPosEta       = new TH2D("hNegEtavsPosEta", "hNegEtavsPosEta; #mu^{+} #eta; #mu^{-} #eta", 200, -10, 10, 200, -10, 10);
	TH2D *hNegPhivsPosPhi       = new TH2D("hNegPhivsPosPhi", "hNegPhivsPosPhi; #mu^{+} #phi; #mu^{-} #phi", 120, -PI, PI, 120, -PI, PI);

	ifstream infile(inFileName.c_str());
	if (!infile.is_open()) { cout << "\t ERROR: I can not open \"" << inFileName << "\"" << endl; return; }

	string temp_string, temp;
	istringstream curstring;

	int nEvts=0; // event_counter
	int evtIdx = 0, NTrk = 0, useless = 0;
	int pdg_id_temp = -1;
	double px_temp = 0, py_temp = 0, pz_temp = 0;

	std::vector<double> px, py, pz, e;
	std::vector<int>    pdg_id;
	TLorentzVector      posFourMom, negFourMom;
	while (getline(infile, temp_string)) 
	{
		curstring.clear(); // needed when using several times istringstream::str(string)
		curstring.str(temp_string);

		if(strstr(temp_string.c_str(), "EVENT")) 
		{
			nEvts++;

			if(nEvts > nAnaEvts) break;
			if(nEvts%100000 == 0) cout<<"Working on "<< nEvts/100000 <<"-th 100k event ..."<<endl;

			// EVENT:          1       2       1
			curstring >> temp >> evtIdx >> NTrk >> useless; //assuming that EVENT always preceeds VERTEX/TRACK so that NTrk is set correctly

			pdg_id.clear();
			px.clear();
			py.clear();
			pz.clear();
			e.clear();

			posFourMom.SetPxPyPzE(0, 0, 0, 0);
			negFourMom.SetPxPyPzE(0, 0, 0, 0);
		}
		else if(strstr(temp_string.c_str(), "TRACK")) 
		{
			curstring >> temp >> useless >> px_temp >> py_temp >> pz_temp >> useless >> useless >> useless >> pdg_id_temp;

			TParticle particle_temp(pdg_id_temp, 0, 0, 0, 0, 0, px_temp, py_temp, pz_temp, 0.0, 0.0, 0.0, 0.0, 0.0);
			double mass_temp = particle_temp.GetMass();
			double e_temp    = TMath::Sqrt(pow(mass_temp, 2) + pow(particle_temp.P(), 2));

			px.push_back(px_temp);
			py.push_back(py_temp);
			pz.push_back(pz_temp);
			e.push_back(e_temp);
			pdg_id.push_back(pdg_id_temp);
		}

		if(NTrk == (int)px.size())
		{
			if(TMath::Abs(pdg_id[0])==13 && TMath::Abs(pdg_id[1])==13)
			{
				hnEvts->Fill(0);

				if(pdg_id[0] < 0){
					posFourMom.SetPxPyPzE(px[0], py[0], pz[0], e[0]);
					negFourMom.SetPxPyPzE(px[1], py[1], pz[1], e[1]);
				}
				else{
					posFourMom.SetPxPyPzE(px[1], py[1], pz[1], e[1]);
					negFourMom.SetPxPyPzE(px[0], py[0], pz[0], e[0]);
				}

				TLorentzVector motherFourMom = posFourMom + negFourMom;

				double pt   = motherFourMom.Pt();
				double mass = motherFourMom.M();
				double y    = motherFourMom.Rapidity();

				if(TMath::Abs(y) < mRapLow || TMath::Abs(y) > mRapHi) continue;

				hnEvts->Fill(1);
				hMvsPtvsRap_STARlight->Fill(y, pt, mass);

				// mMuPtTh, mMuEtaLow, and mMuEtaHi are used to filter muon kinematics during (private/official) MC production
				// see details in https://github.com/mumuustc/starlightLHE/blob/master/scriptForSimulation/Configuration/GenProduction/python/HINPbPbAutumn18GS_STARlight_fragment.py
				// see details in https://twiki.cern.ch/twiki/bin/viewauth/CMS/MCFor2018PbPb5TeV: the LHE-MuonCUT-fragment.py file in "Photo-production of charmonia in ultraperipheral PbPb collisions (Shuai Yang)" request
				if( posFourMom.Pt() >= mMuPtTh && TMath::Abs(posFourMom.Eta()) >= mMuEtaLow && TMath::Abs(posFourMom.Eta()) <= mMuEtaHi
						&& negFourMom.Pt() >= mMuPtTh && TMath::Abs(negFourMom.Eta()) >= mMuEtaLow && TMath::Abs(negFourMom.Eta()) <= mMuEtaHi
				  ){
					if(TMath::Abs(y) >= mRapLow && TMath::Abs(y) <= mRapHi) hnEvts->Fill(2);

					hMvsPtvsRap_muFilter->Fill(y, pt, mass);
					hNegPtvsPosPt->Fill(posFourMom.Pt(), negFourMom.Pt());
					hNegEtavsPosEta->Fill(posFourMom.Eta(), negFourMom.Eta());
					hNegPhivsPosPhi->Fill(posFourMom.Phi(), negFourMom.Phi());
				}
			}
			else{
				cout<<"The first two tracks are not muons !"<<endl;
			}
		}

	} // reading loop of the input file

	infile.close();

	system("mkdir -p mcAccHisto");

	TFile *fOut = new TFile(Form("mcAccHisto/%s.root", parSpec.Data()), "recreate");
	fOut->cd();
	hnEvts->Write();
	hMvsPtvsRap_STARlight->Write();
	hMvsPtvsRap_muFilter->Write();
	hNegPtvsPosPt->Write();
	hNegEtavsPosEta->Write();
	hNegPhivsPosPhi->Write();
	fOut->Close();
}
