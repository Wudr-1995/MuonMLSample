#include "MuonMLSmple.h"
#include "Event/LpmtElecTruthEvent.h"

MuonMLSmple::MuonMLSmple(// TString sf = "simLists",
						 TString cf = "calLists",
						 TString gf = "gdml",
						 TString of = "",
						 TString mo = "exa",
						 int ievt = 0) {
	// SimPaths = sf;
	CalPaths = cf;
	GdPath = gf;
	outPath = of;
	iEvt = ievt;
	mLSRadius = LSR;
	mode = mo;
	Ptab.reserve(NPMT);
	Ptab.resize(NPMT);
	ifstream gdml(GdPath);

	double n, theta, phi;
	for (int i = 0; i < NPMT; i ++) {
		gdml >> n >> theta >> phi;
		TVector3 tmp;
		tmp.SetMagThetaPhi(PMTR, theta, phi);
		Ptab[i].pos = tmp;
		Ptab[i].q = -1;
		Ptab[i].fht = 99999;
		Ptab[i].used = false;
	}

	gdml.close();
}

MuonMLSmple::~MuonMLSmple() {
}

bool MuonMLSmple::ObtainAndDump() {

	// Load files
	// ifstream sps(SimPaths);
	ifstream cps(CalPaths);

	// if (!sps.is_open()) {
	// 	cerr << "The sim paths file dosen's exist." << endl;
	// 	return false;
	// }

	if (!cps.is_open()) {
		cerr << "The calib paths file dosen't exist." << endl;
		return false;
	}

	// string simTmpPath;
	// while (sps.good()) {
	// 	getline(sps, simTmpPath);
	// 	if (simTmpPath.size() == 0)
	// 		continue;
	// 	SimList.push_back(simTmpPath);
	// }

	string calTmpPath;
	while (cps.good()) {
		getline(cps, calTmpPath);
		if (calTmpPath.size() == 0)
			continue;
		CalList.push_back(calTmpPath);
	}

	// sps.close();
	cps.close();

	TChain* simCh = new TChain("Event/Sim/SimEvent");
	for (vector<TString>::iterator it = CalList.begin(); it != CalList.end(); it ++) {
		simCh->Add(*it);
		clog << "File " << *it << "/Event/Sim/SimEvent added." << endl;
	}

	TChain* calCh = new TChain("Event/Calib/CalibEvent");
	for (vector<TString>::iterator it = CalList.begin(); it != CalList.end(); it ++) {
		calCh->Add(*it);
		clog << "File " << *it << "/Event/Calib/CalibEvent added." << endl;
	}

	// Load trees
	TTree* simTree = (TTree*)simCh;
	if (!simTree) {
		std::cerr << "Failed to load sim tree" << std::endl;
		return false;
	}

	TTree* calTree = (TTree*)calCh;
	if (!calTree) {
		cerr << "Failed to load calib tree" << endl;
		return false;
	}

	// Load branches
	JM::SimEvent* simObj = 0;
	simTree->SetBranchAddress("SimEvent", &simObj);
	simTree->GetBranch("SimEvent")->SetAutoDelete(true);

	JM::CalibEvent* calObj = 0;
	calTree->SetBranchAddress("CalibEvent", &calObj);
	calTree->GetBranch("CalibEvent")->SetAutoDelete(true);

	// Analysis
	int nEvts = simTree->GetEntries();
	int nCalEvts = calTree->GetEntries();
	if (nCalEvts != nEvts) {
		cerr << "The number of CalibEvent is not equal to the number of SimEvent." << endl;
		return false;
	}

	clog << nEvts << " events." << endl;

	bool saved = false;
	for (int i = 0; i < nEvts; i ++) {
		clog << "Processing " << i << "th events." << endl;
		simTree->GetEntry(i);
		calTree->GetEntry(i);

		const vector<JM::SimTrack*>& simTrks = simObj->getTracksVec();

		int nMuonTrks = 0;
		int size = simTrks.size();
		clog << "The number of tracks is " << size << endl;
		for (int j = 0; j < size; j ++) {
			if (simTrks[j]->getPDGID() == 13 || simTrks[j]->getPDGID() == - 13)
				nMuonTrks ++;
		}

		clog << "The number of muon tracks is " << nMuonTrks << endl;

		if (!freshPMTData(calObj))
			return false;
		auto q2d = new TH2F("q2d", ";#phi;#theta", 400, - PI, PI, 200, 0, PI);
		auto t2d20 = new TH2F("t2d20inch", ";#phi;#theta", 400, - PI, PI, 200, 0, PI);
		auto nPMTs = new TH2F("nPMTs", ";#phi;#theta", 400, - PI, PI, 200, 0, PI);
		auto n20PMTs = new TH2F("n20PMTs", ";#phi;#theta", 400, - PI, PI, 200, 0, PI);
		auto n3PMTs = new TH2F("n3PMTs", ";#phi;#theta", 400, - PI, PI, 200, 0, PI);
		for (PmtProp& p : Ptab) {
			if (p.res == _PMTINCH20 && p.used) {
				n20PMTs->Fill(p.pos.Phi(), p.pos.Theta());
				nPMTs->Fill(p.pos.Phi(), p.pos.Theta());
			}
			else if (p.res == _PMTINCH3 && p.used) {
				n3PMTs->Fill(p.pos.Phi(), p.pos.Theta());
				nPMTs->Fill(p.pos.Phi(), p.pos.Theta());
			}
			if (p.res == _PMTINCH3 || !p.used)
				continue;
			int biny = p.pos.Theta() * 199 / PI + 1;
			int binx = (p.pos.Phi() + PI) * 399 / 2 / PI + 1;
			q2d->Fill(p.pos.Phi(), p.pos.Theta(), p.q);
			double tmp = t2d20->GetBinContent(binx, biny);
			tmp = (tmp == 0 ? p.fht : (p.fht < tmp ? p.fht : tmp));
			t2d20->SetBinContent(binx, biny, tmp);
		}
		q2d->Scale(1 / (q2d->GetMaximum()));
		t2d20->Scale(1 / (t2d20->GetMaximum()));

		auto t2d3 = new TH2F("t2d3inch", ";#phi;#theta", 400, - PI, PI, 200, 0, PI);
		for (PmtProp& p : Ptab) {
			if (p.res == _PMTINCH20 || !p.used)
				continue;
			int biny = p.pos.Theta() * 199 / PI + 1;
			int binx = (p.pos.Phi() + PI) * 399 / 2 / PI + 1;
			double tmp = t2d3->GetBinContent(binx, biny);
			tmp = (tmp == 0 ? p.fht : (p.fht < tmp ? p.fht : tmp));
			t2d3->SetBinContent(binx, biny, tmp);
		}
		t2d3->Scale(1 / (t2d3->GetMaximum()));

		auto slopes = new TH2F("slopes", ";#phi;#theta", 400, - PI, PI, 200, 0, PI);
		for (PmtProp& p : Ptab) {
			if ((p.res != _PMTINCH3 && p.res != _PMTINCH20) || !p.used || !p.slope)
				continue;
			slopes->Fill(p.pos.Phi(), p.pos.Theta(), p.slope);
			// int biny = p.pos.Theta() * 199 / PI + 1;
			// int binx = (p.pos.Phi() + PI) * 399 / 2 / PI + 1;
			// double tmp = slopes->GetBinContent(binx, biny);
			// slopes->SetBinContent(binx, biny, p.slope / nPMTs->GetBinContent(binx, biny) + tmp);
		}
		slopes->Scale(1 / (slopes->GetMaximum()));
		// slopes->SetMaximum(20);
		// slopes->SetMinimum(-0.01);

		auto slopes20 = new TH2F("slopes20", ";#phi;#theta", 400, - PI, PI, 200, 0, PI);
		for (PmtProp& p : Ptab) {
			if (p.res == _PMTINCH3 || !p.used ||!p.slope)
				continue;
			slopes20->Fill(p.pos.Phi(), p.pos.Theta(), p.slope);
			// int biny = p.pos.Theta() * 199 / PI + 1;
			// int binx = (p.pos.Phi() + PI) * 399 / 2 / PI + 1;
			// double tmp = slopes20->GetBinContent(binx, biny);
			// slopes20->SetBinContent(binx, biny, p.slope / n20PMTs->GetBinContent(binx, biny) + tmp);
		}
		slopes20->Scale(1 / (slopes20->GetMaximum()));

		auto slopes3 = new TH2F("slopes3", ";#phi;#theta", 400, - PI, PI, 200, 0, PI);
		for (PmtProp& p : Ptab) {
			if (p.res == _PMTINCH20 || !p.used || !p.slope)
				continue;
			slopes3->Fill(p.pos.Phi(), p.pos.Theta(), p.slope);
			// int biny = p.pos.Theta() * 199 / PI + 1;
			// int binx = (p.pos.Phi() + PI) * 399 / 2 / PI + 1;
			// double tmp = slopes3->GetBinContent(binx, biny);
			// slopes3->SetBinContent(binx, biny, p.slope / n3PMTs->GetBinContent(binx, biny) + tmp);
		}
		slopes3->Scale(1 / (slopes3->GetMaximum()));

		int nTrks = 0;
		TString samName = outPath + "QTT" + iEvt;
		ofstream sample(samName);
		for (int i = 0; i < simTrks.size(); i ++) {
			if (TMath::Abs(simTrks[i]->getPDGID()) != 13)
				continue;
			TVector3 simInit(simTrks[i]->getInitX(), simTrks[i]->getInitY(), simTrks[i]->getInitZ());
			TVector3 simExit(simTrks[i]->getExitX(), simTrks[i]->getExitY(), simTrks[i]->getExitZ());
			TVector3 simDir(simTrks[i]->getInitPx(), simTrks[i]->getInitPy(), simTrks[i]->getInitPz());

			double dis = TMath::Sqrt(simInit.Mag2() - TMath::Power(simInit * simDir, 2));
			if (dis > LSR || simInit.Mag() < LSR || simExit.Mag() < LSR)
				continue;
			nTrks ++;
		}

		if (!nTrks) {
			delete q2d;
			delete t2d3;
			delete t2d20;
			delete slopes20;
			delete slopes3;
			delete slopes;
			delete nPMTs;
			delete n20PMTs;
			delete n3PMTs;
			continue;
		}

		for (int i = 1; i <= 200; i ++) {
			for (int j = 1; j <= 400; j ++) {
				sample << q2d->GetBinContent(i, j);
				if (j != 400)
					sample << "\t";
			}
			sample << endl;
		}

		for (int i = 1; i <= 200; i ++) {
			for (int j = 1; j <= 400; j ++) {
				sample << t2d20->GetBinContent(i, j);
				if (j != 400)
					sample << "\t";
			}
			sample << endl;
		}

		for (int i = 1; i <= 200; i ++) {
			for (int j = 1; j <= 400; j ++) {
				sample << t2d3->GetBinContent(i, j);
				if (j != 400)
					sample << "\t";
			}
			sample << endl;
		}

		for (int i = 1; i <= 200; i ++) {
			for (int j = 1; j <= 400; j ++) {
				sample << slopes20->GetBinContent(i, j);
				if (j != 400)
					sample << "\t";
			}
			sample << endl;
		}

		for (int i = 1; i <= 200; i ++) {
			for (int j = 1; j <= 400; j ++) {
				sample << slopes3->GetBinContent(i, j);
				if (j != 400)
					sample << "\t";
			}
			sample << endl;
		}

		for (int i = 1; i <= 200; i ++) {
			for (int j = 1; j <= 400; j ++) {
				sample << slopes->GetBinContent(i, j);
				if (j != 400)
					sample << "\t";
			}
			sample << endl;
		}

		for (int i = 0; i < simTrks.size(); i ++) {
			if (TMath::Abs(simTrks[i]->getPDGID()) != 13)
				continue;

			TVector3 simInit(simTrks[i]->getInitX(), simTrks[i]->getInitY(), simTrks[i]->getInitZ());
			TVector3 simExit(simTrks[i]->getExitX(), simTrks[i]->getExitY(), simTrks[i]->getExitZ());
			TVector3 simDir(simTrks[i]->getInitPx(), simTrks[i]->getInitPy(), simTrks[i]->getInitPz());

			double dis = TMath::Sqrt(simInit.Mag2() - TMath::Power(simInit * simDir, 2));
			if (dis > LSR || simInit.Mag() < LSR || simExit.Mag() < LSR)
				continue;

			sample << simInit.Theta() << "\t" << simInit.Phi() << "\t"
				   << simDir.Theta() << "\t" << simDir.Phi() << endl;
		}
		sample << nTrks << endl;

		if (mode == "exa" && !saved && nTrks && iEvt == 1) {
			saved = true;
			auto hits = new TH1F("HitVCharge", ";Hit time / ns;nPE / p.e.", 5000, 0, 1000);
			int nhit = hitTimeExa.size();
			for (int i = 0; i < nhit; i ++) {
				hits->SetBinContent(hitTimeExa[i] * 4999 / 1000 + 1,
									hitChargeExa[i] + hits->GetBinContent(hitTimeExa[i] * 4999 / 1000 + 1));
			}
			DrawExa(q2d, t2d20, t2d3, hits, slopes);
			delete hits;
		}

		sample.close();

		delete q2d;
		delete t2d3;
		delete t2d20;
		delete slopes20;
		delete slopes3;
		delete slopes;
		delete nPMTs;
		delete n20PMTs;
		delete n3PMTs;

		iEvt ++;
	}
	return true;
}

bool MuonMLSmple::freshPMTData(JM::CalibEvent* calEvt) {
	for (int i = 0; i < NPMT; i ++) {
		Ptab[i].used = false;
	}

	const list<JM::CalibPMTChannel*>& chList = calEvt->calibPMTCol();
	list<JM::CalibPMTChannel*>::const_iterator chit = chList.begin();
	bool saved = false;
	while (chit != chList.end()) {
		JM::CalibPMTChannel* ch = *chit ++;
		if (ch->nPE() > 2000 && !saved) {
			saved = true;
			const vector<double>& times = ch->time();
			const vector<double>& charges = ch->charge();
			hitTimeExa.assign(times.begin(), times.end());
			hitChargeExa.assign(charges.begin(), charges.end());
		}
		auto tmph = new TH1F("tmph", "", 2500, 0, 1000);
		const vector<double>& times = ch->time();
		const vector<double>& charges = ch->charge();
		int size = times.size();
		for (int i = 0; i < size; i ++) {
			tmph->SetBinContent(times[i] * 2499 / 1000 + 1,
								charges[i] + tmph->GetBinContent(times[i] * 2499 / 1000 + 1));
		}
		Identifier id = Identifier(ch->pmtId());
		Identifier::value_type value = id.getValue();
		unsigned int pid;
		if (not ((value & 0xFF000000) >> 24 == 0x10)) {
			delete tmph;
			continue;
		}
		pid = CdID::module(id);
		if (pid > NPMT) {
			cerr << "Mis-match: PmtId " << pid << " which exceed the range of the number of PMT." << endl;
			return false;
		}
		Ptab[pid].q = ch->nPE();
		Ptab[pid].fht = ch->firstHitTime();
		double v = tmph->GetMaximumBin() - Ptab[pid].fht;
		if (tmph->GetMaximumBin() - Ptab[pid].fht == 0)
			value = 1;
		Ptab[pid].slope = tmph->GetMaximum() / value;
		if (isinf(Ptab[pid].slope))
			Ptab[pid].slope = 0;
		if (Ptab[pid].q == 0)
			Ptab[pid].slope = 0;
		delete tmph;
		Ptab[pid].loc = 1;
		if (CdID::is3inch(id)) {
			Ptab[pid].res = _PMTINCH3;
			Ptab[pid].used = true;
		}
		else if (CdID::is20inch(id)) {
			Ptab[pid].res = _PMTINCH20;
			Ptab[pid].used = true;
		}
		else {
			Ptab[pid].res = _PMTNULL;
			Ptab[pid].used = false;
		}
	}

	return true;
}

bool MuonMLSmple::DrawExa(TH2F* q2d, TH2F* t2d20, TH2F* t2d3, TH1F* hits, TH2F* slopes20) {

	gStyle->SetOptStat(0000);
	gStyle->SetPalette(1);
	TH2F* q2dc = (TH2F*)q2d->Clone("q2dClone");
	TH2F* t2dc20 = (TH2F*)t2d20->Clone("t2dClone");
	TH2F* t2dc3 = (TH2F*)t2d3->Clone("t2dClone");
	TH1F* hitsc = (TH1F*)hits->Clone("hitsClone");
	TH2F* slopes20c = (TH2F*)slopes20->Clone("slopes20c");

	clog << "Printing......" << endl;
	auto c = new TCanvas("c", "", 800, 500);
	c->SetTopMargin(0.1);
	c->SetLeftMargin(0.2);
	c->SetBottomMargin(0.2);
	c->SetRightMargin(0.15);
	TString pdfPath = "/junofs/users/wudr/Muon_Sim/Real_Muon/Rec_job/MuonMLSmple/test/Example.pdf";
	c->Print(pdfPath + "[");

	c->cd();
	q2dc->GetXaxis()->SetLabelSize(0.06);
	q2dc->GetYaxis()->SetLabelSize(0.06);
	q2dc->GetXaxis()->SetTitleSize(0.08);
	q2dc->GetYaxis()->SetTitleSize(0.08);
	q2dc->Draw("colz");
	c->Print(pdfPath);

	c->cd();
	t2dc20->GetXaxis()->SetLabelSize(0.06);
	t2dc20->GetYaxis()->SetLabelSize(0.06);
	t2dc20->GetXaxis()->SetTitleSize(0.08);
	t2dc20->GetYaxis()->SetTitleSize(0.08);
	t2dc20->Draw("colz");
	c->Print(pdfPath);

	c->cd();
	t2dc3->GetXaxis()->SetLabelSize(0.06);
	t2dc3->GetYaxis()->SetLabelSize(0.06);
	t2dc3->GetXaxis()->SetTitleSize(0.08);
	t2dc3->GetYaxis()->SetTitleSize(0.08);
	t2dc3->Draw("colz");
	c->Print(pdfPath);

	c->cd();
	hitsc->GetXaxis()->SetLabelSize(0.06);
	hitsc->GetYaxis()->SetLabelSize(0.06);
	hitsc->GetXaxis()->SetTitleSize(0.08);
	hitsc->GetYaxis()->SetTitleSize(0.08);
	hitsc->Draw();
	c->Print(pdfPath);

	c->cd();
	slopes20c->GetXaxis()->SetLabelSize(0.06);
	slopes20c->GetYaxis()->SetLabelSize(0.06);
	slopes20c->GetXaxis()->SetTitleSize(0.08);
	slopes20c->GetYaxis()->SetTitleSize(0.08);
	slopes20c->Draw("colz");
	c->Print(pdfPath);

	c->Print(pdfPath + "]");

	return true;
}
