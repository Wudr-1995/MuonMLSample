#ifndef MuonMLSmple_h
#define MuonMLSmple_h

#include <iostream>
#include <map>
#include "TFile.h"
#include "TTree.h"
#include "Event/CalibEvent.h"
#include "Event/ElecEvent.h"
#include "Event/SimEvent.h"
#include "Event/WPRecEvent.h"
#include "Event/CDTrackRecEvent.h"
#include "TString.h"
#include "Identifier/Identifier.h"
#include "Identifier/CdID.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TMath.h"
#include "TH2F.h"
#include <stdio.h>
#include "TStyle.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TArrow.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TChain.h"
#include "TVector3.h"
#include "TStyle.h"
#include <fstream>
#include <vector>
#include <list>

#define UNIT TMath::Pi() / 180
#define PI TMath::Pi()
#define NPMT 43212
#define LSR 17700.
#define PMTR 19500.

using namespace std;

enum Pmttype {
	_PMTNULL,
	_PMTINCH3,
	_PMTINCH20,
	_PMTTT,
};
struct PmtProp{
	TVector3 pos;
	double q;
	double fht;
	double res;
	double trise;
	double amp;
	double integ;
	double slope;
	bool used;
	// Add loc: 1 == CD, 2 == WP, 3 == TT, by wudr
	int loc;
	int value;
	Pmttype type;
};
typedef std::vector<PmtProp> PmtTable;

class MuonMLSmple {
	public:
		MuonMLSmple(TString, TString, TString, TString, int);
		~MuonMLSmple();

		bool ObtainAndDump();

	private:
		bool freshPMTData(JM::CalibEvent*);
		bool DrawExa(TH2F*, TH2F*, TH2F*, TH1F*, TH2F*);
		bool GetSlope(JM::CalibPMTChannel*);

		int iEvt;

		TString SimPaths;
		TString CalPaths;
		TString outPath;
		TString GdPath;
		TString mode;

		PmtTable Ptab;

		vector<TString> SimList;
		vector<TString> CalList;

		vector<double> hitTimeExa;
		vector<double> hitChargeExa;

		double mLSRadius;

		vector<TVector3> SimInits;
		vector<TVector3> SimDirs;

};

#endif
