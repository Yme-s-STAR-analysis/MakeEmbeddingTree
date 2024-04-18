#include "StEmbeddingMaker.h"

#include <TMath.h>

#include <algorithm>
#include <fstream>
#include <vector>
#include <map>

#include "StBTofUtil/tofPathLength.hh"
#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEpdHit.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoPhysicalHelix.h"
#include "StPicoEvent/StPicoMcTrack.h"
#include "StPicoEvent/StPicoMcVertex.h"
#include "StThreeVectorF.hh"
#include "StLorentzVector.hh"
#include "Stiostream.h"
#include "TEfficiency.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TRandom.h"
#include "TTree.h"
#include "phys_constants.h"

#include "StRoot/CentCorrTool/CentCorrTool.h"
#include "StRoot/MeanDcaTool/MeanDcaTool.h"
#include "StRoot/TpcShiftTool/TpcShiftTool.h"
#include "StRoot/TriggerTool/TriggerTool.h"
#include "StRoot/StCFMult/StCFMult.h"

StEmbeddingMaker::StEmbeddingMaker(
	const char* name, 
	StPicoDstMaker* picoMaker,
    const char* outName
) : StMaker(name) {
	mOutputName = outName;
	mPicoDstMaker = picoMaker;
	mPicoDst = 0;
}

StEmbeddingMaker::~StEmbeddingMaker() {}

Int_t StEmbeddingMaker::Init() {
  	mFileOut = new TFile(mOutputName, "recreate");

	// patch 3.0: 
	// - 1) Add a flag 'match', 1 for RC matched, 0 for NOT matched
	// -> thus, we don't need MC information:
	// -> When a track is matched, we will use RC pt, y, eta...
	// -> When a track is NOT matched, we will use MC pt, y, eta...
	// -> And if it is NOT matched, the pt, y and eta of RC tracks would be ones from MC (indeed, there is no RC track for this case)
	// - 2) Some quantities will not be used, like nHitsPoss, dedx, dcaxy/z
	//	Previous quantities:
    //	"cent:vz:"
    //  "pTMc:etaMc:yMc:"
    //  "pTRc:etaRc:yRc:"
    // 	"nHitsFit:nHitsPoss:nHitsRatio:nHitsDedx:dedx:dca:dcaXY:dcaZ";
	
	const char* varList = 
       	"cent:vz:"
		"match:pT:eta:y:"
		"nHitsFit:nHitsDedx:nHitsRatio:dca";

	fDstTree = new TNtuple("fDstTree", "fDstTree", varList);

	// initialize costume modules

	// mean dca tool
	mtDca = new MeanDcaTool();
	mtDca->ReadParams();

	// centrality tool
	mtCent = new CentCorrTool();
	mtCent->EnableIndianMethod(true);
	mtCent->ReadParams();

	// multiplicity and shift tool
	mtShift = new TpcShiftTool();
	mtShift->Init();
	mtMult = new StCFMult();
	mtMult->ImportShiftTool(mtShift);

	// trigger tool
	mtTrg = new TriggerTool();

	nEvents = 0;

	return kStOK;
}

//---------------------------------------------------------
Int_t StEmbeddingMaker::Finish() {
	std::cout << "[LOG] Number of events: " << nEvents << ".\n";
	mFileOut->cd();
	fDstTree->Write();
	mFileOut->Close();
	std::cout << "[LOG] This is the end of this job.\n";
	return kStOK;
}

void StEmbeddingMaker::Clear(Option_t* opt) {}

//---------------------------------------------------------------
Int_t StEmbeddingMaker::Make() {
	if (!mPicoDstMaker) {
		LOG_WARN << " No PicoDstMaker! Skip! " << endm;
		return kStWarn;
	}

	mPicoDst = mPicoDstMaker->picoDst();
	if (!mPicoDst) {
		LOG_WARN << " No PicoDst! Skip! " << endm;
		return kStWarn;
	}

	if (!mPicoDst) {
		return kStOK;
	}

	// Load event
	event = (StPicoEvent*)mPicoDst->event();
	if (!event) {
		cerr << "Error opening picoDst Event, skip!" << endl;
		return kStOK;
	}
	nEvents += 1;

	TVector3 pVtx = event->primaryVertex();
	Double_t vx = pVtx.X();
	Double_t vy = pVtx.Y();
	Double_t vz = pVtx.Z();
	Float_t mB = event->bField();

	if (fabs(vx) < 1.e-5 && 
		fabs(vy) < 1.e-5 &&
		fabs(vz) < 1.e-5) {
		return kStOK;
	}

	// using Ashish's shifted vr cut
	// -> see: https://drupal.star.bnl.gov/STAR/system/files/Vr_xy_N_Vzcut.txt
	vx = vx - 0.0417;
	vy = vy + 0.2715;
	Double_t vr = sqrt(vx * vx + vy * vy);

	if (vr >= 1.0 || fabs(vz) > 50.0) {
		return kStOK;
	}

	Int_t runId = event->runId();
	Int_t trgid = mtTrg->GetTriggerID(event);
	if (trgid < 0) { return kStOK; }

	mtMult->make(mPicoDst);
	Int_t refMult = mtMult->mRefMult;
	Int_t tofMult = mtMult->mTofMult;
	Int_t nTofMatch = mtMult->mNTofMatch;
	Int_t nTofBeta = mtMult->mNTofBeta;

	Int_t refMult3 = mtMult->mRefMult3;
	refMult3 = mtCent->GetIndianRefMult3Corr(
		refMult, refMult3, tofMult, nTofMatch, nTofBeta,
		vz, false
	);
	if (refMult3 < 0) { return kStOK; }
	Int_t cent = mtCent->GetCentrality9(refMult3);
	if (cent < 0 || cent >= 9) { return kStOK; }

	// check DCA
	if (!mtDca->Make(mPicoDst)) { return kStOK; }
	if (mtDca->IsBadMeanDcaZEvent(mPicoDst) || mtDca->IsBadMeanDcaXYEvent(mPicoDst)) {
		return kStOK;
	}

	Int_t numberOfMcTracks = mPicoDst->numberOfMcTracks();
	Int_t numberOfRcTracks = mPicoDst->numberOfTracks();
	Int_t numberOfMcVertices = mPicoDst->numberOfMcVertices();
	if (!numberOfMcVertices || !numberOfMcTracks) { 
		// this event has no MC information, skip it
		return kStOK;
	}
	

	// Reconstructed track loop
	// to construct the map of RcTrack ID -> McTrack ID
	std::map<Int_t, Int_t> mMc2Rc;
	for (Int_t iRcTrk=0; iRcTrk<numberOfRcTracks; iRcTrk++) {
		rcTrack = (StPicoTrack*)mPicoDst->track(iRcTrk);
		if (!rcTrack) {
			continue;
		}
		if (!rcTrack->isPrimary()) {
			continue;
		}

		Int_t idTruth = rcTrack->idTruth(); // index of corresponding MC track
		if (idTruth <= 0 || idTruth > 10000) {
			continue;
		}

		// note that, here idTruth - 1 will be the quantity iMcTrk in later codes
		// not just idTruth. VERY IMPORTANT!
		mMc2Rc.insert(std::pair<Int_t, Int_t>(idTruth - 1, iRcTrk)); 
	}

	// MC track loop
	// to record the MC track and its RC track
	bool flag = false;
	for (Int_t iMcTrk=0; iMcTrk<numberOfMcTracks; iMcTrk++){
		mcTrack = (StPicoMcTrack*)mPicoDst->mcTrack(iMcTrk);
		if (!is_McTrack_from_PV()) {
			continue;
		}
		if (!is_target_particle()) {
			// only fill Geant ID histogram itself when the track is not proton
			continue;
		}
		if (mcTrack->idVtxStart() != 1) {
			// make sure this track is from PV
			continue;
		}

		// don't need this, the flag will be set later
		// if (mMc2Rc.count(iMcTrk) == 0) { // this mc track does not have rc track
		// 	flag = false;
		// }

		// find the reconstruct track
		// Indeed, mcTrack->id() is iMcTrk + 1;
		// cannot just use iRcTrk = mMc2Rc.find(iMcTrk)->second to get the RC track ID
		pair<map<Int_t, Int_t>::iterator, map<Int_t, Int_t>::iterator> ret;
		ret = mMc2Rc.equal_range(iMcTrk);
		map<Int_t, Int_t>::iterator iter;
		Int_t count = 0;
		Int_t iRcTrk = -1;
		Int_t iRcTrk_best = -1;
		Int_t qaTruth_best = -1;
		for (iter=ret.first; iter!=ret.second; iter++, count++) { // loop over the possible reconstructed tracks and find the best one with greatest qaTruth
			iRcTrk = iter->second;
			rcTrack = (StPicoTrack*)mPicoDst->track(iRcTrk);
			if (!rcTrack) {
				continue;
			}
			if (!rcTrack->isPrimary()){
				continue;
			}
			if (rcTrack->qaTruth() > qaTruth_best){
				qaTruth_best = rcTrack->qaTruth();
				iRcTrk_best = iRcTrk;
			}
		} 
		if (count > 0) { // indicates that rc track was found, and set the best one
			rcTrack = (StPicoTrack*)mPicoDst->track(iRcTrk_best);
			if (!rcTrack || !rcTrack->isPrimary()){
				flag = false;
			} else {
				flag = true;
			}
		} else {
			rcTrack = 0;
			flag = false;
		}

		TVector3 rcMom;
		StThreeVector<Float_t> rcMom3;
		StLorentzVector<Float_t> rcMom4;
		Double_t dca = -999.0;
		if (flag) {
			rcMom = rcTrack->pMom();
			rcMom3 = StThreeVector<Float_t>(rcMom.X(), rcMom.Y(), rcMom.Z());
			// calculate rapidity
			Float_t MP = 0.938272;
			Float_t EP = sqrt(rcMom3.mag2() + MP*MP);
			rcMom4 = StLorentzVector<Float_t>(rcMom3, EP);
			StPicoPhysicalHelix helix = rcTrack->helix(mB);
			dca = fabs(helix.geometricSignedDistance(pVtx));
		}

		// prepare the TNtuple
		Float_t array[20] = {0.0};
		Int_t idx = 0;
		// array[idx++] = (Float_t)RefMult3 * 1.0; // now we do not store RefMult3
		array[idx++] = (Float_t)cent * 1.0;
		array[idx++] = (Float_t)vz;
		array[idx++] = flag ? 1.0 : 0.0;
		array[idx++] = flag? (Float_t)rcMom.Pt() : (Float_t)mcTrack->pt();
		array[idx++] = flag? (Float_t)rcMom.PseudoRapidity() : (Float_t)mcTrack->eta();
		array[idx++] = flag? (Float_t)rcMom4.rapidity() : (Float_t)mcTrack->rapidity();
		array[idx++] = (Float_t)flag ? rcTrack->nHitsFit() * 1.0 : 0.0;
		array[idx++] = (Float_t)flag ? rcTrack->nHitsDedx() * 1.0 : 0.0;
		array[idx++] = (Float_t)flag ? (rcTrack->nHitsFit() * 1.0) / (rcTrack->nHitsPoss() * 1.0) : 0.0;
		array[idx++] = (Float_t)flag ? dca : 0.0;

		fDstTree->Fill(array);
	}

	return kStOK;
}

bool StEmbeddingMaker::is_target_particle() {
	return (mcTrack->geantId() == targetID);
}

bool StEmbeddingMaker::is_McTrack_from_PV() {
	if (!mcTrack) {
		return false;
	}
	Int_t idMcVx = mcTrack->idVtxStart();
	while (idMcVx != 1) {
		StPicoMcVertex* mcVertex = (StPicoMcVertex*)mPicoDst->mcVertex(idMcVx - 1);
		Int_t idMcTrack = mcVertex->idOfParentTrack();
		if (!idMcTrack) {
			break;
		}
		StPicoMcTrack* mcTrackP = (StPicoMcTrack*)mPicoDst->mcTrack(idMcTrack - 1);
		idMcVx = mcTrackP->idVtxStart();
		if (!idMcVx) {
			break;
		}
	}
	if (idMcVx != 1) {
		return false;
	}
	return true;
}