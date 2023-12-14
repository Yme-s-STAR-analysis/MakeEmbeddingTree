#ifndef _StEmbeddingMaker_head
#define _StEmbeddingMaker_head
#include "StMaker.h"
#include "StThreeVectorF.hh"
#include "TString.h"
#include "TVector3.h"
#include "StarClassLibrary/StThreeVectorF.hh"
#include "StarClassLibrary/StThreeVectorD.hh"
#include "TNtuple.h"

class StPicoDst;
class StPicoEvent;
class StPicoTrack;
class StPicoMcTrack;
class StPicoDstMaker;
class TH1F;
class TH2F;
class TProfile;
class TTree;
class TH2D;
class TNtuple;

class StRefMultCorr;

class StEmbeddingMaker : public StMaker {
	public:
		StEmbeddingMaker(const char *name, StPicoDstMaker *picoMaker, const char *outName="tofMatchTree.root");
		virtual ~StEmbeddingMaker();

		virtual Int_t Init();
		virtual Int_t Make();
		virtual void  Clear(Option_t *opt="");
		virtual Int_t Finish();

		void set_target_ID(Int_t val) {
			targetID = val;
		}
		bool is_target_particle();
		bool is_McTrack_from_PV();


	private:
		StPicoDstMaker *mPicoDstMaker;
		StPicoDst      *mPicoDst;
		StPicoEvent	   *event;
		StPicoTrack    *rcTrack;
		StPicoMcTrack  *mcTrack;

		Int_t targetID;

		StRefMultCorr* corr;

		TNtuple* fDstTree;

		TString mOutputName;
		TFile* mFileOut;

		Int_t nEvents;

		ClassDef(StEmbeddingMaker, 1)
};


ClassImp(StEmbeddingMaker)

#endif
