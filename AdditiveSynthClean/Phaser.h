/*
 * Phaser.h - Six stage phaser effect with stereo split
 *
 *  Created on: 3/9/19
 *      Author: bcoover
 */

#ifndef PHASER_H_
#define PHASER_H_

#include "FXAlgo.h"

//************PHASE*SHIFTER*************************************
static const int cPhaseAllPassSize	= 3;
static const int cPhaseNumStages = 6;

static const float cPhaseCut1 =	240.0f;
static const float cPhaseCut2 =	240.0f;
static const float cPhaseCut3 = 245.0f;
static const float cPhaseCut4 = 245.0f;
static const float cPhaseCut5 = 245.5f;
static const float cPhaseCut6 = 245.0f;
//Right
static const float cPhaseCut7 = 230.0f;
static const float cPhaseCut8 = 230.0f;
static const float cPhaseCut9 = 245.0f;
static const float cPhaseCut10 = 245.0f;
static const float cPhaseCut11 = 245.0f;
static const float cPhaseCut12 = 245.0f;

static const int cStDlySize = 512;
static const int cStDlyMask = cStDlySize - 1;

class CPhaser: public CFXAlgo{
public:
	CPhaser(float sampleRate, int blockSize, int numChannels);
	virtual ~CPhaser();

	void SetPhaseDepth(float depth);
	void SetDepth(float val);
	void SetPhaseRate(float rate);
	void SetRate(float val);
	void SetMod2(float *mods);
	void ResetPhaseLFO();
	void SetPhaseFeedback(float amount);
	void SetFeedback(float val);
	int ProcessSamples(float **in_audio, float **out_audio, int numSamples);

private:
	void PhaseInit(void);
	void PhaserReset(void);

	float mRate;
	float mRawRate;
	float mRateMod[cBlockSize];
	float mRateFinal;
	float mRawDepth;
	float mDepth;
	float mDepthMod[cBlockSize];
	float mDepthFinal;
	float mLFOStep[2];
	float mLFOAccum;
	float mFeedbackTarget;
	float mRawFeedback;
	float mFeedback;
	float mFeedbackMod[cBlockSize];
	float mRads[2];
	float mMinRads[2];
	float mMaxRads[2];
	// Left
	float mAllPass1 [cPhaseAllPassSize];
	float mAllPass2 [cPhaseAllPassSize];
	float mAllPass3 [cPhaseAllPassSize];
	float mAllPass4 [cPhaseAllPassSize];
	float mAllPass5 [cPhaseAllPassSize];
	float mAllPass6 [cPhaseAllPassSize];
	// Right
	float mAllPass7 [cPhaseAllPassSize];
	float mAllPass8 [cPhaseAllPassSize];
	float mAllPass9 [cPhaseAllPassSize];
	float mAllPass10 [cPhaseAllPassSize];
	float mAllPass11 [cPhaseAllPassSize];
	float mAllPass12 [cPhaseAllPassSize];

	float mLeftShelf [5]; // 1st order hi shelf
	float mRightShelf [5]; // 1st order hi shelf

	// Stereo additions
	// left/right delay split
	int   mLRDelayWrite;
	float mStDelayBuff[cStDlySize];

};

#endif /* PHASER_H_ */

