/*
 * Phaser.cpp - Six stage phaser effect with stereo split
 *
 *  Created on: 3/9/19
 *      Author: bcoover
 */

#include "Phaser.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "UtilityFuncs.h"

using namespace UtilityFuncs;

CPhaser::CPhaser(float sampleRate, int blockSize, int numChannels)
: CFXAlgo(sampleRate, blockSize, numChannels)
{
    
    mLeft = eLEFT;
    mRight = eRIGHT;
    if (mNumChannels == 1)
        mRight = eLEFT;
	PhaseInit();
}

CPhaser::~CPhaser() {

}

void CPhaser::PhaserReset(void)			// called only by PhaseInit
{
	int i;
    // Zero the allpass states

    //Stereo Phaser Additions
    mLRDelayWrite = 0;
    mLFOAccum = 0.0f;
    memset(mStDelayBuff,0,sizeof(float)*cStDlySize);
    memset(mDepthMod, 0, sizeof(float)*cBlockSize);
    memset(mRateMod, 0, sizeof(float)*cBlockSize);
    memset(mFeedbackMod, 0, sizeof(float)*cBlockSize);

    mLeftShelf[3] = 0.0f; mRightShelf[3] = 0.0f;
    mLeftShelf[4] = 0.0f; mRightShelf[4] = 0.0f;

    for (i=0; i<2; i++)
    {
    	// Set the LFO sweep to the low point
    	mRads[i] = mMinRads[i];
    }

	for (i=1; i<cPhaseAllPassSize; i++)
	{
		mAllPass1[i] = 0.0f;
	    mAllPass2[i] = 0.0f;
	    mAllPass3[i] = 0.0f;
	    mAllPass4[i] = 0.0f;
	    mAllPass5[i] = 0.0f;
	    mAllPass6[i] = 0.0f;
	    mAllPass7[i] = 0.0f;
	    mAllPass8[i] = 0.0f;
	    mAllPass9[i] = 0.0f;
	    mAllPass10[i] = 0.0f;
	    mAllPass11[i] = 0.0f;
	    mAllPass12[i] = 0.0f;
	}
}

//**************************************************************************
void CPhaser::PhaseInit(void)	// called only on instantiation
{
	int i;
	float b0, b1, a1, rads;

	rads = cPiF * 4000.0f * mInvSampleRate;
	a1 = -sinf(cQuarterPiF - rads)/sinf(cQuarterPiF + rads);
	// 3.5db of gain
	b0 = 0.5f *(1.496235656094433f + 1.0f - a1*(1.496235656094433f - 1.0f));
	b1 = -0.5f *(1.496235656094433f - 1.0f - a1*(1.496235656094433f + 1.0f));


    // Set the coeffs
	mMinRads[0] = (cPiF * cPhaseCut1 * 2.0f) * mInvSampleRate;
	mMinRads[1] = (cPiF * cPhaseCut3 * 2.0f) * mInvSampleRate;

	PhaserReset();

	for (i=0; i<2; i++)
		mLFOStep[i] = 0.1f;
	mRawDepth = 0.75f;
	SetDepth(mRawDepth); //SetPhaseDepth(15.0f);
	mRawRate = 0.05f;
	SetRate(mRawRate); //SetPhaseRate(1.0f);


	mRawFeedback = 0.545f;
	SetFeedback(mRawFeedback); // Set the feedback target to 0.6
	mFeedback = 0.60f;		// feedback tweak here

	for (i=0; i<2; i++)
	{
		// Set the LFO sweep to the low point
		mRads[i] = mMinRads[i];
	}

	//Hi-shelf filter
	mLeftShelf[0] = b0; mRightShelf[0] = b0;
	mLeftShelf[1] = b1; mRightShelf[1] = b1;
	mLeftShelf[2] = a1; mRightShelf[2] = a1;

	// Stereo Additions
	//mModRate = 0.5f;//C_SP_MOD_RATE;
	//mStep = 0.5f*2.0f*C_SP_MOD_SIZE/SAMPLE_RATE_FLOAT;

	mAllPass2[0] = 0.0f; mAllPass4[0] = 0.0f; mAllPass5[0] = 0.0f; mAllPass6[0] = 0.0f;
	mAllPass7[0] = 0.0f; mAllPass8[0] = 0.0f; mAllPass9[0] = 0.0f; mAllPass10[0] = 0.0f;
	mAllPass11[0] = 0.0f; mAllPass12[0] = 0.0f;

	//DEBUG
	mAllPass1[0] = (1.0f - mRads[0])/(1.0f + mRads[0]);
	mAllPass3[0] = (1.0f - mRads[1])/(1.0f + mRads[1]);
}
//**************************************************************************
void CPhaser::SetPhaseDepth(float depth)	// change the modulation depth of the phase shift
{
	mRawDepth = depth;
}

void CPhaser::SetDepth(float val)
{
	//mDepthFinal = val*18.9f + 1.1;
	mDepthFinal = val*17.9f + 2.1;

	if (mDepthFinal < 0.0f)
		mDepthFinal *= -1.0f;
	mDepthFinal = fmax(mDepthFinal, 1.1f);
	mDepthFinal = fmin(mDepthFinal, 20.0f);

	//DebugPrintf(eDebugPrintfPriority_Notice, "Phaser depth = %f\n", mDepthFinal);

	mMaxRads[0] = (cPiF * cPhaseCut1 * mDepthFinal) * mInvSampleRate;
	mMaxRads[1] = (cPiF * cPhaseCut3 * mDepthFinal) * mInvSampleRate;
}


//**************************************************************************
void CPhaser::SetPhaseRate(float rate)		// 	Change the LFO rate of the phaser
{
	mRawRate = rate;
    printf("Phaser RAW rate = %f\n", mRawRate);
}

void CPhaser::SetRate(float val)
{
	float temp;

	// Cycles per second
	mRateFinal = val*7.95f + 0.05;;
	//if (mRateFinal < 0.0f)
		//mRateFinal *= -1.0f;
	mRateFinal = fmax(mRateFinal, 0.05f);
	mRateFinal = fmin(mRateFinal, 8.0f);

	//printf("Phaser rate = %f\n", mRateFinal);

	temp = (mMaxRads[0] - mMinRads[0])*mRateFinal*2.0f*mInvSampleRate;
	if (mLFOStep[0] < 0.0f)
		mLFOStep[0] = -1.0f*temp;
	else
		mLFOStep[0] = temp;

	temp = (mMaxRads[1] - mMinRads[0])*mRateFinal*2.0f*mInvSampleRate;
	if (mLFOStep[1] < 0.0f)
		mLFOStep[1] = -1.0f*temp;
	else
		mLFOStep[1] = temp;

}

//***************************************************************************
void CPhaser::SetMod2(float *mods)
{
	//mRateMod = modVal * 4.0f;
	int i;
	for (i=0; i<cBlockSize; i++)
		mRateMod[i] = mods[i];
	return;
}

//**************************************************************************
void CPhaser::ResetPhaseLFO()
{
	int i;

	for (i=0; i<2; i++)
	{
		mRads[i] = mMinRads[i];
		if (mLFOStep[i] < 0.0f)
			mLFOStep[i] *= -1.0f;

	}
}

//**************************************************************************
void CPhaser::SetPhaseFeedback(float amount)
{
    mRawFeedback = amount;
    printf("Phaser RAW feedback = %f\n", mRawFeedback);
}

//**************************************************************************
void CPhaser::SetFeedback(float val)
{
    mFeedbackTarget = val*0.55f + 0.3;
    mFeedbackTarget = fmin(mFeedbackTarget, 0.85f);
    mFeedbackTarget = fmax(mFeedbackTarget, 0.3f);

    //printf("Phaser feedback = %f\n", mFeedback);

}


//**************************************************************************
//**************************************************************************
int CPhaser::ProcessSamples(float **in_audio, float **out_audio, int numSamples)
{
	int i, j, index;
	//int mod, count;
    float feedInL, feedInR, outLeft, outRight, phaserOutL;
    //float frac, radianScale, sineVal;

	for (i=0; i<numSamples; i++)
	{
        CFXAlgo::SmoothParameters();
        
        outLeft = outRight = 0.0f;
        SetFeedback(mRawFeedback + mFeedbackMod[i]);
        SetDepth(mRawDepth + mDepthMod[i]);
        SetRate(mRawRate + mRateMod[i]);
		mFeedback = mFeedback + ( (mFeedbackTarget-mFeedback) * cSmoothCoeff);

		// Left
        mAllPass1[0] = (1.0f - mRads[0])/(1.0f + mRads[0]);
        mAllPass3[0] = (1.0f - mRads[1])/(1.0f + mRads[1]);

        feedInL = in_audio[mLeft][i] + mFeedback * mAllPass3[2];
        feedInR = in_audio[mRight][i] + mFeedback * mAllPass9[2];

        //Left
        mAllPass1[2] = mAllPass3[0]*feedInL - mAllPass1[1] + mAllPass3[0]*mAllPass1[2];
        mAllPass2[2] = mAllPass3[0]*mAllPass1[2] - mAllPass2[1] + mAllPass3[0]*mAllPass2[2];
        mAllPass3[2] = mAllPass3[0]*mAllPass2[2] - mAllPass3[1] + mAllPass3[0]*mAllPass3[2];
        mAllPass4[2] = mAllPass3[0]*mAllPass3[2] - mAllPass4[1] + mAllPass3[0]*mAllPass4[2];
        mAllPass5[2] = mAllPass3[0]*mAllPass4[2] - mAllPass5[1] + mAllPass3[0]*mAllPass5[2];
        mAllPass6[2] = mAllPass3[0]*mAllPass5[2] - mAllPass6[1] + mAllPass3[0]*mAllPass6[2];
        
        // Right
        mAllPass7[2] = mAllPass3[0]*feedInR - mAllPass7[1] + mAllPass3[0]*mAllPass7[2];
        mAllPass8[2] = mAllPass3[0]*mAllPass7[2] - mAllPass8[1] + mAllPass3[0]*mAllPass8[2];
        mAllPass9[2] = mAllPass3[0]*mAllPass8[2] - mAllPass9[1] + mAllPass3[0]*mAllPass9[2];
        mAllPass10[2] = mAllPass3[0]*mAllPass9[2] - mAllPass10[1] + mAllPass3[0]*mAllPass10[2];
        mAllPass11[2] = mAllPass3[0]*mAllPass10[2] - mAllPass11[1] + mAllPass3[0]*mAllPass11[2];
        mAllPass12[2] = mAllPass3[0]*mAllPass11[2] - mAllPass12[1] + mAllPass3[0]*mAllPass12[2];
        
        //Set the input state left
        mAllPass1[1] = feedInL;
        mAllPass2[1] = mAllPass1[2];
        mAllPass3[1] = mAllPass2[2];
        mAllPass4[1] = mAllPass3[2];
        mAllPass5[1] = mAllPass4[2];
        mAllPass6[1] = mAllPass5[2];
        //Set the input state right
        mAllPass7[1] = feedInR;
        mAllPass8[1] = mAllPass7[2];
        mAllPass9[1] = mAllPass8[2];
        mAllPass10[1] = mAllPass9[2];
        mAllPass11[1] = mAllPass10[2];
        mAllPass12[1] = mAllPass11[2];

    	// Sweep the frequencies left
        mRads[0] += mLFOStep[0];
        mRads[1] += mLFOStep[1];
        for (j=0; j<2; j++)
        {
        	if (mRads[j] <= mMinRads[j] && mLFOStep[j] < 0.0f)
        		mLFOStep[j] *= -1.0f;
        	if (mRads[j] >= mMaxRads[j] && mLFOStep[j] > 0.0f)
        		mLFOStep[j] *= -1.0f;
        }

        // Shelf
        outLeft = mLeftShelf[0]*mAllPass6[2] + mLeftShelf[1]*mLeftShelf[3] - mLeftShelf[2]*mLeftShelf[4];
        outRight = mRightShelf[0]*mAllPass12[2] + mRightShelf[1]*mRightShelf[3] - mRightShelf[2]*mRightShelf[4];

        mLeftShelf[3] = mAllPass6[2];
        mRightShelf[3] = mAllPass12[2];
        mLeftShelf[4] = outLeft;
        mRightShelf[4] = outRight;

 		// Left
        // Stereo Additions
        mStDelayBuff[mLRDelayWrite++] = outLeft * 0.80f;
        mLRDelayWrite &= cStDlyMask;

        index = floorf((float)mLRDelayWrite);
        index &= cStDlyMask;
        phaserOutL = mStDelayBuff[index];
        //outLeft = fx_wet_mix_current_f[fx_num] * mAllPass6[2] * 1.33f;	// uncomment for no delay
	   	outLeft = mMixDryCurr[i]*in_audio[mLeft][i] + mMixWetCurr[i] * 0.81f * phaserOutL;							// this line enables left channel delay
        
		//Right
		outRight = mMixDryCurr[i]*in_audio[mRight][i] + mMixWetCurr[i] * outRight * 0.65f;
		
		outLeft = fmax (outLeft, -1.0f);
		outLeft = fmin (outLeft, 0.999999f);
		out_audio[eLEFT][i] = outLeft;
        
		outRight = fmax (outRight, -1.0f);
		outRight = fmin (outRight, 0.999999f);
		out_audio[eRIGHT][i] = outRight;

	} // sample number

    return numSamples;
}





