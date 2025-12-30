/*
 * FXAlgo.cpp
 *
 *  Created on: 9/22/18
 *      Author: bcoover
 */

#include "FXAlgo.h"
#include <string.h> // memset
#include <stdlib.h>
#include <math.h>
#include "Tables.h"
#include "UtilityFuncs.h"

using namespace UtilityFuncs;

CFXAlgo::CFXAlgo()
{
    mBlockSize = cBlockSize;
    
    this->mMixWetTarget = (float *) malloc(sizeof(float) * mBlockSize);
    this->mMixWetCurr = (float *) malloc(sizeof(float) * mBlockSize);
    this->mMixDryTarget = (float *) malloc(sizeof(float) * mBlockSize);
    this->mMixDryCurr = (float *) malloc(sizeof(float) * mBlockSize);
    
    memset(mMixWetTarget,0,mBlockSize*sizeof(float));
    memset(mMixWetCurr,0,mBlockSize*sizeof(float));
    memset(mMixDryTarget,1,mBlockSize*sizeof(float));
    memset(mMixDryCurr,1,mBlockSize*sizeof(float));
    
    mMixControl = 0.5f;
    mSampleRate = 44100.0f;
    mInvSampleRate = 1.0f/mSampleRate;
    mNumChannels = 2;
    mLeft = eLEFT;
    mRight = eRIGHT;
}


CFXAlgo::CFXAlgo(float sampleRate, int blockSize, int numChannels)
{
    int i;
    mBlockSize = blockSize;
    
    this->mMixWetTarget = (float *) malloc(sizeof(float) * mBlockSize);
    this->mMixWetCurr = (float *) malloc(sizeof(float) * mBlockSize);
    this->mMixDryTarget = (float *) malloc(sizeof(float) * mBlockSize);
    this->mMixDryCurr = (float *) malloc(sizeof(float) * mBlockSize);
    
    memset(mMixWetTarget,0,mBlockSize*sizeof(float));
    memset(mMixWetCurr,0,mBlockSize*sizeof(float));
    for (i=0; i<mBlockSize; i++)
    {
        mMixDryTarget[i] = 0.707f;
        mMixDryCurr[i] = 0.707f;
    }
    
    mMixControl = 0.0f;
    mSampleRate = sampleRate;
    mInvSampleRate = 1.0f/mSampleRate;
    mNumChannels = numChannels;
    mLeft = eLEFT;
    mRight = eRIGHT;
    if (mNumChannels == 1)
        mRight = eLEFT;
}

CFXAlgo::~CFXAlgo()
{
    if(mMixWetTarget)
        delete mMixWetTarget;
    if(mMixWetCurr)
        delete mMixWetCurr;
    if(mMixDryTarget)
        delete mMixDryTarget;
    if(mMixDryCurr)
        delete mMixDryCurr;
}

void CFXAlgo::SetSampleRate(float sampleRate)
{
    mSampleRate = sampleRate;
    mInvSampleRate = 1.0f/mSampleRate;
}


// This needs to be pure virtual - all effects need to clear their own state
void CFXAlgo::ClearDelayAndState(void)
{
    return;
}


void CFXAlgo::SetMix(float mix)
{
    mMixControl = mix;
}


void CFXAlgo::SmoothParameters()
{
    int i1, i2;
    float frac;
    
    for (int i=0; i<mBlockSize; i++)
    {
        frac = 127.0f * (mMixControl);
        // Clamp it
        if (frac > 127.0f)
            frac = 127.0f;
        else if (frac < 0.0f)
            frac = 0.0f;
        i1 = (int)frac;
        frac = frac - (float)i1;
        i2 = (i1 + 1)&127;
        mMixWetTarget[i] = cPanLawR[i1] + frac*(cPanLawR[i2] - cPanLawR[i1]); // sin pan law
        mMixWetCurr[i] = mMixWetCurr[i] + ( (mMixWetTarget[i] - mMixWetCurr[i]) * cSmopthMixCoeff ); // Smooth it
        mMixDryTarget[i] = cPanLawL[i1] + frac*(cPanLawL[i2] - cPanLawL[i1]); // cos pan law
        mMixDryCurr[i] = mMixDryCurr[i] + ( (mMixDryTarget[i] - mMixDryCurr[i]) * cSmopthMixCoeff ); // Smooth it
    }
}

int CFXAlgo::ProcessSamples(float **in_audio, float **out_audio, int numSamples)
{
    int i;
    for (i=0; i<numSamples; i++)
    {
        //out_audio[LEFT][i] = mMixDryCurr[i] * in_audio[LEFT][i] + 0.0f;
        //out_audio[RIGHT][i] = mMixDryCurr[i] * in_audio[RIGHT][i] + 0.0f;
        out_audio[eLEFT][i] = in_audio[eLEFT][i];
        out_audio[eRIGHT][i] = in_audio[eRIGHT][i];
    }
    
    return numSamples;
}




