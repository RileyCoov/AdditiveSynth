/*
 * FXAlgo.h
 *
 *  Created on: 9/22/18
 *      Author: bcoover
 */

#ifndef FXALGO_H_
#define FXALGO_H_

const int cBlockSize = 512;
const float cSmopthMixCoeff = 0.001f;
const float cSmoothCoeff = 0.001f;
const float c0dBFSRefVoltage = 5.0f;
const float	cControlMin = 0.0;
const float	cControlMid = 5.0;
const float cControlMax = 10.0;
const float cControlScaler = 1.0f/cControlMax;

class CFXAlgo {
public:
    CFXAlgo();
    CFXAlgo(float sampleRate, int blockSize, int numChannels);
    virtual ~CFXAlgo();
    void SetMix(float mix);
    virtual void SmoothParameters();
    virtual int ProcessSamples(float **in_audio, float **out_audio, int numSamples);
    virtual void SetSampleRate(float sampleRate);
    virtual void ClearDelayAndState(void);
    
protected:
    float mSampleRate;
    int mBlockSize;
    int mNumChannels;
    int mLeft;
    int mRight;
    float mInvSampleRate;
    float mMixControl;
    float * mMixWetTarget;
    float * mMixWetCurr;
    float * mMixDryTarget;
    float * mMixDryCurr;
    
};

#endif /* FXALGO_H_ */
