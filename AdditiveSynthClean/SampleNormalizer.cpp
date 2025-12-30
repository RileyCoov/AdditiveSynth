/*=============================================================================================

	SampleNormalizer.cpp

	by Bob Coover

	Copyright:	© 2009 by Bob Coover dba Foxhound Sound, all rights reserved.

	This class will take an interleaved integer sample stream and convert it to
	a normalized, channel separated, float stream. It will also take a channel separated
	float sample stream and convert it to an interleaved integer sample stream.
	

=============================================================================================*/

#include "SampleNormalizer.h"

const int cSixteenFullScale = 32768;
const int cTwentyFourFullScale = 8388608;
const float cSixteenFloatFull = 32786.0f;
const float cTwentyFourFloatFull = 8388608.0f;

CSampleNormalizer::CSampleNormalizer(int bitDepth, int channels)
{
	mBitDepth = bitDepth;
	mChannels = channels;
}

CSampleNormalizer::~CSampleNormalizer()
{
	
}

void CSampleNormalizer::Normalize(float** samples, const short *file, long numSamples)
{
	// Assume sample interleaving if more than one channel of audio
	int ch, i, j;
	for (ch=0; ch<mChannels; ch++)
	{
		for (i=ch, j=0; i<numSamples*mChannels; i+=mChannels, j++)
		{
			samples[ch][j] = (float)file[i] / cSixteenFloatFull;
		}
	}	
}

void CSampleNormalizer::ConvertToShort(short* file, const float* const* samples, long numSamples)
{
	// Assume sample interleaving if more than one channel of audio
	int ch, i, j;
	for (ch=0; ch<mChannels; ch++)
	{
		for (i=ch, j=0; i<numSamples*mChannels; i+=mChannels, j++)
		{
			file[i] = (short)(samples[ch][j]*cSixteenFloatFull);
		}
	}	
}