/*=============================================================================================

	SampleNormalizer.h

	by Bob Coover

	Copyright:	© 2009 by Bob Coover dba Foxhound Sound, all rights reserved.

	This class will take an interleaved integer sample stream and convert it to
	a normalized, channel separated, float stream. It will also take a channel separated
	float sample stream and convert it to an interleaved integer sample stream.
	

=============================================================================================*/

#pragma once

#ifndef CSAMPLE_NORMALIZER_H
#define CSAMPLE_NORMALIZER_H

#include <math.h>

class CSampleNormalizer {

public:
	CSampleNormalizer(int bitDepth, int channels);
	~CSampleNormalizer();

	void Normalize(float** samples, const short* file, long numSamples);
	void ConvertToShort(short *file, const float* const* samples, long numSamples);

private:
	// member variables
	int	mBitDepth;
	int mChannels;
};

#endif