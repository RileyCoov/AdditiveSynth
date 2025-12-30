/*==============================================================================

	RealFFT.h

	by Bob Coover

	Copyright © 2002 by Bob Coover dba Foxhound Sound, all rights reserved. 
	All rights reserved.

	This is a real signal FFT function.

==============================================================================*/

#pragma once

#ifndef REAL_FFT_H
#define REAL_FFT_H

#include "UtilityFuncs.h"
#include <vector>

namespace UtilityFuncs
{


// RealFFT
// ========== 
// Sorensen in-place split-radix FFT for real values
// data: array of floats:
// re(0),re(1),re(2),...,re(size-1)
// 
// output:
// re(0),re(1),re(2),...,re(size/2),im(size/2-1),...,im(1)
// normalized by array length
void RealFFT(float *samples, int len);
void RealFFT(float *input, double *output, int len);
void InvRealFFT(float *samples, int len);
void MagnitudeFFT(float *samples, int len);
void MagnitudeFFTVec(std::vector<float> &samples);

// Miscellaneous
// =============
unsigned int Log2(unsigned int word);

inline
unsigned int Log2(unsigned int word)
{
	if(word == 0)
		throw -1; // WARNING : HACK!!!
		
	unsigned int l2 = 0;
	
	// WARNING : Assumes an unsigned int is 32 bits
	if(word >= (1 << 16))
	{
		word >>= 16;
		l2 += 16;
	};
	if(word >= (1 << 8))
	{
		word >>= 8;
		l2 += 8;
	};
	if(word >= (1 << 4))
	{
		word >>= 4;
		l2 += 4;
	};
	if(word >= (1 << 2))
	{
		word >>= 2;
		l2 += 2;
	};
	if(word >= (1 << 1))
	{
		word >>= 1;
		l2 += 1;
	};
	return l2;
}

} // namespace UtilityFuncs

#endif
