/*==============================================================================

	UtilityFuncs.h

	by Bob Coover

	Copyright � 2006 FoxhoundSound 
	All rights reserved.

	Utility functions, constants, and classes.

==============================================================================*/
#pragma once

#ifndef UTILITY_FUNCS_H
#define UTILITY_FUNCS_H

#include <math.h>

namespace UtilityFuncs
{
    
    enum
    {
        eLEFT = 0,
        eRIGHT,
        eSTEREO
    };
    
    const double	cPi = 3.1415926535897932384626433832795;
	const double	cTwoPi = 6.283185307179586476925286766559;
	const float		cPiF = 3.14159265358979f;
	const float		cTwoPiF = 6.28318530717958f;
	const float		cQuarterPiF = 0.785398163397448f;
	const double	cNoDenromalDouble = 1.0e-33;
	const float		cNoDenormalFloat = 1.0e-20f;
	const double	cDoubleToIntBias = 6755399441055744.0; //(2^52)*1.5
	
	const float		cNeg2PIOverSampleRate = -cTwoPiF/44100.0f;
    const float     cNegPiOverNyquist = -cPiF/(44100.0f*0.5f);


	inline void NoDenormal(double* value)
	{
		if (*value < cNoDenromalDouble && *value > -cNoDenromalDouble)
			*value = 0.0;
	}

	inline void NoDenormalF(float* value)
	{
		if (*value < cNoDenormalFloat && *value > -cNoDenormalFloat)
			*value = 0.0;
	}

	inline int RoundDblToInt(double value)
	{
		value += cDoubleToIntBias;
		return (reinterpret_cast<int*>(&value))[0];
	}

	// circular buffer modulo
	inline int modc(int index, int size)
	{
		while (index<0)
			index = index + size;
		while (index>=size)
			index = index-size;
		return index;
	}
	
	// circular buffer modulo float
	inline float modcf(float index, float size)
	{
		while (index<0.0f)
			index = index + size;
		while (index>=size)
			index = index-size;
		return index;
	}

	template<class FLOAT_TYPE>
	class SmoothingFilter
	{
	public:
		SmoothingFilter() { Reset(); }
		void Reset() {mState = 0.0;}
		void InitState(double state) {mState = state;}
		double GetCurrentState(void) {return mState;}
		void SetSmoothingRate(FLOAT_TYPE freq, FLOAT_TYPE sampleRate) {mB1 = (FLOAT_TYPE)(1.0 - exp(freq*-cTwoPi/sampleRate)); }
        void SetSmoothingCoeff(FLOAT_TYPE coeff) {mB1 = coeff;}
		inline FLOAT_TYPE Smooth(FLOAT_TYPE input) { NoDenormal(&mState); return FLOAT_TYPE(mState += mB1 * (input - mState)); }
	private:
		double mState;
		FLOAT_TYPE mB1;
	};

	inline void HammingWindow(int length, float *hamm_window)
	{
		int half,i;
		half = length>>1;
		for (i=0; i<half; i++)
		{
			hamm_window[half-i-1] = 0.54f + 0.46f*(float)(cos(cTwoPi*i/(length)));
			hamm_window[half+i] = hamm_window[half-1-i];
		}
														
		
		return;
	}

	//-------------------------------------------------------------------------------------------
	//	Helper Function: Create a hanning window in the provided double array
	//-------------------------------------------------------------------------------------------
	inline void hanning (int length, double *han_window) {
		int half,i;

		if (!(length%2))	{		//even length
			half = length/2;
			for (i=0; i<half; i++) {
				han_window[i] = 0.5*(1-cos((cTwoPi*(i+1))/(length+1)));
				han_window[length-1-i] = han_window[i];
			}
		}
		else	{				//odd length		
			half = (length+1)/2;
			for (i=0; i<half; i++)
				han_window[i] = .5*(1-cos((cTwoPi*(i+1))/(length+1)));
			for (i=0; i<half-1; i++)
				han_window[length-1-i] = han_window[i];
		}
		return;
	}

	//-------------------------------------------------------------------------------------------
	//	Helper Function: Create a window the shape of a cosine from -pi/2 to pi/2 in the provided double array
	//-------------------------------------------------------------------------------------------
	inline void halfSine (int length, double *halfSin_window) {

		for (int i=0; i<length; i++) 
			halfSin_window[i] = sin((cPi*i)/(double)(length-1));

		return;
	}

} // end namespace UtilityFuncs

#endif
