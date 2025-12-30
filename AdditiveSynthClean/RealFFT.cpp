/*==============================================================================

	RealFFT.cpp

	by Bob Coover

	Copyright © 2011 by Bob Coover dba Foxhound Sound, all rights reserved. 
	All rights reserved.

	This is an FFT for real signals, (complex signals will not work).

==============================================================================*/

#include "RealFFT.h"

// Standard headers
#include <algorithm>
#include <cmath>


//For Windows disable certain warnings
#if defined( WINDOWS_VERSION )
#pragma warning( disable: 4244 ) // double to float warning
#endif

namespace UtilityFuncs
{

	/////////////////////////////////////////////////////////
	// Sorensen in-place split-radix FFT for real values
	// samples: array of floats:
	// re(0),re(1),re(2),...,re(size-1)
	// 
	// output:
	// re(0),re(1),re(2),...,re(size/2),im(size/2-1),...,im(1)
	// normalized by array length
	//
	// Source: 
	// Sorensen et al: Real-Valued Fast Fourier Transform Algorithms,
	// IEEE Trans. ASSP, ASSP-35, No. 6, June 1987

	void RealFFT(float *samples,int len)
	{
		int i,j,k,i5,i6,i7,i8,i0,id,i1,i2,i3,i4,n2,n4,n8;
		double t1,t2,t3,t4,t5,t6,a3,ss1,ss3,cc1,cc3,a,e,sqrt2;

		sqrt2=sqrt(2.0);
		n4=len-1;

		//data shuffling
		for (i=0,j=0,n2=len/2; i<n4 ; i++){
			if (i<j){
				t1=samples[j];
				samples[j]=samples[i];
				samples[i]=t1;
			}
			k=n2;
			while (k<=j){
				j-=k;
				k>>=1;	
			}
			j+=k;
		}

		/*----------------------*/

		//length two butterflies	
		i0=0;
		id=4;
		do{
			for (; i0<n4; i0+=id){ 
				i1=i0+1;
				t1=samples[i0];
				samples[i0]=t1+samples[i1];
				samples[i1]=t1-samples[i1];
			}
			id<<=1;
			i0=id-2;
			id<<=1;
		} while ( i0<n4 );

		/*----------------------*/
		//L shaped butterflies
		n2=2;
		for(k=len;k>2;k>>=1){  
			n2<<=1;
			n4=n2>>2;
			n8=n2>>3;
			e = 2*cPi/(n2);
			i1=0;
			id=n2<<1;
			do{ 
				for (; i1<len; i1+=id){
					i2=i1+n4;
					i3=i2+n4;
					i4=i3+n4;
					t1=samples[i4]+samples[i3];
					samples[i4]-=samples[i3];
					samples[i3]=samples[i1]-t1;
					samples[i1]+=t1;
					if (n4!=1){
						i0=i1+n8;
						i2+=n8;
						i3+=n8;
						i4+=n8;
						t1=(samples[i3]+samples[i4])/sqrt2;
						t2=(samples[i3]-samples[i4])/sqrt2;
						samples[i4]=samples[i2]-t1;
						samples[i3]=-samples[i2]-t1;
						samples[i2]=samples[i0]-t2;
						samples[i0]+=t2;
					}
				}
				id<<=1;
				i1=id-n2;
				id<<=1;
			} while ( i1<len );
			a=e;
			for (j=2; j<=n8; j++){  
				a3=3*a;
				cc1=cos(a);
				ss1=sin(a);
				cc3=cos(a3);
				ss3=sin(a3);
				a=j*e;
				i=0;
				id=n2<<1;
				do{
					for (; i<len; i+=id){  
						i1=i+j-1;
						i2=i1+n4;
						i3=i2+n4;
						i4=i3+n4;
						i5=i+n4-j+1;
						i6=i5+n4;
						i7=i6+n4;
						i8=i7+n4;
						t1=samples[i3]*cc1+samples[i7]*ss1;
						t2=samples[i7]*cc1-samples[i3]*ss1;
						t3=samples[i4]*cc3+samples[i8]*ss3;
						t4=samples[i8]*cc3-samples[i4]*ss3;
						t5=t1+t3;
						t6=t2+t4;
						t3=t1-t3;
						t4=t2-t4;
						t2=samples[i6]+t6;
						samples[i3]=t6-samples[i6];
						samples[i8]=t2;
						t2=samples[i2]-t3;
						samples[i7]=-samples[i2]-t3;
						samples[i4]=t2;
						t1=samples[i1]+t5;
						samples[i6]=samples[i1]-t5;
						samples[i1]=t1;
						t1=samples[i5]+t4;
						samples[i5]-=t4;
						samples[i2]=t1;
					}
					id<<=1;
					i=id-n2;
					id<<=1;
				} while(i<len);
			}
		}

		//division with array length
		//float normalizer = 1.0/len;
		//for(i=0;i<len;i++) samples[i]*=normalizer;
	}


	/////////////////////////////////////////////////////////
	// Sorensen in-place split-radix FFT for real values
	// samples: array of doubles:
	// re(0),re(1),re(2),...,re(size-1)
	// 
	// output:
	// re(0),re(size/2),re(1),im(1),re(2),im(2),...,re(size/2-1),im(size/2-1)
	// normalized by array length
	//
	// Source: 
	// Source: see the routines it calls ...

	void RealFFT_Unshuffled(float *input,double *output, int len)
	{
		int i,half;

		RealFFT(input,len);
		//unshuffling - not in-place
		half=len/2;
		output[0]=input[0];
		output[1]=input[half];
		for(i=1;i<half;i++) {output[i+i]=input[i];output[i+i+1]=input[len-i];}
	}

	/////////////////////////////////////////////////////////
	// Sorensen in-place inverse split-radix FFT for real values
	// samples: array of doubles:
	// re(0),re(1),re(2),...,re(size/2),im(size/2-1),...,im(1)
	// 
	// output:
	// re(0),re(1),re(2),...,re(size-1)
	// NOT normalized by array length
	//
	// Source: 
	// Sorensen et al: Real-Valued Fast Fourier Transform Algorithms,
	// IEEE Trans. ASSP, ASSP-35, No. 6, June 1987

	void InvRealFFT(float *samples,int len)
	{

		int i,j,k,i5,i6,i7,i8,i0,id,i1,i2,i3,i4,n2,n4,n8,n1;
		float t1,t2,t3,t4,t5,a3,ss1,ss3,cc1,cc3,a,e,sqrt2;

		sqrt2=sqrt(2.0);

		n1=len-1;
		n2=len<<1;
		for(k=len;k>2;k>>=1){  
			id=n2;
			n2>>=1;
			n4=n2>>2;
			n8=n2>>3;
			e = 2*cPi/(n2);
			i1=0;
			do{ 
				for (; i1<len; i1+=id){
					i2=i1+n4;
					i3=i2+n4;
					i4=i3+n4;
					t1=samples[i1]-samples[i3];
					samples[i1]+=samples[i3];
					samples[i2]*=2;
					samples[i3]=t1-2*samples[i4];
					samples[i4]=t1+2*samples[i4];
					if (n4!=1){
						i0=i1+n8;
						i2+=n8;
						i3+=n8;
						i4+=n8;
						t1=(samples[i2]-samples[i0])/sqrt2;
						t2=(samples[i4]+samples[i3])/sqrt2;
						samples[i0]+=samples[i2];
						samples[i2]=samples[i4]-samples[i3];
						samples[i3]=2*(-t2-t1);
						samples[i4]=2*(-t2+t1);
					}
				}
				id<<=1;
				i1=id-n2;
				id<<=1;
			} while ( i1<n1 );
			a=e;
			for (j=2; j<=n8; j++){  
				a3=3*a;
				cc1=cos(a);
				ss1=sin(a);
				cc3=cos(a3);
				ss3=sin(a3);
				a=j*e;
				i=0;
				id=n2<<1;
				do{
					for (; i<len; i+=id){  
						i1=i+j-1;
						i2=i1+n4;
						i3=i2+n4;
						i4=i3+n4;
						i5=i+n4-j+1;
						i6=i5+n4;
						i7=i6+n4;
						i8=i7+n4;
						t1=samples[i1]-samples[i6];
						samples[i1]+=samples[i6];
						t2=samples[i5]-samples[i2];
						samples[i5]+=samples[i2];
						t3=samples[i8]+samples[i3];
						samples[i6]=samples[i8]-samples[i3];
						t4=samples[i4]+samples[i7];
						samples[i2]=samples[i4]-samples[i7];
						t5=t1-t4;
						t1+=t4;
						t4=t2-t3;
						t2+=t3;
						samples[i3]=t5*cc1+t4*ss1;
						samples[i7]=-t4*cc1+t5*ss1;
						samples[i4]=t1*cc3-t2*ss3;
						samples[i8]=t2*cc3+t1*ss3;
			  }
					id<<=1;
					i=id-n2;
					id<<=1;
				} while(i<n1);
			}
		}	

		/*----------------------*/
		i0=0;
		id=4;
		do{
			for (; i0<n1; i0+=id){ 
				i1=i0+1;
				t1=samples[i0];
				samples[i0]=t1+samples[i1];
				samples[i1]=t1-samples[i1];
			}
			id<<=1;
			i0=id-2;
			id<<=1;
		} while ( i0<n1 );

		/*----------------------*/

		//data shuffling
		for (i=0,j=0,n2=len/2; i<n1 ; i++){
			if (i<j){
				t1=samples[j];
				samples[j]=samples[i];
				samples[i]=t1;
			}
			k=n2;
			while (k<=j){
				j-=k;
				k>>=1;	
			}
			j+=k;
		}

		//division with array length
		float normalizer = 1.0/len;
		for(i=0;i<len;i++) samples[i]*=normalizer;

	}

    // In place
    // Magnitude is stored in the first half of the original buffer
    void MagnitudeFFT(float *samples, int len)
    {
        // square DC and Nyquist
        int nyquist = len>>1;
        samples[0] = sqrtf(samples[0] * samples[0]);
        samples[nyquist] = sqrtf(samples[nyquist] * samples[nyquist]);
        
        for (int i=1; i<nyquist; i++)
        {
            // Complex conjugate - cross product cancels
            samples[i] = samples[i] * samples[i];
            samples[len-i] = samples[len-i] * samples[len-i];
            samples[i] += samples[len-i];
            samples[i] = sqrtf(samples[i]);
        }
        return;
    }

    void MagnitudeFFTVec(std::vector<float> &samples)
    {
        int len = (int)samples.size();
        // square DC and Nyquist
        int nyquist = len>>1;
        samples[0] = sqrtf(samples[0] * samples[0]);
        samples[nyquist] = sqrtf(samples[nyquist] * samples[nyquist]);
        
        for (int i=1; i<nyquist; i++)
        {
            // Complex conjugate - cross product cancels
            samples[i] = samples[i] * samples[i];
            samples[len-i] = samples[len-i] * samples[len-i];
            samples[i] += samples[len-i];
            samples[i] = sqrtf(samples[i]);
        }
        return;
    }

	
} // namespace UtilityFuncs
