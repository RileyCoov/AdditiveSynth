//
//  AnalysisInfo.h
//  AdditiveSynthChange
//
//  Created by Riley on 1/30/25.
//
#include <complex>
#include <vector>
#include "RealFFT.h"
#include "UtilityFuncs.h"
using namespace UtilityFuncs;
using namespace std;

class SynthInformation {
public:
    int start;
    int stop;
    int size;
    vector<float> windowApplied;
    bool trans;
    int hop_size;
    
    SynthInformation(int _start, int _stop, int _size, vector<float> _windowApplied, bool _trans, int _hop_size) : start(_start), stop(_stop), size(_size), windowApplied(_windowApplied), trans(_trans), hop_size(_hop_size) {}
};


class STFTAdjustment {
public:
    STFTAdjustment(int longSize, int shortSize, int hopSize);
    
    void processFrames(int num_frames, const vector<float> &audioData, const vector<float> &transientList, int hopSize, const vector<float> &longWindow, const vector<float> &shortWindow, const vector<float> &longToShortWindow, const vector<float> &shortToLongWindow);
    
    const vector<vector<complex<double>>>& getAdjustedSTFT() {return m_readjustedSTFT;}
    const vector<SynthInformation>& getSynthPlacement() {return m_synthPlacement;}
    
private:
    void applyAndProcessWindow(const vector<float>& signal, int start, const vector<float>& window, int size, vector<vector<complex<double>>>& output, int frameNum) {
        vector<float> frame_data(size);
        for (int i = 0; i < size; i++) {
            frame_data[i] = signal[start + i] * window[i];
        }

        RealFFT(frame_data.data(), size);
        vector<complex<double>> fft_result(size / 2);
        fft_result[0] = complex<double>(frame_data[0], 0.0);
        fft_result[(size/2)-1] = complex<double>(frame_data[(size/2)-1], 0.0);
        for (int k = 1; k < size / 2; k++) {
            double re = frame_data[k];
            double im = frame_data[size - k];
            fft_result[k] = complex<double>(re, im);
        }
        output.push_back(fft_result);
    }

    vector<vector<complex<double>>> m_readjustedSTFT;
    vector<SynthInformation> m_synthPlacement;
    int m_longSize;
    int m_shortSize;
    int m_hopSize;
    int m_movementBuffer;
    bool m_appliedShort;
};


//This is the class that will walk through all our frames of ffts, and pull out
//What we want for additive
/*
 class HardCoreAnalysis {
     
 }
 
 */

