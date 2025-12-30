//
//  AnalysisInfo.cpp
//  AdditiveSynthChange
//
//  Created by Riley on 1/30/25.
//
#include <complex>
#include "AnalysisInfo.h"
#include <vector>
#include "RealFFT.h"
#include "UtilityFuncs.h"
using namespace UtilityFuncs;
using namespace std;


STFTAdjustment::STFTAdjustment(int longSize, int shortSize, int hopSize) : m_longSize(longSize), m_shortSize(shortSize), m_hopSize(hopSize) {};

void STFTAdjustment::processFrames(int num_frames, const vector<float> &audioData, const vector<float> &transientList, int hopSize, const vector<float> &longWindow, const vector<float> &shortWindow, const vector<float> &longToShortWindow, const vector<float> &shortToLongWindow) {
    int movementBuffer = (m_longSize/2) + (m_longSize/4 - m_shortSize/4);
    bool appliedShort = false;
    //ERROR HERE
    for (int f = 0; f < num_frames; f++) {
        if (f < num_frames - 1 && (transientList[f] != 0.0f || transientList[f+1] != 0.0f)) {
            // Current frame or next frame has a transient
            bool currTrans = (transientList[f] != 0.0f);
            bool nextTrans = (f < num_frames - 1 && transientList[f+1] != 0.0f);
            if (!appliedShort){
                if (transientList[f] == 0.0f && transientList[f+1] != 0.0f) {
                    // Transition from long to short
                    applyAndProcessWindow(audioData, f * hopSize, longToShortWindow, m_longSize, m_readjustedSTFT,f);
                    SynthInformation hold = SynthInformation(f*hopSize, (f*hopSize)+m_longSize, m_longSize, longToShortWindow, true, movementBuffer); //long_to_short
                    m_synthPlacement.push_back(hold);
                    appliedShort = true;
                }
            } else {
                // Process short windows
                int repetition = m_longSize / m_shortSize;
                for (int i = 0; i < repetition; i++) { // i * SHORT_SIZE / 2 +
                    int movement = i * m_shortSize/2;
                    applyAndProcessWindow(audioData, (f-1) * hopSize +  movementBuffer + movement, shortWindow, m_shortSize, m_readjustedSTFT, f);
                    // i * SHORT_SIZE / 2 +
                    SynthInformation hold = SynthInformation((f-1) * hopSize + movementBuffer + movement, (f-1) * hopSize + movementBuffer+m_shortSize+movement, m_shortSize, shortWindow, false, m_shortSize/2);
                    m_synthPlacement.push_back(hold);
                }
                appliedShort = true;
            }
            
        } else {
            // No transient, use long window
            if (!appliedShort) {
                applyAndProcessWindow(audioData, f * hopSize, longWindow, m_longSize, m_readjustedSTFT, f);
                SynthInformation hold = SynthInformation(f*hopSize, (f*hopSize)+m_longSize, m_longSize, longWindow, false, m_longSize/2);
                m_synthPlacement.push_back(hold);
            } else {
                applyAndProcessWindow(audioData, f * hopSize, shortToLongWindow, m_longSize, m_readjustedSTFT, f);
                SynthInformation hold = SynthInformation(f*hopSize, (f*hopSize)+m_longSize, m_longSize, shortToLongWindow, true, m_shortSize/2); //short_to_long
                m_synthPlacement.push_back(hold);
            }
            appliedShort = false;
        }
    }
}
