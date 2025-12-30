//
//  Windows.cpp
//  AdditiveSynthChange
//
//  Created by Riley on 1/30/25.
//
#include "Windows.h"
#include <vector>
using namespace std;

TransitionWindows::TransitionWindows(int longSize, int shortSize) : mLongSize(longSize), mShortSize(shortSize) {}
NormalWindows::NormalWindows() {}

vector<float> NormalWindows::HanningWindow(int length) {
    vector<float> hanning(length, 0.0f);
    for (int i = 0; i < length; i++) {
        hanning[i] = 0.5f * (1.0f - cos((2.0 * M_PI * i) / (length)));
    }
    return hanning;
}

vector<float> NormalWindows::RectangularWindow(int length) {
    vector<float> rect(length, 1.0f);
    return rect;
}

vector<float> TransitionWindows::createLongToShortWindow(const vector<float>& window_long, const vector<float>& window_short) {
    vector<float> long_to_short(mLongSize, 0.0f);
    
    for (int i = 0; i < mLongSize / 2; i++) {
        long_to_short[i] = window_long[i];
    }
    int length_ones_and_zeros = mLongSize/4 - mShortSize/4;
    
    int new_spot = mLongSize / 2;
    for (int i = new_spot; i < new_spot+length_ones_and_zeros; i++) {
        long_to_short[i] = 1.0f;
    }
    new_spot = new_spot+length_ones_and_zeros;
    for (int i = 0; i < mShortSize / 2; i++) {
        long_to_short[new_spot + i] = window_short[mShortSize / 2 + i];
    }
    return long_to_short;
}

vector<float> TransitionWindows::createShortToLongWindow(const vector<float>& long_to_short) {
    vector<float> short_to_long(mLongSize);
    
    for (int i = 0; i < mLongSize; i++) {
        short_to_long[i] = long_to_short[mLongSize - 1 - i];
    }
    
    return short_to_long;
}

vector<float> TransitionWindows::RectToHann(int length) {
    vector<float> rect_fade_to_hann(length, 0.0f);
    for (int i = 0; i < length/2; i++) {
        rect_fade_to_hann[i] = 1.0f;
    }
    NormalWindows winTemp;
    vector<float> hann = winTemp.HanningWindow(length);
    for (int i = length/2; i < length; i++) {
        rect_fade_to_hann[i] = hann[i];
    }
    
    return rect_fade_to_hann;
}
