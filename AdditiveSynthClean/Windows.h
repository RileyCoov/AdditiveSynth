//
//  Windows.h
//  AdditiveSynthChange
//
//  Created by Riley on 1/30/25.
//
#include <iostream>
#include <vector>
using namespace std;

class TransitionWindows {
public:
    TransitionWindows(int longSize, int shortSize);
    
    vector<float> createLongToShortWindow(const vector<float>& windowLong,
                                          const vector<float>& windowShort);
    
    vector<float> createShortToLongWindow(const vector<float>& longToShort);
    
    vector<float> RectToHann(int length);
    
    int getLongSize() {return mLongSize;}
    int getShortSize() {return mShortSize;}

private:
    const int mLongSize;
    const int mShortSize;
    
    int calculateTransitionLength() {
        return mLongSize/4 - mShortSize/4;
    }
};

class NormalWindows {
public:
    NormalWindows();
    
    vector<float> HanningWindow(int length);
    
    vector<float> RectangularWindow(int length);
    
};
