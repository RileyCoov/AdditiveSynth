//
//  SaveAddititve.h
//  AdditiveSynthChange
//
//  Created by Riley on 2/20/25.
//
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <complex>
#include <fstream> // This is for the binary process
using namespace std;

// The addition of pragma ensures there isn't padding happening which could
//cause more memory than wanted to be used
struct BinaryFrameHeader {
    int32_t frame_index;
    int32_t start;
    int32_t stop;
    int32_t frame_size;
    uint8_t trans;
    int32_t hop_size;
    size_t window_size;
};

class PeakTrack {
public:
    int id;
    double freq_hz;
    double max_db;
    double current_db;
    int peak_bin;
    double phase;
    bool alive;
    bool edit;
    
    PeakTrack(int _id, double _freq, double _mag, int _peak_bin, double _phase) : id(_id), freq_hz(_freq), max_db(_mag), current_db(_mag), peak_bin(_peak_bin), phase(_phase), alive(true), edit(true) {}
};


void save_binary(const std::string& filename, int num_frames,
                 std::vector<SynthInformation> containsSynthPlacement,
                 std::vector<std::vector<PeakTrack>> frames_peaks);
void read_binary(const std::string& filename, int num_frames,
                 std::vector<SynthInformation>& containsSynthPlacement,
                 std::vector<std::vector<PeakTrack>>& frames_peaks);

