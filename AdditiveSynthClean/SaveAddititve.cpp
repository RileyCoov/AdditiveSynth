//
//  SaveAddititve.cpp
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

#include "RealFFT.h"
#include "AnalysisInfo.h"
#include "Phaser.h"
#include "SampleNormalizer.h"
#include "Windows.h"
#include "UtilityFuncs.h"
#include "SaveAddititve.h"
using namespace UtilityFuncs;
using namespace std;

void save_binary(const string& filename, int num_frames, vector<SynthInformation> containsSynthPlacement, vector<vector<PeakTrack>> frames_peaks) {
    //maintain a filebuf object as internal stream buffer
    //Performing input/output operations on the file
    std::ofstream out(filename, std::ios::binary);
    for (int frame_idx = 0; frame_idx < num_frames; frame_idx++) {
        SynthInformation& info = containsSynthPlacement[frame_idx];
        
        BinaryFrameHeader header{ //THis is a struct type which gives me the size nicely
            .frame_index = frame_idx,
            .start = info.start,
            .stop = info.stop,
            .frame_size = info.size,
            .trans = static_cast<uint8_t>(info.trans),
            .hop_size = info.hop_size
        };
        //Just header so literally has the info on containsSynthPlace
        out.write(reinterpret_cast<const char*>(&header), sizeof(header));
        
        //*********NOTE
        // The reason for the const char* casting because of this being in
        //write mode. This allows for the first part to be a pointer to memeory
        //where we are keeping our data
        
        out.write(
                  reinterpret_cast<const char*>(info.windowApplied.data()),
                  info.windowApplied.size() * sizeof(float)
                  ); //Spill the contents of windowApplied
        const auto& peaks = frames_peaks[frame_idx];
        size_t num_peaks = peaks.size();
        out.write(reinterpret_cast<const char*>(&num_peaks), sizeof(size_t));
        
        for (const PeakTrack& peak : peaks) {
            out.write(reinterpret_cast<const char*>(&peak.id), sizeof(int));
            out.write(reinterpret_cast<const char*>(&peak.freq_hz), sizeof(double));
            out.write(reinterpret_cast<const char*>(&peak.max_db), sizeof(double));
            out.write(reinterpret_cast<const char*>(&peak.current_db), sizeof(double));
            out.write(reinterpret_cast<const char*>(&peak.peak_bin), sizeof(int));
            out.write(reinterpret_cast<const char*>(&peak.phase), sizeof(double));
            
            uint8_t alive = static_cast<uint8_t>(peak.alive);
            uint8_t edit = static_cast<uint8_t>(peak.edit);
            out.write(reinterpret_cast<const char*>(&alive), sizeof(uint8_t));
            out.write(reinterpret_cast<const char*>(&edit), sizeof(uint8_t));
        }
        
        //At the end of this function we have the data in teh following order
        //HEADER | WINDOW DATA (Actual window in float) | num_peaks | Peak 1 | Peak 2 ...
        //In total: | Header_0 | window_0 | num_peaks_0 | Peaks_0 | Header_1 | window_1 | num_peaks_1 | Peaks_1 | ... |
    }
}



void read_binary(const std::string& filename, int num_frames,
                 std::vector<SynthInformation>& containsSynthPlacement,
                 std::vector<std::vector<PeakTrack>>& frames_peaks)
{
    std::ifstream in(filename, std::ios::binary);
    if (!in)
    {
        std::cerr << "Error opening file " << filename << std::endl;
        return;
    }

    containsSynthPlacement.clear();
    frames_peaks.clear();

    for (int frame_idx = 0; frame_idx < num_frames; frame_idx++)
    {
        // Read the header.
        BinaryFrameHeader header;
        in.read(reinterpret_cast<char*>(&header), sizeof(header));
        if (!in)
        {
            std::cerr << "Error reading header for frame " << frame_idx << std::endl;
            break;
        }
        
        SynthInformation info(0, 0, 0, {}, false, 0);
        info.start = header.start;
        info.stop = header.stop;
        info.size = header.frame_size;
        info.trans = static_cast<bool>(header.trans);
        info.hop_size = header.hop_size;

        info.windowApplied.resize(header.frame_size);
        in.read(reinterpret_cast<char*>(info.windowApplied.data()),
                header.frame_size * sizeof(float));
        if (!in)
        {
            std::cerr << "Error reading window data for frame " << frame_idx << std::endl;
            break;
        }

        containsSynthPlacement.push_back(info);
        
        // Read the number of peaks.
        size_t num_peaks;
        in.read(reinterpret_cast<char*>(&num_peaks), sizeof(size_t));
        if (!in)
        {
            std::cerr << "Error reading number of peaks for frame " << frame_idx << std::endl;
            break;
        }
        
        // Read each peak and store.
        std::vector<PeakTrack> peaks;
        for (size_t i = 0; i < num_peaks; i++)
        {
            PeakTrack peak(0, 0.0, 0.0, 0, 0.0);
            in.read(reinterpret_cast<char*>(&peak.id), sizeof(int));
            in.read(reinterpret_cast<char*>(&peak.freq_hz), sizeof(double));
            in.read(reinterpret_cast<char*>(&peak.max_db), sizeof(double));
            in.read(reinterpret_cast<char*>(&peak.current_db), sizeof(double));
            in.read(reinterpret_cast<char*>(&peak.peak_bin), sizeof(int));
            in.read(reinterpret_cast<char*>(&peak.phase), sizeof(double));
            
            uint8_t alive, edit;
            in.read(reinterpret_cast<char*>(&alive), sizeof(uint8_t));
            in.read(reinterpret_cast<char*>(&edit), sizeof(uint8_t));
            peak.alive = static_cast<bool>(alive);
            peak.edit = static_cast<bool>(edit);
            
            if (!in)
            {
                std::cerr << "Error reading peak " << i << " for frame " << frame_idx << std::endl;
                break;
            }
            
            peaks.push_back(peak);
        }
        
        frames_peaks.push_back(peaks);
    }
    
    in.close();
}
