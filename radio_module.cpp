#define _USE_MATH_DEFINES
#include <iostream>
#include <vector>
#include <cmath>
#include <cstdint>
#include <tuple>
#include <string>
#include <fstream>
#include <algorithm>

constexpr int16_t NUM_LEVELS = 32767;

using std::cout;
using std::endl;


std::vector<float> create_sin(float ampl, float freq, float phi, size_t num_samples){
    size_t sample_rate = 100;
    std::vector<float> out(num_samples);
    float factor = static_cast<float>(2.0f * M_PI * freq / sample_rate);
    for(size_t n = 0; n < num_samples; ++n){
        out[n] = ampl * std::sin(factor * n + phi);
    }
    return out;
}


std::tuple<std::vector<int16_t>, float> quantize_16(const std::vector<float>& signal, float ampl){
    std::vector<int16_t> out(signal.size());
    float error = 0.0f;
    for(size_t i = 0; i < signal.size(); ++i){
        int32_t value = std::round(signal[i] * NUM_LEVELS / ampl);
        out[i] = static_cast<int16_t>(std::clamp(value, (int32_t)-32768, (int32_t)32767));
        error += std::abs(static_cast<float>(signal[i] - out[i] * ampl / NUM_LEVELS));
    }
    return {out, error};
}


std::vector<float> linear_interpolation(const std::vector<float>& signal){
    std::vector<float> out(2 * signal.size());
    for(size_t i = 0; i < signal.size(); ++i){
        out[2 * i] = signal[i];
    }
    for(size_t i = 0; i < signal.size(); ++i){
        if(i != signal.size() - 1){
            out[2 * i + 1] = (signal[i] + signal[i + 1]) / 2;
        }
        else{
            out[2 * i + 1] = signal[i];
        }
    }
    return out;
}


std::vector<float> catmull_rom_spline(const std::vector<float>& signal){
    std::vector<float> out(2 * signal.size());
    for(size_t i = 0; i < signal.size(); ++i){
        out[2 * i] = signal[i];
        if(i == 0 || i >= signal.size() - 2){
            float x2 = (i + 1 < signal.size()) ? signal[i + 1] : signal[i];
            out[2 * i + 1] = (signal[i] + x2) / 2.0f;
        }
        else{
            out[2 * i + 1] = (-signal[i - 1] + 9 * signal[i] + 9 * signal[i + 1] - signal[i + 2]) / 16.0f; 
        }
    }
    return out;
}


std::vector<int16_t> interpolate_qt_linear(const std::vector<int16_t>& signal_qt) {
    std::vector<int16_t> out(2 * signal_qt.size());
    for(size_t i = 0; i < signal_qt.size(); i++){
        out[2 * i] = signal_qt[i];
        if(i == signal_qt.size() - 1){
            out[2 * i + 1] = (signal_qt[i] + signal_qt[i]) / 2;
        }
        else{
            out[2 * i + 1] = (signal_qt[i] + signal_qt[i + 1]) / 2;
        }
    }
    return out;
}

void write_csv(float freq, float freq_sample, const std::vector<float>& signal, const std::vector<float>& interpolate_signal, const std::vector<int16_t>& fixed_interpolate_signal){
    std::string filename = "signal_f" + std::to_string((int)freq) + "Hz.csv";
    std::ofstream file(filename);
    file << "original,interp_float,interp_fixed\n";
    for(size_t i = 0; i < interpolate_signal.size(); ++i){
        float orig = (i % 2 == 0) ? signal[i / 2]: 0.0f;
        file << orig << "," << interpolate_signal[i] << "," << fixed_interpolate_signal[i] << "\n";
    }
}


int main(){
    size_t num_levels = std::pow(2, sizeof(int16_t) * 8) / 2 - 1;
    cout << num_levels << " " << 2.5 * 4 << endl;
    auto vec = create_sin(8.0f, 45, 0.0f, 30);
    auto [out, error] = quantize_16(vec, 8.0f);
    cout << "Error: " << error << endl;
    auto res = linear_interpolation(vec);
    cout << res.size() << endl;
    return 0;
}