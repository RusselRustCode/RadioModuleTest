#pragma once

#define _USE_MATH_DEFINES
#include <vector>
#include <cmath>
#include <cstdint>
#include <tuple>
#include <string>
#include <fstream>
#include <algorithm>

constexpr int16_t NUM_LEVELS = 32767;
constexpr float   SAMPLE_RATE = 100.0f;

struct SignalData {
    std::vector<float>   original;
    std::vector<int16_t> quantized;
    std::vector<float>   interp_float;
    std::vector<int16_t> interp_fixed;
    float quant_error;
    float freq;
};

std::vector<float> create_sin(float ampl, float freq, float phi,
                               size_t num_samples,
                               float sample_rate = SAMPLE_RATE)
{
    std::vector<float> out(num_samples);
    float factor = 2.0f * M_PI * freq / sample_rate;
    for (size_t n = 0; n < num_samples; ++n)
        out[n] = ampl * std::sin(factor * n + phi);
    return out;
}

std::tuple<std::vector<int16_t>, float> quantize_16(const std::vector<float>& signal,
                                                     float ampl)
{
    std::vector<int16_t> out(signal.size());
    float error = 0.0f;
    for (size_t i = 0; i < signal.size(); ++i) {
        int32_t value = std::round(signal[i] * NUM_LEVELS / ampl);
        out[i] = static_cast<int16_t>(std::clamp(value, (int32_t)-32768, (int32_t)32767));
        error += std::abs(signal[i] - out[i] * ampl / NUM_LEVELS);
    }
    error = std::sqrt(error / signal.size());
    return {out, error};
}

std::vector<float> linear_interpolation(const std::vector<float>& signal)
{
    std::vector<float> out(2 * signal.size());
    for (size_t i = 0; i < signal.size(); ++i)
        out[2 * i] = signal[i];
    for (size_t i = 0; i < signal.size(); ++i) {
        if (i != signal.size() - 1)
            out[2 * i + 1] = (signal[i] + signal[i + 1]) / 2.0f;
        else
            out[2 * i + 1] = signal[i];
    }
    return out;
}

std::vector<float> catmull_rom_spline(const std::vector<float>& signal)
{
    std::vector<float> out(2 * signal.size());
    for (size_t i = 0; i < signal.size(); ++i) {
        out[2 * i] = signal[i];
        if (i == 0 || i >= signal.size() - 2) {
            float x2 = (i + 1 < signal.size()) ? signal[i + 1] : signal[i];
            out[2 * i + 1] = (signal[i] + x2) / 2.0f;
        } else {
            out[2 * i + 1] = (-signal[i-1] + 9*signal[i] + 9*signal[i+1] - signal[i+2]) / 16.0f;
        }
    }
    return out;
}

std::vector<int16_t> interpolate_qt_linear(const std::vector<int16_t>& signal_qt)
{
    std::vector<int16_t> out(2 * signal_qt.size());
    for (size_t i = 0; i < signal_qt.size(); ++i) {
        out[2 * i] = signal_qt[i];
        if (i == signal_qt.size() - 1)
            out[2 * i + 1] = signal_qt[i];
        else
            out[2 * i + 1] = (signal_qt[i] + signal_qt[i + 1]) / 2;
    }
    return out;
}

void write_csv(const SignalData& data)
{
    std::ofstream orig_file("orig_f" + std::to_string((int)data.freq) + "Hz.csv");
    orig_file << "original\n";
    for (auto v : data.original)
        orig_file << v << "\n";

    std::ofstream interp_file("signal_f" + std::to_string((int)data.freq) + "Hz.csv");
    interp_file << "interp_float,interp_fixed\n";
    for (size_t i = 0; i < data.interp_float.size(); ++i)
        interp_file << data.interp_float[i] << "," << data.interp_fixed[i] << "\n";
}

SignalData process_frequency(float freq, float ampl, size_t num_samples)
{
    SignalData data;
    data.freq        = freq;
    data.original    = create_sin(ampl, freq, 0.0f, num_samples);
    auto [q, err]    = quantize_16(data.original, ampl);
    data.quantized   = q;
    data.quant_error = err;
    data.interp_float = catmull_rom_spline(data.original);
    data.interp_fixed = interpolate_qt_linear(data.quantized);
    return data;
}
