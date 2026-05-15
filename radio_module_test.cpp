#include "radio_module.h"
#include <gtest/gtest.h>

// ===== Тесты create_sin =====

TEST(CreateSinTest, OutputSize) {
    auto sig = create_sin(1.0f, 5.0f, 0.0f, 100);
    EXPECT_EQ(sig.size(), 100u);
}

TEST(CreateSinTest, AmplitudeRespected) {
    float ampl = 3.5f;
    auto sig = create_sin(ampl, 5.0f, 0.0f, 1000);
    float max_val = *std::max_element(sig.begin(), sig.end());
    float min_val = *std::min_element(sig.begin(), sig.end());
    EXPECT_NEAR(max_val,  ampl, 0.01f);
    EXPECT_NEAR(min_val, -ampl, 0.01f);
}

TEST(CreateSinTest, ZeroFrequency) {
    auto sig = create_sin(1.0f, 0.0f, 0.0f, 50);
    for (float v : sig)
        EXPECT_NEAR(v, 0.0f, 1e-6f);
}

TEST(CreateSinTest, PhaseShift) {
    float ampl = 1.0f;
    auto sig = create_sin(ampl, 0.0f, M_PI / 2.0f, 10);
    EXPECT_NEAR(sig[0], ampl, 1e-6f);
}

TEST(CreateSinTest, EmptyOutput) {
    auto sig = create_sin(1.0f, 5.0f, 0.0f, 0);
    EXPECT_TRUE(sig.empty());
}

// ===== Тесты quantize_16 =====

TEST(Quantize16Test, OutputSize) {
    auto sig = create_sin(1.0f, 5.0f, 0.0f, 100);
    auto [q, err] = quantize_16(sig, 1.0f);
    EXPECT_EQ(q.size(), 100u);
}

TEST(Quantize16Test, ErrorIsNonNegative) {
    auto sig = create_sin(1.0f, 5.0f, 0.0f, 200);
    auto [q, err] = quantize_16(sig, 1.0f);
    EXPECT_GE(err, 0.0f);
}

TEST(Quantize16Test, ErrorIsSmall) {
    auto sig = create_sin(1.0f, 5.0f, 0.0f, 500);
    auto [q, err] = quantize_16(sig, 1.0f);
    EXPECT_LT(err, 1e-3f);
}

TEST(Quantize16Test, NeverExceedsInt16Range) {
    auto sig = create_sin(1.0f, 5.0f, 0.0f, 1000);
    auto [q, err] = quantize_16(sig, 1.0f);
    for (int16_t v : q) {
        EXPECT_GE(v, -32768);
        EXPECT_LE(v,  32767);
    }
}

TEST(Quantize16Test, DCSignalQuantization) {
    std::vector<float> dc(10, 1.0f);
    auto [q, err] = quantize_16(dc, 1.0f);
    for (int16_t v : q)
        EXPECT_EQ(v, NUM_LEVELS);
}

TEST(Quantize16Test, ClampOverflow) {
    std::vector<float> sig = {2.0f, -2.0f};
    auto [q, err] = quantize_16(sig, 1.0f);
    EXPECT_EQ(q[0], 32767);
    EXPECT_EQ(q[1], -32768);
}

// ===== Тесты linear_interpolation =====

TEST(LinearInterpTest, OutputSizeDoubled) {
    std::vector<float> sig = {1.0f, 2.0f, 3.0f, 4.0f};
    auto out = linear_interpolation(sig);
    EXPECT_EQ(out.size(), 2 * sig.size());
}

TEST(LinearInterpTest, OriginalSamplesPreserved) {
    std::vector<float> sig = {1.0f, 3.0f, 5.0f, 7.0f};
    auto out = linear_interpolation(sig);
    for (size_t i = 0; i < sig.size(); ++i)
        EXPECT_FLOAT_EQ(out[2 * i], sig[i]);
}

TEST(LinearInterpTest, MidpointIsAverage) {
    std::vector<float> sig = {0.0f, 4.0f, 8.0f, 12.0f};
    auto out = linear_interpolation(sig);
    EXPECT_FLOAT_EQ(out[1], 2.0f);
    EXPECT_FLOAT_EQ(out[3], 6.0f);
    EXPECT_FLOAT_EQ(out[5], 10.0f);
}

TEST(LinearInterpTest, LastMidpointEqualsLastSample) {
    std::vector<float> sig = {1.0f, 2.0f, 3.0f};
    auto out = linear_interpolation(sig);
    EXPECT_FLOAT_EQ(out[out.size() - 1], sig.back());
}

// ===== Тесты catmull_rom_spline =====

TEST(CatmullRomTest, OutputSizeDoubled) {
    std::vector<float> sig(10, 1.0f);
    auto out = catmull_rom_spline(sig);
    EXPECT_EQ(out.size(), 2 * sig.size());
}

TEST(CatmullRomTest, OriginalSamplesPreserved) {
    std::vector<float> sig = {1.0f, 2.0f, 4.0f, 3.0f, 5.0f};
    auto out = catmull_rom_spline(sig);
    for (size_t i = 0; i < sig.size(); ++i)
        EXPECT_FLOAT_EQ(out[2 * i], sig[i]);
}

TEST(CatmullRomTest, ConstantSignalUnchanged) {
    std::vector<float> sig(8, 3.0f);
    auto out = catmull_rom_spline(sig);
    for (float v : out)
        EXPECT_NEAR(v, 3.0f, 1e-5f);
}

TEST(CatmullRomTest, BetterThanLinearForSmoothSignal) {
    auto sig       = create_sin(1.0f, 5.0f, 0.0f, 100);
    auto reference = create_sin(1.0f, 5.0f, 0.0f, 200);

    auto linear = linear_interpolation(sig);
    auto spline = catmull_rom_spline(sig);

    float err_linear = 0.0f, err_spline = 0.0f;
    for (size_t i = 0; i < reference.size(); ++i) {
        err_linear += std::abs(reference[i] - linear[i]);
        err_spline += std::abs(reference[i] - spline[i]);
    }
    EXPECT_LT(err_spline, err_linear);
}

// ===== Тесты interpolate_qt_linear =====

TEST(InterpQtLinearTest, OutputSizeDoubled) {
    std::vector<int16_t> sig = {100, 200, 300, 400};
    auto out = interpolate_qt_linear(sig);
    EXPECT_EQ(out.size(), 2 * sig.size());
}

TEST(InterpQtLinearTest, OriginalSamplesPreserved) {
    std::vector<int16_t> sig = {100, 200, 300, 400};
    auto out = interpolate_qt_linear(sig);
    for (size_t i = 0; i < sig.size(); ++i)
        EXPECT_EQ(out[2 * i], sig[i]);
}

TEST(InterpQtLinearTest, MidpointIsAverageOfNeighbours) {
    std::vector<int16_t> sig = {0, 100, 200, 300};
    auto out = interpolate_qt_linear(sig);
    EXPECT_EQ(out[1], 50);
    EXPECT_EQ(out[3], 150);
    EXPECT_EQ(out[5], 250);
}

TEST(InterpQtLinearTest, LastMidpointEqualsLastSample) {
    std::vector<int16_t> sig = {10, 20, 30};
    auto out = interpolate_qt_linear(sig);
    EXPECT_EQ(out[out.size() - 1], sig.back());
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
