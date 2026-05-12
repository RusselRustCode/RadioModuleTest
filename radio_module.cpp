#include <iostream>
#include <random>
#include <vector>
#include <complex>
#include <cmath>
#define _USE_MATH_DEFINES

using std::cout;
using std::endl;
using samples = std::complex<double>;



std::vector<double> create_sin(double ampl, double freq, double phi, size_t num_samples){
    float discret_freq = 1.0f / 100.0f;
    double factor = 2 * M_PI * freq;
    auto sin_func = [=](size_t num_samples) -> double { return ampl * std::sin(factor * num_samples * discret_freq + phi);};

    std::vector<double> output(num_samples);
    
    for(size_t n = 0; n < num_samples; ++n){
        output[n] = sin_func(n);
        
    }

    return output;
    
}

int main(){
    std::random_device rd;
    std::mt19937 gen(rd());

    auto vec = create_sin(5.0, 40, 0.0, 100);

    for(const auto& el: vec){
        cout << el << " "; 
    }
}
