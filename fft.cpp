#define _USE_MATH_DEFINES
#include <iostream>
#include <vector>
#include <complex> 
#include <cmath>
#include <random>
#include <exception>

using samples = std::complex<float>;
using samplesBuff = std::vector<samples>;
using std::cout;
using std::endl;


class fastFourier {
private:
    std::vector<samples> input;
    

public:
    explicit fastFourier(std::vector<samples> data) : input(std::move(data)) {}

    std::vector<samples> dft(){
        if(input.empty()){
            cout << "Вектор пустой" << endl;
            return {};
        }
        samplesBuff out(input.size());

        for(size_t k = 0; k < input.size(); ++k){
            for(size_t i = 0; i < input.size(); ++i){
                float twiddle = 2 * M_PI * k * i / input.size();
                out[k] += input[i] * samples(std::cos(twiddle), -std::sin(twiddle));
            }
        }
        return out;
        
    }

    std::vector<samples> fft(const std::vector<samples>& input){
        if(input.empty()){
           throw std::invalid_argument("fft: fft vector is empty!");
        }

        if(!is_smooth(input.size())){
            throw std::invalid_argument("fft: fft vector size is not 5 smooth: " + std::to_string(input.size()));
        }

        if(input.size() <= 1) return input;
        size_t radix;
        if(input.size() % 2 == 0) radix = 2;
        else if(input.size() % 3 == 0) radix = 3;
        else if(input.size() % 5 == 0) radix = 5;
        else return {};

        size_t new_size = input.size() / radix;
        samplesBuff out(input.size());
        std::vector<samplesBuff> sub(radix, samplesBuff(new_size, samples(0.0f, 0.0f)));
        for(size_t r = 0; r < radix; ++r){
            for(size_t i = 0; i < new_size; ++i){
                sub[r][i] = input[i * radix + r];
            }
        }

        for(size_t r = 0; r < radix; ++r){
            sub[r] = fft(sub[r]);
        }

        if(radix == 2) out = butterfly2(sub, input.size());
        if(radix == 3) out = butterfly3(sub, input.size());
        if(radix == 5) out = butterfly5(sub, input.size());
        
        
        return out;    
    };


    std::vector<samples> ifft(const std::vector<samples>& input){
        if(input.empty()){
            throw std::invalid_argument("ifft: ifft vector is empty!");
        };
        if(!is_smooth(input.size())){
            throw std::invalid_argument("fft: fft vector size is not 5 smooth: " + std::to_string(input.size()));
        }

        auto result = ifft_recursive(input);
        float size = static_cast<float>(input.size());
        for(auto& el: result){
            el /= size;
        }
        return result;
    };

private:
    bool is_smooth(size_t num_of_samples){
        if(num_of_samples == 0) return 0;
        while(num_of_samples % 2 == 0) num_of_samples /= 2;
        while(num_of_samples % 3 == 0) num_of_samples /= 3;
        while(num_of_samples % 5 == 0) num_of_samples /= 5;
        return num_of_samples == 1;
    }

    std::vector<samples> butterfly2(std::vector<samplesBuff>& input, size_t N){
        size_t new_size = N / 2;
        samplesBuff out(N);
        for(size_t k = 0; k < new_size; ++k){
            samples twiddleFactor = samples(
                std::cos(2 * M_PI * k / N),
                -std::sin(2 * M_PI * k / N)
            );

            out[k] = input[0][k] + twiddleFactor * input[1][k];
            out[k + new_size] = input[0][k] - twiddleFactor * input[1][k];
        }
        return out;        
    }

    std::vector<samples> butterfly3(std::vector<samplesBuff>& input, size_t N){
        size_t new_size = N / 3;
        samplesBuff out(N);
        samples w = samples(
            std::cos(2 * M_PI / 3),
            -std::sin(2 * M_PI / 3)
        );
        for(size_t k = 0; k < new_size; k++){
            samples twiddleFactor1 = samples(
                std::cos(2 * M_PI * k / N),
                -std::sin(2 * M_PI * k / N)
            );

            samples twiddleFactor2 = samples(
                std::cos(4 * M_PI * k / N),
                -std::sin(4 * M_PI * k / N)
            );

            out[k] = input[0][k] + twiddleFactor1 * input[1][k] + twiddleFactor2 * input[2][k];
            out[k + new_size] = input[0][k] + w * twiddleFactor1 * input[1][k] + w * w * twiddleFactor2 * input[2][k];
            out[k + 2 * new_size] = input[0][k] + w * w * twiddleFactor1 * input[1][k] + w * twiddleFactor2 * input[2][k];
        }
        return out;
    }

    std::vector<samples> butterfly5(std::vector<samplesBuff>& input, size_t N){
        size_t new_size = N / 5;
        samplesBuff out(N);
        samples w = samples(
            std::cos(2 * M_PI / 5),
            -std::sin(2 * M_PI / 5)
        );

        samples w_2 = w*w;
        samples w_3= w_2 * w;
        samples w_4 = w_2*w_2;
        samples w_8 = w_4*w_4;
        samples w_16 = w_8*w_8;

        for(size_t k = 0; k < new_size; k++){
            samples twiddleFactor1 = samples(
                std::cos(2 * M_PI * k / N),
                -std::sin(2 * M_PI * k / N)
            );
    
            samples twiddleFactor2 = samples(
                std::cos(4 * M_PI * k / N),
                -std::sin(4 * M_PI * k / N)
            );
    
            samples twiddleFactor3 = samples(
                std::cos(6 * M_PI * k / N),
                -std::sin(6 * M_PI * k / N)
            );
    
            samples twiddleFactor4 = samples(
                std::cos(8 * M_PI * k / N),
                -std::sin(8 * M_PI * k / N)
            );

            out[k] = input[0][k] + twiddleFactor1 * input[1][k] + twiddleFactor2 * input[2][k] + twiddleFactor3 * input[3][k] 
            + twiddleFactor4 * input[4][k];
            
            out[k + new_size] = input[0][k] + w * twiddleFactor1 * input[1][k] + w_2 * twiddleFactor2 * input[2][k] + w_3 * twiddleFactor3 * input[3][k]
            + w_4 * twiddleFactor4 * input[4][k];
            
            out[k + 2 * new_size] = input[0][k] + w_2 * twiddleFactor1 * input[1][k] + w_4 * twiddleFactor2 * input[2][k] + w * twiddleFactor3 * input[3][k]
            + w_3 * twiddleFactor4 * input[4][k];
            
            out[k + 3 * new_size] = input[0][k] + w_3 * twiddleFactor1 * input[1][k] + w * twiddleFactor2 * input[2][k] + w_4 * twiddleFactor3 * input[3][k]
            + w_2 * twiddleFactor4 * input[4][k]; 
             
            out[k + 4 * new_size] = input[0][k] + w_4 * twiddleFactor1 * input[1][k] + w_3 * twiddleFactor2 * input[2][k] + w_2 * twiddleFactor3 * input[3][k]
            + w * twiddleFactor4 * input[4][k];
        }
        return out;
    }

    std::vector<samples> ifft_recursive(const std::vector<samples>& input){
        if(input.size() <= 1) return input;
        size_t radix;
        if(input.size() % 2 == 0) radix = 2;
        else if(input.size() % 3 == 0) radix = 3;
        else if(input.size() % 5 == 0) radix = 5;
        else return {};

        size_t new_size = input.size() / radix;
        samplesBuff out(input.size());
        std::vector<samplesBuff> sub(radix, samplesBuff(new_size, samples(0.0f, 0.0f)));
        for(size_t r = 0; r < radix; ++r){
            for(size_t i = 0; i < new_size; ++i){
                sub[r][i] = input[i * radix + r];
            }
        }

        for(size_t r = 0; r < radix; ++r){
            sub[r] = ifft_recursive(sub[r]);
        }

        if(radix == 2) out = inverse_butterfly2(sub, input.size());
        if(radix == 3) out = inverse_butterfly3(sub, input.size());
        if(radix == 5) out = inverse_butterfly5(sub, input.size());
        
        
        return out;    
    };

    std::vector<samples> inverse_butterfly2(std::vector<samplesBuff>& input, size_t N){
        size_t new_size = N / 2;
        samplesBuff out(N);
        for(size_t k = 0; k < new_size; ++k){
            samples twiddleFactor = samples(
                std::cos(2 * M_PI * k / N),
                std::sin(2 * M_PI * k / N)
            );

            out[k] = input[0][k] + twiddleFactor * input[1][k];
            out[k + new_size] = input[0][k] - twiddleFactor * input[1][k];
        }
        return out;        
    }

    std::vector<samples> inverse_butterfly3(std::vector<samplesBuff>& input, size_t N){
        size_t new_size = N / 3;
        samplesBuff out(N);
        samples w = samples(
            std::cos(2 * M_PI / 3),
            std::sin(2 * M_PI / 3)
        );
        for(size_t k = 0; k < new_size; k++){
            samples twiddleFactor1 = samples(
                std::cos(2 * M_PI * k / N),
                std::sin(2 * M_PI * k / N)
            );

            samples twiddleFactor2 = samples(
                std::cos(4 * M_PI * k / N),
                std::sin(4 * M_PI * k / N)
            );

            out[k] = input[0][k] + twiddleFactor1 * input[1][k] + twiddleFactor2 * input[2][k];
            out[k + new_size] = input[0][k] + w * twiddleFactor1 * input[1][k] + w * w * twiddleFactor2 * input[2][k];
            out[k + 2 * new_size] = input[0][k] + w * w * twiddleFactor1 * input[1][k] + w * twiddleFactor2 * input[2][k];
        }
        return out;
    }

    std::vector<samples> inverse_butterfly5(std::vector<samplesBuff>& input, size_t N){
        size_t new_size = N / 5;
        samplesBuff out(N);
        samples w = samples(
            std::cos(2 * M_PI / 5),
            std::sin(2 * M_PI / 5)
        );

        samples w_2 = w*w;
        samples w_3= w_2 * w;
        samples w_4 = w_2*w_2;
        samples w_8 = w_4*w_4;
        samples w_16 = w_8*w_8;

        for(size_t k = 0; k < new_size; k++){
            samples twiddleFactor1 = samples(
                std::cos(2 * M_PI * k / N),
                std::sin(2 * M_PI * k / N)
            );
    
            samples twiddleFactor2 = samples(
                std::cos(4 * M_PI * k / N),
                std::sin(4 * M_PI * k / N)
            );
    
            samples twiddleFactor3 = samples(
                std::cos(6 * M_PI * k / N),
                std::sin(6 * M_PI * k / N)
            );
    
            samples twiddleFactor4 = samples(
                std::cos(8 * M_PI * k / N),
                std::sin(8 * M_PI * k / N)
            );

            out[k] = input[0][k] + twiddleFactor1 * input[1][k] + twiddleFactor2 * input[2][k] + twiddleFactor3 * input[3][k] 
            + twiddleFactor4 * input[4][k];
            
            out[k + new_size] = input[0][k] + w * twiddleFactor1 * input[1][k] + w_2 * twiddleFactor2 * input[2][k] + w_3 * twiddleFactor3 * input[3][k]
            + w_4 * twiddleFactor4 * input[4][k];
            
            out[k + 2 * new_size] = input[0][k] + w_2 * twiddleFactor1 * input[1][k] + w_4 * twiddleFactor2 * input[2][k] + w * twiddleFactor3 * input[3][k]
            + w_3 * twiddleFactor4 * input[4][k];
            
            out[k + 3 * new_size] = input[0][k] + w_3 * twiddleFactor1 * input[1][k] + w * twiddleFactor2 * input[2][k] + w_4 * twiddleFactor3 * input[3][k]
            + w_2 * twiddleFactor4 * input[4][k]; 
             
            out[k + 4 * new_size] = input[0][k] + w_4 * twiddleFactor1 * input[1][k] + w_3 * twiddleFactor2 * input[2][k] + w_2 * twiddleFactor3 * input[3][k]
            + w * twiddleFactor4 * input[4][k];
        }
        return out;
    }
    
};




int main(){
    try{
        std::random_device rd;
        std::mt19937 generator(rd());
        size_t size = 60;
        std::uniform_real_distribution<float> distr (-100.0f, 100.0f);
        std::vector<samples> input(size);
        for(size_t i = 0; i < input.size(); ++i){
            float real = distr(generator);
            float imag = distr(generator);
            input[i] = samples(real, imag);
        }

        for(size_t i = 0; i < 10; i++){
            cout << input[i].real() << " " << input[i].imag() << endl;
        }
        cout<<endl;
        cout<<"Test FFT"<<endl;

        fastFourier ff(input);
        auto result = ff.fft(input);
        for(size_t i = 0; i < 15; ++i){
            cout << result[i] << endl;
        }
        cout << endl;
        cout << "Test Inv FFT" << endl;
        auto inv_res = ff.ifft(result);
        for(auto& el: inv_res){
            cout << el << endl;
        }
        float error = 0.0f;
        for(size_t i = 0; i < size; ++i){
            error += std::abs(input[i] - inv_res[i]);
        }
        cout << endl;
        cout << "Error: " << error << endl;
    }
    catch(const std::invalid_argument& e){
        cout << "Error: " << e.what() << endl; 
    }
    catch(const std::exception& e){
        cout << "Not expected error: " << e.what() << endl;
    }
    return 0;
}
