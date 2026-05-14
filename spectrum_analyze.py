import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def get_spectrum(signal, freq_samples):
    size = len(signal)
    spectrum = np.fft.rfft(signal)
    freqs = np.fft.rfftfreq(size, d=1.0/freq_samples)
    ampl = 2.0 * np.abs(spectrum) / size
    return freqs, ampl

def find_amplitude(freq, ampl, target_freq):
    idx = np.argmin(np.abs(freq - target_freq))
    return ampl[idx]

def analyze_frequency(freq, freq_samples_orig, freq_samples_interp, original, interpolate_float, interpolate_fixed):
    freq_orig, ampl_orig = get_spectrum(original, freq_samples_orig)
    freq_float, ampl_float = get_spectrum(interpolate_float, freq_samples_interp)
    freq_fixed, ampl_fixed = get_spectrum(interpolate_fixed, freq_samples_interp)
    
    A_orig  = find_amplitude(freq_orig,  ampl_orig,  freq)
    A_float = find_amplitude(freq_float, ampl_float, freq)
    A_fixed = find_amplitude(freq_fixed, ampl_fixed, freq)
    
    mirror_freq = freq_samples_orig - freq
    M_float = find_amplitude(freq_float, ampl_float, mirror_freq)
    M_fixed = find_amplitude(freq_fixed, ampl_fixed, mirror_freq)
    
    # затухание основного тона (дБ)
    atten_float = 20 * np.log10(A_float / A_orig + 1e-12)
    atten_fixed = 20 * np.log10(A_fixed / A_orig + 1e-12)

    # подавление зеркала (дБ)
    mirror_float = 20 * np.log10(M_float / A_orig + 1e-12)
    mirror_fixed = 20 * np.log10(M_fixed / A_orig + 1e-12)
    
    return {
        'A_orig':       A_orig,
        'A_float':      A_float,
        'A_fixed':      A_fixed,
        'atten_float':  atten_float,
        'atten_fixed':  atten_fixed,
        'mirror_float': mirror_float,
        'mirror_fixed': mirror_fixed,
    }
    
def plot_spectrum(f, fs_orig, fs_interp, original, interp_float, interp_fixed):
    freqs_orig,  amp_orig  = get_spectrum(original,     fs_orig)
    freqs_float, amp_float = get_spectrum(interp_float, fs_interp)
    freqs_fixed, amp_fixed = get_spectrum(interp_fixed, fs_interp)

    fig, axes = plt.subplots(3, 1, figsize=(12, 8))
    fig.suptitle(f'Спектр сигнала f = {f} Гц', fontsize=14)

    axes[0].stem(freqs_orig, amp_orig, markerfmt='C0o', linefmt='C0-', basefmt='k-')
    axes[0].set_title(f'Исходный сигнал (fs = {fs_orig} Гц)')
    axes[0].set_xlabel('Частота (Гц)')
    axes[0].set_ylabel('Амплитуда')
    axes[0].set_xlim(0, fs_interp / 2)
    axes[0].grid(True)

    axes[1].stem(freqs_float, amp_float, markerfmt='C1o', linefmt='C1-', basefmt='k-')
    axes[1].set_title(f'После интерполяции float (fs = {fs_interp} Гц)')
    axes[1].set_xlabel('Частота (Гц)')
    axes[1].set_ylabel('Амплитуда')
    axes[1].set_xlim(0, fs_interp / 2)
    axes[1].grid(True)

    axes[2].stem(freqs_fixed, amp_fixed, markerfmt='C2o', linefmt='C2-', basefmt='k-')
    axes[2].set_title(f'После интерполяции fixed-point (fs = {fs_interp} Гц)')
    axes[2].set_xlabel('Частота (Гц)')
    axes[2].set_ylabel('Амплитуда')
    axes[2].set_xlim(0, fs_interp / 2)
    axes[2].grid(True)

    plt.tight_layout()
    plt.savefig(f'spectrum_f{f}Hz.png', dpi=150)
    plt.show()
    
def main():
    freq_sample = 100
    freq_sample_inter = 200
    
    frequencies = [5, 10, 20, 40, 49]
    
    all_results = []
    try:
        for f in frequencies:
            filename = f'signal_f{f}Hz.csv'
            df = pd.read_csv(filename)
            
            orig = df["original"]
            interp_float = df['interp_float']
            interp_fixed = df['interp_fixed']
            
            interp_fixed = interp_fixed / 32767.0
            
            result = analyze_frequency(f, freq_sample, freq_sample_inter, orig, interp_float, interp_fixed)
            all_results.append(result)
            
            plot_spectrum(f, freq_sample, freq_sample_inter, orig, interp_float, interp_fixed)
    except FileNotFoundError:
            print(f'Файл {filename} не найден, пропускаем f={f}')
            
if __name__ == '__main__':
    main()