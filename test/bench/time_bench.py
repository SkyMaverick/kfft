#!/usr/bin/env python3

import os, subprocess
import matplotlib.pyplot as plt

from sys import argv

def app_run (app, len_fft):
    print (app, len_fft)
    return subprocess.check_output([app, str(len_fft)], \
            stdin=None, stderr=None, shell=False, universal_newlines=True)

def gen_svg(kfft, fftw, len, step):
    x = range (step, len, step)
    plt.figure(figsize=(40, 20), dpi= 72)

    plt.subplot(2, 1, 1)
    plt.xlabel("seq")  # ось абсцисс
    plt.plot(x, kfft)               # построение графика
#    plt.title("Зависимости: y1 = x, y2 = x^2") # заголовок
    plt.ylabel("kfft") # ось ординат
    plt.grid(True)                # включение отображение сетки
    plt.subplot(2, 1, 2)
    plt.plot(x, fftw)               # построение графика
    plt.xlabel("seq")  # ось абсцисс
    plt.ylabel("fftw") # ось ординат
    plt.grid(True)                # включение отображение сетки
    plt.subplot(2, 1, 1)

    plt.savefig('bench.svg')

def run_analize(**job):
    kfft_app    = job['kfft_app']
    fftw_app    = job['fftw_app']
    len         = int(job['seq_lenght'])
    step        = int(job['seq_step'])
    
    kfft_vals = []
    for i in range (step, len, step):
        kfft_vals.append(app_run(kfft_app, i))
    fftw_vals = []
    for i in range (step, len, step):
        fftw_vals.append(app_run(fftw_app, i))
        
    kfft_vals.sort()
    fftw_vals.sort()

    gen_svg (kfft_vals, fftw_vals, len, step)
    return 0

# Main app
try:
    job = {
        'kfft_app'  : os.path.abspath(argv[1]),
        'fftw_app'  : os.path.abspath(argv[2]),
        'seq_lenght': argv[3],  
        'seq_step'  : argv[4]
    }
except IndexError as e:
    print("Arguments error")
    exit (1)
else:
    exit(run_analize(**job))
