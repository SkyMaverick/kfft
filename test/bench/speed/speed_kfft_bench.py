#!/usr/bin/env python3

import os, subprocess, csv
import matplotlib.pyplot as plt

from sys import argv

def app_run (app, len_fft):
    print ("Testing >>", os.path.basename(app), "size -", len_fft)
    return subprocess.check_output([app, str(len_fft)], \
            stdin=None, stderr=None, shell=False, universal_newlines=True)

def gen_svg(kfft, len, step, out_dir):
    x = range (step, len, step)
    plt.figure(figsize=(40, 20), dpi= 72)

    plt.subplot(1, 1, 1)
    plt.xlabel("seq")
    plt.plot(x, kfft)
    plt.ylabel("kfft")
    plt.grid(True)

    plt.savefig(os.path.join(out_dir,"speed_test.svg"))

def gen_csv(kfft, len, step, out_dir):
    data = [range (step, len, step), kfft]
    with open(os.path.join(out_dir, "speed_test.csv"), "w+") as csv_file:
        writer = csv.writer(csv_file)
        writer.writerows(data)

def run_analize(**job):
    kfft_app    = job['kfft_app']

    len         = int(job['seq_lenght'])
    step        = int(job['seq_step'])
    
    out_dir     = job['result_dir']

    kfft_vals = []
    for i in range (step, len, step):
        kfft_vals.append(float(app_run(kfft_app, i).strip()))

    gen_svg (kfft_vals, len, step, out_dir)
    gen_csv (kfft_vals, len, step, out_dir)
    return 0

# Main app
try:
    job = {
        'kfft_app'  : os.path.abspath(argv[1]),
        'seq_lenght': argv[2],  
        'seq_step'  : argv[3],
        'result_dir': os.path.abspath(argv[4])
    }
except IndexError as e:
    print("Arguments error")
    exit (1)
else:
    exit(run_analize(**job))
