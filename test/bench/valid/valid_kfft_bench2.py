#!/usr/bin/env python3
import os, subprocess, random
from sys import argv

def app_run (app, len_fft):
    print ("Validate >>", os.path.basename(app), "size -", len_fft)
    try:
        return subprocess.call([app, str(len_fft)], \
                stdin=None, stderr=None, shell=False, universal_newlines=True)
    except subprocess.CalledProcessError as e:
        return 3

def run_analize(**job):
    comp_app    = job['comp_app']
    limit       = int(job['seq_limit'])
    count       = int(job['work_count'])

    examples = list(range(1, int(limit) + 1, count))
    print (examples)
    if os.path.isfile(comp_app):
        for size in examples:
            if (app_run(comp_app, size) != 0):
                return 1
    else:
        return 2
    return 0

# Main app
try:
    job = {
        'comp_app'  : os.path.abspath(argv[1]),
        'seq_limit': argv[2],  
        'work_count': argv[3],
    }
except IndexError as e:
    print("Arguments error")
    exit (1)
else:
    exit(run_analize(**job))
