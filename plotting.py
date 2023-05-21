#! /usr/bin/env python

import os, sys, subprocess

# For plotting
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

def CreateData(n):
    """
    Creates data of size 2^n.
    """
    N = 1 << n
    
    periods = [2, 3, 20]
    amplitudes = [2, 5, 10]

    random_shift = [2* np.pi * np.random.random() for _ in periods]

    std = 0.1
    data = std * np.random.randn(N)
    data = data.astype('float32')

    for i in range(len(data)):
        for p, A, r in zip(periods, amplitudes, random_shift):
            data[i] += A * np.sin( r + (i % p) * 2 * np.pi / p)
    
    return data

def CreateStepData(n):
    """
    Creates data of size 2^n.
    """
    N = 1 << n

    period = 100
    
    data = np.zeros(N)
    for i in range(len(data)):
        if i % period >= period//2:
            data[i] = 1

    data = data.astype('float32')
    return data

def MakeInput(data):
    p_input = f"{len(data)} "
    p_input += " ".join([str(num) for num in data])

    return p_input.encode('utf-8')       

def ParseOutput(p_output):
    str_output = p_output.decode('utf-8')
    
    num_elements = int(str_output.split('\n')[0])
    data = np.array(str_output.split('\n')[1].split()).astype('float32')

    return data

def Main():
    FILENAME = "./compressor.exe"
    NUM_FREQUENCIES = 100

    NUM_FREQUENCIES = str(NUM_FREQUENCIES).encode('utf-8')

    data = CreateStepData(10)
    p_input = MakeInput(data)

    process = subprocess.run([FILENAME, NUM_FREQUENCIES], input=p_input, capture_output=True)
    try:
        process.check_returncode()
    except subprocess.CalledProcessError as e:
        print(f"ERROR in file {FILENAME}")
        print(process.stderr.decode('utf-8'))
        raise e

    p_output = process.stdout
    reconstructed_data = ParseOutput(p_output)

    fig, ax = plt.subplots(2, 1)
    sns.lineplot(data, ax=ax[0])
    sns.lineplot(reconstructed_data, ax=ax[1])

    plt.show()

if __name__ == "__main__":
    Main()