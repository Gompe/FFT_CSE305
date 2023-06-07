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
    
    n_periods = np.random.randint(1, 10)

    possible_periods = list(range(1, 21))
    possible_amplitudes = list(range(1, 11))

    periods = np.random.choice(possible_periods, n_periods, replace=False)
    amplitudes = np.random.choice(possible_amplitudes, n_periods, replace=True)

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

def ShowApproximationQuality(data, r_data):
    mae = np.linalg.norm(data - r_data, ord=1) / len(data)
    mrse = np.sqrt(np.linalg.norm(data - r_data, ord=2) / len(data))

    print("Mean Absolute Error:", round(mae, 3))
    print("Mean Root Squared Error:", round(mrse, 3))

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

    MODE = int(input("""Do you wish to create 
        [1] Random Periodic Data (type 1)
        [2] Pulse-like data (type 2)
        Your choice: """))
    
    if MODE not in [1, 2]:
        print("Invalid choice")
        sys.exit()

    NUM_FREQUENCIES = int(input("Choose the number of frequencies that you want to use for the compression: "))
    NUM_FREQUENCIES = str(NUM_FREQUENCIES).encode('utf-8')

    if MODE == 1:
        data = CreateData(10)
    else:
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

    ShowApproximationQuality(data, reconstructed_data)

    fig, ax = plt.subplots(2, 1)
    sns.lineplot(data, ax=ax[0])
    sns.lineplot(reconstructed_data, ax=ax[1])

    print("Showing reconstruction plot:")
    plt.show()

if __name__ == "__main__":
    Main()