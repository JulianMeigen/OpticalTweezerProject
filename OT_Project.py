#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Wed Mar 22 13:40:47 2023

@author: Julian Meigen

This script calculates the significant force fluctuations from given force measurements of an Optical Tweezers (caused by breaking or shifting of actin bundles.)


1. Read and filter important inforamtion:
    a) From the Comment section - get Calibration Data (SENSITIVITY and STIFFNESS) of the Optical Tweezer.
    b) From the Datatable - get xSignal1 in [Volt] - Convert the unit to [pN] using the Sensitivity and Stiffness.
    c) From the Datatable - get Displacement [meter] - Convert the unit to [micro meter]
    d) Only use the FIRST HALF of the Data because the second half is not relevant. (The optical tweezer only returns to the initial position.)
    
2. Find PEAKS and ANTIPEAKS.
     a) Note that not always ANTIPEAK is found to PEAK -> find ANTIPEAK manually
    
3. Calculate the Force-Difference between PEAK and ANTIPEAK:

4. Find Peaks that are not usefull (New Function):
    - integrated in "get_force_and_peak_anti", only returns force differences that are meaningfull and their corresponding peak-antipeak-pairs.
    
5. Plot:
    a) xSignal1 and distance with corresponding PEAKS and ANTIPEAKS.
    b) Histogramm with all force differences.
    
6. Iterate over many different measurements in a folder and plot all the force differences found into a histogram.
    
"""

import re
import os
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import scipy.signal
import glob

# Step 1: Read and filter important inforamtion 
    
# Function that will get File as input and returns the float-values of the sensitivity and stiffness as an tuple (sensitivity, stiffness).
# Note for me: (?: ... ) indicates non-capturing groups. ([\d\.E-]*) will save the value as a string in a list with one entry (->[0]).
def get_calibration(dirName):
    with open(dirName, "r") as file:
        data = file.read()
        sensitivity = float(re.findall(r"(?:\w*_sensitivity: )([\d\.E-]*)", data)[0])
        stiffness = float(re.findall(r"(?:\w*_stiffness: )([\d\.E-]*)", data)[0])
    
    return sensitivity, stiffness


# Function that will get File as input and returns a pandas DataFrame with columnames as columns.#
def get_dataframe_from_txt(dirName):
    
    np_data = np.loadtxt(dirName)
    
    with open(dirName, "r")  as file:
        data_str = file.read()
        column_names = re.findall(r"(?:columns: )(.*)", data_str)[0]
        
        data = pd.DataFrame(np_data, columns=column_names.split())
    return data


# Function that will get two columns from the dataFrame and return only the first half of the measurements as a numpy array.
# Note for me: .to_numpy will return a 2d-array with only one column, array[:,0] will convert it to an 1d-array.
def get_relevant_data_from_dataframe(DataFrame, x_columnname, y_columnname):
    x_column = DataFrame[[x_columnname]].to_numpy()[:int(len(DataFrame)/2), 0]
    y_column = DataFrame[[y_columnname]].to_numpy()[:int(len(DataFrame)/2), 0]
    return x_column, y_column


# Function that will convert the signal (Volt) into piconewton using the sensitivitiy (m/V) and stiffness (N/m) of the optical tweezer.
def convert_signal_to_piconewton(signal, sensitivity, stiffness):
    force = signal * sensitivity * stiffness * 1E12
    return force


def convert_meter_to_micrometer(meter):
    micrometer = meter * 1E6
    return micrometer


# Function that will aplly a IRR-filters to the given signals and remove therefor the noisy-signals.
def filter_signal(signal, Wn=1/250):
    b, a = scipy.signal.iirfilter(4, Wn=Wn, btype="low", ftype="butter")
    yfilt = scipy.signal.filtfilt(b,a, signal)
    return yfilt


# Step 2: Find PEAKS and ANTIPEAKS
 
# Function that will find all peaks and antipeaks in the given data and returns a list of peak-antipeak-tuples.
def peak_antipeak_pairs(signal):
    peaks_found = scipy.signal.find_peaks(signal)[0]
    peaks = np.append(peaks_found, len(signal))
    antipeaks = scipy.signal.find_peaks(-signal)[0]
    
    peak_antipeaks = []
    for i in range(len(peaks)-1):
        peak = peaks[i]
        next_peak = peaks[i+1]
        
        possible_antipeaks = antipeaks[(antipeaks > peak) & (antipeaks < next_peak)]
        if len(possible_antipeaks) > 0: # antipeak with smallest value
            best_antipeak = possible_antipeaks[np.argmin(signal[possible_antipeaks])]
        else: # minimum between two peaks
            best_antipeak = peak + np.argmin(signal[peak:next_peak]) 
        peak_antipeaks.append((peak, best_antipeak))
        
    
    return peak_antipeaks


# Step 3 and Step 4: Calculate the Force-Difference between PEAK and ANTIPEAK and filter Peaks that are not usefull

# Function that will get the forces but also filter them and only returns meaningfull forces and their corresponding peak_antipeak_pairs
def get_force_and_peak_anti(y):
    peak_anti = peak_antipeak_pairs(y)
    force_diff = np.abs(np.diff(y[np.array((peak_anti))]))[:,0]
    
    'Filter Forces: - only get forces between 1% and 95% of the signal range'
    signal_range = max(y)-min(y)
    top_border = 0.95 * signal_range 
    bottom_border = 0.01 * signal_range 
    
    filter_force = force_diff[(force_diff < top_border) & (force_diff > bottom_border)]
    
    'Get corresponding peak_antipeaks'
    idx_lst =  np.where((force_diff < top_border) & (force_diff > bottom_border))[0]
    filter_peak_anti = [peak_anti[i] for i in idx_lst]
    
    return filter_force, filter_peak_anti


# Step 5: Plots

def plot_distance_signal(distance, signal, peak_tuples):
    fig, ax = plt.subplots()
    ax.plot(distance, signal)
    for i in peak_tuples:
        ax.plot(distance[i[0]:i[1]], signal[i[0]:i[1]], c="r")
    ax.set_xlabel("Distanz [$\u03BC$m]")
    ax.set_ylabel("Kraft [pN]")
    return plt

def plot_x_y(x, y, peak_tuples, x_name, y_name):
    fig, ax = plt.subplots()
    ax.plot(x, y)
    for i in peak_tuples:
        ax.plot(x[i[0]:i[1]], y[i[0]:i[1]], c="r")
    ax.set_xlabel(x_name)
    ax.set_ylabel(y_name)
    return plt

def plot_histogram(data):
    fig, ax = plt.subplots()
    n, bins, patches = ax.hist(data, color="b", bins = int(len(data)/4), width=0.5)
    ax.set_xlabel("Kraft [pN]")
    ax.set_ylabel("Häufigkeit")
    ax.set_title("Histogramm der Kraftdifferenzen")
    ax.grid()
    return plt



# Step 6: Check user input and plot one file or iterate over many different measurements in a folder

def CheckFolder(dirName):
    """
    Checks if dirName is:
        - single .txt file (False)
        - Folder containing .txt files (True)
    
    """
    if dirName.endswith(".txt"):
        return False
    else:
        files = glob.glob(f"{dirName}/*.txt")
        if len(files) > 0:
            return True
        else:
            raise Exception(f"No existing .txt files in {dirName}")
            
            
def CheckCalibration(dirName):
    
    if args.calibration is not None:
        sens = int(input("Sensitivity in m/V: "))
        stiff = int(input("Stiffness in N/m: "))
    else:
        sens, stiff = get_calibration(dirName)
    return sens, stiff


def CheckDefault_and_convert(file, x_name, y_name, sens, stiff):
    df = get_dataframe_from_txt(file)
    x,y = get_relevant_data_from_dataframe(df, x_name, y_name)
    
    if (x_name == "distance") &  (y_name == "xSignal1"): # default
        x = convert_meter_to_micrometer(x)
        y = convert_signal_to_piconewton(y, sens, stiff)
        y_filt = filter_signal(y)
    else:
        y_filt = filter_signal(y)
    return x, y_filt





def main(args=None):

    if CheckFolder(args.filename): # Folder with .txt files
        files = glob.glob(f"{args.filename}/*.txt")
        
        if args.calibration is not None:
            sens = float(input("Sensitivity in m/V: "))
            stiff = float(input("Stiffness in N/m: "))
            
        full_force_diff = np.array([])
        for file in files:
            if args.calibration is None:
                sens, stiff = get_calibration(file)
            x, y = CheckDefault_and_convert(file, args.x_axis, args.y_axis, sens, stiff)
            force_diff, peak_anti = get_force_and_peak_anti(y)
            full_force_diff = np.append(full_force_diff, force_diff)
            
            plot_force = plot_distance_signal(x, y, peak_anti)
            
            base = os.path.basename(file)
            filename = os.path.splitext(base)[0]
            plot_force.savefig(f"{filename}.png")
            
        plot_hist = plot_histogram(full_force_diff)
        plot_hist.savefig("Histogramm.png")
        
    
    else: # single .txt file
        sens, stiff = CheckCalibration(args.filename)
        
        x, y = CheckDefault_and_convert(args.filename, args.x_axis, args.y_axis, sens, stiff)
        force_diff, peak_anti = get_force_and_peak_anti(y)
        
        if (args.x_axis == "distance") &  (args.y_axis == "xSignal1"): # default
            plot_force = plot_distance_signal(x, y, peak_anti)
            
            base = os.path.basename(args.filename)
            filename = os.path.splitext(base)[0]
            plot_force.savefig(f"{filename}.png")
        
        else:
            plot_force = plot_x_y(x, y, peak_anti, args.x_axis, args.y_axis)
            
            base = os.path.basename(args.filename)
            filename = os.path.splitext(base)[0]
            plot_force.savefig(f"{filename}.png")

        
        
        


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
                    prog = 'Calculate force fluctuations',
                    description = 'Calculates the significant force fluctuations from given force measurements of an Optical Tweezer.',)
    
    parser.add_argument('-f', '--filename', required=True,
                    help='give Folder Path with multiple Optical Tweezer measurements as .txt files OR single measurement as .txt file')
    parser.add_argument('-x', '--x_axis', required=False, default="distance", help='give different column as x-axis')
    parser.add_argument('-y', '--y_axis', required=False, default="xSignal1", help='give different column as y-axis')
    parser.add_argument('-c', '--calibration', required=False, help='give own calibration measurements')

    args = parser.parse_args()
    
    main(args)