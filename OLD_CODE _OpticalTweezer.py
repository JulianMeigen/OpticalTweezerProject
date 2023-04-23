# -*- coding: utf-8 -*-
"""
Created on Tue Oct  4 23:54:38 2022

@author: julia
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
import glob
from scipy.stats import norm
import matplotlib.mlab as mlab

# Es muss OpticalTweezerKraft() eingegeben werden und dann der Ordnerpfad mit den Messungen sowie die eingerichtete Sensitivität und Steifigkeit der benutzten Falle.
# Mit graph=True können alle Graphen erzeugt und überprüft werden
# Also z.B OpticalTweezerKraft("20220919/processed_curves-20220919", sensitivity = 43, stiffness=0.2094, graph=True)
# Mit steps kann die Anzahl der eingelesenen Daten geändert werden. Für bessere ergebnisse wurde steps auf 10 gestellt. (liest nur jeden 10 Wert ein)
# OPTINAL: Mit var=False kann die Art geändert weren indem die "lows" bestimmt werden (einfach jeder 50 schritt nach einem Peak), kann ggfls für bessere ergebnisse sorgen.

def OpticalTweezerKraft(dirName, var=True, graph=False, steps=10,  sensitivity=42.556, stiffness=0.2094):
    files = glob.glob(f"{dirName}/*.txt")
    
    all_diff= []
    for file in files:
        data = np.loadtxt(file)
        data_signalx_volt = data[:int(len(data)/2):steps,0] # Betrechtet nur die erste Hälfte der Daten.
        data_displacement_m = data[:int(len(data)/2):steps, 13]
        data_signalx = (data_signalx_volt / (sensitivity*0.001))*stiffness
        data_displacement= data_displacement_m*1_000_000
        
        peaks = find_peaks_cwt(data_signalx[:], prominence=3.7) # findet peaks und lows
        lows = find_peaks(-data_signalx[:], prominence=3.7)
        data_peaks = data_signalx[peaks[0]]
        
        
        data_foot =  np.zeros(shape=len(peaks[0]))
        data_foot_idx =  np.zeros(shape=len(peaks[0]), dtype=(int))
        
        if var: # NORMALFALL: Berechnet automatisch mit scipy peaks UND low der Funktion und nimmt (wenn es passt) diesen low-Wert.
        
            if len(peaks[0]) == 0:
                continue
            
            j = 0
            for i in range(len(peaks[0])-1):
                peak_idx = peaks[0][i]
                next_peak_idx = peaks[0][i+1]
                
                if j < len(lows[0]):
                    low_idx = lows[0][j]
                else:
                    low_idx = 0
            
                if peak_idx < low_idx < next_peak_idx:
                    
                    data_foot[i] = data_signalx[low_idx]
                    data_foot_idx[i] = low_idx
                    j += 1
                else:
                    data_range = data_signalx[peak_idx:peak_idx+(int(500 / steps))] # nimmt sich den niedrigesten Werten aus den nächsten 30 werten.
                    data_foot[i] = np.min(data_range)
                    data_foot_idx[i] = np.argmin(data_range) + peak_idx
                
            peak_last_idx = peaks[0][-1]
            data_range_end = data_signalx[peak_last_idx:peak_last_idx+(int(500 / steps))]
            data_foot[-1] = np.min(data_range_end) 
            data_foot_idx[-1] = np.argmin(data_range_end) + peak_last_idx
          
        else:   # nur Zusatz!! durch übergeben von var=False kann man die lows ohne scipy finden. Kann in ausnahmefällen das Ergebnis verbessern.

            
            for i in range(len(peaks[0])):
                peak_idx = peaks[0][i]
                data_range = data_signalx[peak_idx:peak_idx+(int(500 / steps))] # nimmt sich den niedrigesten Werten aus den nächsten 30 werten.
                data_foot[i] = np.min(data_range)
                data_foot_idx[i] = np.argmin(data_range) + peak_idx
                
                
        # Berechnet die Differenz der gefundenen Peaks und Lows. (unabhängig auf von var=True/False)
        data_2d = np.vstack((data_peaks, data_foot))
        data_diff = np.absolute(np.diff(data_2d, axis=0)[0])
        
        # Fügt die berechneten Differenzen der all_diff Liste an, die außerhalb der for Schleife definiert ist und so die Differenzen aller Datein sammelt.
        for value in data_diff:
            all_diff.append(value)
            
        # Fragt ab ob (alle) einzelne Graphen mit den gefundenen peaks und lows ausgegeben werden sollen.
        # Bei graph=True werden Graphen ausgegeben. Dient somit als mögliche Überprüfung der Qualität der erhaltenen Werte.
        if graph:
            print(f"Dies ist der Graph der Datei: {file}")
            fig2, ax = plt.subplots()
            ax.plot(data_displacement, data_signalx)
            ax.scatter(data_displacement[peaks[0]], data_signalx[peaks[0]], color="r")
            ax.scatter(data_displacement[data_foot_idx], data_signalx[data_foot_idx], color="b") # nur als platzhalter, plottet halt nur die gefundenen lows von find_peaks 
            ax.set_xlabel("Distanz [$\u03BC$m]")
            ax.set_ylabel("Kraft [pN]")
    
    
    all_diff_round = np.round(np.array(all_diff)*2)/2
    
    x, y = np.unique(all_diff_round, return_counts=True)
    
    all_diff_array = np.array(all_diff)
    diff_corr = sorted(all_diff_array[all_diff_array < 40])
    
    
    fig, ax = plt.subplots()
    n, bins, patches = ax.hist(diff_corr, bins=22, color="b", width=0.5)
    ax.set_xlabel("Kraft [pN]")
    ax.set_ylabel("Häufigkeit")
    ax.set_title("Histogramm der Kraftdifferenzen")
    ax.grid()
    
    
    mean_all  = np.mean(diff_corr)
    std_all = np.std(diff_corr)
    std_mean = std_all/np.sqrt(len(diff_corr))
    
    start = int(0)
    stop = int(n[0])
    n_array = np.array(n)
    n_int = n_array.astype(int)
    bars_std = []
    bars_mean = []
    #mis_all = []
    for i in range(len(n)-1):
        bar_values = diff_corr[start:stop]
       
        bar_mean = np.mean(bar_values)
        bars_mean.append(bar_mean)
        
        bar_std = np.std(bar_values)
        bars_std.append(bar_std)
        
        # single_mis_raw = []
        # for j in bar_values:
        #     val_minus_mean = (j - mean_all)**2
        #     single_mis_raw.append(val_minus_mean)
            
        # single_mis_sum = np.nansum(np.array(single_mis_raw))
        # single_mis = np.sqrt(single_mis_sum/len(diff_corr))
        # mis_all.append(single_mis)
        
        
        start = start + n_int[i]
        stop = stop + n_int[i+1]
    bars_std.append(np.std(diff_corr[start:]))
    bars_mean.append(np.mean(diff_corr[start:]))
    
    fehler = []
    for i in bars_mean:
        bar_fehler = np.sqrt(((i-mean_all)**2)/len(diff_corr))
        fehler.append(bar_fehler)
    
    var = np.array(fehler)**2
    fehlerfort = np.sqrt(np.nansum(var))
    
    
    #var = np.array(bars_std)**2
    # nansum = np.nansum(sqrt_std)
    #fehler = np.sqrt(np.nansum(var))
   
    (mu, sigma) = norm.fit(diff_corr)
    best_fit = norm.pdf(bins, mu, sigma)
    ax.plot(bins, best_fit, "r--") 
    
    
    
    return(mean_all, std_all, std_mean, fehler, fehlerfort, len(diff_corr), bars_mean)

print(OpticalTweezerKraft("20220919/processed_curves-20220919",steps=10, sensitivity = 43, stiffness=0.2094, graph=True))

