import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from scipy.optimize import curve_fit
import math
import sys

molecules = {"1": "H2O", "2": "H2S"}

for i in molecules: 
    print("SELECT ",i,"FOR ", molecules[i])

x = input()

try:
    molecule = molecules[x]
except KeyError:
    sys.exit("Invalid input")

def data_extractor(molecule):
    """
    Function takes in molecule from user input and returns an array in the order of bond length, bond angle and energy.
    I used os.sep here to avoid issues with different separators in directories between different operating systems
    Data first stored as a list of lists and converted to np array.
    bond length and angle parsed from the the filename(may not be idea as the number of digits might change).
    """
    data = []
    for file in os.listdir(molecule):
        f = open(molecule+os.sep+file,"r")
        for line in f:
            if "SCF Done" in line:
                l = line.split()
                data.append([float(file[5:9]), float(file[14:17]), float(l[4])])
        f.close()
    
    data = np.array(data)
    return data

def surface_plotting(data):
    """
    Takes in data and returns an energy surface plot.
    Used a tri-surface plot here for 3 variables
    """
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.plot_trisurf(data[:,0],data[:,1],data[:,2], cmap = cm.coolwarm)
    ax.set_xlabel('r/Angstroms')
    ax.set_ylabel('Theta/degrees')
    ax.set_zlabel('Energy/Hartrees')
    ax.title.set_text("Potential energy sufarce for "+ molecule)
    plt.savefig("output_"+molecule)
    return

def find_min_values(data):
    """Returns the minimum energy of the molecule and its geometry"""
    index = np.argmin(data[:,2])
    return data[index,:]

def update_units(data):
    """converts to standard units"""
    data = data * np.array([10**-10, math.pi/180, 4.3598*10**-18])
    return data

def vibrational_frequency(data,kr,kt):
    "Takes in (picked) data points, kr, kt and returns energy "
    r = data[:,0]
    t = data[:,1]
    r0,t0,E0 = find_min_values(data)
    energy = E0 + 0.5*kr*((r-r0)**2) + 0.5*kt*((t-t0)**2)
    return energy

def fit_vibrational_frequency(data):
    """
    picks a fixed number of points near the energy mimimum on the surface.
    use curve_fit to fit the picked points to the equation for vibrational energy.
    Returns vibrational frequency in wavenumbers.
    """
    number_of_points = 0
    buffer = 10**-20

    while number_of_points < 50:
        picked_data = []
        buffer = buffer *2
        for row in data:
            if row[2] < find_min_values(data)[2] + buffer:
                picked_data.append(row)
        number_of_points = len(picked_data)
        print(number_of_points)
    
    picked_data = np.array(picked_data)
    picked_data = update_units(picked_data)
  
    r0 = find_min_values(picked_data)[0]

    popt,pcov = curve_fit(vibrational_frequency, picked_data, picked_data[:,2],)
    kr,kt = popt

    mu1 = 2 * 1.6605 * 10**-27 
    mu2 = 0.5 * 1.6605 * 10**-27

    v1 = (1/(2*math.pi))* np.sqrt(kr/mu1) * 3.3356*10**-11
    v2 = (1/(2*math.pi))* np.sqrt(kt/((r0**2)*mu2)) * 3.3356*10**-11

    print("Vibrational frequency 1 is: ", v1,"cm-1")
    print("Vibrational frequency 2 is: ", v2,"cm-1")
    return

data = data_extractor(molecule)
surface_plotting(data)
fit_vibrational_frequency(data)
