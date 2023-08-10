"Addison Powell, July 17th 2023"
from matplotlib import pyplot as plt
import numpy as np
import pickle

def comparison_plotter(exp_names, exp, pcts=5, color_key=None):
    """a function for plotting comparison plots between different RadiationForce objects
    
        Parameters:
            exp_names (list or ndarray): a list of .pkl filenames that store RadiationForce objects
            exp (str): the name of the experiment type to be compared
            pcts (int): 5, the number of unique pcts that were tested at 
                ex: the numbers (10, 15, 20, 25, 30) in the filenames suggest 5 different pcts
            color_key (list or ndarray): None, a list of string containing pyplot colors so you can
                keep the colors well ordered between different comparison_plotter function calls
                the colors correspond to exp_names by index. If None color scheme defaults to pyplots'   
            """

    plt.figure(dpi=150) 
    if color_key is None:
        color_key = [f'C{i}' for i in range(len(exp_names))] 

    #we assume the first experiments in exp_names is the the base object we compare to the others
    with open(exp_names[0], 'rb') as inF:
        base_obj = pickle.load(inF)
    S = base_obj.S
    # plot the mean and std of each expiriment aka each group of self.t files with the same pct and
    # ~the same wattage
    ind = np.searchsorted(base_obj.labels, exp)
    base_obj.plot_single(.1, .7, 1, 1, pcts, ind*base_obj.t*pcts, exp, newlabel=exp_names[0][:-4], c=color_key[0])
    i=1
    
    #we add data from different experiments only if there is a similar sample, use the change sample name func as needed
    for name in exp_names[1:]:
        with open(name, 'rb') as inF:
            obj = pickle.load(inF)

        if exp in obj.samples:
            obj.plot_single(.1, .7, 1, 1, pcts, ind*base_obj.t*pcts, exp, newlabel=exp_names[i][:-4], c=color_key[i])
        i+=1
    
    plt.tight_layout
    handles, labels = plt.gca().get_legend_handles_labels()
    ord = sorted(range(len(labels)), key = lambda k: labels[k])
    plt.legend([handles[i] for i in ord],[labels[i] for i in ord], loc='upper right')
    plt.show()
            
                
                