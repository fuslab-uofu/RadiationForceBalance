"""Addison Powell, May 15 2023

Python 3.11.3
all imported packages are up to date as of 08/08/23

Calculates the attenuation from a radiation force balance experiment and has various functions for plotting and
storing the data.  must follow proper file naming conventions: sample_##pct_T#_##p#W, scrapes termite logs and the 
new directory for writing any new files must be a sub directory of the current working directory.
the root file name (the folder with all the termite logs) must also be a sub directory
"""

import re, os
from scipy.io import savemat
import numpy as np
import copy
from matplotlib import pyplot as plt
import pickle

def save_object(obj, filename):
    """a function for saving python class objects to a file for later access"""
    with open(filename, 'wb') as outp:  # Overwrites any existing file.
        pickle.dump(obj, outp, pickle.HIGHEST_PROTOCOL)

class RadiationForce():
    """A class for calculating, storing, plotting and exporting radiation force balance data

        Attributes:
            self.rootfn (str): rootfn (see parameters)
            self.t (int): trials (see parameters)
            self.samples (list or ndarray): samples
            self.scale (float): transducer
            self.controls (list or ndarray): controls
            self.labels (list or ndarray): keys for the dictionary made from the samples
            and controls
            self.widths (dict): a dictionary mapping the samples to their respective
            widths
            self.raw (dict): the parsed table of data formatted the same as self.data but without the power 
            of sample and alpha columns
            self.data (dict): a dictionary mapping each sample expirimented on to its data
                keys (str or int): filename
                values (list): A list storing the following
                        expiriment    watts   pct     trial#  grams (uniform time)  power of sample  alpha
                    [[str,          float,  float,  float,  list,                 float,           float]
                        # also find water ref trial 1 create to_string method
                there are only 3 columns and the number of rows is determined by the
                number of trials per sample expirimented on.
            self.pw (list): a list of floats each float corresponding to the power of
            ultrasound trough water given a certain pct 
            self.S (ndarray): a matrix of the data from self.data but only the parts with expriment types 
            from self.samples and excludes
            any data that is not type float64
            self.C (ndarray): a matrix of the data from self.data but only the parts with expriment types 
            from self.controls and excludes
            any data that is not type float64

        Functions:
            Most of these have more detailed descriptions in the actual function
            __init__ - initialize a RadiationForce object for an experiment
            parser - ReGeX text parser for reading in .txt files  (automatically used, not really intended for users)
            powerSample - calculates the power of the samples  (automatically used, not really intended for users)
            write - writes a table of the data to a .txt file
            saveMat - saves the results to a .mat file
            plotRaw - plots the results for each file
            plot_single - plot the results of a single file
            plot - plots the primary results of the experiment
            __str__ - a tostring method for printing a RadiationForce object
            find_alpha - a function for calculating attenuation  (automatically used, not really intended for users)
            plotRawSingle - a function for plotting the raw data of a single file
            comparePlot - a function for plotting the results on a single plot
            change_samplename - a function for renaming a sample, this helps with the comparison_plotter.py which
                needs matching sample names
            stats - a function for calculating the mean and standard dev of attenuations
    """
    def __init__(self, rootfn, trials, samples, controls, widths, transducer=0.919, length_correction=1.076238605618438):
        """Parameters:
                rootfn (str): a string containing the folder that has radiation force
                balance data.
                trials (int): an integer number of how many trials were done per 
                different pct values.
                samples (list or 1darray): a list of strings or integers to act as keys,
                each key should act as a label of the different samples that were expiremented on.
                controls (list or 1darray): a list of strings or integers indicating the label
                of control expiriment files.
                widths (list or 1darray): a list of floats that represent the widths of the 
                samples in centimeters, there should be as many widths as there are samples.
                transducer (float): defualt = 0.919, a float representing the scaling coefficient of 
                the transducer that was used. The default is for Allison's breast transducer.
                length_correction (float): default = 1.076238605618438, aa float representing the average
                actual length of a wave propagating through a sample at an angle. The default is for 
                Allison's breast transducer.
        """
        #store the rootfn, create the keys and widths attributes
        self.rootfn = rootfn
        self.t = trials
        self.samples = samples
        self.scale = transducer
        self.lc = length_correction
        self.controls = controls
        self.labels = np.append(samples, controls)
        self.widths = {samples[i]: widths[i] for i in range(len(samples))}

        #skim and calculate the data as well as compute alphas, avg watts, and pw       
        self.raw = self.parser() 
        self.data = copy.deepcopy(self.raw)
        self.powerSample()

        #sort the data by alphabet
        alphabet = list(self.data.keys())
        alphabet.sort()
        self.data = {i: self.data[i] for i in alphabet}
        self.labels = np.sort(self.labels)
        self.pw = self.find_alpha()

        #make a matrix of the sample data for easier access and broadcasting
        samplekeys = [k for k in self.data.keys() if self.data[k][0] in samples]
        M1 = np.array([self.data[i][1:4] for i in samplekeys])
        M2 = np.array([self.data[i][5:] for i in samplekeys])
        self.S = np.concatenate((M1, M2), axis=1)

        #make a matrix of the control data for easier access and broadcasting
        controlkeys = [k for k in self.data.keys() if self.data[k][0] in controls]
        M1= np.array([self.data[i][1:4] for i in controlkeys])
        M2=np.array([self.data[i][5:-1] for i in controlkeys])
        self.C = np.concatenate((M1, M2), axis=1)
        self.samples = list(np.sort(samples))
        self.controls = list(np.sort(controls))

    def parser(self):
        """a ReGeX parser for skimming data from the termite logs in the self.rootfn directory so long as they
        followed proper naming conventions.
        returns:
            dict, the raw data that has been parsed
        """
        #walk trough each file in the given directory
        raw = dict()
        for root, dir, files in os.walk(self.rootfn):
            for file in files:

                #only read files that have data (aka named properly)
                if re.search(r"pct_T\d+_", file) and file.startswith('.') != True:
                    filename = os.path.join(root, file)
                else:
                    continue

                #This code gets the pct watts and trial# from the filename
                pct = float(re.search(r"(\d+)(?=pct)", file).group())
                watts = re.sub(r"(\d+)p?(\d?)(?=W)", r"\1.\2", file)
                watts = float(re.search(r"(\d+).?(\d?)(?=W)", watts).group())
                trial =  int(re.findall('\d+', re.search(r"_T\d+_", file).group())[0])

                #this code gets the name of the sample from the file in question
                ind = [1 if self.labels[i] in file else 0 for i in range(len(self.labels))].index(1)
                extype = self.labels[ind]

                #read in the file
                with open(filename, 'r') as inF:
                    data = [line.strip() for line in inF.readlines() if line.strip()]

                #get rid of the unessecary information (the time because we assume uniform spacing)
                dates = r'^.*?[+-]' # r'^\d{2}:\d{2}:\d{2}\.\d{2}:'
                signs = np.array([re.findall(r'([+\-])', l)[0] for l in data])
                signs = [1 if s == '-' else -1 for s in signs]
                data = [re.sub(dates, '', line) for line in data]
                numbers = np.array([float(re.findall(r'(\d+\.\d+|\d+)', d)[0]) for d in data])
                numbers *= signs

                raw[file] = [extype, watts, pct, trial, numbers]
        return raw
    
    def powerSample(self):
        """Converts weights measured in grams to a single number 
        that represents the power of the sample for that individual expiriment."""
        #calculate mean powers of the samples in each file
        for k in self.raw.keys():
            tol = np.std(self.raw[k][-1])
            zeroish = np.mean(self.raw[k][-1][self.raw[k][-1] < tol])
            weight = np.mean(self.raw[k][-1][self.raw[k][-1] >= tol])
            self.data[k].append((abs(weight - zeroish))*14.715)

    def write(self, filename='experiment_results_data.txt'):
        """writes to an outfile the data and results calculated from the expiriment"""
        with open(filename, 'w') as f:
            f.write(self.__str__()+'\n\nRaw Data:\n' +self.data.__str__() + self.S.__str__() + self.C.__str__())

    def saveMat(self, newdir): ##untested, who uses matlab anyways?
        """a function for saving the parsed data as a .m (matlab) file"""
        #walk through the directory again for the sake of getting each filename,
        for root, dir, files in os.walk(self.rootfn):
            d = {}
            for file in files:
                d[file] = self.raw[file]
            savemat("results", {file:self.raw[file]})
            os.rename("results", os.path.join(newdir, file))

    def plotRaw(self, ymin=0.1, ymax=0.7, pcts=5, show=True):
        """a function for plotting the raw data over a uniform time interval and then plotting the individual
        calculated data
        Parameters:
            ymin (float) = 0.1, the minimum bound for the y axis
            ymax (float) = 0.7, the maximum bound for the x axis
            pcts (int) = 5, the number of distinct pct values the expirement was done at
            show (bool) = true, whether or not to show the plot or leave on stack for show call
        """

        S = self.S
        plt.figure(dpi=200)
        #the subplot dimensions may need to be implemented as def args
        for i, k, l in zip(
            range(len(self.samples)), range(0, len(S), self.t*pcts), self.samples):

            #plot alpha vs watts for each sample
            j = k+self.t*pcts
            plt.subplot(2,3,i+1)
            plt.ylim((ymin, ymax))
            plt.scatter(S[k:j, 0], S[k:j, -1], label=l, s=5)
            plt.xlabel("Watts")
            plt.ylabel("α, nepers/cm")
            plt.title(f"{l}")

        for i, k, l in zip(
            range(len(self.samples)), range(0, len(S), self.t*pcts), self.samples):

            #plot alpha vs pct for each sample
            plt.subplot(2,3,i+4)
            plt.ylim((ymin, ymax))
            j = k+self.t*pcts
            plt.scatter(S[k:j, 1], S[k:j, -1], label=l, s=5)
            plt.xlabel("pct")
            plt.ylabel("α, nepers/cm")
            plt.title(f"{l}")
        plt.tight_layout()
        if show: plt.show()

    def plot_single(self, ymin, ymax, dim1, dim2, pcts, k, l, subplot=None, newlabel=None, c=None):
        """a function for plotting a single experiment
        Parameters:
            ymin (float) the minimum bound for the y axis
            ymax (float) the maximum bound for the x axis
            dim1 (int) the first dimension for the shape of the subplot
            dim2 (int) the second dimension for the shape of the subplot
            pcts (int) the number of distinct pct values the expirement was done at
            subplot (int) = None, the subplot number, defaults to none for a single plot
            newlabel = None, an optional label to overide the default label for the plot legend
        """
        S = self.S
        means = []
        stds = []
        domain = []
        if subplot is not None: plt.subplot(dim1, dim2, subplot+1)
        plt.ylim((ymin, ymax))
        #calculate mean/std for each block of data with matching pct and label
        for x in range(pcts):
            means.append(np.mean(S[(self.t*x)+k:self.t*(x+1)+k, -1])) 
            stds.append(np.std(S[(self.t*x)+k:self.t*(x+1)+k, -1]))
            domain.append(np.mean(S[(self.t*x)+k:self.t*(x+1)+k, 1]))
        plt.errorbar(domain, means, stds, ecolor=c, color=c ,label=l if newlabel is None else newlabel, linestyle="None", marker=".", capsize=3)
        plt.xlabel("pct")
        plt.ylabel("α, nepers/cm")
        plt.title(f"{l}")

    def plot(self, ymin=0.1, ymax=0.7, pcts=5, show=True):
        """a function for plotting some of the data, mainly the alphas and the power output
        Parameters:
            ymin (float) = 0.1, the minimum bound for the y axis
            ymax (float) = 0.7, the maximum bound for the x axis
            pcts (int) = 5, the number of distinct pct values the expirement was done at
            show (bool) = true, whether or not to show the plot or leave on stack for show call
            """
        S = self.S
        plt.figure(dpi=150)  
        
        # plot the mean and std of each expiriment aka each group of self.t files with the same pct and
        # ~the same wattage
        for i, k, l in zip(
            range(len(self.samples)), range(0, len(S), self.t*pcts), self.samples):
            self.plot_single(ymin, ymax, 3, 3, pcts, k, l, i)     
            

        #plot the power output scaled by the transducer efficieny coefficient for the samples
        for i, k, l in zip(
            range(len(self.samples)), range(0, len(S), self.t*pcts), self.samples):
            plt.subplot(3,3,i+4)
            means = []
            stds = []
            domain = []
            plt.ylim((0, 30))
            for x in range(pcts):
                means.append(np.mean(S[(self.t*x)+k:self.t*(x+1)+k, -2]/self.scale)) 
                stds.append(np.std(S[(self.t*x)+k:self.t*(x+1)+k, -2]/self.scale))
                domain.append(np.mean(S[(self.t*x)+k:self.t*(x+1)+k, 0]))
            plt.errorbar(domain, means, stds, label=l, linestyle="None", marker=".", capsize=3)
            plt.xlabel("Power in (W)")
            plt.ylabel("Power out (W)")
            plt.title(f"{l}")
        
        #plot the power output scaled by the transducer efficieny coefficient for the controls
        C=self.C
        for i, k, l in zip(
            range(len(self.controls)), range(0, len(C), self.t*pcts), self.controls):
            plt.subplot(3,3,i+7)
            means = []
            stds = []
            domain = []
            plt.ylim((0, 40))
            for x in range(pcts):
                means.append(np.mean(C[(self.t*x)+k:self.t*(x+1)+k, -1]*(1/self.scale))) 
                stds.append(np.std(C[(self.t*x)+k:self.t*(x+1)+k, -1]*(1/self.scale)))
                domain.append(np.mean(C[(self.t*x)+k:self.t*(x+1)+k, 0]))
            plt.errorbar(domain, means, stds, label=l, linestyle="None", marker=".", capsize=3)
            plt.xlabel("Power in (W)")
            plt.ylabel("Power out (W)")
            plt.title(f"{l}")

        plt.tight_layout()
        if show: plt.show()

    def __str__(self):
        """to string magic method for outputting a table of the data parsed and the data calculated"""
        headers = ["experiment", "watts", "pct", "trial#", "power of sample", "attenuation constant (nepers/cm)"]
        header_line = "filename                    |    " + "    |    ".join(headers)
        separator_line = "-" * len(header_line)

        lines = [header_line, separator_line]

        for key, value in self.data.items():
            line = "{:<33}{:<19}{:<14}{:<14}{:<12}{:<25}{:<14}".format(
                key,value[0],value[1],value[2],value[3],value[5],value[6])
            lines.append(line)

        return "\n".join(lines)
    
    def find_alpha(self):
        """a function for calculating the attenuation constants alpha as well
        as the average wattage and the power of water from one of the control expiriments
        that has been marked in the filename with keyword "water"
        returns:
            power_W (list): a list of floats of the power of water per expiriment        
        """
        datakeys = list(self.data.keys())
        power_W = {}
        f=0
        for label in self.labels:
            #calculate the average power of water per pct conditions
            file_groups = {}
            group_keys = []
            #also create file groupings of self.t len
            while f< len(self.data) and self.data[datakeys[f]][0] == label:
                group = []
                pcts = []
                expect = []
                for i in range(self.t):
                    group.append(datakeys[f+i])
                    pcts.append(self.data[datakeys[f+i]][2])
                    expect.append(self.data[datakeys[f]][2])
                group_keys.append(label+f"_pct{expect[0]}")
                #check that the pct conditions are being met
                if pcts != expect:
                    raise TypeError("wrong expirments being meaned")
                file_groups[label+f"_pct{expect[0]}"] = group
                f+=self.t

            #use the file groupings to find the water reference and calculate the power of the water reference
        for group in group_keys:
            if group.startswith('water') or group.startswith('Water'):
                for k in file_groups[group]:
                    temp = [self.data[i][5] for i in file_groups[group]]
                    power_W[self.data[k][2]] = np.mean(temp)

        #iterate through each file expirimented on
        for key in datakeys:
            if self.data[key][0] not in self.samples:
                self.data[key].append("None")
                continue
            #solve for α = transducer*-ln(Ps/Pw)/(2*d)
            p = power_W[self.data[key][2]]
            self.data[key].append(
                -1*np.log(self.data[key][5]*(p**(-1)))/(self.lc*2*self.widths[self.data[key][0]]))
        
            #returns the pw
        return power_W
    
    def plotRawSingle(self, key, show=True):
        """a function for plotting the raw data of a single file
        Parameters:
            key (str) the filename
            show (bool) = true, whether or not to show the plot or leave on stack for show call
        """
        domain = np.linspace(0,1, len(self.raw[key][-1]))
        plt.scatter(domain, self.raw[key][-1], label="data")
        plt.ylabel("grams")
        plt.xlabel("uniform time")
        plt.hlines(np.std(self.raw[key][-1]), 0, 1, colors='r', label="1 standard dev.")
        plt.legend(bbox_to_anchor=(1.0, 1.0), loc="upper right")
        plt.title(f"{key}, {self.raw[key][:4]}")
        if show: plt.show()

    def comparePlot(self, ymin=0.1, ymax=0.7, pcts=5, show=True):
        """a function for plotting all of the attenuation coeffs on the same plot for reference
        Parameters:
            ymin (float) = 0.1, the minimum bound for the y axis
            ymax (float) = 0.7, the maximum bound for the x axis
            pcts (int) = 5, the number of distinct pct values the expirement was done at
            show (bool) = true, whether or not to show the plot or leave on stack for show call
        """

        k=0
        plt.figure(dpi=120)
        plt.xlabel("Watts")
        plt.ylabel("α, nepers/cm")
        plt.ylim((ymin, ymax)) 
        for i, l in zip(range(len(self.samples)), self.samples):                     
            plt.scatter(self.S[k:k+self.t*pcts, 0], self.S[k:k+self.t*pcts, -1], label=l, s=5)
            k +=self.t*pcts
        plt.title("comparison plot")
        plt.legend(loc="lower right")
        plt.tight_layout()
        if show: plt.show()

    def change_samplename(self, old, new):
        """a function for renaming a sample name, only really updates the plot titles, 
        cannot update the raw data or the readout text file"""

        #change either the sample name or the control name
        try:
            i = self.samples.index(old)
            self.samples[i] = new
            #update the widths dictionary if needed
            self.widths[new] = self.widths[old]
            del self.widths[old]
        except ValueError:
            pass

        try:
            i = self.controls.index(old)
            self.controls[i] = new
        except ValueError:
            pass

        #update the labels attribute
        self.labels = np.append(self.controls, self.labels)
        self.labels = np.sort(self.labels)

    def stats(self, exp_names, show=False):
        """returns and prints the stats (mean, std) of the attenuations of the
        given experiments

        Parameters:
            exp_names (list or ndarray): a list containing strings of different experiment
                names that you would like to be returned
            show (bool): False, indicates whether or not to print the results to the console 

        Returns:
            dict(): labels = experiment names, values = (mean, std)

        Raises:
            ValueError: if exp_names conatins the water reference        
        """
        if 'water' in exp_names or 'Water' in exp_names:
            raise(ValueError("no value for attenuation of water"))
        means = []
        stds = []
        for label in exp_names:
            water = sum([1 if self.data[k][6] == 'None' else 0 for k in self.data.keys()])
            mask = [True if self.data[k][0] == label else False for k in self.data.keys()][:-water]      
            alphas = self.S[:, -1]
            alphas = alphas[mask]
            means.append(np.mean(alphas))
            stds.append(np.std(alphas))
            if show: 
                print('\n', label)
                print("std: ", stds[-1])
                print("mean: ", means[-1])
        return {exp_names[i]: (means[i], stds[i]) for i in range(len(exp_names)) }
                



"""unused code """   
# plot the power of the sample as a function of the power of water
        # err = []
        # for i, k, l, w in zip(
        #     range(len(self.samples)), range(0, len(self.S), self.t*5), self.samples, list(self.widths.values())):
        #     j = k+self.t*5
        #     plt.subplot(3,3,i+7)
        #     domain = [i for i in list(self.pw.values()) for _ in range(3)]
        #     func = np.multiply(domain,np.exp(S[k:j, -1]*w))
        #     plt.plot(domain, func, label="alpha")
        #     plt.plot(domain, S[k:j, -2], label="mean")
        #     err.append(np.linalg.norm(func - S[k:j, -2]))
        #     plt.legend(loc="upper left")
        #     plt.xlabel(r"$P_{w}$")
        #     plt.ylabel(r"$P_{s}$")
        #     plt.title(r'$P_{s}=P_{w}e^{α*d}$' + f" for {l}")