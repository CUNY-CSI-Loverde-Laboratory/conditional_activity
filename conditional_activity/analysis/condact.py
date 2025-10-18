'''
CONDACT --- :mod:`CONDACT.analysis.CONDACT`
===========================================================

This module contains the :class:`CONDACT` class.

This code calculates the Conditional activity (time-correlated transitions) of a degree of freedom.\
    In this case, the first sidechain dihedral angles (chi1)  for selected amino residues in protein\
        (except ALA and GLY) and the sugar-base dihedral angle in nucleotides were used. The code\
            seeks to find the kinetic correlation of amino acids side chains in the 3-dimensional native state\
                of a protein and the sugar-base in a polynucleotide strand. This module was built on MDAnalysis\
                    as a foundation using some functions in the MDAnalysis package.

input: The residue id and names of amino acids and nitrogenous bases of interest.\
        Number of states adopted by selected residues dihedral angle.\
        
output: csv files containing diheral angles of interest,\
        Transition and wait times of each residue in binary file .pkl,\
        Inter-residue conditional activity of the residues of interest in binary file .npy,\
        Dynamical memory of each residue of interest in a txt file,\
        heatmap of inter-residue conditional activity in .png.
'''
from MDAnalysis.core.universe import Universe, AtomGroup
from MDAnalysis.analysis.base import AnalysisBase
import numpy as np
import matplotlib.pyplot as plt
from bisect import bisect, bisect_right
import pickle
import math
from json import dumps, load
import time


begin = time.time()

class CONDACT(AnalysisBase):
    """CONDACT class.
    Performs the Conditional Activity of the dihedral angles on the universe. Trajectory must be long with residue having at least 10 transitions.

        Parameters
        ----------
        universe : Universe
            MDAnalysis Universe object
        selected_resid : str
            Selection string for the residues of interest one 
            wishes to calculate the conditional activity on.
        No_of_peaks_protein :  int
            The number of dihedral states identified for the amino acids of interest, 
            the default is 3
        peak_boundaries_protein : List
            A list of boundaries separating one state from the next.
            The number of items in the list is one less than the number of states identified.
            eg. For 3 state, the will be 2 boundaries. default is [120,240]
            Note: The secondary angle of negative dihedrals angles [-180,0] are used i.e 𝛉 = (𝛉+360)
        states_protein: str
            Each character in the string represents the state of each degree of freedom eg 'XYZ' represent state X,\
            Y and Z. The default is 'XYZ' in which state X has range 0<X<120, state Y has range 120<Y<240, \
            and state Z has range 240<Z<360.
        No_of_peaks_nucleic :  int
            The number of dihedral states identified for the nucleotides of interest.
            the default is 3 in which two states were similar (see "states_nucleic" below)
        peak_boundaries_nucleic : List
            A list of boundaries separating one state from the next.
            The number of items in the list is one less than the number of states identified.
            eg. For 3 state, the will be 2 boundaries. default is [140,330]
            Note: The secondary angle of negative dihedrals angles [-180,0] are used i.e 𝛉 = (𝛉+360)
        states_nucleic: str
            Each character in the string represents the state of each degree of freedom eg 'SAS' represent state S,\
            A and S. The S state is the syn conformation while the A state is the anti-conformation. Based on our probabilty densities\
            of the states, state S had a range 0<S<140 and 330<S<360 while state A had range 140<A<330.
        saving_frequency : int,
            The saving frequency of the simulation in picoseconds. This is NOT the time-step
            used to run the simulation but the time interval for saving the coordinates in 
            the trajectory.
        keep_negative: boolean
            function to convert all negative dihedrals to the positive secondary corresponding
            angles. This makes diheadral range to be from 0 to 360 rather than -180 to 180.

        Note
        ----

        There is no sidechain dihedral for ALA and GLY because the lack CG, OG1, CG1, OG, SG

        .. versionadded:: 2.0.0
            Added 'No_of_peak_protein', 'peak_boundaries_protein' 'No_of_peak_nucleic', 'peak_boundaries_nucleic' keyword
        """

    def __init__(self, universe=None,
                    selected_resid=None, 
                    No_of_peaks_protein=3,
                    peak_boundaries_protein=[120,240],
                    states_protein="XYZ",
                    No_of_peaks_nucleic=3,
                    peak_boundaries_nucleic=[140,330],
                    states_nucleic="SAS",
                    saving_frequency=None,
                    keep_negative=False):
        
        self.universe = universe
        self.selected_resid = selected_resid
        self.No_of_peaks_protein = No_of_peaks_protein
        self.peak_boundaries_protein = peak_boundaries_protein
        self.states_protein = states_protein
        self.No_of_peaks_nucleic = No_of_peaks_nucleic
        self.peak_boundaries_nucleic = peak_boundaries_nucleic
        self.states_nucleic = states_nucleic
        self.saving_frequency = saving_frequency
        self.keep_negative = keep_negative

        self.states = dict()
        self.transition_times = dict()
        self.wait_times = dict()
        
        # Dictionary of all selected nucleic acid residues of interest and their sugar-base dihedral angle residues
        self.dihedral_atom_map_nucleic = {
            # DNA
            "DA": ["O4'", "C1'", "N9", "C4"],
            "DG": ["O4'", "C1'", "N9", "C4"],
            "DT": ["O4'", "C1'", "N1", "C2"],
            "DC": ["O4'", "C1'", "N1", "C2"],
            "DA5": ["O4'", "C1'", "N9", "C4"],
            "DG5": ["O4'", "C1'", "N9", "C4"],
            "DT5": ["O4'", "C1'", "N1", "C2"],
            "DC5": ["O4'", "C1'", "N1", "C2"],
            "DA3": ["O4'", "C1'", "N9", "C4"],
            "DG3": ["O4'", "C1'", "N9", "C4"],
            "DT3": ["O4'", "C1'", "N1", "C2"],
            "DC3": ["O4'", "C1'", "N1", "C2"],
            
            # RNA
            "A":  ["O4'", "C1'", "N9", "C4"],
            "G":  ["O4'", "C1'", "N9", "C4"],
            "U":  ["O4'", "C1'", "N1", "C2"],
            "C":  ["O4'", "C1'", "N1", "C2"],
            "A5": ["O4'", "C1'", "N9", "C4"],
            "G5": ["O4'", "C1'", "N9", "C4"],
            "U5": ["O4'", "C1'", "N1", "C2"],
            "C5": ["O4'", "C1'", "N1", "C2"],
            "A3": ["O4'", "C1'", "N9", "C4"],
            "G3": ["O4'", "C1'", "N9", "C4"],
            "U3": ["O4'", "C1'", "N1", "C2"],
            "C3": ["O4'", "C1'", "N1", "C2"]
        }

    def _check_universe(self):
        """Ensure that the MDAnalysis Universe object is provided."""
        if self.universe is None:
            raise ValueError("MDAnalysis Universe object must be provided.")

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # defining different states of the amino acids in the protein
    def define_state_protein(self, angle):
        # defining different states
        if angle < 0:
            raise ValueError ("Angle must be positive because secondary angles used")
        i=bisect(self.peak_boundaries_protein, angle)
        return self.states_protein[i]
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # defining different states of the nucleotides in the DNA and RNA
    def define_state_nucleic(self, angle):
        # defining different states
        if angle < 0:
            raise ValueError ("Angle must be positive because secondary angles used")
        i=bisect(self.peak_boundaries_nucleic, angle)
        return self.states_nucleic[i]
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    def find_gt(self, a, x):
        # Bisection search through transition times to calculate adjacent transitions
        'Find leftmost value greater than x'
        i = bisect_right(a, x)
        return a[i] if i != len(a) else 0
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    def residues(self):
        # Selecting residues of interest to compute dihedral angle 
        if self.selected_resid is None:
            raise ValueError("please select residues of interest in format a-b or a b ...")
        if self.No_of_peaks_protein - 1 != len(self.peak_boundaries_protein):
            raise ValueError(f"If there are {self.No_of_peaks_protein} peaks, the list of boundaries separating the amino acid peaks in protein in should have {self.No_of_peaks_protein - 1} numbers. For example, if number of peaks is 3, the peak boundaries can be the list of 2 elements i.e [120,240]")
        if self.No_of_peaks_nucleic - 1 != len(self.peak_boundaries_nucleic):
            raise ValueError(f"If there are {self.No_of_peaks_nucleic} peaks, the list of boundaries separating the nucleic peaks in nucleotides in should have {self.No_of_peaks_nucleic - 1} numbers. For example, if number of peaks is 3, the peak boundaries can be the list of 2 elements i.e [140,330]")

        interested = self.universe.select_atoms(f'resid {self.selected_resid} and not(resname ALA or resname GLY)')
        residue_type = interested.residues
        residue_names = list(residue_type.resnames)
        return interested, residue_names, residue_type
   
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    def run(self):
        OUTFILE = 'Dihedral_Angle.csv'
        # Run the residue function---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        interested, residue_names, residue_type = self.residues()
        with open(OUTFILE, 'w') as output:
            # create header in ouput file ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            header = "frame, time(ps)"
            for i in range(len(residue_type)):
                resid_key = f"{residue_names[i]} {residue_type[i].resid}"
                header += f",{resid_key}"
                self.transition_times[resid_key] = []
                self.wait_times[resid_key] = []
                self.states[resid_key] = -1
            output.write(header)

             # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            
            for ts in self.universe.trajectory: #[0:500]
                """calculate current time in picoseconds and add frame number and current time to output file---"""
                current_time = self.saving_frequency * self.universe.trajectory.frame
                output.write(f'\n {self.universe.trajectory.frame}, {current_time}')
                """calculating dihedrals for each residue"""

                for residue in residue_type:
                    resid_key = f"{residue.resname} {residue.resid}"
                    resname = residue.resname
                    #get atoms from the current residue to calculate the base-sugar dihedral angles 𝛉 in nucleic #--------------------------------------------------------------------------------------------------------------
                    if resname in self.dihedral_atom_map_nucleic:
                        # Nucleic acids (DNA or RNA)
                        atom_names = self.dihedral_atom_map_nucleic[resname]
                        atoms = [self.universe.select_atoms(f"resid {residue.resid} and name {name}") for name in atom_names]
                        angle = (sum(atoms)).dihedral.value()
                        if not self.keep_negative:
                            angle = angle if angle >= 0 else angle + 360
                        output.write(f', {angle:.2f}')

                        if self.states[resid_key] != self.define_state_nucleic(angle):
                            prev = self.transition_times[resid_key][-1] if self.wait_times[resid_key] else 0
                            self.transition_times[resid_key].append(self.universe.trajectory.frame * self.saving_frequency)
                            if self.universe.trajectory.frame != 0:
                                self.wait_times[resid_key].append(self.universe.trajectory.frame * self.saving_frequency - prev)
                            self.states[resid_key] = self.define_state_nucleic(angle)

                    else:
                        # get proper atoms from the current residue for calculating dihedral angles chi1 (𝞆1) in amino acids---------------------------------------------------------------------------------------------------------------------
                        try:
                            N = interested.select_atoms(f"resid {residue.resid} and name N")
                            CA = interested.select_atoms(f"resid {residue.resid} and name CA")
                            CB = interested.select_atoms(f"resid {residue.resid} and name CB")
                            CG = interested.select_atoms(f"resid {residue.resid} and (name CG or name OG1 or name CG1 or name OG or name SG)")
                            angle = (sum([N, CA, CB, CG])).dihedral.value()
                            if not self.keep_negative:
                                angle = angle if angle >= 0 else angle + 360
                            output.write(f', {angle:.2f}')

                            if self.states[resid_key] != self.define_state_protein(angle):
                                prev = self.transition_times[resid_key][-1] if self.wait_times[resid_key] else 0
                                self.transition_times[resid_key].append(self.universe.trajectory.frame * self.saving_frequency)
                                if self.universe.trajectory.frame != 0:
                                    self.wait_times[resid_key].append(self.universe.trajectory.frame * self.saving_frequency - prev)
                                self.states[resid_key] = self.define_state_protein(angle)
                        except Exception as e:
                            # In case atom group missingt in current residue, write NaN #-----------------------------------------------------------------------------------------------------------------
                            output.write(', NaN')

        with open('transition_times.pkl', 'wb') as f:
            pickle.dump(self.transition_times, f)
        with open('wait_times.pkl', 'wb') as t:
            pickle.dump(self.wait_times, t)
        return self.transition_times, self.wait_times
        
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    def mean_persistent_time(self):
        # Run the residue and loop function----------------------------------------------------------------------------------------------------------------------------------------
        with open('wait_times.pkl', 'rb') as f:
            wait_times= pickle.load(f)
        interested, residue_names, residue_type = self.residues()
        observation_time = len(self.universe.trajectory) * self.saving_frequency
        persistence_times = dict()
        list_persistence_times = []
        '''Calculating the persistent time for each residue. The persistence time of each amino acid is the\
            average time in which the amino acid's degree of freedom (in this case dihedral angle) remains\
            in a specific state. It uses the waiting time W(X) of each amino acids as shown below:
            Tp[X] = 1/2𝜏 * sum{1,N(X)} W(X, T(X, i))^2'''
            
        for res in wait_times.values():
            # Ignoring all residues with less than 10 transitions, hence less than 9 wait times
            if len(res) < 10:
                res = [0]
            # start the summation: wait-time * p(wait-time), and p(wait-time) is just wait-time i.e (wait-time ^ 2)---
            sum_ = 0
            for time in range(len(res)):
                # The sequential addition of the product of the square of each wait-time  i.e [A1*A1] + [A2*A2]+...[An*An] 
                sum_ += res[time] * res[time]
            residue_pers_times = sum_ / (2 * observation_time)
            list_persistence_times.append(residue_pers_times)
        for i in range(len(residue_type)):
            persistence_times.update({f"{residue_names[i]} {residue_type[i].resid}": list_persistence_times[i]})
        return persistence_times, list_persistence_times
    
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    def mean_exchange_time(self):
        '''Calculating the persistent time for each residue. The exchange time of amino acid X with\
            respect to Y is the average time taken for the degree of freedom (dihedral angle)  of Y\
            to change after a change in the degree of freedom of X. 
                It uses the waiting time W(X) of each amino acids as shown below:
            Tp[X] = 1/𝜏 * sum{1,N(Y)-1} W(X,T(Y,i+1))W(Y,T(Y,i))'''
         # Call over trajectory to get transitions ---------------------------------------------------------------
        with open('transition_times.pkl', 'rb') as f:
            transition_times= pickle.load(f)
        observation_time = len(self.universe.trajectory) * self.saving_frequency
        exchange_times = []
        # for loops to iterate keys of dictionary
        for res1 in transition_times:
            # Ignoring all residues with less than 10 transitions
            if len(transition_times[res1]) < 11:
                transition_times[res1] = [0]
            exchange_times_ = []
            for res2 in transition_times:
                # Ignoring all residues with less than 10 transitions
                if len(transition_times[res2]) < 11:
                    transition_times[res2] = [0]
                sum_ = 0 # reset the sum
                
                # for loops to iterate arrays within dictionary
                for i in range(1, len(transition_times[res1])):
                    # ensures we only get the first
                    if (self.find_gt(transition_times[res2], transition_times[res1][i])) < transition_times[res1][i]:
                        break 
                    else:
                        sum_ += ((self.find_gt(transition_times[res2], transition_times[res1][i]) - transition_times[res1][i]) * (transition_times[res1][i] - transition_times[res1][i - 1]))
                exchange_times_.append(sum_ / observation_time)# get avg
            exchange_times.append(exchange_times_) # add to the array
        return exchange_times
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    def mean_conditional_activity(self):
        # using exchange-time / persistent-time to find conditional activiy ------------------------------------
        interested, residue_names, residue_type = self.residues()
        transition_times, wait_times = self.run()
        persistence_times, list_persistence_times = self.mean_persistent_time()
        exchange_times = self.mean_exchange_time()
        conditional_activity = []
        for res in range(len(exchange_times)):
            conditional_activity_ = []
            for time_ in range(len(exchange_times[res])):
                if (exchange_times[res][time_] == 0.0) or (list_persistence_times[res]== 0.0):
                    print(f"{residue_names[res]} {residue_type[res].resid} or {residue_names[time_]} {residue_type[time_].resid} had statistically fewer transitions hence Conditional Activity A[{residue_names[res]} {residue_type[res].resid}][{residue_names[time_]} {residue_type[time_].resid}] could not be determine.")
                    CA = math.nan
                    conditional_activity_.append(CA)
                else:
                    CA = -np.log(exchange_times[res][time_] / list_persistence_times[res])
                    conditional_activity_.append(CA if (CA != np.inf and CA != -np.inf) else (math.nan))
            conditional_activity.append(conditional_activity_)
        np.save("Conditional_Activity_list.npy",conditional_activity)

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        # Extract dynamic memory from conditional activity matrix -----------------------------------------------------------------------------
        Dynamic_Memory = np.diag(np.array(conditional_activity))
        np.save('Dynamic_Memory.npy', Dynamic_Memory)

        # Write conditional activity matrix from results and save to binary file---------------------------------------------------------------------
        conditional_activity_matrix = np.array(conditional_activity)
        np.save('Conditional_Activity_matrix.npy', conditional_activity_matrix)
        return conditional_activity, persistence_times, exchange_times

# # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# # total time taken 
end = time.time()
print(f"Total runtime of the program is {end - begin}") 