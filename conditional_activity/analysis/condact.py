'''
CONDACT --- :mod:`CONDACT.analysis.CONDACT`
===========================================================

This module contains the :class:`CONDACT` class.

This code calculates the Conditional activity (time-correlated transitions) of a degree of freedom.\
    In this case, the first sidechain dihedral angles (chi1)  for selected amino residues in protein\
        (except ALA and GLY) and the sugar-phosphate backbone dihedral angle in DNA were used. The code\
            seeks to find the kinetic correlation of amino acids side chains in the 3-dimensional native state\
                of a protein and the sugar-phosphate backbone in a DNA strand. This module was built on MDAnalysis\
                    as a foundation using some functions in the MDAnalysis package.

input: The residue id and names of amino acids and nitrogenous bases of interest.\
        Number of states adopted by selected residues dihedral angle.\
        
output: csv files containing diheral angles of interest,\
        Transition and wait times of each residue,\
        Inter-residue conditional activity of the residues of interest,\
        Dynamical memory of each residue of interest,\
        heatmap of inter-residue conditional activity.
'''
from MDAnalysis.core.universe import Universe, AtomGroup
import numpy as np
from bisect import bisect
import pickle
import math
from json import dumps, load
from bisect import bisect_right
import time

begin = time.time()

States = dict()
transition_times = dict()
wait_times = dict()
OUTFILE = 'Dihedral_Angle.csv'


# from typing import Union, TYPE_CHECKING
from MDAnalysis.analysis.base import AnalysisBase
# if TYPE_CHECKING:
#     from MDAnalysis.core.universe import Universe, AtomGroup

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
        No_of_peaks_DNA :  int
            The number of dihedral states identified for the DNA residue of interest.
            the default is 3 in which two states were similar (see "states_DNA" below)
        peak_boundaries_DNA : List
            A list of boundaries separating one state from the next.
            The number of items in the list is one less than the number of states identified.
            eg. For 3 state, the will be 2 boundaries. default is [90,270]
            Note: The secondary angle of negative dihedrals angles [-180,0] are used i.e 𝛉 = (𝛉+360)
        states_DNA: str
            Each character in the string represents the state of each degree of freedom eg 'SAS' represent state S,\
            A and S. The S state is the syn conformation while the A state is the anti-conformation. Based on our probabilty densities\
            of the states, state S had a range 0<S<90 and 270<S<360 while state A had range 90<A<270.
        saving_frequency : int,
            The saving frequency of the simulation in picoseconds. This is NOT the time-step
            used to run the simulation but the time interval for saving the coordinates in 
            the trajectory.
        keep_negative: boolean
            function to convert all negative dihedrals to the positive secondary corresponding
            angles. This makes diheadral range to be from 0 to 360 rather than -180 to 180.

        Note
        ----

        There is no side chain dihedral for ALA and GLY because the lack CG, OG1, CG1, OG, SG

        .. versionadded:: 2.0.0
            Added 'No_of_peak_protein', 'peak_boundaries_protein' 'No_of_peak_DNA', 'peak_boundaries_DNA' keyword
        """

    def __init__(self, universe,
                    selected_resid=None, 
                    No_of_peaks_protein=3,
                    peak_boundaries_protein=[120,240],
                    states_protein="XYZ",
                    No_of_peaks_DNA=3,
                    peak_boundaries_DNA=[140,330],
                    states_DNA="SAS",
                    saving_frequency=None,
                    keep_negative=False):
        
        self.u = universe
        self.selected_resid = selected_resid
        self.No_of_peaks_protein = No_of_peaks_protein
        self.peak_boundaries_protein = peak_boundaries_protein
        self.states_protein = states_protein
        self.No_of_peaks_DNA = No_of_peaks_DNA
        self.peak_boundaries_DNA = peak_boundaries_DNA
        self.states_DNA = states_DNA
        self.saving_frequency = saving_frequency
        self.keep_negative = keep_negative

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # defining different states
    def define_state_protein(self, angle):
        # defining different states
        if angle < 0:
            raise ValueError ("Angle must be positive because secondary angles used")
        i=bisect(self.peak_boundaries_protein, angle)
        return self.states_protein[i]
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # defining different states
    def define_state_DNA(self, angle):
        # defining different states
        if angle < 0:
            raise ValueError ("Angle must be positive because secondary angles used")
        i=bisect(self.peak_boundaries_DNA, angle)
        return self.states_DNA[i]
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    def find_gt(self, a, x):
        # Bisection search through transition times to calculate adjacent transitions
        'Find leftmost value greater than x'
        i = bisect_right(a, x)
        if i != len(a):
            return a[i]
        else:
            return 0
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    def residues(self):
        # Selecting residues of interest to compute dihedral angle 
        if self.selected_resid is None:
            raise ValueError("please select residues of interest in format a-b or a b ...")
        if self.No_of_peaks_protein - 1 != len(self.peak_boundaries_protein):
            raise ValueError(f"If there are {self.No_of_peaks_protein} peaks, the list of boundaries separating the amino acid peaks in protein in should have {self.No_of_peaks_protein - 1} numbers. For example, if number of peaks is 3, the peak boundaries can be the list of 2 elements i.e [120,240]")
        if self.No_of_peaks_DNA - 1 != len(self.peak_boundaries_DNA):
            raise ValueError(f"If there are {self.No_of_peaks_DNA} peaks, the list of boundaries separating the DNA peaks in nucleotides in should have {self.No_of_peaks_DNA - 1} numbers. For example, if number of peaks is 3, the peak boundaries can be the list of 2 elements i.e [120,240]")

        interested = universe.select_atoms(f'resid {self.selected_resid} and not(resname ALA or resname GLY)')
        residue_type = interested.residues
        residue_names = list(residue_type.resnames)
        return interested, residue_names, residue_type
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    def run(self):

        # Run the residue function---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        interested, residue_names, residue_type = self.residues()
        with open(OUTFILE, 'w') as output:
            # create header in ouput file ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            header = "frame, time(ps)"
            for i in range(len(residue_type)):
                header += f",{residue_names[i]} {residue_type[i].resid}"
                transition_times[f"{residue_names[i]} {residue_type[i].resid}"] = []
                wait_times[f"{residue_names[i]} {residue_type[i].resid}"] = []
                States[f"{residue_names[i]} {residue_type[i].resid}"] = -1
            output.write(header)
             # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            for ts in universe.trajectory: #[0:500]
                """calculate current time in picoseconds and add frame number and current time to output file---"""
                current_time = self.saving_frequency * universe.trajectory.frame
                output.write(f'\n {universe.trajectory.frame}, {current_time}')
                """calculating dihedrals for each residue"""
                for residue in residue_type:
                    residue_names = residue.resname
                #get atoms from the current residue to calculate the base-sugar dihedral angles 𝛉 in DNA #---------------------------------------------------------------------------------------------
                    if residue_names == "DA5":
                    	O4 = universe.select_atoms(f"resid {residue.resid} and (name O4')")
                    	C1 = universe.select_atoms(f"resid {residue.resid} and (name C1')")
                    	N9 = universe.select_atoms(f"resid {residue.resid} and (name N9)")
                    	C4 = universe.select_atoms(f"resid {residue.resid} and (name C4)")
                    	angle = (sum([O4, C1, N9, C4])).dihedral.value()
                    	if not self.keep_negative:
                            angle = angle if angle >= 0 else angle + 360
                    	output.write(f', { angle }')

                    	if States[f"{residue.resname} {residue.resid}"] != self.define_state_DNA(angle):
                            prev = transition_times[f"{residue.resname} {residue.resid}"][-1] if wait_times[f"{residue.resname} {residue.resid}"] != [] else 0
                            transition_times[f"{residue.resname} {residue.resid}"].append(universe.trajectory.frame * self.saving_frequency)
                            # The waiting time W(X) is the time between two consecutive transition times i.e W(X)=T2 - T1
                            if universe.trajectory.frame != 0:
                                wait_times[f"{residue.resname} {residue.resid}"].append(universe.trajectory.frame * self.saving_frequency - prev)
                            States[f"{residue.resname} {residue.resid}"] = self.define_state_DNA(angle)
                    # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                    elif residue_names == "DG5":
                    	O4 = universe.select_atoms(f"resid {residue.resid} and (name O4')")
                    	C1 = universe.select_atoms(f"resid {residue.resid} and (name C1')")
                    	N9 = universe.select_atoms(f"resid {residue.resid} and (name N9)")
                    	C4 = universe.select_atoms(f"resid {residue.resid} and (name C4)")
                    	angle = (sum([O4, C1, N9, C4])).dihedral.value()
                    	if not self.keep_negative:
                            angle = angle if angle >= 0 else angle + 360
                    	output.write(f', { angle }')

                    	if States[f"{residue.resname} {residue.resid}"] != self.define_state_DNA(angle):
                            prev = transition_times[f"{residue.resname} {residue.resid}"][-1] if wait_times[f"{residue.resname} {residue.resid}"] != [] else 0
                            transition_times[f"{residue.resname} {residue.resid}"].append(universe.trajectory.frame * self.saving_frequency)
                            # The waiting time W(X) is the time between two consecutive transition times i.e W(X)=T2 - T1
                            if universe.trajectory.frame != 0:
                                wait_times[f"{residue.resname} {residue.resid}"].append(universe.trajectory.frame * self.saving_frequency - prev)
                            States[f"{residue.resname} {residue.resid}"] = self.define_state_DNA(angle)
                    # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                    elif residue_names == "DT5":
                    	O4 = universe.select_atoms(f"resid {residue.resid} and (name O4')")
                    	C1 = universe.select_atoms(f"resid {residue.resid} and (name C1')")
                    	N1 = universe.select_atoms(f"resid {residue.resid} and (name N1)")
                    	C2 = universe.select_atoms(f"resid {residue.resid} and (name C2)")
                    	angle = (sum([O4, C1, N1, C2])).dihedral.value()
                    	if not self.keep_negative:
                            angle = angle if angle >= 0 else angle + 360
                    	output.write(f', { angle }')

                    	if States[f"{residue.resname} {residue.resid}"] != self.define_state_DNA(angle):
                            prev = transition_times[f"{residue.resname} {residue.resid}"][-1] if wait_times[f"{residue.resname} {residue.resid}"] != [] else 0
                            transition_times[f"{residue.resname} {residue.resid}"].append(universe.trajectory.frame * self.saving_frequency)
                            # The waiting time W(X) is the time between two consecutive transition times i.e W(X)=T2 - T1
                            if universe.trajectory.frame != 0:
                                wait_times[f"{residue.resname} {residue.resid}"].append(universe.trajectory.frame * self.saving_frequency - prev)
                            States[f"{residue.resname} {residue.resid}"] = self.define_state_DNA(angle)
                    # ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                    elif residue_names == "DC5":
                    	O4 = universe.select_atoms(f"resid {residue.resid} and (name O4')")
                    	C1 = universe.select_atoms(f"resid {residue.resid} and (name C1')")
                    	N1 = universe.select_atoms(f"resid {residue.resid} and (name N1)")
                    	C2 = universe.select_atoms(f"resid {residue.resid} and (name C2)")
                    	angle = (sum([O4, C1, N1, C2])).dihedral.value()
                    	if not self.keep_negative:
                            angle = angle if angle >= 0 else angle + 360
                    	output.write(f', { angle }')

                    	if States[f"{residue.resname} {residue.resid}"] != self.define_state_DNA(angle):
                            prev = transition_times[f"{residue.resname} {residue.resid}"][-1] if wait_times[f"{residue.resname} {residue.resid}"] != [] else 0
                            transition_times[f"{residue.resname} {residue.resid}"].append(universe.trajectory.frame * self.saving_frequency)
                            # The waiting time W(X) is the time between two consecutive transition times i.e W(X)=T2 - T1
                            if universe.trajectory.frame != 0:
                                wait_times[f"{residue.resname} {residue.resid}"].append(universe.trajectory.frame * self.saving_frequency - prev)
                            States[f"{residue.resname} {residue.resid}"] = self.define_state_DNA(angle)
                    # --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                    elif residue_names == "DA":
                    	O4 = universe.select_atoms(f"resid {residue.resid} and (name O4')")
                    	C1 = universe.select_atoms(f"resid {residue.resid} and (name C1')")
                    	N9 = universe.select_atoms(f"resid {residue.resid} and (name N9)")
                    	C4 = universe.select_atoms(f"resid {residue.resid} and (name C4)")
                    	angle = (sum([O4, C1, N9, C4])).dihedral.value()
                    	if not self.keep_negative:
                            angle = angle if angle >= 0 else angle + 360
                    	output.write(f', { angle }')

                    	if States[f"{residue.resname} {residue.resid}"] != self.define_state_DNA(angle):
                            prev = transition_times[f"{residue.resname} {residue.resid}"][-1] if wait_times[f"{residue.resname} {residue.resid}"] != [] else 0
                            transition_times[f"{residue.resname} {residue.resid}"].append(universe.trajectory.frame * self.saving_frequency)
                            # The waiting time W(X) is the time between two consecutive transition times i.e W(X)=T2 - T1
                            if universe.trajectory.frame != 0:
                                wait_times[f"{residue.resname} {residue.resid}"].append(universe.trajectory.frame * self.saving_frequency - prev)
                            States[f"{residue.resname} {residue.resid}"] = self.define_state_DNA(angle)
                    # ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                    elif residue_names == "DG":
                    	O4 = universe.select_atoms(f"resid {residue.resid} and (name O4')")
                    	C1 = universe.select_atoms(f"resid {residue.resid} and (name C1')")
                    	N9 = universe.select_atoms(f"resid {residue.resid} and (name N9)")
                    	C4 = universe.select_atoms(f"resid {residue.resid} and (name C4)")
                    	angle = (sum([O4, C1, N9, C4])).dihedral.value()
                    	if not self.keep_negative:
                            angle = angle if angle >= 0 else angle + 360
                    	output.write(f', { angle }')

                    	if States[f"{residue.resname} {residue.resid}"] != self.define_state_DNA(angle):
                            prev = transition_times[f"{residue.resname} {residue.resid}"][-1] if wait_times[f"{residue.resname} {residue.resid}"] != [] else 0
                            transition_times[f"{residue.resname} {residue.resid}"].append(universe.trajectory.frame * self.saving_frequency)
                            # The waiting time W(X) is the time between two consecutive transition times i.e W(X)=T2 - T1
                            if universe.trajectory.frame != 0:
                                wait_times[f"{residue.resname} {residue.resid}"].append(universe.trajectory.frame * self.saving_frequency - prev)
                            States[f"{residue.resname} {residue.resid}"] = self.define_state_DNA(angle)
                    # -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                    elif residue_names == "DT":
                    	O4 = universe.select_atoms(f"resid {residue.resid} and (name O4')")
                    	C1 = universe.select_atoms(f"resid {residue.resid} and (name C1')")
                    	N1 = universe.select_atoms(f"resid {residue.resid} and (name N1)")
                    	C2 = universe.select_atoms(f"resid {residue.resid} and (name C2)")
                    	angle = (sum([O4, C1, N1, C2])).dihedral.value()
                    	if not self.keep_negative:
                            angle = angle if angle >= 0 else angle + 360
                    	output.write(f', { angle }')

                    	if States[f"{residue.resname} {residue.resid}"] != self.define_state_DNA(angle):
                            prev = transition_times[f"{residue.resname} {residue.resid}"][-1] if wait_times[f"{residue.resname} {residue.resid}"] != [] else 0
                            transition_times[f"{residue.resname} {residue.resid}"].append(universe.trajectory.frame * self.saving_frequency)
                            # The waiting time W(X) is the time between two consecutive transition times i.e W(X)=T2 - T1
                            if universe.trajectory.frame != 0:
                                wait_times[f"{residue.resname} {residue.resid}"].append(universe.trajectory.frame * self.saving_frequency - prev)
                            States[f"{residue.resname} {residue.resid}"] = self.define_state_DNA(angle)
                    # -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                    elif residue_names == "DC":
                    	O4 = universe.select_atoms(f"resid {residue.resid} and (name O4')")
                    	C1 = universe.select_atoms(f"resid {residue.resid} and (name C1')")
                    	N1 = universe.select_atoms(f"resid {residue.resid} and (name N1)")
                    	C2 = universe.select_atoms(f"resid {residue.resid} and (name C2)")
                    	angle = (sum([O4, C1, N1, C2])).dihedral.value()
                    	if not self.keep_negative:
                            angle = angle if angle >= 0 else angle + 360
                    	output.write(f', { angle }')

                    	if States[f"{residue.resname} {residue.resid}"] != self.define_state_DNA(angle):
                            prev = transition_times[f"{residue.resname} {residue.resid}"][-1] if wait_times[f"{residue.resname} {residue.resid}"] != [] else 0
                            transition_times[f"{residue.resname} {residue.resid}"].append(universe.trajectory.frame * self.saving_frequency)
                            # The waiting time W(X) is the time between two consecutive transition times i.e W(X)=T2 - T1
                            if universe.trajectory.frame != 0:
                                wait_times[f"{residue.resname} {residue.resid}"].append(universe.trajectory.frame * self.saving_frequency - prev)
                            States[f"{residue.resname} {residue.resid}"] = self.define_state_DNA(angle)
                    # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                    elif residue_names == "DA3":
                     	O4 = universe.select_atoms(f"resid {residue.resid} and (name O4')")
                     	C1 = universe.select_atoms(f"resid {residue.resid} and (name C1')")
                     	N9 = universe.select_atoms(f"resid {residue.resid} and (name N9)")
                     	C4 = universe.select_atoms(f"resid {residue.resid} and (name C4)")
                     	angle = (sum([O4, C1, N9, C4])).dihedral.value()
                     	if not self.keep_negative:
                     	    angle = angle if angle >= 0 else angle + 360
                     	output.write(f', { angle }')
                     	
                     	if States[f"{residue.resname} {residue.resid}"] != self.define_state_DNA(angle):
                     	    prev = transition_times[f"{residue.resname} {residue.resid}"][-1] if wait_times[f"{residue.resname} {residue.resid}"] != [] else 0
                     	    transition_times[f"{residue.resname} {residue.resid}"].append(universe.trajectory.frame * self.saving_frequency)
                     	    # The waiting time W(X) is the time between two consecutive transition times i.e W(X)=T2 - T1
                     	    if universe.trajectory.frame != 0:
                     	        wait_times[f"{residue.resname} {residue.resid}"].append(universe.trajectory.frame * self.saving_frequency - prev)
                     	    States[f"{residue.resname} {residue.resid}"] = self.define_state_DNA(angle)
                    # ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                    elif residue_names == "DG3":
                    	O4 = universe.select_atoms(f"resid {residue.resid} and (name O4')")
                    	C1 = universe.select_atoms(f"resid {residue.resid} and (name C1')")
                    	N9 = universe.select_atoms(f"resid {residue.resid} and (name N9)")
                    	C4 = universe.select_atoms(f"resid {residue.resid} and (name C4)")
                    	angle = (sum([O4, C1, N9, C4])).dihedral.value()
                    	if not self.keep_negative:
                            angle = angle if angle >= 0 else angle + 360
                    	output.write(f', { angle }')

                    	if States[f"{residue.resname} {residue.resid}"] != self.define_state_DNA(angle):
                            prev = transition_times[f"{residue.resname} {residue.resid}"][-1] if wait_times[f"{residue.resname} {residue.resid}"] != [] else 0
                            transition_times[f"{residue.resname} {residue.resid}"].append(universe.trajectory.frame * self.saving_frequency)
                            # The waiting time W(X) is the time between two consecutive transition times i.e W(X)=T2 - T1
                            if universe.trajectory.frame != 0:
                                wait_times[f"{residue.resname} {residue.resid}"].append(universe.trajectory.frame * self.saving_frequency - prev)
                            States[f"{residue.resname} {residue.resid}"] = self.define_state_DNA(angle)
                    # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                    elif residue_names == "DT3":
                    	O4 = universe.select_atoms(f"resid {residue.resid} and (name O4')")
                    	C1 = universe.select_atoms(f"resid {residue.resid} and (name C1')")
                    	N1 = universe.select_atoms(f"resid {residue.resid} and (name N1)")
                    	C2 = universe.select_atoms(f"resid {residue.resid} and (name C2)")
                    	angle = (sum([O4, C1, N1, C2])).dihedral.value()
                    	if not self.keep_negative:
                            angle = angle if angle >= 0 else angle + 360
                    	output.write(f', { angle }')

                    	if States[f"{residue.resname} {residue.resid}"] != self.define_state_DNA(angle):
                            prev = transition_times[f"{residue.resname} {residue.resid}"][-1] if wait_times[f"{residue.resname} {residue.resid}"] != [] else 0
                            transition_times[f"{residue.resname} {residue.resid}"].append(universe.trajectory.frame * self.saving_frequency)
                            # The waiting time W(X) is the time between two consecutive transition times i.e W(X)=T2 - T1
                            if universe.trajectory.frame != 0:
                                wait_times[f"{residue.resname} {residue.resid}"].append(universe.trajectory.frame * self.saving_frequency - prev)
                            States[f"{residue.resname} {residue.resid}"] = self.define_state_DNA(angle)
                    # -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                    elif residue_names == "DC3":
                    	O4 = universe.select_atoms(f"resid {residue.resid} and (name O4')")
                    	C1 = universe.select_atoms(f"resid {residue.resid} and (name C1')")
                    	N1 = universe.select_atoms(f"resid {residue.resid} and (name N1)")
                    	C2 = universe.select_atoms(f"resid {residue.resid} and (name C2)")
                    	angle = (sum([O4, C1, N1, C2])).dihedral.value()
                    	if not self.keep_negative:
                            angle = angle if angle >= 0 else angle + 360
                    	output.write(f', { angle }')

                    	if States[f"{residue.resname} {residue.resid}"] != self.define_state_DNA(angle):
                            prev = transition_times[f"{residue.resname} {residue.resid}"][-1] if wait_times[f"{residue.resname} {residue.resid}"] != [] else 0
                            transition_times[f"{residue.resname} {residue.resid}"].append(universe.trajectory.frame * self.saving_frequency)
                            # The waiting time W(X) is the time between two consecutive transition times i.e W(X)=T2 - T1
                            if universe.trajectory.frame != 0:
                                wait_times[f"{residue.resname} {residue.resid}"].append(universe.trajectory.frame * self.saving_frequency - prev)
                            States[f"{residue.resname} {residue.resid}"] = self.define_state_DNA(angle)

                # get proper atoms from the current residue for calculating dihedral angles chi1 (𝞆1) in amino acids---------------------------------------------------------------------------------------------------------------------            
                    elif residue_names not in ["DA5", "DG5", "DT5", "DC5", "DA", "DG", "DT", "DC","DA3", "DG3", "DT3", "DC3"]:
                        N = interested.select_atoms(f'resid {residue.resid} and name N')
                        CA = interested.select_atoms(f'resid {residue.resid} and name CA')
                        CB = interested.select_atoms(f'resid {residue.resid} and name CB')
                        CG = interested.select_atoms(f'resid {residue.resid} and (name CG or name OG1 or name CG1 or name OG or name SG)')
                        angle = (sum([N, CA, CB, CG])).dihedral.value()
                        if not self.keep_negative:
                            angle = angle if angle >= 0 else angle + 360
                        output.write(f', { angle }')

                        if States[f"{residue.resname} {residue.resid}"] != self.define_state_protein(angle):
                            prev = transition_times[f"{residue.resname} {residue.resid}"][-1] if wait_times[f"{residue.resname} {residue.resid}"] != [] else 0
                            transition_times[f"{residue.resname} {residue.resid}"].append(universe.trajectory.frame * self.saving_frequency)
                            # The waiting time W(X) is the time between two consecutive transition times i.e W(X)=T2 - T1
                            if universe.trajectory.frame != 0:
                                wait_times[f"{residue.resname} {residue.resid}"].append(universe.trajectory.frame * self.saving_frequency - prev)
                            States[f"{residue.resname} {residue.resid}"] = self.define_state_protein(angle)
        
        with open('transition_times.pkl', 'wb') as f:
            pickle.dump(transition_times, f)
        with open('wait_times.pkl', 'wb') as t:
            pickle.dump(wait_times, t)
        # print(transition_times)
        # print(wait_times)
        return transition_times, wait_times
        

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    def mean_persistent_time(self):
        # Run the residue and loop function----------------------------------------------------------------------------------------------------------------------------------------
        with open('wait_times.pkl', 'rb') as f:
            wait_times= pickle.load(f)
        interested, residue_names, residue_type = self.residues()
        observation_time = len(universe.trajectory) * self.saving_frequency
        persistence_times = dict()
        list_persistence_times = []
        '''Calculating the persistent time for each residue. The persistence time of each amino acid is the\
            average time in which the amino acid's degree of freedom (in this case dihedral angle) remains\
            in a specific state. It uses the waiting time W(X) of each amino acids as shown below:
            Tp[X] = 1/2𝜏 * sum{1,N(X)} W(X, T(X, i))^2'''
            
        for res in wait_times.values():
            # Ignoring all residues with less than 10 transitions, hence less than 9 wait times
            res = res[1:]
            if len(res) < 9:
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
         # Call over trajectory to get transitions --------------------------------------------------------
        with open('transition_times.pkl', 'rb') as f:
            transition_times= pickle.load(f)
        observation_time = len(universe.trajectory) * self.saving_frequency
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
                for i in range(2, len(transition_times[res1])):
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
                    print(f"{residue_names[res]} {residue_type[res].resid} or {residue_names[time_]} {residue_type[time_].resid} had statistically fewer transitions hence Conditional Activity CA[{residue_names[res]} {residue_type[res].resid}][{residue_names[time_]} {residue_type[time_].resid}] could not be determine.")
                    CA = math.nan
                    conditional_activity_.append(CA)
                else:
                    CA = -np.log(exchange_times[res][time_] / list_persistence_times[res])
                    conditional_activity_.append(CA if (CA != np.inf and CA != -np.inf) else (math.nan))
            conditional_activity.append(conditional_activity_)
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        # Extract dynamic memory from conditional activity matrix -----------------------------------------------------------------------------
        Dynamic_Memory = np.diag(np.array(conditional_activity))
        np.savetxt('Dynamic_Memory.txt', Dynamic_Memory, fmt='%.7f')

        # Plot heatmap for conditional activity matrix ----------------------------------------------------------------------------------------------
        import matplotlib.pyplot as plt
        fig, ax=plt.subplots(figsize=(11, 8))
        orig_map=plt.get_cmap('hot')
        reversed_map = orig_map.reversed()
        im=ax.imshow(conditional_activity,interpolation='nearest', origin='lower', cmap = reversed_map, aspect='auto', vmin=0)
        cbar2 = fig.colorbar(im)
        cbar2.ax.set_ylabel('$CA$', fontsize = 10.0)
        plt.savefig("Conditional_Activity.png",format='png', dpi=300)

        # Write conditional activity matrix from results and save to text file---------------------------------------------------------------------
        conditional_activity_matrix = np.array(conditional_activity)
        print(conditional_activity_matrix)
        np.savetxt('Conditional_Activity.txt', conditional_activity_matrix, fmt='%.7f', delimiter=' ')
        return conditional_activity, persistence_times, exchange_times
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# universe = Universe("/Users/augustineonyema/MolSSI/Conditional_Activity/Codes_Cond_Act/dry_Lys.prmtop", 
#                     "/Users/augustineonyema/MolSSI/Conditional_Activity/Codes_Cond_Act/test_dry_Lys.xtc")
universe = Universe("pytest_LYS.prmtop", 
                    "pytest_LYS.xtc")

Study = CONDACT(universe,
                    selected_resid='1-3', #1-1268
                    No_of_peaks_protein=3,
                    peak_boundaries_protein=[120,240],
                    states_protein = "XYZ",
                    peak_boundaries_DNA =[140,330],
                    No_of_peaks_DNA=3,
                    states_DNA = "SAS",
                    saving_frequency=10,
                    keep_negative=False) # 33 34 35 36 41 46 52 53 58 59 63 83 98 101 107 108 109
conditional_activity, persistence_times, exchange_times = Study.mean_conditional_activity()

print(f"persistence_times{persistence_times}")
print(f"exchange_times{exchange_times}")
print(conditional_activity)

# total time taken 
end = time.time() 
# print(f"Total runtime of the program is {end - begin}") 
