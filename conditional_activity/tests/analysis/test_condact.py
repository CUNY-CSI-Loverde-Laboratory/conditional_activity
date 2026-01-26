"""
Unit and regression test for the conditional_activity package.
"""

# Import package, test suite, and other packages as needed
import sys
import os
import pytest
import pickle
from MDAnalysis import Universe
from conditional_activity.analysis.condact import CONDACT
# from conditional_activity.tests.utils import make_Universe
import math


# create universe for pytest
universe = Universe("pytest_LYS.prmtop", "pytest_LYS.xtc")
OUTFILE = 'Dihedral_Angle.csv'

def test_conditional_activity_imported():
    """Test that the module imports correctly."""
    assert "conditional_activity.analysis.condact" in sys.modules


def test_check_universe_raises_value_error():
    # Create CONDACT instance without providing a universe
    condact = CONDACT(universe=None)

    # Expect ValueError when _check_universe is called
    with pytest.raises(ValueError, match="MDAnalysis Universe object must be provided."):
        condact._check_universe()


# Check default values and initialized parameters
def test_initialization_with_defaults():
    condact = CONDACT(universe)
    assert condact.No_of_peaks_protein == 3
    assert condact.peak_boundaries_protein == [120, 240]
    assert condact.states_protein == "XYZ"
    assert condact.No_of_peaks_nucleic == 3
    assert condact.peak_boundaries_nucleic == [140, 330]
    assert condact.states_nucleic == "SAS"
    assert condact.saving_frequency is None
    assert condact.keep_negative is False


# Check custom values
def test_initialization_with_custom_values():
    condact = CONDACT(
        universe,
        selected_resid="1-5",
        No_of_peaks_protein=5,
        peak_boundaries_protein=[50, 150, 250, 310],
        states_protein="ABCDE",
        No_of_peaks_nucleic=4,
        peak_boundaries_nucleic=[100, 200, 300],
        states_nucleic="WXYZ",
        saving_frequency=20,
        keep_negative=True
    )
    assert condact.selected_resid == "1-5"
    assert condact.No_of_peaks_protein == 5
    assert condact.peak_boundaries_protein == [50, 150, 250, 310]
    assert condact.states_protein == "ABCDE"
    assert condact.No_of_peaks_nucleic == 4
    assert condact.peak_boundaries_nucleic == [100, 200, 300]
    assert condact.states_nucleic == "WXYZ"
    assert condact.saving_frequency == 20
    assert condact.keep_negative is True


# Test the state of the dihedral angles of nucleotides and amino acids
def test_define_state_protein_and_nucleic():
    condact = CONDACT(universe)

    # protein
    assert condact.define_state_protein(10) == "X"
    assert condact.define_state_protein(130) == "Y"
    assert condact.define_state_protein(260) == "Z"
    assert condact.define_state_protein(260) != "X"
    assert condact.define_state_protein(130) != "Z"
    with pytest.raises(ValueError):
        condact.define_state_protein(-40)

    # nucleic
    assert condact.define_state_nucleic(10) == "S"
    assert condact.define_state_nucleic(200) == "A"
    assert condact.define_state_nucleic(340) == "S"
    assert condact.define_state_nucleic(340) != "A"
    with pytest.raises(ValueError):
        condact.define_state_nucleic(-50)
    
    with pytest.raises(ValueError) as excinfo:
        condact.define_state_nucleic(-20)
    print(f"Error message: {str(excinfo.value)}")
    # Assert the error message
    assert str(excinfo.value) == "Angle must be positive because secondary angles used"


# Test the greater than functions to find time greater than a specific transition time
def test_find_gt():
    condact = CONDACT(universe)
    transition_times = {'LYS 33': [0, 7910, 7920, 9570, 11600, 12610, 15140],
            'PHE 34': [0,7420,7520,7590,8130,9330,9690,10010,10840,11410,11490,11800,13360,13420,14750]}
    assert condact.find_gt(transition_times['PHE 34'], 9570) == 9690
    with pytest.raises(AssertionError):
        assert condact.find_gt(transition_times['PHE 34'], 14750) == 12610
    with pytest.raises(AssertionError):
        assert condact.find_gt(transition_times['LYS 33'], 7420) == 0
    assert condact.find_gt(transition_times['LYS 33'], 7420) == 7910


# Test selected nucleic acid residues
def test_selected_resid_none_error_message():
    """Ensure residues() raises the exact ValueError when selected_resid is None."""
    condact = CONDACT(universe, selected_resid=None)
    with pytest.raises(ValueError) as excinfo:
        condact.residues()
    assert str(excinfo.value) == "please select residues of interest in format a-b or a b ..."

def test_peak_boundaries_protein_length_error_message():
    """Check exact wording for protein peak boundary mismatch."""
    condact = CONDACT(
        universe,
        selected_resid="1-3",
        No_of_peaks_protein=3,
        peak_boundaries_protein=[120],  # check length of boundary list
    )
    with pytest.raises(ValueError) as excinfo:
        condact.residues()
    assert "If there are 3 peaks, the list of boundaries separating the amino acid peaks in protein in should have 2 numbers." in str(excinfo.value)

def test_peak_boundaries_nucleic_length_error_message():
    """Check exact wording for nucleic peak boundary mismatch."""
    condact = CONDACT(
        universe,
        selected_resid="1-3",
        No_of_peaks_nucleic=3,
        peak_boundaries_nucleic=[140],  # # check length of boundary list
    )
    with pytest.raises(ValueError) as excinfo:
        condact.residues()
    assert "If there are 3 peaks, the list of boundaries separating the nucleic peaks in nucleotides in should have 2 numbers." in str(excinfo.value)


# Test selected amino acid residues
def test_residues_selection_excludes_ala_gly():
    condact = CONDACT(universe, selected_resid="1-5")
    ag, names, residues = condact.residues()
    assert isinstance(names, list)
    assert 'ALA' not in names
    assert 'GLY' not in names
    assert len(names) > 0
    assert len(ag) > 0

# Test the run function if csv file is created and contains postive dihedral angles
def test_run():
    # Create an instance of ConditionalActivity for testing
    condact = CONDACT(universe,
                              selected_resid="1-3",
                              No_of_peaks_protein=3,
                              peak_boundaries_protein=[120, 240] ,
                              saving_frequency=10,
                               states_protein= "XYZ")
    
    transition_times, wait_times = condact.run()
    # Perform assertions on the output
    assert isinstance(transition_times, dict), "transition_times should be a dictionary"
    assert isinstance(wait_times, dict), "wait_times should be a dictionary"

    # Length mismatch between transition_times and wait_times for residue
    for key in transition_times:
        assert len(transition_times.keys()) == len(wait_times.keys())

    # Check if the output file exists
    assert os.path.exists(OUTFILE), f"Output file '{OUTFILE}' does not exist"

    # Read the output file and verify the dihedral angle calculations
    with open(OUTFILE, 'r') as f:
        lines = f.readlines()
    assert lines, f"Output file '{OUTFILE}' is not empty"

    # Ignore the header which starts with 'frame, time(ps), residue ...'
    for line in lines[1:]:  # Skip the header line
        values = line.strip().split(',')[2:]  # Skip 'frame, time(ps)'
        for value in values:
            try:
                angle = float(value.strip())
                assert angle >= 0, f"Found negative angle value: {angle}"
            except ValueError:
                pytest.fail(f"Failed to convert value '{value}' to float")  


def test_mean_persistent_time():
    # Create an instance of Conditional Activity for testing  
    condact = CONDACT(universe,
                        selected_resid="1-3",
                        saving_frequency=10)
    
    persistence_times, list_persistence_times = condact.mean_persistent_time()
    assert persistence_times == {'GLN 1': 947.1316666666667, 'ASN 2': 0.0, 'ASP 3': 3805.0616666666665}
    assert list_persistence_times == [947.1316666666667, 0.0, 3805.0616666666665]

def test_mean_exchange_time():
    # Create an instance of ConditionalActivity for testing
    condact = CONDACT(universe,
                              selected_resid="1-3",
                              saving_frequency=10)
    # call mean_exchange_times function
    exchange_time = condact.mean_exchange_time()
    assert exchange_time == [[319.9166666666667, 0.0, 103.28333333333333], [0.0, 0.0, 0.0], [3056.0233333333335, 0.0, 776.2666666666667]]

def test_mean_conditional_activity():
    condact = CONDACT(universe,
                        selected_resid="1-3",
                        No_of_peaks_protein=3,
                        peak_boundaries_protein= [120, 240] ,
                        saving_frequency=10,
                        states_protein= "XYZ")
    conditional_activity, persistence_times, exchange_times = condact.mean_conditional_activity()
    assert conditional_activity == [[1.0853775738595968, math.nan, 2.2159620983323873], [math.nan, math.nan, math.nan], [0.21921769139965167, math.nan, 1.5895913731466649]]

    # Check if the output file exists
    final = "Conditional_Activity_list.npy"
    assert os.path.exists(final), f"Output file '{final}' does not exist"