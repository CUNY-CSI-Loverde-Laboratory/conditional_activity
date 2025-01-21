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
from conditional_activity.tests.utils import make_Universe
import math
# import conditional_activity



# import MDAnalysis
# from MDAnalysis.analysis.base import AnalysisBase

# from numpy.testing import assert_allclose



# universe = Universe("/Users/augustineonyema/MolSSI/Protein_Conditional_Activity/Protein_Conditional_Activity/Protein_Conditional_Activity/tests/pytest_LYS.prmtop", 
#                     "/Users/augustineonyema/MolSSI/Protein_Conditional_Activity/Protein_Conditional_Activity/Protein_Conditional_Activity/tests/pytest_LYS.xtc")

universe = Universe("pytest_LYS.prmtop", 
                    "pytest_LYS.xtc")


OUTFILE = 'Dihedral_Angle.csv'

def test_conditional_activity_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "conditional_activity" in sys.modules


def test_mdanalysis_logo_length(mdanalysis_logo_text):
    """Example test using a fixture defined in conftest.py"""
    logo_lines = mdanalysis_logo_text.split("\n")
    assert len(logo_lines) == 46, "Logo file does not have 46 lines!"

@pytest.fixture
def sample_universe():
    # You can set up a mock universe object if needed or just pass None for testing purposes.
    return CONDACT(universe)


def test_initialization_with_defaults(sample_universe):
    condact = sample_universe
    # Check default values
    assert condact.No_of_peaks_protein == 3
    assert condact.peak_boundaries_protein == [120, 240]
    assert condact.states_protein == "XYZ"
    assert condact.No_of_peaks_DNA == 3
    assert condact.peak_boundaries_DNA == [140, 330]
    assert condact.states_DNA == "SAS"
    assert condact.saving_frequency is None
    assert condact.keep_negative is False


def test_initialization_with_custom_values(sample_universe):
    # universe = None  # Replace with actual universe object if necessary
    condact = CONDACT(
        sample_universe,
        selected_resid=10,
        No_of_peaks_protein=5,
        peak_boundaries_protein=[100, 200],
        states_protein="ABC",
        No_of_peaks_DNA=4,
        peak_boundaries_DNA=[150, 350],
        states_DNA="DEF",
        saving_frequency=10,
        keep_negative=True
    )
    
    # Check custom values
    assert condact.selected_resid == 10
    assert condact.No_of_peaks_protein == 5
    assert condact.peak_boundaries_protein == [100, 200]
    assert condact.states_protein == "ABC"
    assert condact.No_of_peaks_DNA == 4
    assert condact.peak_boundaries_DNA == [150, 350]
    assert condact.states_DNA == "DEF"
    assert condact.saving_frequency == 10
    assert condact.keep_negative is True


def test_define_state_protein(sample_universe):
    # Create an instance of CONDACT for testing
    states_protein = "XYZ"
    activity = CONDACT(sample_universe,
                       states_protein=states_protein)
    
    # Test define_state function
    angle = 45
    stateX = activity.define_state_protein(angle)
    # Assert the expected result
    assert stateX == "X"

    angle = 200
    stateY = activity.define_state_protein(angle)
    # Assert the expected result
    assert stateY == "Y"

    angle = 300
    stateZ = activity.define_state_protein(angle)
    # Assert the expected result
    assert stateZ == "Z"

    angle = 20
    invalid_stateX = activity.define_state_protein(angle)
    # Assert the expected result
    assert invalid_stateX != "Z"

    angle = 130
    invalid_stateY = activity.define_state_protein(angle)
    # Assert the expected result
    assert invalid_stateY != "X"

    angle = 260
    invalid_stateZ = activity.define_state_protein(angle)
    # Assert the expected result
    assert invalid_stateZ != "Y"

    angle = -20
    with pytest.raises(ValueError) as excinfo:
        activity.define_state_protein(angle)
    print(f"Error message: {str(excinfo.value)}")
    # Assert the error message
    assert str(excinfo.value) == "Angle must be positive because secondary angles used"


def test_find_gt(sample_universe):
    activity = CONDACT(sample_universe)
    transition_times = {'LYS 33': [0, 7910, 7920, 9570, 11600, 12610, 15140],
            'PHE 34': [0,7420,7520,7590,8130,9330,9690,10010,10840,11410,11490,11800,13360,13420,14750]}
    assert activity.find_gt(transition_times['PHE 34'], 9570) == 9690
    with pytest.raises(AssertionError):
        assert activity.find_gt(transition_times['PHE 34'], 14750) == 12610
    with pytest.raises(AssertionError):
        assert activity.find_gt(transition_times['LYS 33'], 7420) == 0
    assert activity.find_gt(transition_times['LYS 33'], 7420) == 7910


def test_selected_resid_none(sample_universe):
    selected_resid = None
    activity = CONDACT(sample_universe,
                       selected_resid=selected_resid)
    
    with pytest.raises(ValueError, match=r"please select residues of interest in format a-b or a b ..."):
        activity.residues()

def test_peak_boundaries_length_DNA(sample_universe):
    selected_resid = "1-3"
    No_of_peaks_DNA = 3
    peak_boundaries_DNA = [120]  # Only 1 boundary, which is incorrect for 3 peaks
    
    # Test case for peak_boundaries length not matching No_of_peaks - 1
    activity = CONDACT(sample_universe,
                       selected_resid=selected_resid,
                       No_of_peaks_protein=No_of_peaks_DNA,
                       peak_boundaries_protein=peak_boundaries_DNA)
    
    # Update the expected error message based on No_of_peaks_protein = 3
    with pytest.raises(ValueError, match=r"If there are 3 peaks, the list of boundaries separating the amino acid peaks in protein in should have 2 numbers."):
        activity.residues()


def test_peak_boundaries_length(sample_universe):
    selected_resid = "1-3"
    No_of_peaks_protein = 3
    peak_boundaries_protein = [120]  # Only 1 boundary, which is incorrect for 3 peaks
    
    # Test case for peak_boundaries length not matching No_of_peaks - 1
    activity = CONDACT(sample_universe,
                       selected_resid=selected_resid,
                       No_of_peaks_protein=No_of_peaks_protein,
                       peak_boundaries_protein=peak_boundaries_protein)
    
    # Update the expected error message based on No_of_peaks_protein = 3
    with pytest.raises(ValueError, match=r"If there are 3 peaks, the list of boundaries separating the amino acid peaks in protein in should have 2 numbers."):
        activity.residues()


def test_residues(sample_universe):
    # Create an instance of ConditionalActivity for testing
    selected_resid = "1-3"
    activity = CONDACT(sample_universe,
                       selected_resid=selected_resid)
    
    interested, residue_names, residue_type = activity.residues()
    assert 'ALA' not in residue_names
    assert 'GLY' not in residue_names
    assert len(residue_names) > 0
    assert len(interested) > 0
    assert len(residue_type) > 0
    assert isinstance(residue_names, list)


def test_run(sample_universe):
    # Create an instance of ConditionalActivity for testing
    selected_resid = "1-3"
    No_of_peaks_protein = 3
    peak_boundaries_protein = [120, 240]  
    states_protein = "XYZ"
    saving_frequency = 10

    activity = CONDACT(sample_universe,
                              selected_resid=selected_resid,
                              No_of_peaks_protein=No_of_peaks_protein,
                              peak_boundaries_protein=peak_boundaries_protein,
                              saving_frequency=saving_frequency,
                               states_protein= states_protein)
    
    transition_times, wait_times = activity.run()
    # Perform assertions on the output
    assert isinstance(transition_times, dict), "transition_times should be a dictionary"
    assert isinstance(wait_times, dict), "wait_times should be a dictionary"
    
    # Length mismatch between transition_times and wait_times for residue
    for key in transition_times:
        assert len(transition_times.keys()) == len(wait_times.keys())

    # Run the method that generates the output file
    activity.run()

    # Check if the output file exists
    assert os.path.exists(OUTFILE), f"Output file '{OUTFILE}' does not exist"

    # Read the output file and verify the dihedral angle calculations
    with open(OUTFILE, 'r') as f:
        lines = f.readlines()
    assert lines, f"Output file '{OUTFILE}' is empty"

    # Ignore the header which starts with 'frame, time(ps), residue ...'
    for line in lines[1:]:  # Skip the header line
        values = line.strip().split(',')[2:]  # Skip 'frame, time(ps)'
        for value in values:
            try:
                angle = float(value.strip())
                assert angle >= 0, f"Found negative angle value: {angle}"
            except ValueError:
                pytest.fail(f"Failed to convert value '{value}' to float")  


def test_mean_persistent_time(sample_universe):
    # Create an instance of Conditional Activity for testing  
    selected_resid = "1-3"
    saving_frequency = 10
    activity = CONDACT(sample_universe,
                        selected_resid=selected_resid,
                        saving_frequency=saving_frequency)
    
    persistence_times, list_persistence_times = activity.mean_persistent_time()
    assert persistence_times == {'GLN 1': 928.405, 'ASN 2': 0, 'ASP 3': 848.0216666666666}
    assert list_persistence_times == [928.405, 0.0, 848.0216666666666]



def test_mean_exchange_time(sample_universe):
    # Create an instance of ConditionalActivity for testing
    selected_resid = "1-3"
    saving_frequency = 10
    activity = CONDACT(sample_universe,
                              selected_resid=selected_resid,
                              saving_frequency=saving_frequency)
    # call mean_exchange_times function
    exchange_time = activity.mean_exchange_time()
    assert exchange_time == [[307.9033333333333, 0.0, 2622.8366666666666], [0.0, 0.0, 0.0], [85.52333333333333, 0.0, 554.2666666666667]]


def test_mean_conditional_activity(sample_universe):
    selected_resid = "1-3"
    No_of_peaks_protein = 3
    peak_boundaries_protein = [120, 240]  
    states_protein = "XYZ"
    saving_frequency = 10

    activity = CONDACT(sample_universe,
                        selected_resid=selected_resid,
                        No_of_peaks_protein=No_of_peaks_protein,
                        peak_boundaries_protein=peak_boundaries_protein,
                        saving_frequency=saving_frequency,
                        states_protein= states_protein)
    conditional_activity, persistence_times, exchange_times = activity.mean_conditional_activity()
    assert conditional_activity == [[1.103682179097811, math.nan, -1.0385436482668298], [math.nan, math.nan, math.nan], [2.2941169425854744, math.nan, 0.425260267069592]]
