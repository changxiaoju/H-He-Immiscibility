
# Atomic Miscibility Analysis Based on Reweighted Conditional Probability (AMA-RCP)

## Objective
This project provides two scripts to calculate the Pcond values from a molecular trajectory in XYZ format and perform re-weighting of the results based on specific criteria.

## Dependencies
To run this project, you need the following Python packages:
- Python 3.6+
- numpy
- pandas
- matplotlib
- scipy

You can install all required dependencies using the following command:

```
pip install numpy pandas matplotlibi scipy
```

In addition, the code uses standard Python libraries:
- os
- glob
- re
- math (log, exp)
- time

These libraries are included in the Python standard library and do not require separate installation.

## Installation and Usage
0. **Clone the MolecularModTools repository**:
   -git clone https://github.com/changxiaoju/MolecularModTools

1. **Prepare your working directory**:
   - Download `00.get_Pcond_from_xyz_traj.py` and  `01.re-weight_Pcond_get_x1_x2.py` to your working directory.

1. **Prepare your XYZ trajectory file**:
   - The trajectory file should be in XYZ format (e.g., `example.xyz`), where each frame contains the atomic coordinates.
   - Put the trajectory file in work directory.

2. **Modify the scripts**:
   - Open the scripts `00.get_Pcond_from_xyz_traj.py` and `01.re-weight_Pcond_get_x1_x2.py`.
   - Update the following variables as needed:
     - **ADD_PATH_TO**: The path to MolecularModTools.
     - **PATH_TO_XYZ_FILE**: Ensure the paths to your XYZ file and other required data are correct.
     - **NAME_OF_XYZ_FILE**: Specify the name of the results should be saved.
     - **Parameters**: Modify number of neighbors `nNB`, abundance values `xA`, ratio of trajectory to use `ratio` according to your data.

3. **Run the scripts**:
   - First, calculate the Pcond values:
     ```
     python 00.get_Pcond_from_xyz_traj.py
     ```
   - Next, perform the re-weighting of the calculated Pcond values:
     ```
     python 01.re-weight_Pcond_get_x1_x2.py
     ```

## Output
- **Pcond Results**:
  - The first script generates 2 files containing the calculated Pcond values based on the input trajectory.
  
- **Re-weighted Results**:
  - The second script prints the abundance diverge results. You can modify these two scripts for batch processing.

