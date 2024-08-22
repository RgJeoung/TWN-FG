# TWN-FG (TWN-based Fragment Growing)
Source code for TWN-FG method.
![graphical_abs](https://github.com/user-attachments/assets/de6784d7-d8f6-4b6b-9c51-62c38db82a9d)
TWN-FG is conducted through the following procedure.
## Installation
This code was tested in windows with Python 3.9.

A yaml file containing all requirements is provided.

You can set up a virtual environment in conda as follows:

    conda env create -f TWN-Region-Analysis.yaml
    conda activate TWN-Region-Analysis

## Preset
To run this code, you need the following three preset files.
1. Molecular dynamics (MD) simulation trajectories for target protein
2. Topological water networks (TWN)
3. Boundary file

###### MD trajectories
We used GROMACS version 5.1.4 and the CHARMM27 forcefield to get MD trajecotries. For this analysis, you can use input files available in DATA directory or you can prepare your input files in similar way.

###### TWNs
TWNs can be obtained by running the TWN_Extractor_v1.exe file.
The input file must contain a_input in the upper PDB directory (for example, 1NVR in the trajectory directory) and must be entered in the format of "trajectory number.pdb" (for example, 0.pdb).

For example:

    TWN_Extractor_v1.exe -path ./DATA/trajectory/1NVR

The output file is created as "a_output" in the same pdb directory by default. In detail, you can set the output directory name and whether to include duplicates of rings of larger size ("True" by default which means include duplicates).

For example:

    TWN_Extractor_v1.exe -path ./DATA/trajectory/1NVR -opname "output directory name" -dup False
    
After the extraction process is completed, you must change the rfour directory created in the output directory to the "protein-name_pdb-code_ring-type" format to create a TWN directory like "CHK1_1NVR_R4" in the DATA/TWN directory.
###### Boundary file
Boundary file can be obtained by running the Boundary_file_maker.py file.
The input file must be same as the input of TWN_Extractor_v1.exe.

For example:

    python boundary_file_maker.py -path ./DATA/trajectory/1NVR

You can choose options after typing this at the prompt (for example, centering methods and sizes). After running you can get a ".bd" format file in the directory name of "boundary".

## Start TWN region analysis
After preparing all required presets, you are ready to run the main code, "TWN-Region-Analysis.py". Please make sure to prepare all directories like the example "DATA".

Script for running code is following:

    python TWN-Region-Analysis.py -d ./DATA

After running the code, you can get three directories.
1. TWN Pattern
2. TWN Region
3. logs

TWN patterns indicate independent positions and shapes of water networks within protein binding pocket. TWN regions are identified based on the TWN patterns.

For better understanding of output files, you can read our paper.

## Contact (Questions/Bugs/Requests)
Questions : Please ask our professor <nskang@cnu.ac.kr>

Bugs/Requests : Please submit a GitHub issue or contact me <ray971125@naver.com>
