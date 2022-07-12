# helix

Input files are:
 - Amino acid letter + h (or + b): Ah, Rh, Lb... are the energy files of each amino acid in human (h) or bacterial (b) membrane. Charged amino acids: K+b, K0b...
 - input is the initial set of sequences from which the algorithm should evolve. It is a random set of sequences.
 - init_pos20 are the coordinates for a 20 amino acid long alpha helix (used to set the initial coordinates)
 - all.gro contains the coordinates of the atoms of all amino acids --> used to calculate the center of mass of each amino acid
 
