"""Amino-acid property scales and pKa values used by the descriptor module."""

STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")

# Kyte-Doolittle hydropathy, used for GRAVY and sliding-window hydrophobic patches.
KYTE_DOOLITTLE = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5,
    "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5,
    "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6,
    "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}

# Eisenberg consensus scale, used for the CDR hydrophobicity index.
EISENBERG = {
    "A": 0.62, "R": -2.53, "N": -0.78, "D": -0.90, "C": 0.29,
    "Q": -0.85, "E": -0.74, "G": 0.48, "H": -0.40, "I": 1.38,
    "L": 1.06, "K": -1.50, "M": 0.64, "F": 1.19, "P": 0.12,
    "S": -0.18, "T": -0.05, "W": 0.81, "Y": 0.26, "V": 1.08,
}

# Side-chain pKa values (EMBOSS set) for Henderson-Hasselbalch charge curves.
PKA_BASIC = {"K": 10.8, "R": 12.5, "H": 6.5}
PKA_ACIDIC = {"D": 3.9, "E": 4.1, "C": 8.5, "Y": 10.1}
PKA_N_TERMINUS = 8.6
PKA_C_TERMINUS = 3.6

# Formal charge sign per ionisable residue, for simple counting descriptors.
POSITIVE = "KR"
NEGATIVE = "DE"
