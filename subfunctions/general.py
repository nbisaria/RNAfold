####### This script has all the parameters and helper functions that the rest of the library generation will use i think it will wor not sure


import RNA
import subprocess
import pandas as pd
import imp
from random import random
imp.load_source("walkerrandom",("/lab/bartel1_ata/nbisaria/RNAfold/subfunctions/walkerrandom.py"))
from walkerrandom import Walkerrandom
from collections import defaultdict




def GetMotifstoFilter():
    splicesites = ['GGUAAGU']
    polyA = ['AAUAA']
    seqs = splicesites + polyA
    return seqs

class Parameters():
    def __init__(self, concentrations=None):
        # 
        self.RT = 0.616  # kcal/mol at 37 degrees celsius
        self.min_deltaG = -12
        self.max_deltaG = -3
        self.tempK = 310


def GetdG_dotbracket(seq):
    output_ = subprocess.check_output('echo \'' + seq +'\' | RNAfold -d 0', shell=True).split('\n')
    dG = float(output_[1].split(' (')[1].split(')')[0])
    dotbracket = output_[1].split(' (')[0]
    return dG, dotbracket

def Get_MS2():
    MS2_seq = 'gcgcACATGAGGATCACCCATGTgc' #extended hp, from zalatan et al.
    dG, dotbracket = GetdG_dotbracket(MS2_seq)
    print(dG, dotbracket, MS2_seq)
    return MS2_seq
def Get_PP7():
    PP7 = 'aacaTAAGGAGTTTATATGGAAACCCTTAtg' #from zalatan et al.
    dG, dotbracket = GetdG_dotbracket(PP7)
    print(dG, dotbracket, PP7)
    return PP7

def get_rc(seq,rna=False):
    """Gives the reverse complementary seqence of DNA.

        Args:
            seq: The sequence to be reverse complemented, in either DNA or RNA
            form.

        Returns:
            The DNA reverse comlementary sequence.
    """
    # Builds dictionary mapping each sequence to its complement.
    nucleotide_complement_map = {
    "C" : "G",
    "G" : "C",
    "T" : "A",
    "U" : "A"  
    }
    # Assigns A to U or T depending on 'rna' flag.
    if rna:
        nucleotide_complement_map["A"] = "U"
    else:
        nucleotide_complement_map["A"] = "T"
    # Constructs a list of the characters substituted and reversed in order.
    rev_complement_list = [nucleotide_complement_map[i] for i in seq][::-1]
    rev_complement = "".join(rev_complement_list)
    return(rev_complement)

def Get_PairedEndPrimersSeqs():
    # read 5' to 3' 
    TrueSeqR2 = 'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT' # Tm 67C
    TrueSeqR1 = 'GATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT' # Tm 67C

    return TrueSeqR2, TrueSeqR1
