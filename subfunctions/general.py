####### This script has all the parameters and helper functions that the rest of the library generation will use i think it will wor not sure


import RNA
import subprocess
from random import random
from walkerrandom import Walkerrandom









class Parameters():
    def __init__(self, concentrations=None):
        # 
        self.RT = 0.616  # kcal/mol at 37 degrees celsius
        self.min_deltaG = -12
        self.max_deltaG = -3
        self.tempK = 310

def GetScaffoldSeqs():

    #scaffold sequences need to x nt in length with certain amount/varying amount of AU content and stability
    #two sequences are the 5' and 3' sequences that will be threaded together to make the final sequence
    f_seqs = {
    "scaff_1" : ["UGGAAUGUAA,AGAAGUAUGUAU"],
    "scaff_2" : "UGGAAUGUAA-AGAAGUAUGUAU",
    "scaff_3" : "UUAAUGCUAA-UCGUGAUAGGGGU",
    "scaff_4" : "UAAGGCACGC-GGUGAAUGCCAA",
    "scaff_5" : "UGAGGUAGUA-GGUUGUAUAGUU",
	}
	return f_seqs


def GetLinkerSeqs():
    l_seqs = {
    "link_1" : ['AAU', 'UAU', 'AUU', 'UUA'],
    "link_2" : "UGGAAUGUAAAGAAGUAUGUAU",
    "link_3" : "UUAAUGCUAAUCGUGAUAGGGGU",
    "link_4" : "UAAGGCACGCGGUGAAUGCCAA",
    "link_5" : "UGAGGUAGUAGGUUGUAUAGUU",
	}

def GetdG_dotbracket(seq):
    output_ = subprocess.check_output('echo \'' + seq +'\' | RNAfold -d 0', shell=True)
    dG = float(output_[1].split(' (')[1].split(')')[0])
    dotbracket = output_[1].split(' (')[0]
    return dG, dotbracket

def Get_MS2():
    MS2_seq = ''
    return MS2_seq

def GenerateRandomScaffAndLinker():
    per_AU_scaff = 0.6
    per_AU_linker = 0.85 
    len_scaff = 10
    len_linker = 3
    numhits = 0
    seed(1)
    while numhits <10 :

    

    # simulate sequence





abcd = dict( A=3, U=3, C=2, G=2 )
  # keys can be any immutables: 2d points, colors, atoms ...
wrand = Walkerrandom( abcd.values(), abcd.keys() )
wrand.random()


