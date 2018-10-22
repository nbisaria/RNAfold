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
    output_ = subprocess.check_output('echo \'' + seq +'\' | RNAfold -d 0', shell=True).split('\n')
    dG = float(output_[1].split(' (')[1].split(')')[0])
    dotbracket = output_[1].split(' (')[0]
    return dG, dotbracket

def Get_MS2():
    MS2_seq = ''
    return MS2_seq

def GenerateRandomScaffAndLinker():
    per_AU_scaff = 0.65
    per_AU_linker = 0.85 
    len_scaff = 10
    len_linker = 3
    Scaffprob = dict( A=per_AU_scaff/2, U=per_AU_scaff/2, C=(1-per_AU_scaff)/2, G=(1-per_AU_scaff)/2 )
    srand = Walkerrandom( Scaffprob.values(), Scaffprob.keys() )
    Linkprob = dict( A=per_AU_linker/2, U=per_AU_linker/2, C=(1-per_AU_linker)/2, G=(1-per_AU_linker)/2 )
    lrand = Walkerrandom( Linkprob.values(), Linkprob.keys())
    numhits = 10
    def makeSeqwProb(l,wrand):
        s = ''
        for i in range(0,l):
            s = s + wrand.random()
        return s
    seqs_final = defaultdict()
    while numhits < 10:
        l1 = makeSeqwProb(len_linker, lrand)
        l2 = makeSeqwProb(len_linker, lrand)
        l3 = makeSeqwProb(len_linker, lrand)
        l4 = makeSeqwProb(len_linker, lrand)
        s1 = makeSeqwProb(len_scaff, srand)
        s2 = makeSeqwProb(len_scaff, srand)
        seq_final =  l1+s1 + l2 + l3 + s2 + l4 
        dG, dotbracket = GetdG_dotbracket(seq_final)
        if (dG > -0.1) and ('(' not in dotbracket):
            print(seq_final)
            print(dG, dotbracket)
            numhits = numhits + 1

