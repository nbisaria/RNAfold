####### This script has all the parameters and helper functions that the rest of the library generation will use i think it will wor not sure


import RNA
import subprocess
from random import random
from walkerrandom import Walkerrandom
from collections import defaultdict




def GetMotifstoFilter():
    splicesites = ['AAUAA','GGUAAGU']
    seqs = splicesites
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
    MS2_seq = 'gcgcACATGAGGATCACCCATGTgc' #extended hp
    return MS2_seq

def GenerateRandomScaffAndLinker(n):
    per_AU_scaff = 0.65
    per_AU_linker = 0.85 
    len_scaff = 10
    len_linker = 3
    Scaffprob = dict( A=per_AU_scaff/2, U=per_AU_scaff/2, C=(1-per_AU_scaff)/2, G=(1-per_AU_scaff)/2 )
    srand = Walkerrandom( Scaffprob.values(), Scaffprob.keys() )
    Linkprob = dict( A=per_AU_linker/2, U=per_AU_linker/2, C=(1-per_AU_linker)/2, G=(1-per_AU_linker)/2 )
    lrand = Walkerrandom( Linkprob.values(), Linkprob.keys())
    numhits = 0
    contexts_final = defaultdict()
    seqsfilter = GetMotifstoFilter()
    def makeSeqwProb(l,wrand):
        s = ''
        for i in range(0,l):
            s = s + wrand.random()
        return s
    
    while numhits < n:
        l1 = makeSeqwProb(len_linker, lrand)
        l2 = makeSeqwProb(len_linker, lrand)
        l3 = makeSeqwProb(len_linker, lrand)
        l4 = makeSeqwProb(len_linker, lrand)
        s1 = makeSeqwProb(len_scaff, srand)
        s2 = makeSeqwProb(len_scaff, srand)
        MS2 = Get_MS2()
        seq_final =  l1+MS2 + l2 +s1 + s2 + l3 + MS2 + l4 
        dG, db = GetdG_dotbracket(seq_final)
        constrain = db[0:3] + db[11:15] + db[22:25] + db[-3:] + db[-26:-22] + db[-14:-10]
        constrain = db[0:3] + db[11:15] + db[-3:] + db[-14:-10]
        # print(dG,dotbracket)
        # print(db)
        if (dG > -14.8) and ('(' not in constrain) and (')' not in constrain) and not any(motif in seq_final  for motif in seqsfilter):
            print(seq_final)
            print(dG, db)
            numhits = numhits + 1
            contexts_final['c' + str(numhits)] = [l1,l2, MS2, s1, s2, MS2, l3, l4]
    return contexts_final

def main():
    n=10
    c_final = GenerateRandomScaffAndLinker(n)

if __name__ == "__main__":
    main()


