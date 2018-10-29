import RNA
import subprocess
import pandas as pd
import imp
from random import random
imp.load_source("walkerrandom",("/lab/bartel1_ata/nbisaria/RNAfold/subfunctions/walkerrandom.py"))
from walkerrandom import Walkerrandom
from collections import defaultdict
imp.load_source("general",("/lab/bartel1_ata/nbisaria/RNAfold/subfunctions/general.py"))
from general import *

def GenerateRandomScaffAndLinker(n):
    per_AU_scaff = 0.65
    per_AU_linker = 0.85 
    len_scaff = 10
    len_linker = 2
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
    MS2 = Get_MS2()
    PP7 = Get_PP7()
    R2, R1 = Get_PairedEndPrimersSeqs()
    db_final = '....(((((((.((((......)))))))))))..........................(((((((.((....)))))))))..'
    while numhits < n:
        l1 = makeSeqwProb(len_linker, lrand)
        l2 = makeSeqwProb(len_linker, lrand)
        l3 = makeSeqwProb(len_linker, lrand)
        l4 = makeSeqwProb(len_linker, lrand)
        s1 = makeSeqwProb(len_scaff, srand)
        s2 = makeSeqwProb(len_scaff, srand)
        # right now R2 and revcomp(R1) makes a giant hairpin:
        #'...(((....(((..(((((((((((((((((....)))))))))))).)))))..)))...))).......'
        # this is using the R2 + get_rc(R1)
        # perhaps need to add a constant new sequence + BstX1 sites and use some dark cycles
        seq_final =  R2 + l1+PP7 + l2 +s1 + s2 + l3 + MS2 + l4 + get_rc(R1)
        dG, db = GetdG_dotbracket(seq_final)
        # constrain = db[0:3] + db[11:15] + db[22:25] + db[-3:] + db[-26:-22] + db[-14:-10]
        # constrain = db[0:3] + db[11:15] + db[-3:] + db[-14:-10]
        print(dG,db)
        # print(db)
        if (dG > -21) and (db[len(R2):] == db_final) and not any(motif in seq_final  for motif in seqsfilter):
            print(seq_final)
            print(dG, db)
            numhits = numhits + 1
            contexts_final['c' + str(numhits)] = [l1,PP7, l2, s1, s2, l3, MS2, l4]
    return contexts_final

def main():
    n=10
    c_final = GenerateRandomScaffAndLinker(n)
    print(c_final)

if __name__ == "__main__":
    main()

