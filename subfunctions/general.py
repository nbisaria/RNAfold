####### This script has all the parameters and helper functions that the rest of the library generation will use i think it will wor not sure













class Parameters():
    def __init__(self, concentrations=None):
        # 
        self.RT = 0.616  # kcal/mol at 37 degrees celsius
        self.min_deltaG = -12
        self.max_deltaG = -3
        self.tempK = 310

def GetScaffoldSeqs():
    f_seqs = {
    "scaff_1" : "UGGAAUGUAAAGAAGUAUGUAU",
    "scaff_2" : "UGGAAUGUAAAGAAGUAUGUAU",
    "scaff_3" : "UUAAUGCUAAUCGUGAUAGGGGU",
    "scaff_4" : "UAAGGCACGCGGUGAAUGCCAA",
    "scaff_5" : "UGAGGUAGUAGGUUGUAUAGUU",
	}
	return f_seqs


def GetLinkerSeqs():
    l_seqs = {
    "link_1" : "UGGAAUGUAAAGAAGUAUGUAU",
    "link_2" : "UGGAAUGUAAAGAAGUAUGUAU",
    "link_3" : "UUAAUGCUAAUCGUGAUAGGGGU",
    "link_4" : "UAAGGCACGCGGUGAAUGCCAA",
    "link_5" : "UGAGGUAGUAGGUUGUAUAGUU",
	}
