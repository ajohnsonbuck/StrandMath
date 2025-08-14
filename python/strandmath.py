# -*- coding: utf-8 -*-
"""
Alex Johnson-Buck, 2025
The University of Michigan

Repository and documentation found at:
    https://github.com/ajohnsonbuck/StrandMath
"""

from typing import Union, Iterable, Optional
import copy
import re
import random
import numpy as np
import pandas as pd

class Strand:
    Modlist = ["+","b","r"]
    Nucleotides = ["A","C","G","T","U","a","c","g","t","u"]
    
    def __init__(self, sequence, name = ""):
        if isinstance(sequence,Strand): # If already a sequence, return a copy of the original Strand object
            self.name = copy.copy(sequence.name)
            self.sequence = copy.copy(sequence.sequence)
        elif isinstance(sequence,str):                    # if sequence is provided as string, interpret as single sequence
            self.sequence = [self.fromString(sequence)]  
            if isinstance(name,str):
                self.name = [name]
            else:
                raise TypeError("For single sequence, name must be a single string")
        elif isinstance(sequence,list):                 # if sequence is provided as list, interpret as multiple sequences
            if all(isinstance(seq, str) for seq in sequence): 
                self.sequence = []
                for seq in sequence:
                    self.sequence.append(self.fromString(seq))
            self.name = []
            if isinstance(name,list):
                if all(isinstance(item,str) for item in name):
                    for item in name:
                        self.name.append(item)
            elif name == "":
                for ind in range(len(sequence)):
                    self.name.append("")
            else:
                raise TypeError("For multiple input sequences, name must be provided as a list of strings")
        
    @property
    def name(self):
        return self._name
    
    @name.setter
    def name(self,value):
        if isinstance(value, str):
            self._name = [value] 
        elif isinstance(value, list): 
            self._name = value 
        else:
            raise TypeError("name must be a string or a list of strings")      

    def __getitem__(self,ind): # Indexing returns a new Strand object containing the corresponding indexed strands and names
        name = self.name[ind]
        seq = self.sequence[ind]
        if isinstance(name,str):
            name = [name]
            seq = [seq]
        newStrand = copy.copy(self)
        newStrand.name = name
        newStrand.sequence = seq
        return newStrand
    
    def __setitem__(self,ind,value):
        if isinstance(value,Strand):
            if value.numel() == 1:
                self.name[ind] = value.name[0]
                self.sequence[ind] = value.sequence[0]
            else:
                raise ValueError("Strand item assignment operand must be a Strand object with one sequence")
        else:
            raise TypeError("Strand only supports item assignment for Strand operands")

    def __repr__(self):
        if self.numel()==1:
            str1 = "".join(self.sequence[0])
            return f"Strand({str1!r})"
        else:
            return f"Strand({self.numel()} Sequences)"

    def __add__(self, other: Union['Strand',str]): # Addition of Strands means concatenating their respective sequences
        if isinstance(other, str): 
            other = Strand(other)
        if self.numel()==1:
            newStrand = copy.copy(other)
            for ind in range(other.numel()):
                newStrand.sequence[ind] = self.sequence[0] + other.sequence[ind]
                if len(self.name[0])>0:
                    newStrand.name[ind] = f"{self.name[0]} + {other.name[ind]}"
            return newStrand
        elif other.numel()==1:
            newStrand = copy.copy(self)
            for ind in range(self.numel()):
                newStrand.sequence[ind] = self.sequence[ind] + other.sequence[0]
                if len(other.name[0])>0:
                    newStrand.name[ind] = f"{self.name[ind]} + {other.name[0]}"
            return newStrand
        return NotImplemented
    
    def __radd__(self, other: str): # other + self = Strand(other) + self
        other = Strand(other)
        return other + self
    
    def __neg__(self): # Negative = reverse sequence 5'-to-3'
        return copy.copy(self).reverse()
    
    def __sub__(self,other): # self - other = self + other.reverse()
        if isinstance(other,str):
            other = Strand(other)
        return self + (-other)
    
    def __rsub__(self,other): # other - self = other + self.reverse()
        return other + (-self)
    
    def __invert__(self): # Invert operator ~ yields reverse complement
        return self.reverseComplement()
    
    def __eq__(self,other): # Equals operator: return True if two Strands have same sequences in the same order, and False otherwise
        if self.sequence == other.sequence:
            return True
        return False
    
    def __mul__(self,other):
        if isinstance(other,int): # Multiplying a Strand array by a constant b concatenates the sequence b times
            seq = self.string() * other
            return Strand(seq)
        elif isinstance(other,Strand): # Multiplying two Strand arrays of size m and n results in a Duplex array of size m*n
            return NotImplemented
        else:
            raise TypeError("* operator can only multiply a Strand by an int")
            
    def __rmul__(self,other):
        if isinstance(other,str):
            other = Strand(other)
        return self * other
    
    def isSymmetric(self): # Return True for each sequence that it is its own reverse complement (regardless of type of sugar), and False otherwise
        out = []
        for ind in range(self.numel()):
            if self[ind].toDNA() == ~(self[ind].toDNA()):
                out.append(True)
            else:
                out.append(False)
        if len(out)==1:
            out = out[0]
        return out
    
    def removeDuplicates(self): # Remove any duplicate sequences; keep only the first instance of any sequence + its name
        out = copy.copy(self)
        seqnew = []
        namenew = []
        seen = {}
        for seqrow, namerow in zip(self.sequence, self.name):
            key = tuple(seqrow)
            if key not in seen:
                seen[key]=True
                seqnew.append(seqrow)
                namenew.append(namerow)
        out.sequence = seqnew
        out.name = namenew
        return out 
    
    def scramble(self): # Scramble sequence(s) in Strand object
        out = copy.copy(self)
        for ind in range(out.numel()):
            random.shuffle(out.sequence[ind])
            out.name[ind] = out.name[ind]+'_scrambled'
        return out
    
    def fromString(self,str1: str): # Extract sequence from string input
        str1 = Strand.cleanString(str1) # Remove empty spaces and termini markers
        sequence = []
        c = 0
        m = 0
        while m < len(str1):
            n = m
            while not(str1[n] in Strand.Nucleotides):
                n += 1 
            sequence.append(str1[m:n] + str1[n].upper())
            m = n+1 
            c += 1 
        return sequence
        
    def string(self): # String representation of sequence
        strlist = []
        for item in self.sequence:
            strlist.append("".join(item))
        if len(strlist)==1:
            strlist = strlist[0]
        return strlist
        
    def bareString(self): # String representation of sequence(s), stripped of modifications
        strlist = []
        string = self.string()
        if isinstance(string,str):
            string = [string]
        for item in string:
            strlist.append(Strand.removeMods(item))
        if len(strlist) == 1:
            strlist = strlist[0]
        return strlist
    
    def bareSequence(self): # List representation of sequence(s), stripped of modifications
        bareSeq = copy.copy(self.sequence)
        for seq in range(len(bareSeq)):
            for nt in range(len(bareSeq[seq])):
                bareSeq[seq][nt] = Strand.removeMods(bareSeq[seq][nt])
        return bareSeq
            
    def removeMods(str1: str): # Remove modifications from string representation of sequence
        pattern = '|'.join(map(re.escape,Strand.Modlist))
        return re.sub(pattern,'', str1)
    
    def cleanString(str1: str): # Clean up sequence to remove termini markers and empty spaces
        pattern = '|'.join(map(re.escape,[" ", "5'-","-3'","5-","-3","5'","3'"])) 
        return re.sub(pattern,'',str1)
    
    def numel(self): # Number of sequences in Strand object
        return len(self.sequence)
    
    def toDNA(self): # Convert to DNA sequence
        out = copy.copy(self)
        out.sequence = out.bareSequence()
        out.sequence = [[s.replace('U','T') for s in row] for row in out.sequence]
        return out
        
    def toRNA(self): # Convert to RNA sequence
        out = copy.copy(self)
        out.sequence = out.bareSequence()
        out.sequence = [[s.replace('T','U') for s in row] for row in out.sequence]
        out.sequence = [['r' + s for s in row] for row in out.sequence]
        return out
    
    def toLNA(self): # Convert to LNA sequence
        out = copy.copy(self)
        out.sequence = out.bareSequence()
        out.sequence = [['+' + s for s in row] for row in out.sequence]
        return out
    
    def toBNA(self): # Convert to BNA sequence
        out = copy.copy(self)
        out.sequence = out.bareSequence()
        out.sequence = [['b' + s for s in row] for row in out.sequence]
        return out
    
    def len(self): # Number of nucleotides in each sequence
        L = []
        for item in self.sequence:
            L.append(len(item))
        if len(L)==1:
            L = L[0]
        return L
    
    def reverse(self): # reverse sequence(s)
        out = copy.copy(self)
        out.sequence = [row[::-1] for row in out.sequence]
        out.name = [name + '_reverse' for name in out.name]
        return out
    
    def reverseComplement(self): # Create reverse complement of all sequences in Strand
        out = copy.copy(self)
        out = out.reverse()
        for i, seq in enumerate(out.sequence):
            for j, nt in enumerate(seq):
                if nt[-1]=="C":
                    out.sequence[i][j] = out.sequence[i][j].replace("C","G")
                elif nt[-1]=="G":
                    out.sequence[i][j] = out.sequence[i][j].replace("G","C")
                elif nt[-1]=="T":
                    out.sequence[i][j] = out.sequence[i][j].replace("T","A")
                elif nt[-1]=="U":
                    out.sequence[i][j] = out.sequence[i][j].replace("U","A")
                elif nt[-2:]=="rA":
                    out.sequence[i][j] = out.sequence[i][j].replace("rA","rU")
                elif nt[-1]=="A":
                    out.sequence[i][j] = out.sequence[i][j].replace("A","T")
        out.name = [name + '_complement' for name in out.name]
        return out                    

    def crop(self,ind: list=(0,'end')): # Crop sequence(s) to the specified nucleotide range
        if any(not(isinstance(item, int)) and not(item=='end') for item in ind):
            raise TypeError("Indices for Strand.crop must be a list of integers or 'end'")
        if ind[1]=='end':
            ind[1] = np.inf
        if not(ind[1]>ind[0]):
            raise ValueError("Second index of Strand.crop must be larger than the first index")
        out = copy.copy(self)
        for i, row in enumerate(out.sequence):
            out.sequence[i] = out.sequence[i][max([0,ind[0]]):min([len(row),ind[1]])]
        return out
    
    def gcContent(self): # Fraction G and C bases for each sequence in Strand object
        fGC = [] 
        bareSeq = self.bareString()
        if isinstance(bareSeq,str):
            bareSeq = [bareSeq]
        for item in bareSeq:
            fGC.append((item.count('G') + item.count('C'))/len(item))
        if len(fGC)==1: # return float scalar if only one sequence in Strand; list otherwise
            fGC = fGC[0]
        return fGC
    
    def gc(self):
        return self.gcContent()
    
    def print(self): # Display names and sequences of all items in Strand object
        strlist = self.string()
        if isinstance(strlist,str):
            strlist = [strlist]
        for (stritem,nameitem) in zip(strlist, self.name):
            if nameitem != "":
                print(nameitem)
            print("5'-" + stritem + "-3'")
    
    def random(length: int, seqtype: str='DNA', gcContent: float=0.5,name: str=''): # Generate random sequence
        if name == '':
            name = f"Random_{seqtype}_{length}nt_{gcContent}GC"
        nGC = int(np.round(gcContent*length))
        str1 = "".join(random.choice(["C","G"]) for _ in range(nGC))
        str2 = "".join(random.choice(["A","T"]) for _ in range(length-nGC))
        seq = list(str1+str2)
        random.shuffle(seq)
        out = Strand("".join(seq),name=name)
        if seqtype=='RNA':
            out = out.toRNA()
        elif seqtype=='LNA':
            out = out.toLNA()
        elif seqtype=='BNA':
            out = out.toBNA()
        return out
    
    def polyN(nt: str, n: int):
        # Create homopolymer sequence of nucleotide nt with length n
        if nt[0]=='d': # Remove 'd' prefix for DNA, if present
            nt = nt.replace("d","")
        name = 'dT10' # f"({nt}){n}"
        out = n*Strand(nt)
        out.name = name
        return out
    
    def polyT(n: int):
        out = Strand.polyN('T',n)
        out.name = '(dT)' + str(n)
        return out 

class Multistrand:
    def __init__(self,strand0="",strand1=""):
        self.Strands = Strand(["",""])
        strand0 = Strand(strand0)
        strand1 = Strand(strand1)
        if len(strand1.string()) == 0:
            if len(strand0.string()) > 0:
                strand1 = ~strand0
        self.Strands[0] = strand0
        self.Strands[1] = strand1
        
        # if isinstance(strands,list):
        #     if len(strands) > 1 and all(isinstance(seq,str) for seq in strands):
        #         self.Strands = Strand(strands)
        #     else:
        #         raise ValueError("Multistrand must be initialized with two or more nucleic acid sequences, formatted as Strand or a list of str")
    
class Duplex:
    PARAMETERS = pd.read_csv("../NN_Parameters.csv")
    Strands = Strand(["",""]) # Strand array containing two interacting nucleic acids
    Schema = [[],[]] # 2xN List showing register of two sequences in interaction
    PairingState = [] # 1xN List showing pairing state ('p'=paired, 'w'=wobble,''=mismatch,'d'= dangling/overhang)
    NearestNeighbors = [] # List of codes for all nearest-neighbor interactions within the Duplex
    Nbp = np.array([]) # Number of base pairs in interaction (initialize as empty unless provided)
    Length = np.array([]) # ength of interaction, including mismatches and overhangs
    fGC = -np.inf # GC content of interaction
    dS0 = -np.inf # Entropy of hybridization at standard conditions
    dH0 = np.inf # Enthalpy of hybridization at standard conditions
    dG0 = np.inf # Free energy of hybridization at standard conditions
    

# # Example usage:
# A = Strand(["rArGrCrU","GACCTA"],name=["Sequence1","Sequence2"]) # Create a Strand object with two sequences
# A.print() # Show all sequences

# B = Strand('TTTTT',name="T5") # Create a strand object with a single sequence (5T linker)

# (B+A).print() # Concatenation

# print(f"GC content of second sequence in Strand object A is {A[1].gc()*100} %") 

# (~A).print() # Reverse complement

# Strand.random(20,seqtype='RNA',gcContent=0.25).print() # Create and print a random RNA sequence of length 20 with 25% GC content

# M = Multistrand('AGGC','TTAGTG')

# M.Strands.print()

Strand.polyN('rA',5).print()