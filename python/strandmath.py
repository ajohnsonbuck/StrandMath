# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

from typing import Union, Iterable, Optional
import copy
import re
import random
import numpy as np

class Strand:
    Modlist = ["+","b","r"]
    Nucleotides = ["A","C","G","T","U","a","c","g","t","u"]

    def __init__(self, sequence, name = ""):
        if isinstance(sequence,str):                    # if sequence is provided as string, interpret as single sequence
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

    def __getitem__(self,ind):
        name = self.name[ind]
        seq = self.sequence[ind]
        if isinstance(name,str):
            name = [name]
            seq = [seq]
        newStrand = copy.deepcopy(self)
        newStrand.name = name
        newStrand.sequence = seq
        return newStrand

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
            newStrand = copy.deepcopy(other)
            for ind in range(other.numel()):
                newStrand.sequence[ind] = self.sequence[0] + other.sequence[ind]
                if len(self.name[0])>0:
                    newStrand.name[ind] = f"{self.name[0]} + {other.name[ind]}"
            return newStrand
        elif other.numel()==1:
            newStrand = copy.deepcopy(self)
            for ind in range(self.numel()):
                newStrand.sequence[ind] = self.sequence[ind] + other.sequence[0]
                if len(other.name[0])>0:
                    newStrand.name[ind] = f"{self.name[ind]} + {other.name[0]}"
            return newStrand
        return NotImplemented
    
    def __radd__(self, other: str):
        other = Strand(other)
        return other + self
    
    def __neg__(self):
        return copy.deepcopy(self).reverse()
    
    def __sub__(self,other):
        if isinstance(other,str):
            other = Strand(other)
        return self + (-other)
    
    def __rsub__(self,other):
        return other + (-self)
    
    def __invert__(self):
        return self.reverseComplement()
    
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
        
    def bareString(self): # String representation, stripped of modifications
        strlist = []
        string = self.string()
        if isinstance(string,str):
            string = [string]
        for item in string:
            strlist.append(Strand.removeMods(item))
        if len(strlist) == 1:
            strlist = strlist[0]
        return strlist
    
    def bareSequence(self):
        bareSeq = copy.deepcopy(self.sequence)
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
        out = copy.deepcopy(self)
        out.sequence = out.bareSequence()
        out.sequence = [[s.replace('U','T') for s in row] for row in out.sequence]
        return out
        
    def toRNA(self): # Convert to DNA sequence
        out = copy.deepcopy(self)
        out.sequence = out.bareSequence()
        out.sequence = [[s.replace('T','U') for s in row] for row in out.sequence]
        out.sequence = [['r' + s for s in row] for row in out.sequence]
        return out
    
    def toLNA(self): # Convert to LNA sequence
        out = copy.deepcopy(self)
        out.sequence = out.bareSequence()
        out.sequence = [['+' + s for s in row] for row in out.sequence]
        return out
    
    def toBNA(self): # Convert to BNA sequence
        out = copy.deepcopy(self)
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
        out = copy.deepcopy(self)
        out.sequence = [row[::-1] for row in out.sequence]
        out.name = [name + '_reverse' for name in out.name]
        return out
    
    def reverseComplement(self): # Create reverse complement of all sequences in Strand
        out = copy.deepcopy(self)
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

    def crop(self,ind: list=(0,'end')):
        if any(not(isinstance(item, int)) and not(item=='end') for item in ind):
            raise TypeError("Indices for Strand.crop must be a list of integers or 'end'")
        if ind[1]=='end':
            ind[1] = np.inf
        if not(ind[1]>ind[0]):
            raise ValueError("Second index of Strand.crop must be larger than the first index")
        out = copy.deepcopy(self)
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
    
    def random(length: int, seqtype: str='DNA', gcContent: float=0.5,name: str=''):
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

# Example usage:
A = Strand(["rArGrCrU","GACCTA"],name=["Sequence1","Sequence2"]) # Create a Strand object with two sequences
A.print() # Show all sequences

B = Strand('TTTTT',name="T5") # Create a strand object with a single sequence (5T linker)

(B+A).print() # Concatenation

print(f"GC content of second sequence in Strand object A is {A[1].gc()[0]*100} %") 

(~A).print() # Reverse complement

Strand.random(20,seqtype='RNA',gcContent=0.25).print() # Create and print a random RNA sequence of length 20 with 25% GC content
