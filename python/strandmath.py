# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

from typing import Union, Iterable, Optional
import copy
import re

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

    def __repr__(self):
        if self.numel()==1:
            str1 = "".join(self.sequence[0])
            return f"Strand({str1!r})"
        else:
            return f"Strand({self.numel()} Sequences)"

    def __add__(self, other: Union['Strand',str]):
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
        for item in string:
            strlist.append(Strand.removeMods(item))
        if len(strlist) == 1:
            strlist = strlist[0]
        return strlist
            
    def removeMods(str1: str): # Remove modifications from string representation of sequence
        pattern = '|'.join(map(re.escape,Strand.Modlist))
        return re.sub(pattern,'', str1)
    
    def cleanString(str1: str): # Clean up sequence to remove termini markers and empty spaces
        pattern = '|'.join(map(re.escape,[" ", "5'-","-3'","5-","-3","5'","3'"])) 
        return re.sub(pattern,'',str1)
    
    def numel(self): # Number of sequences in Strand object
        return len(self.sequence)
    
    def len(self): # Number of nucleotides in each sequence
        L = []
        for item in self.sequence:
            L.append(len(item))
        if len(L)==1:
            L = L[0]
        return L      
        
    def print(self): # Display names and sequences of all items in Strand object
        strlist = self.string()
        if isinstance(strlist,str):
            strlist = [strlist]
        for (stritem,nameitem) in zip(strlist, self.name):
            if nameitem != "":
                print(nameitem)
            print("5'-" + stritem + "-3'")

# Example usage:
# A = Strand("rArGrCrT",name="strand1")
A = Strand(["rArGrCrT","GACCTA"],name=["strand1","strand2"])
A.print()

B = Strand('TTTT',name="T4")

C = B + A
C.print()

D = A + B
D.print()

# b = Strand("CGA")
# print(a + b)  # Strand('ATGCGA')
# print(a + [b, Strand("TTT")])  # [Strand('ATGCGA'), Strand('ATGTTT')]