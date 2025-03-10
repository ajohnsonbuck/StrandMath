# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

from typing import Union, Iterable, Optional
import re

class Strand:
    Modlist = ["+","b","r"]
    Nucleotides = ["A","C","G","T","U","a","c","g","t","u"]
    
    def __init__(self, sequence, name: Optional[str] = None):
        if isinstance(sequence,str):                    # if sequence is provided as string
            self.sequence = self.fromString(sequence) 
        elif isinstance(sequence,list):
            if all(isinstance(nt,str) for nt in sequence): # If sequence is provided as list of nucleotides
                self.sequence = self.fromString(sequence)

    def __repr__(self):
        str1 = "".join(self.sequence)
        return f"Strand({str1!r})"

    def __add__(self, other: Union['Strand', list]):
        if isinstance(other, Strand):
            return Strand(self.sequence + other.sequence)
        elif isinstance(other, list):
            return [self + item for item in other]
        return NotImplemented

    def __iter__(self):
        yield self  # Treat individual object as an iterable of one
    
    def fromString(self,str1: str):
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
        
    def string(self):
        return "".join(self.sequence)
        
    def bareString(self):
        return Strand.removeMods(self.string())            
    
    def removeMods(str1: str):
        pattern = '|'.join(map(re.escape,Strand.Modlist))
        return re.sub(pattern,'', str1)
    
    def cleanString(str1: str):
        pattern = '|'.join(map(re.escape,[" ", "5'-","-3'","5-","-3","5'","3'"])) 
        return re.sub(pattern,'',str1)
        
    def print(self):
        str1 = self.string()
        print("5'-" + str1 + "-3'")

# Example usage:
# A = Strand(["A","C","G","T"])
a = Strand("rArTrG")
# b = Strand("CGA")
# print(a + b)  # Strand('ATGCGA')
# print(a + [b, Strand("TTT")])  # [Strand('ATGCGA'), Strand('ATGTTT')]