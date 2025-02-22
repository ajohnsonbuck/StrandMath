# symboligo

## Description
Library of classes for the convenient, high-throughput handling of 
oligonucleotides including DNA, RNA, and LNA.  Unlike other publicly 
available tools, symboligo allows the intuitive handling of **thousands
of sequences in parallel** through built-in class methods, implemented as an intuitive
symbolic operations.

Current functionality includes:
1. Prediction of Tm, free energy,enthalpy, and entropy of hybridization
2. Thermodynamic parameters for DNA, RNA, DNA/RNA, LNA/DNA, and (through approximation) LNA/RNA duplexes
3. Displaying the longest predicted duplex between two sequences
4. Interconversion between DNA, RNA, and LNA sequences
5. Generation of complements, reverse sequences, and repeats/concatenated sequences
6. Random sequence generation with desired length and GC content
7. and more

## Installation
1. Clone or download the repository
2. Place the repository folder and its subfolders on your Matlab path

## Example usage
### Random sequence generation
     A = Strand('random', 20); % Generate a random DNA 20-mer

### Concatenation and reverse
     B = A + polydT(10); % Concatenate a (dT)10 sequence to the 3'-end of sequence A
     D = -B; % Flip sequence B 5'-to-3'

### Reverse complement and type conversion
     C = B'; % Create a new strand C that is the reverse complement of B
     C = B.reverseComplement; % Another equivalent way to get the reverse complement
     C = C.toRNA % Convert C to RNA

### Hybridize
     P = B * C; % Find base-pairing between B and C
     P.print; % Show base-pairing between B and C as well as the standard Gibbs free energy of hybridization

### Tm estimation
     Tm = P.estimateTm; % Estimate Tm for pair P at 1 M Na+
     Tm = estimate_Tm('ATAGCGCCTAAT','Na',0.1,'conc',1E-6); % Estimate Tm for specified sequence and its reverse complement at 100 mM Na+ and 1 uM oligo

### High-throughput hybridization
    load('validation_Sugimoto_etal_1995.mat'); % Load DNA/RNA validation set
    seqs = SugimotoTable2.seqs; % Store sequences from Table 2 as an array of Strand objects
    pairs = seqs .* (seqs.toDNA)'; % Perform element-wise hybridization of each sequence in Sugimoto et al Table 2 to its DNA reverse complement
    pairs.print; % Show all duplexes and their predicted standard Gibbs free energies
    Tm = pairs(5).Tm; % get the Tm of the fifth pair

## Further Description
The basic class of symboligo is the Strand.  

Strand objects can contain a single
oligonucleotide sequence or an array of sequences, allowing for high-throughput operations.

They can be combined and allowed to interact through
symbolic operations like + (concatenation), - (adding the reverse sequence), * 
(combinatorial hybridization), .* (element-wise hybridization), and ' (reverse complement).

Sequences can also be concatenated in head-to-tail repeats by multiplying them by a positive integer.
(For example, Strand('AAG')*20 generates 20 repeats of the sequence AAG)

Strand objects must be initially created with a Strand() call, but afterwards can be manipulated symbolically.
Furthermore, if an operation is called involving a Strand and a string representing a nucleotide sequence, 
the string will be understood to represent a second Strand.  For example:

     A = Strand('AAT+CGG+CTG').toLNA;
     B = 'TTTTT' + A;



