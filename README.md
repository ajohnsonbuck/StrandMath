# symboligo

## Description
Library of classes for the convenient, high-throughput handling of 
oligonucleotides including DNA, RNA, and LNA.  Unlike other publicly 
available tools, symboligo allows the rapid manipulation of **thousands
of sequences in parallel** through built-in class methods, implemented as intuitive
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

Strand objects can contain a single oligonucleotide sequence or an array of sequences, allowing for high-throughput operations.

They can be combined and allowed to interact through
symbolic operations like `+` (concatenation), `-` (reverse sequence), `*` 
(combinatorial hybridization), `.*` (element-wise hybridization), and `'` (reverse complement).

Sequences can also be concatenated in head-to-tail repeats by multiplying them by a positive integer.
(For example, `Strand('AAG')*20` generates 20 repeats of the DNA sequence AAG)

Strand objects must be initially created with a Strand() call, but afterwards can be manipulated symbolically.
Furthermore, if an operation is called involving a Strand and a string representing a nucleotide sequence, 
the string will be understood to represent a second Strand.  For example:

     A = Strand('AAT+CGG+CTG').toLNA;
     B = 'TTTTT' + A;

Strand objects can also be created by passing an `n x 1` cell array of char or string variables as arguments to Strand().  For example,

     N = Strand({'ATTG';'rCrArGrA';'GGAATTC'});

creates a 3-element Strand array containing the three DNA sequences specified.  Note that RNA nucleotides are specified with the `r` prefix.

Strand sequences are stored internally as cell arrays of char variables representing individual nucleotides, 
and can be accessed as the Strand.Sequence property (for instance, `N(2).Sequence{3}` would return the char `'rG'` in the above example).
Alternatively, Strand sequences can be printed into a more convenient format for viewing with the .print() method:

     N.print; % Prints the three sequences contained in N in the form 5'-NNN...N-3'

By passing the argument `'bare'` to the print method, any modification prefixes are omitted.  For instance:

     N(2).print('bare');

will return `5'-CAGA` in the above example.

## Sources of thermodynamic nearest-neighbor parameters
Predictions of Tm, free energy, enthalpy, and energy of hybridization use nearest-neighbor models and parameters published in the following papers:
- **RNA/DNA hybrids:** Sugimoto,N. et al., Biochemistry, 34, 11211
- **RNA/RNA hybrids:** Xia et al. Biochemistry 1998, 37, 42, 14719–14735
- **RNA dangling ends:** Serra, M.J. and Turner, D.H. (1995) Predicting Thermodynamic Properties of RNA. Methods Enzymol., 259, 242-261.
- **RNA mismatches:** Mathews et al. JMB Volume 288, Issue 5, 21 May 1999, Pages 911-940.; Xia, T., Mathews, D.H. and Turner, D.H. (1999) In Söll, D. G., Nishimura, S. and Moore, P. B. (eds.), Prebiotic Chemistry, Molecular Fossils, Nucleosides, and RNA. Elsevier, New York, pp. 21-47.; Mathews et al. (2004) Proc. Natl. Acad. Sci. USA, 101, 7287-7292 & Lu, et al. (2006) Nucleic Acids Res., 34 4912 - 4924.; 
- **DNA/DNA hybrids:** Allawi,H., SantaLucia,J.,Jr., Biochemistry, 36, 10581
- **DNA/DNA mismatches:** Peyret et al. Biochemistry 1999, 38, 12, 3468–3477; Hatim T. Allawi and John SantaLucia, Biochemistry 1997 36 (34), 10581-10594; Hatim T. Allawi, John SantaLucia. Volume 26, Issue 11, 1 June 1998, Pages 2694–2701; Allawi and SantaLucia, Biochemistry 1998, 37, 26, 9435–9444
- **DNA dangling ends:** Bommarito et al. (2000), Nucleic Acids Research 28(9), 1929-1934.
- **LNA/DNA hybrids:** McTigue et al. Biochemistry 2004, 43, 18, 5388–5405; Owczarzy,R et al., Biochemistry, 50, 9352

Parameters for LNA/RNA duplexes are approximated, as these parameters are unavailable empirically.
