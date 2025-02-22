# symboligo

## Description
Library for the high-throughput handling of 
oligonucleotides including DNA, RNA, and LNA.  

Unlike other freely available tools, symboligo allows the rapid manipulation of thousands
of sequences in parallel through built-in class methods, implemented as intuitive
symbolic operations.

Current functionality includes:
1. Prediction of Tm, standard Gibbs free energy, enthalpy, and entropy of hybridization
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
### Creating an oligo (Strand)
     N = Strand('ACAGAGATTAGAAACCCAG');      % Create a nucleic acid Strand object
     N.print;      % Show result

          5'-ACAGAGATTAGAAACCCAG-3'

### Random sequence generation
     A = Strand('random', 'length', 20, 'GCcontent', 0.4);      % Generate a random DNA 20-mer with 40% GC content
     A.print;

          5'-TTCATTTTCTCCAAGGAGCT-3'

### Concatenation
     B = A + polydT(10);      % Concatenate a (dT)10 sequence to the 3'-end of sequence A
     B.print;

          5'-TTCATTTTCTCCAAGGAGCTTTTTTTTTTT-3'

### Reverse complement
     C = A';      % Create a new strand C that is the reverse complement of A
     C = A.reverseComplement;      % Another way to get the same result
     H = toRNA(A + 'GAAA' + A');      % Create a hairpin-forming sequence with the RNA equivalent of A and its reverse complement separated by a GAAA tetraloop
     H.print;

          5'-rUrUrCrArUrUrUrUrCrUrCrCrArArGrGrArGrCrUrGrArArArArGrCrUrCrCrUrUrGrGrArGrArArArArUrGrArA-3'

### Hybridize
     P = B * C;      % Find base-pairing between B and C
     P.print;      % Show base-pairing between B and C as well as the standard Gibbs free energy of hybridization

          5'- T T C A T T T T C T C C A A G G A G C T T T T T T T T T T T-3'
              | | | | | | | | | | | | | | | | | | | |                    
          3'- A A G T A A A A G A G G T T C C T C G A                    -5'

          dG0 = -23.4 kcal/mol

### Tm estimation
Standard Tm prediction conditions are 0.2 μM oligo, 1 M Na+, and 0 M Mg2+.  However, other conditions can be specified.

     Tm = P.estimateTm      % Estimate Tm for pair P at 1 M Na+ and 0.2 μM oligo.

          Tm = 67.4881
          
     Tm = estimate_Tm('ATAGCGCCTAAT','Na',0.1,'conc',1E-6)      % Estimate Tm, in degrees Celsius, for specified sequence and its reverse complement at 100 mM Na+ and 1 μM oligo

          Tm = 40.3000

### High-throughput hybridization and thermodynamics
    load('validation_Sugimoto_etal_1995.mat');      % Load DNA/RNA validation set
    N = SugimotoTable2.seqs;      % Store sequences from Table 2 as an array of Strand objects
    P= N .* (N.toDNA)';      % Perform element-wise hybridization of each of the 64 sequences in Sugimoto et al Table 2 to its DNA reverse complement
    Tm = P.estimateTm('conc', 100-6);      % Estimate Tm for all 64 duplexes at 100 μM oligo and 1 M Na+
    P(1:4).print;      % Show the first four RNA/DNA duplexes and their predicted standard Gibbs free energies

     5'- C G G C T-3'
         | | | | |
     3'-rGrCrCrGrA-5'
   

     dG0 = -5.2 kcal/mol

   
     5'- A G C C G-3'
         | | | | |
     3'-rUrCrGrGrC-5'
   

      dG0 = -5.1 kcal/mol

   
     5'- C C A C C-3'
         | | | | |
     3'-rGrGrUrGrG-5'
   

      dG0 = -5.4 kcal/mol

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

     A = Strand('AAT+CGG+CTG');
     B = 'TTTTT' + A;      % 'TTTTT' is understood as Strand('TTTTT') since it is being added to another Strand
     print(B)

          5'-TTTTTAAT+CGG+CTG-3'          

Strand objects can also be created by passing an `n x 1` cell array of char or string variables as arguments to Strand().  For example,

     N = Strand({'ATTG';'rCrArGrA';'GGAA+TTC'},'name',{'Sequence 1', 'Sequence 2', 'Sequence 3'});

creates a 3-element Strand array containing the DNA sequence ATTG, the RNA sequence CAGA, and the DNA sequence containing one LNA residue GGAA+TTC. The strands
are also named by passing 'name' as an argument to Strand. Strand names can be accessed, for instance, by:

     N(2).Name

          ans =

              'Sequence 2'

Note that RNA nucleotides are specified with the `r` prefix, and LNA residues with the `+` prefix.

Strand sequences are stored internally as cell arrays of char variables representing individual nucleotides, 
and can be accessed as the Strand.Sequence property (for instance, `N(2).Sequence{3}` would return the char `'rG'` in the above example).
Alternatively, Strand sequences can be printed into a more convenient format for viewing with the .print() method:

     N.print;

          Sequence 1
          5'-ATTG-3'

          Sequence 2
          5'-rCrArGrA-3'

          Sequence 3
          5'-GGAA+TTC-3'
     

By passing the argument `'bare'` to the print method, any modification prefixes are omitted.  For instance:

     N(2).print('bare');

          Sequence 2
          5'-CAGA-3'

## Sources of thermodynamic nearest-neighbor parameters
Predictions of Tm, free energy, enthalpy, and energy of hybridization use nearest-neighbor models and parameters published in the following papers:
- **RNA/DNA hybrids:** Sugimoto, N. et al. (1995) Biochemistry, 34, 11211.
- **RNA/RNA hybrids:** Xia et al. (1998) Biochemistry, 37, 42, 14719–14735.
- **RNA dangling ends:** Serra, M.J. and Turner, D.H. (1995) Methods Enzymol., 259, 242-261.
- **RNA mismatches:** Mathews et al. (1999) JMB Volume 288, Issue 5, 911-940; Xia, T., Mathews, D.H. and Turner, D.H. (1999) In Söll, D. G., Nishimura, S. and Moore, P. B. (eds.), Prebiotic Chemistry, Molecular Fossils, Nucleosides, and RNA. Elsevier, New York, pp. 21-47.; Mathews et al. (2004) Proc. Natl. Acad. Sci. USA, 101, 7287-7292; Lu, et al. (2006) Nucleic Acids Res., 34 4912-4924.
- **DNA/DNA hybrids:** Allawi, H., SantaLucia, J., Jr. (1997) Biochemistry, 36, 10581.
- **DNA/DNA mismatches:** Peyret et al. (1999) Biochemistry, 38, 12, 3468–3477; Allawi, H., SantaLucia, J., Jr. (1997) Biochemistry 36 (34), 10581-10594; Allawi, H., SantaLucia, J., Jr. (1998) Nucleic Acids Res. 26 (11), 2694–2701; Allawi, H., SantaLucia, J., Jr., Biochemistry 1998, 37, 26, 9435–9444.
- **DNA dangling ends:** Bommarito et al. (2000), Nucleic Acids Res. 28(9), 1929-1934.
- **LNA/DNA hybrids:** McTigue et al. Biochemistry 2004, 43(18), 5388–5405; Owczarzy, R. et al. (2011) Biochemistry, 50, 9352.

Parameters for LNA/RNA duplexes are approximated, as these parameters are unavailable empirically.

### Validation
By running the included script `validate_NN_parameters.m`, you can confirm that the predictions of the model match those shown in Sugimoto et al. 1995 (RNA/DNA), Oczarzy et al. 2011 (DNA/DNA, DNA mismatches, LNA, and LNA mismatches) and Xia et al. 1998 (RNA/RNA) to within 0.1 degree Celsius.
