function N = walter_lab_miRNAs()
    seqs = cell(7,1);
    seq_names = cell(size(seqs));
    seq_names{1} = 'hsa-miR-141';
    seqs{1} = rna2dna('UAACACUGUCUGGUAAAGAUGG');
    seq_names{2} = 'hsa-miR-375';
    seqs{2} = rna2dna('UUUGUUCGUUCGGCUCGCGUGA');
    seq_names{3} = 'cel-miR-39';
    seqs{3} = rna2dna('UCACCGGGUGUAAAUCAGCUUG');
    seq_names{4} = 'hsa-let-7a';
    seqs{4} = rna2dna('UGAGGUAGUAGGUUGUAUAGUU');
    seq_names{5} = 'hsa-miR-29';
    seqs{5} = rna2dna('UAGCACCAUCUGAAAUCGGUUA');
    seq_names{6} = 'hsa-miR-16';
    seqs{6} = rna2dna('UAGCAGCACGUAAAUAUUGGCG');
    seq_names{7} = 'hsa-miR-21';
    seqs{7} = rna2dna('UAGCUUAUCAGACUGAUGUUGA');
    

    N = NucleicAcid(seqs,'name',seq_names).toRNA;

end