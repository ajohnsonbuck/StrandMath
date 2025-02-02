% Given a small set of probes (e.g., 16-20), find the smallest set that can
% cover the 7 miRNAs in our possession unambiguously

[seqs_orig, seq_names] = walter_lab_miRNAs();

% Probe list as of 2025-01-29
subseq_list = {'GTTG',	'GGAG',	'ATCG',	'GTAT',	'TCTG',	'ACTC',	'CTGG',	'GTGT',	'TCAA',	'AAGG',	'TAAC',	'TCCA',	'TTGC',	'TCCC',	'GCAC',	'TGGC'};

% subseq_list = {'TGGC', 'CTGA', 'TTGC', 'AATC'};

options.required_sequences = {'TGGC', 'GTTG', 'ATCG'};

options.probe_size = 4; % Probe size to use
options.max_trials = 10000; % Maximum number of probe sets to try
options.min_representation = 0; % Minimum number of target sequences for a probe to bind to
options.max_fraction_ambiguous = 0.5; % Maximum fraction of target sequences that cannot be distinguished for a probe set to be valid
options.max_fraction_missed = 0.5; % Maximum fraction of target sequences that can be missed for a probe set to be valid
options.starting_set = length(subseq_list); % Starting size of probe sequence set
options.target_set_size = 4; % Desired size of probe sequence set
options.find_subsubseqs = 1; % Look for N-1 complementary regions within N-residue sub-sequences
options.gc_content_range = [0 1]; % inclusive range of allowed gc content for probe sequences
options.remove_self_complementary_probes = 1; % Remove self-complementary probes from candidate pool
options.stop_at_theoretical_minimum = 0; % Stop trying to reduce probe set size at theoretical minimum size for ternary code
options.eliminate_duplicate_targets = 1; % Eliminate any duplicate target sequences in the target set
options.max_self_complementarity = 2; % Maximum allowed self-complementarity

% options.mask = 'nnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnn';
options.mask = '------------nnnnnnnnnnnnnnnnnnnnnnnnn'; % g13-g22 only

seqs = seqs_orig;
for n = 1:length(seqs)
    seqs{n} = rna2dna(seqs{n}); % Convert RNA sequences to DNA equivalents
end

seqs = apply_mask(seqs,options); % Mask unused parts of target sequences

is_subset = determine_subsequence_presence(seqs,subseq_list,options); % Construct matrix of presence/absence of probe-binding sequences/subsequences

ind = 1:length(subseq_list); % Include all subsequences

options = mark_required_sequences(options,subseq_list,ind); % Mark positions of any required probe target sequences

theoretical_minimum = calculate_theoretical_minimum(seqs,options); % Calculate theoretical minimum number of probes

% Decrease probe number as much as feasible
ntrials = 0;
results = initialize_results(options);
while ntrials<options.max_trials && results.min_set > max([theoretical_minimum, options.target_set_size])
    [ntrials, results] = test_probe_set(ntrials,subseq_list,is_subset,ind,options,results,'y');
end

% Optimize results with lowest feasible number of probes
ntrials = 0;
while ntrials<options.max_trials && results.nambiguous_best > 0 && results.nmissed_best > 0
    [ntrials, results] = test_probe_set(ntrials,subseq_list,is_subset,ind,options,results,'n');
end

disp('Done.');

table1 = results_table_rna(seqs_orig,seq_names,results);