% Find deltaG for all probes binding to all miRNA targets in a defined list

% NEED TO INCORPORATE MISMATCHES INTO CODE, INCLUDING TERMINAL MISMATCHES

T = 37;

[seqs_orig, seq_names] = walter_lab_miRNAs();

seqs = seqs_orig;
for n=1:length(seqs)
   seqs{n} = seqs{n}(12:end); % only g12-end 
   % seqs{n} = seqs{n}(13:end); % only g13-end 
end

% probe_targets = {'ATCG','TCAT', 'GGCT', 'TCAA', 'GGAC', 'GAAG', 'CCTC', 'GCAA', 'TGGC', 'ACCG', 'GTTG', 'GTAT', 'GTCC', 'TAGT', 'CTGC', 'TGTA'}; %2025-01-29 98.8% unambiguous
probe_targets = {'TGGC','GTTG','ATCG'};

probes = probe_targets;
for n = 1:length(probe_targets)
    probes{n} = reverse_complement(probe_targets{n});
end

dG = zeros(length(seq_names),length(probes));

for m = 1:length(seqs_orig)
    for n = 1:length(probes)
        [~, ~, dG(m,n)] = parameters_NA(add_modifications(probes{n},{'+','+','+','+'}),'target',seqs{m},'type','RNA','temperature',T,'model_overhangs','true'); % LNA parameters
    end
end

% imshow(-dG/max(max(-dG)),'InitialMagnification',2000);

dG = dG/1000;