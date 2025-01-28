function is_subset = determine_subsequence_presence(seqs,subseq_list,options)

nsubseqs = length(subseq_list);

is_subset = zeros(size(seqs,1),length(subseq_list)); % Initialize array of subsequence presence/absence in each sequence

for n = 1:nsubseqs
    if mod(n,100)==0
        fprintf(1,'Determined presence of %4u of %4u subsequences in sequence set.\n',n,nsubseqs);
    end
    for m = 1:length(seqs)
        is_subset(m,n) = contains(seqs{m,1},subseq_list{n});
    end
end

% Find sub-subsequences (optional)
if options.find_subsubseqs == 1
    len = options.subsubseq_size;
    m = length(subseq_list{1})-len+1; % number of subsubsequences
    fprintf(1,'Finding which sub-subsequences are present in each target sequence...\n');
    % Make second pass to identify partial complementarity
    for n = 1:nsubseqs
        % determine subsubsequences
        for p = 1:m
            subsubseq = subseq_list{n}(p:p+len-1);
            % Check if subsubseq is a substring of each target; if so, label
            % it a 2
            for q = 1:size(seqs,1)
                if is_subset(q,n) == 0
                    if contains(seqs{q,1},subsubseq)==1
                        is_subset(q,n) = 2;
                    end
                end
            end
        end
    end
end

end