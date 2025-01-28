function [ntrials, results] = test_probe_set(ntrials,subseq_list,is_subset,ind,options,results,decrease_set)
    if mod(ntrials,1000)==0
        fprintf(1,'Trials completed: %3u \n',ntrials);
    end
    ambiguous = zeros(size(is_subset,1),1);
    if strcmp(decrease_set,'y')
        set_size = results.min_set-1;
    else
        set_size = results.min_set;
    end
    ind2 = randsample(length(ind),set_size,false);
    ind2 = enforce_required_sequences(ind2,options.required_sequence_ind); % Enforce required sequences, if any
    is_subset_2 = is_subset(:,ind(ind2));
    fail = 0;
    nmissed = sum(sum(is_subset_2,2)==0);
    [is_subset_3,~,idc] = unique(is_subset_2,'rows'); % Find unique rows
    for n = 1:length(idc)
        if sum(idc==idc(n))>1
            ambiguous(n)=1; % flag non-unique rows
        end
    end
    % nambiguous = size(is_subset_2,1)-size(is_subset_3,1);
    nambiguous = sum(ambiguous);
    if strcmp(decrease_set,'y')
        if nambiguous/size(is_subset_2,1) > options.max_fraction_ambiguous || nmissed/size(is_subset_2,1) > options.max_fraction_missed
            fail = 1;
        end
    else
        if nambiguous >= results.nambiguous_best || nmissed > results.nmissed_best
            fail = 1;
        end
    end
    if fail == 0
        results.min_set = set_size;
        results.subseqs_best = subseq_list(ind(ind2));
        results.nambiguous_best = nambiguous;
        results.ambiguous_best = ambiguous;
        results.nmissed_best = nmissed;
        results.is_subset_best = is_subset_2;
        fprintf(1,'Set of %3u subsequences found satisfying requirements. %3u ambiguous sequences (%.1f %%). %3u missed sequences (%.1f %%).\n',results.min_set,nambiguous,nambiguous/n*100,nmissed,nmissed/n*100);
    end
    results.fail = fail;
    ntrials = ntrials+1;
end