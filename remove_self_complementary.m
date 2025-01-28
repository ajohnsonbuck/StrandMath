function seqs = remove_self_complementary(seqs,len)
    N = length(seqs);
    nsubseqs = length(seqs{1})-len+1;
    for m = 1:N
        n = N-m+1;
        seq = seqs{n};
        rc = reverse_complement(seq);
        removed = 0;
        for p = 1:nsubseqs
            for q = 1:nsubseqs
                if sum(seq(q:q+len-1)==rc(p:p+len-1))>=len && removed==0
                    % seq(q:q+len-1)
                    % rc(p:p+len-1)
                    % pause;
                    seqs(n)=[];
                    removed = 1;
                end
            end
        end
    end
end