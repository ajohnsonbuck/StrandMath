function seqs = apply_gc_cutoffs(seqs,gc_content_range)
    N = length(seqs);
    for m = 1:N
        n = N-m+1;
        gc = (sum(seqs{n}=='C') + sum(seqs{n}=='G'))/length(seqs{n});
        if  gc < gc_content_range(1) || gc > gc_content_range(2)
            seqs(n) = [];
        end
    end
end