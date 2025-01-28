function rc = reverse_complement(seq)
    rc = seq;
    for m = 1:length(seq)
        n = length(seq)-m+1;
        base = seq(n);
        if strcmp(base,'C')
            comp = 'G';
        elseif strcmp(base,'G')
            comp = 'C';
        elseif strcmp(base,'T')
            comp = 'A';
        elseif strcmp(base, 'U')
            comp = 'A';
        elseif strcmp(base,'A')
            comp = 'T';
        end
        rc(m)=comp;
    end
end