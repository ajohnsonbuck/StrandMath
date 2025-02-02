function rc = reverse_complement(seq)
    rc = seq;
    for m = 1:length(seq)
        n = length(seq)-m+1;
        base = seq(n);
        if strcmpi(base,'C')
            comp = 'G';
        elseif strcmpi(base,'G')
            comp = 'C';
        elseif strcmpi(base,'T')
            comp = 'A';
        elseif strcmpi(base, 'U')
            comp = 'A';
        elseif strcmpi(base,'A')
            comp = 'T';
        end
        rc(m)=comp;
    end
end