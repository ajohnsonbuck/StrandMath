function rc = reverse_complement(seq)
modlist = {'b','r','+'};
cellInput = true;
if ~isa(seq,'cell')
    seq = {seq};
    cellInput = false;
end
rc = seq;
for p = 1:length(seq)
    entry = seq{p};
    entry = erase(entry,modlist);
    seqentry = entry;
    for m = 1:length(seqentry)
        n = length(seqentry)-m+1;
        base = seqentry(n);
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
        else
            comp = '-';
        end
            entry(m)=comp;
    end
    rc{p} = entry;
end
if cellInput == false %Convert back to char in case of single input
    rc = rc{1};
end
end