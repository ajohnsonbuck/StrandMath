function fGC = gc_content(seq)
    nGC = 0;
    for n = 1:length(seq)
        if strcmpi(seq(n),'G') || strcmpi(seq(n),'C')
            nGC = nGC + 1;
        end
    end
    fGC = nGC/length(seq);
end