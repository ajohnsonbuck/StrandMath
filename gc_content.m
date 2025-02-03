function fGC = gc_content(seq)
    % Calculate fractional GC content of a sequence

    seq = parse_modifications(seq); % Remove modification prefixes from sequence

    nGC = 0;
    for n = 1:length(seq)
        if strcmpi(seq(n),'G') || strcmpi(seq(n),'C')
            nGC = nGC + 1;
        end
    end
    fGC = nGC/length(seq);
end