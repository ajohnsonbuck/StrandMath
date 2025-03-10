function N = polyN(nt,n)
% Create homopolymer sequence of nucleotide nt with length n
    if isa(nt,'char') || isa(nt,'string')
        nt = erase(nt,'d'); % Remove 'd' prefix for DNA if provided
        N = n*Strand(char(nt));
        N.Name = ['(', nt,')',num2str(n)];
    end
end