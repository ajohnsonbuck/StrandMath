function y = construct_nucleic_acid_subseqs(K)
    x = 'ATGC';
    %// Sample data

%// Create all possible permutations (with repetition) of letters stored in x
C = cell(K, 1);             %// Preallocate a cell array
[C{:}] = ndgrid(x);         %// Create K grids of values
y = cellfun(@(x){x(:)}, C); %// Convert grids to column vectors
y = [y{:}];                 %// Obtain all permutations
y  = cellstr(y);
end