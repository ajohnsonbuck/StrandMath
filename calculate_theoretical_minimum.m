function theoretical_minimum = calculate_theoretical_minimum(seqs,options)

if options.stop_at_theoretical_minimum == 1
    theoretical_minimum = ceil(logN(size(seqs,1)+1,3)); % Theoretical minimum size of set, assuming ternary code
else
    theoretical_minimum = 0;
end

end

function y = logN(x,N)
    y = log10(x)/log10(N);
end
