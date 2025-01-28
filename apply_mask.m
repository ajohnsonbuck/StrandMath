function seqs = apply_mask(seqs,options)

for n = 1:length(seqs)
    for m=1:length(seqs{n})
        if options.mask(m)=='-'
              seqs{n}(m)='-';
        end
    end
end

end