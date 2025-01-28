function sequences = human_only(sequences)
    ind = zeros(0,1);
    for n = 1:size(sequences,1)
        if contains(sequences{n,1},'Homo sapiens') || contains(sequences{n,1},'hsa-')
            ind = cat(1,ind,n);
        end
    end
    sequences = sequences(ind,:);
end