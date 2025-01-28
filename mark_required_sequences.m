function options = mark_required_sequences(options,subseq_list,ind)

options.required_sequence_ind = zeros(0,1);
for i = 1:length(options.required_sequences)
    for j = 1:length(ind)
        if strcmp(subseq_list{ind(j)},options.required_sequences{i})
            options.required_sequence_ind = cat(1,options.required_sequence_ind,j);
        end
    end
end

end