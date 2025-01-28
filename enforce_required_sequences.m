function ind2 = enforce_required_sequences(ind,required_sequence_ind)
    ind2 = ind;
    % Enforce required sequences
    if isempty(required_sequence_ind) == 0
        for m = 1:length(required_sequence_ind)
            if sum(ind2==required_sequence_ind(m))==0
                replaced=0;
                while replaced==0
                    n=unidrnd(length(ind2));
                    if sum(required_sequence_ind==ind2(n))==0
                        ind2(n)=required_sequence_ind(m);
                        replaced = 1;
                    end
                end
            end
        end
    end
end