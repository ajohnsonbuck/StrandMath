function [seq2] = add_modifications(seq,mods)
    seq2 = '';
    for n = 1:length(seq)
        if ~isempty(mods{n})
            seq2 = strcat(seq2,mods{n});
        end
        seq2 = strcat(seq2,seq(n));
    end
end