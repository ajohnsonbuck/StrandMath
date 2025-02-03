function [seq2,mods] = parse_modifications(seq)
    modlist = {'+','b','r'};
    seq2 = erase(seq,modlist); % Remove modification prefixes from sequences
    mods = cell(1,length(seq2));
    c = 1;
    for n = 1:length(seq)
        if ismember(seq(n),modlist)
            mods{c} = seq(n);
        else
            c = c+1;
        end
    end
end