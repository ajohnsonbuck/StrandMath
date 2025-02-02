function len = NA_length(seq)
    seq = erase(seq,{'+','b'});
    len = length(seq);
end