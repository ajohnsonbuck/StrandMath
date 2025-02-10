function seqs_out = generate_reverse_motif(seqs)
    % Flips nearest-neighbor motif (top strand becomes bottom and
    % vice-versa). Currently only works for DNA.
    seqs2 = cell(size(seqs,1),2);
    seqs3 = cell(size(seqs,1),2);
    seqs_out = cell(size(seqs));
    for n = 1:length(seqs)
        seqs2(n,:) = strsplit(seqs{n,1},'/');
    end
    for n = 1:length(seqs)
        seqs3{n,1} = reverse(seqs2{n,2});
        seqs3{n,2} = reverse(seqs2{n,1});
    end
    for n = 1:length(seqs)
        seqs_out{n,1} = [seqs3{n,1},'/',seqs3{n,2}];
    end
end