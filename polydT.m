function N = polydT(n)
% Create poly-T DNA sequence of length n
    N = n*NucleicAcid('T');
    N.Name = ['(dT)',num2str(n)];
end