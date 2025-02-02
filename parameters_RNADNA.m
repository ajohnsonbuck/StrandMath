function [dH0, dS0, dG] = parameters_RNADNA(rseq,T)
    % From Sugimoto,N. et al., Biochemistry, 34, 11211, per IDT Oligoanalyzer
    % 2025-01-31: Confirmed identical results for several oligo sequences in Sugimoto
    % et al. Table 2.

    % Parameter table
    seqs = {'AA','AC','AG','AU','CA','CC','CG','CU','GA','GC','GG','GU','UA','UC','UG','UU'};
    dH0s = -[7.8, 5.9, 9.1, 8.3, 9.0, 9.3, 16.3, 7.0, 5.5, 8.0, 12.8, 7.8, 7.8, 8.6, 10.4, 11.5]; % in kcal/mol
    dS0s = -[21.9, 12.3, 23.5, 23.9, 26.1, 23.2, 47.1, 19.7, 13.5, 17.1, 31.9, 21.6, 23.2, 22.9, 28.4, 36.4]; % in cal/mol/K

    dH0 = 1.9; % Initiation
    dS0 = -3.9; % Initiation
    for n = 1:length(rseq)-1 % Iterate through all nearest neighbor sequences and add to deltaH and deltaS
        nn = rseq(n:n+1);
        dH0 = dH0 + dH0s(strcmp(seqs,nn));
        dS0 = dS0 + dS0s(strcmp(seqs,nn));
    end
    dH0 = dH0*1000; % Convert to cal/mol

    dG = dH0 - (273.15+T)*dS0;
end