function [dH0, dS0, dG] = parameters_NA(seq,type,T)
    % seq = PROBE sequence (e.g., DNA if RNA-DNA duplex, or BNA if RNA-BNA
    % duplex); target sequence is implied
    % type = target type (DNA or RNA)

    if strcmp(type,'RNA')
        % RNA-DNA parameters from From Sugimoto,N. et al., Biochemistry, 34, 11211, per IDT Oligoanalyzer
        % 2025-01-31: Confirmed identical results for several oligo sequences in Sugimoto
        % et al. Table 2.
        seqs = {'AA',   'AC',	'AG',   'AT',   'CA',   'CC',   'CG',   'CT',	'GA',   'GC',   'GG',   'GT',   'TA',   'TC',   'TG',   'TT'};
        dH0s = -[11.5,	7.8,    7,	    8.3,    10.4,   12.8,   16.3,   9.1,    8.6,    8,	    9.3,    5.9,    7.8,    5.5,	9,     7.8];
        dS0s = -[36.4,  21.6,   19.7,   23.9,   28.4,   31.9,   47.1,   23.5,   22.9,   17.1,   23.2,   12.3,   23.2,   13.5,   26.1,  21.9];
        dH0 = 1.9; % Initiation
        dS0 = -3.9; % Initiation
    elseif strcmp(type, 'DNA')
        % RNA-DNA parameters from Allawi,H., SantaLucia,J.,Jr., Biochemistry, 36, 10581
        seqs = {'AA',   'AC',   'AG',   'AT',  'CA',   'CC',   'CG',   'CT',   'GA',   'GC',   'GG',   'GT',   'TA',   'TC',   'TG',   'TT'};
        dH0s = -[7.9,   8.4,    7.8,    7.2,    8.5,    8.0,    10.6,   7.8,    8.2,    9.8,    8.0,    8.4,    7.2,    8.2,    8.5,    7.9]; % in kcal/mol
        dS0s = -[22.2,  22.4,   21.0,   20.4,   22.7,   19.9,   27.2,   21.0,   22.2,   24.4,   19.9,   22.4,   21.3,   22.2,   22.7,   22.2]; % in cal/mol/K
        % Isolated LNAs
        seqs = horzcat(seqs,{'+AA', '+AC',	'+AG',	'+AT',	'+CA',	'+CC',	'+CG',	'+CT',	'+GA',	'+GC',	'+GG',	'+GT',	'+TA',	'+TC',	'+TG',	'+TT', 'A+A',	'A+C',	'A+G',	'A+T',	'C+A',	'C+C',	'C+G',	'C+T',	'G+A',	'G+C',	'G+G',	'G+T',	'T+A',	'T+C',	'T+G',	'T+T'});
        % Isolated LNAs
        ddH1 = dH0s + [0.707,	1.131,	0.264,	2.282,	1.049,	2.096,	0.785,	0.708,	3.162,	-0.36,	-2.844,	-0.212,	-0.046,	1.893,	-1.54,	1.528];
        ddS1 = dS0s + [2.477,	4.064,	2.613,	7.457,	4.32,	7.966,	3.709,	4.175,	10.544,	-0.251,	-6.68,	0.073,	1.562,	6.685,	-3.044,	5.298];
        ddH2 = dH0s + [0.992,	2.89,	-1.2,	1.816,	1.358,	2.063,	-0.276,	-1.671,	0.444,	-0.925,	-0.943,	-0.635,	1.591,	0.609,	2.165,	2.326];
        ddS2 = dS0s + [4.065,	10.576,	-1.826,	6.863,	4.367,	7.565,	-0.718,	-4.07,	2.898,	-1.111,	-0.933,	-0.342,	5.281,	3.169,	7.163,	8.051];
        dH0s = cat(2,dH0s,ddH1,ddH2);
        dS0s = cat(2,dS0s,ddS1,ddS2);

        % Initiation and symmetry
        [dH0, dS0] = DNA_initiation(seq);
    end

    n = 1;
    while n < length(seq) % Iterate through all nearest neighbor sequences and add to deltaH and deltaS
        nn = nearest_neighbor(seq,n);
        if NA_length(nn) == 2
        dH0 = dH0 + dH0s(strcmpi(seqs,nn));
        dS0 = dS0 + dS0s(strcmpi(seqs,nn));
        end
        if strcmp(nn(1),'+') || strcmp(nn(1),'b')
            n = n+2;
        else
            n = n+1;
        end
    end
    dH0 = dH0*1000; % Convert to cal/mol

    % if strcmp(type, 'RNA-BNA')
    %     dS0 = dS0+4.9*length(seq);  % Assume dS0 decreases by 3.33 for each BNA modification -- based on estimate that deltaG drops by ~1000 cal/mol for each increase in ~6 deg C to the melting temp
    % end

    dG = dH0 - (273.15+T)*dS0;
end

function [dH0, dS0] = DNA_initiation(seq)
    seq2 = erase(seq,{'+','b'});
    if strcmpi(seq2(1),'G') || strcmpi(seq2(1),'C')
        dH0 = 0.1;
        dS0 = -2.8;
    else
        dH0 = 2.3;
        dS0 = 4.1;
    end
    if strcmpi(seq2(end),'G') || strcmpi(seq2(end),'C')
        dH0 = dH0 + 0.1;
        dS0 = dS0 - 2.8;
    else
        dH0 = dH0 + 2.3;
        dS0 = dS0 + 4.1;
    end
    if strcmpi(seq, reverse_complement(seq))
    dS0 = dS0 - 1.4; % Symmetry correction for self-complementary sequences 
    end
end

function nn = nearest_neighbor(seq,m)
    n = m+1;
    nn = seq(m:n);
    while NA_length(nn) < 2 && n<length(seq)
        n = n+1;
        nn = seq(m:n);
    end
end