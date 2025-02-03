function [dH0, dS0, dG] = parameters_NA(seq, varargin)
    % seq = PROBE sequence (e.g., DNA if RNA-DNA duplex, or BNA if RNA-BNA
    % duplex); target sequence is implied
    % type = target type (DNA or RNA)
    %
    % 2025-02-02 - To do:
    % Incorporate explicit dangling end treatment
    % Incorporate explicit mismatch treatment
    % Refine and test RNA-LNA/RNA-BNA parameters

    args = varargin;

    [type, T, mismatches, overhangs, startpos, endpos] = parse_input(seq,args);

    if strcmp(type,'RNA')
        % RNA-DNA parameters from From Sugimoto,N. et al., Biochemistry, 34, 11211, per IDT Oligoanalyzer
        % 2025-01-31: Confirmed identical results for several oligo sequences in Sugimoto et al. Table 2.
        seqs = {'AA',   'AC',	'AG',   'AT',   'CA',   'CC',   'CG',   'CT',	'GA',   'GC',   'GG',   'GT',   'TA',   'TC',   'TG',   'TT'};
        dH0s = -[11.5,	7.8,    7,	    8.3,    10.4,   12.8,   16.3,   9.1,    8.6,    8,	    9.3,    5.9,    7.8,    5.5,	9,     7.8];
        dS0s = -[36.4,  21.6,   19.7,   23.9,   28.4,   31.9,   47.1,   23.5,   22.9,   17.1,   23.2,   12.3,   23.2,   13.5,   26.1,  21.9];
        % Estimated BNA parameters, based on Abdur Rahman, S. M.; Seki, S.; Obika, S.; Yoshikawa, H.; Miyashita, K.; Imanishi, T. Design, Synthesis, and Properties of 2′,4′-BNANC: A Bridged Nucleic Acid Analogue. J. Am. Chem. Soc. 2008, 130 (14), 4886−4896.
        % Assume only entropy changes (prob. not accurate); adjusted dS0 to maximize agreement with Tms in Table 2 in Rahman et al. at 100 mM Na and 1 uM oligo
        seqs = horzcat(seqs, {'bAA',   'bAC',	'bAG',   'bAT',   'bCA',   'bCC',   'bCG',   'bCT',	'bGA',   'bGC',   'bGG',   'bGT',   'bTA',   'bTC',   'bTG',   'bTT', 'AbA',   'AbC',	'AbG',   'AbT',   'CbA',   'CbC',   'CbG',   'CbT',	'GbA',   'GbC',   'GbG',   'GbT',   'TbA',   'TbC',   'TbG',   'TbT', 'bAbA',   'bAbC',	'bAbG',   'bAbT',   'bCbA',   'bCbC',   'bCbG',   'bCbT',	'bGbA',   'bGbC',   'bGbG',   'bGbT',   'bTbA',   'bTbC',   'bTbG',   'bTbT'});
        ddS1 = dS0s + ones(1,16)*2.1;
        ddS2 = dS0s + ones(1,16)*2.1;
        ddS3 = dS0s + ones(1,16)*4.9;
        % dS0s = cat(2,dS0s,ddS1,ddS2,ddS3);
        % dH0s = cat(2,dH0s,dH0s,dH0s,dH0s);
        % Assume identical parameters for RNA-LNA hybridization as for RNA-BNA (reasonable based on Table 2 of Rahman et al.)
        seqs = horzcat(seqs, {'+AA',   '+AC',	'+AG',   '+AT',   '+CA',   '+CC',   '+CG',   '+CT',	'+GA',   '+GC',   '+GG',   '+GT',   '+TA',   '+TC',   '+TG',   '+TT', 'A+A',   'A+C',	'A+G',   'A+T',   'C+A',   'C+C',   'C+G',   'C+T',	'G+A',   'G+C',   'G+G',   'G+T',   'T+A',   'T+C',   'T+G',   'T+T', '+A+A',   '+A+C',	'+A+G',   '+A+T',   '+C+A',   '+C+C',   '+C+G',   '+C+T',	'+G+A',   '+G+C',   '+G+G',   '+G+T',   '+T+A',   '+T+C',   '+T+G',   '+T+T'});
        ddS4 = dS0s + ones(1,16)*2.1;
        ddS5 = dS0s + ones(1,16)*2.1;
        ddS6 = dS0s + ones(1,16)*4.9;
        dS0s = cat(2,dS0s,ddS1,ddS2,ddS3,ddS4,ddS5,ddS6);
        dH0s = cat(2,dH0s,dH0s,dH0s,dH0s,dH0s,dH0s,dH0s);
        dH0 = 1.9; % Initiation
        dS0 = -3.9; % Initiation
    elseif strcmp(type, 'DNA')
        % DNA-DNA parameters from Allawi,H., SantaLucia,J.,Jr., Biochemistry, 36, 10581
        seqs = {'AA',   'AC',   'AG',   'AT',  'CA',   'CC',   'CG',   'CT',   'GA',   'GC',   'GG',   'GT',   'TA',   'TC',   'TG',   'TT'};
        dH0s = -[7.9,   8.4,    7.8,    7.2,    8.5,    8.0,    10.6,   7.8,    8.2,    9.8,    8.0,    8.4,    7.2,    8.2,    8.5,    7.9]; % in kcal/mol
        dS0s = -[22.2,  22.4,   21.0,   20.4,   22.7,   19.9,   27.2,   21.0,   22.2,   24.4,   19.9,   22.4,   21.3,   22.2,   22.7,   22.2]; % in cal/mol/K
        % Parameters for isolated LNA residues from McTigue et al. Biochemistry 2004, 43, 18, 5388–5405
        seqs = horzcat(seqs,{'+AA', '+AC',	'+AG',	'+AT',	'+CA',	'+CC',	'+CG',	'+CT',	'+GA',	'+GC',	'+GG',	'+GT',	'+TA',	'+TC',	'+TG',	'+TT', 'A+A',	'A+C',	'A+G',	'A+T',	'C+A',	'C+C',	'C+G',	'C+T',	'G+A',	'G+C',	'G+G',	'G+T',	'T+A',	'T+C',	'T+G',	'T+T'});
        ddH1 = dH0s + [0.707,	1.131,	0.264,	2.282,	1.049,	2.096,	0.785,	0.708,	3.162,	-0.36,	-2.844,	-0.212,	-0.046,	1.893,	-1.54,	1.528];
        ddS1 = dS0s + [2.477,	4.064,	2.613,	7.457,	4.32,	7.966,	3.709,	4.175,	10.544,	-0.251,	-6.68,	0.073,	1.562,	6.685,	-3.044,	5.298];
        ddH2 = dH0s + [0.992,	2.89,	-1.2,	1.816,	1.358,	2.063,	-0.276,	-1.671,	0.444,	-0.925,	-0.943,	-0.635,	1.591,	0.609,	2.165,	2.326];
        ddS2 = dS0s + [4.065,	10.576,	-1.826,	6.863,	4.367,	7.565,	-0.718,	-4.07,	2.898,	-1.111,	-0.933,	-0.342,	5.281,	3.169,	7.163,	8.051];
        % Parameters for multiple consecutive LNA residues from Owczarzy,R et al., Biochemistry, 50, 9352
        % Tested several perfectly matched DNA-only and LNA-containing probes from Table 4 of Owczarzy et al., and found agreement within 0.1 deg C.
        seqs = horzcat(seqs,{'+A+A',   '+A+C',	'+A+G',   '+A+T',   '+C+A',   '+C+C',   '+C+G',   '+C+T',	'+G+A',   '+G+C',   '+G+G',   '+G+T',   '+T+A',   '+T+C',   '+T+G',   '+T+T'});
        ddH3 = dH0s - [2.091, 2.989, 4.993, 7.503, 5.677, 7.399, 3.958, 7.937, 5.759, 6.309, 5.022, 8.961, 3.118, 0.966, 1.546, 2.519];
        ddS3 = dS0s + [-4.975, -6.563, -10.607, -20.350, -12.798, -16.475, -8.039, -20.218, -12.897, -16.338, -9.773, -23.458, -4.808, 0.665, 0.109, -5.483]; 
        dH0s = cat(2,dH0s,ddH1,ddH2, ddH3);
        dS0s = cat(2,dS0s,ddS1,ddS2, ddS3);
        % Initiation and symmetry
        [dH0, dS0] = DNA_initiation(seq,startpos,endpos);
    end

    %Overhangs
    % 5'overhangs
    oseqs5 =   {'AA', 'AC',	   'AG',   'AT',  'CA',  'CC',   'CG',   'CT',	 'GA',   'GC',    'GG',    'GT',   'TA',  'TC',   'TG',   'TT'};
    o5dH0s =   [0.2,  -6.3,    -3.7,  -2.9,   0.6,   -4.4,   -4.0,   -4.1,   -1.1,   -5.1,    -3.9,    -4.2,   -6.9,  -4.0,   -4.9,   -0.2]; % in kcal/mol
    o5dS0s =   [2.3,  -17.1,   -10.0, -7.6,   3.3,   -12.6,  -11.9,  -13.0,  -1.6,   -14.0,  -10.9,  -15.0,   -20.0,  -10.9,  -13.8,  -0.5]; % in cal/mol/K
    % 3'overhangs
    oseqs3 =   {'TA',   'GA',   'CA',   'AA',   'TC',   'GC',   'CC',   'AC',   'TG',   'GG',   'CG',   'AG',   'TT',   'GT',   'CT',   'AT'}; % 3'overhangs (listed in 5'-3' order for the strand with the overhang)
    o3dH0s =   [-0.7,   -2.1,   -5.9,   -0.5,    4.4,   -0.2,   -2.6,   4.7,    -1.6,   -3.9,   -3.2,   -4.1,   2.9,    -4.4,   -5.2,   -3.8]; % in kcal/mol
    o3dS0s =   [-0.8,   -3.9,   -16.5,  -1.1,   14.9,   -0.1,   -7.4,   14.2,   -3.6,   -11.2,  -10.4,  -13.1,  10.4,   -13.1,  -15.0,  -12.6]; % in cal/mol/K

    % Iterate through all nearest neighbor sequences and add to deltaH and deltaS
    [seq2, mods] = parse_modifications(seq);
    for n = startpos:endpos-1
        nn = add_modifications(seq2(n:n+1),mods(n:n+1));
        dH0 = dH0 + dH0s(strcmpi(seqs,nn));
        dS0 = dS0 + dS0s(strcmpi(seqs,nn));
    end

    dS0
    % Apply dangling end (overhang) corrections
    if ~isempty(overhangs.seq{1})
        dH0 = dH0 + o5dH0s(strcmpi(oseqs5,overhangs.seq{1}));
        dS0 = dS0 + o5dS0s(strcmpi(oseqs5,overhangs.seq{1}));
    end
    if ~isempty(overhangs.seq{2})
        dH0 = dH0 + o3dH0s(strcmpi(oseqs3,overhangs.seq{2}));
        dS0 = dS0 + o3dS0s(strcmpi(oseqs3,overhangs.seq{2}));
    end
    if ~isempty(overhangs.target{1})
        dH0 = dH0 + o5dH0s(strcmpi(oseqs5,overhangs.target{1}));
        dS0 = dS0 + o5dS0s(strcmpi(oseqs5,overhangs.target{1}));
    end
    if ~isempty(overhangs.target{2})
        dH0 = dH0 + o3dH0s(strcmpi(oseqs3,overhangs.target{2}));
        dS0 = dS0 + o3dS0s(strcmpi(oseqs3,overhangs.target{2}));
    end
    dS0

    dH0 = dH0*1000; % Convert to cal/mol

    dG = dH0 - (273.15+T)*dS0;
end

function [dH0, dS0] = DNA_initiation(seq,startpos,endpos)
    seq2 = parse_modifications(seq);
    if strcmpi(seq2(startpos),'G') || strcmpi(seq2(startpos),'C')
        dH0 = 0.1;
        dS0 = -2.8;
    else
        dH0 = 2.3;
        dS0 = 4.1;
    end
    if strcmpi(seq2(endpos),'G') || strcmpi(seq2(endpos),'C')
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

function [type, T, mismatches, overhangs, startpos, endpos] = parse_input(seq,args)

    type = 'DNA'; % default to DNA
    T = 37; % default temp 37 C
    
    target = '';
    
    for n = 1:2:length(args)
        if strcmpi(args{n}, 'type')
            type = args{n+1};
        elseif strcmpi(args{n}, 'temperature')
            T = args{n+1};
        elseif strcmpi(args{n}, 'target')
            target = args{n+1};
            % may want to check that target sequence matches type
        end
    end
    
    if isempty(target)
        mismatches = {};
        overhangs.seq = {'', ''};
        overhangs.target = {'', ''};
        startpos = 1;
        endpos = length(parse_modifications(seq));
    else
        [mismatches, overhangs, startpos, endpos] = find_complementarity(seq, target);
    end
end