function Tm = estimate_Tm(seq,varargin)
    
    % Parse input arguments
    args = varargin;
    [type, conc, Na, Mg, target] = parse_input(seq,args);

    % Calculate thermodynamic parameters
    [dH0, dS0, ~] = parameters_NA(seq,'type',type,'concentration',conc,'Na',Na,'Mg',Mg,'temperature',37,'target',target);

    % Calculate Tm under standard conditions
    Tm = Tm_1MNa(dH0,dS0,conc);

    % Apply salt correction
    Tm = salt_correction(Tm,length(parse_modifications(seq)),gc_content(seq),Na,Mg);
end

function [type, conc, Na, Mg, target] = parse_input(seq,args)
    type = 'DNA';
    conc = 0.2E-6;
    Na = 1;
    Mg = 0;
    target = ''; % default is no target specified (fully complementary)
    for n = 1:2:length(args)
        if strcmpi(args{n}, 'type')
            type = args{n+1};
        elseif strcmpi(args{n}, 'concentration')
            conc = args{n+1};
        elseif strcmpi(args{n}, 'Na')
            Na = args{n+1};
        elseif strcmpi(args{n}, 'Mg')
            Mg = args{n+1};
        elseif strcmpi(args{n}, 'target')
            target = args{n+1};
        end
    end
end