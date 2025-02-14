function Tm = estimate_Tm(probe,varargin)
    % Estimate Tm of longest duplex between two provided sequences.
    % If only one argument is provided, it is assumed to be the probe, and
    % the target is assumed to be fully complementary.
    % ------------------------
    % First argument:
    % Probe sequence (type: NucleicAcid object, char, string, or cell array)
    % Name-value pair arguments:
    % 'Target'
    % Target sequence (type: NucleicAcid object, char, string, or cell array)
    % 'concentration' or 'conc'
    % Probe concentration, in mol/L (type: float)
    % 'Na'
    % Sodium concentration, in mol/L (type: float)
    % 'Mg'
    % Sodium concentration, in mol/L (type: float)

    % Parse input arguments
    if ~isa(probe,'NucleicAcid')
        probe = NucleicAcid(probe);
    end
    args = varargin;
    [type, conc, Na, Mg, target] = parse_input(probe, args);

    pair = NucleicAcidPair(probe,target); % Generate nucleic acid pair

    Tm = pair.longestDuplex.estimateTm('concentration',conc,'Na',Na,'Mg',Mg); % estimate Tm of longest duplex

    Tm = round(Tm,1);

end

function [type, conc, Na, Mg, target] = parse_input(probe, args)
    type = 'DNA';
    conc = 0.2E-6;
    Na = 1;
    Mg = 0;
    target = ''; % default is no target specified (fully complementary)
    for n = 1:2:length(args)
        if strcmpi(args{n}, 'type')
            type = args{n+1};
        elseif strcmpi(args{n}, 'concentration') || strcmpi(args{n}, 'conc')
            conc = args{n+1};
        elseif strcmpi(args{n}, 'Na')
            Na = args{n+1};
        elseif strcmpi(args{n}, 'Mg')
            Mg = args{n+1};
        elseif strcmpi(args{n}, 'target')
            target = args{n+1};
            if ~isa(target,'NucleicAcid')
                target = NucleicAcid(target);
            end
        else
            fprintf(1,'Warning: estimate_Tm() did not recognize argument "%s". Ignored this argument and any that immediately follow.  Please re-run function without this argument.', num2str(args{n}));
        end
    end
    if isempty(target)
        target = probe.reverseComplement();
    end
end