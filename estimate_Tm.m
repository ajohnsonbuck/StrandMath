function Tm = estimate_Tm(probe,varargin)
    
    % Parse input arguments
    if ~isa(probe,'NucleicAcid')
        probe = NucleicAcid(probe);
    end
    args = varargin;
    [type, conc, Na, Mg, target] = parse_input(probe, args);

    pair = NucleicAcidPair(probe,target); % Generate nucleic acid pair
    pair = pair.findLongestDuplex(); % Find longest duplex
    duplex = pair.Duplexes{1}.estimateTm('concentration',conc,'Na',Na,'Mg',Mg); % estimate Tm
    Tm = duplex.Tm;

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
        elseif strcmpi(args{n}, 'concentration')
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