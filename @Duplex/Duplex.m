classdef Duplex
    % Add mismatch treatment - only create paired Nearest Neighbor code if n and n+1 of
    % PairingState are 'p'; add codes for terminal mismatches
    properties
        Sequences = {Strand(); Strand()}; % Cell array of two Strand objects
        Schema = cell(2,0); % 2xN cell array showing register of two sequences in interaction
        PairingState = {}; % 1xN cell array showing pairing state ('p'=paired, 'w'=wobble,''=mismatch,'d'= dangling/overhang)
        NearestNeighbors = {};
        Nbp = []; % Number of base pairs in interaction
        Length = []; % Length of interaction, including mismatches and overhangs
        fGC = -Inf;
        dS0 = -Inf; % Entropy of hybridization at standard conditions
        dH0 = Inf; % Enthalpy of hybridization at standard conditions
        dG0 = Inf; % Free energy of hybridization at standard conditions
    end
    properties (Constant)
        parameters = readtable("NN_Parameters.csv"); % Load nearest neighbor parameters;
    end
    methods
        function obj = Duplex(schema,varargin) % Constructor
            obj.Schema = schema;
            if ~isempty(varargin) % Parse initialization arguments
                for n = 1:2:length(varargin)
                    if strcmpi(varargin{n},'Sequences')
                        obj.Sequences = varargin{n+1};
                    elseif strcmpi(varargin{n},'PairingState')
                        obj.PairingState = varargin{n+1};
                    elseif strcmpi(varargin{n}, 'Nbp')
                        obj.Nbp = varargin{n+1};
                    end
                end
            end
            if isempty(obj.PairingState) % Determine pairing state if not provided
                obj = determinePairingState(obj);
            end
            if isempty(obj.Nbp) % Determine number of base pairs if not provided
                nOverhangs = sum(strcmp(obj.PairingState,'d'));
                nMismatches = sum(strcmp(obj.PairingState,'-'));
                obj.Nbp = length(obj.PairingState)-nOverhangs-nMismatches;
            end
            % parameters = readtable(obj.ParametersFile); % Load nearest neighbor parameters
            obj = obj.findNearestNeighbors(); % Find nearest neighbors
            obj = determineSymmetryAndInitiation(obj); % Determine symmetry and initiation factors and add them to obj.NearestNeighbors property
            obj = obj.estimateThermodynamics(Duplex.parameters); % Estimate deltaS, deltaH, deltaG
            obj = obj.gcContent();
            % obj = obj.estimateTm(); % Don't estimate by default - user may have specific conditions to request for Tm estimation
            obj.Length = length(obj.PairingState); % 
        end
        function obj = determinePairingState(obj)
            obj.PairingState = cell(1,size(obj.Schema,2));
            for n = 1:length(obj.PairingState)
                if contains(obj.Schema{1,n},'A') || contains(obj.Schema{1,n},'a')
                    if contains(obj.Schema{2,n},'T') || contains(obj.Schema{2,n},'t') || contains(obj.Schema{2,n},'U') || contains(obj.Schema{2,n},'u')
                        obj.PairingState{n} = 'p'; % paired
                    elseif isempty(obj.Schema{2,n})
                        obj.PairingState{n} = 'd'; % dangling end/overhang
                    else
                        obj.PairingState{n} = '-'; % unpaired
                    end
                elseif contains(obj.Schema{1,n},'C') || contains(obj.Schema{1,n},'c')
                    if contains(obj.Schema{2,n},'G') || contains(obj.Schema{2,n},'g')
                        obj.PairingState{n} = 'p';
                    elseif isempty(obj.Schema{2,n})
                        obj.PairingState{n} = 'd';
                    else
                        obj.PairingState{n} = '-';
                    end
                elseif contains(obj.Schema{1,n},'G') || contains(obj.Schema{1,n},'g')
                    if contains(obj.Schema{2,n},'C') || contains(obj.Schema{2,n},'c')
                        obj.PairingState{n} = 'p';
                    elseif contains(obj.Schema{2,n},'U') || contains(obj.Schema{2,n},'u')
                        obj.PairingState{n} = 'w'; % wobble
                    elseif isempty(obj.Schema{2,n})
                        obj.PairingState{n} = 'd';
                    else
                        obj.PairingState{n} = '-';
                    end
                elseif contains(obj.Schema{1,n},'T') || contains(obj.Schema{1,n},'t') || contains(obj.Schema{1,n},'U') || contains(obj.Schema{1,n},'u')
                    if contains(obj.Schema{2,n},'A') || contains(obj.Schema{2,n},'a')
                        obj.PairingState{n} = 'p';
                    elseif contains(obj.Schema{2,n},'G') || contains(obj.Schema{2,n},'g')
                        obj.PairingState{n} = 'w';
                    elseif isempty(obj.Schema{2,n})
                        obj.PairingState{n} = 'd';
                    else
                        obj.PairingState{n} = '-';
                    end
                elseif isempty(obj.Schema{1,n})
                    if ~isempty(obj.Schema{2,n})
                        obj.PairingState{n} = 'd';
                    end
                end
            end
        end
        function obj = findNearestNeighbors(obj)
            [schema,pairingState] = trimSchema(obj);
            nn = {};
            for n = 1:size(schema,2)-1
                subcell = schema(1:2,n:n+1);
                char1 = '';
                if n == 1 && contains(pairingState{n},{'-','d'}) % Mark terminal mismatches and overhangs as 3' or 5'
                    if isempty(subcell{1,1})
                        subcell = rot90(subcell,2); % If overhang is on bottom strand, flip overhang for lookup
                        char1 = strcat(char1,'3t');
                    else
                        char1 = strcat(char1,'5t');
                    end
                elseif n == size(schema,2)-1 && contains(pairingState{n+1},{'-','d'}) % Mark terminal mismatches and overhangs as 3' or 5'
                    if isempty(subcell{1,2})
                        subcell = rot90(subcell,2); % If overhang is on bottom strand, flip overhang for lookup
                        char1 = strcat(char1,'5t');
                    else
                        char1 = strcat(char1,'3t');
                    end
                end
                char1 = strcat(char1, subcell{1,1}, subcell{1,2},'/', subcell{2,1}, subcell{2,2});
                nn = [nn, char1];
            end
            % Handle special cases with more than 2 nearest neighbors
            % Internal RNA mismatches and loops
            for n = 1:numel(nn)-1
                if contains(schema{1,n+1},{'r','+','b'}) && contains(schema{2,n+1},{'r','+','b'}) && strcmp(pairingState{n+1},'-') && strcmp(pairingState{n},'p') % if RNA internal mismatch
                    p = 1;
                    while strcmp(pairingState{n+p},'-')
                        p = p+1;
                    end
                    nn{n} = strcat(erase(char(join(schema(1,n:n+p))),' '),'/',erase(char(join(schema(2,n:n+p))),' '));
                    for q = 1:p-1
                        nn{n+q}=''; % eliminate duplicate nn entries covered by the larger entry
                    end
                end
            end
            nn(cellfun(@isempty,nn))=[]; % Clear any empty cells of nn
            obj.NearestNeighbors = nn;
        end
        function [schema,pairingState] = trimSchema(obj)
            % Trim schema beyond first unpaired nucleotides at termini
            schema = obj.Schema;
            schema = replace(schema,'-',''); % Replace masked cell elements with empty char
            pairingState = obj.PairingState;
            comp = cellfun(@(x) any(ismember(x, {'p','w'})), pairingState); % identify 
            ind = movmax(comp,3); % Maximal region
            ind2 = any(~cellfun(@isempty,schema),1); % Region occupied by at least one nucleotide from maximal region
            ind = ind & ind2; % Trim region
            startpos = find(ind,1,'first');
            endpos = find(ind,1,'last');
            schema = schema(:, startpos:endpos); % trim
            pairingState = pairingState(startpos:endpos);
        end
        function obj = determineSymmetryAndInitiation(obj)
            % Determine type
            rBases = sum(contains(obj.Schema,'r'),2);
            rBases = rBases>0;
            rBases = sum(rBases);
            dBases = sum(~cellfun(@isempty,obj.Schema) & ~contains(obj.Schema,{'r','+','b'}),2);
            dBases = dBases>0;
            dBases = sum(dBases);
            if rBases==0
                % DNA/DNA
                ind = find(contains(obj.PairingState,'p'));
                ind = [min(ind), max(ind)];
                ind = contains(obj.Schema(1,ind),'G','IgnoreCase',true) | contains(obj.Schema(1,ind),'C','IgnoreCase',true);
                nGCterm = sum(ind);
                nATterm = 2 - nGCterm; % Note: this nATterm also includes GU wobbles, as these should be penalized the same as terminal AUs
                if nGCterm > 0
                    for n = 1:nGCterm
                        obj.NearestNeighbors = ['DNAinitGC', obj.NearestNeighbors];
                    end
                end
                if nATterm > 0
                    for n = 1:nATterm
                        obj.NearestNeighbors = ['DNAinitAT', obj.NearestNeighbors];
                    end
                end
                % Check for and apply symmetry
                ind = find(contains(obj.PairingState,'p'));
                if sum(~strcmp(obj.Schema(1,ind),fliplr(obj.Schema(2,ind))))==0
                    obj.NearestNeighbors = ['DNAsymmetry', obj.NearestNeighbors];
                end
            elseif rBases==1 && dBases==1
                % RNA/DNA
                obj.NearestNeighbors = ['DNARNAinit', obj.NearestNeighbors];
            else
                % RNA/RNA, LNA/RNA, BNA/RNA
                ind = find(contains(obj.PairingState,'p'));
                ind = [min(ind), max(ind)];
                ind = contains(obj.Schema(1,ind),'A','IgnoreCase',true) | contains(obj.Schema(1,ind),'U','IgnoreCase',true);
                nAUterm = sum(ind);
                obj.NearestNeighbors = ['RNAinit', obj.NearestNeighbors];
                if nAUterm > 0
                    for n = 1:nAUterm
                        obj.NearestNeighbors = ['RNAinitAU', obj.NearestNeighbors];
                    end
                end
                % Check for and apply symmetry
                ind = find(contains(obj.PairingState,'p'));
                if sum(~strcmp(obj.Schema(1,ind),fliplr(obj.Schema(2,ind))))==0
                    obj.NearestNeighbors = ['RNAsymmetry', obj.NearestNeighbors];
                end
            end
        end
        function obj = estimateThermodynamics(obj,parameters)
            codes = parameters.NearestNeighbor;
            dH0s = parameters.dH0;
            dS0s = parameters.dS0;
            % Apply parameters
            dH0 = 0;
            dS0 = 0;
            for n = 1:length(obj.NearestNeighbors)
                ind = strcmpi(obj.NearestNeighbors{n},codes);
                if sum(ind) == 0
                    % Try approximating LNA and BNA as RNA
                    nn = replace(obj.NearestNeighbors{n},{'b','+'},'r');
                    nn = replace(nn,'T','U');
                    ind = strcmpi(nn,codes);
                    ind = strcmpi(nn,codes);
                    if sum(ind) == 0
                        % Try approximating RNA, LNA and BNA as DNA
                        nn = erase(obj.NearestNeighbors{n},{'r','b','+'});
                        ind = strcmpi(nn,codes);
                        if sum(ind) == 0
                            % Try approximating terminal mismatches as
                            % internal mismatches
                            nn = erase(obj.NearestNeighbors{n},{'5t','3t'});
                            ind = strcmpi(nn,codes);
                            if sum(ind) == 0
                                % Try approximating RNA as DNA
                                nn = erase(obj.NearestNeighbors{n},'r');
                                nn = replace(nn,'U','T');
                                if sum(ind) == 0
                                    fprintf(1,'Warning: code %s was not found in the lookup table of Parameters.  Ignoring this motif.\n',obj.NearestNeighbors{n});
                                else
                                    dH0 = dH0 + dH0s(ind);
                                    dS0 = dS0 + dS0s(ind);
                                    fprintf(1, 'Code %s was not found in the lookup table of Parameters.  Approximating as %s.\n', obj.NearestNeighbors{n}, char(codes(ind)));
                                end
                            else
                                dH0 = dH0 + dH0s(ind);
                                dS0 = dS0 + dS0s(ind);
                                fprintf(1, 'Code %s was not found in the lookup table of Parameters.  Approximating as %s.\n', obj.NearestNeighbors{n}, char(codes(ind)));
                            end
                        else
                            dH0 = dH0 + dH0s(ind);
                            dS0 = dS0 + dS0s(ind);
                            fprintf(1, 'Code %s was not found in the lookup table of Parameters.  Approximating as %s.\n', obj.NearestNeighbors{n}, char(codes(ind)));
                        end
                    else
                        dH0 = dH0 + dH0s(ind);
                        dS0 = dS0 + dS0s(ind);
                        fprintf(1, 'Code %s was not found in the lookup table of Parameters.  Approximating as %s.\n', obj.NearestNeighbors{n}, char(codes(ind)));
                    end
                else
                    dH0 = dH0 + dH0s(ind);
                    dS0 = dS0 + dS0s(ind);
                end
            end
            obj.dH0 = dH0*1000; % Convert to cal/mol before setting enthapy property
            obj.dS0 = dS0; % Set entropy property
            % Calculate deltaG
            obj.dG0 = obj.dH0 - 310.15*obj.dS0;
        end
        function Tm = estimateTm(objArray, varargin)
            Tm = zeros(numel(objArray),1);
            conc = 0.2E-6;
            Na = 1;
            Mg = 0;
            for n = 1:2:length(varargin)
                if strcmpi(varargin{n}, 'Na')
                    Na = varargin{n+1};
                elseif strcmpi(varargin{n}, 'Mg')
                    Mg = varargin{n+1};
                elseif strcmpi(varargin{n}, 'concentration') || strcmpi(varargin{n}, 'conc')
                    conc = varargin{n+1};
                else
                    fprintf(1,'Warning: Duplex.estimateTm() did not recognize argument "%s". Ignored this argument and any that immediately follow.  Please re-run function without this argument.', num2str(varargin{n}));
                end
            end
            R = 1.987204258; % Gas constant, cal/(mol K)
            for m = 1:numel(objArray)
                % Check for symmetry
                if sum(contains(objArray(m).NearestNeighbors,'symm'))>0
                    a = 1; % for self-complementary duplexes
                else
                    a = 4; % for non-self-complementary duplexes
                end
                % Calculate Tm from entropy, enthalpy, and concentration
                Tm(m) = objArray(m).dH0/(objArray(m).dS0 + R*log(conc/a))-273.15;
                % Apply salt correction
                Tm(m) = objArray(m).salt_correction(Tm(m),Na,Mg);
            end
        end
        function Tm = salt_correction(obj, Tm_1MNa, Na, Mg)
            % Convert from Celsius to Kelvin
            Tm_1MNa = Tm_1MNa + 273.15;
            R = sqrt(Mg)/Na;
            if R < 0.22
                Tm = 1/(1/Tm_1MNa + ((4.29*obj.fGC - 3.95)*log(Na) + 0.940*(log(Na))^2)*(1E-5)); % Monovalent correction
            elseif R < 6
                a = 3.92*(0.843-0.352*sqrt(Na)*log(Na));
                d = 1.42*(1.279-4.03E-3*log(Na)-8.03E-3*(log(Na))^2);
                g = 8.31*(0.486-0.258*log(Na)+5.25E-3*(log(Na))^3);
                Tm = 1/(1/Tm_1MNa + (a - 0.911*log(Mg)+obj.fGC*(6.26+d*log(Mg))+1/(2*(obj.Nbp-1))*(-48.2+52.5*log(Mg)+g*(log(Mg))^2))*1E-5);
            else
                a = 3.92;
                d = 1.42;
                g = 8.31;
                Tm = 1/(1/Tm_1MNa + (a - 0.911*log(Mg)+fGC*(6.26+d*log(Mg)+1/(2*(obj.Nbp-1))*(-48.2+52.5*log(Mg)+g*log(Mg))))*1E-5);
            end
            % Convert from Kelvin back to Celsius
            Tm = Tm-273.15;
        end
        function dG = estimateDeltaG(objArray,varargin)
            R = 1.987204258;% Gas constant, cal/(mol K)
            T = 37; % Temperature default 37 C
            c = 1; % Concentration default 1 M
            if ~isempty(varargin)
                for n = 1:2:length(varargin)
                    if strcmpi(varargin{n},'temperature') || strcmpi(varargin{n},'T')
                        T = varargin{n+1};
                    elseif strcmpi(varargin{n},'concentration') || strcmpi(varargin{n}, 'conc')
                        c = varargin{n+1};
                    else
                        fprintf(1,'Warning: Duplex.estimateDeltaG() did not recognize argument "%s". Ignored this argument and any that immediately follow.  Please re-run function without this argument.', num2str(varargin{n}));
                    end
                end
            end
            T = T+273.15;
            dG = zeros(size(objArray));
            for n = 1:numel(objArray)
                if sum(contains(objArray(n).NearestNeighbors,'symm'))>0
                    a = 1; % for self-complementary duplexes
                else
                    a = 4; % for non-self-complementary duplexes
                end
                dG(n) = objArray(n).dH0 - T*objArray(n).dS0 - R*T*log(c/a);
            end
        end
        function obj = gcContent(obj)
            % Calculate GC content of interaction
            nGC = 0;
            for n = 1:size(obj.Schema,2)
                if contains(obj.Schema{1,n},'G') || contains(obj.Schema{1,n},'g')
                    if contains(obj.Schema{2,n},'C') || contains(obj.Schema{2,n},'c')
                        nGC = nGC+1;
                    end
                elseif contains(obj.Schema{1,n},'C') || contains(obj.Schema{1,n},'c')
                    if contains(obj.Schema{2,n},'G') || contains(obj.Schema{2,n},'g')
                        nGC = nGC+1;
                    end
                end
            end
            obj.fGC = nGC/obj.Nbp;
        end
        function NN = nn(objArray)
            for n = 1:numel(objArray)
                NN{n} = objArray(n).NearestNeighbors;
            end
            if numel(objArray)==1
                NN = NN{:};
            end
        end
        function print(objArray,varargin)
            flipSequences=false; % put bottom sequence on top
            if numel(varargin)>0
                for n = 1:numel(varargin)
                    if strcmp(varargin{n},'flip')
                        flipSequences = true;
                    end
                end
            end
            for m = 1:numel(objArray)
                if flipSequences == true
                    objArray(m).Sequences = flipud(objArray(m).Sequences);
                    objArray(m).Schema = rot90(objArray(m).Schema,2);
                    objArray(m).PairingState = fliplr(objArray(m).PairingState);
                end
                % Find positions to place sequence names
                name1start = find(~cellfun(@isempty,objArray(m).Schema(1,:)),1,'first');
                name2start = find(~cellfun(@isempty,objArray(m).Schema(2,:)),1,'first');
                line0 = '\n ';
                line1 = char("\n5'-");
                line2 = '\n   ';
                line3 = char("\n3'-");
                line4 = '\n ';
                line5 = '\n dG0 = ';
                for n = 1:size(objArray(m).Schema,2)
                    if length(objArray(m).Schema{1,n})==2
                        line1 = [line1,objArray(m).Schema{1,n}];
                    elseif length(objArray(m).Schema{1,n})==1
                        line1 = [line1,' ',objArray(m).Schema{1,n}];
                    else
                        line1 = [line1,'  ',objArray(m).Schema{1,n}];
                    end
                    if length(objArray(m).Schema{2,n})==2
                        line3 = [line3,objArray(m).Schema{2,n}];
                    elseif length(objArray(m).Schema{2,n})==1
                        line3 = [line3,' ',objArray(m).Schema{2,n}];
                    else
                        line3 = [line3,'  ',objArray(m).Schema{2,n}];
                    end
                    if objArray(m).PairingState{n}=='p'
                        line2 = [line2,' |'];
                    elseif objArray(m).PairingState{n}=='w'
                        line2 = [line2,' o'];
                    else
                        line2 = [line2,'  '];
                    end
                end
                line1 = [line1,char("-3'")];
                line3 = [line3,char("-5'")];
                for n = 1:name1start
                    line0 = [line0, '  '];
                end
                for n = 1:name2start
                    line4 = [line4, '  '];
                end
                fprintf(1,[line0, '%s'],objArray(m).Sequences{1}.Name);
                fprintf(1,[line1,line2,line3],'\n');
                fprintf(1,[line4,'%s\n'],objArray(m).Sequences{2}.Name);
                fprintf(1,[line5,'%.1f kcal/mol\n'],objArray(m).dG0/1000);
            end
            fprintf(1,'\n')
        end
    end
end