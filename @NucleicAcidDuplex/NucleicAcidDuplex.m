classdef NucleicAcidDuplex
    % Add mismatch treatment - only create paired Nearest Neighbor code if n and n+1 of
    % PairingState are 'p'; add codes for terminal mismatches
    properties
        Sequences = cell(2,1); % Cell array of strings or char with two (full) strand sequences
        SequenceNames = cell(2,1); % Cell array of strings or char with two strand names
        Schema = cell(2,0); % 2xN cell array showing register of two sequences in interaction
        PairingState = {}; % 1xN cell array showing pairing state ('p'=paired, 'w'=wobble,''=mismatch,'d'= dangling/overhang)
        NearestNeighbors = {};
        Nbp = []; % Number of base pairs in interaction
        Length = []; % Length of interaction, including mismatches and overhangs
        Tm = -Inf; % Melting temperature
        dS0 = -Inf; % Entropy of hybridization at standard conditions
        dH0 = Inf; % Enthalpy of hybridization at standard conditions
        dG0 = Inf; % Free energy of hybridization at standard conditions
        ParametersFile = "NN_Parameters.csv";
    end
    methods
        function obj = NucleicAcidDuplex(schema,varargin) % Constructor
            obj.Schema = schema;
            if ~isempty(varargin) % Parse initialization arguments
                for n = 1:2:length(varargin)
                    if strcmpi(varargin{n},'DuplexName')
                        obj.DuplexName = varargin{n+1};
                    elseif strcmpi(varargin{n},'Sequences')
                        obj.Sequences = varargin{n+1};
                    elseif strcmpi(varargin{n},'SequenceNames')
                        obj.SequenceNames = varargin{n+1};
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
            parameters = readtable(obj.ParametersFile); % Load nearest neighbor parameters
            obj = obj.findNearestNeighbors(); % Find nearest neighbors
            obj = determineSymmetryAndInitiation(obj); % Determine symmetry and initiation factors and add them to obj.NearestNeighbor property
            obj = obj.estimateThermodynamics(parameters); % Estimate deltaS, deltaH, deltaG
            obj = obj.estimateTm();
            obj.Length = length(obj.PairingState); % 
        end
        function obj = determinePairingState(obj)
            obj.PairingState = cell(1,size(obj.Schema,2));
            for n = 1:length(obj.PairingState)
                if contains(obj.Schema{1,n},'A') || contains(obj.Schema{1,n},'a')
                    if contains(obj.Schema{2,n},'T') || contains(obj.Schema{2,n},'t') || contains(obj.Schema{2,n},'U') || contains(obj.Schema{2,n},'u')
                        obj.PairingState{n} = 'p';
                    elseif isempty(obj.Schema{2,n})
                        obj.PairingState{n} = 'd';
                    else
                        obj.PairingState{n} = '-';
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
                        obj.PairingState{n} = 'w';
                    elseif isempty(obj.Schema{2,n})
                        obj.PairingState{n} = 'd';
                    else
                        obj.PairingState{n} = '-';
                    end
                elseif contains(obj.Schema{1,n},'T') || contains(obj.Schema{1,n},'t') || contains(obj.Schema{1,n},'U') || contains(obj.Schema{1,n},'u')
                    if contains(obj.Schema{2,n},'A') || contains(obj.Schema{2,n},'a')
                        obj.PairingState{n} = 'p';
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
            nn = {};
            for n = 1:size(obj.Schema,2)-1
                subcell = obj.Schema(1:2,n:n+1);
                char1 = '';
                if n == 1 && contains(obj.PairingState{n},{'-','d'}) % Mark terminal mismatches and overhangs as 3' or 5'
                    if isempty(subcell{1,1})
                        subcell = rot90(subcell,2); % If overhang is on bottom strand, flip overhang for lookup
                        char1 = strcat(char1,'3t');
                    else
                        char1 = strcat(char1,'5t');
                    end
                elseif n == size(obj.Schema,2)-1 && contains(obj.PairingState{n+1},{'-','d'}) % Mark terminal mismatches and overhangs as 3' or 5'
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
            obj.NearestNeighbors = nn;
        end
        function obj = determineSymmetryAndInitiation(obj)
            % Determine type
            rBases = sum(contains(obj.Schema,'r'),2);
            rBases = rBases>0;
            rBases = sum(rBases);
            if rBases==0
                % DNA/DNA
                ind = find(contains(obj.PairingState,'p'));
                ind = [min(ind), max(ind)];
                ind = contains(obj.Schema(1,ind),'G','IgnoreCase',true) | contains(obj.Schema(1,ind),'C','IgnoreCase',true);
                nGCterm = sum(ind);
                nATterm = 2 - nGCterm;
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
            elseif rBases==1
                % RNA/DNA
                obj.NearestNeighbors = ['DNARNAinit', obj.NearestNeighbors];
            elseif rBases==2
                % RNA/RNA
                disp('Warning: Initiation and symmetry parameters for RNA/RNA duplexes not yet in database. Initiation and symmetry will be ignored.')
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
                if sum(ind) > 0
                dH0 = dH0 + dH0s(ind);
                dS0 = dS0 + dS0s(ind);
                else
                    fprintf(1,'Error: code %s was not found in the lookup table of Parameters.  Ignoring this motif.\n',obj.NearestNeighbors{n});
                end
            end
            obj.dH0 = dH0*1000; % Convert to cal/mol before setting enthapy property
            obj.dS0 = dS0; % Set entropy property
            % Calculate deltaG
            obj.dG0 = obj.dH0 - 310.15*obj.dS0;
        end
        function obj = estimateTm(obj, varargin)
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
                    fprintf(1,'Warning: NucleicAcidDuplex.estimateTm() did not recognize argument "%s". Ignored this argument and any that immediately follow.  Please re-run function without this argument.', num2str(varargin{n}));
                end
            end
            % Calculate Tm under standard conditions
            R = 1.987204258; % Gas constant, cal/(mol K)
            % Calculate Tm from entropy, enthalpy, and concentration
            obj.Tm = obj.dH0/(obj.dS0 + R*log(conc/4))-273.15;
            % Apply salt correction
            obj.Tm = obj.salt_correction(obj.Tm,Na,Mg);
        end
        function Tm = salt_correction(obj,Tm_1MNa, Na, Mg)
            % Convert from Celsius to Kelvin
            Tm_1MNa = Tm_1MNa + 273.15;
            fGC = obj.gcContent();
            R = sqrt(Mg)/Na;
            if R < 0.22
                Tm = 1/(1/Tm_1MNa + ((4.29*fGC - 3.95)*log(Na) + 0.940*(log(Na))^2)*(1E-5)); % Monovalent correction
            elseif R < 6
                a = 3.92*(0.843-0.352*sqrt(Na)*log(Na));
                d = 1.42*(1.279-4.03E-3*log(Na)-8.03E-3*(log(Na))^2);
                g = 8.31*(0.486-0.258*log(Na)+5.25E-3*(log(Na))^3);
                Tm = 1/(1/Tm_1MNa + (a - 0.911*log(Mg)+fGC*(6.26+d*log(Mg))+1/(2*(obj.Nbp-1))*(-48.2+52.5*log(Mg)+g*(log(Mg))^2))*1E-5);
            else
                a = 3.92;
                d = 1.42;
                g = 8.31;
                Tm = 1/(1/Tm_1MNa + (a - 0.911*log(Mg)+fGC*(6.26+d*log(Mg)+1/(2*(obj.Nbp-1))*(-48.2+52.5*log(Mg)+g*log(Mg))))*1E-5);
            end
            % Convert from Kelvin back to Celsius
            Tm = Tm-273.15;
        end
        function dG = estimateDeltaG(obj,varargin)
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
                        fprintf(1,'Warning: NucleicAcidDuplex.estimateDeltaG() did not recognize argument "%s". Ignored this argument and any that immediately follow.  Please re-run function without this argument.', num2str(varargin{n}));
                    end
                end
            end
            T = T+273.15;
            dG = obj.dH0 - T*obj.dS0 - R*T*log(c/4);
        end
        function fGC = gcContent(obj)
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
            fGC = nGC/obj.Nbp;
        end
    end
end