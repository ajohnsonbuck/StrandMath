classdef NucleicAcidPair
    properties
        Sequences = cell(2,1); % Cell array of nucleic acid sequence objects
        Duplexes = {}; % Cell array of NucleicAcidDuplex objects describing one or more duplexes formed by the pair. By default, the first is the longest.
    end
    methods
        function obj = NucleicAcidPair(varargin) % Constructor
            if length(varargin) == 1
               if isa(varargin{1},'NucleicAcid')
                    obj.Sequences{1} = varargin{1};
                    obj.Sequences{2} = NucleicAcid(obj.Sequences{1}.reverseComplement('string')); % Create reverse complement
               else
                   disp('Error: Input must be one or two NucleicAcid objects');
               end
            elseif length(varargin) == 2
               if isa(varargin{1},'NucleicAcid') && isa(varargin{2},'NucleicAcid')
                   obj.Sequences{1} = varargin{1};
                   obj.Sequences{2} = varargin{2};
               else
                   disp('Error: Input must be one or two NucleicAcid objects');
               end
            end
            if sum(contains(obj.Sequences{2}.Sequence,'r'))<sum(contains(obj.Sequences{1}.Sequence,'r'))
                obj.Sequences = flipud(obj.Sequences); % If either sequence has RNA, ensure Sequence 2 has more RNA residues
            end
            if sum(contains(obj.Sequences{2}.Sequence,'+'))>sum(contains(obj.Sequences{2}.Sequence,'+'))
                obj.Sequences = flipud(obj.Sequences); % If either sequence has LNA, ensure Sequence 1 has more LNA residues
            end
            if sum(contains(obj.Sequences{2}.Sequence,'b'))>sum(contains(obj.Sequences{2}.Sequence,'b'))
                obj.Sequences = flipud(obj.Sequences); % If either sequence has BNA, ensure Sequence 1 has more BNA residues
            end
            obj = findLongestDuplex(obj);
        end
        function obj = findLongestDuplex(obj) % Find duplex with largest number of base pairs
            % Create schema with padding (empty cells) for all possible registers
            schema = cell(2,obj.Sequences{2}.len + (obj.Sequences{1}.len-1)*2);
            schema(2,obj.Sequences{1}.len:obj.Sequences{1}.len+obj.Sequences{2}.len-1) = obj.Sequences{2}.toDNA.reverseComplement.BareSequence; % Reverse complement of bare DNA version of first sequence
            seq1 = obj.Sequences{1}.toDNA.BareSequence; % first sequence to be slid across second sequence and compared
            ncomp_best = 0; % highest number of complementary base pairs
            comp_best = zeros(1,size(schema,2)); % matrix of complementary base pairs
            % Determine register of schema with most base pairs
            for n=1:size(schema,2)-obj.Sequences{1}.len+1
                schema(1,:) = cell(1,size(schema,2)); % empty first row
                schema(1,n:n+obj.Sequences{1}.len-1) = seq1;
                comp = cellfun(@strcmp,schema(1,:),schema(2,:));
                ncomp = sum(comp);
                if ncomp > ncomp_best
                    ncomp_best = ncomp;
                    comp_best = comp;
                    nbest = n;
                end
            end
            % Reconstruct schema with largest number of base pairs
            schema = cell(2,obj.Sequences{2}.len + (obj.Sequences{1}.len-1)*2);
            schema(2,obj.Sequences{1}.len:obj.Sequences{1}.len+obj.Sequences{2}.len-1) = obj.Sequences{2}.reverse().Sequence;
            schema(1,nbest:nbest+obj.Sequences{1}.len-1) = obj.Sequences{1}.Sequence;
            % Trim schema beyond first unpaired nucleotides at termini
            ind = movmax(comp_best,3); % Maximal region
            ind2 = any(~cellfun(@isempty,schema),1); % Region occupied by at least one nucleotide from maximal region
            ind = ind & ind2; % Trim region
            startpos = find(ind,1,'first');
            endpos = find(ind,1,'last');
            schema = schema(:, startpos:endpos); % trim
            schema(cellfun(@isempty,schema))={''}; % Replace empty cell elements with empty char
            % Create duplex object
            obj.Duplexes{1} = NucleicAcidDuplex(schema);
        end
        function duplex = longestDuplex(obj)
            duplex = obj.Duplexes{1};
        end
        function list(obj) % List nucleic acid sequences in pair as strings
            for n = 1:2
                fprintf(1,'Sequence %d: %s\n',n,obj.Sequences{n}.String);
            end
        end
    end
end