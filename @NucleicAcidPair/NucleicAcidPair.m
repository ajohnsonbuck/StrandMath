classdef NucleicAcidPair
    properties
        Sequences = {NucleicAcid(); NucleicAcid()}; % Cell array of two NucleicAcid objects
        Duplexes = {}; % Cell array of NucleicAcidDuplex objects describing one or more duplexes formed by the pair. By default, the first is the longest.
    end
    methods
        function obj = NucleicAcidPair(varargin) % Constructor
            if length(varargin) == 1
                if isa(varargin{1},'NucleicAcid')
                    obj.Sequences{1} = varargin{1};
                elseif isa(varargin{1},'string') || isa(varargin{1},'char') || isa(varargin{1},'cell')
                    obj.Sequences{1} = NucleicAcid(varargin{1});
                else
                    disp('Error: Input must be one or two NucleicAcid objects, chars, strings, or sequences');
                end
                obj.Sequences{2} = obj.Sequences{1}.reverseComplement; % Create reverse complement
            elseif length(varargin) == 2
                for n = 1:2
                    if isa(varargin{n},'NucleicAcid')
                        obj.Sequences{n} = varargin{n};
                    elseif isa(varargin{n},'string') || isa(varargin{n},'char') || isa(varargin{n},'cell')
                        obj.Sequences{n} = NucleicAcid(varargin{n});
                    else
                        disp('Error: Input must be one or two NucleicAcid objects, chars, strings, or sequences');
                    end
                end
            end
            if sum(contains(obj.Sequences{2}.Sequence,'r'))<sum(contains(obj.Sequences{1}.Sequence,'r'))
                obj.Sequences = flipud(obj.Sequences); % If either sequence has RNA, ensure Sequence 2 has more RNA residues
            end
            if sum(contains(obj.Sequences{2}.Sequence,'+'))>sum(contains(obj.Sequences{1}.Sequence,'+'))
                obj.Sequences = flipud(obj.Sequences); % If either sequence has LNA, ensure Sequence 1 has more LNA residues
            end
            if sum(contains(obj.Sequences{2}.Sequence,'b'))>sum(contains(obj.Sequences{1}.Sequence,'b'))
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
            % Trim schema of any padding
            ind = any(~cellfun(@isempty,schema),1);
            startpos = find(ind,1,'first');
            endpos = find(ind,1,'last');
            schema = schema(:, startpos:endpos); % trim
            schema(cellfun(@isempty,schema))={''}; % Replace empty cell elements with empty char
            % Create duplex object
            obj.Duplexes{1} = NucleicAcidDuplex(schema,'Sequences',obj.Sequences);
        end
        function duplex = longestDuplex(obj)
            duplex = obj.Duplexes{1};
        end
        function list(obj) % List nucleic acid sequences in pair as strings
            for n = 1:2
                fprintf(1,'Sequence %d: %s\n',n,obj.Sequences{n}.String);
            end
        end
        function Tm = estimateTm(obj,varargin)
            args = varargin;
            duplex = obj.longestDuplex();
            if ~isempty(varargin)
                Tm = duplex.estimateTm(args{:});
            else
                Tm = duplex.estimateTm();
            end
        end
        function print(objArray)
            for m = 1:numel(objArray)
                for n = 1:numel(objArray(m).Sequences)
                    fprintf(1,'\n Sequence %d: %s',n,objArray(m).Sequences{n}.Name)
                    fprintf(1,['\n5-',objArray(m).Sequences{n}.String,'-3\n']);
                end
            end
            fprintf(1,'\n');
        end
    end
end