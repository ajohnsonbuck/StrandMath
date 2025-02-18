classdef Multistrand
    properties
        Strands(2,1) = Strand(); % 2x1 Strand array
        Duplexes = {}; % Cell array of Duplex objects describing one or more duplexes formed by the pair. By default, the first is the longest.
    end
    methods
        function objArray = Multistrand(varargin) % Constructor
            args = varargin;
            if length(args) == 1
                if isa(args{1},'Strand')
                    objArray(1,numel(args{1})) = Multistrand();
                    for n = 1:numel(args{1})
                        objArray(n).Strands(1) = args{1}(n);
                    end
                elseif isa(args{1},'string') || isa(args{1},'char') || isa(args{1},'cell')
                    objArray(1).Strands(1) = Strand(args{1});
                else
                    error('Input must be one or two Strand objects, chars, strings, or sequences');
                end
                for n = 1:numel(objArray)
                    objArray(n).Strands(2) = objArray(n).Strands(1).reverseComplement; % Create reverse complement
                    if ~isempty(objArray(n).Strands(1).Name)
                        objArray(n).Strands(2).Name = [objArray(n).Strands(1).Name,'_reverseComplement'];
                    end
                end
            elseif length(args) == 2
                for n = 1:2
                    if isa(args{n},'Strand')
                        for p = 1:numel(args{1})
                            objArray(p).Strands(n) = args{n}(p);
                        end
                    elseif isa(args{n},'string') || isa(args{n},'char') || isa(args{n},'cell')
                        objArray(1).Strands(n) = Strand(args{n});
                    else
                        error('Input must be one or two Strand objects, chars, strings, or sequences');
                    end
                end
            end
            if numel(objArray) > 100
                wb1 = waitbar(0,'Creating multi-strand objects and calculating duplexes...');
            end
            for n = 1:numel(objArray)
                if numel(objArray) > 100 && mod(n,100)==0
                waitbar(n/numel(objArray),wb1,['Created multi-strand object ',num2str(n),' of ',num2str(numel(objArray))]);
                end
                if sum(contains(objArray(n).Strands(2).Sequence,'r'))<sum(contains(objArray(n).Strands(1).Sequence,'r'))
                    objArray(n).Strands = flipud(objArray(n).Strands); % If either sequence has RNA, ensure Sequence 2 has more RNA residues
                end
                if sum(contains(objArray(n).Strands(2).Sequence,'+'))>sum(contains(objArray(n).Strands(1).Sequence,'+'))
                    objArray(n).Strands = flipud(objArray(n).Strands); % If either sequence has LNA, ensure Sequence 1 has more LNA residues
                end
                if sum(contains(objArray(n).Strands(2).Sequence,'b'))>sum(contains(objArray(n).Strands(1).Sequence,'b'))
                    objArray(n).Strands = flipud(objArray(n).Strands); % If either sequence has BNA, ensure Sequence 1 has more BNA residues
                end
                if ~isempty(objArray(n).Strands(1).String)
                    objArray(n) = findLongestDuplex(objArray(n));
                end
            end
        end
        function a = findLongestDuplex(a) % Find duplex with largest number of base pairs
            objArray = a;
            for m = 1:numel(objArray)
                objArray(m) = applyMask(objArray(m));
                % Create schema with padding (empty cells) for all possible registers
                schema = cell(2,objArray(m).Strands(2).len + (objArray(m).Strands(1).len-1)*2);
                encodedSchema = ones(size(schema)); % Initialize with 1 (code for empty position)
                encodedSchema(2,objArray(m).Strands(1).len:objArray(m).Strands(1).len+objArray(m).Strands(2).len-1) = Multistrand.encodeSequence(objArray(m).Strands(2).reverse.bareSequence); % Encoded bare version of first sequence
                seq1 = Multistrand.encodeSequence(objArray(m).Strands(1).bareSequence); % encoded first sequence to be slid across second sequence and compared
                nbest = objArray(m).Strands(1).len;
                score_best = 0; % highest complementarity score
                % Determine register of schema with most base pairs
                for n=1:size(schema,2)-objArray(m).Strands(1).len+1
                    encodedSchema(1,:) = 1; % Empty first row
                    encodedSchema(1,n:n+objArray(m).Strands(1).len-1) = seq1; % place encoded Sequence{1} into first row of encodedSchema at position n
                    score = Multistrand.scoreBasePairs(encodedSchema); % score base pairs of encodedSchema
                    if score > score_best
                        score_best = score;
                        nbest = n;
                    end
                end
                % Reconstruct schema with largest number of base pairs
                schema = cell(2,objArray(m).Strands(2).len + (objArray(m).Strands(1).len-1)*2);
                schema(2,objArray(m).Strands(1).len:objArray(m).Strands(1).len+objArray(m).Strands(2).len-1) = objArray(m).Strands(2).reverse().Sequence;
                schema(1,nbest:nbest+objArray(m).Strands(1).len-1) = objArray(m).Strands(1).Sequence;
                % Trim schema of any padding
                ind = any(~cellfun(@isempty,schema),1);
                startpos = find(ind,1,'first');
                endpos = find(ind,1,'last');
                schema = schema(:, startpos:endpos); % trim
                schema(cellfun(@isempty,schema))={''}; % Replace empty cell elements with empty char
                % Create duplex object and place in original Multistrand array
                a(m).Duplexes{1} = Duplex(schema,'Strands',objArray(m).Strands);
            end
        end
        function duplex = longestDuplex(objArray)
            for n = 1:numel(objArray) 
                duplex(n) = objArray(n).Duplexes{1};
            end
        end
        function list(obj) % List nucleic acid sequences in pair as strings
            for n = 1:2
                fprintf(1,'Sequence %d: %s\n',n,obj.Strands(n).String);
            end
        end
        function Tm = estimateTm(objArray,varargin)
            args = varargin;
            Tm = zeros(numel(objArray),1);
            for n = 1:numel(objArray)
                duplex = objArray(n).longestDuplex();
                if ~isempty(varargin)
                    Tm(n) = duplex.estimateTm(args{:});
                else
                    Tm(n) = duplex.estimateTm();
                end
            end
        end
        function objArray = applyMask(objArray)
            for m = 1:numel(objArray)
                for n = 1:numel(objArray(m).Strands)
                    mask = objArray(m).Strands(n).Mask;
                    if isempty(mask)
                        mask = repmat('n',1,objArray(m).Strands(n).len);
                    end
                    str1 = objArray(m).Strands(n).String;
                    for p = 1:objArray(m).Strands(n).len
                        if strcmp(mask(p),'-')
                            objArray(m).Strands(n).Sequence{p}='-';
                        end
                    end
                    objArray(m).Strands(n) = objArray(m).Strands(n).fromSequence;
                    objArray(m).Strands(n).UnmaskedString = str1;
                end
            end
        end
        function print(objArray, varargin)
            if numel(varargin)==0
                objArray.longestDuplex.print;
            else
                if strcmpi(varargin{1}, 'longestDuplex')
                    if numel(varargin)==2
                        if strcmpi(varargin{2},'flip')
                            objArray.longestDuplex.print('flip');
                        else
                            error(['Unknown second argument ', num2string(varargin{2}), ' passed to Multistrand.print.']);
                        end
                    else
                        objArray.longestDuplex.print;
                    end
                elseif strcmpi(varargin{1},'flip')
                    A = objArray.longestDuplex;
                    A.print('flip');
                elseif strcmpi(varargin{1},'strands')
                    for m = 1:numel(objArray)
                        for n = 1:numel(objArray(m).Strands)
                            fprintf(1,'\n Sequence %d: %s',n,objArray(m).Strands(n).Name)
                            fprintf(1,[char("\n5'-"),objArray(m).Strands(n).String,char("-3'\n")]);
                        end
                    end
                    fprintf(1,'\n');
                else 
                    error("Unknown argument passed to Multistrand.print; allowed arguments are 'strands' or 'longestDuplex'");
                end
            end
        end
    end
    methods (Static)
        function score = scoreBasePairs(encodedSchema) 
            persistent scoreMat
            if isempty(scoreMat)
                scoreMat = [0,	0,	0,	0,	0,	0;... % Rows and cols are: (empty), A, C, G, T, U; Score: G-C = 6, A-U/T = 4, G-U = 3, G-T = 2
                            0,	0,	0,	0,	4,	4;...
                            0,	0,	0,	6,	0,	0;...
                            0,	0,	6,	0,	2,	3;...
                            0,	4,	0,	2,	0,	0;...
                            0,	4,	0,	3,	0,	0];
                scoreMat = int8(scoreMat);
            end
            score = zeros(width(encodedSchema),1);
            for n = 1:width(encodedSchema)
                score(n,1) = scoreMat(encodedSchema(1,n),encodedSchema(2,n));
            end
            score = sum(score);
        end
        function encoded = encodeSequence(seq)
            % Initialize array; any positions not occupied by a nucleobase
            % are 1
            encoded = ones(size(seq));
            % Encode bases as numbers
            encoded(strcmp(seq,'A'))=2;
            encoded(strcmp(seq,'C'))=3;
            encoded(strcmp(seq,'G'))=4;
            encoded(strcmp(seq,'T'))=5;
            encoded(strcmp(seq,'U'))=6;
        end
    end
end