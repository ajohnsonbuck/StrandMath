classdef NucleicAcidPair
    properties
        Sequences = cell(2,1); % Cell array of nucleic acid sequence objects
        SequenceNames = cell(2,1); % Cell array of strings or char with two strand names
        Duplexes = {}; % Cell array of NucleicAcidDuplex objects describing one or more duplexes formed by the pair
    end
    methods
        function obj = NucleicAcidPair(varargin) % Constructor
            if length(varargin) == 1
               if isa(varargin{1},'NucleicAcid')
                    obj.Sequences{1} = varargin{1};
                    obj.Sequences{2} = NucleicAcid(obj.Sequences{1}.reverseComplement()); % Create reverse complement
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
        end
        function obj = findLongestDuplex(obj)
            % Assume fully complementary for now
            schema = [obj.Sequences{1}.Sequence;fliplr(obj.Sequences{2}.Sequence)];
            obj.Duplexes{1} = NucleicAcidDuplex(schema);
        end
    end
end