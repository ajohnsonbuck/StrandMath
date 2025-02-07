classdef NucleicAcid
    properties
        String = '';
        Sequence = {};
        BareString = '';
        BareSequence = {};
        Modifications = {};
    end
    properties (SetAccess = private)
        Modlist = {'+','b','r'};
    end
    methods
        function obj = NucleicAcid(varargin)
            if nargin == 1
                str1 = varargin{1};
                obj = fromString(obj, str1);
            else
                for n = 1:length(varargin)
                    if strcmpi(varargin{n},'String')
                       str1 = varargin{n+1};
                       obj = fromString(obj, str1);
                    elseif strcmpi(varargin{n},'Sequence')
                       obj.Sequence = varargin{n+1};
                       obj = fromSequence(obj);
                    end
                end
            end
        end
        function obj = fromSequence(obj) % Populate object from input cell array of nucleotides
            obj.String = '';
            for n = 1:length(obj.Sequence)
                obj.String = strcat(obj.String,obj.Sequence{n});
            end
            seq = obj.String;
            obj.BareString = erase(seq,obj.Modlist); % Remove modification prefixes from sequences
            obj.BareSequence = cell(1,length(obj.BareString));
            obj.Modifications = cell(1,length(obj.BareString));
            c = 1;
            for n = 1:length(seq)
                if ismember(seq(n),obj.Modlist)
                    obj.Modifications{c} = seq(n);
                else
                    c = c+1;
                end
            end
            for n = 1:length(obj.BareString)
                obj.BareSequence{n} = obj.BareString(n);
            end
        end
        function obj = fromString(obj,str1) % Populat object from input string
            obj.String = str1;
            seq = obj.String;
            obj.BareString = erase(seq,obj.Modlist); % Remove modification prefixes from sequences
            obj.Sequence = cell(1,length(obj.BareString));
            obj.BareSequence = cell(1,length(obj.BareString));
            for n = 1:length(obj.Sequence)
                obj.Sequence{n} = '';
            end
            obj.Modifications = cell(1,length(obj.BareString));
            c = 1;
            for n = 1:length(seq)
                if ismember(seq(n),obj.Modlist)
                    obj.Modifications{c} = seq(n);
                    obj.Sequence{c} = strcat(obj.Sequence{c},seq(n)); 
                else
                    obj.Sequence{c} = strcat(obj.Sequence{c},seq(n));
                    c = c+1;
                end
            end
            for n = 1:length(obj.BareString)
                obj.BareSequence{n} = obj.BareString(n);
            end
        end
        function obj = stripModifications(obj)
            
        end
        function obj = toDNA(obj)
            str1 = replace(obj.BareString, 'U', 'T');
            obj = fromString(obj,str1);
        end
        function obj = toRNA(obj)
            seq1 = replace(obj.BareSequence, 'T', 'U');
            for n = 1:length(seq1)
                seq1{n} = strcat('r',seq1{n});
            end
            obj.Sequence = seq1;
            obj = fromSequence(obj);
        end 
        function rc = reverseComplement(obj, varargin)
            type = 'DNA';
            if ~isempty(varargin)
                if strcmpi(varargin{1},'RNA')
                    type = 'RNA';
                end
            end
            rc = obj.BareString;
            for m = 1:length(rc)
                n = length(rc)-m+1;
                base = obj.BareString(n);
                if strcmpi(base,'C')
                    comp = 'G';
                elseif strcmpi(base,'G')
                    comp = 'C';
                elseif strcmpi(base,'T')
                    comp = 'A';
                elseif strcmpi(base, 'U')
                    comp = 'A';
                elseif strcmpi(base,'A')
                    if strcmpi(type,'RNA')
                        comp = 'U';
                    else
                        comp = 'T';
                    end
                end
                rc(m)=comp;
            end
        end
        function len = length(obj)
           len = length(obj.BareString);
        end
    end
end