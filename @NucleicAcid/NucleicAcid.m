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
                seq = varargin{1};
                if isa(seq,'char') || isa(seq,'string')
                    obj = fromString(obj, seq);
                else
                    obj = fromSequence(obj, seq);
                end
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
        function obj = fromSequence(obj,varargin) % Populate object from input cell array of nucleotides
            if ~isempty(varargin)
                obj.Sequence = varargin{1};
            end
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
        function obj = fromString(obj,str1) % Populate object from input string
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
        function objArray = toDNA(objArray)
            for n = 1:numel(objArray)
                str1 = replace(objArray(n).BareString, 'U', 'T');
                objArray(n) = fromString(objArray(n),str1);
            end
        end
        function objArray = toRNA(objArray)
            for m = 1:numel(objArray)
                seq1 = replace(objArray(m).BareSequence, 'T', 'U');
                for n = 1:length(seq1)
                    seq1{n} = strcat('r',seq1{n});
                end
                objArray(m).Sequence = seq1;
                objArray(m) = objArray(m).fromSequence();
            end
        end
        function r = reverse(objArray,varargin)
            outputType = 'NucleicAcid'; % Default: provide reverse complement as NucleicAcid unless otherwise provided as argument
            r = cell(size(objArray));
            if ~isempty(varargin)
                for n = 1:length(varargin)
                    if strcmpi(varargin{n},'char') || strcmpi(varargin{n},'string')
                        outputType = 'char';
                    elseif strcmpi(varargin{n},'sequence') || strcmpi(varargin{n},'cell')
                        outputType = 'sequence';
                    end
                end
            end
            if strcmp(outputType,'NucleicAcid')
                rNA = objArray; % copy object array initially
            end
            for n = 1:numel(objArray)
                seq = fliplr(objArray(n).Sequence);
                if strcmpi(outputType,'char')
                    r{n} = horzcat(seq{:});
                elseif strcmpi(outputType,'sequence')
                    r{n} = seq;
                elseif strcmp(outputType,'NucleicAcid')
                    rNA(n) = NucleicAcid(seq);
                end
            end
            if strcmpi(outputType,'NucleicAcid')
                r = rNA;
            end
            if isa(r,'cell') && numel(objArray)==1
                r = r{1};
            end
        end
        function rc = reverseComplement(objArray, varargin)
            type = 'DNA';
            outputType = 'NucleicAcid'; % Default: provide reverse complement as NucleicAcid unless otherwise provided as argument
            rc = cell(size(objArray));
            if ~isempty(varargin)
                for n = 1:length(varargin)
                    if strcmpi(varargin{n},'RNA')
                        type = 'RNA';
                    elseif strcmpi(varargin{n},'char') || strcmpi(varargin{n},'string')
                        outputType = 'char';
                    elseif strcmpi(varargin{n},'sequence') || strcmpi(varargin{n},'cell')
                        outputType = 'sequence';
                    end
                end
            end
            if strcmp(outputType,'NucleicAcid')
                rcNA = objArray; % copy object array initially
            end
            for j = 1:numel(objArray)
                rc{j} = objArray(j).BareSequence;
                for m = 1:length(rc{j})
                    n = length(rc{j})-m+1;
                    base = objArray(j).BareSequence{n};
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
                    rc{j}{m}=comp;
                end
                if strcmpi(type,'RNA')
                    for m = 1:length(rc{j})
                        rc{j}{m} = ['r',rc{j}{m}];
                    end
                end
                if strcmpi(outputType,'char')
                    rc{j} = horzcat(rc{j}{:}); % Convert to string (actually 1D char)
                elseif strcmpi(outputType,'NucleicAcid')
                    rcNA(j) = NucleicAcid(rc{j});
                end
            end
            if strcmpi(outputType,'NucleicAcid')
                rc = rcNA;
            end
            if isa(rc,'cell') && numel(objArray)==1
                rc = rc{1};
            end
        end
        function L = len(objArray)
            L = zeros(size(objArray));
            for n = 1:numel(objArray)
                L(n) = length(objArray(n).BareString);
            end
        end
        function fGC = gcContent(obj)
            nGC = 0;
            for n = 1:length(obj.BareSequence)
                if strcmpi(obj.BareSequence(n),'G') || strcmpi(obj.BareSequence(n),'C')
                    nGC = nGC + 1;
                end
            end
            fGC = nGC/length(obj.BareSequence);
        end
    end
end