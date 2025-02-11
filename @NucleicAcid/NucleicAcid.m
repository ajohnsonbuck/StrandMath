classdef NucleicAcid
    properties
        Name = '';
        String = '';
        Sequence = {};
    end
    properties (Hidden)
        BareString = '';
        BareSequence = {};
        Modifications = {};
        Mask = '';
        UnmaskedString = '';
    end
    properties (Hidden, SetAccess = private)
        Modlist = {'+','b','r'};
    end
    methods
        function obj = NucleicAcid(varargin) % Constructor
            if numel(varargin)>0
                if strcmpi(varargin{1},'random')
                    L = 10;
                    fGC = 0.5;
                    if length(varargin)>1
                        for n = 2:2:length(varargin)
                            if strcmpi(varargin{n},'length') || strcmpi(varargin{n},'size')
                                L = varargin{n+1};
                            elseif strcmpi(varargin{n},'GCcontent') || strcmpi(varargin{n},'GC_content') || strcmpi(varargin{n},'fGC')
                                fGC = varargin{n+1};
                            end
                        end
                    end
                    seq = randomSequence(obj, L, fGC);
                    obj = fromString(obj, seq);
                else
                    seq = varargin{1};
                    if isa(seq,'char') || isa(seq,'string')
                        obj = fromString(obj, seq);
                    else
                        obj = fromSequence(obj, seq);
                    end
                end
            end
            if length(varargin)>1
                for n = 2:2:length(varargin)
                    if strcmpi(varargin{n},'Mask')
                       obj.Mask = varargin{n+1};
                    elseif strcmpi(varargin{n},'Name')
                       obj.Name = varargin{n+1};
                    end
                end
            end
            if isempty(obj.Mask)
                for n = 1:obj.len
                    obj.Mask = [obj.Mask, 'n'];
                end
            end
            obj = applyMask(obj,obj.Mask);
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
        function seq = randomSequence(obj,L,fGC) % Generate random DNA sequence of length L with fractional GC content fGC
            AT = {'A','T'};
            GC = {'G','C'};
            nGC = round(fGC*L);
            nAT = L-nGC;
            ind = unidrnd(2,nGC,1);
            seq = '';
            for n = 1:nGC % Add specified number of Gs and/or Cs
                seq = [seq, GC{ind(n)}];
            end
            ind = unidrnd(2,nAT,1);
            for n = 1:nAT % Populate remainder of sequence with As and Ts
                seq = [seq, AT{ind(n)}];
            end
            seq = seq(randperm(length(seq))); % Shuffle sequence
        end
        function obj = fromString(obj,str1) % Populate object from input char or string
            obj.String = char(str1); % Convert string to char
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
        function objArray = toLNA(objArray)
            for n = 1:numel(objArray)
                seq1 = objArray(n).BareSequence;
                for p = 1:length(seq1)
                    seq1{p} = ['+',seq1{p}];
                end
                objArray(n) = fromSequence(objArray(n),seq1);
            end
        end
        function objArray = toBNA(objArray)
            for n = 1:numel(objArray)
                seq1 = objArray(n).BareSequence;
                for p = 1:length(seq1)
                    seq1{p} = ['b',seq1{p}];
                end
                objArray(n) = fromSequence(objArray(n),seq1);
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
        function fGC = gc(obj) % Alias for gcContent()
            fGC = gcContent(obj);
        end
        function objArray = applyMask(objArray,mask)
            for n = 1:numel(objArray)
                str1 = objArray(n).String;
                for p = 1:objArray(n).len
                    if strcmp(mask(p),'-')
                       objArray(n).Sequence{p}='-';
                    end
                end
                objArray(n) = objArray(n).fromSequence;
                objArray(n).UnmaskedString = str1; % Retain unmasked sequence for reference
            end
        end
        function print(objArray)
            for n = 1:numel(objArray)
                if ~isempty(objArray(n).Name)
                    fprintf(1,'\n%s',objArray(n).Name);
                end
                fprintf(1,['\n5-',objArray(n).String,'-3\n']);
            end
            fprintf(1,'\n');
        end
    end
end