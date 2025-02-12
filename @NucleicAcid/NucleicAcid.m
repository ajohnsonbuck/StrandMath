classdef NucleicAcid
    properties
        Name = '';
        String = '';
        Sequence = {};
    end
    properties (Hidden)
        Mask = '';
        UnmaskedString = '';
    end
    properties (Hidden, SetAccess = private)
        Modlist = {'+','b','r'};
    end
    methods
        function objArray = NucleicAcid(varargin) % Constructor
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
                    seq = NucleicAcid.randomSequence(L, fGC);
                    objArray(1) = fromString(objArray(1), seq);
                else
                    seq = varargin{1};
                    if isa(seq,'char') || isa(seq,'string')
                        objArray(1) = fromString(objArray(1), seq);
                    elseif isa(seq,'cell') && size(seq,1)>1 && size(seq,2)==1 % if input argumeent is a vertical cell, assume those cells contain sequences
                        objArray(1,numel(seq)) = NucleicAcid(); % preallocate object array
                        for n = 1:numel(seq)
                            objArray(n) = fromString(objArray(1), seq{n,1});
                        end
                    else
                        objArray(1) = fromSequence(objArray(1), seq);
                    end
                end
            end
            if length(varargin)>1
                for n = 2:2:length(varargin)
                    if strcmpi(varargin{n},'Mask')
                        for p = 1:numel(objArray)
                            objArray(p).Mask = varargin{n+1};
                        end
                    elseif strcmpi(varargin{n},'Name')
                        if numel(objArray)==1
                            objArray(1).Name = varargin{n+1};
                        else
                        for p = 1:numel(objArray)
                            objArray(p).Name = varargin{n+1}{p};
                        end
                        end
                    end
                end
            end
            for m = 1:numel(objArray)
                if isempty(objArray(m).Mask)
                    for n = 1:objArray(m).len
                        objArray(m).Mask = [objArray(m).Mask, 'n'];
                    end
                end
            end
            for p = 1:numel(objArray)
                objArray(p) = applyMask(objArray(p),objArray(p).Mask);
            end
        end
        function obj = fromString(obj,str1) % Populate object from input char or string
            obj.String = erase(char(str1),{' ',char("5'-"),char("-3'"),'5-','-3',char("5'"),char("3'")}); % Convert string to char and remove empty spaces and termini
            obj.Sequence = cell(1,length(obj.bareString));
            mods = obj.modifications;
            str1 = obj.bareString;
            for n = 1:length(obj.Sequence)
                obj.Sequence{n} = [obj.modifications{n}, str1(n)];
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
        end
        function str = bareString(objArray)
            str = cell(numel(objArray),1);
            for n = 1:numel(objArray)
                str{n} = erase(objArray(n).String,objArray(n).Modlist); % Remove modification prefixes from sequences
            end
            if numel(str)==1
                str = str{:};
            end
        end
        function bareSeqs = bareSequence(objArray)
            bareSeqs = cell(numel(objArray),1);
            for n = 1:numel(objArray)
                bareSeqs{n} = cell(1,objArray(n).len);
                str = objArray(n).bareString;
                for p = 1:length(str)
                    bareSeqs{n}{1,p} = str(p);
                end
            end
            if numel(objArray) == 1
                bareSeqs = bareSeqs{1};
            end
        end
        function mods = modifications(objArray)
            mods = cell(numel(objArray),1);
            for n = 1:numel(objArray)
                mods{n} = cell(1,objArray(n).len);
                c = 1;
                for p = 1:length(objArray(n).String)
                    if ismember(objArray(n).String(p),objArray(n).Modlist)
                        mods{n}{c} = objArray(n).String(p);
                    else
                        c = c+1;
                    end
                end
            end
            if numel(objArray) == 1
                mods = mods{1};
            end
        end
        function objArray = toDNA(objArray)
            for n = 1:numel(objArray)
                str1 = replace(objArray(n).bareString, 'U', 'T');
                objArray(n) = fromString(objArray(n),str1);
            end
        end
        function objArray = toLNA(objArray)
            for n = 1:numel(objArray)
                seq1 = objArray(n).bareSequence;
                for p = 1:length(seq1)
                    seq1{p} = ['+',seq1{p}];
                end
                objArray(n) = fromSequence(objArray(n),seq1);
            end
        end
        function objArray = toBNA(objArray)
            for n = 1:numel(objArray)
                seq1 = objArray(n).bareSequence;
                for p = 1:length(seq1)
                    seq1{p} = ['b',seq1{p}];
                end
                objArray(n) = fromSequence(objArray(n),seq1);
            end
        end
        function objArray = toRNA(objArray)
            for m = 1:numel(objArray)
                seq1 = replace(objArray(m).bareSequence, 'T', 'U');
                for n = 1:length(seq1)
                    if contains(seq1{n}(end),{'A','T','U','G','C'})
                        seq1{n} = strcat('r',seq1{n});
                    end
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
            keepName = false; % Default: wipe name of reverse complement unless 'keepname' specified
            rc = cell(size(objArray));
            if ~isempty(varargin)
                for n = 1:length(varargin)
                    if strcmpi(varargin{n},'RNA')
                        type = 'RNA';
                    elseif strcmpi(varargin{n},'char') || strcmpi(varargin{n},'string')
                        outputType = 'char';
                    elseif strcmpi(varargin{n},'sequence') || strcmpi(varargin{n},'cell')
                        outputType = 'sequence';
                    elseif strcmpi(varargin{n},'keepName')
                        keepName = true;
                    end
                end
            end
            if strcmp(outputType,'NucleicAcid')
                rcNA = objArray; % copy object array initially
            end
            for j = 1:numel(objArray)
                rc{j} = objArray(j).bareSequence;
                for m = 1:length(rc{j})
                    n = length(rc{j})-m+1;
                    base = objArray(j).bareSequence{n};
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
                    if keepName
                        name = rcNA(j).Name;
                    else
                        name = [rcNA(j).Name,'_reverseComplement'];
                    end
                    rcNA(j) = NucleicAcid(rc{j},'name',name);
                end
            end
            if strcmpi(outputType,'NucleicAcid')
                rc = rcNA;
            end
            if isa(rc,'cell') && numel(objArray)==1
                rc = rc{1};
            end
        end
        function objArray = scramble(objArray)
            for n = 1:numel(objArray)
                seq = objArray(n).Sequence;
                name = objArray(n).Name;
                objArray(n) = NucleicAcid(seq(randperm(numel(objArray(n).Sequence))));
                objArray(n).Name = strcat(name,'_scrambled');
            end
        end
        function L = len(objArray)
            L = zeros(size(objArray));
            for n = 1:numel(objArray)
                L(n) = length(objArray(n).bareString);
            end
        end
        function fGC = gcContent(obj)
            nGC = 0;
            for n = 1:length(obj.bareSequence)
                if strcmpi(obj.bareSequence(n),'G') || strcmpi(obj.bareSequence(n),'C')
                    nGC = nGC + 1;
                end
            end
            fGC = nGC/length(obj.bareSequence);
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
                objArray(n).Mask = mask;
                objArray(n).UnmaskedString = str1; % Retain unmasked sequence for reference
            end
        end
        function objArray = unmask(objArray)
            for n = 1:numel(objArray)
                str1 = objArray(n).UnmaskedString;
                objArray(n) = objArray(n).fromString(str1);
                objArray(n).Mask = replace(objArray(n).Mask,'-','n');
            end
        end
        function print(objArray,varargin)
            printBare = false;
            if numel(varargin)>0
                for n = 1:numel(varargin)
                    if strcmpi(varargin{n},'bare')
                        printBare = true;
                    end
                end
            end
            for n = 1:numel(objArray)
                if ~isempty(objArray(n).Name)
                    fprintf(1,'\n%s',objArray(n).Name);
                end
                if printBare
                    str1 = objArray(n).bareString;
                else
                    str1 = objArray(n).String;
                end
                fprintf(1,[char("\n5'-"),str1,char("-3'\n")]);
            end
            fprintf(1,'\n');
        end
        function obj = plus(obj1,obj2)
            if isa(obj1,'NucleicAcid') && isa(obj2,'NucleicAcid')
                obj = NucleicAcid();
                str1 = strcat(obj1.String,obj2.String);
                obj = obj.fromString(str1);
                if ~isempty(obj1.Name) && ~isempty(obj2.Name)
                    obj.Name = [obj1.Name,' + ', obj2.Name];
                else
                    obj.Name = obj1.Name;
                end
            end
        end
        function obj = mtimes(obj1,c)
            obj = NucleicAcid();
            if isa(obj1,'NucleicAcid') && isnumeric(c) && (round(c,0)==c) && c > 0
                str = '';
                for n = 1:c
                    str = [str,obj1.String];
                end
                obj = obj.fromString(str);
                obj.Name = [obj1.Name,' x ',num2str(c)];
            end
        end
        function ans = isSymmetric(objArray)
           for n = 1:numel(objArray)
                ans(n) = false;
                if cellfun(@strcmp,objArray(n).toDNA.Sequence,objArray(n).reverseComplement.toDNA.Sequence)
                    ans(n) = true;
                end
           end
        end
    end
    methods (Static)
        function seq = randomSequence(L,fGC) % Generate random DNA sequence of length L with fractional GC content fGC

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
    end
end