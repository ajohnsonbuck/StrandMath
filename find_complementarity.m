function [mismatches, overhangs, startpos, endpos] = find_complementarity(seq, target)
    trc = reverse_complement(target);
    [seq2,mods] = parse_modifications(seq);
    if ~isempty(target)
        %find register of greatest complementarity
        r = -Inf; % Initialize register of greatest complementarity
        comp.seq = zeros(1,length(seq2));
        comp.target = zeros(1,length(target));
        ncomp = 0;
        for o = 1-length(seq2):length(trc)-1 % Range of registers to probe
            s1 = max([1,-o+1]);
            s2 = min([length(seq2), length(target)-o]);
            t1 = length(target)-min([length(target), length(seq2)+o])+1;
            t2 = length(target)-max([1,1+o])+1;
            seq2comp = zeros(1,length(seq2));
            seq2comp(s1:s2) = seq2(s1:s2)==reverse_complement(target(t1:t2));
            targetcomp = zeros(1,length(target));
            targetcomp(t1:t2) = fliplr(reverse_complement(target(t1:t2))==seq2(s1:s2));
            if sum(seq2comp) > sum(comp.seq)
                comp.seq = seq2comp;
                comp.target = targetcomp;
                ncomp = sum(seq2comp);
                r = o;
            end
        end
        % determine overhangs and mismatches
        overhangs.seq = {'',''};
        overhangs.target = {'',''};
        % Find probe overhangs
        for n = 2:length(seq2)
            if comp.seq(n) == 1
                if comp.seq(n-1) == 0
                    overhangs.seq{1} = seq2(n-1:n);
                else
                    break
                end
            end
        end
        for n = 2:length(seq2)
            p = length(seq2)-n+1;
            if comp.seq(p) == 1
                if comp.seq(p+1) == 0
                    overhangs.seq{2} = seq2(p:p+1);
                else
                    break
                end
            end
        end
        % Find target overhangs
        for n = 2:length(target)
            if comp.target(n) == 1
                if comp.target(n-1) == 0
                    overhangs.target{1} = target(n-1:n);
                else
                    break
                end
            end
        end
        for n = 2:length(target)
            p = length(target)-n+1;
            if comp.target(p) == 1
                if comp.target(p+1) == 0
                    overhangs.target{2} = target(p:p+1);
                else
                    break
                end
            end
        end
        mismatches = {};
        % find first and last complementary nucleotide of probe
        startpos = 0;
        n = 1;
        while startpos == 0
            if comp.seq(n) == 1
                startpos = n;
            else
                n=n+1;
            end
        end
        n = length(seq2);
        endpos = Inf;
        while endpos == Inf
            if comp.seq(n) == 1
                endpos = n;
            else
                n = n-1;
            end
        end
        for n = startpos:endpos
            if comp.seq(n) == 0
                mismatch = '';
                for p = n:n+1
                if ~isempty(mods{p})
                    mismatch = strcat(mismatch,mods{p});
                end
                mismatch = strcat(mismatch,seq2(p));
                end
                mismatch = strcat(mismatch,'/');
                mismatch = strcat(mismatch,target(length(target)-r-n+1));
                mismatch = strcat(mismatch,target(length(target)-r-n));
                mismatches = horzcat(mismatches,mismatch);
            end
        end
    end
end