function [mismatches, overhangs, startpos, endpos] = find_complementarity(probe, target)
    target = target.toDNA(); % Convert target to DNA to simplify assessing complementarity
    % Should repace with a Duplex or Interaction class
    if ~isempty(target)
        %find register of greatest complementarity
        r = -Inf; % Initialize register of greatest complementarity
        comp.probe = zeros(1,probe.length());
        comp.target = zeros(1,target.length());
        for o = 1-probe.length():target.length()-1 % Range of registers to probe
            s1 = max([1,-o+1]);
            s2 = min([probe.length(), target.length()-o]);
            t1 = target.length()-min([target.length(), probe.length()+o])+1;
            t2 = target.length()-max([1,1+o])+1;
            probecomp = zeros(1,probe.length());
            probecomp(s1:s2) = probe.BareString(s1:s2)==reverse_complement(target.BareString(t1:t2));
            targetcomp = zeros(1,target.length());
            targetcomp(t1:t2) = fliplr(reverse_complement(target.BareString(t1:t2))==probe.BareString(s1:s2));
            if sum(probecomp) > sum(comp.probe)
                comp.probe = probecomp;
                comp.target = targetcomp;
                r = o;
            end
        end
        % determine overhangs and mismatches
        overhangs.probe = {'',''};
        overhangs.target = {'',''};
        % Find probe overhangs
        for n = 2:probe.length()
            if comp.probe(n) == 1
                if comp.probe(n-1) == 0
                    overhangs.probe{1} = probe.BareString(n-1:n);
                else
                    break
                end
            end
        end
        for n = 2:probe.length()
            p = probe.length()-n+1;
            if comp.probe(p) == 1
                if comp.probe(p+1) == 0
                    overhangs.probe{2} = probe.BareString(p:p+1);
                else
                    break
                end
            end
        end
        % Find target overhangs
        for n = 2:target.length()
            if comp.target(n) == 1
                if comp.target(n-1) == 0
                    overhangs.target{1} = target.BareString(n-1:n);
                else
                    break
                end
            end
        end
        for n = 2:target.length()
            p = target.length()-n+1;
            if comp.target(p) == 1
                if comp.target(p+1) == 0
                    overhangs.target{2} = target.BareString(p:p+1);
                else
                    break
                end
            end
        end
        mismatches = {};
        % Check for terminal mismatches
        if ~isempty(overhangs.target{2}) && ~isempty(overhangs.probe{1}) % terminal mismatch on 5' end of probe
            mismatches = horzcat(mismatches,strcat('5t',overhangs.probe{2},'x',overhangs.target{1}(1)));
            overhangs.target{2} = '';
            overhangs.probe{1} = '';
        end
        if ~isempty(overhangs.target{1}) && ~isempty(overhangs.probe{2}) % terminal mismatch on 3' end of probe
            mismatches = horzcat(mismatches,strcat('3t',overhangs.probe{2},'x',overhangs.target{1}(1)));
            overhangs.target{1} = '';
            overhangs.probe{2} = '';
        end
        % find first and last complementary nucleotide of probe
        startpos = 0;
        n = 1;
        while startpos == 0 && n < length(comp.probe)
            if comp.probe(n) == 1
                startpos = n;
            else
                n=n+1;
            end
        end
        n = probe.length();
        endpos = Inf;
        while endpos == Inf && n > 0
            if comp.probe(n) == 1
                endpos = n;
            else
                n = n-1;
            end
        end
        if startpos>0 && endpos < Inf
        for n = startpos:endpos
            if comp.probe(n) == 0
                mismatch = '';
                for p = n:n+1
                    mismatch = strcat(mismatch,probe.Sequence{p});
                end
                mismatch = strcat(mismatch,'/');
                mismatch = strcat(mismatch,target.Sequence{target.length()-r-n+1});
                mismatch = strcat(mismatch,target.Sequence{target.length()-r-n});
                mismatches = horzcat(mismatches,mismatch);
            end
        end
        end
    end
end