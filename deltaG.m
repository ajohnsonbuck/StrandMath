function dG = deltaG(dH0,dS0,varargin)
    % Gas constant, cal/(mol K)
    R = 1.987204258; 
    % Temperature default 37 C
    T = 37;
    % Concentration default 1 M
    c = 1;
    if ~isempty(varargin)
        for n = 1:2:length(varargin)
            if strcmpi(varargin{n},'temperature') || strcmpi(varargin{n},'T')
                T = varargin{n+1};
            elseif strcmpi(varargin{n},'concentration') || strcmpi(varargin{n}, 'conc')
                c = varargin{n+1};
            else
                fprintf(1,'Warning: deltaG() did not recognize argument "%s". Ignored this argument and any that immediately follow.  Please re-run function without this argument.', num2str(varargin{n}));
            end
        end
    end
    T = T+273.15;
    dG = dH0 - T*dS0 - R*T*log(c/4);
end