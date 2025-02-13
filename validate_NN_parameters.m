function rmsd = validate_NN_parameters()
    rmsd = zeros(3,1) + Inf;

    % Load datasets
    load('validation_Sugimoto_etal_1995.mat'); %DNA/RNA
    load('validation_Owczarzy_etal_2011.mat'); %DNA/DNA and LNA/DNA
    load('validation_Xia_etal_1998.mat'); % RNA/RNA

    for n = 1:numel(SugimotoTable2.seqs)
        SugimotoTable2.Tm_pred(n) = estimate_Tm(SugimotoTable2.seqs(n),'target',toDNA(SugimotoTable2.seqs(n)'),'conc',100E-6);
    end

    for n = 1:numel(OwczarzyTable4.seqs)
        OwczarzyTable4.Tm_pred(n) = estimate_Tm(OwczarzyTable4.seqs(n),'target',OwczarzyTable4.comp(n),'conc',2E-6);
    end

    for n = 1:numel(XiaTable1.seqs)
        if XiaTable1.seqs(n).isSymmetric
            conc = 100E-6; % for symmetric sequences
        else
            conc = 200E-6; % for non-symmetric sequences
        end
        XiaTable1.Tm_pred(n,1) = estimate_Tm(XiaTable1.seqs(n),'target',XiaTable1.seqs(n).reverseComplement.toRNA,'conc',conc);
    end

    rmsd(1) = sqrt(mean((SugimotoTable2.Tm_pred - SugimotoTable2.Tms).^2));
    rmsd(2) = sqrt(mean((OwczarzyTable4.Tm_pred - OwczarzyTable4.Tms).^2));
    rmsd(3) = sqrt(mean((XiaTable1.Tm_pred - XiaTable1.Tms).^2));

    uitable("RowName",{'Sugimoto et al. 1995','Owczarzy et al. 2011', 'Xia et al. 1998'},...
        "ColumnName", strcat('RMSD (',char(176),'C)'),"Data",rmsd,'Units','normalized','Position',[0.05 0.05 0.9 0.9]);

    % Results 2025-02-12:
    %           RMSD (deg C)
    % Sugimoto  0.0250
    % Owczarzy  0.0189
    % Xia       0.0365
    
end