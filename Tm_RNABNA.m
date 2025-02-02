function Tm = Tm_RNABNA(rseq,conc,Na,Mg)
    [dH0, dS0, dG0] = parameters_RNABNA(rseq);

    % % Assume dS0 decreases by 3.33 for each BNA modification -- based on
    % % estimate that deltaG drops by ~1000 cal/mol for each increase in ~6
    % % deg C to the melting temp
    % dS0 = dS0+4.9*length(rseq);

    Tm = Tm_1MNa(dH0,dS0,conc);

    Tm = salt_correction(Tm,length(rseq),gc_content(rseq),Na,Mg);

    % Tm = Tm + 6*length(rseq); % Assumes each BNA adds 6 deg C to Tm, per Abdur Rahman, S. M.; Seki, S.; Obika, S.; Yoshikawa, H.; Miyashita, K.; Imanishi, T. Design, Synthesis, and Properties of 2′,4′-BNANC: A Bridged Nucleic Acid Analogue. J. Am. Chem. Soc. 2008, 130 (14), 4886−4896.
end