function Tm = Tm_RNADNA(rseq,conc,Na,Mg)
    [dH0, dS0, ~] = parameters_RNADNA(rseq,37);

    Tm = Tm_1MNa(dH0,dS0,conc);

    Tm = salt_correction(Tm,length(rseq),gc_content(rseq),Na,Mg);
end