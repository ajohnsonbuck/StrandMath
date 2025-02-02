function Tm = Tm_NA(seq,type,conc,Na,Mg)
    [dH0, dS0, ~] = parameters_NA(seq,type,37);

    Tm = Tm_1MNa(dH0,dS0,conc);

    Tm = salt_correction(Tm,NA_length(seq),gc_content(seq),Na,Mg);
end