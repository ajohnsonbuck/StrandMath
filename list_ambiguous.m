function list_ambiguous(seqs_orig,ia,seq_names,degenerate_sequences)
    shown = [];
    seqs_orig = seqs_orig(ia);
    count = 0;
    for n = 1:length(seqs_orig)
        if isempty(degenerate_sequences{n})==0 && sum(shown==n)==0
            count = count+1;
            fprintf(1,"Showing degenerate sequence group %3u:\n",count)
            fprintf(1,strcat(seq_names{n},"  ",(seqs_orig{n}),"\n"));
            shown = cat(1,shown,n);
            ind = degenerate_sequences{n};
            for p = 1:length(ind)
                fprintf(1,strcat(seq_names{ind(p)},"  ",(seqs_orig{ind(p)}),"\n"));
                shown = cat(1,shown,ind(p));
            end
            fprintf(1,"\n");
            usr1 = input("q to quit, any other key for next group:","s");
            if strcmp(usr1,'q')
                break
            end
        end
    end
end