function table1 = results_table_rna(seqs_orig,seq_names,results)
    header = horzcat({'Sequence Name', 'Sequence'},dna2rna(reshape(results.subseqs_best,1,length(results.subseqs_best))));
    data = horzcat(reshape(seq_names,length(seq_names),1),dna2rna(reshape(seqs_orig,length(seqs_orig),1)),num2cell(results.is_subset_best));
    table1 = uitable("ColumnName",header,"Data",data,'Units','normalized','Position',...
    [0.05 0.05 0.9 0.9],'ColumnWidth',{100,200,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100});
end