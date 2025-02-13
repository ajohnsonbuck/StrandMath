% terminal = {'rC','rC','rG','rG'};
% 
% bases = {'rA','rC','rG','rU'};
% c=1;
% seqs = cell(16,1);
% for x = 1:4
%     for y = 1:4
%         seqs{c} = [terminal{1},bases{x},terminal{2},'/',terminal{3},bases{y},terminal{4}];
%         c = c+1;
%     end
% end
% 
% n = reshape(m',16,1);
% 
% out = horzcat(seqs,num2cell(n));


% Semi-automated

vert_offset = 3;
interval = 16;
vert_end = height(sheet);

bases = {'rA','rC','rG','rU'};

out = cell(0,2);

for start = vert_offset:interval:vert_end

    terminal = cell(1,4);
    row1 = erase(sheet{start,3},' ');
    row2 = erase(sheet{start+1,3},' ');
    terminal{1} = strcat('r', row1(1));
    terminal{2} = strcat('r', row1(2));
    terminal{3} = strcat('r', row2(1));
    terminal{4} = strcat('r', row2(2));

    m = sheet(start+8:start+11,2:5);

    c=1;
    seqs = cell(16,1);
    for x = 1:4
        for y = 1:4
            seqs{c} = [terminal{1},bases{x},terminal{2},'/',terminal{3},bases{y},terminal{4}];
            c = c+1;
        end
    end

    n = reshape(m',16,1);

    out = vertcat(out,(horzcat(seqs,n)));

end