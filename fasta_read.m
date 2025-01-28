function fastaData = fasta_read(fastaFile)
% Read fasta file and place contents into cell array

% Open the file for reading
fid = fopen(fastaFile, 'r');
if fid == -1
    error('Could not open the file %s. Check if it exists and is accessible.', fastaFile);
end

% Initialize cell array to store sequences
fastaData = {};
currentSeq = '';
header = '';

% Read the file line by line
while ~feof(fid)
    line = fgetl(fid);
    
    % Check if the line is a header line
    if startsWith(line, '>')
        % If there is an existing sequence, store it
        if ~isempty(currentSeq)
            fastaData{end+1, 1} = header; % Add the header
            fastaData{end, 2} = currentSeq; % Add the sequence
        end
        
        % Reset for the new sequence
        header = line; % Store the new header
        currentSeq = ''; % Reset the sequence
    else
        % Concatenate the current line to the sequence
        currentSeq = [currentSeq, line];
    end
end

% Add the last sequence if any
if ~isempty(currentSeq)
    fastaData{end+1, 1} = header;
    fastaData{end, 2} = currentSeq;
end

% Close the file
fclose(fid);

end