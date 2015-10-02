function write_directory = generate_write_directories(Data)

if isempty(Data.write_directory) || size(Data.write_directory) == 1   % if write directories need to be created
    if isempty(Data.write_directory)                             % if no write directory is specified, write into parent directory
        root = Data.parent_directory;
    else                                                         % if write directory is specified, use it
        root = Data.write_directory;
    end
    write_directory = cell(size(Data.filter_to_run));           % make write directory with entries for each filter
    for i=1:size(Data.filter_to_run)                               % generate subdirectories for each filter
        append = Data.filter_to_run{i};                          % use filter name as subdirectory name
        write_directory(i) = cellstr(fullfile(root,append,filesep));
        if Data.additive == 1                                    % if filters are additive, cascade file structure
            root = write_directory{i};
        end
    end
else if size(Data.write_directory)~=size(Data.filter_to_run)
    error('Size of write_directory does not match number of filters applied')
    else
        write_directory = Data.write_directory;
    end
end

% Make new directories if non-existent
for i=1:size(write_directory)
    if not(exist(write_directory{i},'dir'))
        % This creates the output directory
        mkdir(write_directory{i});
    end
end