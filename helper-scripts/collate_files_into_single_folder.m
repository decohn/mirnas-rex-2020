% Input: source_folder is the char name of a folder that contains many 
% subfolders, each of which may contain one or more files.

% When run, this will delete all subfolders of source_folder, and move
% any files contained within the subfolders into source_folder itself.

% Returns: null.

function [] = collate_files_into_single_folder(source_folder)
    
    % when initially called, init a tracker for current_folder (index 2)
    if ~isa(source_folder, 'cell')
        source_folder = {source_folder, source_folder};
    end
    
    subfolders = dir(source_folder{2});
    
    for i = 3:size(subfolders, 1)
        if ~subfolders(i).isdir
            % base case
            if ~isequal(source_folder{1}, source_folder{2})
                movefile([source_folder{2}, filesep, subfolders(i).name], [source_folder{1}, filesep, subfolders(i).name]);
            end
        else
            % recursive case
            collate_files_into_single_folder({source_folder{1}, [source_folder{2}, filesep, subfolders(i).name]});
            rmdir([source_folder{2}, filesep, subfolders(i).name]);
        end
    end
end