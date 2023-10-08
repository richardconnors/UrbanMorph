
% Set the location of JULIA code and data files
% you can just do 
% dirLocation = "C:\blah\blah\blah"

% For me this depends on whether I am using my home PC or work laptop so I
% check the name of the computer and set the location accordingly
[~, thisPC] = system ('hostname');
if strcmpi(strtrim(thisPC),"PP0695") % this is my laptop
  % dirLocation = "C:\Users\richard.connors\Dropbox\_LuxWork\MEVERST\code\TY_Flexbus\";
else
  dirLocation = "Q:\REPOS\flexbus_v0.6\";
end

% so now I have the location of the "main folder" I concatenate the filename:
% the main JULIA file
code_main = append(replace(dirLocation, "\","\\"), 'example.jl'); % join strings with no space

% the directory with all the data FOLDERS I want to cycle through
data_dir = join([dirLocation, 'data\'],"");

% Use the 'dir' function to list all files and folders in the specified directory
listing = dir(data_dir);
% Initialize an empty cell array to store subdirectories matching the pattern
dataFolders = {};
% Define the regular expression pattern to match directory names
pattern = 'P\d+S\d+C\d+DP\d+\.\d+';
% Iterate through the 'listing' and identify subdirectories matching the pattern
for i = 1:length(listing)
    % Check if the item is a directory, not "." or "..", and matches the pattern
    if listing(i).isdir && ~ismember(listing(i).name, {'.', '..'}) && ~isempty(regexp(listing(i).name, pattern, 'once'))
        % Add the full path of the matching subdirectory to the cell array
        dataFolders{end+1} = fullfile(data_dir, listing(i).name);
    end
end
% 'dataFolders' now contains a list of full paths to subdirectories matching the pattern
disp(dataFolders');
nDataFolders = length(dataFolders);

for ff = 1%:nDataFolders
  % need the path with double \\ for JULIA
  F1 = append(replace(dataFolders{1}, "\","\\"),"\\");  
  % run JULIA via the windows command
  cd(dirLocation) % because inside Julia it assumes we are in the REPO
  [status,cmdout]=system(join(["julia", code_main, F1]))

  %  Here are some hard-coded filenames if needed to check things are
  %  working - you need to set the path for your computer of course!
  %   [status,cmdout]=system("julia Q:\Dropbox\_LuxWork\MEVERST\code\Flexbus\main.jl " + ...
  %     "Q:\\Dropbox\\_LuxWork\\MEVERST\\code\\Flexbus\\data\\converted_Peak\\c10.txt " + ...
  %     "Q:\\Dropbox\\_LuxWork\\MEVERST\\code\\Flexbus\\data\\converted_Peak\\Timetable.csv")

end


