
% Set the location of JULIA code and data files
% you can just do
% dirLocation = "C:\blah\blah\blah"

% For me this depends on whether I am using my home PC or work laptop so I
% check the name of the computer and set the location accordingly
[~, thisPC] = system ('hostname');
if strcmpi(strtrim(thisPC),"PP0695") % this is my laptop
  dirLocation = "C:\Users\richard.connors\Documents\REPOS\flexbus_v0.6\";
else
  dirLocation = "Q:\REPOS\flexbus_v0.6\";
end

% so now I have the location of the "main folder" I concatenate the
% filename:
% the main JULIA file
code_main = join([dirLocation, 'example.jl'],""); % join strings with no space

% the directory with all the data files I want to cycle through
data_dir = join([dirLocation, 'data\'],"");

% how many scenario drectories are in this folder?
% Specify the directory you want to search (replace with your desired path)
subdirectories = dir(data_dir);
subdirectories = subdirectories([subdirectories.isdir]);

% Define the pattern for matching directory names
pattern = 'P(\d+)S(\d+)C(\d+)DP([\d.]+)';
% Initialize an empty cell array to store matching directory names
matchingDirectories = {};
% Loop through the subdirectories and check if they match the pattern
for i = 1:length(subdirectories)
  dirname = subdirectories(i).name;
  if regexp(dirname, pattern)
    matchingDirectories{end+1} = dirname;
  end
end

% Display the list of matching directories
disp('Matching Directories:');
disp(matchingDirectories');
nDir = length(matchingDirectories);

data_dir = replace(data_dir, "\","\\");


for ff = 1%:nDir % for each .txt file in the data_dir
  % need the path with double \\ for JULIA
  F1 = join([data_dir,matchingDirectories{ff},'\\'],''); % join with no whitespace
  % run JULIA via the windows command
  RUN_STRING = join(["julia", code_main, F1])
  cd(dirLocation) % because inside Julia it assumes we are in the REPO
  [status,cmdout]=system(RUN_STRING)

  %  Here are some hard-coded filenames if needed to check things are
  %  working - you need to set the path for your computer of course!
  %   [status,cmdout]=system("julia Q:\Dropbox\_LuxWork\MEVERST\code\Flexbus\main.jl " + ...
  %     "Q:\\Dropbox\\_LuxWork\\MEVERST\\code\\Flexbus\\data\\converted_Peak\\c10.txt " + ...
  %     "Q:\\Dropbox\\_LuxWork\\MEVERST\\code\\Flexbus\\data\\converted_Peak\\Timetable.csv")

end


