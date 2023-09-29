
% Set the location of JULIA code and data files
% you can just do 
% dirLocation = "C:\blah\blah\blah"

% For me this depends on whether I am using my home PC or work laptop so I
% check the name of the computer and set the location accordingly
[~, thisPC] = system ('hostname');
if strcmpi(strtrim(thisPC),"PP0695") % this is my laptop
  dirLocation = "C:\Users\richard.connors\Dropbox\_LuxWork\MEVERST\code\TY_Flexbus\";
else
  dirLocation = "Q:\REPOS\flexbus_v0.6\";
end

% so now I have the location of the "main folder" I concatenate the
% filename:
% the main JULIA file
code_main = join([dirLocation, 'example.jl'],""); % join strings with no space

% the directory with all the data files I want to cycle through
data_dir = join([dirLocation, 'data\non_peak\'],"");
% grab all the .txt files in this data_dir
all_data = dir(join([data_dir '*.txt'],""));
% grab what should be the ONLY .csv file in the folder
data_timetable = dir(join([data_dir '*.csv'],""));

% NOTE: the for loop should run TY's JULIA code for each .txt file
% However I guess each time this would overwrite the results from the
% previous run!
% Right now this is not a problem since it doesn't even run once 🤣
% TO DO: update results .csv files to include name of calling .txt file

for ff = 1%:length(all_data) % for each .txt file in the data_dir
  % need the path with double \\ for JULIA
  F1 = join(['"',replace(data_dir, "\","\\"),all_data(ff).name,'"'],"");
  F2 = join(['"',replace(data_dir, "\","\\"),data_timetable(1).name,'"'],"");
  % run JULIA via the windows command
  setpwd(dirLocation) % because inside Julia it assumes we are in the REPO
  [status,cmdout]=system(join(["julia", code_main, F1, F2]))

  %  Here are some hard-coded filenames if needed to check things are
  %  working - you need to set the path for your computer of course!
  %   [status,cmdout]=system("julia Q:\Dropbox\_LuxWork\MEVERST\code\Flexbus\main.jl " + ...
  %     "Q:\\Dropbox\\_LuxWork\\MEVERST\\code\\Flexbus\\data\\converted_Peak\\c10.txt " + ...
  %     "Q:\\Dropbox\\_LuxWork\\MEVERST\\code\\Flexbus\\data\\converted_Peak\\Timetable.csv")

end


