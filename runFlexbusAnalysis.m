rng('default') % for reproducability

% sort out folder and file names
repo_folder = 'Q:\REPOS\';
[~, hostname] = system('hostname');
if contains(hostname, 'PP0695'),   repo_folder = 'C:\Users\richard.connors\Documents\REPOS\'; end
repo_flexbus = [repo_folder, 'Flexbus3_v0.7\'];
saveData_UrbanMorph = [repo_folder, 'UrbanMorph\data\'];
data_repo = [repo_flexbus, 'data\P100_10_150\'];
load([data_repo,'scenarioWorkspace'])
data_repo = [repo_flexbus, 'data\P100_10_150\'];

T = T_Passenger;
dFromO= sqrt(T.passenger_X.^2 + T.passenger_Y.^2);
T = [T, table(dFromO, 'VariableName', {'DistanceFromStation'})];

% get list of all folders we need to iterate through
items = dir(data_repo); nDay = 1;
for i = 1:numel(items)
  % Check if the item is a directory and not '.' or '..'
  if items(i).isdir && ~ismember(items(i).name, {'.', '..'})
    % load V_DF and cus_experience files from this folder

    thisFolder = [data_repo, items(i).name, filesep];
    folderItems = dir(thisFolder);
    allFiles = string({folderItems.name});

    V_DF_file = fullfile(thisFolder, allFiles(find(contains(allFiles,'V_DF'))));
    cus_file = fullfile(thisFolder, allFiles(find(contains(allFiles,'cus_'))));
    if ~isempty(V_DF_file) && ~isempty(cus_file)
      vdf = readtable(V_DF_file);  % You can use readmatrix for newer MATLAB versions
      this_pax = vdf(strcmp(vdf.v_type,'customer'),:); % has columns of x and y data
      cus_data = readtable(cus_file);  % You can use readmatrix for newer MATLAB versions

      originalPaxID = pSchedule(:,nDay);
      this_data = nan(size(T,1),1);
      this_data(originalPaxID) = cus_data.ride_time;
      T = [T, table(this_data, 'VariableName', {items(i).name})];
      nDay = nDay + 1;
    else
      disp([items(i).name,' no VDF_ file']);
    end
  end
end
disp(i)
disp(nDay-1)
% what do we want to analyse?
[~,ind] = sort(T.DistanceFromStation);
figure
boxplot(T{ind,7:end}'); ylabel('Ride Distance')
yyaxis right
h = scatter(1:100,T{ind,6}); ylabel('Distance to station')
title('Distribution of ride time with distance fromm station')


figure;
h_bs = scatter(T_busStop.busStop_X,T_busStop.busStop_Y,'+');
h_bs.MarkerEdgeColor = 0.8*ones(3,1);
hold on
h_pax = scatter(T_Passenger.passenger_X,T_Passenger.passenger_Y,'r.');
for i=1:10
  h_pax = scatter(T_Passenger.passenger_X(pSchedule(:,i)),T_Passenger.passenger_Y(pSchedule(:,i)),'bo');
  title(sprintf('Day %02d',i));
  pause(0.5)
  delete(h_pax)
end
