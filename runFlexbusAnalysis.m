rng('default') % for reproducability

% sort out folder and file names
repo_folder = 'Q:\REPOS\';
[~, hostname] = system('hostname');
if contains(hostname, 'PP0695'),   repo_folder = 'C:\Users\richard.connors\Documents\REPOS\'; end
repo_flexbus = [repo_folder, 'Flexbus3_v0.7\'];
saveData_UrbanMorph = [repo_folder, 'UrbanMorph\data\'];

% %=======================
% thisTEST = 'data\P100_10_150\';
% thisTEST = 'data\P100_50_30\';
% thisTEST = 'data\P100_50D30_P\';
% thisTEST = 'data\test\';
% %=======================

data_repo = [repo_flexbus, 'data\P100_50D30_P\'];
load([data_repo,'scenarioWorkspace'])
fprintf('nPassengers = %d \n', nPassengers);
fprintf('Pax_maxRadius = %.1f \n', Pax_maxRadius)
fprintf('Pax_minRadius = %.1f \n', Pax_minRadius)
fprintf('demandPeakness = %.1f \n', demandPeakness)
fprintf('nSample = %d \n', nSample)
fprintf('nDays = %d \n', nDays)

data_repo = [repo_flexbus, 'data\P100_50D30_P\'];

T = T_Passenger;
dFromO= sqrt(T.passenger_X.^2 + T.passenger_Y.^2);
T = [T, table(dFromO, 'VariableName', {'DistanceFromStation'})];
busOccupancyData = [];
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
    route_files =  allFiles(startsWith(allFiles,'route_detail_'));
    if ~isempty(V_DF_file) && ~isempty(cus_file)
      vdf = readtable(V_DF_file);  % You can use readmatrix for newer MATLAB versions
      this_pax = vdf(strcmp(vdf.v_type,'customer'),:); % has columns of x and y data
      cus_data = readtable(cus_file);  % You can use readmatrix for newer MATLAB versions
      allCus = cus_data.cus;
      originalPaxID = pSchedule(:,i-2); % i = 1,2 are '.' and '..' then count through file
      this_data = nan(size(T,1),1);
      this_rideTime = cus_data.ride_time; this_rideTime(~this_rideTime) = NaN;
      this_data(originalPaxID) = this_rideTime;
      T = [T, table(this_data, 'VariableName', {items(i).name})]; %#ok<AGROW>

      servedCus = [];
      for ff = 1:length(route_files) % one for each bus?
        this_route = readtable(fullfile(thisFolder,route_files{ff}));  % You can use readmatrix for newer MATLAB versions
        busOccupancyData = [busOccupancyData; this_route.num_passenger]; %#ok<AGROW>
        tour = [this_route.x, this_route.y]
        C = string(this_route.cus_set); C = C(~startsWith(C,"Int64"));
        pattern = '\d+';
        % Loop through each element of the string array
        for i = 1:numel(C)
          % Find all numerical values in the current string
          matches = regexp(C(i), pattern, 'match');
          % Convert the matched values to a numeric array and append them to the vector
          servedCus = [servedCus, str2double(matches{1})];
        end

        % Display the resulting vector of numerical
      end
      nDay = nDay + 1;
    else
      disp([items(i).name,' no VDF and cus file']);
    end
  end
end
% fprintf('nPassengers = %d \n', nPassengers);
fprintf('num days failed = %d \n', nDays - nDay+1);

% what do we want to analyse?
[~,ind] = sort(T.DistanceFromStation);
figure('pos',[ 162,524,1511,420])
boxplot(T{ind,7:end}'); ylabel('Ride Distance')
yyaxis right
h = scatter(1:100,T{ind,6}); ylabel('Distance to station')
title('Distribution of ride time with distance from station')

figure;
histogram(busOccupancyData(busOccupancyData>0))
title('Bus Occupancy'); ylabel('Number of Trips')

return

figure;
h_bs = scatter(T_busStop.busStop_X,T_busStop.busStop_Y,'+');
h_bs.MarkerEdgeColor = 0.8*ones(3,1);
hold on
h_st = scatter(0,0,'gs'); h_st.MarkerFaceColor = 'g';
h_pax = scatter(T_Passenger.passenger_X,T_Passenger.passenger_Y,'r.');
for i=1:10
  h_pax = scatter(T_Passenger.passenger_X(pSchedule(:,i)),T_Passenger.passenger_Y(pSchedule(:,i)),'bo');
  title(sprintf('Day %02d',i));
  pause(0.5)
  delete(h_pax)
end
