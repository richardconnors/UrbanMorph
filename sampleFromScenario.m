function [T_busFleet, T_Passenger,T_busStop,T_Charger,T_Station,allStationDeps,T_depot] = sampleFromScenario(nSample, p, T_busFleet, T_Passenger,T_busStop,T_Charger,T_Station,allStationDeps,T_depot)


% to ensure always get nSample passengers
thisPax = randperm(p.nPax);
thisPax = sort(thisPax(1:nSample));
T_Passenger = T_Passenger(thisPax,:);

% ====== BUS FLEET
% for each type of bus downsize fleet by ratio of nSample to nPax



pSchedule = zeros(nSample,nDays);
this_Pax = sort(randperm(nPassengers,nSample));
pSchedule(:,i) = this_Pax ;
T_Passenger(this_Pax,:);


%===========================================================================
return

thisFolder = 'C:\Users\richard.connors\Documents\REPOS\UrbanMorph\test_data\P150_S1_R40.0_W1.0_MP01.4_SS03.0_DP2.0_r04';
folderItems = dir(thisFolder); allFiles = string({folderItems.name});
cus_file = allFiles(contains(allFiles,'cus_'));
route_files =  allFiles(startsWith(allFiles,'route_detail_'));
% load all the usual scenario tables and params
params_file = fullfile(thisFolder, allFiles(contains(allFiles,'parameters'))); p = readtable(params_file);
filename = fullfile(thisFolder, allFiles(contains(allFiles,'passenger'))); T_Passenger = readtable(filename);
filename = fullfile(thisFolder, allFiles(contains(allFiles,'busFleet'))); T_busFleet = readtable(filename);
filename = fullfile(thisFolder, allFiles(contains(allFiles,'busStopXY'))); T_busStop = readtable(filename);
filename = fullfile(thisFolder, allFiles(contains(allFiles,'chargerXY')));  T_Charger = readtable(filename);
filename = fullfile(thisFolder, allFiles(contains(allFiles,'stationXY'))); T_Station = readtable(filename);
filename = fullfile(thisFolder, allFiles(contains(allFiles,'transitTimetable'))); allStationDeps = readtable(filename);
filename = fullfile(thisFolder, allFiles(contains(allFiles,'depotXY'))); T_depot = readtable(filename);








