
% generateScenario allows the inputs listed below to be specified.
% You can set values here, and run this script to call the function
% Outputs are written to comma delimited text files (easy to change this)
% Output text file names currently hardcoded as below, saved to PWD
% % transitTimetable.txt
% % stationXY.txt
% % chargerXY.txt
% % paxData.txt
% % busStopXY.txt

% % All distances in KM
% % Timetable of departures for each station is hardcoded
% % For each passenger departure time is randomly assigned
% % Demand distribution is controlled and can go from uniform -> peaked
% nTransitStations = 2;
% Station_separation = 5; % distance between multiple stations
% nPassengers = 100; % number passengers to generate
% % Passengers will be generated uniformly distributed in an annulus around each station
% % minRadius and maxRadius can have single value = used for all stations
% % OR can be a row vector size [1 x nTransitStations]
% Pax_minRadius = 0.5; % min radius away from station for passenger locations
% Pax_maxRadius = 2.0; % max radius around station for passenger locations
% paxSeparation = 0.05; % min distance between passengers
% BS_separation = 1; % spacing for grid of potential bus stop locations
% maxWalkingDist = BS_separation*1.05/sqrt(2); %
% nCharger  = 4; % how many chargers PER TOWN/TRANSIT STATION.
% charger_radius = 1.5; % chargers located on circle around each town centre
% demandPeakness = 0; % 0 = uniform. 1 = peaked in middle
% PLOTFLAG = 1; % plots the network data


repo_folder = 'Q:\REPOS\';
[~, hostname] = system('hostname');
if contains(hostname, 'PP0695')
  repo_folder = 'C:\Users\richard.connors\Documents\REPOS\';
end
repo_flexbus = [repo_folder, 'Flexbus3_v0.7\']; 
saveData_UrbanMorph = [repo_folder, 'UrbanMorph\data\']; % on Kuzuri


nTransitStations = 1;
Station_separation = 1; % distance between multiple stations
nPassengers = 100; % number passengers to generate
% Passengers will be generated uniformly distributed in an annulus around each station
% minRadius and maxRadius can have single value = used for all stations
% OR can be a row vector size [1 x nTransitStations]
Pax_minRadius = 1.5; % min radius away from station for passenger locations
Pax_maxRadius = 10.0; % max radius around station for passenger locations
paxSeparation = 0.05; % min distance between passengers
BS_separation = 1; % spacing for grid of potential bus stop locations
maxWalkingDist = 1.2; % BS_separation*sqrt(2)/2; %
nCharger  = 4; % how many chargers PER TOWN/TRANSIT STATION.
charger_radius = 2; % chargers located on circle around each town centre
demandPeakness = 1; % 0 = uniform. 1 = peaked in middle
PLOTFLAG = 0; % plots the network data


% here generate the total scenario with all passengers
[T_busFleet, T_Passenger,T_busStop,T_Charger,T_Station,allStationDeps,T_depot] = generateScenario(saveData_UrbanMorph,nTransitStations, Station_separation, nPassengers, Pax_minRadius, Pax_maxRadius,...
  paxSeparation, maxWalkingDist, BS_separation, nCharger, charger_radius, demandPeakness, PLOTFLAG);

% ====== 
% ====== 
% now subsample the T_passenger and write to yumeng style files.
% will need to generate different filenames?

nDays = 30; nSample = 75;

% ====== BUS FLEET
% we will have too many buses since nBus computed on total population
busType = 1; maxPax = 10; maxKWH = 35.8; % kWh - does this correspond to 100% SOC or maxSOC?
minSOC = 10; maxSOC = 80; SOC = 20; consumption = 0.24; %kW.km
nBus = ceil(0.5*nSample*0.7/maxPax);
busFleet1 = table(busType, maxPax, maxKWH, minSOC, maxSOC, SOC, consumption); busFleet1 = repmat(busFleet1,nBus,1);
% distribute SOC amongst buses
busFleet1.SOC = linspace(20,80,nBus)'; % even spacing of bus SOC from 20 -> 80

busType = 2; maxPax = 20; maxKWH = 53.7; % kWh - does this correspond to 100% SOC or maxSOC?
minSoC = 10; maxSoC = 80; SoC = 20; consumption = 0.29; %kW.km
nBus = ceil(0.5*nSample*0.7/maxPax);
busFleet2 = table(busType,maxPax,maxKWH,minSOC,maxSOC,SOC,consumption); busFleet2 = repmat(busFleet2,nBus,1);
busFleet2.SOC = linspace(20,80,nBus)'; % even spacing of bus SOC from 20 -> 80
this_busFleet = [busFleet1;busFleet2];


data_repo = [repo_flexbus, 'data\P100_75_30\'];
pSchedule = zeros(nSample,nDays);
rng('default') % for reproducability
for i = 1:nDays
  this_Pax = sort(randperm(nPassengers,nSample));
  pSchedule(:,i) = this_Pax ;

  T_thisPax = T_Passenger(this_Pax,:);
  % do I need to re-number passengers to be consecutive?

  saveFolder = [data_repo, sprintf('P%d_%dD%d_%03d',nPassengers,nSample,nDays,i)];
  if ~isfolder(saveFolder), mkdir(saveFolder); end
  saveToYumengFormat(saveFolder,this_busFleet, T_thisPax ,T_busStop,T_Charger,T_Station,allStationDeps,T_depot,maxWalkingDist)
end
save([data_repo,'scenarioWorkspace'])

% runFlexbusAnalysis(data_repo)


