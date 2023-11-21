
% generateScenario allows the inputs listed below to be specified.
% You can set values here, and run this script to call the function
% Outputs are written to comma delimited text files (easy to change this)
% Output text file names currently hardcoded as below, saved to PWD
% transitTimetable.txt
% stationXY.txt
% chargerXY.txt
% paxData.txt
% busStopXY.txt
% ========================================================================
% All distances in KM
% Timetable of departures for each station is hardcoded
% For each passenger departure time is randomly assigned
% Demand distribution is controlled and can go from uniform -> peaked
% nTransitStations = 2;
% Station_separation = 5; % distance between multiple stations
% nPassengers = 100; % number passengers to generate
% Passengers will be generated uniformly distributed in an annulus around each station
% minRadius and maxRadius can have single value = used for all stations
% OR can be a row vector size [1 x nTransitStations]
% Pax_minRadius = 0.5; % min radius away from station for passenger locations
% Pax_maxRadius = 2.0; % max radius around station for passenger locations
% paxSeparation = 0.05; % min distance between passengers
% BS_separation = 1; % spacing for grid of potential bus stop locations
% maxWalkingDist = BS_separation*1.05/sqrt(2); %
% nCharger  = 4; % how many chargers PER TOWN/TRANSIT STATION.
% charger_radius = 1.5; % chargers located on circle around each town centre
% demandPeakness = 0; % 0 = uniform. 1 = peaked in middle
% PLOTFLAG = 1; % plots the network data

repo_flexbus = [get_repo_folder, 'Flexbus3_v0.8.4\'];
repo_urbanMorph = [get_repo_folder, 'UrbanMorph\'];
saveData_UrbanMorph = [repo_urbanMorph, 'data_bicentric\'];

% ====== PARAMETERS for saving
p.maxWalkingDist = 1.0; % BS_separation*sqrt(2)/2; %
p.walkingSpeed = 0.085; % km/minute
p.busSpeed = 0.83; % km/minute
p.busStopServiceTime = 0.5; % minutes
p.maxTransitWaitingTime = 15; % minutes
p.nTotalPop = 200;
p.nPax = 100;
p.Pax_minRadius = 1.5;
p.Pax_maxRadius = 5.0;
p.paxSeparation = 0.05;
p.BS_separation = 1.4;
p.nStation = 2;
p.nCharger = 1; % how many chargers PER TOWN/TRANSIT STATION
p.chargerRadius = 3;
p.demandPeakness = 2;
p.stationSeparation = 3;
p = orderfields(p);

rng('shuffle')
cityDiameter = [5,15]; nC = numel(cityDiameter);
cityGap= [1,2,5,10]; nGap = numel(cityGap);

p.nPax = 100; p.nTotalPop = 200;
nRuns = 20; % with random shuffle
for i = 1:nC % each city size
  p.Pax_maxRadius = cityDiameter(i); % max radius around station for passenger locations
  for j = 1:nGap % each gap size
    p.stationSeparation = cityGap(j);

    %===================================================================
    % Create one population and then sample from it
    %==================================================================
    p.nPax = 200; % to do generation
    [T_busFleet, T_Passenger,T_busStop,T_Charger,T_Station,allStationDeps,T_depot] = generateScenario(p);
    p.nPax = 100;
    for rr = 1:nRuns
      % to ensure always get nSample passengers
      thisPax = randperm(p.nTotalPop); % reorder all pax ids
      thisPax = sort(thisPax(1:p.nPax)); % take first few ids and re-sort
      thisT_Passenger = T_Passenger(thisPax,:);


      % calculate bus stops we need to use for this case
      BS_XY = [T_busStop.busStop_X,T_busStop.busStop_Y];
      pax_XY = [thisT_Passenger.passenger_X,thisT_Passenger.passenger_Y];
      % remove bus stops more than walking distance away from any pax
      [~,dist] = dsearchn(pax_XY,BS_XY);
      BSinReach = dist<=p.maxWalkingDist;
      nBS = sum(BSinReach);
      busStop_ID = (1:nBS)';
      busStop_X = BS_XY(BSinReach,1);
      busStop_Y = BS_XY(BSinReach,2);
      thisT_busStop = table(busStop_ID,busStop_X,busStop_Y);

      this_instance_folder = [scenarioName(p),sprintf('_r%02d%',rr)]; % append name with _r05 for run number
      saveFolder = [saveData_UrbanMorph, this_instance_folder]; % put this scenario in this folder
      scenarioSave(saveFolder, p, T_busFleet , thisT_Passenger,thisT_busStop,T_Charger,T_Station,allStationDeps,T_depot)
    end

    %===========================================================================
  end
end
