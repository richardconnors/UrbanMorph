
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

repo_folder = 'Q:\REPOS\'; % default
[~, hostname] = system('hostname');
if contains(hostname, 'PP0695')
  repo_folder = 'C:\Users\richard.connors\Documents\REPOS\';
end
repo_flexbus = [repo_folder, 'Flexbus3_v0.7\'];
save_Yumengdata_to = [repo_flexbus, 'data\yumengData\'];
% saveData_UrbanMorph = [repo_folder, 'UrbanMorph\data\'];
saveData_UrbanMorph = [repo_flexbus, 'data\urbanMorph\'];

nStations = 2;
Station_separation = 3; % distance between multiple stations
nPassengers = 100; % number passengers to generate
% Passengers will be generated uniformly distributed in an annulus around each station
% minRadius and maxRadius can have single value = used for all stations
% OR can be a row vector size [1 x nTransitStations]
Pax_minRadius = 1.3; % min radius away from station for passenger locations
Pax_maxRadius = 10.0; % max radius around station for passenger locations
paxSeparation = 0.05; % min distance between passengers
BS_separation = 1; % spacing for grid of potential bus stop locations
maxWalkingDist = 1.2; % BS_separation*sqrt(2)/2; %
nCharger  = 4; % how many chargers PER TOWN/TRANSIT STATION.
charger_radius = 3; % chargers located on circle around each town centre
demandPeakness = 0; % 0 = uniform. 1 = peaked in middle
PLOTFLAG = 0; % plots the network data

% ====== PARAMETERS for saving
params.maxWalkingDist = round(maxWalkingDist,5);
params.walkingSpeed = 0.085; % km/minute
params.busSpeed = 0.83; % km/minute
params.busStopServiceTime = 0.5; % minutes
params.maxTransitWaitingTime = 15; % minutes
params.nPax = nPassengers;
params.Pax_minRadius = Pax_minRadius;
params.Pax_maxRadius = Pax_maxRadius;
params.paxSeparation = paxSeparation;
params.BS_separation = BS_separation;
params.nStation = nStations;
params.nCharger = nCharger;
params.chargerRadius = charger_radius;
params.demandPeakness = demandPeakness;
params = orderfields(params);

SCENARIO = 'two_cities'; % 'sample_from_population'  'expanding_radius' 'mp_density' 'two_cities'
SAVE_TO_YUMENG = 0;
switch SCENARIO
  case 'sample_from_population'
    %===========================================================================
    % Create one population and then sample from it
    %===========================================================================
    % here generate the total scenario with all passengers
    [T_busFleet, T_Passenger,T_busStop,T_Charger,T_Station,allStationDeps,T_depot] =...
      generateScenario(saveData_UrbanMorph,nStations, Station_separation, nPassengers, Pax_minRadius, Pax_maxRadius,...
      paxSeparation, maxWalkingDist, BS_separation, nCharger, charger_radius, demandPeakness, PLOTFLAG);

    % now subsample the T_passenger and write to yumeng style files.
    nDays = 30; nSample = 50;
    % ====== BUS FLEET
    % (need to amend bus fleet to be smaller for sampled daily customer numbers)
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


    pSchedule = zeros(nSample,nDays);
    rng('default') % for reproducability
    for i = 1:nDays
      this_Pax = sort(randperm(nPassengers,nSample));
      pSchedule(:,i) = this_Pax ;
      T_thisPax = T_Passenger(this_Pax,:);
      if SAVE_TO_YUMENG
        saveFolder = [save_Yumengdata_to, sprintf('P%d_%dD%d_%03d',nPassengers,nSample,nDays,i)]; %#ok<*UNRCH>
        if ~isfolder(saveFolder), mkdir(saveFolder); end
        saveToYumengFormat(saveFolder,this_busFleet, T_thisPax ,T_busStop,T_Charger,T_Station,allStationDeps,T_depot,maxWalkingDist)
      else
        saveFolder = [saveData_UrbanMorph, sprintf('P%d_%dD%d_%03d',nPassengers,nSample,nDays,i)];
        if ~isfolder(saveFolder), mkdir(saveFolder); end

        writetable(T_busFleet,[saveFolder '\busFleet.csv'],'Delimiter',',');
        writetable(T_Passenger,[saveFolder '\passengerData.csv'],'Delimiter',',')
        writetable(T_busStop,[saveFolder '\busStopXY.csv'],'Delimiter',',')
        writetable(T_Charger,[saveFolder '\chargerXY.csv'],'Delimiter',',')
        writetable(T_Station,[saveFolder '\stationXY.csv'],'Delimiter',',')
        writetable(allStationDeps,[saveFolder '\transitTimetable.csv'],'Delimiter',',');
        writetable(T_depot,[saveFolder '\depotXY.csv'],'Delimiter',',')
        writetable(struct2table(params),[saveFolder '\parameters.csv'],'Delimiter',',');
      end
      save([saveFolder,'\ml_wkspce']) % save all setup parameters etc

    end
    %===========================================================================

  case 'expanding_radius'
    %===========================================================================
    % create a sequence of cities and run each one
    %===========================================================================
    rng('default') % for reproducability
    N = 10; cityDiameter = linspace(2,30,N);
    for i = 1:N
      Pax_maxRadius = cityDiameter(i); % max radius around station for passenger locations
      % here generate the total scenario with all passengers
      [T_busFleet, T_Passenger,T_busStop,T_Charger,T_Station,allStationDeps,T_depot] = ...
        generateScenario(saveData_UrbanMorph,nStations, Station_separation, nPassengers, Pax_minRadius, Pax_maxRadius,...
        paxSeparation, maxWalkingDist, BS_separation, nCharger, charger_radius, demandPeakness, PLOTFLAG);

      this_instance_folder = sprintf('P%d_R%04.1f',nPassengers,Pax_maxRadius);
      if SAVE_TO_YUMENG
        saveFolder = [save_Yumengdata_to, this_instance_folder ]; % put this scenario in this folder
        if ~isfolder(saveFolder), mkdir(saveFolder); end
        saveToYumengFormat(saveFolder,T_busFleet, T_Passenger ,T_busStop,T_Charger,T_Station,allStationDeps,T_depot,maxWalkingDist)
      else
        saveFolder = [saveData_UrbanMorph, this_instance_folder ]; % put this scenario in this folder
        if ~isfolder(saveFolder), mkdir(saveFolder); end

        writetable(T_busFleet,[saveFolder '\busFleet.csv'],'Delimiter',',');
        writetable(T_Passenger,[saveFolder '\passengerData.csv'],'Delimiter',',')
        writetable(T_busStop,[saveFolder '\busStopXY.csv'],'Delimiter',',')
        writetable(T_Charger,[saveFolder '\chargerXY.csv'],'Delimiter',',')
        writetable(T_Station,[saveFolder '\stationXY.csv'],'Delimiter',',')
        writetable(allStationDeps,[saveFolder '\transitTimetable.csv'],'Delimiter',',');
        writetable(T_depot,[saveFolder '\depotXY.csv'],'Delimiter',',')
        writetable(struct2table(params),[saveFolder '\parameters.csv'],'Delimiter',',');
      end
      save([saveFolder,'\ml_wkspce']) % save all setup parameters etc
    end

  
  case 'two_cities'
    %===========================================================================
    % create a sequence of city separations
    %===========================================================================
    rng('default') % for reproducability
    N = 10; cityseparation = linspace(3,25,N);
    for i = 1:N
      Station_separation = cityseparation(i); % distance between multiple stations
      % here generate the total scenario with all passengers
      [T_busFleet, T_Passenger,T_busStop,T_Charger,T_Station,allStationDeps,T_depot,fh] = ...
        generateScenario(saveData_UrbanMorph,nStations, Station_separation, nPassengers, Pax_minRadius, Pax_maxRadius,...
        paxSeparation, maxWalkingDist, BS_separation, nCharger, charger_radius, demandPeakness, PLOTFLAG);

      this_instance_folder = sprintf('P%dR%04.1f_Sep%04.1f',nPassengers,Pax_maxRadius, Station_separation);
      if SAVE_TO_YUMENG
        saveFolder = [save_Yumengdata_to, this_instance_folder]; % put this scenario in this folder
        if ~isfolder(saveFolder), mkdir(saveFolder); end
        saveToYumengFormat(saveFolder,T_busFleet, T_Passenger ,T_busStop,T_Charger,T_Station,allStationDeps,T_depot,maxWalkingDist)
      else
        saveFolder = [saveData_UrbanMorph, this_instance_folder]; % put this scenario in this folder
        if ~isfolder(saveFolder), mkdir(saveFolder); end

        writetable(T_busFleet,[saveFolder '\busFleet.csv'],'Delimiter',',');
        writetable(T_Passenger,[saveFolder '\passengerData.csv'],'Delimiter',',')
        writetable(T_busStop,[saveFolder '\busStopXY.csv'],'Delimiter',',')
        writetable(T_Charger,[saveFolder '\chargerXY.csv'],'Delimiter',',')
        writetable(T_Station,[saveFolder '\stationXY.csv'],'Delimiter',',')
        writetable(allStationDeps,[saveFolder '\transitTimetable.csv'],'Delimiter',',');
        writetable(T_depot,[saveFolder '\depotXY.csv'],'Delimiter',',')
        writetable(struct2table(params),[saveFolder '\parameters.csv'],'Delimiter',',');
        savefig(fh, [saveFolder,'\',this_instance_folder,'.fig']);
        if ~PLOTFLAG, close(fh);  end

      end
      save([saveFolder,'\ml_wkspce']) % save all setup parameters etc
    end

  
  case 'mp_density'
    %===========================================================================
    % create a sequence of cities and run each one
    %===========================================================================
    rng('default') % for reproducability
    maxBS = sqrt(2)*maxWalkingDist;
    N = 5; bs_distance= linspace(0.1,maxBS,N);
    for i = 1:N
      nPassengers = 50; % number passengers to generate
      Pax_minRadius = 1.5; % min radius away from station for passenger locations
      Pax_maxRadius = 10; % max radius around station for passenger locations
      BS_separation = bs_distance(i);
      % here generate the total scenario with all passengers
      [T_busFleet, T_Passenger,T_busStop,T_Charger,T_Station,allStationDeps,T_depot] = ...
        generateScenario(saveData_UrbanMorph,nStations, Station_separation, nPassengers, Pax_minRadius, Pax_maxRadius,...
        paxSeparation, maxWalkingDist, BS_separation, nCharger, charger_radius, demandPeakness, PLOTFLAG);


      if SAVE_TO_YUMENG
        saveFolder = [save_Yumengdata_to, sprintf('P%d_BS%04.1f',nPassengers,BS_separation)]; % put this scenario in this folder
        if ~isfolder(saveFolder), mkdir(saveFolder); end
        saveToYumengFormat(saveFolder,T_busFleet, T_Passenger ,T_busStop,T_Charger,T_Station,allStationDeps,T_depot,maxWalkingDist)
      else
        saveFolder = [saveData_UrbanMorph, sprintf('P%d_BS%04.1f',nPassengers,BS_separation)]; % put this scenario in this folder
        if ~isfolder(saveFolder), mkdir(saveFolder); end

        writetable(T_busFleet,[saveFolder '\busFleet.csv'],'Delimiter',',');
        writetable(T_Passenger,[saveFolder '\passengerData.csv'],'Delimiter',',')
        writetable(T_busStop,[saveFolder '\busStopXY.csv'],'Delimiter',',')
        writetable(T_Charger,[saveFolder '\chargerXY.csv'],'Delimiter',',')
        writetable(T_Station,[saveFolder '\stationXY.csv'],'Delimiter',',')
        writetable(allStationDeps,[saveFolder '\transitTimetable.csv'],'Delimiter',',');
        writetable(T_depot,[saveFolder '\depotXY.csv'],'Delimiter',',')
        writetable(struct2table(params),[saveFolder '\parameters.csv'],'Delimiter',',');
      end
      save([saveFolder,'\ml_wkspce']) % save all setup parameters etc
    end
end




