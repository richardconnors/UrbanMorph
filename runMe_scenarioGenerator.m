
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
saveData_UrbanMorph = [repo_urbanMorph, 'sampling_data_high\'];

% ====== PARAMETERS for saving
p.maxWalkingDist = 1.0; % BS_separation*sqrt(2)/2; %
p.walkingSpeed = 0.085; % km/minute
p.busSpeed = 0.83; % km/minute
p.busStopServiceTime = 0.5; % minutes
p.maxTransitWaitingTime = 15; % minutes
p.nTotalPop = 50;
p.nPax = 20;
p.Pax_minRadius = 1.5;
p.Pax_maxRadius = 5.0;
p.paxSeparation = 0.05;
p.BS_separation = 1.4;
p.nStation = 1;
p.nCharger = 1; % how many chargers PER TOWN/TRANSIT STATION
p.chargerRadius = 3;
p.demandPeakness = 2;
p.stationSeparation = 3;
p = orderfields(p);

rng('shuffle') 
% nC = 7; cityDiameter = linspace(3,30,nC); 
% nP = 4;  cityCustomers = ceil(linspace(20,100,nP)); % customers to pick up each day
cityDiameter = [3,5,10,15,20,30]; nC = numel(cityDiameter);
%cityCustomers = [20,50,100,150]; nP = numel(cityCustomers);  
cityCustomers = [250,500]; nP = numel(cityCustomers);  
nRuns = 20; % with random shuffle
for j = 1:nP % each population size
  p.nPax = cityCustomers(j); p.nTotalPop = 2*p.nPax;
  for i = 1:nC % each city size
    p.Pax_maxRadius = cityDiameter(i); % max radius around station for passenger locations
    
    %===================================================================
    % Create one population and then sample from it
    %==================================================================
    p.nPax = p.nTotalPop; % to do generation
    [T_busFleet, T_Passenger,T_busStop,T_Charger,T_Station,allStationDeps,T_depot] = generateScenario(p);
    p.nPax = cityCustomers(j);

    % % downsize bus fleet?? - will be different for each instance of i,j loops
    % thisT_busFleet = T_busFleet([],:);
    % nBusTypes = max(T_busFleet.busType);
    % for bb = 1:nBusTypes
    %   ind = find(T_busFleet.busType==bb);
    %   nThisBus = numel(ind); newNumBus = ceil(nThisBus*p.nPax./p.nTotalPop);
    %   thisT_busFleet = [thisT_busFleet; T_busFleet(ind(1:newNumBus),:)]; %#ok<AGROW>
    % end
    thisT_busFleet = T_busFleet;


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
      scenarioSave(saveFolder, p, thisT_busFleet , thisT_Passenger,thisT_busStop,T_Charger,T_Station,allStationDeps,T_depot)
    end

    %===========================================================================
  end
end










return
% 
% switch SCENARIO
%   case 'makeFigure'
%     % here generate the total scenario with all passengers
%     [T_busFleet, T_Passenger,T_busStop,T_Charger,T_Station,allStationDeps,T_depot] = generateScenario(p);
%     this_instance = scenarioName(p);
%     fh = plotScenario(T_busFleet, T_Passenger,T_busStop,T_Charger,T_Station,allStationDeps,T_depot);
%     print(fh, [saveFolder,'\',this_instance_folder,'.jpg'], '-djpeg', '-r300');
%     close(fh);
% 
%   case 'radius_and_population_grid'
%     %===========================================================================
%     % create a sequence of cities and run each one
%     %===========================================================================
%     nC = 8; nP = 7;
%     cityDiameter = linspace(3,40,nC);
%     pop = ceil(linspace(10,150,nP));
%     for j = 1:nP
%       p.nPax = pop(j);
%       for i = 1:nC
%         p.Pax_maxRadius = cityDiameter(i); % max radius around station for passenger locations
%         % here generate the total scenario with all passengers
%         [T_busFleet, T_Passenger,T_busStop,T_Charger,T_Station,allStationDeps,T_depot] = generateScenario(p);
% 
%         this_instance_folder = [scenarioName(p),sprintf('_r%02d%',rr)]; % append name with _r05 for run number
%         saveFolder = [saveData_UrbanMorph, this_instance_folder]; % put this scenario in this folder
%         scenarioSave(saveFolder, p, T_busFleet, T_Passenger,T_busStop,T_Charger,T_Station,allStationDeps,T_depot)
%       end
%     end
%   case 'expanding_radius_fixedPop'
%     %===========================================================================
%     % create a sequence of cities and run each one
%     %===========================================================================
%     nC = 10; cityDiameter = linspace(2,30,nC);
%     p.nPax = 100;
%     for i = 1:nC
%       p.Pax_maxRadius = cityDiameter(i); % max radius around station for passenger locations
%       % here generate the total scenario with all passengers
%       [T_busFleet, T_Passenger,T_busStop,T_Charger,T_Station,allStationDeps,T_depot] = generateScenario(p);
% 
%       this_instance_folder = scenarioName(p);
%       saveFolder = [saveData_UrbanMorph, this_instance_folder]; % put this scenario in this folder
%       scenarioSave(saveFolder, p, T_busFleet, T_Passenger,T_busStop,T_Charger,T_Station,allStationDeps,T_depot)
%     end
% 
%   case 'expanding_radius_fixedDensity'
%     %===========================================================================
%     % create a sequence of cities and run each one
%     %===========================================================================
%     rng('default') % for reproducability
%     nC = 10; cityDiameter = linspace(5,30,nC);
%     baseArea = pi*(10.^2) - pi*(p.Pax_minRadius.^2); ppkm2 = 25/baseArea;
%     for i = 1:nC
%       p.Pax_maxRadius = cityDiameter(i); % max radius around station for passenger locations
%       area = pi*(p.Pax_maxRadius.^2) - pi*(p.Pax_minRadius.^2);
%       p.nPax = ceil(ppkm2*area);
%       % here generate the total scenario with all passengers
%       [T_busFleet, T_Passenger,T_busStop,T_Charger,T_Station,allStationDeps,T_depot] = generateScenario(p);
% 
%       this_instance_folder = scenarioName(p);
%       saveFolder = [saveData_UrbanMorph, this_instance_folder]; % put this scenario in this folder
%       scenarioSave(saveFolder, p, T_busFleet, T_Passenger,T_busStop,T_Charger,T_Station,allStationDeps,T_depot)
%     end
% 
%   case 'mp_density'
%     %===========================================================================
%     % create a sequence of cities and run each one
%     %===========================================================================
%     nC = 5; bs_distance= linspace(0.2,3,nC);
%     for i = 1:nC
%       p.nPax = 100; % number passengers to generate
%       p.Pax_minRadius = 1.5; % min radius away from station for passenger locations
%       p.Pax_maxRadius = 7.5; % max radius around station for passenger locations
%       p.BS_separation = bs_distance(i);
%       p.maxWalkingDist = 1.05*bs_distance(i)*sqrt(2)/2;
%       % here generate the total scenario with all passengers
%       [T_busFleet, T_Passenger,T_busStop,T_Charger,T_Station,allStationDeps,T_depot] = generateScenario(p);
%       height(T_busStop)
%       this_instance_folder = scenarioName(p);
%       saveFolder = [saveData_UrbanMorph, this_instance_folder]; % put this scenario in this folder
%       scenarioSave(saveFolder, p, T_busFleet, T_Passenger,T_busStop,T_Charger,T_Station,allStationDeps,T_depot)
%     end
% 
%   case 'two_cities'
%     %===========================================================================
%     % create a sequence of city separations
%     %===========================================================================
%     nC = 10; cityseparation = linspace(3,25,nC);
%     for i = 1:nC
%       p.stationSeparation = cityseparation(i); % distance between multiple stations
%       % here generate the total scenario with all passengers
%       [T_busFleet, T_Passenger,T_busStop,T_Charger,T_Station,allStationDeps,T_depot] = generateScenario(p);
% 
%       this_instance_folder = scenarioName(p);
%       saveFolder = [saveData_UrbanMorph, this_instance_folder]; % put this scenario in this folder
%       scenarioSave(saveFolder, p, T_busFleet, T_Passenger,T_busStop,T_Charger,T_Station,allStationDeps,T_depot)
%     end
% 
% end
% 
