
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
% save_Yumengdata_to = [repo_flexbus, 'data\yumengData\'];
% saveData_UrbanMorph = [repo_flexbus, 'data\urbanMorph_toRun\'];
saveData_UrbanMorph = [repo_urbanMorph, 'data\'];

% ====== PARAMETERS for saving
p.maxWalkingDist = 1.0; % BS_separation*sqrt(2)/2; %
p.walkingSpeed = 0.085; % km/minute
p.busSpeed = 0.83; % km/minute
p.busStopServiceTime = 0.5; % minutes
p.maxTransitWaitingTime = 15; % minutes
p.nPax = 100;
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

PLOTFLAG = 1; % plots the network data
SCENARIO = 'radius_and_population_grid'; % 'sample_from_population'  'expanding_radius' 'mp_density' 'two_cities'
SAVE_TO_YUMENG = 0;
SAVE_TO_URBANMORPH = 1;

rng('shuffle') % for reproducability
nRuns = 10; % with random shuffle
for rr = 1:nRuns
  switch SCENARIO
    case 'makeFigure'
      % here generate the total scenario with all passengers
      [T_busFleet, T_Passenger,T_busStop,T_Charger,T_Station,allStationDeps,T_depot] = generateScenario(p);
      this_instance = scenarioName(p);
      fh = plotScenario(T_busFleet, T_Passenger,T_busStop,T_Charger,T_Station,allStationDeps,T_depot);
      print(fh, [saveFolder,'\',this_instance_folder,'.jpg'], '-djpeg', '-r300');
      close(fh);
    case 'sample_from_population'
      %===========================================================================
      % Create one population and then sample from it
      %===========================================================================
      % here generate the total scenario with all passengers
      [T_busFleet, T_Passenger,T_busStop,T_Charger,T_Station,allStationDeps,T_depot] = generateScenario(p);

      % now subsample the T_passenger
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
        end
        if SAVE_TO_URBANMORPH
          saveFolder = [saveData_UrbanMorph, sprintf('P%d_%dD%d_%03d',nPassengers,nSample,nDays,i)];
          if ~isfolder(saveFolder), mkdir(saveFolder); end

          writetable(this_busFleet,[saveFolder '\busFleet.csv'],'Delimiter',',');
          writetable(T_thisPax,[saveFolder '\passengerData.csv'],'Delimiter',',')
          writetable(T_busStop,[saveFolder '\busStopXY.csv'],'Delimiter',',')
          writetable(T_Charger,[saveFolder '\chargerXY.csv'],'Delimiter',',')
          writetable(T_Station,[saveFolder '\stationXY.csv'],'Delimiter',',')
          writetable(allStationDeps,[saveFolder '\transitTimetable.csv'],'Delimiter',',');
          writetable(T_depot,[saveFolder '\depotXY.csv'],'Delimiter',',')
          writetable(struct2table(p),[saveFolder '\parameters.csv'],'Delimiter',',');
        end
        if SAVE_TO_YUMENG || SAVE_TO_URBANMORPH
          variables = whos;
          non_graphics_variables = variables(~startsWith({variables.class},'matlab.ui'));
          variable_names = {non_graphics_variables.name};
          save([saveFolder,'\ml_wkspce'],variable_names{:}) % save all setup parameters etc
        end

      end
      %===========================================================================

    case 'radius_and_population_grid'
      %===========================================================================
      % create a sequence of cities and run each one
      %===========================================================================
      nC = 8; nP = 7;
      cityDiameter = linspace(3,40,nC);
      pop = ceil(linspace(10,150,nP));
      for j = 1:nP
        p.nPax = pop(j);
        for i = 1:nC
          p.Pax_maxRadius = cityDiameter(i); % max radius around station for passenger locations
          % here generate the total scenario with all passengers
          [T_busFleet, T_Passenger,T_busStop,T_Charger,T_Station,allStationDeps,T_depot] = generateScenario(p);

          this_instance_folder = [scenarioName(p),sprintf('_r%02d%',rr)]; % append name with _r05 for run number
          saveFolder = [saveData_UrbanMorph, this_instance_folder]; % put this scenario in this folder
          scenarioSave(saveFolder, p, T_busFleet, T_Passenger,T_busStop,T_Charger,T_Station,allStationDeps,T_depot)
        end
      end
    case 'expanding_radius_fixedPop'
      %===========================================================================
      % create a sequence of cities and run each one
      %===========================================================================
      nC = 10; cityDiameter = linspace(2,30,nC);
      p.nPax = 100;
      for i = 1:nC
        p.Pax_maxRadius = cityDiameter(i); % max radius around station for passenger locations
        % here generate the total scenario with all passengers
        [T_busFleet, T_Passenger,T_busStop,T_Charger,T_Station,allStationDeps,T_depot] = generateScenario(p);

        this_instance_folder = scenarioName(p);
        saveFolder = [saveData_UrbanMorph, this_instance_folder]; % put this scenario in this folder
        scenarioSave(saveFolder, p, T_busFleet, T_Passenger,T_busStop,T_Charger,T_Station,allStationDeps,T_depot)
      end

    case 'expanding_radius_fixedDensity'
      %===========================================================================
      % create a sequence of cities and run each one
      %===========================================================================
      rng('default') % for reproducability
      nC = 10; cityDiameter = linspace(5,30,nC);
      baseArea = pi*(10.^2) - pi*(p.Pax_minRadius.^2); ppkm2 = 25/baseArea;
      for i = 1:nC
        p.Pax_maxRadius = cityDiameter(i); % max radius around station for passenger locations
        area = pi*(p.Pax_maxRadius.^2) - pi*(p.Pax_minRadius.^2);
        p.nPax = ceil(ppkm2*area);
        % here generate the total scenario with all passengers
        [T_busFleet, T_Passenger,T_busStop,T_Charger,T_Station,allStationDeps,T_depot] = generateScenario(p);

        this_instance_folder = scenarioName(p);
        saveFolder = [saveData_UrbanMorph, this_instance_folder]; % put this scenario in this folder
        scenarioSave(saveFolder, p, T_busFleet, T_Passenger,T_busStop,T_Charger,T_Station,allStationDeps,T_depot)
      end

    case 'mp_density'
      %===========================================================================
      % create a sequence of cities and run each one
      %===========================================================================
      nC = 5; bs_distance= linspace(0.2,3,nC);
      for i = 1:nC
        p.nPax = 100; % number passengers to generate
        p.Pax_minRadius = 1.5; % min radius away from station for passenger locations
        p.Pax_maxRadius = 7.5; % max radius around station for passenger locations
        p.BS_separation = bs_distance(i);
        p.maxWalkingDist = 1.05*bs_distance(i)*sqrt(2)/2;
        % here generate the total scenario with all passengers
        [T_busFleet, T_Passenger,T_busStop,T_Charger,T_Station,allStationDeps,T_depot] = generateScenario(p);
        height(T_busStop)
        this_instance_folder = scenarioName(p);
        saveFolder = [saveData_UrbanMorph, this_instance_folder]; % put this scenario in this folder
        scenarioSave(saveFolder, p, T_busFleet, T_Passenger,T_busStop,T_Charger,T_Station,allStationDeps,T_depot)
      end

    case 'two_cities'
      %===========================================================================
      % create a sequence of city separations
      %===========================================================================
      nC = 10; cityseparation = linspace(3,25,nC);
      for i = 1:nC
        p.stationSeparation = cityseparation(i); % distance between multiple stations
        % here generate the total scenario with all passengers
        [T_busFleet, T_Passenger,T_busStop,T_Charger,T_Station,allStationDeps,T_depot] = generateScenario(p);

        this_instance_folder = scenarioName(p);
        saveFolder = [saveData_UrbanMorph, this_instance_folder]; % put this scenario in this folder
        scenarioSave(saveFolder, p, T_busFleet, T_Passenger,T_busStop,T_Charger,T_Station,allStationDeps,T_depot)
      end

  end


end

