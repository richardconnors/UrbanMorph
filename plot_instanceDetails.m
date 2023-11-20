%
%
% LOAD DATA
% repo_flexbus = [get_repo_folder, 'Flexbus3_v0.8.4\'];
% data_to_load = [repo_flexbus, 'results\grid_results_v2\'];
% data_to_load = 'Q:\REPOS\UrbanMorph\data\';

baseFolderName = 'P050_S1_R30.0_W1.0_MP01.4_SS03.0_DP2.0_r';
for i=1:20
  thisFolderName = [baseFolderName, sprintf('%02d',i)];

  data_to_load = [get_repo_folder, 'UrbanMorph\sampling_data\'];
  thisFolder = [data_to_load, thisFolderName,filesep];
  folderItems = dir(thisFolder);
  allFiles = string({folderItems.name});

  % load all the usual scenario tables and params
  params_file = fullfile(thisFolder, allFiles(contains(allFiles,'parameters'))); p = readtable(params_file);
  filename = fullfile(thisFolder, allFiles(contains(allFiles,'passenger'))); T_Pax = readtable(filename);
  filename = fullfile(thisFolder, allFiles(contains(allFiles,'busStopXY'))); T_busStop = readtable(filename);
  filename = fullfile(thisFolder, allFiles(contains(allFiles,'chargerXY'))); T_Charger = readtable(filename);
  filename = fullfile(thisFolder, allFiles(contains(allFiles,'stationXY'))); T_Station = readtable(filename);
  filename = fullfile(thisFolder, allFiles(contains(allFiles,'depotXY'))); T_depot = readtable(filename);
  filename = fullfile(thisFolder, allFiles(contains(allFiles,'busFleet'))); T_busFleet = readtable(filename);

  fh = plotScenario(T_Pax,T_busStop,T_Charger,T_Station,T_depot); hold on;

  cus_file = allFiles(contains(allFiles,'cus_'));
  T_C = readtable(fullfile(thisFolder, cus_file));
  T_Pax = [T_Pax, T_C]; paxOK = T_Pax.ride_time>0; % these are the served customers

  route_files =  allFiles(startsWith(allFiles,'route_detail_'));
  this_fleetSize = numel(route_files);
  % new instance so we zero all the totals
  total_veh_kms = 0; total_empty_kms = 0; total_cus_kms = 0; total_direct_kms = 0; allCus = [];

  cmap = colormap_generator(this_fleetSize);
  for ff = 1:this_fleetSize % one for each bus. Could do multiple loops
    R = readtable(fullfile(thisFolder,route_files{ff}));  % You can use readmatrix for newer MATLAB versions
    this_leg_kms = sqrt(diff(R.x).^2 + diff(R.y).^2);
    this_leg_occ = R.num_passenger(1:end-1);
    this_cus_kms = this_leg_kms'*this_leg_occ; % empty bus automatically not counted
    this_empty_kms = this_leg_kms'*~this_leg_occ; % empty bus automatically not counted

    C = string(R.cus_set); C = C(~startsWith(C,"Int64")); % extract customer numbers on this bus
    this_bus_cus = cell2mat(cellfun(@str2double, regexp(C(:)', '\d+', 'match'), 'UniformOutput', false));
    allCus = [allCus;this_bus_cus(:)];
    this_direct_kms = sum(T_Pax.DistanceFromStation(ismember(T_Pax.cus, this_bus_cus)));
    total_veh_kms = total_veh_kms + sum(this_leg_kms);
    total_empty_kms = total_empty_kms + this_empty_kms;
    total_cus_kms = total_cus_kms + this_cus_kms;
    total_direct_kms = total_direct_kms + this_direct_kms;
    % R_bus(end+1,:) = {i,p.nPax,p.Pax_maxRadius,p.BS_separation,p.maxWalkingDist,sum(this_leg_kms),this_empty_kms,this_direct_kms,this_cus_kms,max(this_leg_occ)}; %#ok<*SAGROW>

    figure(fh)
    % lh = line(R.x,R.y); lh.Color = cmap(ff,:); lh.HandleVisibility = 'off';
    hold on;
    % quiver(R.x(1:end-1), R.y(1:end-1), diff(R.x), diff(R.y), 0, 'MaxHeadSize', 0.5, 'Color', 'red');
    axis('manual')
    ah = arrow([R.x(1:end-1),R.y(1:end-1)],[R.x(2:end),R.y(2:end)],'EdgeColor',cmap(ff,:),'FaceColor',cmap(ff,:));
    for hhh=1:numel(ah), ah(hhh).HandleVisibility = 'off';  end



  end % loop bus routes
  % allCus 


  % do we want text numbers?
  th = text(T_Pax.passenger_X,T_Pax.passenger_Y, string(T_Pax.passenger_ID), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');

end

