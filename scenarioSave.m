function scenarioSave(saveFolder, params, T_busFleet, T_Passenger,T_busStop,T_Charger,T_Station,allStationDeps,T_depot)

this_instance_folder = scenarioName(params);

if ~isfolder(saveFolder), mkdir(saveFolder); end

writetable(T_busFleet,[saveFolder '\busFleet.csv'],'Delimiter',',');
writetable(T_Passenger,[saveFolder '\passengerData.csv'],'Delimiter',',')
writetable(T_busStop,[saveFolder '\busStopXY.csv'],'Delimiter',',')
writetable(T_Charger,[saveFolder '\chargerXY.csv'],'Delimiter',',')
writetable(T_Station,[saveFolder '\stationXY.csv'],'Delimiter',',')
writetable(allStationDeps,[saveFolder '\transitTimetable.csv'],'Delimiter',',');
writetable(T_depot,[saveFolder '\depotXY.csv'],'Delimiter',',')
writetable(struct2table(params),[saveFolder '\parameters.csv'],'Delimiter',',');

fh = plotScenario(T_Passenger,T_busStop,T_Charger,T_Station,T_depot);
fh.Visible = 'off';
print(fh, [saveFolder,'\',this_instance_folder,'.jpg'], '-djpeg', '-r300');
close(fh);

% variables = evalin('base', 'who');
% non_graphics_variables = variables(~startsWith({variables.class},'matlab.ui'));
% % evalin('base', ['save(''' filename ''')']);
% variable_names = {non_graphics_variables.name};
% save([saveFolder,'\ml_wkspce'],variable_names{:}) % save all setup parameters etc


% if SAVE_TO_YUMENG
%   saveFolder = [save_Yumengdata_to, this_instance_folder ]; % put this scenario in this folder
%   if ~isfolder(saveFolder), mkdir(saveFolder); end
%   saveToYumengFormat(saveFolder,T_busFleet, T_Passenger ,T_busStop,T_Charger,T_Station,allStationDeps,T_depot,maxWalkingDist)
% end
