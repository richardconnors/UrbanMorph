function saveToYumengFormat(saveFolder,T_busFleet, T_Passenger,T_busStop,T_Charger,T_Station,allStationDeps,T_depot,maxWalkingDist)
%=========================================================================
% % ====== SAVE TO YUMENG FORMAT 2 FILES

nPax = size(T_Passenger,1);
nBS = max(T_busStop.busStop_ID);
nStation = max(allStationDeps.StationID);
nCharger = size(T_Charger,1);
busType = max(T_busFleet.busType);

yumeng_filename = sprintf('c-%d-bs-%d.txt',nPax,nBS);
yumeng_file = [saveFolder,'\',yumeng_filename];

tService = 0.5; busSpeed = 50/60; walkSpeed = 5.1/60;
fileID = fopen(yumeng_file, 'w'); % 'w' means overwrite any existing content
fprintf(fileID, '%d %d %d %d %.1f %.1f %.2f %.2f', nPax, nBS, nStation, busType, maxWalkingDist, tService, busSpeed, walkSpeed);
for k = 1:busType
  thisBus = T_busFleet(T_busFleet.busType==k,:);
  nThisBus = sum(T_busFleet.busType==k);
  thismin = thisBus.maxKWH(1).*thisBus.minSOC(1)./100; % use 1st row in subtable
  thismax = thisBus.maxKWH(1).*thisBus.maxSOC(1)./100;
  fprintf(fileID, '\n');
  fprintf(fileID, '%d %d %.2f %.2f %.2f', nThisBus, thisBus.maxPax(1), thisBus.consumption(1), thismax, thismin);
end
fprintf(fileID, '\n');
nDummies = 3; % something used for exact solution
fprintf(fileID, '%d %.2f %d',nCharger, T_Charger.charger_rate(1), nDummies); 
fprintf(fileID, '\n');
fprintf(fileID, '0 %.1f %.1f',T_depot.depot_X,T_depot.depot_Y); % what is 3 indicating?????
fclose(fileID);

fileID = fopen(yumeng_file, 'a'); % 'w' means overwrite any existing content
Y_Passenger = T_Passenger(:,[1,2,3,5]);
writetable(Y_Passenger, yumeng_file, 'Delimiter', ' ', 'WriteVariableNames', false,'WriteMode','append');
Y_busStop = T_busStop; Y_busStop.busStop_ID = Y_busStop.busStop_ID+nPax;
writetable(Y_busStop, yumeng_file, 'Delimiter', ' ', 'WriteVariableNames', false, 'WriteMode', 'append');

nRows = Y_busStop.busStop_ID(end); 
rowNums = table((nRows+1:nRows+nStation)');
Y_Station = [rowNums,T_Station(:,[2,3])];
writetable(Y_Station, yumeng_file, 'Delimiter', ' ', 'WriteVariableNames', false, 'WriteMode', 'append');

nRows = nRows+1+nStation;
rowNums = table((nRows+1:nRows+nCharger)');
Y_Charger = [rowNums,T_Charger(:,[2,3])];
writetable(Y_Charger, yumeng_file, 'Delimiter', ' ', 'WriteVariableNames', false, 'WriteMode', 'append');
fclose(fileID);

YumengTimetable = allStationDeps(:,[1,3,4,5]);
YumengTimetable.Properties.VariableNames = {'No',	'E',	'L',	'no_layer'};

yumeng_timetable_csv = [saveFolder,'\', sprintf('c-%d-bs-%d.csv',nPax,nBS)];
writetable(YumengTimetable,yumeng_timetable_csv,'Delimiter',',');

% ====== end SAVE TO YUMENG
%=========================================================================