% ====== START network geometry
% Total population from which travellers are sampled each day
nPax = 1000;

rng default % make randomness repeatable
station_XY = [0,0];

% ====== generate passenger pax_XY
maxIter = 1000; % for algorithm ensuring pax min separation
paxSeparation = 0.05;
this_Pax_minRadius = 0.5;
this_Pax_maxRadius = 10.0;
pax_XY = generatePassengers(station_XY(1), station_XY(2), this_Pax_minRadius, this_Pax_maxRadius, nPax, paxSeparation, maxIter);

% how many people will use the bus each day?
N = 10; % number of pax sampled from whol population for today
% how many days/samples/tours will we run?
nEnsemble = 1000;

TThistory = zeros(nPax,nEnsemble);
for k = 1:nEnsemble
  % randomly sample from the whole population
  idx = randperm(nPax,N); 
  this_XY = [0,0; pax_XY(idx,:)]; %include the station at the origin.
  
  % SOLVE TSP on todays demand
  [p,L] = tspsearch(this_XY,10);
  % put the station at the start and end of the tour
  start = find(p==1); tour = [p(start:end), p(1:start)];
  
  % order of coordinates and leg-distances
  tour_XY = this_XY(tour,:);
  D = sqrt(sum(diff(tour_XY).^2, 2));
  
  % PLOT the tour here if you want to
  % tspplot(p,this_XY); hold on; plot(0,0,'rs')
  
  % calculate the experienced cost for each individual
  % NOTE arbitrary which direction the tour is driven
  onboard = []; 
  for j=2:N % put people on board as we go through the tour
    onboard = [onboard;idx(tour(j)-1)]; % by population ID
    TThistory(onboard,k) = TThistory(onboard,k) + D(j)';
  end
  % and now with all onboard, go to station to end.
  TThistory(onboard,k) = TThistory(onboard,k) + D(end);
end

% reorder pax by distance from the origin for plotting
DfromO = sqrt(pax_XY(:,1).^2+pax_XY(:,2).^2);
[~,distOrder] = sort(DfromO);

TThistory(TThistory==0)=NaN; % dont want to plot the zeros
boxplot(TThistory(distOrder,:)','plotstyle','compact')
