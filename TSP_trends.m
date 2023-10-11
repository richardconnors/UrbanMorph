% ====== START network geometry
% Total population from which travellers are sampled each day
nPax = 50;

rng('shuffle') % make randomness repeatable
station_XY = [0,0];

% ====== generate passenger pax_XY
maxIter = 1000; % for algorithm ensuring pax min separation
paxSeparation = 0.05;
this_Pax_minRadius = 0.5;
this_Pax_maxRadius = 10.0;
pax_XY = generatePassengers(station_XY(1), station_XY(2), this_Pax_minRadius, this_Pax_maxRadius, nPax, paxSeparation, maxIter);

% distance of pax from the origin
DfromO = sqrt(pax_XY(:,1).^2+pax_XY(:,2).^2);
[~,distOrder] = sort(DfromO);

%======= how many people will use the bus each day?
M = 5; % number of pax sampled from whole population for today
%=======

% how many days/samples/tours will we run?
nEnsemble = 5000;  L = zeros(nEnsemble,1);
TThistory = zeros(nPax,nEnsemble);
for k = 1:nEnsemble
  % randomly sample from the whole population
  idx = randperm(nPax,M); 
  this_XY = [0,0; pax_XY(idx,:)]; %include the station at the origin.
  
  % SOLVE TSP on todays demand
  [p,L(k)] = tspsearch(this_XY,min(M,10));
  % put the station at the start and end of the tour
  start = find(p==1); tour = [p(start:end), p(1:start)];

  % most efficient direction is to do longest leg empty
  if norm(this_XY(tour(end-1),:)) > norm(this_XY(tour(2),:))
    tour = fliplr(tour);
  end
  % order of coordinates and leg-distances
  tour_XY = this_XY(tour,:);
  D = sqrt(sum(diff(tour_XY).^2, 2));
  
  % PLOT the tour here if you want to
%   tspplot(p,this_XY); hold on; 
%   sh = plot(0,0,'gs'); sh.MarkerSize = 10; sh.MarkerFaceColor = 'g'; 
%   drawnow
%   keyboard
  
  % calculate the experienced cost for each individual
  % NOTE arbitrary which direction the tour is driven
  onboard = []; 
  for j=2:M % put people on board as we go through the tour
    onboard = [onboard;idx(tour(j)-1)]; % by population ID
    TThistory(onboard,k) = TThistory(onboard,k) + D(j)';
  end
  % and now with all onboard, go to station to end.
  TThistory(onboard,k) = TThistory(onboard,k) + D(end);
end


TThistory(TThistory==0)=NaN; % dont want to plot the zeros
figure
subplot(1,2,1)
boxplot(TThistory(distOrder,:)','plotstyle','compact')
title(sprintf('Distribution of individual IVD. Mean Tour %.1f',mean(L))); xlabel('Passengers: ordered by distance from station')
subplot(1,2,2)
boxplot(TThistory(distOrder,:)'./repmat(DfromO(distOrder),1,nEnsemble)','plotstyle','compact')
title('Distribution of IVD multiplier'); xlabel('Passengers: ordered by distance from station')
ax = gca; ax.YLim = [0,30];

