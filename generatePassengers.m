function X = generatePassengers(x0, y0, minRadius, maxRadius, nPax, minPaxDist, maxIter)
% around centre [x0,y0] normally distributed within (rough) extent = radius
% (radius = 3*std dev of normal)
exponentialDecayFlag = 1;
% check of inputs mostly for testing
if ~nargin
  x0 = 0; y0 = 0; % length of square sides in m
  minRadius = 500;
  maxRadius = 1500; % roughly max dist away in from xy_centre in m
  nPax = 50; % number of people
  minPaxDist = 25; % minimum distance between people
  maxIter = 1000; % max iterations of algorithm
end

% Generate 2D epanechnikov-distribution
% X = [x0,y0] + (sin(asin(2*rand(nPax,2)-1)/3))*xy_centre;
% X = [x0,y0] + 0.3*maxRadius*randn(nPax,2);
r = sqrt(minRadius^2+(maxRadius^2-minRadius^2)*rand(nPax,1)); % Using square root here ensures distribution uniformity by area

% Exponential decay of radial density
if exponentialDecayFlag
  decay = 1;
  r = -log(exp(-decay*minRadius)-(exp(-decay*minRadius)-exp(-decay*maxRadius)).*rand(nPax,1))/decay;
  % to see the distribtuion of radii you could plot the histogram
  % histogram(r,linspace(minRadius,maxRadius,ceil(nPax/10)))
end 


t = 2*pi*rand(nPax,1);  
X = [x0 + r.*cos(t), y0 + r.*sin(t)];
D = squareform(pdist(X)); D(1:nPax+1:end) = Inf; % set diagonal to be inf

% check if separation between passenger locations violated.
nConflicts = sum(sum(D<minPaxDist,1)>0);
[maxConflicts,i]=max(sum(D<minPaxDist,1));

count = 0;
while nConflicts>0 && count<maxIter
  r = sqrt(minRadius^2+(maxRadius^2-minRadius^2)*rand); 
  if exponentialDecayFlag
    r = -log(exp(-decay*minRadius)-(exp(-decay*minRadius)-exp(-decay*maxRadius)).*rand)/decay;
  end
  t = 2*pi*rand;
  X(i(1),:) = [x0 + r.*cos(t), y0 + r.*sin(t)];
  %   X(i(1),:) = [x0,y0] + 0.3*maxRadius*randn(1,2);
  D = squareform(pdist(X)); D(1:nPax+1:end) = Inf; % set diagonal to be inf
  nConflictsEach = sum(D<minPaxDist,1);
  nConflicts = sum(nConflictsEach>0);
  [maxConflicts,i]=max(nConflictsEach);
  count = count + 1;
end

DISPLAY = 1;

if DISPLAY
  fprintf('GeneratePax: [Iterations,Conflicts] = [%d,%d] \n', count, nConflicts);
%   cla; scatter(X(:,1),X(:,2)); axis square
  
end




 