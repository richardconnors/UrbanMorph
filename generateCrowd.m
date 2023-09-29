function X = generateCrowd(sideLength, nPax, minDist, nIter)

% a comment from Richard  to overwrite Jos's commemnt

if ~nargin
  sideLength = 1000; % <-- Choose length of square sides
  nPax = 50; % <-- Choose number of points
  minDist = 100;
  nIter = 500;
end

x0 = sideLength/2; y0 = sideLength/2; % <-- Choose center of square


% Generate 2D epanechnikov-distribution
X = [x0,y0] + (sin(asin(2*rand(nPax,2)-1)/3))*sideLength;

XYR = [x0,y0]+[[-1;1;1;-1;-1],[-1;-1;1;1;-1]]*sideLength/2;
XB = interp1((0:4)'*sideLength,XYR,linspace(0,4*sideLength,200));
XB(end,:) = [];

% Repulsion of seeds to avoid them to be too close to each other
nPax = size(X,1);
Xmin = [x0-sideLength/2,y0-sideLength/2];
Xmax = [x0+sideLength/2,y0+sideLength/2];

% Point on boundary
XR = x0+[-1,1,1,-1,-1]*sideLength/2;
YR = y0+[-1,-1,1,1,-1]*sideLength/2;

cla;
hold on
plot(XR,YR,'r-');
h = plot(X(:,1),X(:,2),'b.');
axis equal

d2min = minDist*minDist;
beta = 0.5;
for k = 1:nIter
  XALL = [X; XB];
  DT = delaunayTriangulation(XALL);
  T = DT.ConnectivityList;
  containX = ismember(T,1:nPax);
  b = any(containX,2);
  TX = T(b,:);
  [r,i0] = find(containX(b,:));
  i = mod(i0+(-1:1),3)+1;
  i = TX(r + (i-1)*size(TX,1));
  T = accumarray([i(:,1);i(:,1)],[i(:,2);i(:,3)],[nPax 1],@(x) {x});
  maxd2 = 0;
  R = zeros(nPax,2);
  move = false(nPax,1);
  for i=1:nPax
    Ti = T{i};
    P = X(i,:) - XALL(Ti,:);
    nP2 = sum(P.^2,2);
    if any(nP2<2*d2min)
      move(i) = true;
      move(Ti(Ti<=nPax)) = true;
    end
    maxd2 = maxd2 + mean(nP2);
    b = Ti > nPax;
    nP2(b) = nP2(b)*5; % reduce repulsion from each point of the border
    R(i,:) = sum(P./max((nP2-d2min),1e-3),1);
  end
  if ~any(move)
    break
  end
  if k==1
    v0 = (sideLength*5e-3)/sqrt(maxd2/nPax);
  end
  R = R(move,:);
  v = v0/sqrt(max(sum(R.^2,2)));
  X(move,:) = X(move,:) + v*R;

  % Project back if points falling outside the rectangle
  X = min(max(X,Xmin),Xmax);

  set(h,'XData',X(:,1),'YData',X(:,2));
%   pause(0.01);
end

theta = linspace(0,2*pi,65);
xc = minDist/2*sin(theta);
yc = minDist/2*cos(theta);
% plot circles of diameter dmin around random points
for i=1:nPax
  plot(X(i,1)+xc,X(i,2)+yc,'k');
end
