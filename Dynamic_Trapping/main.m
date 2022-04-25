clear; close all; clc;

%% Set-up
% System parameters (SI units)
p.w = 0.06; % Diameter of pheromone production
p.v0 = 0.04; % RAnt speed
p.Ls = 0.01; % Sensor distance RAnts
p.kp = 0.1; % pheromone production rate
p.km = 2*p.kp; % pheromone decay rate
p.G = 0.1;

p.nr = 5; % number of RAnts
p.sr = 0.04; % RAnt size
p.tE = 10; % terminal time

% Problem parameters
p.L = 0.5; % size of square-shaped arena

% numerical parameters
p.nx = 210; % Grid point number in x
p.ny = p.nx; % Grid point number in y
p.dt = 1e-2; % Time step
p.cx = repmat(linspace(0,p.L,p.nx),[p.ny 1]); % x location of pheromone matrix entries
p.cy = repmat(linspace(p.L,0,p.ny)',[1 p.nx]); % y location of pheromone matrix entries

% Initial conditions
rng(2); % Randiom initial seed
x0 = randInit(p);

%% Integration
[T,X] = euler(@(t,x) dynamics(t,x,p),[0 p.tE],x0,p);

%% Plot results
repSpeed = 1; % Replay speed
p.vidFlag = 0; % set to 1 if video is to be recorded
p.vidDir = '../video/trapping.avi';
animation(T,X,repSpeed,p) % animation

%% Sub-functions
% Coupled RAnts and photormone dynamics
function dx = dynamics(~,x,p)
% unpack
c = x.c;

% RAnt dynamics
dx = rantRules(x,p);

% Photormone dynamics
cMask = findCMask(x,p);
dx.c = p.kp*cMask - p.km*c;
end

% RAnt local rules
function dx = rantRules(x,p)
% unpack
r = x.r; c = x.c; th = x.th;

pv = [cos(th),sin(th)];
cDiff = findDiff(x,p); % find difference in sensor signals

dx.r = p.v0*pv;
dx.th = p.G*cDiff/p.Ls;
end

% Find difference in left and right sensor signals of RAnt
function cDiff = findDiff(x,p)
% unpack
c = x.c; r = x.r; th = x.th;

n = [-sin(th),cos(th)]; % normal vector to RAnt orientation

rL = r + n*p.Ls/2; % left sensor position
rR = r - n*p.Ls/2; % right sensor position

cL = zeros(p.nr,1); cR = zeros(p.nr,1);
for i=1:p.nr
    cL(i) = getVal(rL(i,:),c,p);
    cR(i) = getVal(rR(i,:),c,p);
end

cDiff = cL-cR;
end

% get interpolated value of photormone field c at location r
function cr = getVal(r,c,p)
it = [p.ny-r(2)/p.L*p.ny,r(1)/p.L*p.nx];
flv = [mod(p.ny-floor(r(2)/p.L*p.ny)-2,p.ny)+1,mod(floor(r(1)/p.L*p.nx)-1,p.nx)+1];
clv = [mod(flv(1),p.ny)+1,mod(flv(2),p.nx)+1];
A = [c(flv(1),flv(2)),c(flv(1),clv(2));
    c(clv(1),flv(2)),c(clv(1),clv(2))];
wty = [1-(it(1)-floor(it(1))),it(1)-floor(it(1))];
wtx = [1-(it(2)-floor(it(2)));it(2)-floor(it(2))];
cr = wty*A*wtx;
end

% Find location of photormone production
function cMask = findCMask(x,p)
% unpack
c = x.c; r = x.r;

cMask = zeros(p.ny*p.nx,1);
% loop over all RAnts
for i=1:p.nr
    % apply periodic boundary condition
    cx = p.cx;
    cy = p.cy;
    
    % find distance matrix
    dist = sqrt((cx-r(i,1)).^2+(cy-r(i,2)).^2);
    it = find(dist<p.w/2);
    dpx = sqrt(2)*p.L/p.nx;
    itc = find(dist>p.w/2 & dist<p.w/2+dpx);
    
    % update photormone production mask
    cMask(it) = cMask(it) + 1;
    cMask(itc) = cMask(itc) + (p.w/2+dpx-dist(itc))/dpx; % pixel correction
end

cMask = reshape(cMask,[p.ny p.nx]);
end

% Compute random initial position of RAnts
function x0 = randInit(p)
% Find solution for non-overlapping RAnt position
D = 1;
while D>0
    r = p.sr+rand(p.nr,2)*(p.L-2*p.sr);        
    dist = eye(p.nr)+sqrt((r(:,1)-r(:,1)').^2+(r(:,2)-r(:,2)').^2);
    D = sum(sum(triu(dist<p.sr)));
end
x0.r = r;
x0.th = 2*pi*rand(p.nr,1);
x0.c = zeros(p.ny,p.nx); % initial pheromone concentration
end

% compute contact dynamics. Rants cannot penetrate each other nor boundary
function [x,dx] = contact(x,dx,p)
r = x.r;
dr = dx.r;

rM = r+dr*p.dt/2; % midpoint displacement

% Rant distances
dist = eye(p.nr)+sqrt((rM(:,1)-rM(:,1)').^2+(rM(:,2)-rM(:,2)').^2);
lDist = triu(dist<p.sr);

for i=1:p.nr
    % Rant-Rant interactions
    idx = find(lDist(i,:)==1);
    for j=idx
        n = r(j,:)-r(i,:); n = n/norm(n); t = [-n(2);n(1)];
        if dr(i,:)*n'>0
            dr(i,:) = (dr(i,:)*t)*t';
        end
        if dr(j,:)*n'<0
            dr(j,:) = (dr(j,:)*t)*t';
        end
    end
    
    % Rant-boundary interactions
    if rM(i,1) < p.sr
        t = [0;1];
        dr(i,:) = (dr(i,:)*t)*t';
    elseif rM(i,1) > p.L-p.sr
        t = [0;1];
        dr(i,:) = (dr(i,:)*t)*t';
    end
    
    if rM(i,2) < p.sr
        t = [1;0];
        dr(i,:) = (dr(i,:)*t)*t';
    elseif rM(i,2) > p.L-p.sr
        t = [1;0];
        dr(i,:) = (dr(i,:)*t)*t';
    end
end

dx.r = dr;
end

% Explicit first order Euler integration scheme
function [T,X] = euler(f,tSpan,x0,p)
N = round((tSpan(2)-tSpan(1))/p.dt);
T = linspace(tSpan(1),tSpan(2),N);

X.r = zeros(N,2*p.nr); X.th = zeros(N,p.nr); X.c = zeros(N,p.nx*p.ny);
X.r(1,:) = x0.r(:); X.th(1,:) = x0.th(:); X.c(1,:) = x0.c(:);
for i=1:N-1
    x.r = reshape(X.r(i,:),[p.nr 2]);
    x.th = reshape(X.th(i,:),[p.nr 1]);
    x.c = reshape(X.c(i,:),[p.ny p.nx]);
    
    [x,dx] = contact(x,f(0,x),p); % apply collision constraints
    rNew = x.r+dx.r*p.dt;
    X.r(i+1,:) = rNew(:);
    X.th(i+1,:) = x.th(:) + dx.th(:)*p.dt;
    X.c(i+1,:) = x.c(:) + dx.c(:)*p.dt;
    outputFcn(T(i),0,0,tSpan);
end
end

% display integration progress
function status = outputFcn(t,~,~,tSpan)
status = 0;
if numel(t)==1
    t0 = tSpan(1); tF = tSpan(end);
    pct = ceil((t(1)-t0)/(tF-t0)*100);
    strOut = ['Progress: ',num2str(pct),'/100'];
    strCR = repmat('\b',1,length(strOut));
    fprintf([strCR,strOut]);
end
end