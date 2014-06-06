%% Flocks with nearest neighbours
%function [elong, groupdir] = couzin(N, T, w, p, a, rho, plotornot)
N = 10;
T = 1000;
w = 0.5;
p = 0.5;
a = 0.01;
rho = 0.1;
plotornot = 1;

%{ 
	N - Number of individuals
	T - Number of timesteps
	w - Weight of influence of preferred direction
	p - Proportion of indiviuals with knowledge of preferred direction
	a - Minimum distance between particles
%}

N = 10;
T = 1000;
w = 0.5;
p = 0.5;
a = 0.01;
rho = 0.1;
plotornot = 1;
%Set up movie

if plotornot
	fig=figure;
end
dt = 0.1; % time step
s = 0.5; % constant speed
L=20; %L is the size of the domain on which the particles can move

c=zeros(2,N,T+1); % Position of all individuals at all timesteps
v=zeros(2,N,T+1); % Direction vector of all individuals at all timesteps
d=zeros(2,N,T+1); % Desired direction of all individuals at all timesteps
g=[zeros(1,round(N*p)); ones(1,round(N*p))]; % Preferred direction of informed individuals
cDiff = zeros(2,N,N);
c(:,:,1) = rand(2,N); % Random initial position
v(:,:,1) = rand(2,N); % Random initial direction

elong = zeros(1,T);
groupdir = zeros(2,T);
% T(i,j) gives the angle with the x axis of the direction of motion of the ith
% particle at time j
%T=zeros(N,J+1);
%T(:,1)=2*pi*rand(N,1); %define initial direction of all particles
%fi = zeros(1,J+1);
%fi(1) = abs(1/N*sqrt(sum(cos(T(:,1)))^2+sum(sin(T(:,1)))^2)); % Global alignment

sigma = 0.01;
% Use normrnd(0, sigma) for random angle

%For all time steps
for t=1:T
	
	D = squareform(pdist(c(:,:,t)'),'tomatrix');
	for i = 1:N
		for j = 1:N
			cDiff(:,i,j) = c(:,i,t)-c(:,j,t);
		end
	end
	for i = 1:N
		% Find particles that are too close
		tooCloseIdx = find(D(i,:)<a);
		tooCloseIdx(tooCloseIdx==i) = []; % Do not count itself
		closeEnoughIdx = find(D(i,:)<rho);
		%closeEnoughIdx(closeEnoughIdx==i) = [];
		% Calculate desired direction
		if tooCloseIdx
			d(:,i,t) = -sum(bsxfun(@rdivide,cDiff(:,tooCloseIdx,i),...
                D(i,tooCloseIdx)),2);
		elseif closeEnoughIdx
			norms = arrayfun(@(idx) norm(v(:,idx,t)),...
                1:size(v(:,closeEnoughIdx,t),2));
			d(:,i,t) = sum(bsxfun(@rdivide,...
                cDiff(:,closeEnoughIdx(closeEnoughIdx~=i),i),...
                D(i,closeEnoughIdx(closeEnoughIdx~=i))),2) + ...
				sum(bsxfun(@rdivide,v(:,closeEnoughIdx,t),norms),2);
			norms = arrayfun(@(idx) norm(d(:,idx,t)), 1:size(d(:,:,t),2));
			d(:,:,t)=d(:,:,t)./[norms; norms]; % Normalize
		end
	end
    
	d(:,1:round(N*p),t) = d(:,1:round(N*p),t) + w*g;
    
    % Add a random (gaussian) angle to direction
    randAngles = normrnd(0,0.01,1,size(d(:,:,t),2));
    V = atan2(d(2,:,t),d(1,:,t)) + randAngles;
    d(:,:,t) = [cos(V) ; sin(V)]; % d is now normalized
    groupdir(:,t) = mean(d(:,:,t),2);
	
	v(:,:,t+1) = s*d(:,:,t);
	c(:,:,t+1) = c(:,:,t) + dt*v(:,:,t+1);
	
	[bb, e] = orientedBoundingBox(c(:,:,t), groupdir(:,t));
	elong(t) = e;
	
	if plotornot
		qscale = 1;
		quiver(c(1,:,t),c(2,:,t),d(1,:,t)*qscale,d(2,:,t)*qscale,...
			'AutoScaleFactor',0.1,'Marker','.','Markersize',10,'Color','k');
		hold on
		plot(bb(1,[1:end 1]), bb(2,[1:end 1]),'r-')
		hold off
		xlabel('X position')
		ylabel('Y position')
        xlim([0 5])
        ylim([0 5])
		drawnow
		pause(0.1)
	end

	
end



%movien = close(movien); %finishes the movie



