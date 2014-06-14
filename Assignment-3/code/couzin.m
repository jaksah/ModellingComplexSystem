%% Flocks with nearest neighbours
function [elong, groupdir] = couzin(N, T, w, p, a, rho, gam, plotornot)

%{ 
	N - Number of individuals
	T - Number of timesteps
	w - Weight of influence of preferred direction
	p - Proportion of indiviuals with knowledge of preferred direction
	a - Minimum distance between particles
    rho - Attraction distance
    gam - Standard deviation of random angle
%}

%Set up plot
if plotornot
	fig=figure;
end
dt = 0.2; % time step
s = a; % constant speed
L=20; %L is the size of the domain on which the particles can move
theta = 2;

c=zeros(2,N,T+1); % Position of all individuals at all timesteps
v=zeros(2,N,T+1); % Direction vector of all individuals at all timesteps
d=zeros(2,N,T+1); % Desired direction of all individuals at all timesteps
nInformed = round(N*p);
g=ones(2,nInformed)/sqrt(2); % Preferred direction of informed individuals
cDiff = zeros(2,N,N);
c(:,:,1) = rand(2,N); % Random initial position
v(:,:,1) = rand(2,N); % Random initial direction

elong = zeros(1,T);
groupdir = zeros(2,T);
groupdir(:,1) = mean(d(:,:,1),2);
% T(i,j) gives the angle with the x axis of the direction of motion of the ith
% particle at time j
%T=zeros(N,J+1);
%T(:,1)=2*pi*rand(N,1); %define initial direction of all particles
%fi = zeros(1,J+1);
%fi(1) = abs(1/N*sqrt(sum(cos(T(:,1)))^2+sum(sin(T(:,1)))^2)); % Global alignment


%For all time steps
for t=1:T
	d(:,:,t+1)=d(:,:,t);
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
			repel = -sum(bsxfun(@rdivide,cDiff(:,tooCloseIdx,i),...
                D(i,tooCloseIdx)),2);
			d(:,i,t+1) = d(:,i,t) + repel;
			%d(:,i,t+1) = repel;
		elseif closeEnoughIdx
			norms = arrayfun(@(idx) norm(v(:,idx,t)),...
                1:size(v(:,closeEnoughIdx,t),2));
			d(:,i,t+1) = sum(bsxfun(@rdivide,...
                cDiff(:,closeEnoughIdx(closeEnoughIdx~=i),i),...
                D(i,closeEnoughIdx(closeEnoughIdx~=i))),2) + ...
				sum(bsxfun(@rdivide,v(:,closeEnoughIdx,t),norms),2);
		end
	end
    norms = arrayfun(@(idx) norm(d(:,idx,t+1)), 1:size(d(:,:,t+1),2));
	d(:,:,t+1)=d(:,:,t+1)./[norms; norms]; % Normalize
	
	d(:,1:nInformed,t+1) = d(:,1:nInformed,t+1) + w*g;
    
    % Add a random (gaussian) angle to direction
    randAngles = normrnd(0,0.01,1,size(d(:,:,t+1),2));
    V = atan2(d(2,:,t+1),d(1,:,t+1)) + randAngles;
    d(:,:,t+1) = [cos(V) ; sin(V)]; % d is now normalized
    groupdir(:,t+1) = mean(d(:,:,t+1),2);
	d_angles = atan2(d(2,:,t+1),d(1,:,t+1));
	v_angles = atan2(v(2,:,t),v(1,:,t));
	for i=1:N
		ang_diff = d_angles(i)-v_angles(i);
		if(abs(ang_diff) < theta*dt)
			v(:,i,t+1) = s*d(:,i,t+1);
		else
			new_v = v_angles(i) + sign(ang_diff)*theta*dt;
			v(:,i,t+1) = [cos(new_v); sin(new_v)];
		end
	end
	
	v(:,:,t+1) = s*d(:,:,t+1);
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
        xlim([-10 40])
        ylim([-10 40])
		drawnow
		pause(0.001)
    end
end