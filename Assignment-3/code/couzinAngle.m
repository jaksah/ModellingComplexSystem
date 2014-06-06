%% Flocks with nearest neighbours
function [elong, groupdir] = couzin(N, T, w, prop, a, rho, gam, plotornot)
%{ 
	N - Number of individuals
	T - Number of timesteps
	w - Weight of influence of preferred direction
	p - Proportion of indiviuals with knowledge of preferred direction
	a - Minimum distance between particles
    rho - Attraction distance

N=100;
T=1000;
w=0.25;
prop=0.15;
a = 0.5;
rho = 6;
gam=0.0;
plotornot=1;
%}
if plotornot
	fig=figure;
end

dt = 0.2; % time step
s = a; % constant speed
L=20; %L is the size of the domain on which the particles can move

p=zeros(2,N,T+1); % Position of all individuals at all timesteps
d=zeros(N,T+1); % Desired direction of all individuals at all timesteps
d(:,1) = 2*pi*rand(1,N);
nInformed = round(N*prop);
g=pi/4*ones(nInformed,1); % Preferred direction of informed individuals
g = g + normrnd(0,gam,nInformed,1);
pDiff = zeros(2,N,N);
p(:,:,1) =rand(2,N); % Random initial position

elong = zeros(1,T);
groupdir = zeros(1,T);
% T(i,j) gives the angle with the x axis of the direction of motion of the ith
% particle at time j
%T=zeros(N,J+1);
%T(:,1)=2*pi*rand(N,1); %define initial direction of all particles
%fi = zeros(1,J+1);
%fi(1) = abs(1/N*sqrt(sum(cos(T(:,1)))^2+sum(sin(T(:,1)))^2)); % Global alignment

%For all time steps
for t=1:T
	d(:,t+1) = d(:,t);
	D = squareform(pdist(p(:,:,t)'),'tomatrix');
	for i = 1:N
		for j = 1:N
			pDiff(:,i,j) = p(:,i,t)-p(:,j,t);
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
            repel = -sum(bsxfun(@rdivide,pDiff(:,tooCloseIdx,i),...
                D(i,tooCloseIdx)),2);
			%d(i,t+1) = (d(i,t) + atan2(repel(2),repel(1)))/2;
			d(i,t+1) = atan2(repel(2),repel(1));
		elseif closeEnoughIdx
			attract = sum(bsxfun(@rdivide,...
                pDiff(:,closeEnoughIdx(closeEnoughIdx~=i),i),...
                D(i,closeEnoughIdx(closeEnoughIdx~=i))),2);
			d(i,t+1) = (atan2(attract(2),attract(1)) +...
                mean(d(closeEnoughIdx,t)))/2;
		end
	end
    
    % Turn towards preferred direction
	d(1:nInformed,t+1) = d(1:nInformed,t+1) + w*(g-d(1:nInformed,t+1));
    
    % Add a random (gaussian) angle to direction
    randAngles = normrnd(0,0.01,1,size(d(:,t+1),2));
    d(:,t+1) = d(:,t+1) + randAngles;
    
    % Normalize angles between 0 and 2*pi
    d(:,t+1) = mod(d(:,t+1),2*pi);
    
    % Calculate mean direction of the particles
    groupdir(t) = mean(d(:,t+1));
	
    % Update positions
	p(:,:,t+1) = p(:,:,t) + dt*s*[cos(d(:,t+1)'); sin(d(:,t+1)')];
	
    % Calculate elongation
	[bb, e] = orientedBoundingBox(p(:,:,t), groupdir(t));
	elong(t) = e;
	
	if plotornot
		qscale = 1;
		quiver(p(1,:,t),p(2,:,t),cos(d(:,t)'*qscale),sin(d(:,t)'*qscale),...
			'AutoScaleFactor',0.1,'Marker','.','Markersize',10,'Color','k');
		hold on
		plot(bb(1,[1:end 1]), bb(2,[1:end 1]),'r-')
		hold off
		xlabel('X position')
		ylabel('Y position')
        %xlim([0 40])
        %ylim([0 40])
		drawnow
		pause(0.001)
	end
end