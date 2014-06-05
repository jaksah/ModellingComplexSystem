%% Flocks with nearest neighbours
%function fi = couzin(N, T, w, p, a, rho, plotornot)
N = 100;
T = 1000;
w = 0.5;
p = 0.5;
a = 0.05;
rho = 0.5;
plotornot = 1;

%{ 
	N - Number of individuals
	T - Number of timesteps
	w - Weight of influence of preferred direction
	p - Proportion of indiviuals with knowledge of preferred direction
	a - Minimum distance between particles
%}
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
    
	v(:,:,t+1) = s*d(:,:,t);
	c(:,:,t+1) = c(:,:,t) + dt*v(:,:,t+1);
	
	if plotornot
		qscale = 1;
		quiver(c(1,:,t),c(2,:,t),d(1,:,t)*qscale,d(2,:,t)*qscale,...
			'AutoScaleFactor',0.1,'Marker','.','Markersize',10,'Color','k');
% 		plot([x(i,j), x(i,j+1)] ,[y(i,j),y(i,j+1)],'k-','markersize',4) %plots the first half of the particles in black
% 		axis([0 L 0 L]);
% 		hold on
% 		plot(x(i,j+1) ,y(i,j+1),'k.','markersize',10)
		xlabel('X position')
		ylabel('Y position')
        xlim([0 10])
        ylim([0 10])
		drawnow
		pause(0.1)
	end
% 	
% 	%For each particle
% 	for i=1:N
% 		%finds how many particles are in the interaction radius of each
% 		%particle
% % 		xDiff(:,1) = x(i,j)-x(:,j);
% % 		xDiff(:,2) = x(i,j)-x(:,j)-L;
% % 		xDiff(:,3) = x(i,j)-x(:,j)+L;
% % 		
% % 		yDiff(:,1) = y(i,j)-y(:,j);
% % 		yDiff(:,2) = y(i,j)-y(:,j)-L;
% % 		yDiff(:,3) = y(i,j)-y(:,j)+L;
% % 		
% % 		for xIdx = 1:size(xDiff,2)
% % 			for yIdx = 1:size(yDiff,2)
% % 				A(:,(yIdx-1)*3+xIdx) = ((xDiff(:,xIdx)).^2+(yDiff(:,yIdx)).^2).^0.5;
% % 			end
% % 		end
% 
% 		D = min(A,[],2); % Calculate smallest distance to particle i
% 		[sortedD, sortedIdxD] = sort(D);
% 		nNearestIdx = sortedIdxD(1:n+1); % 1 is the particle itself.
% 		
% 		% Get those particles that are too near particle i
% 		tooNear = ((sortedD(2:end) < sep)).*sortedIdxD(2:end);
% 		tooNear(tooNear==0)=[];
% 		
% 		repDirX = 0;
% 		repDirY = 0;
% 		% Calculate the repulsion term and weight it with 1/d^2
% 		if(useRepulsion && ~isempty(tooNear))
% 			diffX = min(xDiff(tooNear),[],2);
% 			diffY = min(yDiff(tooNear),[],2);
% 			vLen = (diffX.^2+diffY.^2).^0.5;
% 			diffX = diffX./vLen./D(tooNear).^2;
% 			diffY = diffY./vLen./D(tooNear).^2;
% 			repDirX = sum(diffX);
% 			repDirY = sum(diffY);
% 			repDirX = repDirX./length(tooNear);
% 			repDirY = repDirY./length(tooNear);
% 		end
% 		
% 		sc=(alpha)*sum(cos(T(nNearestIdx,j)))/(n+1) + (1-alpha)*repDirX;
% 		ss=(alpha)*sum(sin(T(nNearestIdx,j)))/(n+1) + (1-alpha)*repDirY;
% 		S=atan2(ss,sc);
% 		
% 		
% 		T(i,j+1)=S+e*(rand-0.5); %adds noise to the measured angle
% 		x(i,j+1)=x(i,j)+v*cos(T(i,j+1)); %updates the particles' x-coordinates
% 		y(i,j+1)=y(i,j)+v*sin(T(i,j+1)); %updates the particles' y coordinates
% 		% Jumps from the right of the box to the left or vice versa
% 		x(i,j+1)=mod(x(i,j+1),L);
% 		%Jumps from the top of the box to the bottom or vice versa
% 		y(i,j+1)=mod(y(i,j+1),L);
% 		%Plot particles
% 		if makemovie
% 			if abs(x(i,j)-x(i,j+1))<v & abs(y(i,j)-y(i,j+1))<v
% 				plot([x(i,j), x(i,j+1)] ,[y(i,j),y(i,j+1)],'k-','markersize',4) %plots the first half of the particles in black
% 				axis([0 L 0 L]);
% 				hold on
% 				plot(x(i,j+1) ,y(i,j+1),'k.','markersize',10)
% 				xlabel('X position')
% 				ylabel('Y position')
% 			end
% 		end
% 	end
% 	
% 	fi(j+1) = abs(1/N*sqrt(sum(cos(T(:,j+1)))^2+sum(sin(T(:,j+1)))^2)); % Global alignment
% 	if makemovie
% 		hold off
% 		M(j)=getframe; %makes a movie fram from the plot
% 		%movien = addframe(movien,M(j)); %adds this movie fram to the movie
% 	end
	
end


%movien = close(movien); %finishes the movie



