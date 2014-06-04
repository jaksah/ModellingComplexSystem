%% Flocks with nearest neighbours
function fi = couzin(nParticles,tTime,eta,nNeighbours,sep,useRepulsion, plotornot)
%Set up movie
makemovie=plotornot;
if makemovie
	fig=figure;
end
%movien = avifile('Vicsekmovie','FPS',3,'compression','none')
J=tTime; %Number of timestep t0 be used
UJ=1;   %Rate at which film is updated
t=1/J; %Size of one time step
N=nParticles; %Number of particles
e=eta; %e is eta the noise parameter, whose maximum value is 2*pi
%r=1 %The radius of influence of a particle
n=nNeighbours; %Number of nearest neighbours, exclusive itself
alpha = 0.5; % Parameter to weight the effect from neighbours vs repulsion [0 1]
L=20; %L is the size of the domain on which the particles can move
v=0.5; %v is the speed at which the particles move
% x(i,j) gives the x coordinate of the ith particle at time j
x=zeros(N,J+1);
x(:,1)=L*rand(N,1); %define initial x coordiantes of all particles
% y(i,j) gives the y coordinate of the ith particle at time j
y=zeros(N,J+1);
y(:,1)=L*rand(N,1); %define initial y coordiantes of all particles
% T(i,j) gives the angle with the x axis of the direction of motion of the ith
% particle at time j
T=zeros(N,J+1);
T(:,1)=2*pi*rand(N,1); %define initial direction of all particles
fi = zeros(1,J+1);
fi(1) = abs(1/N*sqrt(sum(cos(T(:,1)))^2+sum(sin(T(:,1)))^2)); % Global alignment
%For all time steps
for j=1:J
	%For each particle
	for i=1:N
		%finds how many particles are in the interaction radius of each
		%particle
		xDiff(:,1) = x(i,j)-x(:,j);
		xDiff(:,2) = x(i,j)-x(:,j)-L;
		xDiff(:,3) = x(i,j)-x(:,j)+L;
		
		yDiff(:,1) = y(i,j)-y(:,j);
		yDiff(:,2) = y(i,j)-y(:,j)-L;
		yDiff(:,3) = y(i,j)-y(:,j)+L;
		
		for xIdx = 1:size(xDiff,2)
			for yIdx = 1:size(yDiff,2)
				A(:,(yIdx-1)*3+xIdx) = ((xDiff(:,xIdx)).^2+(yDiff(:,yIdx)).^2).^0.5;
			end
		end
		
		D = min(A,[],2); % Calculate smallest distance to particle i
		[sortedD, sortedIdxD] = sort(D);
		nNearestIdx = sortedIdxD(1:n+1); % 1 is the particle itself.
		
		% Get those particles that are too near particle i
		tooNear = ((sortedD(2:end) < sep)).*sortedIdxD(2:end);
		tooNear(tooNear==0)=[];
		
		repDirX = 0;
		repDirY = 0;
		% Calculate the repulsion term and weight it with 1/d^2
		if(useRepulsion && ~isempty(tooNear))
			diffX = min(xDiff(tooNear),[],2);
			diffY = min(yDiff(tooNear),[],2);
			vLen = (diffX.^2+diffY.^2).^0.5;
			diffX = diffX./vLen./D(tooNear).^2;
			diffY = diffY./vLen./D(tooNear).^2;
			repDirX = sum(diffX);
			repDirY = sum(diffY);
			repDirX = repDirX./length(tooNear);
			repDirY = repDirY./length(tooNear);
		end
		
		sc=(alpha)*sum(cos(T(nNearestIdx,j)))/(n+1) + (1-alpha)*repDirX;
		ss=(alpha)*sum(sin(T(nNearestIdx,j)))/(n+1) + (1-alpha)*repDirY;
		S=atan2(ss,sc);
		
		
		T(i,j+1)=S+e*(rand-0.5); %adds noise to the measured angle
		x(i,j+1)=x(i,j)+v*cos(T(i,j+1)); %updates the particles' x-coordinates
		y(i,j+1)=y(i,j)+v*sin(T(i,j+1)); %updates the particles' y coordinates
		% Jumps from the right of the box to the left or vice versa
		x(i,j+1)=mod(x(i,j+1),L);
		%Jumps from the top of the box to the bottom or vice versa
		y(i,j+1)=mod(y(i,j+1),L);
		%Plot particles
		if makemovie
			if abs(x(i,j)-x(i,j+1))<v & abs(y(i,j)-y(i,j+1))<v
				plot([x(i,j), x(i,j+1)] ,[y(i,j),y(i,j+1)],'k-','markersize',4) %plots the first half of the particles in black
				axis([0 L 0 L]);
				hold on
				plot(x(i,j+1) ,y(i,j+1),'k.','markersize',10)
				xlabel('X position')
				ylabel('Y position')
			end
		end
	end
	
	fi(j+1) = abs(1/N*sqrt(sum(cos(T(:,j+1)))^2+sum(sin(T(:,j+1)))^2)); % Global alignment
	if makemovie
		hold off
		M(j)=getframe; %makes a movie fram from the plot
		%movien = addframe(movien,M(j)); %adds this movie fram to the movie
	end
	
end


%movien = close(movien); %finishes the movie



