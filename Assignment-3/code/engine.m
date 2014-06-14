clear all
%N = 30; % Number of individuals
T = 1000; % Number of timesteps
w = 0.5; % Weight of influence of preferred direction
%p = 0.1; % Proportion of indiviuals with knowledge of preferred direction
a = 1; % Minimum distance between particles
rho = 6; % Radius of particles to interact with
gam = 0; % Uncertainty of preferred direction
plotornot = 0; 
reps = 5;
g = [1;0];
gAngle = atan2(g(2),g(1));

%Nvals = [10 30 50 100 200];
Nvals = [30];
pvals = 0:0.1:1;
elongvals = zeros(length(pvals),length(Nvals));
accuracy = zeros(length(pvals),length(Nvals));
groupdirvals = zeros(length(pvals),length(Nvals));

tic;
for Nidx=1:length(Nvals)
	N = Nvals(Nidx)
	for pidx=1:length(pvals)
		p = pvals(pidx)
		angles = zeros(1,reps);
		for rep=1:reps
			rep
			toc
			[e, cent] = couzin(N, T, w, p, a, rho, gam, g, plotornot);
			vec = cent(:,end) - cent(:,end-50);
			angles(rep) = abs(atan2(vec(2),vec(1)))/pi;
			elongvals(pidx,Nidx) = elongvals(pidx,Nidx) + mean(e(end-10:end));
		end
		accuracy(pidx, Nidx) = mean((angles-gAngle).^2);
		elongvals(pidx,Nidx) = elongvals(pidx,Nidx)/reps;
        groupdirvals(pidx,Nidx) = groupdirvals(pidx,Nidx)/reps;
	end
end

%% 
%pvals = 0:0.1:1;
plot(pvals,accuracy(:,1),'o-')
xlabel('p')
ylabel('Accuracy')
hold on
%%

plot(pvals,elongvals(:,1),'rx-')
xlabel('p')
ylabel('Elongation')
hold on
%pvals = 0:0.01:1;

y = (N./(((sqrt(N.*pvals)-ones(size(pvals)))*a).^2));
%y = (2./pvals);
plot(pvals, y,'k.:')
set(gca, 'YLim', [0 5])
legend('Simulation', 'Approximation','Master Equation')





