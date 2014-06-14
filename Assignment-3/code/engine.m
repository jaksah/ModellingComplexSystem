
%N = 30; % Number of individuals
T = 1500; % Number of timesteps
w = 0.5; % Weight of influence of preferred direction
%p = 0.1; % Proportion of indiviuals with knowledge of preferred direction
a = 1; % Minimum distance between particles
rho = 6; % Radius of particles to interact with
gam = 0; % Uncertainty of preferred direction
plotornot = 0; 
reps = 1;
g = [1;1];
gAngle = atan2(g(2),g(1));

%Nvals = [10 30 50 100 200];
Nvals = [30];
pvals = 0:0.05:1;
elongvals = zeros(length(pvals),length(Nvals));
accuracy = zeros(length(pvals),length(Nvals));
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
			angles(rep) = atan2(vec(2),vec(1));
			elongvals(pidx,Nidx) = elongvals(pidx,Nidx) + mean(e(end-100:end));
		end
		accuracy(pidx, Nidx) = mean(sum((angles-gangle).^2));
		elongvals(pidx,Nidx) = elongvals(pidx,Nidx)/reps;
	end
end
%%
load dataN50
plot(pvals,elongvals(:,1),'o-')
xlabel('p')
ylabel('Elongitude')
hold on
pvals = 0:0.01:1;

y = (N./(((sqrt(N.*pvals)-ones(size(pvals)))*a).^2));
%y = (2./pvals);
plot(pvals, y,'.')
set(gca, 'YLim', [0 10])
legend('Simulation', 'Approximation')

%% 
plot(pvals,accuracy,'o-')
xlabel('p')
ylabel('Accuracy')
hold on
