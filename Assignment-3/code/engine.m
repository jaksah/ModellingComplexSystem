
N = 10; % Number of individuals
T = 100; % Number of timesteps
w = 1; % Weight of influence of preferred direction
p = 0.1; % Proportion of indiviuals with knowledge of preferred direction
a = 1; % Minimum distance between particles
rho = 6; % Radius of particles to interact with
gamma = 0; % Uncertainty of preferred direction
plotornot = 0; 
reps = 10;

Nvals = [10 30 50 100 200];
pvals = 0:0.1:1;

for N=Nvals
	for p=pvals
		for rep=1:reps
			[elong, groupdir] = couzin(N, T, w, p, a, rho, plotornot);
		end
	end
end