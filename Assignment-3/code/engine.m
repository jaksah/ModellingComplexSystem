
%N = 30; % Number of individuals
T = 1500; % Number of timesteps
w = 0.5; % Weight of influence of preferred direction
%p = 0.1; % Proportion of indiviuals with knowledge of preferred direction
a = 1; % Minimum distance between particles
rho = 6; % Radius of particles to interact with
gam = 0; % Uncertainty of preferred direction
plotornot = 0; 
reps = 1;

%Nvals = [10 30 50 100 200];
Nvals = [30];
pvals = 0:0.05:1;
elongvals = zeros(length(pvals),length(Nvals));
groupdirvals = zeros(length(pvals),length(Nvals));
tic;
for Nidx=1:length(Nvals)
	N = Nvals(Nidx)
	for pidx=1:length(pvals)
		p = pvals(pidx)
		for rep=1:reps
			rep
			toc
			[e, groupdir] = couzin(N, T, w, p, a, rho, gam, plotornot);
			elongvals(pidx,Nidx) = elongvals(pidx,Nidx) + mean(e(end-10:end));
		end
		elongvals(pidx,Nidx) = elongvals(pidx,Nidx)/reps;
        groupdirvals(pidx,Nidx) = groupdirvals(pidx,Nidx)/reps;
	end
end
%%

%%
figure()
plot(pvals,elongvals(:,1),'o-')
xlabel('p')
ylabel('Group elongation')

figure()
plot(pvals,accuracy(:,1),'o-')
xlabel('p')
ylabel('Group accuracy')