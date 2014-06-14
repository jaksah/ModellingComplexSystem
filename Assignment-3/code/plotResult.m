clear all
close all

load accuracy10
tot_acc(1,:) = accuracy; 
load y10
tot_y(1,:) = y;
load elongvals10
tot_elong(1,:) = elongvals;

load Vals30
tot_y(2,:) = y30;
tot_acc(2,:) = accuracy30;
tot_elong(2,:) = elongvals30;

load y50
tot_y(3,:) = y;
load accuracy50
tot_acc(3,:) = accuracy;
load elongvals50
tot_elong(3,:) = elongvals;

load result100
tot_y(4,:) = y100;
tot_acc(4,:) = accuracy100;
tot_elong(4,:) = elongvals100;

p = linspace(0,1, size(tot_y,2));
cPicker = 'ko-';
cPicker2 = 'k:';
Ns = [10 30 50 100];
w = 0.5;
a = 1;
for i = 1:size(tot_y,1)
	N = Ns(i);
	%y_approx(i,:) = (N./(((sqrt(w*N.*p)-1)*a).*((sqrt(N.*p)-1)*a)));
	%y_approx(i,:) = (((N./sqrt(N*p*w)-1)*a)./(((sqrt(w*N.*p)-1)*a)));
	%y_approx(i,:) = y_approx(i,:)./(1+(y_approx(i,:)).*exp(-p*10000));
	y_approx(i,:) = 1./(p.*w);
end

for i = 1:size(tot_y,1)
	figure;
	plot(p,tot_elong(i,:), cPicker)
	hold on
	plot(p,y_approx(i,:),cPicker2);
	title(['N = ', num2str(Ns(i))]);
	ax1 = ylabel('Elongation');
	ax2 = xlabel('Proportion of informed individuals');
	legend('Simulation', 'Approximation')
	set(gca, 'FontSize', 16, 'YLim', [0 10])
	set(ax1, 'FontSize', 16);
	set(ax2, 'FontSize', 16);
	mystr = sprintf('elongation%d.tikz',Ns(i));
	matlab2tikz_short(mystr);
end
%ph = linspace(0,1,101);


%legend('N=10', 'N=30', 'N=50', 'N=100')





