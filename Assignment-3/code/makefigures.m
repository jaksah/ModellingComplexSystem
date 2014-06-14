load results_tot
acc = (0.5-accuracy)*2;
styles = {'rv-','b^-','gs-','ko-'};
figure()
hold on
for i = 1:length(acc(:,1));
    plot(pvals,acc(i,:),styles{i},'LineWidth',2);
end
hold off
ax1 = xlabel('Proportion of informed individuals');
ax2 = ylabel('Accuracy');
legend('N = 10','N = 30','N = 50','N = 100')
set(gca, 'YLim', [0, 1.1], 'FontSize',16);
set(ax1,'FontSize',16);
set(ax2,'FontSize',16);
matlab2tikz_short('accuracy.tikz');


%{
figure()
hold on
for i = 1:length(elongvals(:,1));
    plot(pvals,elongvals(i,:),styles{i},'LineWidth',2);
end
hold off
xlabel('Proportion of informed individuals')
ylabel('Elongation')
legend('N = 10','N = 30','N = 50','N = 100')
%}