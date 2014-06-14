load results_tot
acc = (0.5-accuracy)*2;
styles = {'rv-','b^-','gs-','ko-'};
figure()
hold on
for i = 1:length(acc(:,1));
    plot(pvals,acc(i,:),styles{i},'LineWidth',2);
end
hold off
xlabel('Proportion of informed individuals')
ylabel('Accuracy')
legend('N = 10','N = 30','N = 50','N = 100')

figure()
hold on
for i = 1:length(elongvals(:,1));
    plot(pvals,elongvals(i,:),styles{i},'LineWidth',2);
end
hold off
xlabel('Proportion of informed individuals')
ylabel('Elongation')
legend('N = 10','N = 30','N = 50','N = 100')