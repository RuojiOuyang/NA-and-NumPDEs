x = -5 : 0.01 : 5;
p2 = -0.0384615* x.^2 + 1;
p4 = 0.00530504* x.^4 + -0.171088 * x.^2 + 1;
p6 = -0.000840633* x.^6 + 8.67362e-19 * x.^5 + 0.0335319 * x.^4 + -0.351364 * x.^2 + 2.77556e-17 * x.^1 + 1;
p8 = 0.000137445* x.^8 + -0.00658016 * x.^6 + 0.0981875 * x.^4 + -6.93889e-17 * x.^3 + -0.528121 * x.^2 + 2.77556e-17 * x.^1 + 1;
y = 1./(1+x.^2);
fig = figure;
plot(x,p2,'Linewidth',2)
hold on
plot(x,p4,'Linewidth',2)
plot(x,p6,'Linewidth',2)
plot(x,p8,'Linewidth',2)
plot(x,y,'Linewidth',4)
set(gca,'FontSize',20)
legend = legend('{\it n=2}','{\it n=4}','{\it n=6}','{\it n=8}','{\it f(x)}');
set(legend,'FontSize',20,'FontName', 'Times New Roman', 'Location','NorthEastOutside')
saveas(gcf, '../../output/images/B.png')
