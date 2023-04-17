x =0 : 0.01 : 13;
tick = 0:1:13;
p0 = -0.000182013* x.^8 + 0.00832472 * x.^7 + -0.15313 * x.^6 + 1.45825 * x.^5 + -7.69148 * x.^4 + 22.0325 * x.^3 + -30.2859 * x.^2 + 14.3238 * x.^1 + 75;
y = 81 + 0 .* x;
plot(x,p0,'Linewidth',2)
hold on
plot(x,y,'Linewidth',2)
set(gca,'FontSize',24)
set(gca,'xtick',tick)
legend('{\it y=p''(x)}','{\it y=81}')
set(legend,'FontSize',20,'FontName', 'Times New Roman', 'Location','NorthEastOutside')
saveas(gcf, '../../output/images/D.png')
