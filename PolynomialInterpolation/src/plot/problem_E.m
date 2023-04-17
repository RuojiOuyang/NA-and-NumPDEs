x =0 : 1 : 28;
p1 = 4.1477e-05* x.^6 + -0.00371557 * x.^5 + 0.128281 * x.^4 + -2.11512 * x.^3 + 16.2855 * x.^2 + -43.0127 * x.^1 + 6.67;
p2 = 8.6768e-06* x.^6 + -0.000777473 * x.^5 + 0.0265858 * x.^4 + -0.424283 * x.^3 + 2.98227 * x.^2 + -5.85018 * x.^1 + 6.67;
fig1 = figure;
plot(x,p1,'Linewidth',2)
hold on
plot(x,p2,'Linewidth',2)
set(gca,'FontSize',20)
legend('{\it sp1}','{\it sp2}')
set(legend,'FontSize',20,'FontName', 'Times New Roman', 'Location','NorthEastOutside')
saveas(fig1, '../../output/images/E1.png')

fig2 = figure;
x =0 : 1 : 43;
p1 = 4.1477e-05* x.^6 + -0.00371557 * x.^5 + 0.128281 * x.^4 + -2.11512 * x.^3 + 16.2855 * x.^2 + -43.0127 * x.^1 + 6.67;
p2 = 8.6768e-06* x.^6 + -0.000777473 * x.^5 + 0.0265858 * x.^4 + -0.424283 * x.^3 + 2.98227 * x.^2 + -5.85018 * x.^1 + 6.67;
plot(x,p1,'Linewidth',2)
hold on
plot(x,p2,'Linewidth',2)
set(gca,'FontSize',20)
legend('{\it sp1}','{\it sp2}')
set(legend,'FontSize',20,'FontName', 'Times New Roman', 'Location','NorthEastOutside')
saveas(fig2, '../../output/images/E2.png')
