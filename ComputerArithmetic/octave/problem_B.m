hold on
norm = [0 0.5 -0.5 0.75 -0.75 0.625 -0.625 0.875 -0.875 1 -1 1.5 -1.5 1.25 -1.25 1.75 -1.75 2 -2 3 -3 2.5 -2.5 3.5 -3.5];
subnorm = [0.25 -0.25 0.125 -0.125 0.375 -0.375];
extended = [norm, subnorm];

lnorm = plot(0,1,'ko','MarkerSize',6,'MarkerFaceColor','k');
lsub = plot(0.25,1,'kd','MarkerSize',6,'MarkerFaceColor','k');
norm_axis = plot([-4.1 4.1],[0 0],'-r','LineWidth',4);
extend_axis = plot([-4.1 4.1],[1 1],'-b','LineWidth',4);

plot(norm,0,'ko','MarkerSize',6,'MarkerFaceColor','k');
plot(norm,1,'ko','MarkerSize',6,'MarkerFaceColor','k');
plot(subnorm,1,'kd','MarkerSize',6,'MarkerFaceColor','k');

legend('normalized numbers', 'subnormalized numbers', 'normalized FPN', 'extended FPN')
legend('FontSize',15,'Location','EastOutside')
% plot(1,0,'bo','MarkerSize',5,'MarkerFaceColor','b')

%legend(lnorm, 'normalized numbers', lsub, 'subnormalized numbers');
daspect([1 1 1])

height = -0.5;
fs = 20;
text(-4,height,'-4','Horiz','center','Vert','bottom','Fontsize',fs)
text(-3,height,'-3','Horiz','center','Vert','bottom','Fontsize',fs)
text(-2,height,'-2','Horiz','center','Vert','bottom','Fontsize',fs)
text(-1,height,'-1','Horiz','center','Vert','bottom','Fontsize',fs)
text(0,height,'0','Horiz','center','Vert','bottom','Fontsize',fs)
text(1,height,'1','Horiz','center','Vert','bottom','Fontsize',fs)
text(2,height,'2','Horiz','center','Vert','bottom','Fontsize',fs)
text(3,height,'3','Horiz','center','Vert','bottom','Fontsize',fs)
text(4,height,'4','Horiz','center','Vert','bottom','Fontsize',fs)

axis off

scrsz=get(0,'ScreenSize');
set(gcf,'Position',scrsz);
saveas(gcf, './images/B.png')