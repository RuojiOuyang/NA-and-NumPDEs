%% For n = 1
figure(1)
syms x
y11 = (1 + sign(1-x))./2.*(1-x);
y12 = (1 + sign(2-x))./2.*(2-x); y22 = y12-y11;
y13 = (1 + sign(3-x))./2.*(3-x); y23 = y13-y12; y33 = y23 - y22;
%%
subplot(3,3,1)
x = 0:0.01:1;
plot(x,subs(y11,x));
xlim([0,3])
set(gca,'XTick',[0 1 2 3]);
set(gca,'XTicklabel',{'','t_{i-1}','t_{i}','t_{i+1}'},'FontSize',15)
daspect([1 1 1])
%%
subplot(3,3,4)
x = 1:0.01:2;
plot(x,subs(y12,x));
xlim([0,3])
set(gca,'XTick',[0 1 2 3]);
set(gca,'XTicklabel',{'','t_{i-1}','t_{i}','t_{i+1}'},'FontSize',15)
daspect([1 1 1])
%%
subplot(3,3,7)
x = 2:0.01:3;
plot(x,subs(y13,x));
xlim([0,3])
set(gca,'XTick',[0 1 2 3]);
set(gca,'XTicklabel',{'','t_{i-1}','t_{i}','t_{i+1}'},'FontSize',15)
daspect([1 1 1])
%%
subplot(3,3,5)
x = 0:0.01:3;
plot(x,subs(y22,x));
xlim([0,3])
set(gca,'XTick',[0 1 2 3]);
set(gca,'XTicklabel',{'','t_{i-1}','t_{i}','t_{i+1}'},'FontSize',15)
daspect([1 1 1])
%%
subplot(3,3,8)
x = 0:0.01:3;
plot(x,subs(y23,x));
xlim([0,3])
set(gca,'XTick',[0 1 2 3]);
set(gca,'XTicklabel',{'','t_{i-1}','t_{i}','t_{i+1}'},'FontSize',15)
daspect([1 1 1])
%%
subplot(3,3,9)
x = 0:0.01:3;
plot(x,subs(y33,x));
xlim([0,3])
ylim([0,1])
set(gca,'XTick',[0 1 2 3]);
set(gca,'XTicklabel',{'','t_{i-1}','t_{i}','t_{i+1}'},'FontSize',15)
daspect([1 1 1])
sgtitle('Linear B-spline','FontSize',25)
scrsz=get(0,'ScreenSize');
set(gcf,'Position',scrsz);
saveas(gcf, './images/G1.png')
%% For N = 2
figure(2)
syms x
y11 = (1 + sign(1-x))./2.*(1-x).^2;
y12 = (1 + sign(2-x))./2.*(2-x).^2; y22 = (y12 - y11);
y13 = (1 + sign(3-x))./2.*(3-x).^2; y23 = (y13 - y12); y33 = (y23 - y22)/2;
y14 = (1 + sign(4-x))./2.*(4-x).^2; y24 = (y14 - y13); y34 = (y24 - y23)/2; y44 = (y34 - y33);
%%
subplot(4,4,1)
x = 0:0.01:1;
plot(x,subs(y11,x));
xlim([0,4])
set(gca,'XTick',[0 1 2 3 4]);
%set(gca,'YTick',[]);
set(gca,'XTicklabel',{'','t_{i-1}','t_{i}','t_{i+1}','t_{i+2}'},'FontSize',15)
%%
subplot(4,4,5)
x = 1:0.01:2;
plot(x,subs(y12,x));
xlim([0,4])
set(gca,'XTick',[0 1 2 3 4]);
%set(gca,'YTick',[]);
set(gca,'XTicklabel',{'','t_{i-1}','t_{i}','t_{i+1}','t_{i+2}'},'FontSize',15)
%%
subplot(4,4,9)
x = 2:0.01:3;
plot(x,subs(y13,x));
xlim([0,4])
set(gca,'XTick',[0 1 2 3 4]);
%set(gca,'YTick',[]);
set(gca,'XTicklabel',{'','t_{i-1}','t_{i}','t_{i+1}','t_{i+2}'},'FontSize',15)
%%
subplot(4,4,13)
x = 3:0.01:4;
plot(x,subs(y14,x));
xlim([0,4])
set(gca,'XTick',[0 1 2 3 4]);
%set(gca,'YTick',[]);
set(gca,'XTicklabel',{'','t_{i-1}','t_{i}','t_{i+1}','t_{i+2}'},'FontSize',15)
%%
subplot(4,4,6)
x = 0:0.01:4;
plot(x, subs(y22,x));
xlim([0,4])
ylim([0,7])
set(gca,'XTick',[0 1 2 3 4]);
%set(gca,'YTick',[]);
set(gca,'XTicklabel',{'','t_{i-1}','t_{i}','t_{i+1}','t_{i+2}'},'FontSize',15)
%%
subplot(4,4,10)
x = 0:0.01:4;
plot(x, subs(y23,x));
xlim([0,4])
ylim([0,7])
set(gca,'XTick',[0 1 2 3 4]);
%set(gca,'YTick',[]);
set(gca,'XTicklabel',{'','t_{i-1}','t_{i}','t_{i+1}','t_{i+2}'},'FontSize',15)
%%
subplot(4,4,14)
x = 0:0.01:4;
plot(x, subs(y24,x));
xlim([0,4])
ylim([0,7])
set(gca,'XTick',[0 1 2 3 4]);
%set(gca,'YTick',[]);
set(gca,'XTicklabel',{'','t_{i-1}','t_{i}','t_{i+1}','t_{i+2}'},'FontSize',15)
%%
subplot(4,4,11)
x = 0:0.01:4;
plot(x, subs(y33,x));
xlim([0,4])
set(gca,'XTick',[0 1 2 3 4]);
%set(gca,'YTick',[]);
set(gca,'XTicklabel',{'','t_{i-1}','t_{i}','t_{i+1}','t_{i+2}'},'FontSize',15)
%%
subplot(4,4,15)
x = 0:0.01:4;
plot(x, subs(y34,x));
xlim([0,4])
set(gca,'XTick',[0 1 2 3 4]);
%set(gca,'YTick',[]);
set(gca,'XTicklabel',{'','t_{i-1}','t_{i}','t_{i+1}','t_{i+2}'},'FontSize',15)
%%
subplot(4,4,16)
x = 0:0.01:4;
plot(x, subs(y44,x));
xlim([0,4])
set(gca,'XTick',[0 1 2 3 4]);
%set(gca,'YTick',[]);
set(gca,'XTicklabel',{'','t_{i-1}','t_{i}','t_{i+1}','t_{i+2}'},'FontSize',15)
sgtitle('Quadratic B-spline','FontSize',25)
scrsz=get(0,'ScreenSize');
set(gcf,'Position',scrsz);
saveas(gcf, './images/G2.png')