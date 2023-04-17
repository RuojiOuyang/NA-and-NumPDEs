theta = linspace(0,2*pi,1000);
subplot(2,2,1);
plot(exp(1j.*theta)-1);
axis equal;
set(gca,'Xlim',[-2.5 0.5],'XTick',[-2.5:0.5:0.5]);
set(gca,'Ylim',[-1.5 1.5],'YTick',[-1.5:0.5:1.5]);
hold on;
h = 0.02;
c = [0.8,1.6,2,2.4];
m = 50;
A1 = diag(repmat([1],1,m))-c(1)/2.*(diag(repmat([3],1,m))+diag(repmat([-4],1,m-1),-1) ...
    +diag(repmat([1],1,m-2),-2)+diag(repmat([1],1,2),m-2)+diag(repmat([-4],1,1),m-1)) ...
    +c(1)^2/2.*(diag(repmat([1],1,m))+diag(repmat([-2],1,m-1),-1)+diag(repmat([1],1,m-2),-2) ...
    +diag(repmat([1],1,2),m-2)+diag(repmat([-2],1,1),m-1));
[lam1] = eig(A1);
plot(lam1-1,'.');
hold off;

subplot(2,2,2);
plot(exp(1j.*theta)-1);
axis equal;
set(gca,'Xlim',[-2.5 0.5],'XTick',[-2.5:0.5:0.5]);
set(gca,'Ylim',[-1.5 1.5],'YTick',[-1.5:0.5:1.5]);
hold on;
h = 0.02;
c = [0.8,1.6,2,2.4];
m = 50;
A1 = diag(repmat([1],1,m))-c(2)/2.*(diag(repmat([3],1,m))+diag(repmat([-4],1,m-1),-1) ...
    +diag(repmat([1],1,m-2),-2)+diag(repmat([1],1,2),m-2)+diag(repmat([-4],1,1),m-1)) ...
    +c(2)^2/2.*(diag(repmat([1],1,m))+diag(repmat([-2],1,m-1),-1)+diag(repmat([1],1,m-2),-2) ...
    +diag(repmat([1],1,2),m-2)+diag(repmat([-2],1,1),m-1));
[lam1] = eig(A1);
plot(lam1-1,'.');
hold off;

subplot(2,2,3);
plot(exp(1j.*theta)-1);
axis equal;
set(gca,'Xlim',[-2.5 0.5],'XTick',[-2.5:0.5:0.5]);
set(gca,'Ylim',[-1.5 1.5],'YTick',[-1.5:0.5:1.5]);
hold on;
h = 0.02;
c = [0.8,1.6,2,2.4];
m = 50;
A1 = diag(repmat([1],1,m))-c(3)/2.*(diag(repmat([3],1,m))+diag(repmat([-4],1,m-1),-1) ...
    +diag(repmat([1],1,m-2),-2)+diag(repmat([1],1,2),m-2)+diag(repmat([-4],1,1),m-1)) ...
    +c(3)^2/2.*(diag(repmat([1],1,m))+diag(repmat([-2],1,m-1),-1)+diag(repmat([1],1,m-2),-2) ...
    +diag(repmat([1],1,2),m-2)+diag(repmat([-2],1,1),m-1));
[lam1] = eig(A1);
plot(lam1-1,'.');
hold off;

subplot(2,2,4);
plot(exp(1j.*theta)-1);
axis equal;
set(gca,'Xlim',[-2.5 0.5],'XTick',[-2.5:0.5:0.5]);
set(gca,'Ylim',[-1.5 1.5],'YTick',[-1.5:0.5:1.5]);
hold on;
h = 0.02;
c = [0.8,1.6,2,2.4];
m = 50;
A1 = diag(repmat([1],1,m))-c(4)/2.*(diag(repmat([3],1,m))+diag(repmat([-4],1,m-1),-1) ...
    +diag(repmat([1],1,m-2),-2)+diag(repmat([1],1,2),m-2)+diag(repmat([-4],1,1),m-1)) ...
    +c(4)^2/2.*(diag(repmat([1],1,m))+diag(repmat([-2],1,m-1),-1)+diag(repmat([1],1,m-2),-2) ...
    +diag(repmat([1],1,2),m-2)+diag(repmat([-2],1,1),m-1));
[lam1] = eig(A1);
plot(lam1-1,'.');
hold off;