clc;clear;
x = [-2:0.05:25]';
y1 = exp(-20*(x-19).^2)+exp(-(x-22).^2);
y2 = exp(-20*(x-2).^2)+exp(-(x-5).^2);
plot(x,y1,':');
axis([15 25 -0.4 1]);
hold on;
A1 = diag(repmat([1],1,540),-1)+diag(repmat([1],1,540),1) ... 
    +diag(repmat([1],1,1),-540)+diag((repmat([1],1,1)),540);
A2 = diag(repmat([-1],1,540),-1)+diag(repmat([1],1,540),1) ... 
    +diag(repmat([1],1,1),-540)+diag((repmat([-1],1,1)),540);
A = 0.5.*A1-0.4.*A2;
for i = 0:0.04:17
    y2 = A*y2;
end
plot(x,y2,'.-');
title('Lax-Friedrichs')
hold off;