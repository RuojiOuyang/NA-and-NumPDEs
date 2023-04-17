clc;clear;
x1 = [-2:0.05:25]';
k = 0.05;
y1 = exp(-20*(x1-19).^2)+exp(-(x1-22).^2);
y2 = exp(-20*(x1-2).^2)+exp(-(x1-5).^2);
y3 = exp(-20*(x1-2-k).^2)+exp(-(x1-5-k).^2);
plot(x1,y1,':');
axis([15 25 -0.4 1]);
hold on;
A = diag(repmat([-1],1,540),-1)+diag(repmat([1],1,540),1) ... 
    +diag(repmat([1],1,1),-540)+diag((repmat([-1],1,1)),540);
for i = 0.1:0.05:17
    y = y2-A*y3;
    y2 = y3;
    y3 = y;
end
plot(x1,y3,'.-');
title('leapfrog')
hold off;