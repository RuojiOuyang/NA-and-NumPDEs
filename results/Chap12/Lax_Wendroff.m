clc;clear;
x1 = [-2:0.05:25]';
y1 = exp(-20*(x1-19).^2)+exp(-(x1-22).^2);
plot(x1,y1,':');
axis([15 25 -0.4 1]);
hold on;
y2 = exp(-20*(x1-2).^2)+exp(-(x1-5).^2);
A1 = diag(repmat([1],1,541));
A2 = diag(repmat([-1],1,540),-1)+diag(repmat([1],1,540),1) ... 
    +diag(repmat([1],1,1),-540)+diag((repmat([-1],1,1)),540);
A3 = diag(repmat([1],1,540),-1)+diag(repmat([1],1,540),1) ... 
    +diag(repmat([1],1,1),-540)+diag((repmat([1],1,1)),540) ...
    +diag(repmat([-2],1,541));
A = A1-0.4.*A2+0.8*0.4.*A3;
for i = 0:0.04:17
    y2 = A*y2;
end
plot(x1,y2,'.-');
title('Lax-Wendroff');
hold off;