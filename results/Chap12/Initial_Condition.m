x = -2:0.05:8;
y = exp(-20*(x-2).^2)+exp(-(x-5).^2);
plot(x,y,'.-');
axis([-2 8 -0.4 1]);
title('initial condition')