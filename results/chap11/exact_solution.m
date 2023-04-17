x = linspace(0,1,101);
y1 = zeros(1,101);
y2 = zeros(1,101);
y3 = zeros(1,101);
t_0 = 1/400;
for k = 1:1000
    t = t_0;
    A = 40/(k*k*pi*pi)*(-sin(9*k*pi/20)+2*sin(k*pi/2)-sin(11*k*pi/20));
    y1 = y1+A*exp(-k*k*pi*pi*t).*sin(k*pi*x);
    y2 = y2+A*exp(-k*k*pi*pi*2*t).*sin(k*pi*x);
    y3 = y3+A*exp(-k*k*pi*pi*10*t).*sin(k*pi*x);
end

subplot(1,3,1);
plot(x,y1,'-');
axis([0 1 0 0.5]);
subplot(1,3,2);
plot(x,y2,'-');
axis([0 1 0 0.5]);
subplot(1,3,3);
plot(x,y3,'-');
axis([0 1 0 0.5]);
