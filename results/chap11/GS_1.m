x = linspace(0,1,21);
y0 = zeros(19,1);
y1 = zeros(19,1);
y2 = zeros(19,1);
y3 = zeros(19,1);
y0(10) = 1;
A = -2*diag(diag(ones(19,19)))+diag(diag(ones(18,18)),1)+diag(diag(ones(18,18)),-1);
I = diag(diag(ones(19,19)));
for i = 1:10
    y0 = (2.*I+A)/(2.*I-A)*y0;
    if i == 1
        y1 = y0;
    end
    if i == 2
        y2 = y0;
    end
    if i == 10
        y3 = y0;
    end
end
y1 = [0;y1;0];
y2 = [0;y2;0];
y3 = [0;y3;0];

subplot(1,3,1);
plot(x,y1,'.-');
axis([0 1 0 0.5]);
subplot(1,3,2);
plot(x,y2,'.-');
axis([0 1 0 0.5]);
subplot(1,3,3);
plot(x,y3,'.-');
axis([0 1 0 0.5]);