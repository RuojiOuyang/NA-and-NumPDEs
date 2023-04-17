eta1 = 1;
eta2 = 1.5;
K = [0.2, 0.1, 0.05];
E = zeros(6,2);
lam = -1e+6;
y1 = zeros(1,30);
y2 = zeros(1,30);
for j = 1:3
    u0 = eta1;
    k = K(j);
    t = 0;
    u1 = u0;
    u2 = u0;
    s = 1;
    while (t<3)
        u1 = (u1-k*lam*cos(t+k)-k*sin(t+k))/(1-k*lam);
        u2 = (u2+k/2*lam*u2-k/2*(lam*(cos(t)+cos(t+k))+sin(t)+sin(t+k)))/(1-k*lam/2);
        if (j == 2)
            y1(s) = u2;
            s = s+1;
        end
        t = t+k;
    end
    E(j,1) = abs(u1-(u0-1)*exp(3*lam)-cos(3));
    E(j,2) = abs(u2-(u0-1)*exp(3*lam)-cos(3));
end

for j = 1:3
    u0 = eta2;
    k = K(j);
    t = 0;
    u1 = u0;
    u2 = u0;
    s = 1;
    while (t<3)
        u1 = (u1-k*lam*cos(t+k)-k*sin(t+k))/(1-k*lam);
        u2 = (u2+k/2*lam*u2-k/2*(lam*(cos(t)+cos(t+k))+sin(t)+sin(t+k)))/(1-k*lam/2);
        if (j == 2)
            y2(s) = u2;
            s = s+1;
        end
        t = t+k;
    end
    E(3+j,1) = abs(u1-(u0-1)*exp(3*lam)-cos(3));
    E(3+j,2) = abs(u2-(u0-1)*exp(3*lam)-cos(3));
end
disp(E);
x = linspace(0,3,31);
y1 = [1,y1];
y2 = [1.5,y2];
plot(x,y1,'k--o',x,cos(x),'k-');
title('eta = 1');
plot(x,y2,'k--o',x,0.5*exp(lam*x)+cos(x),'k-');
title('eta = 1.5');