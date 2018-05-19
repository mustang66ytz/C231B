function xNext = discreteModel(x, u)
m = 2237;
Jz = 5112;
a = 1.46;
b = 1.55;

[Ff, Fr] = tireModel(u, x);
xDot = [x(2)*x(3) - 2/m*Ff*sin(u(2)) + u(1)/m;
    -x(1)*x(3) + 2/m*(Ff*cos(u(2)) + Fr);
    2/Jz*(a*Ff*cos(u(2)) - b*Fr);
    x(1)*cos(x(6)) - x(2)*sin(x(6));
    x(1)*sin(x(6)) + x(2)*cos(x(6));
    x(3)];
xNext = x + xDot*0.01;
end