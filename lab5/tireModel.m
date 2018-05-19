function [Ff, Fr] = tireModel(u, z)
m = 2237;
a = 1.46;
b = 1.55;
A = -6.8357;
B = 0.0325;
C = 238.9874;

alphaf = atan((z(2)+a*z(3))/z(1) - u(2));
alphar = atan((z(2)-b*z(3))/z(1));
Ff = m*a*A*sin(C*atan(B*alphaf))/(a+b);
Fr = m*a*A*sin(C*atan(B*alphar))/(a+b);

end