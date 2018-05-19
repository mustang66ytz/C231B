function dzdt = nonlinearVehicleModel(z, u)

m = 2237;
Jz = 5112;
a = 1.46;
b = 1.55;

[Ff, Fr] = tireModel(u, z);
dzdt = [z(2)*z(3)-(2/m)*Ff*sin(u(2))+u(1)/m;
    -z(1)*z(3)+(2/m)*(Ff*cos(u(2)) + Fr);
    (2/Jz)*(a*Ff*cos(u(2)) - b*Fr);
    z(1)*cos(z(6)) - z(2)*sin(z(6));
    z(1)*sin(z(6)) + z(2)*cos(z(6));
    z(3)];


end