function y = nonlinearVehicleModelMeasure(z, u)

m = 2237;

[Ff, Fr] = tireModel(u, z);
y = [z(1);
     2 / m * (Ff * cos(u(2)) + Fr);
     z(3)];
 
end