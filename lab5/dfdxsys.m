function result = dfdxsys(zVal, uVal)
syms vx vy wz X Y psi_ Fx delta;
z = [vx; vy; wz; X; Y; psi_];
u = [Fx; delta];

result = jacobian(discreteModel(z, u),z);
result = subs(result,{vx, vy, wz, X, Y, psi_, Fx, delta},...
    {zVal(1), zVal(2), zVal(3), zVal(4), zVal(5), zVal(6), uVal(1), uVal(2)});
result = double(result);
end