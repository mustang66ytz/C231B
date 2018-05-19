function dxdt = kinematicsModel(x, delta, vx, L)
dxdt = [vx*cos(x(3)); vx*sin(x(3)); vx*tan(delta)/L];
end