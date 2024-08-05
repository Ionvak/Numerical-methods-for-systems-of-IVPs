% This function returns the numerical solution to a given system of IVPs
% using the 4th order Runge-Kutta method.
function z_RK = RK4(f,dt,T,Z)

timespan = 0:dt:T;

z_RK = zeros(12,size(timespan,2));
z_RK(:,1) = Z;

for k = 2:size(timespan,2)
    f1 = f(0,z_RK(:,k-1));
    f2 = f(0,z_RK(:,k-1)+(dt/2)*f1);
    f3 = f(0,z_RK(:,k-1)+(dt/2)*f2);
    f4 = f(0,z_RK(:,k-1)+dt*f3);
    z_RK(:,k) = z_RK(:,k-1) + (dt/6)*(f1+2*f2+2*f3+f4);
end