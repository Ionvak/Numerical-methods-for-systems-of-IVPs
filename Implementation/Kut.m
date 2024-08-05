% This function returns the numerical solution to a given system of IVPs
% using the Kutta method.
function z_K = Kut(f,dt,T,Z)

timespan = 0:dt:T;

z_K = zeros(12,size(timespan,2));
z_K(:,1) = Z;

for j = 2:size(timespan,2)
    f1 = f(0,z_K(:,j-1));
    f2 = f(0,z_K(:,j-1)+(dt/2)*f1);
    f3 = f(0,z_K(:,j-1)-dt*f1+2*dt*f2);
    z_K(:,j) = z_K(:,j-1) + (dt/6)*(f1+4*f2+f3);
end