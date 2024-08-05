% This function returns the numerical solution to a given system of IVPs
% using the 2nd order Adam-bashford method.
function z_AB = AB2(f,dt,T,Z)

timespan = 0:dt:T;

z_AB = zeros(12,size(timespan,2));
z_AB(:,1) = Z;
z_AB(:,2) = Z + dt*f(0,Z);

for i = 3:size(timespan,2)
    z_AB(:,i) = z_AB(:,i-1) + dt*((3/2)*f(0,z_AB(:,i-1)) - (1/2)*f(0,z_AB(:,i-2)));
end