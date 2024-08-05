% This function is for use in finding the approximation criterion of the
% numerical solution given the vector of parameters P.
function j = J(P) 

% Importing data from the file "data_13.csv".
data13 = readtable("Documents/MATLAB/enume-2024-mahmoud-elshekh-ali-323930/Assignment C/data_13.csv");
data13 = data13{:,:};
timespan = 0:0.02:5;

% Defining the parameters used in the planar three-body problem.
G = 1;
m1 = P(1);
m2 = P(2);
m3 = P(3);

% Defining the initial conditions for the system of IVPs used in the planar
% three-body problem.
z = [
    data13(1,2);
    data13(1,3);
    (data13(2,2) - data13(1,2))/0.02;
    (data13(2,3) - data13(1,3))/0.02;
    data13(1,4);
    data13(1,5);
    (data13(2,4) - data13(1,4))/0.02;
    (data13(2,5) - data13(1,5))/0.02;
    data13(1,6);
    data13(1,7);
    (data13(2,6) - data13(1,6))/0.02;
    (data13(2,7) - data13(1,7))/0.02;
    ];

% Defining the functions describing the distances between the bodies
% for use in the system of IVPs used in the planar three-body problem.
r12 = @(Z) sqrt((Z(1) - Z(5))^2 + (Z(2) - Z(6))^2);
r23 = @(Z) sqrt((Z(5) - Z(9))^2 + (Z(6) - Z(10))^2);
r13 = @(Z) sqrt((Z(1) - Z(9))^2 + (Z(2) - Z(10))^2);

% Defining the set of IVPs which make up the planar three-body
% problem.
f = @(t,Z) [
    Z(3);
    Z(4);
    -G*m2*(Z(1) - Z(5))/r12(Z)^3 - G*m3*(Z(1) - Z(9))/r13(Z)^3;
    -G*m2*(Z(2) - Z(6))/r12(Z)^3 - G*m3*(Z(2) - Z(10))/r13(Z)^3;
    Z(7);
    Z(8);
    -G*m3*(Z(5) - Z(9))/r23(Z)^3 - G*m1*(Z(5) - Z(1))/r12(Z)^3;
    -G*m3*(Z(6) - Z(10))/r23(Z)^3 - G*m1*(Z(6) - Z(2))/r12(Z)^3;
    Z(11);
    Z(12);
    -G*m1*(Z(9) - Z(1))/r13(Z)^3 - G*m2*(Z(9) - Z(5))/r23(Z)^3;
    -G*m1*(Z(10) - Z(2))/r13(Z)^3 - G*m2*(Z(10) - Z(6))/r23(Z)^3;
    ];

% Using ode45 to obtain a numerical solution using the parameters P
% to compare with the data imported from the file "data_13.csv".
[ts,R] = ode45(f,timespan,z);

% Finding the sought after approximation criterion.
j = sum( sqrt( (R(:,1) - data13(:,2)).^2 ) + sqrt( (R(:,2) - data13(:,3)).^2 ));
j = j + sum( sqrt( (R(:,5) - data13(:,4)).^2 ) + sqrt( (R(:,6) - data13(:,5)).^2 ) );
j = j + sum(sqrt( (R(:,9) - data13(:,6)).^2) + sqrt( (R(:,10) - data13(:,7)).^2 ) );

end