
function dxdt = ode_orbitdynamics_nonlinear(t,x)
  
%Normalized variables
%mu = M_earth/M;
mu = 3.040357143 * 10^-6;
omega = 1;

%Centered at barycentre
r = 1;
x_1 = -mu;
x_2 = 1- mu;

pos_sun = [x_1, 0, 0];
pos_earth = [x_2, 0, 0];

v_r = [x(1), x(3), x(5)];

r1 = norm(v_r - pos_sun);
r2 = norm(v_r - pos_earth);

%2. Define new variables: y=x(1), ydot=x(2), etc
x1 = x(1); %x
x2 = x(2); %xdot
y1 = x(3); %y
y2 = x(4); %ydot
z1 = x(5); %z
z2 = x(6); %zdot

%3. Write out your first order system of ODE’s:
dxdt = zeros(size(x));
dxdt(1) = x2;
dxdt(2) = 2*y2 + (1 - (1-mu)/r1^3 - mu/r2^3)* x1 + mu*(1-mu)/r2^3 - mu*(1-mu)/r1^3;
dxdt(3) = y2;
dxdt(4) = -2*x2 + (1 - (1-mu)/r1^3 - mu/r2^3)* y1;
dxdt(5) = z2;
dxdt(6) = (- (1-mu)/r1^3 - mu/r2^3)* z1;

end


