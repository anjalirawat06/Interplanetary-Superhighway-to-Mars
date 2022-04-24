%Normalized variables
%mu = M_earth/M;
mu = 3.040357143 * 10^-6;
omega = 1;

%correct normalization, identification of constants for a system is extremely imp

%Centered at barycentre
r = 1;
x_1 = -mu;
x_2 = 1- mu;

L_1 = 0.989986;    
x0 = 0.9888383910739;
y0 = 0;
z0 = 0.0008152222855;
x0_dot = 0;
y0_dot =0.0089606022073; 
z0_dot = 0;

pos_sun = [x_1, 0, 0];
pos_earth = [x_2, 0, 0];


for repeat=1:50
    
    %define tstart, tend and the number of time steps you want to take, n:
    tstart = 0;
    tend = 3.5;
    n = 0.00001;
    tspan = linspace(tstart,tend,tend/n);
    %define the initial conditions making sure to use the right ordering
    xinit = [x0, x0_dot,  y0, y0_dot, z0, z0_dot];
    
    OPTIONS = odeset('RelTol',3e-14,'AbsTol',1e-16);
    %Get x. I have left options as the default value so it’s not included in the function call below
    [t,x] = ode113(@ode_orbitdynamics_nonlinear, tspan, xinit, OPTIONS);
    %define the output variables:
    X = x(:,1);
    X_dot = x(:,2);
    Y = x(:,3);
    Y_dot = x(:,4);
    Z = x(:,5);
    Z_dot = x(:,6);

    index = find(Y < 0, 1, 'first');
    T_half = n*(index-1); %check


    v_r = [X, Y, Z];
    r_1 = v_r - pos_sun;
    r_2 = v_r - pos_earth;
    r1 = (r_1(:,1).^2 + r_1(:,2).^2 + r_1(:,3).^2).^0.5 ;
    r2 = (r_2(:,1).^2 + r_2(:,2).^2 + r_2(:,3).^2).^0.5 ;
    

    
    N=6;
    phi_int(1:N^2) 	   = reshape(eye(N),N^2,1);
    phi_int(1+N^2:N+N^2) = xinit;
    h = 0.00001;
    tspan1 = linspace(0,index*h,index);
    
    RelTol = 3.e-14 ; 
    AbsTol = 1e-16;
    OPTIONS = odeset('RelTol',RelTol,'AbsTol',AbsTol);

    [t,phi] = ode113(@varEq,tspan1,phi_int, OPTIONS);
    
    phi = phi(index,:);
    phi  = reshape(phi(1:36),6,6);
    
    
    x_half = [X(index),X_dot(index), Y(index), Y_dot(index), Z(index), Z_dot(index)];
    x_desired = [0, 0]; %x_dot z_dot
    del_xf = x_desired - [x_half(2), x_half(6)]; %error in final state
   
    xdotdot = X(index) + 2*Y_dot(index) - (1-mu)*(X(index) +mu)/r1(index)^3 - mu*(X(index) - (1-mu))/r2(index)^3 ; %small formula mistake
    zdotdot = -(1-mu)*Z(index)/r1(index)^3 - mu*Z(index)/r2(index)^3 ;

    PHI11 = phi(4,1) - phi(2,1)*(xdotdot/Y_dot(index));
    PHI12 = phi(4,5) - phi(2,5)*(xdotdot/Y_dot(index));
    PHI21 = phi(6,1) - phi(2,1)*(zdotdot/Y_dot(index));
    PHI22 = phi(6,5) - phi(2,5)*(zdotdot/Y_dot(index));
    PHI = [PHI11  PHI12;
           PHI21  PHI22];
 
    
    PHI_inv = inv(PHI);

    del_x0 = PHI_inv*del_xf.'; %del_x0 = [del_x0, del, del_y_dot]
    
    %changing initial state
    x0 = x0 + del_x0(1);
    y0_dot = y0_dot + del_x0(2);

end

tstart = 0;
tend = 6;
n = 0.00001;
tspan = linspace(tstart,tend,tend/n);
%define the initial conditions making sure to use the right ordering
xinit = [x0, x0_dot,  y0, y0_dot, z0, z0_dot];

OPTIONS = odeset('RelTol',3e-14,'AbsTol',1e-14);
%Get x. I have left options as the default value so it’s not included in the function call below
[t,x] = ode113(@ode_orbitdynamics_nonlinear, tspan, xinit,OPTIONS);
%define the output variables:
X = x(:,1);
X_dot = x(:,2);
Y = x(:,3);
Y_dot = x(:,4);
Z = x(:,5);
Z_dot = x(:,6);


Vel = (X_dot.^2 + Y_dot.^2 + Z_dot.^2).^0.5;

figure();
plot(t,Vel); % use this is you are plotting a 3D plot title(’ ’)
xlabel('time'); 
ylabel('velocity'); 
grid on;


figure();
plot3(X,Y,Z); % use this is you are plotting a 3D plot title(’ ’)
xlabel('X Label'); 
ylabel('Y Label'); 
zlabel('Z label'); 
grid on;

xlim([0.9 1.1]);
hold on
%scatter(pos_sun,0,1000, 'filled','r');
scatter(pos_earth,0,50, 'filled','b');
scatter(L_1,0,10, 'filled','g');
hold off

