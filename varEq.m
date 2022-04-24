function PHIdot = varEq(t,phi)

mu = 3.040357143 * 10^-6 ; 

x(1:6) = phi(37:42);
phi  = reshape(phi(1:36),6,6);

X=x(1);
Xdot = x(2);
Y=x(3);
Ydot = x(4);
Z=x(5);
Zdot = x(6);


pos_sun = [-mu, 0, 0];
pos_earth = [1-mu, 0, 0];

v_r = [X, Y, Z];
r1 = norm(v_r - pos_sun);
r2 = norm(v_r - pos_earth);

% The following are three double partial derivatives of the effective potential U(x,y)

Uxx = 1 - (1-mu)*(1 - 3*((X + mu)^2)/(r1^2))/(r1^3) - mu*(1 - 3*((X - (1-mu))^2/(r2^2))/(r2^3)) ; %small formula mistake here
Uyy = 1 - (1-mu)*(1 - 3*(Y^2)/(r1^2))/(r1^3) - mu*(1 - 3*(Y^2)/(r2^2))/(r2^3) ; 
Uzz = -(1-mu)*(1 - 3*Z^2/r1^2)/r1^3 - mu*(1 - 3*Z^2/r2^2)/r2^3 ;  %small formula mistake here
Uxy = 3*Y*(1 - mu)*(X + mu)/r1^5 + 3*mu*Y*(X - (1-mu))/r2^5 ; 
Uyx = Uxy;
Uxz = 3*Z*(1 - mu)*(X + mu)/r1^5 + 3*mu*Z*(X - (1-mu))/r2^5 ;
Uzx = Uxz;
Uyz = 3*Z*(1 - mu)*Y/r1^5 + 3*mu*Z*Y/r2^5 ;
Uzy = Uyz; 


% The following is the Jacobian matrix

        A =   [0     0      0    1  0  0;
               0     0      0    0  1  0;
               0     0      0    0  0  1;
               Uxx   Uxy   Uxz   0  2  0;
               Uyx   Uyy   Uyz  -2  0  0;
               Uzx   Uzy   Uzz   0  0  0];


phidot = A * phi; % variational equation

PHIdot        = zeros(42,1);
PHIdot(1:36)  = reshape(phidot,36,1);
PHIdot(37)    = x(2);
PHIdot(38)    = 2*Ydot + (1 - (1-mu)/r1^3 - mu/r2^3)* X + mu*(1-mu)/r2^3 - mu*(1-mu)/r1^3;
PHIdot(39)    = x(4);
PHIdot(40)    = -2*Xdot + (1 - (1-mu)/r1^3 - mu/r2^3)* Y;
PHIdot(41)    = x(6);
PHIdot(42)    = (- (1-mu)/r1^3 - mu/r2^3)* Z;


end
