function ret=vertical_plane_model(t,y)


%% Import glider vairables

load('glider_vars.mat');

%% Define terms

v1 = y(1);                      % v1
v2 = 0;                         % v2
v3 = y(2);                      % v3

theta = y(3);                   % Pitch

Omega1 = 0;                     % Angular velocity about e1 axis
Omega2 = y(4);                  % Angular velocity about e2 axis
Omega3 = 0;                     % Angular velocity about e3 axis

rp1 = y(5);                     % rp1
rp2 = 0;                        % rp2
rp3 = y(6);                     % rp3

rb1 = 0;                        % rb1
rb2 = 0;                        % rb2
rb3 = 0;                        % rb3

Pp1 = y(7);                     % Pp1
Pp2 = 0;                        % Pp2
Pp3 = y(8);                     % Pp3

x = y(9);                       % x-position in earth frame
z = y(10);                      % z-position (depth) in earth frame

mb = y(11);                     % Ballast mass

alpha = atan(v3/v1);            % Angle of attack

%% Check on ballast mass

if glide_dir == 'U'
    if mb <= mb_d
        mb = mb_d;              % Turn off ballast when desired mass is obtained
    end    
elseif glide_dir == 'D'
    if mb >= mb_d
        mb = mb_d;              % Turn off ballast when desired mass is obtained
    end  
end

%% Position Vectors

% Limit the longitudinal position of primary moving mass
if glide_dir == 'U'             % During upward glide primary moving mass moves in negative e1 direction
    if rp1 <= rp1_d
        rp1 = rp1_d;
    end
elseif glide_dir == 'D'         % During downward glide primary moving mass moves in positive e1 direction
    if rp1 >= rp1_d
        rp1 = rp1_d;
    end    
end

rp = [rp1 rp2 rp3]';
rp_c = [ 0   -rp3  rp2
         rp3  0   -rp1
        -rp2  rp1  0  ];

rb = [rb1 rb2 rb3]';
rb_c = [ 0   -rb3  rb2
         rb3  0   -rb1
        -rb2  rb1  0  ]; 
    
%% Vehicle Velocity Vectors

v = [v1 v2 v3]';                    % Linear velocity
Omega = [Omega1 Omega2 Omega3]';    % Angular velocity

%% Momentum Matrices

% Pp = mbar*(v + cross(Omega,rp) + rp_dot);
% Pb = mb*(v + cross(Omega,rb) + rb_dot);

Pp = [Pp1 Pp2 Pp3]';                % Linear momentum of movable point mass
Pb = mb*v;                          % Linear momentum of ballast mass

%% Transformation Matrix

RT = [cos(theta)        0       -sin(theta)
      0                 1       0
      sin(theta)        0       cos(theta)];

%% Force and Torque terms

L = (KL0 + KL*((alpha)))*(v1^2 + v3^2);                      % Lift force                         
D = (KD0 + KD*((alpha))^2)*(v1^2 + v3^2);                    % Drag force
MDL = (KM0 + KM*((alpha)))*(v1^2 + v3^2) ...                 % Moment
      + KOmega1_2*Omega2 + KOmega2_2*(Omega2^2);

Fext = [-D 0 -L]';                                           % Total external force
Text = [0 MDL 0]';                                           % Total external mmoment

%% Internal Mass Velocity Vectors

% Limit the rate of movement of primary moving mass
if glide_dir == 'U'
    rp1_dot = (1/mbar)*Pp1 - v1 - rp3*Omega2;                % rp1_dot in upward glide
    if rp1 <= rp1_d                                          % When desired value reached, set rate = 0
        rp1_dot = 0;
    end
elseif glide_dir == 'D'
    rp1_dot = (1/mbar)*Pp1 + v1 - rp3*Omega2;                % rp1_dot in downward glide
    if rp1 >= rp1_d                                          % When desired value reached, set rate = 0
        rp1_dot = 0;
    end
end

rp2_dot = 0;                                                 % rp2_dot
rp3_dot = 0;                                                 % rp3_dot

rp_dot = [rp1_dot rp2_dot rp3_dot]';                          % rp_dot
rb_dot = [0 0 0]';                                            % rb_dot

%% Control Inputs

% Set acceleration input of primary moving mass
if glide_dir == 'U'
    wp1 = -0.02;            % wp1 in upward glide (m/s^2)
elseif glide_dir == 'D' 
    wp1 = 0.02;             % wp1 in downward glide (m/s^2)
end
wp2 = 0;                    % wp2
wp3 = 0;                    % wp3

wp = [wp1                   % (w1) Accleration input of movable point mass mbar along body axis e1
      wp2                   % (w2) Accleration input of movable point mass mbar along body axis e2
      wp3];                 % (w3) Accleration input of movable point mass mbar along body axis e3

% Accleration input of ballast point mass mb in body coordinates
wb = 0;                     

% Ballast rate (kg/s)
u4 = ballast_rate;          

% Check on ballast rate
if glide_dir == 'U'
    if mb <= mb_d
        u4 = 0;             % Turn off ballast when desired mass is obtained
    end    
elseif glide_dir == 'D'
    if mb >= mb_d
        u4 = 0;             % Turn off ballast when desired mass is obtained
    end  
end

%% Calculate net buoyancy

m0 = mh + mw + mb + mbar - m;

 %% Control Transformation

Minv = inv(M);
Jinv = inv(J);

F = [ (Minv - rp_c*Jinv*rp_c + (I3/mbar))       (Minv - rp_c*Jinv*rb_c)               
      (Minv - rb_c*Jinv*rp_c)                   (Minv - rb_c*Jinv*rb_c + (I3/mb))];    

H = inv(F);

Zp = -Minv*[cross((M*v + Pp + Pb),Omega) + m0*g*RT*k + Fext] - cross(Omega,rp_dot) ...
      -Jinv*cross([cross((J*Omega + rp_c*Pp + rb_c*Pb),Omega) + cross((M*v),v) + Text ...
      + cross(cross(Omega,rp),Pp) + cross(cross(Omega,rb),Pb) ...
      + (mbar*rp_c + mb*rb_c)*g*RT*k],rp);

Zb = -Minv*[cross((M*v + Pp + Pb),Omega) + m0*g*RT*k + Fext] - cross(Omega,rb_dot) ...
     -Jinv*cross([cross((J*Omega + rp_c*Pp + rb_c*Pb),Omega) + cross((M*v),v) + Text ...
     + cross(cross(Omega,rp),Pp) + cross(cross(Omega,rb),Pb) ...
     + (mbar*rp_c + mb*rb_c)*g*RT*k],rb);
  
 u = H*[(-Zp + wp)
        (-Zb + wb)];
 
 %% Control Inputs
 
 u1 = u(1);
 u3 = 0; 
 
 if glide_dir == 'U'
    if rp1 <= rp1_d
        u1 = 0;
    end  
 elseif glide_dir == 'D'
    if rp1 >= rp1_d
        u1 = 0;
    end
 end
 
%% Equations of Motion

x_dot       = v1*cos(theta) + v3*sin(theta);
z_dot       = -v1*sin(theta) + v3*cos(theta);
theta_dot   = Omega2;
Omega2_dot  = (1/j2)*((m3-m1)*v1*v3 - (rp1*Pp1 + rp3*Pp3)*Omega2 - mbar*g*(rp1*cos(theta) + rp3*sin(theta)) + MDL - rp3*u1 + rp1*u3);
v1_dot      = (1/m1)*(-m3*v3*Omega2 - Pp3*Omega2 - m0*g*sin(theta) + L*sin(alpha) - D*cos(alpha) - u1);
v3_dot      = (1/m3)*(m1*v1*Omega2 + Pp1*Omega2 + m0*g*cos(theta) - L*cos(alpha) - D*sin(alpha) - u3);
Pp1_dot     = u1;
Pp3_dot     = u3;
mb_dot      = u4;

%% return
ret = [v1_dot v3_dot theta_dot Omega2_dot rp1_dot rp3_dot Pp1_dot Pp3_dot x_dot z_dot mb_dot]';
 
%% Miscellaneous

% rw1 = 0;
% rw2 = 0;
% rw3 = 0;
% rw = [rw1;rw2;rw3];
% rw_c = [ 0   -rw3  rw2
%          rw3  0   -rw1
%         -rw2  rw1  0  ]; 
   
   
% I = [ (M + (mbar + mb + mw)*I3)   (-mbar*rp_c - mb*rb_c)               mbar*I3      mb*I3       mw*I3
%       (mbar*rp_c + mb*rb_c)       (J - mbar*rp_c*rp_c - mb*rb_c*rb_c)  mbar*rp_c    mb*rb_c     mw*rw_c
%       mbar*I3                     -mbar*rp_c                           mbar*I3      Z3          Z3
%       mb*I3                       -mb*rb_c                             Z3           mb*I3       Z3
%       mw*I3                       -mw*rw_c                             Z3           Z3          mw*I3];

% Iinv = inv(I);
