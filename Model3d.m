function ret=Model3d(t,y)


%% Import glider vairables

load('glider_vars.mat');

%% Define terms

b       = y(1:3);                 % Position in earth frame
Omega   = y(4:6);                 % Angular rate
v       = y(7:9);                 % Velocity
rp      = y(10:12);               % Position of primary moving mass
rs      = y(13:15);               % Position of secondary moving mass
rb      = y(16:18);               % Position of ballast mass
Pp      = y(19:21);               % Linear momentum of primary moving mass
Ps      = y(22:24);               % Linear momentum of secondary moving mass
%Pb     = y(25:27);               % Linear momentum of ballast moving mass
mb      = y(28);                  % Ballast mass
theta   = y(29);                  % Pitch
alpha = atan(v(3)/v(1));          % Angle of attack

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

% Limit the longitudinal position of moving masses
if glide_dir == 'U'             % During upward glide primary moving mass moves in negative e1 direction
    if rp(1) <= rp1_d
        rp(1) = rp1_d;
    end
    if rs(1) <= rs1_d
        rs(1) = rs1_d;
    end
elseif glide_dir == 'D'         % During downward glide primary moving mass moves in positive e1 direction
    if rp(1) >= rp1_d
        rp(1) = rp1_d;
    end 
    if rs(1) >= rs1_d
        rs(1) = rs1_d;
    end 
end

% Primary moving mass
rp_c = [ 0      -rp(3)  rp(2)
         rp(3)  0       -rp(1)
        -rp(2)  rp(1)   0  ];
    
% Secondary moving mass
rs_c = [ 0      -rs(3)  rs(2)
         rs(3)  0       -rs(1)
        -rs(2)  rs(1)   0  ];

% Ballast mass
rb_c = [ 0      -rb(3)  rb(2)
         rb(3)  0       -rb(1)
        -rb(2)  rb(1)   0  ]; 
    
    
%% Transformation matrix

R = [cos(yaw)*cos(theta)    (-sin(yaw)*cos(phi) + cos(yaw)*sin(theta)*sin(phi)) (sin(yaw)*sin(phi) + cos(yaw)*cos(phi)*sin(theta))
     sin(yaw)*cos(theta)    (cos(yaw)*cos(phi) + sin(phi)*sin(theta)*sin(yaw))  (-cos(yaw)*sin(phi) + sin(theta)*sin(yaw)*cos(phi))
     -sin(theta)            cos(theta)*sin(phi)                                 cos(theta)*cos(phi)];

RT = R';

%% Force and Torque terms

L = (KL0 + KL*((alpha)))*(v(1)^2 + v(3)^2);                      % Lift force                         
D = (KD0 + KD*((alpha))^2)*(v(1)^2 + v(3)^2);                    % Drag force
MDL = (KM0 + KM*((alpha)))*(v(1)^2 + v(3)^2) ...                 % Moment
      + KOmega1_2*Omega(2) + KOmega2_2*(Omega(2)^2);

F_ext = [-D 0 -L]';                                          % Total external force
T_ext = [0 MDL 0]';                                          % Total external mmoment

%% Internal Mass Velocity Vectors

% Limit the rate of movement of primary moving mass
if glide_dir == 'U'
    rp_dot = (1/mbar)*Pp - v - cross(Omega, rp);             % rp_dot in upward glide
    if ms > 0                                                
       rs_dot = (1/ms)*Ps - v - cross(Omega, rs);            % rs_dot in upward glide
    else
       rs_dot = [0 0 0]'; 
    end
    
    if rp(1) <= rp1_d                                        % When desired value reached, set rate = 0
        rp_dot = [0 0 0]';
    end
    if rs(1) <= rs1_d                                        % When desired value reached, set rate = 0
        rs_dot = [0 0 0]';
    end
elseif glide_dir == 'D'
    rp_dot = (1/mbar)*Pp + v - cross(Omega, rp);             % rp_dot in upward glide
    if ms > 0                                                
       rs_dot = (1/ms)*Ps + v - cross(Omega, rs);            % rs_dot in upward glide
    else
       rs_dot = [0 0 0]'; 
    end
    
    if rp(1) >= rp1_d                                        % When desired value reached, set rate = 0
        rp_dot = [0 0 0]';
    end
    if rs(1) >= rs1_d                                        % When desired value reached, set rate = 0
        rs_dot = [0 0 0]';
    end
end

rb_dot = [0 0 0]';                                           % rb_dot

%% Momentum Matrix

Pb = mb*(v + cross(Omega,rb) + rb_dot);                     % Linear momentum of ballast mass

%% Calculate net buoyancy

m0 = mh + mw + mb + ms + mbar - m;

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
 
% Check on ballast rate
if glide_dir == 'U'
    if mb <= mb_d
        u_ballast_rate = 0;             % Turn off ballast when desired mass is obtained
    end    
elseif glide_dir == 'D'
    if mb >= mb_d
        u_ballast_rate = 0;             % Turn off ballast when desired mass is obtained
    end  
end

%% Control Transformation

F = [ (Minv - rp_c*Jinv*rp_c + (I3/mbar))       (Minv - rp_c*Jinv*rb_c)               
      (Minv - rb_c*Jinv*rp_c)                   (Minv - rb_c*Jinv*rb_c + (I3/mb))];    

H = inv(F);

Zp = -Minv*[cross((M*v + Pp + Pb),Omega) + m0*g*RT*k + F_ext] - cross(Omega,rp_dot) ...
      -Jinv*cross([cross((J*Omega + rp_c*Pp + rb_c*Pb),Omega) + cross((M*v),v) + T_ext ...
      + cross(cross(Omega,rp),Pp) + cross(cross(Omega,rb),Pb) ...
      + (mbar*rp_c + mb*rb_c)*g*RT*k],rp);

Zb = -Minv*[cross((M*v + Pp + Pb),Omega) + m0*g*RT*k + F_ext] - cross(Omega,rb_dot) ...
     -Jinv*cross([cross((J*Omega + rp_c*Pp + rb_c*Pb),Omega) + cross((M*v),v) + T_ext ...
     + cross(cross(Omega,rp),Pp) + cross(cross(Omega,rb),Pb) ...
     + (mbar*rp_c + mb*rb_c)*g*RT*k],rb);
  
 u = H*[(-Zp + wp)
        (-Zb + wb)];
    
 %% Control Inputs
 
 if glide_dir == 'U'
     
    if rp(1) <= rp1_d
        u_bar = [0 0 0]';
    else
        u_bar = [-0.18 0 0]';
        %u_bar = [u(1) 0 0]';
    end 
    
    if rs(1) <= rs1_d
        us    = [0 0 0]';
    else
        us    = [-0.18 0 0]';
    end
     
    if ms == 0                              % If secondary mass doed not exist
        us    = [0 0 0]';
    end
    
 elseif glide_dir == 'D'
    
    if rp(1) >= rp1_d
        u_bar = [0 0 0]';
    else
        u_bar = [0.18 0 0]';
        %u_bar = [u(1) 0 0]';
    end 
    
    if rs(1) >= rs1_d
        us    = [0 0 0]';
    else
        us    = [0.18 0 0]';
    end
    
    if ms == 0                              % If secondary mass doed not exist
        us    = [0 0 0]';
    end
 end
 
 ub = [0 0 0]';

%% Equations of motion

T_bar      = cross((J*Omega + rp_c*Pp + rs_c*Ps + rb_c*Pb), Omega) + cross((M*v),v) ...
             + cross(cross(Omega, rp), Pp) + + cross(cross(Omega, rs), Ps) + + cross(cross(Omega, rb), Pb) + (mbar*rp_c + ms*rs_c + mb*rb_c)*g*RT*k ...
             + T_ext - rp_c*u_bar - rs_c*us - rb_c*ub;

F_bar      = cross((M*v + Pp + Ps + Pb), Omega) + m0*g*RT*k + F_ext - u_bar - us - ub;

b_dot      = R*v;
Omega_dot  = Jinv*T_bar;
v_dot      = Minv*F_bar;
% rp_dot     = (1/mbar)*Pp - v - cross(Omega, rp);
% rs_dot     = (1/ms)*Pd - v - cross(Omega, rs);
% rb_dot     = (1/mb)*Pb - v - cross(Omega, rb);
Pp_dot     = u_bar;
Ps_dot     = us;
Pb_dot     = ub;
mb_dot     = u_ballast_rate;
theta_dot  = Omega(2);

%% Return

ret = [b_dot' Omega_dot' v_dot' rp_dot' rs_dot' rb_dot' Pp_dot' Ps_dot' Pb_dot' mb_dot theta_dot]';
