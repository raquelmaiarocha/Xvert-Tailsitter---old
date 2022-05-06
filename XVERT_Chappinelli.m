% X-Vert Model - Romain Chiappinelli
clear all
clc


%% TO DO's %%
% https://www.mathworks.com/help/matlab/matlab_prog/add-reminders-to-files.html

% wind model (Vwb) --> Appendix B

% DONE
% use slipstream velocity in aerodynamics (wing and rod)
% use slipstream velocity in thrusters (inflow velocity)



%% NOMENCLATURE & DATA %%

% AIRCRAFT DATA
xvert = aircraftdata;

% AIR DATA
global rho Vb Wb Vwb q R
rho = 1.2; % [kg/m^3]

% Vb = [u, v, w]'  --> velocity expressed in the body frame [m/s]
Vb = [10;0;2];

% Wb = [p, q, r]'  --> body angular rates in the body frame [rad/s]
Wb = [0;0;0];

% wind velocity in the body frame [m/s] (todo)
Vwb = [0;0;0];              


%% ROTATION MATRIX

% p = [N,E,D]'   --> inertial position
% Φ = [ϕ, θ, ψ]'  --> orientation 
roll = 0; % degree
pitch = 10; % degree
yaw = 0; % degree

% q = [q0, q1, q2, q3]'  --> orientation
q = quaternion([yaw pitch roll],"eulerd","zyx","frame");

R = rotmat(q,"frame");
% R     --> ned2body
% R'    --> body2ned

%% GRAVITY FORCE Fg
FgNED = [0 0 xvert.m*9.8]; % [m/s^2]
Fg = R*FgNED'; % Body frame
Mg = [0,0,0]';

%% CONTROL STATE SPACE
% sr = 0; % elevon right signal
% sl = 0; % elevon left signal
% 
% der = sr*xvert.dfmax; % right elevon defletion (linear model) --> eq 2.51 Tese 0
% del = sl*xvert.dfmax; % left elevon defletion (linear model) --> eq 2.51 Tese 0
% 
% dthr = 0.7; % throttle right signal
% dthl = 0.7; % throttle left signal

% de = (der + del)/2; % elevator (pitch)
% da = (der - del)/2; % aileron (roll)
% dt = (dthr + dthl)/2; % thrusters (thrust)
% dr = (dthr - dthl)/2; % rudder (yaw)

%% AIRCRAFT EQUATIONS OF MOTION

dp = R'*Vb;                   % --> inertial velocity --> inertial position 
dq = 0.5 * [4x4] * q             % --> quaternion orientation rates --> quaternion orientation --> LER APENDICE B (Small Unmaned ... Beard)

dVb = Fb/xvert.m - Wb*Vb;        % --> acceleration -> velocity expressed in the body frame [m/s]
dWb = xvert.I_1 * (Mb-cross(Wb,xvert.I*Wb));  % --> body angular rates derivative --> body angular rates in the body frame


%% CENARIOS

cases = [% sr     sl     dthr   dthl    description
           0.3    0.3    0.8    0.8  ;   %'forward flight'
           0.5    -0.5   0.7    0.7  ;   %'roll right'
           0.5    0.5    0.7    0.7  ;   %'pitch up'
           0      0      0.5    0.9  ;   %'yaw left'
           0      0      0.9    0.9      %'increase velocity'
    ];

names = {'forward flight'; 'roll right'; 'pitch up'; 'yaw left'; 'increase velocity'};

%% BODY FORCES AND MOMENTS

for i=1:5

    [sr, sl, dthr, dthl] = matsplit(cases(i,:),1);

    der = sr*xvert.dfmax;
    del = sl*xvert.dfmax;

    [Fth, Mth, Tr, Tl] = thrusters(Vb, Wb, Vwb, rho, xvert, dthr, dthl);
    
    [Fa, Ma] = wing_aerodynamics(xvert,Vb, Wb, Vwb, rho, der, del, Tr, Tl);
    [Fr, Mr] = structural_aerodynamics(Vb, Wb, Vwb, rho, xvert, Tr, Tl);
    [Fc, Mc] = ground_contact(Vb, Wb, Vwb, rho, xvert, R);
    
    
    Fb = Fg + Fth + Fa + Fr + Fc;      % net body forces exerted at the CM [N]
    Mb = Mth + Ma + Mr + Mc;           % net body moments exerted at the CM [Nm]
    
    % body moments Mb = [L, M, N]' and the forward force F
    [L, M, N, F] = simplified_model(xvert, Vb, Wb, Vwb, rho, dthl, dthr, del, der);

    % print results
    name = names(i)

    %forces  = [Fg, Fth, Fa, Fr, Fc, Fb, [F,0,0]']
    %moments = [Mth, Ma, Mr, Mc, Mb, [L,M,N]']
    
    %names = {'Thrust'; 'Aerodynamics'; 'Drag'; 'Ground'; 'Total'; 'Simplified'}
    FT = table(Fg, Fth, Fa, Fr, Fc, Fb, [F,0,0]')
    MT = table(Mg, Mth, Ma, Mr, Mc, Mb, [L,M,N]')
end

% PLOT GRAPHS
plot_aerdynamics()
plot_propeller(xvert)


%% SIMPLIFIED ANALYTICAL AERODYNAMIC MODEL (pag 40 Tese 0)
function [L, M, N, F] = simplified_model(xvert, Vb, Wb, Vwb, rho, dthl, dthr, del, der)
    
    l = 0.145; % [m] thrusters lateral position

    [u, ~, w, ~, aa, ~] = airspeed(Vb, Wb, Vwb, xvert.CM, xvert);

    [Tr, Qr, ~] = propeller(u,rho,xvert,dthr);
    [Tl, Ql, ~] = propeller(u,rho,xvert,dthl);

    Cm = (5.18e-4*aa^5 - 1.03e-3*aa^3 + 2.72e-2*aa)*(aa-pi)*(aa+pi); % eq D.2 Tese 0
    km = 0.8; %  gain used in the controller
    M0 = km*0.5*rho*(u^2+w^2)*xvert.Sw*xvert.cw*Cm; % eq D.1 Tese 0
    
    % table H.1 Tese 0
    cx = 9.91e-4; % [m^3/rad]
    cy = 4.75e-4; % [m^3/rad]
    bx = 9.37e-4; % [m^3/rad]
    by = 3.48e-4; % [m^3/rad]

    Pd = 0.5*rho+(u^2+w^2); % XZ dynamic pressure

    % eq 2.41 Tese 0
    L = cx/(pi*xvert.rp^2)*(del*Tl-der*Tr) +  Pd*bx*(del-der) + (Qr-Ql);
    M = -cy/(pi*xvert.rp^2)*(del*Tl+der*Tr) - Pd*(cy+by)*(del+der) + M0;
    N = l*(Tl-Tr);
    F = Tl+Tr;

end

%% WING AERODYNAMICS (Fa Ma)
function [Fa, Ma] = wing_aerodynamics(xvert, Vb, Wb, Vwb, rho, der, del, Tr, Tl)
    
    % Realocate variables
    [Si, bi, ci, cfi] = matsplit(xvert.wing,1);
    riG = xvert.wing(:,[5:7]);

    Fi = zeros(11,3);
    Mi = zeros(11,3);
    Mi_ac = zeros(11,3);
    
    % Horizontal Wing Segments
    for i=1:9
        [ui, vi, wi, Vi, aai, bbi] = airspeed(Vb, Wb, Vwb, [riG(i,:)]', xvert);
        ui = slipstream(Vb, Wb, Vwb, rho, xvert, [riG(i,:)]', Tr, Tl);
        
        switch i
            case {1, 5, 9}
                [CLwi,CDwi,Cmwi] = aerodynamics(aai,0);
            case {2, 3, 4}
                [CLwi,CDwi,Cmwi] = aerodynamics(aai,del);
            case {6, 7, 8}
                [CLwi,CDwi,Cmwi] = aerodynamics(aai,der);
        end

        Fi(i,1) = 0.5*rho*bi(i)*ci(i) * (ui^2+wi^2) * (CLwi*sin(aai)-CDwi*cos(aai));
        Fi(i,3) = 0.5*rho*bi(i)*ci(i) * (ui^2+wi^2) * (-CLwi*cos(aai)-CDwi*sin(aai));

        Mi(i,2) = 0.5*rho*bi(i)*ci(i)^2 * (ui^2+wi^2) * Cmwi;
        
        Mi_ac_aux = Mi(i,:)' + cross(riG(i,:)',Fi(i,:)');
        [Mi_ac(i,1), Mi_ac(i,2), Mi_ac(i,3)] = matsplit(Mi_ac_aux);
    end
    
    % Vertical Wing Segments
    for i=10:11
        [ui, vi, wi, Vi, aai, bbi] = airspeed(Vb, Wb, [riG(i,:)]', Vwb, xvert);
        ui = slipstream(Vb, Wb, Vwb, rho, xvert, [riG(i,:)]', Tr, Tl);

        [CLwi,CDwi,Cmwi] = aerodynamics(bbi,0);

        Fi(i,1) = 0.5*rho*bi(i)*ci(i) * (ui^2+vi^2) * (CLwi*sin(bbi)-CDwi*cos(bbi));
        Fi(i,2) = 0.5*rho*bi(i)*ci(i) * (ui^2+vi^2) * (-CLwi*cos(bbi)-CDwi*sin(bbi));

        Mi(i,3) = 0.5*rho*bi(i)*ci(i)^2 * (ui^2+vi^2) * Cmwi;

        Mi_ac_aux = Mi(i,:)' + cross(riG(i,:)',Fi(i,:)');
        [Mi_ac(i,1), Mi_ac(i,2), Mi_ac(i,3)] = matsplit(Mi_ac_aux);
    end
    
    Fa = sum(Fi)'; % Sum of Matrix Rows
    Ma = sum(Mi_ac)'; % Sum of Matrix Rows
end

%% STRUCTURAL SEGMENTS AERODYNAMICS (Fr Mr)
function [Fr, Mr] = structural_aerodynamics(Vb, Wb, Vwb, rho, xvert, Tr, Tl)

    if Vb == [0;0;0]
        % set rod drag to zero when inflow velocity is zero --> pag 33 Tese 0
        Fr = [0;0;0];
        Mr = [0;0;0];
    else
        % 21 right side segments coordinates
        %   - 12 propeller protector circle [3 mm]
        %   - 3 propeller protector radial [3 mm]
        %   - 6 landing gear rod seg [7 mm]
        
        rc = [%r1Gx r1Gy∗ r1Gz r2Gx r2Gy r2Gz dj
                46.5 252.0 0.0 19.0 252.0 -66.0 7
                19.0 252.0 -66.0 -10.0 252.0 -66.0 7
                -10.0 252.0 -66.0 0.0 252.0 25.0 7
                0.0 252.0 25.0 -10.0 252.0 69.0 7
                -10.0 252.0 69.0 14.5 252.0 69.0 7
                14.5 252.0 69.0 81.0 252.0 0.0 7
                177.0 207.8 36.2 177.0 182.3 62.8 3
                177.0 181.3 62.8 177.0 146.0 72.5 3
                177.0 145.0 72.5 177.0 109.8 62.8 3
                177.0 108.8 62.8 177.0 83.2 36.2 3
                177.0 82.2 36.2 177.0 73.5 0.0 3
                177.0 72.5 0.0 177.0 83.2 -36.2 3
                177.0 82.2 -36.2 177.0 109.7 -62.8 3
                177.0 108.7 -62.8 177.0 146.0 -72.5 3
                177.0 145.0 -72.5 177.0 182.3 -62.8 3
                177.0 181.3 -62.8 177.0 208.8 -36.3 3
                177.0 207.8 -36.3 177.0 218.5 -0.0 3
                177.0 217.5 -0.0 177.0 208.8 36.3 3
                172.0 140.5 9.5 175.0 114.5 54.6 3
                172.0 140.5 -9.5 175.0 114.5 -54.6 3
                172.0 157.0 -0.0 175.0 209.0 -0.0 3
        ]*10^-3;
        
        % 21 right side segments coordinates
        lc = [rc(:,1) -rc(:,2) rc(:,3) rc(:,4) -rc(:,5) rc(:,6) rc(:,7)];
        
        % coordinates, lenght, diameter and center position of all segments
        c = [rc; lc];
        l = zeros(length(c),3);
        ln = zeros(length(c),1);
        d = zeros(length(c),1);
        r = zeros(length(c),3);
        v = zeros(length(c),3);
        vn = zeros(length(c),1);
        u = zeros(length(c),3);
        tt = zeros(length(c),1);
        N = zeros(length(c),1);
        
        for j=1:length(c)
            
            % vector length (lx, ly, lz)
            l(j,1) = c(j,4)-c(j,1);
            l(j,2) = c(j,5)-c(j,2);
            l(j,3) = c(j,6)-c(j,3);
    
            ln(j) = norm(l(j,:)); % norm of the vector
            
            % diameter
            d(j)   = c(j,7);
            
            % center position
            r(j,1) = (c(j,4)+c(j,1))/2;
            r(j,2) = (c(j,5)+c(j,2))/2;
            r(j,3) = (c(j,6)+c(j,3))/2;
    
            % velocity
            [v(j,1), v(j,2), v(j,3), vn(j), ~, ~] = airspeed(Vb, Wb, Vwb, r(j,:)', xvert);
            v(j,1) = slipstream(Vb, Wb, Vwb, rho, xvert, [r(j,:)]', Tr, Tl);
    
            vn(j) = norm(v(j,:)); % norm of the vector
    
            % unit vector
            if dot(v(j,:)',cross(cross(v(j,:)',l(j,:)'),l(j,:)')) >= 0
                %[u(j,1), u(j,2), u(j,3)]' = cross(cross(v(j,:)',l(j,:)'),l(j,:)')/(abs(vn(j))*ln(j)^2);
                u(j,:) = -cross(cross(v(j,:)',l(j,:)'),l(j,:)')/(abs(vn(j))*ln(j)^2);
            else 
                %[u(j,1); u(j,2); u(j,3)]' = -cross(cross(v(j,:)',l(j,:)'),l(j,:)')/(abs(vn(j))*ln(j)^2);
                u(j,:) = cross(cross(v(j,:)',l(j,:)'),l(j,:)')/(abs(vn(j))*ln(j)^2);     
            end
    
            % theta
            tt(j) = asin(norm(cross(v(j,:)',l(j,:)'))/(abs(vn(j))*abs(ln(j))));
    
            % drag force
            CDrods = 1.1;
            N(j) = 0.5*rho*vn(j)^2*ln(j)*d(j)*CDrods*sin(tt(j))^2;
            
        end
        
        Fr = sum(N.*u)'; % Sum of Matrix Rows
        Mr = sum(cross(r,N.*u))'; % Sum of Matrix Rows
    end
end

%% THRUSTERS FORCES AND MOMENTS (Fth Mth)
function [Fth, Mth, Tr, Tl] = thrusters(Vb,Wb,Vwb,rho,xvert,dthr,dthl)
    
    l = 0.145; % [m] thrusters lateral position
    Ith = 1.626e-6; % [kg/m^2] thruster rotational inertia

    rthr = [0.175; l; 0]; % [m] right propeller position --> Table 2-3 Tese 0
    rthl = [0.175; -l; 0]; % [m] left propeller position --> Table 2-3 Tese 0

    [uir, vir, wir, Vir, aair, bbir] = airspeed(Vb, Wb, Vwb, rthr, xvert);
    [uil, vil, wil, Vil, aail, bbil] = airspeed(Vb, Wb, Vwb, rthl, xvert);

    [Tr, Qr, wr] = propeller(uir,rho,xvert,dthr);
    [Tl, Ql, wl] = propeller(uil,rho,xvert,dthl);

    %dthr_aux = inverse_propeller(Vb,Wb,Vwb,rho,xvert,Tr,rthr)
    %dthl_aux = inverse_propeller(Vb,Wb,Vwb,rho,xvert,Tl,rthl)

    Fth = [Tl+Tr; 0; 0];
    
    Mt = [Qr-Ql; 0; l*(Tl-Tr)]; % moments generated by the thrusters
    Mgyro = Ith*(wl-wr)+[0; -Wb(2); Wb(3)];

    Mth = Mt + Mgyro;
end

function [T, Q, omega] = propeller(ui,rho,xvert,dth)
   
    % ESC wp (Vbatt, tau) --> eq 2.43 Tese 0
    % Simplified empirical model (Nahon)
    omega = xvert.Vbatt^(0.8)*(xvert.wp2*dth^2 + xvert.wp1*dth + xvert.wp0);
    
    % Propeller
    J = (pi*ui)/(omega*xvert.rp);
    
    % Experimental coefficients --> eq 2.47 Tese 0
    CT = (xvert.CT2*(J^2)) + (xvert.CT1*(J)) + (xvert.CT0);
    CP = (xvert.CP2*(J^2)) + (xvert.CP1*(J)) + (xvert.CP0);
    
    T = (4/(pi^2))*rho*(omega^2)*(xvert.rp^4)*CT;
    Q = (4/(pi^3))*rho*(omega^2)*(xvert.rp^5)*CP;

end

%% INVERSE THRUSTER MODEL
function dth = inverse_propeller(Vb,Wb,Vwb,rho,xvert,T,rth)
    
    [ui, ~, ~, ~, ~, ~] = airspeed(Vb, Wb, Vwb, rth, xvert);
    omega = pi/(2*xvert.CT0*xvert.rp)*(sqrt(xvert.CT1^2*ui^2 + 4*xvert.CT0*(T/(4*rho*xvert.rp^2)-xvert.CT2*ui^2))-xvert.CT1*ui);
    
    dth = (-xvert.wp1 + sqrt(xvert.wp1^2-4*xvert.wp2*(xvert.wp0-omega/xvert.Vbatt^0.8)))/(2*xvert.wp2);

end

%% GROUND CONTACT DYNAMICS (Fc Mc)
function [Fc, Mc] = ground_contact(Vb, Wb, Vwb, rho, xvert, R)
    
    rk = [% rGx     rGy     rGz 
            177.0   145.0   72.5
            177.0   217.5   0.0
            177.0   145.0  -72.5
            177.0  -145.0   72.5
            177.0  -217.5   0.0
            177.0  -145.0  -72.5
            249.4   0.0     0.0
            118.0   251.0   0.0
            118.0  -251.0   0.0
            -15.0   251.0   71.0
            -15.0   251.0  -71.0
            -15.0  -251.0   71.0
            -15.0  -251.0  -71.0 
    ]*10^-3; % Reference Points for Contact Force Calculations, given in the geometric reference frame G.

    
    kp = 100; % [1/s^2] --> eq 2.35 Tese 0
    kv = 5; % [1/s] --> eq 2.35 Tese 0

    dk = zeros(13,1); % each point k penetrating the ground at a depth dk
    VIk = zeros(13,3);

    for k=1:13 

        %VIk(k,:) = (R'*(Vb+cross(Wb,rk(k,:)')))'; % Inertial velocity ???
        
        FIk(k,:) = [0, 0, -xvert.m*kp*dk(k)]' - xvert.m*kv*VIk(k); % Inertial force at its respective location
        Fk(k,:) = R*FIk(k,:)'; 
    end

    Fc = sum(Fk)'; % Sum of Matrix Rows
    Mc = sum(cross(rk,Fk))'; % Sum of Matrix Rows
end

%% VELOCITY AND ANGLE OF ATTACK
% Vi = [ui, vi, wi]'    --> velocity relative to the air, at any point i, located at a position ri on the aircraft
% aai                   --> angle of attack
% bbi                   --> sideslip angle

% This equation holds for any point outside of the slipstream
function [ui, vi, wi, Vi, aai, bbi] = airspeed(Vb, Wb, Vwb, rG, xvert)
    ri = rG-xvert.CM;
    Vi_v = Vb + cross(Wb,ri) - Vwb;
    [ui, vi, wi] = matsplit(Vi_v);
    
    Vi = norm([ui, vi, wi]); % [m/s]
    aai = atan2(wi,ui); % [rad]
    bbi = atan2(vi,ui); % [rad]
end

% Inside the slipstream (Momentum Theory)
% function [Vxl, Vxr, Vdxl, Vdxr] = slipstream(Vb, Wb, Vwb, xvert)
% 
%     r_thl = [0.177, 0.145, 0]; % right thruster position [m] --> table 2-3 and H.1
%     r_thr = [0.177, -0.145, 0]; % left thruster position [m] --> table 2-3 and H.1
%    
%     [Vinl,~,~] = airspeed(Vb, Wb, Vwb, r_thl, xvert); % inflow velocity at the left propeller
%     [Vinr,~,~] = airspeed(Vb, Wb, Vwb, r_thr, xvert); % inflow velocity at the right propeller
%     
%     Vxl = sqrt(Vinl^2 + 2*Tl/(rho*pi*xvert.rp^2)); % downstream velocity at the left propeller
%     Vxr = sqrt(Vinr^2 + 2*Tr/(rho*pi*xvert.rp^2)); % downstream velocity at the right propeller
%     
%     Vdxl = (Vinl+Vxl)/2; % velocity at the left propeller disc
%     Vdxr = (Vinr+Vxr)/2; % velocity at the right propeller disc
% end

function ui = slipstream(Vb, Wb, Vwb, rho, xvert, rG, Tr, Tl)
    
    [ui,~,~,~,~,~] = airspeed(Vb, Wb, Vwb, rG, xvert);

    if (0.1008<rG(2)) && (rG(2)<0.1892) % slipstrem right side
        
        r_thr = [0.177, 0.145, 0]' - xvert.CM; % left thruster position [m] --> table 2-3 and H.1

        Vin = Vb + cross(Wb,r_thr) - Vwb; % inflow velocity at the right propeller

        Vx = sqrt(Vin(1)^2 + 2*Tr/(rho*pi*xvert.rp^2));

        if rG(1)<0.148 % downstream velocity at the right propeller
             
            ui=Vx; 
        
        elseif (0.148<rG(1) && rG(1)<0.175) % velocity at the right propeller disc

            ui = (Vin(1)+Vx)/2;
        else
            [ui,~,~,~,~,~] = airspeed(Vb, Wb, Vwb, rG, xvert);
        end
        
    elseif (-0.1892<rG(2)) && (rG(2)<-0.1008) % slipstrem left side
        
        r_thl = [0.177, -0.145, 0]' - xvert.CM; % right thruster position [m] --> table 2-3 and H.1

        Vin = Vb + cross(Wb,r_thl) - Vwb; % inflow velocity at the left propeller
        
        Vx = sqrt(Vin(1)^2 + 2*Tl/(rho*pi*xvert.rp^2)); 

        if rG(1)<0.148 % downstream velocity at the left propeller

            ui=Vx;
        
        elseif (0.148<rG(1) && rG(1)<0.175) % velocity at the left propeller disc

            ui = (Vin(1)+Vx)/2;
        else
            [ui,~,~,~,~,~] = airspeed(Vb, Wb, Vwb, rG, xvert);
        end
    end
end

%% AUXILIAR FUNCTIONS
function varargout = matsplit(A,dim)
    
    if nargin==1
        varargout = num2cell(A);
    else
        varargout = num2cell(A,dim);
    end
end
