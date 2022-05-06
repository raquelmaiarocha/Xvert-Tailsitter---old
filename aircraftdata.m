function xvert = aircraftdata
    
    % Body
    xvert.m = 0.2; % mass of the vehicle [kg]
    % inertia matrix about the body frame axes centered at the center of mass (CM) (eq A.2)
    xvert.I = [  2.986344   -0.022259   -0.012405   ;
                -0.022259    0.637364    0.001210   ;
                -0.012405    0.001210    3.487448   
              ] *10^-3;
    xvert.I_1 = inv([  2.986344   -0.022259   -0.012405   ;
                -0.022259    0.637364    0.001210   ;
                -0.012405    0.001210    3.487448   
              ] *10^-3);
    
    xvert.CM = [0.12834; 0.00110; -0.00344]; % center of mass (CM) in the geometric frame (eq A.1)
    xvert.AC = [0.1264; 0; 0]; % aerodynamic center (AC) in the geometric frame (eq 2.15)

    % Wing
    xvert.cw=0.17;  % mean cord [m]
    xvert.bw=0.5;  % span [m]
    xvert.Sw=0.08;  % area [m^2]
    xvert.delta=deg2rad(19.8); % sweep angle [rad] (eq 2.16)
    xvert.AR=3.13; % aspect ratio [-] (eq 2.16)
    xvert.k0=0.87; % oswald s efficiency factor [-] (Table H.1)
    
    % Additional coefficients
    xvert.CLq=3.1851;
    xvert.Cmq=-2.4487;
    xvert.Cmdaa=-0.4897;
    
    % Flaps
    xvert.dfmax = deg2rad(39); % max flap deflection
    xvert.xflap=0.768;  % [m]
    xvert.yflap=0.750;  % [m]
    xvert.cflap=0.062;  % [m]
    xvert.bflap=0.190;  % [m]
    xvert.Sflap=xvert.bflap*xvert.cflap;  % [m^2]
    
    % Propellers --> eq 2.47 Tese 0
    xvert.rp=0.0625;  % propeller radius [m] --> table H.1
    xvert.ds = sqrt(2)*0.0625; % diameter of the slipstream when the airplane is static (i.e. no inflow velocity to the thrusters)
    xvert.CT0=0.1342;
    xvert.CT1=-0.1196;
    xvert.CT2=-0.1281;
    xvert.CP0=0.0522;
    xvert.CP1=-0.0146;
    xvert.CP2=-0.0602;
    
    % Motor --> eq 2.43 Tese 0
    xvert.wp0 = -4.27;
    xvert.wp1 = 356.34;
    xvert.wp2 = -84.75;
    
    % Electric motor
    xvert.Vbatt=7.4;
    % xvert.Vbatt=14.8;
    % xvert.Lm=3.6*(10^(-4));
    % xvert.Rm=0.032;
    % xvert.Ke=9.68*(10^(-2));
    % xvert.Jm=1.4*(10^(-5));
    % xvert.Kt=2.47*(10^(-2));
    % xvert.Bm=2.9*(10^(-4));
    
    % Wing Segments (Table 2-1 and 2-2)
    xvert.wing=[% Si(1)      bi(2)   ci(3)    cfi(4)  riG(5,6,7)
    
                % Horizontal Segments
                  2.408      22.5    110.0    0       [82.5,   -239.3,  0];
                  4.785      38.3    125.0    63      [93.75,  -208,    0];
                  13.051     88.4    148.2    63      [111.15, -143.4,  0];
                  11.100     63.3    175.9    63      [131.9,  -68.4,   0];
                  17.042     75.0    231.4    0       [173,     0,      0];
                  11.100     63.3    175.9    63      [131.9,   68.4,   0];
                  13.051     88.4    148.2    63      [111.15,  143.4,  0];
                  4.785      38.3    125.0    63      [93.75,   208,    0];
                  2.408      22.5    110.0    0       [82.5,    239.3,  0];
    
                % Vertical Segments
                  3.120      26      120     0       [85, -252, 0];
                  3.120      26      120     0       [85, 252, 0];
    
                ]*10^-3;
end