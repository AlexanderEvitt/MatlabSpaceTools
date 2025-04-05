function U = ICRS_to_ITRS(JD,TT_UTC,UT1_UTC,xp,yp)
    % Returns matrix for converting celestial to terrestrial
    % Inputs:
    % JD: Julian date at timestep in UT1
    % TT_UTC: Seconds difference, TT - UTC
    % UT1_UTC: Seconds difference, UT1 - UTC
    % xp: polar motion, arcseconds
    % yp: polar motion, arcseconds

    % Time
    TT = JD - ((-UT1_UTC + TT_UTC)/(24*60*60)); % convert to TT first
    T = (TT - 2451545.0)/36525.0;
    
    % Precession term
    % Angles all in arcseconds
    zeta = (2306.2181)*T + (0.30188)*T^2 + (0.017998)*T^2;
    nu = (2004.3109)*T - (0.42665)*T^2 - (0.041833)*T^3;
    z = zeta + (0.79280)*T^2 + (0.000205)*T^3;
    
    % Convert arcseconds -> radians and calculate by eq. 5.46
    P = Rz(deg2rad(z/3600))*Ry(deg2rad(-nu/3600))*Rz(deg2rad(zeta/3600));
    
    % Nutation term
    epsilon = 23.4392911 - (46.8150/3600)*T - (0.00059/3600)*T^2 + (0.001813/3600)*T^3;
    [dpsi,deps] = IAU1980Nutation(T);
    N = Rx(deg2rad(epsilon+deps))*Rz(deg2rad(dpsi))*Rx(deg2rad(-epsilon));
    
    % Rotation term
    T0 = (floor(JD)+0.5 - 2451545)/36525; % in UT1
    T = (JD - (UT1_UTC/(24*60*60)) - 2451545)/36525; % in UT1
    GMST = (24110.54841) + (8640184.812866)*T0 + (0.093104)*T^2 - (0.0000062)*T^3 + 1.002737909350795*(86400*(JD - floor(JD)));
    GAST = (2*pi*GMST/86400) + deg2rad(dpsi)*cosd(epsilon);
    Theta = Rz(-GAST);
    
    % Polar motion term
    xp = deg2rad(xp/3600); % convert to radians
    yp = deg2rad(yp/3600);
    Pi = Ry(xp)*Rx(yp);
    
    % Finally:
    U = Pi*Theta*N*P;
end

%% Helper functions

function mat = Rx(angle)
    % Returns rotation matrix of angle (radians) about x-axis
    mat = [1 0 0; 0 cos(angle) -sin(angle); 0 sin(angle) cos(angle)];
end

function mat = Ry(angle)
    % Returns rotation matrix of angle (radians) about y-axis
    mat = [cos(angle) 0 sin(angle); 0 1 0; -sin(angle) 0 cos(angle)];
end

function mat = Rz(angle)
    % Returns rotation matrix of angle (radians) about z-axis
    mat = [cos(angle) -sin(angle) 0; sin(angle) cos(angle) 0; 0 0 1];
end

function deg = to_deg(deg,min,sec)
    deg = deg + (min/60) + (sec/3600);
end