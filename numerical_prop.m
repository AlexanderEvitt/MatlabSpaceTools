% Alexander Evitt, 26 Oct 2024
clear; close all; clc;

mu = 3.986e5;
n = 10000;
line2 = '2 33153 0.0103 212.4984 0002853 338.7926 217.7793 1.00269427 52146';

% Extract initial conditions from TLE
[a,e,i,RAAN,argp,theta] = extract_TLE(line2,mu);
[r0, v0] = orbital_elements_to_rv(a,e,i,RAAN,argp,theta, mu);
ts = linspace(0,86164.0905,n);

% Numerically integrate with ode45
options = odeset('AbsTol',1e-14,'RelTol',1e-14);
[t_out, r_out] = ode45(@orbdyn, ts, [r0;v0], options);

% Plot Earth on ECI plot
figure(1)
subplot(3,1,1)
axis equal;
npanels=20;
erad = 6378.1; % equatorial radius (km)
prad = 6356.8; % polar radius (km)
view(9,26);
hold on;
axis vis3d;
[ xx, yy, zz ] = ellipsoid(0, 0, 0, erad, erad, prad, npanels);
globe = surf(xx, yy, -zz, "FaceColor", "none", "EdgeColor", 0.5*[1 1 1]);
xlabel("x-ECI (km)")
ylabel("y-ECI (km)")
zlabel("z-ECI (km)")

% Plot orbit
plot3(r_out(:,1),r_out(:,2),r_out(:,3))

final_pos = r_out(end,1:3);
final_vel = r_out(end,4:6);
plot3(final_pos(1),final_pos(2),final_pos(3),"rx")

% Plot ECEF plot
subplot(3,1,2)
axis equal;
npanels=20;
erad = 6378.1; % equatorial radius (km)
prad = 6356.8; % polar radius (km)
view(9,26);
hold on;
axis vis3d;
[ xx, yy, zz ] = ellipsoid(0, 0, 0, erad, erad, prad, npanels);
globe = surf(xx, yy, -zz, "FaceColor", "none", "EdgeColor", 0.5*[1 1 1]);

A = @(alpha) [cos(alpha) -sin(alpha) 0; sin(alpha) cos(alpha) 0; 0 0 1];

ECEFs = zeros(3,n);
for i = 1:n
    trans = A(-2*pi*double(i/n))*(r_out(i,1:3).');
    ECEFs(:,i) = trans;
end

plot3(ECEFs(1,:),ECEFs(2,:),ECEFs(3,:))
hold on
plot3(ECEFs(1,end),ECEFs(2,end),ECEFs(3,end),"rx")
xlabel("x-ECEF (km)")
ylabel("y-ECEF (km)")
zlabel("z-ECEF (km)")

% Plot ground track
subplot(3,1,3)
LLAs = ecef2lla(1000.*ECEFs.');
scatter(LLAs(:,2),LLAs(:,1), "b.")
xlabel("Longitude (deg)")
ylabel("Latitude (deg)")
xlim([-180,180])









function dydt = orbdyn(t,y)
    mu = 3.986004418e5; % km^3 / s^2
    r = y(1:3);
    v = y(4:6);
    
    a = (-mu*r)/(norm(r)^3);
    a = a + -(3/2)*0.0010827*((mu/norm(r)^2)*(6371/norm(r))^2)*[(1 - 5*(r(3)/norm(r))^2)*r(1)/norm(r);(1 - 5*(r(3)/norm(r))^2)*r(2)/norm(r);(3 - 5*(r(3)/norm(r))^2)*r(3)/norm(r)];
    
    dydt(1:3,1) = v;
    dydt(4:6,1) = a;
end