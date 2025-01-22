% Alexander Evitt - 30 Aug 2024

syms gamma beta alpha
Ax = [1 0 0; 0 cos(gamma) -sin(gamma); 0 sin(gamma) cos(gamma)];
Ay = [cos(beta) 0 sin(beta); 0 1 0; -sin(beta) 0 cos(beta)];
Az = [cos(alpha) -sin(alpha) 0; sin(alpha) cos(alpha) 0; 0 0 1];