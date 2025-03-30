function acc = nonspherical_accel(CS, r, order, min_order, Rp, mu)
    % Nonspherical gravity calculator
    % Alexander Evitt, 26 March 2025
    
    % CS ([order]x[order]): table of spherical harmonic coefficients
    % r (3x1): point where acceleration is evaluated, relative to planet
    % order (int): maximum n
    % min_order (int): minimum n
    % Rp (float): planet radius
    % mu (float): planet gravitational parameter

    % Set up constants
    x = r(1); y = r(2); z = r(3);
    r = norm(r);
    Rpovr2 = Rp/(r^2); % precalculate these for speed
    muovRp2 = mu/(Rp^2);

    % Initialize V and W matrices (up to n+1)
    V = zeros(order+2, order+2);
    W = zeros(order+2, order+2);
    
    % Set the seed values
    V(1,1) = Rp/r; 
    W(1,1) = 0; 

    % Compute values for V and W using recursion up to n+1
    % First, compute the zonal terms Vn0 (m = 0)
    for n = 1:order+1
        V(n+1,1) = down(0,n,z,V,Rp,r);
        W(n+1,1) = down(0,n,z,W,Rp,r);
    end
    
    % Compute tesseral terms Vnm and Wnm
    for m = 1:order % starts at 1 because we know the 0,0 element
        % First, compute diagonal elements using Eqn (3.29)
        V(m+1, m+1) = (2*m - 1) * ((x * Rpovr2) * V(m, m) - (y * Rpovr2) * W(m, m));
        W(m+1, m+1) = (2*m - 1) * ((x * Rpovr2) * W(m, m) + (y * Rpovr2) * V(m, m));
        
        for n = m+1:order+1
            % Compute off-diagonal elements using Eqn (3.30)
            V(n+1, m+1) = down(m,n,z,V,Rp,r);
            W(n+1, m+1) = down(m,n,z,W,Rp,r);
        end
    end

    % Summation for acceleration
    xdd = 0;
    ydd = 0;
    zdd = 0;
    
    for n = min_order:order
        % Start summations with m = 0 terms
        xdd = xdd + (mu/(Rp^2))*(-CS(n+1,1)*V(n+2,2));
        ydd = ydd + (mu/(Rp^2))*(-CS(n+1,1)*W(n+2,2));
        zdd = zdd + (mu/(Rp^2))*(n-0+1)*(-CS(n+1,1)*V(n+2,0+1) - 0*W(n+2,1));
    
        % Iterate through remaining terms
        for m = 1:n
            fr = factorial_ratio(n-m+2,n-m);
            Cnm = CS(n+1,m+1);
            if m == 0
                Snm = 0;
            else
                Snm = CS(m,n+1);
            end
    
            xdd = xdd + muovRp2*0.5*((-Cnm*V(n+2,m+2) - Snm*W(n+2,m+2)) + fr*(Cnm*V(n+2,m) + Snm*W(n+2,m)));
            ydd = ydd + muovRp2*0.5*((-Cnm*W(n+2,m+2) + Snm*V(n+2,m+2)) + fr*(-Cnm*W(n+2,m) + Snm*V(n+2,m)));
            zdd = zdd + muovRp2*(n-m+1)*(-Cnm*V(n+2,m+1) - Snm*W(n+2,m+1));
        end
    end

    % Assign output
    acc = [xdd;ydd;zdd];
end

function out = down(m,n,z,M,Rp,r)
    if n > 1
        out = ((2*n - 1) / (n - m)) * (z * Rp / r^2) * M(n, m + 1) - ((n + m - 1) / (n - m)) * (Rp^2 / r^2 * M(n-2+1, m+1));
    else
        out = ((2*n - 1) / (n - m)) * (z * Rp / r^2) * M(n, m + 1);
    end
end

function out = factorial_ratio(num,dem)
    % Returns num!/dem! given num>dem
    out = prod(dem+1:num);
end