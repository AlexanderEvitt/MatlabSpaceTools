function [v1,v2] = piter(r1,r2,df,t,mu)
    % Total constants
    alpha = 0.01;
    tol = 1e-8;

    % Get klm constants
    k = norm(r1)*norm(r2)*(1 - cos(df));
    l = norm(r1) + norm(r2);
    m = norm(r1)*norm(r2)*(1 + cos(df));

    % Initial guess for p
    p_i = k/(l + sqrt(2*m));
    p_ii = k/(l - sqrt(2*m));
    p = (p_i+p_ii)/2;
    diff = 1;

    % iterate with p-iteration
    while diff >= tol
        % Get F and G from p
        F = 1 - (norm(r2)/p)*(1- cos(df));
        G = (norm(r1)*norm(r2)*sin(df))/sqrt(mu*p);
        Fdot = sqrt(mu/p)*tan(df/2)*(((1 - cos(df))/p) - (1/norm(r1)) - (1/norm(r2)));
        Gdot = 1 - (norm(r1)/p)*(1- cos(df));

        % Split by a to determine if ellipse or hyperbola
        a = (m*k*p)/((2*m - l^2)*p^2 + 2*k*l*p - k^2);

        if a > 0
            % Calculate dE
            dE = wrapTo2Pi(atan2( (-norm(r1)*norm(r2)*Fdot)/(sqrt(mu*a)), 1 - (norm(r1)/a)*(1 - F)));

            % Iterating function and its derivative
            iter = G + sqrt((a^3)/mu) * (dE - sin(dE));
            diter = (-G/(2*p)) - (3/2)*a*(iter-G)*((k^2 + (2*m - l^2)*p^2)/(m*k*p^2)) + sqrt((a^3)/mu)*(2*k*sin(dE))/(p*(k - l*p));
            
            % Update p
            p = p + (alpha*(t - iter)/diter);


        elseif a < 0
            % Calculate dF (diff in hyperbolic eccentric anomaly)
            dF = acosh(1 - (norm(r1)/a)*(1-F));

            % Iterating function and its derivative
            iter = G + sqrt(((-a)^3)/mu) * (dF - sinh(dF));
            diter = (-G/(2*p)) - (3/2)*a*(itert-G)*((k^2 + (2*m - l^2)*p^2)/(m*k*p^2)) + sqrt(((-a)^3)/mu)*(2*k*sinh(dF))/(p*(k - l*p));
            
            % Update p
            p = p + alpha*(t - iter)/diter;

        else
            p = NaN;

        end

        % Calculate diff
        diff = (t - iter);
    end
    
    % Recalculate F and G from final p
    F = 1 - (norm(r2)/p)*(1- cos(df));
    G = (norm(r1)*norm(r2)*sin(df))/sqrt(mu*p);
    Gdot = 1 - (norm(r1)/p)*(1- cos(df));

    % Convert F and G into final answers
    v1 = (r2 - F*r1)/G;
    v2 = (Gdot*r2 - r1)/G;
end