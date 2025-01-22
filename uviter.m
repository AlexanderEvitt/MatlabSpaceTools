function [v1,v2] = uviter(r1,r2,df,t,mu)
    % Total constants
    alpha = 0.01;
    tol = 1e-8;

    % Initial guess for z
    z = 0;
    diff = 1;

    % iterate with p-iteration
    while diff >= tol
        % Get Stumpff constants
        C = stumpffC(z);
        S = stumpffS(z);

        % Find parameter A by long way (-1) or short way (1)
        way = 1;
        A = way*sqrt(norm(r1)*norm(r2)*(1 + cos(df)));

        % More parameter definitions
        y = norm(r1) + norm(r2) - A*(1 - z*S)/sqrt(C);
        x = sqrt(y/C);

        % Define the iterative function
        iter = S*x^3 + A*sqrt(y);
        iter = iter/sqrt(mu);

        % Define the derivatives
        if abs(z) < 1e-6
            dC = -(1/factorial(4)) + (2*z/factorial(6));
            dS = -(1/factorial(5)) + (2*z/factorial(7));
        else
            dC = (1/(2*z))*(1 - S*z - 2*C);
            dS = (1/(2*z))*(C-3*S);
        end
        diter = (dS - 3*S*dC/(2*C))*x^3 + (A/8)*((3*S*sqrt(y)/C) + (A/x));
        diter = diter/sqrt(mu);

        % Update z
        z = z + alpha*(t - iter)/diter;

        % Calculate diff
        diff = (t - iter);
    end
    
    % Recalculate F and G from final z
    F = 1 - (x^2 / norm(r1))*C;
    G = t - (x^3 / sqrt(mu))*S;
    Fdot = -(sqrt(mu)/(norm(r1)*norm(r2)))*x*(1 - z*S);
    Gdot = 1 - (x^2 / norm(r2))*C;

    % Convert F and G into final answers
    v1 = (r2 - F*r1)/G;
    v2 = (Gdot*r2 - r1)/G;
end

function C = stumpffC(z)
    if z > 0
        C = (1 - cos(sqrt(z)))/z;
    elseif z == 0
        C = 0.5;
    elseif z < 0
        C = (cosh(sqrt(-z)) - 1)/(-z);
    else
        C = NaN;
    end
end

function S = stumpffS(z)
    if z > 0
        S = (sqrt(z) - sin(sqrt(z)))/(sqrt(z)^3);
    elseif z == 0
        S = 1/6;
    elseif z < 0
        S = (sinh(sqrt(-z)) - sqrt(-z))/(sqrt(-z)^3);
    else
        S = NaN;
    end
end