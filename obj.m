% Objective Function
% 
% Inputs:
% Initial Guess: xk
% 
% Outputs:
% Function Evaluation: f
% Gradient: delf
% 
% Ian Pylkkanen
% December 16, 2020

function [f, delf] = obj(xk)
% Step-size for Complex Step
h = 1e-60;

% Solve function at initial guess xk
f = subobj(xk);

% Define gradient matrix
delf = zeros(size(xk,1),1);

% Implement Complex Step to find gradient
for i=1:size(xk)
    xc = xk;
    xc(i) = xk(i) + complex(0,h);
    delf(i) = imag(subobj(xc))/h;
end

    % Subfunction to solve for f
    function [f] = subobj(dv)
        % Define variables
        x1 = dv(1);
        x2 = dv(2); % Uncomment of 2D problem
        x3 = dv(3); % Uncomment for 3D problem
        x4 = dv(4); % Uncomment of 4D problem
        
        % Define function f and solve
        %f = x1.^2;
        %f = cos(x1);
        %f = (x1+2*x2-7)^2 + (2*x1+x2-5)^2;
        f = 100*(x1^2-x2)^2 + (x1-1)^2 + (x3-1)^2 + 90*(x3^2-x4)^2 + 10.1*((x2-1)^2 + (x4-1)^2) + 19.8*(x2-1)*(x4-1);
    end
end