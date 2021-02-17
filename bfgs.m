% BFGS Algorithm
% 
% Adopted from: Nocedal, Wright "Numercial Optimization", 2006
% Algorithm 6.1
% 
% Inputs:
% nital guess: xk
% Tolerance: eps
% Objective Function: obj
%
% Outputs:
% Local Minimizer: f(xk), xk
% 
% Ian Pylkkanen
% December 16, 2020

function bfgs(xk,eps)

% Define Initial Guess of Inverse Hessian
Hk = eye(length(xk));

% Call Objective Function ant Initial Guess
[fk, delfk] = obj(xk);

while norm(delfk) > eps
    pk = -Hk*delfk; % Define Search Direction
    pk = pk/norm(pk);
    alpha1 = 1; % Initial Guess for Step-Length
    alphak = linesearch(xk, alpha1, pk); % Define Step Length from Line-Search
    xk1 = xk + alphak*pk; % Define new x-value
    sk = xk1 - xk; % Displacement
    [fk1, delfk1] = obj(xk1); % Solve objective at new xk1
    yk = delfk1 - delfk; % Change in Gradient
    rhok = 1./(yk'*sk); % 
    I = eye(length(xk));
    Hk = (I - rhok*sk*yk')*Hk*(I - rhok*yk*sk') + rhok*sk*sk'; % Solve for Inverse Hessian
    xk = xk1;
    [fk, delfk] = obj(xk);
end
% Define Final Values
delfk = norm(delfk); % Gradient
xk; % x-values
fk = norm(fk); %Function Evaluation
end