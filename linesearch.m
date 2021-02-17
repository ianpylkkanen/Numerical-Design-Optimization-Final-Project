% Line-Search Algorithm
% 
% Adopted from: Nocedal, Wright "Numercial Optimization", 2006
% Algorithm 3.5
% 
% Inputs:
% Initial Guess: xk
% Initial Step Length: alpha1
% Search Direction: pk
%
% Outputs:
% Step Length Satisfying Wolfe Conditions: alphastar
% 
% Ian Pylkkanen
% December 16, 2020

function [alphastar] = linesearch(xk, alpha1, pk)
% Define Initial Parameters
alpha0 = 0;
alphamax = 5*alpha1;
xk0 = xk + alpha0*pk;
c1 = 10^-4;
c2 = .9;
i = 1;
% Call Objective Function to be minimized
[phi0, phip0] = obj(xk0);
% Define Objective at Initial Guess
phip0 = phip0'*pk;
[phi_0, phip_0] = obj(xk);
phip_0 = phip_0'*pk;

while 1
    xk1 = xk + alpha1*pk;
    [phi1, phip1] = obj(xk1);
    phip1 = phip1'*pk;
    % If alpha fails sufficient decrease, send to Zoom
    if phi1 > phi_0 + c1*alpha1*phip_0 || (phi1 >= phi0 && i>1) 
        alphastar = zoom(alpha0, alpha1);
        break
    end
    % If alpha passes sufficient decrease, set alphastar=alpha1
    if abs(phip1) <= -c2*phip_0
        alphastar = alpha1;
        break
    end
    % If alpha fails curvature condition, send to Zoom
    if phip1 >= 0
        alphastar = zoom(alpha1, alpha0);
        break
    end
    alpha0 = alpha1;
    alpha1 = (alpha1+alphamax)/2;
    i = i+1;
end

% Zoom Function

% Adopted from: Nocedal, Wright "Numercial Optimization", 2006
% Algorithm 3.6

% Inputs:
% Lower and Upper Bounds: [alphalo, alphahi]

% Outputs:
% Step Length Satisfying Wolfe Conditions: alphastar

% Zoom Function looks to find alphastar to satisfy Wolfe Conditions
function [alphastar] = zoom(alphalo, alphahi)
while 1
    alphaj = (alphahi+alphalo)/2; % Define Trial Step Length
    xj = xk + pk*alphaj;
    [phij, phipj] = obj(xj);
    phipj = phipj'*pk;
    [philo, phiplo] = obj(xk + pk*alphalo);
    phiplo = phiplo'*pk;
    if phij > phi0 + c1*alphaj*phip0 || phij >= philo
        alphahi = alphaj;        
    else
        if abs(phipj) <= -c2*phip0
            alphastar = alphaj;
            return
        end
        if phipj*(alphahi-alphalo) >= 0
            alphahi = alphalo;
        end
        alphalo = alphaj;
    end
end
end
end