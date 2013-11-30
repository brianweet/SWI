% Function to be optimized:
% Author: K. Passino, Version: 5/23/00
function J=nutrientsfunc1(theta,flag)

if flag==1  % Test to see if main program indicated nutrients
	J=0;
	return
end


% To test concentric pattern of groups of bacteria

%J=-1000+1000*exp(-0.07*((theta(1,1)-15)^2+(theta(2,1)-15)^2));

% Sphere 
%J = sum(theta.^2);

% ackley
x = theta';
[ps,n] = size(x);
  J = 20 - 20 * exp(-0.2 * sqrt(sum(x .^ 2, 2) / n)) ...
        - exp(sum(cos(2 * pi * x), 2) / n) + exp(1);
