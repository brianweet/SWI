
function fposition=Live_fn(x)

%     p=0;q=0;
%     for k=1:5
%         p=p+k*cos((k+1)*x(1)+k);
%         q=q+k*cos((k+1)*x(2)+k);
%     end
% 
% fposition=p*q+(x(1)+1.42513)^2+(x(2)+.80032)^2;

fposition = sum(x.^2);

% fposition = 100 * sum((x(1:end-1).^2 - x(2:end)).^2) + sum((x(1:end-1)-1).^2);

% [ps,n] = size(x);
%   fposition = 20 - 20 * exp(-0.2 * sqrt(sum(x .^ 2, 2) / n)) ...
%         - exp(sum(cos(2 * pi * x), 2) / n) + exp(1);
