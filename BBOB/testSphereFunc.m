function fposition = testSphereFunc( x )
%TESTSPHEREFUNC Summary of this function goes here
%   Detailed explanation goes here

% % Custom: Sphere function
%fposition = sum(x.^2);

% % Custom: Rosenbrock function
fposition = 100 * sum((x(1:end-1).^2 - x(2:end)).^2) + sum((x(1:end-1)-1).^2);


%     p=0;q=0;
%     for k=1:5
%         p=p+k*cos((k+1)*x(1)+k);
%         q=q+k*cos((k+1)*x(2)+k);
%     end
% 
% fposition=p*q+(x(1)+1.42513)^2+(x(2)+.80032)^2;

end

