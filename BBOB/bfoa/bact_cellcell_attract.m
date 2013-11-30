%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bacterial cell-to-cell attraction function
% Author: K. Passino
% Version: 5/16/00
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x=0:31/200:30;   % For our function the range of values we are considering
y=x;

% Compute the function that we are trying to find the minimum of.

theta(:,1)=[15; 15]; % Just put at middle
theta(:,2)=[15; 20]; % Put next to it


for jj=1:length(x)
	for ii=1:length(y)
		z(ii,jj)=bact_cellcell_attract_func([x(jj);y(ii)],theta,2,0);
	end
end

figure(1)
clf
surf(x,y,z);
colormap(jet)
% Use next line for generating plots to put in black and white documents.
colormap(white);
view(-53,78)
xlabel('x=\theta_1');
ylabel('y=\theta_2');
zlabel('z=J');
title('Cohesion-repulsion function');
zoom(2)
%rotate3d on


figure(2)
clf
contour(x,y,z,25)
colormap(jet)
xlabel('x=\theta_1');
ylabel('y=\theta_2');
title('Cohesion-repulsion function for two cells');