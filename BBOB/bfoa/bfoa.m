function xbest = bfoa(FUN, DIM, ftarget, maxfunevals)
%TODO BwE: use ftarget and maxfunevals

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   E. coli Bacterial Swarm Foraging for Optimization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Kevin Passino
%   Version: 6/1/00
% 
% This program simulates the minimization of a simple function with
% chemotaxis, swarming, reproduction, and elimination/dispersal of a 
% E. coli bacterial population.
%
% To change the nutrientsfunc search on it.  For
% example, change it to nutrientsfunc1 to study another type of swarm behavior.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%clear all            % Initialize memory
close all

p=DIM;                 % Dimension of the search space

S=50;	% The number of bacteria in the population (for convenience, require this to be an
        % an even number)

Nc=100; % Number of chemotactic steps per bacteria lifetime (between reproduction steps), assumed
        % for convenience to be the same for every bacteria in the population
Ns=4;   % Limits the length of a swim when it is on a gradient


Nre=4;	 % The number of reproduction steps (right now the plotting is designed for Nre=4)
Sr=S/2;	 % The number of bacteria reproductions (splits) per generation (this
		 % choice keeps the number of bacteria constant)
		 

Ned=2; % The number of elimination-dispersal events (Nre reproduction steps in between each event)

ped=0.25; % The probabilty that each bacteria will be eliminated/dispersed (assume that 
          % elimination/dipersal events occur at a frequency such that there can be 
		  % several generations of bacteria before an elimination/dispersal event but
		  % for convenience make the elimination/dispersal events occur immediately after
		  % reproduction)

flag=2; % If flag=0 indicates that will have nutrients and cell-cell attraction
        % If flag=1 indicates that will have no (zero) nutrients and only cell-cell attraction
		% If flag=2 indicates that will have nutrients and no cell-cell attraction
		
% Initial population

P(:,:,:,:,:)=0*ones(p,S,Nc,Nre,Ned);  % First, allocate needed memory

% Initialize locations of bacteria all at the center (for studying one swarming case - when nutrientsfunc1 is used)
% for m=1:S
% 	P(:,m,1,1,1)=[15;15];
% end


% % Another initialization possibility: Randomly place on domain:
for m=1:S
	%P(:,m,1,1,1)=(15*((2*round(rand(p,1))-1).*rand(p,1))+[15;15]);
	%BwE Because of bbob: initialize domain uniformly randomly in [-5,5]^DIM
	%P(:,m,1,1,1)= 10 * rand(DIM, popsize) - 5;
    P(:,m,1,1,1)=(5^DIM*((2*round(rand(p,1))-1).*rand(p,1)));
end

% Next, initialize the parameters of the bacteria that govern
% part of the chemotactic behavior 

C=0*ones(S,Nre); % Allocate memory

% Set the basic run length step (one step is taken of this size if it does not go up a gradient
% but if it does go up then it can take as many as Ns such steps).  C(i,k) is the step size for
% the ith bacteria at the kth reproduction step.  For now, the step size is assumed to be constant since 
% we assume that perfect copies of bacteria are made, but later you can add evolutionary effects to modify 
% this and perhaps Ns and Nc.

runlengthunit=0.1;
C(:,1)=runlengthunit*ones(S,1);


% Allocate memory for cost function:

J=0*ones(S,Nc,Nre,Ned);
Jhealth=0*ones(S,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%---------------------------------
% Elimination-dispersal loop: 
%---------------------------------

for ell=1:Ned

%---------------------------------
% Reproduction loop: 
%---------------------------------

for k=1:Nre

%---------------------------------
% Swim/tumble (chemotaxis) loop:
%---------------------------------

for j=1:Nc

	for i=1:S  % For each bacterium
		
		% Compute the nutrient concentration at the current location of each bacterium
		
% 		J(i,j,k,ell)=nutrientsfunc(P(:,i,j,k,ell),flag);
%		J(i,j,k,ell)=nutrientsfunc1(P(:,i,j,k,ell),flag);
		%BwE Because of bbob: use the FUN that is given as parameter
		J(i,j,k,ell)=feval(FUN, P(:,i,j,k,ell));

		% Next, add on cell-cell attraction effect:
		
		J(i,j,k,ell)=J(i,j,k,ell)+bact_cellcell_attract_func(P(:,i,j,k,ell),P(:,:,j,k,ell),S,flag);
		
		%-----------
		% Tumble:
		%-----------
		
		Jlast=J(i,j,k,ell); % Initialize the nutrient concentration to be the one at the tumble
							% (to be used below when test if going up gradient so a run should take place)

		% First, generate a random direction

		Delta(:,i)=(2*round(rand(p,1))-1).*rand(p,1);

		% Next, move all the bacteria by a small amount in the direction that the tumble resulted in
		% (this implements the "searching" behavior in a homogeneous medium)
		
		P(:,i,j+1,k,ell)=P(:,i,j,k,ell)+C(i,k)*Delta(:,i)/sqrt(Delta(:,i)'*Delta(:,i));
										% This adds a unit vector in the random direction, scaled
										% by the step size C(i,k)
		
		%---------------------------------------------------------------------
		% Swim (for bacteria that seem to be headed in the right direction):
		%---------------------------------------------------------------------

%		J(i,j+1,k,ell)=nutrientsfunc(P(:,i,j+1,k,ell),flag); % Nutrient concentration for each bacterium after
%		J(i,j+1,k,ell)=nutrientsfunc1(P(:,i,j+1,k,ell),flag); % Nutrient concentration for each bacterium after
															% a small step (used by the bacterium to
															% decide if it should keep swimming)
		%BwE Because of bbob: use the FUN that is given as parameter
		J(i,j+1,k,ell)=feval(FUN, P(:,i,j+1,k,ell));
		
		% Next, add on cell-cell attraction effect:
		
		J(i,j+1,k,ell)=J(i,j+1,k,ell)+bact_cellcell_attract_func(P(:,i,j+1,k,ell),P(:,:,j+1,k,ell),S,flag);
															
		m=0; % Initialize counter for swim length 
		
		while m<Ns  % While climbing a gradient but have not swam too long...
			
			m=m+1;
			
			if J(i,j+1,k,ell)<Jlast  % Test if moving up a nutrient gradient.  If it is then move further in
				                     % same direction
				Jlast=J(i,j+1,k,ell); % First, save the nutrient concentration at current location
									  % to later use to see if moves up gradient at next step
									  
				% Next, extend the run in the same direction since it climbed at the last step
				
				P(:,i,j+1,k,ell)=P(:,i,j+1,k,ell)+C(i,k)*Delta(:,i)/sqrt(Delta(:,i)'*Delta(:,i));
				
%				J(i,j+1,k,ell)=nutrientsfunc(P(:,i,j+1,k,ell),flag); % Find concentration at where
%				J(i,j+1,k,ell)=nutrientsfunc1(P(:,i,j+1,k,ell),flag); % Find concentration at where
																	% it swam to and give it new cost value
				%BwE Because of bbob: use the FUN that is given as parameter
				J(i,j+1,k,ell)=feval(FUN, P(:,i,j+1,k,ell));													
																	
				% Next, add on cell-cell attraction effect:
		
				J(i,j+1,k,ell)=J(i,j+1,k,ell)+bact_cellcell_attract_func(P(:,i,j+1,k,ell),P(:,:,j+1,k,ell),S,flag);
																	
			else  % It did not move up the gradient so stop the run for this bacterium
				m=Ns;
			end
		
		end	% Test if should end run for bacterium
	end  % Go to next bacterium

	
%---------------------------------
end  % j=1:Nc
%---------------------------------

	% Reproduction
	
	Jhealth=sum(J(:,:,k,ell),2);  % Set the health of each of the S bacteria.
	                                     % There are many ways to define this; here, we sum
										 % the nutrient concentrations over the lifetime of the
										 % bacterium.

	% Sort cost and population to determine who can reproduce (ones that were in best nutrient
	% concentrations over their life-time reproduce)
	
	[Jhealth,sortind]=sort(Jhealth); % Sorts the nutrient concentration in order 
									% of ascending cost in the first dimension (bacteria number)
									% sortind are the indices in the new order
		
	P(:,:,1,k+1,ell)=P(:,sortind,Nc+1,k,ell); % Sorts the population in order of ascending Jhealth (the
											% ones that got the most nutrients were the ones with 
											% the lowest Jhealth values)
	
	C(:,k+1)=C(sortind,k); % And keeps the chemotaxis parameters with each bacterium at the next generation
		
	% Split the bacteria (reproduction)
	
	for i=1:Sr
		P(:,i+Sr,1,k+1,ell)=P(:,i,1,k+1,ell); % The least fit do not reproduce, the most 
		 									% fit ones split into two identical copies 
		C(i+Sr,k+1)=C(i,k+1); 	% and they get the same parameters as for their mother
	end

	% Evolution can be added here (can add random modifications to C(i,k), Ns, Nc, etc)
	

%---------------------------------
end  % k=1:Nre
%---------------------------------
	
	% Eliminate and disperse (on domain for our function) - keep same parameters C
	
	for m=1:S
		if ped>rand  % Generate random number and if ped bigger than it then eliminate/disperse
			P(:,m,1,1,ell+1)=(15*((2*round(rand(p,1))-1).*rand(p,1))+[15;15]);
		else
			P(:,m,1,1,ell+1)=P(:,m,1,Nre+1,ell);  % Bacteria that are not dispersed
		end
	end

%---------------------------------
end  % ell=1:Ned
%---------------------------------


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the function we are seeking the minimum of:
%BwE range is [-5,5]^DIM

x=-5^DIM:5^DIM/100:5^DIM;   % For our function the range of values we are considering
y=x;

% Compute the function that we are trying to find the minimum of.

for jj=1:length(x)
	for ii=1:length(y)
%		z(ii,jj)=nutrientsfunc([x(jj);y(ii)],flag);
%		z(ii,jj)=nutrientsfunc1([x(jj);y(ii)],flag);
		%BwE Because of bbob: use the FUN that is given as parameter
		%This is only for plotting so this should be removed.
		z(ii,jj)=feval(FUN, [x(jj);y(ii)]);
		
	end
end

% First, show the actual function to be maximized

figure(1)
clf
surf(x,y,z);
colormap(jet)
% Use next line for generating plots to put in black and white documents.
colormap(white);
xlabel('x=\theta_1');
ylabel('y=\theta_2');
zlabel('z=J');
title('Nutrient concentration (valleys=food, peaks=noxious)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Next, provide some plots of the results of the simulation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


t=0:Nc;  % For use in plotting (makes t=0 correspond to the 1 index and plots to Nc+1)

% As is Figure 2 (4) shows parameter trajectories for 4 generations, then there is an elimination/dispersal event
% Figure 3 (5) shows parameter trajectories for the following 4 generations (after the elim/disp event)

for kk=1:Ned
figure(kk+1) 
clf
for mm=1:Nre
subplot(2,2,mm)
for nn=1:S % Plot all bacteria trajectories for generation mm
plot(t,squeeze(P(1,nn,:,mm,kk)),t,squeeze(P(2,nn,:,mm,kk)))
%plot(t,squeeze(P(1,nn,:,mm,kk)),'-',t,squeeze(P(2,nn,:,mm,kk)),'-')
axis([min(t) max(t) -5^DIM 5^DIM])
hold on
end
T=num2str(mm);
T=strcat('Bacteria trajectories, Generation=',T);
title(T)
xlabel('Iteration, j')
ylabel('\theta_1, \theta_2')
hold off
end
end

for kk=1:Ned
figure(Ned+kk+1) 
clf
for mm=1:Nre
subplot(2,2,mm)
contour(x,y,z,25)
colormap(jet)
for nn=1:S  % Plot all bacteria trajectories for generation mm
hold on
plot(squeeze(P(1,nn,:,mm,kk)),squeeze(P(2,nn,:,mm,kk)))
%plot(squeeze(P(1,nn,:,mm,kk)),squeeze(P(2,nn,:,mm,kk)),'-')
axis([-5^DIM 5^DIM -5^DIM 5^DIM])
end
T=num2str(mm);
T=strcat('Bacteria trajectories, Generation=',T);
title(T)
% Use next line for generating plots to put in black and white documents.
%colormap(gray);
xlabel('\theta_1');
ylabel('\theta_2');
hold off
end
end

%%%%%%%%%%%%
%pause % Can leave this in if want to avoid movie (then hit control-C)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Next, show a movie of the chemotactic steps:

figure(2*Ned+2) 
clf
		contour(x,y,z,25)
		colormap(jet)
		axis([-5^DIM,5^DIM,-5^DIM,5^DIM]);
		xlabel('\theta_1');
		ylabel('\theta_2');
		title('Bacteria movements');
hold on

M = moviein(Nc);
	for j=1:Nc;
% Can change generation step and elimination-dispersal step on next line.
% Currently for 1,1 - the first generation in the first elimination dispersal step
        for i=1:S
		v=plot(squeeze(P(1,i,j:j+1,1,1)),squeeze(P(2,i,j:j+1,1,1)),'-');
		set(v,'MarkerSize',3);
		end
        M(:,j)=getframe;
    end;
	
%movie(M,0)
%save bacteria_swarm


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%