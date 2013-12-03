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
max_generation = 1000;

% Set the basic run length step (one step is taken of this size if it does not go up a gradient
C=0*ones(S,max_generation); % Allocate memory

runlengthunit=0.1;
C(:,1)=runlengthunit*ones(S,1);

e=zeros(1,max_generation);
e(1,1)=100;

n = 10;
alpha = 10;
beta = 10;


%Nc=50; % Number of chemotactic steps per bacteria lifetime (between reproduction steps), assumed
        % for convenience to be the same for every bacteria in the population
Ns=4;   % Limits the length of a swim when it is on a gradient


%Nre=4;	 % The number of reproduction steps (right now the plotting is designed for Nre=4)
Sr=S/2;	 % The number of bacteria reproductions (splits) per generation (this
		 % choice keeps the number of bacteria constant)
		 

%Ned=10; % The number of elimination-dispersal events (Nre reproduction steps in between each event)

ped=0.25; % The probabilty that each bacteria will be eliminated/dispersed (assume that 
          % elimination/dipersal events occur at a frequency such that there can be 
		  % several generations of bacteria before an elimination/dispersal event but
		  % for convenience make the elimination/dispersal events occur immediately after
		  % reproduction)

% Initial population

%P(:,:,:,:,:)=0*ones(p,S,Nc,Nre,Ned);  % First, allocate needed memory
P(:,:,:)=zeros(p,S,max_generation);  % First, allocate needed memory

% BBOB search space bounds (lb = -ub)
ub = 4;

% % Another initialization possibility: Randomly place on domain:
for m=1:S
	%BwE Because of bbob: initialize domain uniformly randomly in [-4,4]
    P(:,m,1)=(ub*((2*round(rand(p,1))-1).*rand(p,1)));
end




% Allocate memory for cost function:

J=0*ones(S,max_generation);
%JbestX=cell(S,max_generation);
JbestX=zeros(S,max_generation,p);
JbestJ=inf(S,max_generation);
Jhealth=0*ones(S,1);

% BBOB variables
fbest = inf;
xbest = zeros(2,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%---------------------------------
% Elimination-dispersal loop: 
%---------------------------------

%for ell=1:Ned

%---------------------------------
% Reproduction loop: 
%---------------------------------

%for k=1:Nre

%---------------------------------
% Swim/tumble (chemotaxis) loop:
%---------------------------------

%for j=1:Nc
for t = 1:max_generation
	for i=1:S  % For each bacterium
        %BwE: keep track of best position for this bacterium for this generation
		sBestJ=inf;
        sBestX=[0;0];
        sasas = P(:,i,t);
        
		% Compute the nutrient concentration at the current location of this bacterium
		[J(i,t), fbest, xbest, sBestJ, sBestX] = ...
            evaluate_function(FUN, P(:,i,t), fbest, xbest, sBestJ, sBestX);
        %BwE Because of bbob: use the FUN that is given as parameter and quit if fbest < ftarget
        if fbest < ftarget
            return
        end
		
		%-----------
		% Tumble:
		%-----------
		
		Jlast=J(i,t); % Initialize the nutrient concentration to be the one at the tumble
							% (to be used below when test if going up gradient so a run should take place)

		% First, generate a random direction

		Delta(:,i)=(2*round(rand(p,1))-1).*rand(p,1);

		% Next, move all the bacteria by a small amount in the direction that the tumble resulted in
		% (this implements the "searching" behavior in a homogeneous medium)
		
		P(:,i,t+1)=P(:,i,t)+C(i,t)*Delta(:,i)/sqrt(Delta(:,i)'*Delta(:,i));
										% This adds a unit vector in the random direction, scaled
										% by the step size C(i,k)
		
		%---------------------------------------------------------------------
		% Swim (for bacteria that seem to be headed in the right direction):
		%---------------------------------------------------------------------
        % Nutrient concentration for each bacterium after
        % a small step (used by the bacterium to decide if it should keep swimming)
        
		%BwE Because of bbob: use the FUN that is given as parameter and quit if fbest < ftarget
        [J(i,t+1), fbest, xbest, sBestJ, sBestX] = ...
            evaluate_function(FUN, P(:,i,t+1), fbest, xbest, sBestJ, sBestX);
        if fbest < ftarget
            return
        end
															
		m=0; % Initialize counter for swim length 
		
		while m<Ns  % While climbing a gradient but have not swam too long...
			
			m=m+1;
			
			if J(i,t+1)<Jlast  % Test if moving up a nutrient gradient.  If it is then move further in
				                     % same direction
				Jlast=J(i,t+1); % First, save the nutrient concentration at current location
									  % to later use to see if moves up gradient at next step
									  
				% Next, extend the run in the same direction since it climbed at the last step
				
				P(:,i,t+1)=P(:,i,t+1)+C(i,t)*Delta(:,i)/sqrt(Delta(:,i)'*Delta(:,i));
				
                % Find concentration at where it swam to and give it new cost value
                %BwE Because of bbob: use the FUN that is given as parameter and quit if fbest < ftarget
                [J(i,t+1), fbest, xbest, sBestJ, sBestX] = ...
                    evaluate_function(FUN, P(:,i,t+1), fbest, xbest, sBestJ, sBestX);
                if fbest < ftarget
                    return
                end								
			else  % It did not move up the gradient so stop the run for this bacterium
				m=Ns;
			end
		
		end	% Test if should end run for bacterium
        
        %BwE: Save best position of this run for this bacteria
        JbestX(i,t,:)=sBestX;
        JbestJ(i,t)=sBestJ;
        
        
        
	end  % Go to next bacterium
    
        x = P(1,:,t);
        y = P(2,:,t);
        clf    
        plot(x, y , 'h')   
        axis([-5 5 -5 5]);
        pause(.01)
        drawnow
	
%---------------------------------
%end  % j=1:Nc
%---------------------------------

	% Reproduction
	
    %BwE dont need to sum JHealth according to paper 3
    %We also simplify the “bacterial health” by the bacterium’s current fitness.
        
    % Set the health of each of the S bacteria.
    % There are many ways to define this; here, we sum
	% the nutrient concentrations over the lifetime of the bacterium.
	%J health=sum(J(:,:,k,ell),2);         
                                         
	% Sort cost and population to determine who can reproduce (ones that were in best nutrient
	% concentrations over their life-time reproduce)
	 
	[Jhealth,sortind]=sort(J(:,t+1)); % Sorts the nutrient concentration in order 
									% of ascending cost in the first dimension (bacteria number)
									% sortind are the indices in the new order
		
	P(:,:,t+1)=P(:,sortind,t+1); % Sorts the population in order of ascending Jhealth (the
											% ones that got the most nutrients were the ones with 
											% the lowest Jhealth values)
    JbestX(:,:,:)=JbestX(sortind,:,:);
    JbestJ(:,:)=JbestJ(sortind,:);
	
    %BwE I removed this, think we dont need it
	%C(:,k+1)=C(sortind,k); % And keeps the chemotaxis parameters with each bacterium at the next generation
		
	% Split the bacteria (reproduction)
	
	for b=1:Sr
		P(:,b+Sr,t+1)=P(:,b,t+1); % The least fit do not reproduce, the most 
		 									% fit ones split into two identical copies 
		%C(i+Sr,t+1)=C(i,t); 	% and they get the same parameters as for their mother
        JbestX(b+Sr,t,:)=JbestX(b,t,:);
        JbestJ(b+Sr,t)=JbestJ(b,t);
	end

	% Evolution can be added here (can add random modifications to C(i,k), Ns, Nc, etc)
	

%---------------------------------
%end  % k=1:Nre
%---------------------------------
	
	% Eliminate and disperse (on domain for our function) - keep same parameters C
	
	for m=1:S
		if ped>rand  % Generate random number and if ped bigger than it then eliminate/disperse
            %BwE change random number produced min -4 max 4
            P(:,m,t+1)=(ub*((2*round(rand(p,1))-1).*rand(p,1)));
		else
			P(:,m,t+1)=P(:,m,t+1);  % Bacteria that are not dispersed
		end
    end
    
    
    
    % TODO BwE: Implement this correctly
    % understand when and why it creates a small stepsize fbest < e(t)
    
    % Adaptation
    % Where t is the current generation number
    % fbest is the best fitness value among all the bacteria in the colony
    % etis the required precision in the current generation
    % and n, a, and ß are user-defined constants
%     if mod(t,n) == 0
%         if fbest < E(t)
%             C(:,t) = C(:,t-n)/alpha;
%             E(t) = E(t)/beta;
%         else
%             C(:,t) = C(:,t+1-n);
%             E(t+1) = E(t+1-n);
%         end
%     else
%         C(:,t+1) = C(:,t);
%         E(t+1) = E(t);
%     end
    

     if mod(t,n) == 0
         
        %Reinitialize colony position from the best-so-far position found by each bacterium    
        for b=1:S
            [~,idx]=sort(JbestJ(b,:));
            P(:,b,t+1)= JbestX(b,idx(1),:);
        end
        
        %if (min(JbestJ(:,t))-ftarget) < e(t)
        if (fbest-ftarget) < e(t)
            start = t-n+1;
            if start < 1
                start = 1;
            end                
            C(:,t+1) = C(:,start)/alpha;
            display(['smaller:' num2str(C(1,t+1)) ]) 
            e(t+1) = e(t)/beta;
        else
            start = t-n+1;
            if start < 1
                start = 1;
            end
            C(:,t+1) = C(:,start);
            display(['keep:' num2str(C(1,t+1)) ]) 
            e(t+1) = e(start);
        end
    else
        C(:,t+1) = C(:,t);
        e(t+1) = e(t);
     end    
%---------------------------------
end
%---------------------------------
    display(['fbest:' num2str(fbest) ])
%---------------------------------
%end  % ell=1:Ned
%---------------------------------


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Plot the function we are seeking the minimum of:
% %BwE range is [-5,5]
% 
% x=-5:5/100:5;   % For our function the range of values we are considering
% y=x;
% 
% % Compute the function that we are trying to find the minimum of.
% 
% for jj=1:length(x)
% 	for ii=1:length(y)
% %		z(ii,jj)=nutrientsfunc([x(jj);y(ii)],flag);
% %		z(ii,jj)=nutrientsfunc1([x(jj);y(ii)],flag);
% 		%BwE Because of bbob: use the FUN that is given as parameter
% 		%This is only for plotting so this should be removed.
% 		z(ii,jj)=feval(FUN, [x(jj);y(ii)]);
% 		
% 	end
% end
% 
% % First, show the actual function to be maximized
% 
% figure(1)
% clf
% surf(x,y,z);
% colormap(jet)
% % Use next line for generating plots to put in black and white documents.
% colormap(white);
% xlabel('x=\theta_1');
% ylabel('y=\theta_2');
% zlabel('z=J');
% title('Nutrient concentration (valleys=food, peaks=noxious)');
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Next, provide some plots of the results of the simulation.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% t=0:Nc;  % For use in plotting (makes t=0 correspond to the 1 index and plots to Nc+1)
% 
% % As is Figure 2 (4) shows parameter trajectories for 4 generations, then there is an elimination/dispersal event
% % Figure 3 (5) shows parameter trajectories for the following 4 generations (after the elim/disp event)
% 
% for kk=1:Ned
% figure(kk+1) 
% clf
% for mm=1:Nre
% subplot(2,2,mm)
% for nn=1:S % Plot all bacteria trajectories for generation mm
% plot(t,squeeze(P(1,nn,:,mm,kk)),t,squeeze(P(2,nn,:,mm,kk)))
% %plot(t,squeeze(P(1,nn,:,mm,kk)),'-',t,squeeze(P(2,nn,:,mm,kk)),'-')
% axis([min(t) max(t) -5 5])
% hold on
% end
% T=num2str(mm);
% T=strcat('Bacteria trajectories, Generation=',T);
% title(T)
% xlabel('Iteration, j')
% ylabel('\theta_1, \theta_2')
% hold off
% end
% end
% 
% for kk=1:Ned
% figure(Ned+kk+1) 
% clf
% for mm=1:Nre
% subplot(2,2,mm)
% contour(x,y,z,25)
% colormap(jet)
% for nn=1:S  % Plot all bacteria trajectories for generation mm
% hold on
% plot(squeeze(P(1,nn,:,mm,kk)),squeeze(P(2,nn,:,mm,kk)))
% %plot(squeeze(P(1,nn,:,mm,kk)),squeeze(P(2,nn,:,mm,kk)),'-')
% axis([-5 5 -5 5])
% end
% T=num2str(mm);
% T=strcat('Bacteria trajectories, Generation=',T);
% title(T)
% % Use next line for generating plots to put in black and white documents.
% %colormap(gray);
% xlabel('\theta_1');
% ylabel('\theta_2');
% hold off
% end
% end
% 
% %%%%%%%%%%%%
% %pause % Can leave this in if want to avoid movie (then hit control-C)
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Next, show a movie of the chemotactic steps:
% 
% figure(2*Ned+2) 
% clf
% 		contour(x,y,z,25)
% 		colormap(jet)
% 		axis([-5,5,-5,5]);
% 		xlabel('\theta_1');
% 		ylabel('\theta_2');
% 		title('Bacteria movements');
% hold on
% 
% M = moviein(Nc);
% 	for j=1:Nc;
% % Can change generation step and elimination-dispersal step on next line.
% % Currently for 1,1 - the first generation in the first elimination dispersal step
%         for i=1:S
% 		v=plot(squeeze(P(1,i,j:j+1,1,1)),squeeze(P(2,i,j:j+1,1,1)),'-');
% 		set(v,'MarkerSize',3);
% 		end
%         M(:,j)=getframe;
%     end;
% 	
% %movie(M,0)
% %save bacteria_swarm


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%