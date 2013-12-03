function xbest = bfoa2(FUN, DIM, ftarget, maxfunevals,n,alpha,beta)
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
max_gen=1000;

S=100;	
Ns=4;
Sr=S/2;
ped=0.25;

e=zeros(1,max_gen);
e(1,1)=100;
% n = 10;
% alpha = 2;
% beta = 5;

% Initial population
P(:,:,:)=0*ones(p,S,max_gen);  % First, allocate needed memory

% BBOB search space bounds
ub = 4;

% Randomly place on domain:
for m=1:S
    P(:,m,1)=(ub*((2*round(rand(p,1))-1).*rand(p,1)));
end
C=0*ones(S,max_gen);
runlengthunit=0.1;
C(:,1)=runlengthunit*ones(S,1);


% Allocate memory for cost function:
J=0*ones(S,max_gen);
Jhealth=0*ones(S,1);

% BBOB variables
fbest = inf;
xbest = zeros(2,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(1) 
% clf
% 		%contour(x,y,z,25)
% 		%colormap(jet)
% 		%axis([1,max_gen,-5,5]);
% 		
% hold on

%---------------------------------
% Swim/tumble (chemotaxis) loop:
%---------------------------------
Ja(:,:)=inf(S,max_gen+1);
Xa(:,:,:)=zeros(S,max_gen+1,p);
for t=1:max_gen

	for i=1:S  % For each bacterium
		
		% Compute the nutrient concentration at the current location of each bacterium
		
		%BwE Because of bbob: use the FUN that is given as parameter and quit if fbest < ftarget
        [J(i,t), fbest, xbest] = evaluate_function2(FUN, P(:,i,t), fbest, xbest);
        if fbest < ftarget
            return
        end        
        if J(i,t) < Ja(i,t)
            Ja(i,t)=J(i,t);
            Xa(i,t,:)=P(:,i,t);
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
        [J(i,t+1), fbest, xbest] = evaluate_function2(FUN, P(:,i,t+1), fbest, xbest);
        if fbest < ftarget
            return
        end
        if J(i,t+1) < Ja(i,t+1)
            Ja(i,t+1)=J(i,t+1);
            Xa(i,t+1,:)=P(:,i,t+1);
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
                [J(i,t+1), fbest, xbest] = evaluate_function2(FUN, P(:,i,t+1), fbest, xbest);
                if fbest < ftarget
                    return
                end	
                if J(i,t+1) < Ja(i,t+1)
                    Ja(i,t+1)=J(i,t+1);
                    Xa(i,t+1,:)=P(:,i,t+1);
                end
			else  % It did not move up the gradient so stop the run for this bacterium
				m=Ns;
			end
		
		end	% Test if should end run for bacterium
	end  % Go to next bacterium

	% Reproduction
	
    % Set the health of each of the S bacteria.
    % There are many ways to define this; here, we sum
	% the nutrient concentrations over the lifetime of the bacterium.
	Jhealth=sum(J(:,:),2);  
    
    [Jhealth,sortind]=sort(J(:,t+1));
                                         
	% Sort cost and population to determine who can reproduce (ones that were in best nutrient
	% concentrations over their life-time reproduce)
	
	%[Jhealth,sortind]=sort(Jhealth); % Sorts the nutrient concentration in order 
									% of ascending cost in the first dimension (bacteria number)
									% sortind are the indices in the new order
		
	P(:,:,t+1)=P(:,sortind,t+1); % Sorts the population in order of ascending Jhealth (the
											% ones that got the most nutrients were the ones with 
											% the lowest Jhealth values)
		
	% Split the bacteria (reproduction)
	for i=1:Sr
		P(:,i+Sr,t+1)=P(:,i,t+1); % The least fit do not reproduce, the most 
		 									% fit ones split into two identical copies 
	end

	% Eliminate and disperse (on domain for our function) - keep same parameters C
    for m=1:S
        if ped>rand  % Generate random number and if ped bigger than it then eliminate/disperse
            P(:,m,t+1)=(ub*((2*round(rand(p,1))-1).*rand(p,1)));
        else
			%P(:,m,1,1,ell+1)=P(:,m,1);  % Bacteria that are not dispersed
        end
    end
    
    
    % Adaptation
    % Where t is the current generation number
    % fbest is the best fitness value among all the bacteria in the colony
    % etis the required precision in the current generation
    % and n, a, and ß are user-defined constants
    ta = t+1;
    if mod(t,n) == 0
        
        %Reinitialize colony position from the best-so-far position found by each bacterium
        for b=1:S
            [~,idx]=sort(Ja(b,:));
            P(:,b,t+1)= Xa(b,idx(1),:);
        end
        
        if fbest < e(t)
            C(:,t+1) = C(:,ta-n)/alpha;
            e(t+1) = e(ta-n)/beta;
%             display(['low:' num2str(C(1,t+1))])
        else
            C(:,t+1) = C(:,ta-n);
            e(t+1) = e(ta-n);
%             display(['keep:' num2str(C(1,t+1))])
        end
    else
        C(:,t+1) = C(:,t);
        e(t+1) = e(t);
    end
        
    
% v=semilogy(t,fbest,'-');
% set(v,'MarkerSize',3);
% drawnow   
%---------------------------------
end  % t=1:max_gen
%---------------------------------
display(['End:' num2str(fbest) ' stepsize:' num2str(C(1,max_gen)) ...
    ' n:' num2str(n) ' a:' num2str(alpha) ' b:' num2str(beta)])
	


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