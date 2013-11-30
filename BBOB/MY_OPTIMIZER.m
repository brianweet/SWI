function xbest = MY_OPTIMIZER(FUN, DIM, ftarget, maxfunevals)
% MY_OPTIMIZER(FUN, DIM, ftarget, maxfunevals)
% samples new points uniformly randomly in [-5,5]^DIM
% and evaluates them on FUN until ftarget of maxfunevals
% is reached, or until 1e8 * DIM fevals are conducted. 

%   maxfunevals = min(1e8 * DIM, maxfunevals); 
%   popsize = min(maxfunevals, 200);
%   fbest = inf;
%   for iter = 1:ceil(maxfunevals/popsize)
%     xpop = 10 * rand(DIM, popsize) - 5;      % new solutions
%     [fvalues, idx] = sort(feval(FUN, xpop)); % evaluate
%     if fbest > fvalues(1)                    % keep best
%       fbest = fvalues(1);
%       xbest = xpop(:,1);
%     end
%     if feval(FUN, 'fbest') < ftarget         % COCO-task achieved
%       break;                                 % (works also for noisy functions)
%     end
%   end 

p=DIM;                       % dimension of search space 
%s=26;                        % The number of bacteria 
s=20;
Nc=50;                       % Number of chemotactic steps 
Ns=4;                        % Limits the length of a swim 
Nre=4;                       % The number of reproduction steps 
Ned=2;                       % The number of elimination-dispersal events 
Sr=s/2;                      % The number of bacteria reproductions (splits) per generation 
Ped=0.25;                    % The probabilty that each bacteria will be eliminated/dispersed 
c(:,1)=0.05*ones(s,1);       % the run length  
for m=1:s                    % the initital posistions 
    P(1,:,1,1,1)= 5*rand(s,1)';
    P(2,:,1,1,1)= 5*rand(s,1)';
end                                                                  
     
%%
%Main loop 
    

%Elimination and dispersal loop 
for ell=1:Ned
    

%Reprodution loop


    for K=1:Nre    

%  swim/tumble(chemotaxis)loop   

        for j=1:Nc
            
            for i=1:s        
                J(i,j,K,ell) = feval(FUN, P(:,i,j,K,ell));
                %Live_fn(P(:,i,j,K,ell));         

% Tumble

                        
                Jlast=J(i,j,K,ell);   
                Delta(:,i)=(2*round(rand(p,1))-1).*rand(p,1); 	             	
                P(:,i,j+1,K,ell)=P(:,i,j,K,ell)+c(i,K)*Delta(:,i)/sqrt(Delta(:,i)'*Delta(:,i)); % This adds a unit vector in the random direction            
 
% Swim (for bacteria that seem to be headed in the right direction)     
                
                J(i,j+1,K,ell)=...
                    feval(FUN, P(:,i,j+1,K,ell));
                %Live_fn(P(:,i,j+1,K,ell));  
                m=0;         % Initialize counter for swim length 
                    while m<Ns     
                          m=m+1;
                          if J(i,j+1,K,ell)<Jlast  
                             Jlast=J(i,j+1,K,ell);    
                             P(:,i,j+1,K,ell)=P(:,i,j+1,K,ell)+c(i,K)*Delta(:,i)/sqrt(Delta(:,i)'*Delta(:,i)) ;  
                             J(i,j+1,K,ell)= ...
                                 feval(FUN, P(:,i,j+1,K,ell));
                             % Live_fn(P(:,i,j+1,K,ell));  
                          else       
                             m=Ns ;     
                          end        
                    
                    end 
                J(i,j,K,ell)=Jlast;
                sprintf('The value of interation i %3.0f ,j = %3.0f  , K= %3.0f, ell= %3.0f' , i, j, K ,ell );
                   
            end % Go to next bacterium
            
            x = P(1,:,j,K,ell);
            y = P(2,:,j,K,ell);
            clf    
            plot(x, y , 'h')   
            axis([-5 5 -5 5]);
            pause(.1)
        end  % Go to the next chemotactic    

                 
%Reprodution                                              
        Jhealth=sum(J(:,:,K,ell),2);              % Set the health of each of the S bacteria
        [Jhealth,sortind]=sort(Jhealth);          % Sorts the nutrient concentration in order of ascending 
        P(:,:,1,K+1,ell)=P(:,sortind,Nc+1,K,ell); 
        c(:,K+1)=c(sortind,K);                    % And keeps the chemotaxis parameters with each bacterium at the next generation
                                     

%Split the bacteria (reproduction)                             
            for i=1:Sr
                P(:,i+Sr,1,K+1,ell)=P(:,i,1,K+1,ell); % The least fit do not reproduce, the most fit ones split into two identical copies  
                c(i+Sr,K+1)=c(i,K+1);                 
            end   
        end %  Go to next reproduction    


%Eliminatoin and dispersal
        for m=1:s 
            if  Ped>rand % % Generate random number 
                P(1,:,1,1,1)= 50*rand(s,1)';
                P(2,:,1,1,1)= .2*rand(s,1)';  
            else 
                P(:,m,1,1,ell+1)=P(:,m,1,Nre+1,ell); % Bacteria that are not dispersed
            end        
        end 
        
        reproduction = J(:,1:Nc,Nre,ell);
        [jlastreproduction,O] = min(reproduction,[],2);  % min cost function for each bacterial 
        [~,I] = min(jlastreproduction);
        xbest=P(:,I,O(I,:),K,ell);
        
        if xbest < ftarget
            break;
        end  
    end % Go to next elimination and disperstal 
  
