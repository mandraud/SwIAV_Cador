function [popPigsTime,popSowsTime,popFarrowTime] = Sows_SIRS_MainProgram(betaS1, betaS2, betaV,...
      betAir,betweenBatch,weaningAge,batchRearing,litterSize,Ns,Np,...
      mu,durCycle,Initial,MaxTime,nBatch)
% Setup and start the event driven simulator

% Sets up default parameters if necessary.
if nargin == 0; 
  Parameters;
end  
n=0;
EventTable;
ChangeTable;
    
% The main iteration 

    [popPigsTime,popSowsTime,popFarrowTime] = Stoch_Iteration(MaxTime,Initial,betaS1, betaS2, betaV,...
        betAir,betweenBatch,weaningAge,batchRearing,litterSize,Ns,Np,mu,durCycle,...
        ChangeSows,ChangeFarrow,ChangePiglets,tabEvSows,tabEvPigs,indBatch,...
        nBatch,transitions_piglets,betaNoAir,...
        transitions_sows_rep,transitions_farrowing,betaFarrowST1,betaFarrowST2,numInfectedSeeds,debVacc);


% STOCH ITERATION    
% Do the iterations using the full event driven stochastic methodology
% This is a relatively general version of the event-driven methodology know
% as Gillespie's Direct Algorithm
function [popPigsTime,popSowsTime,popFarrowTime]=Stoch_Iteration(MaxTime,Initial,betaS1, betaS2, betaV,...
      betAir,betweenBatch,weaningAge,batchRearing,...
    litterSize,Ns,Np,mu,durCycle,ChangeSows,ChangeFarrow,ChangePiglets,tabEvSows,tabEvPigs,indBatch,...
    nBatch,transitions_piglets,betaNoAir,...
    transitions_sows_rep,transitions_farrowing,betaFarrowST1,betaFarrowST2,numInfectedSeeds,debVacc)

   popPigsTime = zeros(28,2,nBatch,MaxTime);
   popSowsTime = zeros(53,3,batchRearing,MaxTime);
   popFarrowTime = zeros(28+53,nBatch,MaxTime);
   alphaBatch=zeros(1,3,nBatch);
   
    indexPiglets=0; 
%     totalBatch=0;

    popSowsTime(:,:,:,1) = Initial;
    oldPop=Initial;

    loop=1;
    currentTime = 0;
    popSowsTime_ind = 1;
   
while (currentTime < MaxTime)  
   
    [step,newPop]=Iterate(oldPop,betaS1,betaS2,betaV,Ns,ChangeSows,transitions_sows_rep);
    
    % if step == 0, go to MaxTime 
     if(step == 0),
         step = MaxTime - currentTime;
     end
     
    loop=loop+1;
    previousTime = currentTime;
    currentTime  = min(currentTime+step,MaxTime);
       
    % Get first event from previousTime till currentTime
    IsEvent=(previousTime<tabEvSows(1,:) & currentTime>tabEvSows(1,:));
    if any(IsEvent)
        ind=find(IsEvent==1,1);
   
        numBatch=tabEvSows(2,ind);
%         numBatch
        currentTime=tabEvSows(1,ind);
        Event=tabEvSows(3,ind);
        
        % start from previous population
        newPop = oldPop;
        
        switch Event
         %%% Introduction of new sows in the service room
         case 1  
            newPop(1,1,numBatch) = Ns;

         %%% Entrance in gestating room
         case 2  
            newPop(:,2,numBatch)=oldPop(:,1,numBatch); 
            newPop(:,1,numBatch)=0; 

         %%% Entrance in farrowing room
         case 3
         %   Infectious pressure from the piglets on susceptible sows passing through the corridor  
            infectPress=CorridorPressure(Np,popPigsTime(:,:,:,floor(currentTime)));
            if any(infectPress>0)
                SI=sum(rand(oldPop(1,2,numBatch),1) <= 1-exp(-betAir*sum(infectPress)));
                I1=round(SI*infectPress(1)/sum(infectPress));
                I2=SI-I1;
                Transfer=[I1;I2;...
                         sum(rand(oldPop(12,2,numBatch),1)<=1-exp(-betAir*infectPress(1)));...
                         sum(rand(oldPop(13,2,numBatch),1)<=1-exp(-betAir*infectPress(1)));...
                         sum(rand(oldPop(14,2,numBatch),1)<=1-exp(-betAir*infectPress(1)));...
                         sum(rand(oldPop(15,2,numBatch),1)<=1-exp(-betAir*infectPress(1)));...
                         sum(rand(oldPop(16,2,numBatch),1)<=1-exp(-betAir*infectPress(1)));...
                         sum(rand(oldPop(17,2,numBatch),1)<=1-exp(-betAir*infectPress(1)));...
                         sum(rand(oldPop(18,2,numBatch),1)<=1-exp(-betAir*infectPress(1)));...
                         sum(rand(oldPop(19,2,numBatch),1)<=1-exp(-betAir*infectPress(1)));...
                         sum(rand(oldPop(2,2,numBatch),1)<=1-exp(-betAir*infectPress(2)));...
                         sum(rand(oldPop(3,2,numBatch),1)<=1-exp(-betAir*infectPress(2)));...
                         sum(rand(oldPop(4,2,numBatch),1)<=1-exp(-betAir*infectPress(2)));...
                         sum(rand(oldPop(5,2,numBatch),1)<=1-exp(-betAir*infectPress(2)));...
                         sum(rand(oldPop(6,2,numBatch),1)<=1-exp(-betAir*infectPress(2)));...
                         sum(rand(oldPop(7,2,numBatch),1)<=1-exp(-betAir*infectPress(2)));...
                         sum(rand(oldPop(8,2,numBatch),1)<=1-exp(-betAir*infectPress(2)));...
                         sum(rand(oldPop(9,2,numBatch),1)<=1-exp(-betAir*infectPress(2)))]; 
                 
            newPop([1 12 13:19 2 3:9],2,numBatch)=newPop([1 12 13:19 2 3:9],2,numBatch)-[SI;Transfer(3:18)];
            newPop([2 12 20 21 10 11],2,numBatch)=newPop([2 12 20 21 10 11],2,numBatch)+[Transfer(1:3);sum(Transfer(4:10));Transfer(11);sum(Transfer(12:18))];
            end

            newPop(:,3,numBatch)=newPop(:,2,numBatch); 
            newPop(:,2,numBatch)=0; 

         %%% Birth of the piglets
         case 33
            indexPiglets=indexPiglets+1;
            Time=[currentTime min(MaxTime,currentTime+weaningAge)];
            newPop(:,3,numBatch)=0;        

            [popFarrowTime,alpha] = Farrowing_Room_Program(oldPop(:,3,numBatch),litterSize,ChangeFarrow,Time,...
                indexPiglets,transitions_farrowing,popFarrowTime,betaFarrowST1,betaFarrowST2,Ns,Np);
            
    
        %%% Weaning
        %%% Piglets go to the nursery
         case 4
            % Increment the total number of piglet batches
            indice=find((indBatch(1,:)==currentTime),1);
            totalBatch=indBatch(2,indice);
                
            % Calculation of the infectious pressure from the piglets on 
            % susceptible sows passing through the corridor
            infectPress=CorridorPressure(Np,popPigsTime(:,:,:,currentTime));


            % copy piglets into piglet matrix
%             ExportEvent=(currentTime<debVacc+(nExportBatch*betweenBatch) & currentTime>debVacc);
%             if any(ExportEvent)
%                 popFarrowTime(1:28,totalBatch,currentTime)=0;
%                 n=n+1;
%             end
%             ExportEvent=(currentTime<debVacc+182+(1*betweenBatch) & currentTime>debVacc+182);
%             if any(ExportEvent)
%                 popFarrowTime(1:28,totalBatch,currentTime)=0;
%             end
            
            popPigsTime(:,1,totalBatch,currentTime) = popFarrowTime(1:28,totalBatch,currentTime);
            alphaBatch(1,1:3,totalBatch)=alpha(1,1:3);
            popFarrowTime(1:28,totalBatch,currentTime)=0;
            
       propM=sum(popPigsTime(1:7,1,totalBatch,currentTime))/sum(popPigsTime(:,1,totalBatch,currentTime));
%             % batch to batch vaccination at weaning
%        if currentTime >= debVacc && currentTime <= (debVacc+(5*betweenBatch))
%             popPigsTime(18,1,totalBatch,currentTime)=sum(popPigsTime([7 8 10 14 17],1,totalBatch,currentTime),1);
%             popPigsTime(19:24,1,totalBatch,currentTime)=popPigsTime(6:-1:1,1,totalBatch,currentTime);
% %             popPigsTime(22,1,totalBatch)=popPigsTime(6,1,totalBatch);
% %             popPigsTime(23,1,totalBatch)=popPigsTime(5,1,totalBatch);
% %             popPigsTime(24,1,totalBatch)=popPigsTime(4,1,totalBatch);
%             popPigsTime([1:7 8 10 14 17],1,totalBatch)=0;
%       end

            % process piglet events during 2 batches
            Time=[currentTime min(MaxTime,currentTime+betweenBatch)]; 
            [popPigsTime] = Piglets_MSIRS_Program(durCycle,numInfectedSeeds,betAir,Np,ChangePiglets,...
                Time,totalBatch,tabEvPigs,popPigsTime,transitions_piglets,betaNoAir,alphaBatch,debVacc,propM);

            % copy the sow population from popFarrow back into popSow
            newPop(:,3,numBatch)=popFarrowTime(29:81,totalBatch,Time(1));
            popFarrowTime(29:81,totalBatch,Time(1))=0;
            
            % Vaccination au retour des truies de Maternité
            VaccEvent=(Time(1)-weaningAge<tabEvSows(1,tabEvSows(3,:)==6) & Time(1)>tabEvSows(1,tabEvSows(3,:)==6));
            if any(VaccEvent)
                newPop(29,3,numBatch)=sum(newPop([1 3:9 13:19 22:28 38:44 45:51],3,numBatch),1);
                newPop([1 3:9 13:19 22:28 38:44 45:51],3,numBatch)=0;
            end
          
            % Sows get back to service room: infection process  
            %We should also include the FOI exerted on (and by) vaccinated sows
            if any(infectPress>0)
                    SI=sum(rand(newPop(1,3,numBatch),1) <= 1-exp(-betAir*sum(infectPress)));
                    I1=round(SI*infectPress(1)/sum(infectPress));
                    I2=SI-I1;
            Transfer=[I1;I2;...
                     sum(rand(newPop(12,3,numBatch),1)<=1-exp(-betAir*infectPress(1)));...
                     sum(rand(newPop(13,3,numBatch),1)<=1-exp(-betAir*infectPress(1)));...
                     sum(rand(newPop(14,3,numBatch),1)<=1-exp(-betAir*infectPress(1)));...
                     sum(rand(newPop(15,3,numBatch),1)<=1-exp(-betAir*infectPress(1)));...
                     sum(rand(newPop(16,3,numBatch),1)<=1-exp(-betAir*infectPress(1)));...
                     sum(rand(newPop(17,3,numBatch),1)<=1-exp(-betAir*infectPress(1)));...
                     sum(rand(newPop(18,3,numBatch),1)<=1-exp(-betAir*infectPress(1)));...
                     sum(rand(newPop(19,3,numBatch),1)<=1-exp(-betAir*infectPress(1)));...
                     sum(rand(newPop(2,3,numBatch),1)<=1-exp(-betAir*infectPress(2)));...
                     sum(rand(newPop(3,3,numBatch),1)<=1-exp(-betAir*infectPress(2)));...
                     sum(rand(newPop(4,3,numBatch),1)<=1-exp(-betAir*infectPress(2)));...
                     sum(rand(newPop(5,3,numBatch),1)<=1-exp(-betAir*infectPress(2)));...
                     sum(rand(newPop(6,3,numBatch),1)<=1-exp(-betAir*infectPress(2)));...
                     sum(rand(newPop(7,3,numBatch),1)<=1-exp(-betAir*infectPress(2)));...
                     sum(rand(newPop(8,3,numBatch),1)<=1-exp(-betAir*infectPress(2)));...
                     sum(rand(newPop(9,3,numBatch),1)<=1-exp(-betAir*infectPress(2)))]; 
                 
            newPop([1 12 13:19 2 3:9],3,numBatch)=newPop([1 12 13:19 2 3:9],3,numBatch)-[SI;Transfer(3:18)];
            newPop([2 12 20 21 10 11],3,numBatch)=newPop([2 12 20 21 10 11],3,numBatch)+[Transfer(1:3);sum(Transfer(4:10));Transfer(11);sum(Transfer(12:18))];
            end
             
            % Switch locations + replacement process
            Cull=zeros(53,1);
            for k=1:53
                Cull(k)=sum(rand(newPop(k,3,numBatch),1)<=mu);
            end
           
            if sum(Cull)==0
                q=find(newPop(:,3,numBatch)>0);
                Cull(q(randi(size(q,1))))=1; %insure the introduction of at least one gilt
            end
            newPop(:,1,numBatch)=newPop(:,3,numBatch)-Cull; 
            newPop(:,3,numBatch)=0; 
            newPop(1,1,numBatch)=newPop(1,1,numBatch)+sum(Cull);
            
            % Seed 1 infected sow in the first batch on 'durCycle' (ST1)
            if (currentTime>=3*durCycle && numBatch==1 && numInfectedSeeds(1)>0)
                newPop(1,1,numBatch)=newPop(1,1,numBatch)-numInfectedSeeds(1);
                newPop(2,1,numBatch)=newPop(2,1,numBatch)+numInfectedSeeds(1);
                numInfectedSeeds(1)=0;
            end
            if (currentTime>=4*durCycle && numBatch==1 && numInfectedSeeds(2)>0)
                newPop(1,1,numBatch)=newPop(1,1,numBatch)-numInfectedSeeds(2);
                newPop(12,1,numBatch)=newPop(12,1,numBatch)+numInfectedSeeds(2);
                numInfectedSeeds(2)=0;
            end
            
         % Vaccination bande à bande
         case 5
            newPop(29,2,numBatch)=sum(newPop([1 3:9 13:19 22:28 38:44 45:51],2,numBatch),1);
            newPop([1 3:9 13:19 22:28 38:44 45:51],2,numBatch)=0;
         
         % Vaccination de masse
         case 6
            newPop(29,1:3,1:batchRearing)=sum(newPop([1 3:9 13:19 22:28 38:44 45:51],1:3,1:batchRearing),1);
            newPop([1 3:9 13:19 22:28 38:44 45:51],1:3,1:batchRearing)=0;
        end
    end
    
%     % prepare new step
%     oldPop = newPop;

     % logging
    if currentTime > popSowsTime_ind,

        if(currentTime > popSowsTime_ind+1),
            % time step > 1
            num_steps = floor(currentTime) - popSowsTime_ind+1;
            popSowsTime(:,:,:,popSowsTime_ind:(floor(currentTime)-1)) = repmat(oldPop,[1,1,1,num_steps-1]);
            popSowsTime(:,:,:,floor(currentTime)) = newPop;
            popSowsTime_ind = popSowsTime_ind + num_steps;
        else 
            % time step == 1
            popSowsTime(:,:,:,popSowsTime_ind) = newPop;
            popSowsTime_ind = popSowsTime_ind + 1;
            
        end
     end
    % prepare new step
    oldPop = newPop;
    
end

%% Do the actual iteration step
function [step, new_value]=Iterate(oldPop,betaS1,betaS2,betaV,Ns,ChangeSows,transitions_sows_rep)

    Nsows = sum(sum(oldPop,3));
    lambda=zeros(2,3);
    % lambda in service rooms
    % betawS/2 because animals are housed in individual rearing system >> less contacts between animals
    lambda(1,1)=(betaS1/2*sum(sum(oldPop([2 21 53 10 20],1,:),3))  + betaV/2*sum(oldPop(36,1,:)))/Nsows(1);
    lambda(2,1)=(betaS2/2*sum(sum(oldPop([11 12 52 10 20],1,:),3)) + betaV/2*sum(oldPop(37,1,:)))/Nsows(1);
    
    %lambda in gestating room
    if Nsows(2)>0,
       lambda(1,2)=(betaS1*sum(sum(oldPop([2 21 53 10 20],2,:),3)) + betaV*sum(oldPop(36,2,:)))/Nsows(2);
       lambda(2,2)=(betaS2*sum(sum(oldPop([11 12 52 10 20],2,:),3)) +  betaV*sum(oldPop(37,2,:)))/Nsows(2); 
    end

    % lambda in farrowing room => individual batches?!
     lambda(1,3)=(betaS1*sum(sum(oldPop([2 21 53 10 20],3,:),3)) + betaV*sum(oldPop(36,3,:))) /Ns;
     lambda(2,3)=(betaS2*sum(sum(oldPop([11 12 52 10 20],3,:),3))+ betaV*sum(oldPop(37,3,:)))/Ns;

    % adapt transitions with force of infection
    transitions_sows_rep([1 29 31:37 53:59],1,:) = transitions_sows_rep([1 29 31:37 53:59],1,:)* lambda(1,1);
    transitions_sows_rep([10 12:18 20 60:66],1,:) = transitions_sows_rep([10 12:18 20 60:66],1,:)* lambda(2,1);
    transitions_sows_rep([1 29 31:37 53:59],2,:) = transitions_sows_rep([1 29 31:37 53:59],1,:)* lambda(1,2);
    transitions_sows_rep([10 12:18 20 60:66],2,:) = transitions_sows_rep([10 12:18 20 60:66],1,:)* lambda(2,2);

    transitions_sows_rep([1 29 31:37 53:59],3,:) = transitions_sows_rep([1 29 31:37 53:59],1,:)* lambda(1,3);
    transitions_sows_rep([10 12:18 20 60:66],3,:) = transitions_sows_rep([10 12:18 20 60:66],1,:)* lambda(2,3);
%     if any(lambda>0)
%         display('stop')
%     end
    states=[1:9,2,10,3:9,11,1,12:19,12,20,13:19,21,22:28,29:35,29:35,29:35,36,37,38:44,45:51,52,53];

    % apply transition rates
    Rate = transitions_sows_rep(:,:,1).*sum(oldPop(states,:,:),3);
    
    new_value=oldPop;
    [step,j,k] = GetEventSows(Rate);

    if step >0, 
%         a=sum(oldPop); 
        b=find(new_value(states(j),k,:)>0); %4 5 6
        c=size(b,1); %3
        l=b(randi(c));  %b(1) ou b(2) ou b(3)   
        new_value(:,k,l)=new_value(:,k,l) + ChangeSows(j,:)';
    end
    

     
 
    
   
    