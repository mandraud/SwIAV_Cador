function [popPigsTime] = Piglets_MSIRS_Program(durCycle,numInfectedSeeds,betAir,Np,ChangePiglets,Time,totalBatch,tabEvPigs,...
    popPigsTime,transitions_piglets,betaNoAir,alphaBatch,debVacc,propM)
% Do the iterations using the full event driven stochastic methodology
% This is a relatively general version of the event-driven methodology know
% as Gillespie's Direct Algorithm
 
%     firstActiveBatch = find(any(any(popPigsTime(:,:,:,Time(1)))),1);
    active_batches = find(any(any(popPigsTime(:,:,:,Time(1)))),totalBatch); %firstActiveBatch:totalBatch;
    oldPop = popPigsTime(:,:,active_batches,Time(1));

    % reduce the 'totalBatch' to the total number of current batches
    totalBatch = size(oldPop,3);
    
    loop=1;
    currentTime = Time(1);

    popPigsTime_ind = floor(Time(1));   
 
    % FIX: if Time(1) from tableEvSows is also present in tableEvPigs
    if any(tabEvPigs(1,:) == Time(1)),
        currentTime = currentTime-0.0001;
    end
                   
    % prepare transition values
    transitions_rep = repmat(transitions_piglets,[1,2,totalBatch]);
    Npigs= sum(max(oldPop>0,[],1),3)*Np;
    Npigs = max(repmat(Npigs-Np,[1,1,totalBatch]),0);
    betaAir_rep = betAir./Npigs;
    betaNoAir_rep = repmat(betaNoAir,[1,2,totalBatch]);
   

    
    infectPress=CorridorPressure(Np,oldPop);%popPigsTime(:,:,:,floor(currentTime)));
                %infections S->I1 ou S->I2
               if any(infectPress>0)
                    MI=zeros(7,1);
                    SI=sum(rand(oldPop(8,1,totalBatch),1) <= 1-exp(-betAir*sum(infectPress)));
                    for i=1:7
                        MI(i)=sum(rand(oldPop(i,1,totalBatch),1) <= 1-exp(-betAir*sum(infectPress)));
                    end
                    MI1=round(MI*infectPress(1)/sum(infectPress)).*(oldPop(1:7,1,totalBatch)>0);
                    MI2=MI-MI1;
                    I1=round(SI*infectPress(1)/sum(infectPress)).*(oldPop(8,1,totalBatch)>0);
                    I2=SI-I1;
                    I1=I1+sum(MI1);
                    I2=I2+sum(MI2);
                    
                    transfer= [I1;I2;...
                         sum(rand(oldPop(9,1,totalBatch),1) <= 1-exp(-betAir*infectPress(2)));...
                         sum(rand(oldPop(10,1,totalBatch),1) <= 1-exp(-betAir*infectPress(2)));...
                         sum(rand(oldPop(13,1,totalBatch),1) <= 1-exp(-betAir*infectPress(1)));...
                         sum(rand((oldPop(14,1,totalBatch)),1) <= 1-exp(-betAir*infectPress(1)))];

                    oldPop([1:7 8 9 10 13 14],1,totalBatch)= oldPop([1:7 8 9 10 13 14],1,totalBatch)-[MI;SI;transfer(3:6)];
                    oldPop([9 13 11 27 15 28],1,totalBatch)= oldPop([9 13 11 27 15 28],1,totalBatch)+transfer;
  
    popPigsTime(:,:,active_batches,Time(1))=oldPop;
    

end
    % perform time loop
    while (currentTime<Time(2)) 

        [step,newPop]=Iterate(oldPop,ChangePiglets,totalBatch,transitions_rep,betaAir_rep,betaNoAir_rep,active_batches,alphaBatch,propM,Np);

        % if step == 0, go to MaxTime 
        if(step == 0),
            step = Time(2) - currentTime;
        end
      
        loop=loop+1;
        previousTime = currentTime;
        currentTime = min(currentTime+step,Time(2));

        IsEvent=(previousTime<tabEvPigs(1,:) & currentTime>tabEvPigs(1,:));
        if any(IsEvent)
            ind=find(IsEvent==1,1);
            batch=find(active_batches==tabEvPigs(2,ind),1); % - (firstActiveBatch-1);
            currentTime=tabEvPigs(1,ind);
            Event=tabEvPigs(3,ind);

            % start from previous population
            newPop = oldPop;
            
            switch Event
                %%% Exit from nursery to finishing room with corridor pressure
                case 1  
                newPop(:,2,batch)= newPop(:,1,batch); 
                newPop(:,1,batch)=0; 

                %corridor pressure
                infectPress=CorridorPressure(Np,oldPop);%popPigsTime(:,:,:,floor(currentTime)));
                %infections S->I1 ou S->I2
               if any(infectPress>0)
                    SI=sum(rand(oldPop(8,1,batch),1) <= 1-exp(-betAir*sum(infectPress)));
                    I1=round(SI*infectPress(1)/sum(infectPress));
                    I2=SI-I1;
                    transfer= [I1;I2;...
                         sum(rand(oldPop(9,1,batch),1) <= 1-exp(-betAir*infectPress(2)));...
                         sum(rand(oldPop(10,1,batch),1) <= 1-exp(-betAir*infectPress(2)));...
                         sum(rand(oldPop(13,1,batch),1) <= 1-exp(-betAir*infectPress(1)));...
                         sum(rand((oldPop(14,1,batch)),1) <= 1-exp(-betAir*infectPress(1)))];

                    newPop([8 9 10 13 14],2,batch)= newPop([8 9 10 13 14],2,batch)-[SI;transfer(3:6)];
                    newPop([9 13 11 27 15 28],2,batch)= newPop([9 13 11 27 15 28],2,batch)+transfer;
                end
%                 transfertST2= [sum(rand(oldPop(8,1,batch),1) <= 1-exp(-betAir*infectPress(2)));...
%                     sum(rand(oldPop(9,1,batch),1) <= 1-exp(-betAir*infectPress(1)));...
%                     sum(rand(oldPop(10,1,batch),1) <= 1-exp(-betAir*infectPress(2)))];
%                 newPop([8 9 10],2,batch)= newPop([8 9 10],2,batch)-transfertST2.*(newPop([8 9 10],2,batch)>0);
%                 newPop([13 11 12],2,batch)= newPop([13 12],2,batch)+transfertST2.*(newPop([8 9 10],2,batch)>0);
         
%                 % batch to batch vaccination
%             if currentTime>=debVacc
%                 newPop(18,2,batch)=sum(newPop([8 10 14 17],2,batch),1);
%                 newPop(21,2,batch)=newPop(7,2,batch);
%                 newPop(22,2,batch)=newPop(6,2,batch);
%                 newPop(23,2,batch)=newPop(5,2,batch);
%                 newPop(24,2,batch)=newPop(4,2,batch);
%                 newPop([4:7 8 10 14 17],2,batch)=0;
%             end
                % update Npigs and betaAir_rep
                Npigs = sum(max(newPop>0,[],1),3)*Np;
                Npigs = repmat(Npigs-Np,[1,1,totalBatch]);
                betaAir_rep = betAir/Npigs;

                %%% Exit from finishing room to slaughterhouse
                case 2
                newPop(:,2,batch)=0;

                % update Npigs and betaAir_rep
                Npigs= sum(max(newPop>0,[],1),3)*Np;
                Npigs = repmat(Npigs-Np,[1,1,totalBatch]);
                betaAir_rep = betAir/Npigs;
                
                
%                 % Vaccination de masse
%                 case 3
%                     newPop(18,1:2,:)=sum(newPop([8 10 14 17],1:2,:),1);
%                     newPop(21,1:2,:)=newPop(7,1:2,:);
%                     newPop(22,1:2,:)=newPop(6,1:2,:);
%                     newPop(23,1:2,:)=newPop(5,1:2,:);
%                     newPop(24,1:2,:)=newPop(4,1:2,:);
%                     newPop([4:7 8 10 14 17],1:2,:)=0;
            end
        end

        oldPop = newPop;
          if any(oldPop(:,:,:)<0)
                   totalBatch;
               end             
        % logging
        if currentTime > popPigsTime_ind,

            if(currentTime > popPigsTime_ind+1),
                % time step > 1
                num_steps = floor(currentTime) - popPigsTime_ind+1;
                popPigsTime(:,:,active_batches,popPigsTime_ind:floor(currentTime)) = repmat(newPop,[1,1,1,num_steps]);
                popPigsTime_ind = popPigsTime_ind + num_steps;
            else 
                % time step == 1
                popPigsTime(:,:,active_batches,popPigsTime_ind) = newPop;
                popPigsTime_ind = popPigsTime_ind + 1;
            end
        end

    end
    
    % last time step
    popPigsTime(:,:,active_batches,floor(currentTime)) = newPop;



function [step, new_value]=Iterate(oldp,ChangePiglets,totalBatch,transitions_rep,betaAir_rep,betaNoAir_rep,active_batches,alphaBatch,propM,Np)
%% Do the actual iteration step

    lambda=FOI(oldp,totalBatch,betaAir_rep,betaNoAir_rep,propM,Np);
%     transitions_rep(8:14,:,:)
% if any(lambda(2,:)>0)
%     display(lambda);
% end
    transitions_rep(8:14,:,:) = (repmat((alphaBatch(1,1,active_batches)+alphaBatch(1,2,active_batches)),7,2).* ...
        transitions_rep(8:14,:,:)+ repmat(alphaBatch(1,3,active_batches),7,2));
%     transitions_rep(8:14,:,:)
    transitions_rep(15:21,:,:) = (repmat((alphaBatch(1,3,active_batches)+alphaBatch(1,2,active_batches)),7,2).* ...
        transitions_rep(15:21,:,:)+ repmat(alphaBatch(1,1,active_batches),7,2));
    
    transitions_rep([8:14 22 30 31 41:47],:,:) = ...
        transitions_rep([8:14 22 30 31 41:47],:,:).* repmat(lambda(1,:,:),[17,1,1]); %Subtype n°1
    transitions_rep([15:21 24 25 28 48:54],:,:) = ...
        transitions_rep([15:21 24 25 28 48:54],:,:).* repmat(lambda(2,:,:),[17,1,1]); %Subtype n°2
%  if any(lambda)>0
%      lambda
%  end
    states=[1:7,1:7,1:7,8:10,9,11:12,8,13,14,13,15,16,18:24,18:24,18:24,25,26,27,28];
    Rate = transitions_rep .* oldp(states,:,:);
    
    new_value=oldp;
    [step, j,k,l] = GetEvent(Rate);  
    if step>0
        new_value(:,k,l) = new_value(:,k,l) + ChangePiglets(j,:)'; 
    end
  
 
%% Calculate the force of infection    
function lambda=FOI(oldp,totalBatch,betaAir_rep,betaNoAir_rep,propM,Np)
betaNoAir_rep([11 12],:,:)=1.1*(3.6*propM+(1-propM))/Np;
lambdaNoAirST1 = sum(betaNoAir_rep([1:5 11],:,:).*oldp([9 11 15 16 28 25],:,:),1);
lambdaNoAirST2 = sum(betaNoAir_rep([6:10 12],:,:).*oldp([13 11 15 12 27 26],:,:),1);

lambdaAirST1  = max(betaAir_rep.*sum(oldp([9 11 15 16 25 28],:,:),1),0); % total pressure ST1
if any(oldp(25,:,:)>0)
oldp(25,:,:)
end
% a1=lambdaAirST1;
lambdaAirST1(:,1,:)=sum(lambdaAirST1(:,1,:),3)-lambdaAirST1(:,1,:);
lambdaAirST1(:,2,:)=sum(lambdaAirST1(:,2,:),3)-lambdaAirST1(:,2,:);
% lambdaAirST1  = repmat(sum(lambdaAirST1,3),[1,1,totalBatch]) - lambdaAirST1; % remove own batch pressure


lambdaAirST2  = max(betaAir_rep.*sum(oldp([13 11 15 12 26 27],:,:),1),0); % total pressure ST2
lambdaAirST2(:,1,:)=sum(lambdaAirST2(:,1,:),3)-lambdaAirST2(:,1,:);
lambdaAirST2(:,2,:)=sum(lambdaAirST2(:,2,:),3)-lambdaAirST2(:,2,:);
% lambdaAirST2  = repmat(sum(lambdaAirST2,3),[1,1,totalBatch]) - lambdaAirST2; % remove own batch pressure

lambda = [lambdaNoAirST1+lambdaAirST1;lambdaNoAirST2+lambdaAirST2];

% extra check required in case Npigs == 0
lambda(isnan(lambda)) = 0;

