 function [popFarrowTime,alpha] = Farrowing_Room_Program(popSows,litterSize,ChangeFarrow,Time,indexPiglets,...
     transitions_farrowing,popFarrowTime,betaFarrowST1,betaFarrowST2,Ns,Np)
% Do the iterations using the full event driven stochastic methodology
% This is a relatively general version of the event-driven methodology 
% known as Gillespie's Direct Algorithm

   % Initiate piglet vector
   initPiglets=zeros(28,1);
   % recovered and non-primary infected sows => maternal immunity
%    initPiglets(1)=round(sum(popSows([3:9 11 13:19 21:28 29:35]))*litterSize);

   %make correspondance between sows' immunity waning and MDA delivery
   initPiglets(1)=round(sum(popSows([3 11 13 21 22 29 38 45 52 53]))*litterSize); 
   initPiglets(2)=round(sum(popSows([4    14    23 30 39 46]))*litterSize);
   initPiglets(3)=round(sum(popSows([5    15    24 31 40 47]))*litterSize);
   initPiglets(4)=round(sum(popSows([6    16    25 32 41 48]))*litterSize);
   initPiglets(5)=round(sum(popSows([7    17    26 33 42 49]))*litterSize);
   initPiglets(6)=round(sum(popSows([8    18    27 34 43 50]))*litterSize);
   initPiglets(7)=round(sum(popSows([9    19    28 35 44 51]))*litterSize);

   % susceptible and primary infected sows => susceptible piglets
   initPiglets(8)=round(sum(popSows([1:2 10 12 20 36 37]))*litterSize);
   
   alpha=zeros(1,3);
   RSows= sum(popSows([3:9 11 13:19 21:28 29:35 38:44 45:51 52 53]));
   if RSows>0
       alpha(1)=sum(popSows([3:9 11 52]))/RSows;           %proportion of Sows immune to serotype 1
       alpha(2)=sum(popSows([22:28 29:35 38:44 45:51]))/RSows;      %proportion of Sows immune to both serotypes
       alpha(3)=sum(popSows([13:19 21 53]))/RSows;         %proportion of Sows immune to serotype 2 
   
   end
   
    propM=sum(initPiglets(2:7))/sum(initPiglets); 

   % add sows
   oldPop=[initPiglets; popSows];
   loop=1;
   current_time = Time(1);
   popFarrowTime_ind = floor(current_time);
 
    while (current_time<Time(2))  
        [step,newPop]=Iterate(oldPop,ChangeFarrow,transitions_farrowing,betaFarrowST1,betaFarrowST2,alpha,propM,Ns,Np);
        
        loop=loop+1;
        current_time=min(current_time+step,Time(2));
 
%         if current_time==Time(2)
%             newPop = oldPop;
%         end
        
        % prepare for new iteration
        oldPop=newPop;

        % logging
        if current_time > popFarrowTime_ind,
                % time step > 1
            if(current_time > popFarrowTime_ind+1),
                num_steps = floor(current_time) - popFarrowTime_ind+1;
                popFarrowTime(:,indexPiglets,popFarrowTime_ind:floor(current_time)) = repmat(newPop,[1,1,num_steps]);
                popFarrowTime_ind = popFarrowTime_ind + num_steps;
            else
                % time step == 1
                popFarrowTime(:,indexPiglets,popFarrowTime_ind) = newPop;
                popFarrowTime_ind = popFarrowTime_ind + 1;
            end
        end
    end

    % last time step
    popFarrowTime(:,indexPiglets,floor(current_time)) = newPop;

%     sum(popFarrowTime(1:26,indexPiglets,floor(current_time)) )
%% Do the actual iteration step
function [step, new_pop]=Iterate(old_pop,ChangeFarrow,transitions_farrowing,betaFarrowST1,betaFarrowST2,alpha,propM,Ns,Np)
    
%    betaFarrowST1([5,11])=0.9*3.8*propM+0.9*(1-propM)/(Ns+Np);
 %   betaFarrowST2([5,11])=0.9*3.8*propM+0.9*(1-propM)/(Ns+Np);

    lambda=zeros(1,2);
    lambda(1) = sum(betaFarrowST1 .* old_pop([9,11,15,16,25,28,28+[2,10,20,21,36,53]]));
    lambda(2) = sum(betaFarrowST2 .* old_pop([13,11,15,12,26,27,28+[10,11,12,20,37,52]]));
if lambda(2)>0
    lambda;
end;
    transitions_farrowing(8:14) = ((alpha(1)+alpha(2)).* ...
        transitions_farrowing(8:14)+alpha(3));
    transitions_farrowing(15:21) = ((alpha(3)+alpha(2)).* ...
        transitions_farrowing(15:21)+alpha(1));

    K1 = [8:14 22 30 31 41:47 58+[1 29 31:37 53:59]]; %58 transitions chez les porcs, 84 chez les truies
    K2 = [15:21 24 25 28 48:54 58+[10 12:18 20 60:66]]; %58 transitions chez les porcs, 84 chez les truies
        
    transitions_farrowing(K1) = transitions_farrowing(K1) * lambda(1);
    transitions_farrowing(K2) = transitions_farrowing(K2) * lambda(2);
    states=[1:7,1:7,1:7,8:10,9,11,12,8,13:14,13,15,16,18:24,18:24,18:24,25,26,27,28,...
        28+[1:9,2,10,3:9,11,1,12:19,12,20,13:19,21,22:28,29:35,29:35,29:35,36,37,38:44,45:51,52,53]];
    Rate = transitions_farrowing .* old_pop(states);
  
    sumRates=sum(Rate);
        
    step = -log(rand(1))/sumRates;
if step<0
    display(step);
end
    % find which event to do 
    m=find(cumsum(Rate)>=rand(1)*sumRates, 1);

    new_pop=old_pop;
    if sumRates>0
        new_pop=new_pop + ChangeFarrow(m,:)';
    end