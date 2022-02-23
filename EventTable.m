    %% Sows
    %Evt 1, introduction of the sows in the beeding herd
    tIntroduction=((1:batchRearing)-1)*betweenBatch;
    batchNumberINTRO=mod(tIntroduction,durCycle)/betweenBatch+1;
    INTRO=[tIntroduction;batchNumberINTRO;ones(1,batchRearing)];
    %Evt 2, entrance in gestating room
    tGestatingRoom=durServ+((1:floor(MaxTime/betweenBatch))-1)*betweenBatch;
    batchNumberGEST=mod(0:length(tGestatingRoom)-1,batchRearing)+1;
    GEST=[tGestatingRoom;batchNumberGEST;2*ones(1,nBatch)];
    %Evt 3, entrance in farrowing room
    tFarrowingRoom=durGest+durServ+((1:floor(MaxTime/betweenBatch))-1)*betweenBatch;
    batchNumberFARROW=mod(0:length(tFarrowingRoom)-1,batchRearing)+1;
    FARROW=[tFarrowingRoom;batchNumberFARROW;3*ones(1,nBatch)];
    %Evt 4, return to service room
    tServiceRoom=durGest+durServ+durFarrow+((1:floor(MaxTime/betweenBatch))-1)*betweenBatch;
    batchNumberSERV=mod(0:length(tServiceRoom)-1,batchRearing)+1;
    SERV=[tServiceRoom;batchNumberSERV;4*ones(1,nBatch)];
    %Evt 33, piglet's birth
    tBirth=durGest+durServ+durFarrow-weaningAge+((1:floor(MaxTime/betweenBatch))-1)*betweenBatch;
    batchNumberBIRTH=mod(0:length(tBirth)-1,batchRearing)+1;
    BIRTH=[tBirth+0.5;batchNumberBIRTH;33*ones(1,nBatch)];
    %BIRTH=[tBirth;batchNumberBIRTH;33*ones(1,nBatch)];
    
    %Evt 5
    % Batch to batch vaccination VACC
     VACC=[tFarrowingRoom(tFarrowingRoom>debVacc)-7.9;batchNumberFARROW(tFarrowingRoom>debVacc);5*ones(1,length(tFarrowingRoom(tFarrowingRoom>debVacc)))];
    %Evt 6
    % Mass vaccination VACC
        
    MASS=[(debVacc:FrequenceVacc:MaxTime)+0.1;ones(1,length(debVacc:FrequenceVacc:MaxTime));6*ones(1,length(debVacc:FrequenceVacc:MaxTime))];
    
    tabEvSows=sortrows([INTRO GEST FARROW BIRTH SERV]',1)'; 
 
    
    %% Piglets
    %Evt 1, exit from Nursery
    tNursery=[tBirth+exitNursery;1:nBatch;ones(1,nBatch)];
            %Bandes de porcs entrant en PS
            indBatch=[tNursery(1,:)-exitNursery+weaningAge;tNursery(2,:)];
    
    %Evt 2, exit from Finishing Room
    tFinishingRoom=[BIRTH(1,:)+exitFinishingRoom;1:nBatch;2*ones(1,nBatch)];
    
%     %Evt 3, mass vaccination
%     MASS_piglets=[debVacc:FrequenceVacc:1600;ones(1,length(debVacc:FrequenceVacc:1600));3*ones(1,length(debVacc:FrequenceVacc:1600))];
   
    % enlever les bandes exportées des évènements porcelets
        Export=[find(tServiceRoom<debVacc+(1*betweenBatch) & tServiceRoom>debVacc),...
             find(tServiceRoom<debVacc+exitFinishingRoom+(1*betweenBatch) & tServiceRoom>debVacc+exitFinishingRoom)];
        tNursery(:,Export)=[];   
        tFinishingRoom(:,Export)=[];
        
    tabEvPigs=sortrows([tNursery tFinishingRoom ]',1)'; %MASS_piglets
    
