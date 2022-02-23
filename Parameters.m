 % Model parameters   
 
batchRearing=7;

FrequenceVacc=90;

betAir=0.1; %0.01 0.5 or 1
sigmaS=1/600; %120, 240 or 360
sigmaS12=1/90;
sigmaP=sigmaS;

betaS1=2.43; %1.15; %Rose et al., 2013 106-127 jours à l'infection pour H1avN1
betaS2=2.43; %3.22; %betaS1; %1.52; %Rose et al., 2013 38-56 jours à l'infection pour H1huN2=0.57 quand H1avN1=0.43 donc rapport=1.33
%betaS12=0.54; %Rose et al., 2013 42-50 jours à l'infection pour H1avN1*H1huN2*rH1avN2

betaP1=2.43; %2.43; %Essai
betaP2=2.43; %3.22; % 2.43 * rapport 1.33

gammaS1=1/8.5;%6.1; 
gammaS2=gammaS1;%6.1; %1/8.5; %Rose et al., 2013 38-56 jours à l'infection pour H1huN2=10.4 quand H1avN1=7.5 => rapport=1.39
gammaS12=gammaS1;%6.1; %1/7.6; %Rose et al., 2013 42-50 jours à l'infection pour H1avN1*H1huN2*rH1avN2
gammaRI=1/6;
gammaP1=gammaS1;
gammaP2=gammaS2;
gammaP12=gammaS12;
nExportBatch=0;
sigmaM=1/70.0;

%% Vaccination
betaV=.28; %Romagosa: de 0.09 à 0,53
sigmaV=1/105; %Truies essai animalerie
gammaV=1/6.0; %Romagosa= 3.5 jours pour les vaccinés, 4.5 pour les non-vaccinés


eps(1:7,1,1)=[0.024,0.13,0.33,0.49,0.72,0.87,0.94]; 

mu=0;%.385/2.35;

if batchRearing==4
   weaningAge=21;
   Ns=round(192/batchRearing); %192 Nb of sows in the total herd
   betweenBatch=35;
   durFarrow=weaningAge+7;
   durServ=33;
   durGest=79; %115-28-5
   
   exitNursery=82; 
   exitFinishingRoom=182; 
end
if batchRearing==5
   weaningAge=21;
    Ns=round(247/batchRearing); 
   betweenBatch=28;
   durFarrow=weaningAge+6;
   durServ=34;
   durGest=79;
   
   exitNursery=75; 
   exitFinishingRoom=182; 
end
if batchRearing==7
   weaningAge=28;
    Ns=round(203/batchRearing);
   betweenBatch=21;
   durFarrow=weaningAge+7;
   durServ=34;
   durGest=78;
   
   exitNursery=68; 
   exitFinishingRoom=182; 
end
if batchRearing==10
   weaningAge=21;
    Ns=round(433/batchRearing);
   betweenBatch=14;
   durFarrow=weaningAge+6;
   durServ=34;
   durGest=79;
   
   exitNursery=75; 
   exitFinishingRoom=182;
end
if batchRearing==20
   weaningAge=21;
    Ns=round(615/batchRearing);
   betweenBatch=7;
   durFarrow=weaningAge+6;
   durServ=33;
   durGest=80;

   exitNursery=68; 
   exitFinishingRoom=182;
end

% durCycle=durFarrow+durServ+durGest;

% durFarrow=weaningAge+5;
% durServ=32;
% durGest=82; %115-28-5
    %115-28-5
   
% exitNursery=77; 
% exitFinishingRoom=182;  

if weaningAge==21
   durCycle=140;
else
   durCycle=147;
end


debVacc=4*durCycle+90;


numInfectedSeeds = [1 0];

litterSize=12; 
Np=round(Ns*litterSize);

MaxTime=2440;

Initial=zeros(53,3,batchRearing);
Initial(1,1,1)=Ns;

nBatch=floor(MaxTime/betweenBatch); 

%% Pre-process the transition rates.
n = 7;
Eps=repmat(eps,[1 2]);

% PIGLETS
transitions_piglets = ones(58,1);

% waning maternally immunity M => M_i => S
transitions_piglets(1:7)  = sigmaM*n;
% infection rate 
transitions_piglets(8:14) = eps(1:7); % .lambda M_i => I1
transitions_piglets(15:21) = eps(1:7); % .lambda M_i => I2

% infection rate S => I1
transitions_piglets(22) = 1; % .lambda // Sensibilité du porc=1
% I1 => R1
transitions_piglets(23) = gammaP1;
% infection rate R1 => R1I2
transitions_piglets(24) = 1; % .lambda // Sensibilité du porc=1
% infection rate I1 => I12
transitions_piglets(25) = 1; % .lambda // Sensibilité du porc=1
% recovery I12 => R1I2
transitions_piglets(26) = gammaP2*2; % I2 partagé entre les transitions 26 et 27
% recovery R1I2 => R12
transitions_piglets(27) = gammaP2*2;

% infection rate S => I2
transitions_piglets(28) = 1; % .lambda // Sensibilité du porc=1
% I2 => R2
transitions_piglets(29) = gammaP2;
% R2 => R2I1
transitions_piglets(30) = 1; % .lambda // Sensibilité du porc=1
% infection rate I2 => I21
transitions_piglets(31) = 1; % .lambda // Sensibilité du porc=1
% infection rate I21 => R2I1
transitions_piglets(32) = gammaP1*2; % I1 partagé entre les transitions 32 et 33
% recovery R2I1 => R12
transitions_piglets(33) = gammaP1*2;

% vaccine immunity waning
transitions_piglets(34:40)  = sigmaV*n;
% infection rate while having vaccine
transitions_piglets(41:47) = 1; % .lambda V_i => I1
transitions_piglets(48:54) = 1; % .lambda V_i => I2
% recovery rate
transitions_piglets(55:56) = gammaV;
transitions_piglets(57:58) = [gammaP1 gammaP2];
% SOWS
transitions_sows  = ones(84,1);

% infection rate S => I1
transitions_sows(1) = 1; % .lambda // Sensibilité de la truie=1
% I1 => R1
transitions_sows(2) = gammaS1;
% waning immunity R1 => R1_i => S
transitions_sows(3:9) = sigmaS*n;
% infection rate I1 => I12
transitions_sows(10) = 1; % .lambda // Sensibilité de la truie=1
% I12 => R1I2 => R12
transitions_sows(11) = gammaS2*2; % I2 partagé entre les transitions 11 et 83
transitions_sows(83) = gammaS2*2; % I2 partagé entre les transitions 11 et 83

transitions_sows(19) = gammaS2;
% infection rate R1_i => R1I2
transitions_sows(12:18) = 1; % .lambda // Sensibilité de la truie=1

% infection rate S => I2
transitions_sows(20) = 1; % .lambda // Sensibilité de la truie=1
% I2 => R2
transitions_sows(21) = gammaS2;
% waning immunity R2 => R2_i => S
transitions_sows(22:28) = sigmaS*n;
% infection rate I2 => I21
transitions_sows(29) = 1; % .lambda // Sensibilité de la truie=1
% I21 => R2I1 => R12
transitions_sows(30) = gammaS1*2; % I1 partagé entre les transitions 30 et 84
transitions_sows(84) = gammaS1*2; % I1 partagé entre les transitions 30 et 84

transitions_sows(38) = gammaS1;
% infection rate R2_i => R2I1
transitions_sows(31:37) = 1; % .lambda // Sensibilité de la truie=1

% waning immunity R12 => R12_i => S
transitions_sows(39:45)=sigmaS*n;

% vaccine immunity waning
transitions_sows(46:52)  = sigmaV*n;
% infection rate while having vaccine
transitions_sows(53:59) = 1; % .lambda V_i => I1
transitions_sows(60:66) = 1; % .lambda V_i => I2
% Recovery rate
transitions_sows(67:68) = gammaV;

% waning immunity R12 => R12_i => R2_4
transitions_sows(69:75)=sigmaS12*n;
% waning immunity R12 => R12_i => R1_4
transitions_sows(76:82)=sigmaS12*n;

transitions_sows_rep = repmat(transitions_sows,[1,3,batchRearing]);


transitions_farrowing = [transitions_piglets ; transitions_sows];

%% Pre-process the beta values.

% beta no air: PIGLETS
betaNoAir = [betaP1 ; betaP1 ; betaP1 ; betaP1; betaP1 ; betaP2 ; betaP2 ; betaP2 ; betaP2; betaP2; betaV; betaV] / Np; 


% beta FARROW
betaFarrowST1 = [betaP1; betaP1; betaP1;  betaP1; betaV; betaP1;
    betaS1; betaS1; betaS1; betaS1; betaV; betaS1]/(Ns+Np);
betaFarrowST2 = [betaP2; betaP2; betaP2; betaP2; betaV; betaP2;
    betaS2; betaS2; betaS2; betaS2; betaV; betaS2]/(Ns+Np);

%% Pre-process for proportion of M piglets
alphaBatch=zeros(1,3,nBatch);
