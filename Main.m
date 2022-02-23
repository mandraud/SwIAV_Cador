% Main program.

clear all;

Parameters; 

%%
tic
nbIter = 20 ;
for j=1:nbIter

    rng('shuffle');
    [popPigsTime,popSowsTime,popFarrowTime] = Sows_SIRS_MainProgram;
    j2=num2str(j);
    %Outputs;
    
    v3 = strcat('outputs/Mu0Sigma600/SIVAC',j2,'.mat');
    save(v3);
end
toc

