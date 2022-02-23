function infectPress=CorridorPressure(Np,popPigsTime_step)
 % Infection pressure from the piglets on susceptibles passing through the corridor 

 % Number of piglet batches at each animal passing in the corridor
 NP = sum(sum(any(popPigsTime_step),3));
 % Number of infectious piglets at that time
 infectedPigs = [sum(sum(sum(popPigsTime_step([9 11 15 16 25 28],:,:))));...
     sum(sum(sum(popPigsTime_step([11 12 13 15 26 27],:,:))))];
 % infection pressure
 infectPress =infectedPigs/(NP*Np);
 if isnan(infectPress(1))
     infectPress(1)=0;
 end
 if isnan(infectPress(2))
     infectPress(2)=0;
 end