# SwIAV_Cador
The model is made of three component: 
  - the driving one is the dynamics of the breeding herd, dersigned in Sows_SIRS_MainProgram
  - From here, piglets are generated and evolve in two environments: 
  the farrowing sector in Farrowing_Room_Program, followed by the post-weaning/farrowing sector governed by Piglets_MSIRS_Program
 
 All events related to batch-farrowing systems are set deteministic in EventTable and parameters are set in Parameters
 Stochastic epidemiological events are defined using Gillespie algorithm with transition matrix defined in ChangeTable

To run: execute Main 
