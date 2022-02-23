% Transition changes for the Sows, Piglets and Farrowing room

ChangeSows=zeros(84,53);
ChangeSows(1:9,1:9)=diag(-ones(9,1),0)+diag(ones(8,1),1); % S=> I1 => R1{1:7}
ChangeSows(9,1)=1; % R1{7} => S
ChangeSows(10,2)=-1; ChangeSows(10,10)=1; % I1 => I12
ChangeSows(11,10)=-1; ChangeSows(11,52)=1; % I12 =>R1I2
ChangeSows(12:18,3:9)=diag(-ones(7,1),0); % R1{1:7} => R1I2
ChangeSows(12:18,11)=1; % R1{1:7} => R1I2
ChangeSows(19,11)=-1; ChangeSows(19,38)=1; % R1I2 => R12
ChangeSows(20,1)=-1; ChangeSows(20,12)=1; % S => I2
ChangeSows(21:28,12:19)=diag(-ones(8,1),0)+diag(ones(7,1),1); % I2 => R2{1:7}
ChangeSows(28,1)=1; % R2{7} => S
ChangeSows(29,12)=-1; ChangeSows(29,20)=1; % I1 => I21
ChangeSows(30,20)=-1; ChangeSows(30,53)=1; % I21 => R2I1
ChangeSows(31:37,13:19)=diag(-ones(7,1),0); % R2{1:7} => R2I1
ChangeSows(31:37,21)=1; % R2{1:7} => R2I1
ChangeSows(38,21)=-1; ChangeSows(38,45)=1; % R1I2 => R12
ChangeSows(39:45,22:28)=diag(-ones(7,1),0)+diag(ones(6,1),1); % R12{1:7} => S
ChangeSows(45,1)=1; % R12{7} => S
ChangeSows(46:52,29:35)=diag(-ones(7,1),0)+diag(ones(6,1),1); % V12{1:7} => S
ChangeSows(52,1)=1; % V12{1:7} => S
ChangeSows(53:59,29:35)=diag(-ones(7,1),0); % V12{1:7} => V12I1
ChangeSows(53:59,36)=1; % V12{1:7} => V12I1
ChangeSows(60:66,29:35)=diag(-ones(7,1),0); % V12{1:7} => V12I2
ChangeSows(60:66,37)=1; % V12{1:7} => V12I2
ChangeSows(67,36)=-1; ChangeSows(67,3)=1; % V12I1 => R12
ChangeSows(68,37)=-1; ChangeSows(68,22)=1; % V12I2 => R12
ChangeSows(69:75,38:44)=diag(-ones(7,1),0)+diag(ones(6,1),1); % R12{1:7} => R2{4}
ChangeSows(75,16)=1; % R12{7} => R2{4}
ChangeSows(76:82,45:51)=diag(-ones(7,1),0)+diag(ones(6,1),1); % R12{1:7} => R1{4}
ChangeSows(82,6)=1; % R12{7} => R1{4}

ChangeSows(83,52)=-1;
ChangeSows(84,53)=-1;
ChangeSows(83,3)=1;
ChangeSows(84,22)=1;



ChangePiglets=zeros(58,28);
ChangePiglets(1:7,1:7)=diag(-ones(7,1),0)+diag(ones(6,1),1); % M{1:7} => S
ChangePiglets(7,8)=1; % M7 => S
ChangePiglets(8:14,1:7)=diag(-ones(7,1),0); % M{1:7} => I1
ChangePiglets(15:21,1:7)=diag(-ones(7,1),0); % M{1:7} => I2
ChangePiglets(8:14,9)=1; % M{1:7} => I1
ChangePiglets(15:21,13)=1; % M{1:7} => I2
ChangePiglets(22:24,8:10)=diag(-ones(3,1),0)+diag(ones(2,1),1); % S => I1 => R1 => R1I2
ChangePiglets(24,27)=1; % R1 => R1I2
ChangePiglets(25,9)=-1; ChangePiglets(25,11)=1; % I1 => I12
ChangePiglets(26,11)=-1; ChangePiglets(26,12)=1; % I12 => R1I2
ChangePiglets(27,12)=-1; ChangePiglets(27,17)=1; % R1I2 => R12
ChangePiglets(28,8)=-1; ChangePiglets(28,13)=1; % S => I2
ChangePiglets(29:30,13:14)=diag(-ones(2,1),0)+diag(ones(1,1),1); % I2 => R2 => R2I1
ChangePiglets(30,28)=1; % R2 => R2I1
ChangePiglets(31,13)=-1; ChangePiglets(31,15)=1; % I2 => I21
ChangePiglets(32,15)=-1; ChangePiglets(32,16)=1; % I21 => R2I1
ChangePiglets(33,16:17)=[-1 1]; % R2I1 => R12
ChangePiglets(34:40,18:24)=diag(-ones(7,1),0)+diag(ones(6,1),1); % V12{1:7} => S
ChangePiglets(40,8)=1; % V12 => S
ChangePiglets(41:47,18:24)=diag(-ones(7,1),0); % V12{1:7} => V12I1
ChangePiglets(41:47,25)=1; % V12{1:7} => V12I1
ChangePiglets(48:54,18:24)=diag(-ones(7,1),0); % V12{1:7} => V12I2
ChangePiglets(48:54,26)=1; % V12{1:7} => V12I2
ChangePiglets(55,25)=-1; ChangePiglets(55,10)=1; % V12I1 => R12
ChangePiglets(56,26)=-1; ChangePiglets(56,17)=1; % V12I2 => R12
ChangePiglets(57,27)=-1; ChangePiglets(57,17)=1; % V12I2 => R12
ChangePiglets(58,28)=-1; ChangePiglets(58,17)=1; % V12I2 => R12

ChangeFarrow=[ChangePiglets zeros(58,53);zeros(84,28) ChangeSows];
