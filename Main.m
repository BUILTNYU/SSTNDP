experiment = 2; % 1: 5-by-5 grid, 2: NYC PUMA

if experiment == 1
    % PART FOR GRID EXPERIMENT
    load('Exp1_Input.mat');
    M = 1;      % # of trials per segment expansion
    Run = 10;  % # of simulations
    policy = 300;   % 100: MAB, 200: KG, 300: KGCB
    nReqRoute = 3;  % # of routes to be operated
    maxL = 4;       % # of nodes 
    d = 2;      % level of demand (1-10)
    u = 1;      % urban environment setting (1: random, 2: monocentric, 3: polycentric)
    
    [OptChc1,Result1,TotFlow1,T1,ErrRate1]=SSTNDP_Exp1_Grid(M,Run,policy,nReqRoute,maxL,PreDefTruth,GridInput,d,u);
else
    % PART FOR NYC PUMA EXPERIMENT
    DataList = {'Exp2_Input_1.mat','Exp2_Input_2.mat','Exp2_Input_3.mat','Exp2_Input_4.mat','Exp2_Input_5.mat'};
    Data = DataList{1}; % choose one dataset from DataList
    M = 1;
    Run = 10;
    level = 1;      % variability level (1: low, 2: medium, 3: high)
    policy = 400;    % 100: MAB, 200: KG, 300: KGCB, 400: Greedy
    nReqRoute = 5;
    maxL = 8;
    
    [OptChc2,Result2,TotFlow2,T2,ErrRate2] = SSTNDP_Exp2_NYC(M,Run,level,policy,nReqRoute,maxL,Data);
end