function [OptChc,Result,TotFlow,T,ErrRate]=SSTNDP_Exp2_NYC(M,Run,level,policy,nReqRoute,maxL,Data) % with saved pilot
load(Data);

if level==1
    TrCorrCov=TrueCorrPairCovL;
    TrCorrMean=TrueCorrPairMeanL;
    TrUncorr=TrueUncorrPairL;
    FlowCluster=FlowClusterL;
    PilotPre=PilotPreL;
elseif level==2
    TrCorrCov=TrueCorrPairCovM;
    TrCorrMean=TrueCorrPairMeanM;
    TrUncorr=TrueUncorrPairM;
    FlowCluster=FlowClusterM;
    PilotPre=PilotPreM;
else
    TrCorrCov=TrueCorrPairCovH;
    TrCorrMean=TrueCorrPairMeanH;
    TrUncorr=TrueUncorrPairH;
    FlowCluster=FlowClusterH;
    PilotPre=PilotPreH;
end

T=zeros(Run,1); % time calculation
Result=struct;  % prepare the space for archiving results

% Truth for 2 flow groups (uncorrelated vs correlated)
% TrueUncorrFlow: mean and std.dev
% TrueCorrFlowMean: mean
% TrueCorrFlowCov: covariance

%% Separate PUMAinfo by whether OD is correlated or not
m=0; n=0;
for i=1:nOD
    if ismember(PUMAinfo(i,1),TrUncorr(:,1))
        m=m+1;
        PUMAinfoUncorr(m,:)=PUMAinfo(i,:);
    else
        n=n+1;
        PUMAinfoCorr(n,:)=PUMAinfo(i,:);
    end
end


PairperNode=cell(nNode,2);
for i=1:size(PUMAinfoUncorr,1)
    PairperNode{PUMAinfoUncorr(i,2),1}(end+1)=PUMAinfoUncorr(i,1);
    PairperNode{PUMAinfoUncorr(i,3),1}(end+1)=PUMAinfoUncorr(i,1);
end
for i=1:size(PUMAinfoCorr,1)
    PairperNode{PUMAinfoCorr(i,2),2}(end+1)=PUMAinfoCorr(i,1);
    PairperNode{PUMAinfoCorr(i,3),2}(end+1)=PUMAinfoCorr(i,1);
end

%%
mabsetting=0; kggeneral=0; kgcb=0; grdy=0;
if policy==100
    mabsetting=1;	% 1 for Multi-armed bandit algorithm, 0 otherwise
elseif policy==200
    kggeneral=1;    % 1 for Knowledge gradient, 0 otherwise
elseif policy==300
    kgcb=1;         % 1 for Knowledge gradient with correlated belief, 0 otherwise
else % no policy: greedy algorithm
    grdy=1;
end
assert((mabsetting + kggeneral + kgcb + grdy)==1);
    
%% REPEAT EXPERIMENTS
%{ 
STEPS
1. Designate N_11, the initial node of the 1st route R1, from the
network node set       
2. Extend link by link from N_11
3. When reaching maxL, finalize R1
For 1<k<=nReqRoute
4. Designate N_k1, the initial node of the k-th route Rk, from the
system node set (to increase transferability)
5. Extend link by link from N_k1
6. When reaching maxL, finalize Rk
7. k=k+1
8. When k>nReqRoute, stop.
%}
for i=1:Run
    Result(i).Err=[0,0];
end

nSim=0; % # of experiments conducted
while nSim<Run
    try
    tic
    %% PREPARE SPACE FOR RESULT AND ARCHIVES
    EstMeanODFlow=[];   % Estimated mean OD flow of all OD pairs
    
    % LINK ARCHIVE % Use "IncludedLink" instead
    IncludedLink=[];        % space for links included in the system (empty at the beginning)
    RemainingLink=EdgeInfo; % links remaining uncovered (all links remain at the beginning)
    
    % NODE ARCHIVE
    NodeRoute = []; % archive nodes included in route during each route installation
    NodeSystemNet = []; % archive all nodes included (without grouping)
    
    % OD PAIR ARCHIVE
    ODSystemNet=[]; % ODs served by routes
    
    % LIST OF ROUTES VISITING A NODE
    NodeRouteList=cell(nNode,2);    % routes that serves a node
    for i=1:nNode
        NodeRouteList{i,1}=i; % 1st column: node ID, 2nd column: routes visiting i-th node
    end
    
    % INDICATOR OF TRANSFERABILITY
    RTintsct=eye(nReqRoute);	% if RTintsct(i,j)=1: i and j are transferable
    RTGroups=cell(nReqRoute,2); % group of routes mutually transferable
    for i=1:nReqRoute
        RTGroups{i,1}=i; % 1st column: route ID, 2nd column: transferable routes from i-th route
    end
    
    % OPPORTUNITY COST
    % Opportunity cost can be observed by subtracting the value of chosen
    % link from that of the best link 
    % (Only under oracle setting. In reality, it's impossible.)
    OPPCOST={};	% archive for opportunity costs per extension
    OppCost=cell(maxL,5);	% empty cell for opportunity cost [route,chosen node,choice made,true best,opportunity cost]
    
    %% RECALL GENERATED PRIORS
    PrUncorr=PilotPre.PrUncorrFlow;
    PrCorrMean=PilotPre.PrCorrFlowMean;
    PrCorrCov=PilotPre.PrCorrFlowCov;
    ObsOD=PilotPre.ObsOD;               % observed OD flows during pilots
    PrecisU=PilotPre.PrecisU;
    PrecisC=PilotPre.PrecisC;
    nUncorrFlow=size(PrUncorr,1);
    nCorrFlow=size(PrCorrMean,1);
    
    %% DESIGNATE THE INITIAL NODE OF THE 1ST ROUTE AND BEGIN THE EXTENSION
    RT=1; % begin with the 1st route
    OppCost{1,1}=RT;    % opportunity cost of the 1st route will be archived
    
    % LARGEST IN- AND OUTFLOW TO ALL OTHER NODES (FROM PILOT)
    NodeNetFlow=zeros(1,nNode);	% aggregated in- and outflow of all nodes
    for i=1:nUncorrFlow
        af=PrUncorr(PrUncorr(:,1)==PUMAinfoUncorr(i,1),2);
        NodeNetFlow(PUMAinfoUncorr(i,2))=NodeNetFlow(PUMAinfoUncorr(i,2))+af;
        NodeNetFlow(PUMAinfoUncorr(i,3))=NodeNetFlow(PUMAinfoUncorr(i,3))+af;
    end
    for i=1:nCorrFlow
        af=PrCorrMean(PrCorrMean(:,1)==PUMAinfoCorr(i,1),2);
        NodeNetFlow(PUMAinfoCorr(i,2))=NodeNetFlow(PUMAinfoCorr(i,2))+af;
        NodeNetFlow(PUMAinfoCorr(i,3))=NodeNetFlow(PUMAinfoCorr(i,3))+af;
    end
    [~,start]=max(NodeNetFlow(:));	% node with maximum total flow
    if policy<0 % choose starting point randomly
        start=randi(nNode);
    end
    NodeRoute=start;                % add this node to 1st route
    OppCost{1,2}=start;             % 1st choice is "start"
    NodeSystemNet=start;            % list of nodes covered by system
    NodeRouteList{start,2}=RT;      % route "RT" covers node "start"
        
    % IDENTIFY AVAILABLE LINKS FROM THE INITIAL NODE
    LinkChoice1=[find(EdgeInfo(2,:)==start),find(EdgeInfo(3,:)==start)]; % find links connected to the node
    nCandLink=size(LinkChoice1,2);          % # of available links
    ODChoice=cell(1,nCandLink);             % cell for total OD covered if extending to a link
    ODChoiceIndUncorr=zeros(nUncorrFlow,nCandLink);
    ODChoiceIndCorr=zeros(nCorrFlow,nCandLink);

    % PREPARE TRUTH AND PRIOR OF MEAN AND STD.DEV OF REWARD FROM CHOICES
    Mu=zeros(nCandLink,1);      % mean flow based on truth
    Mu0=zeros(nCandLink,1);     % mean flow based on prior
    Sigma=zeros(nCandLink,1);	% std.dev of flow based on truth
    Sigma0=zeros(nCandLink,1);	% std.dev of flow based on prior
    for i=1:nCandLink
        % Since there are only two nodes, one OD pair exists
        ODTemp=EdgeInfo(2:3,LinkChoice1(i))'; 
        % Convert (o,d) to OD flow ID
        ODFLOWIDTemp=zeros(1,1); 
        [~,ODFLOWIDTemp(1,1)]=ismember(ODTemp(1,:),PUMAinfo(:,2:3),'rows'); % bring flow ID of given (o,d)
        ODFLOWIDTemp(ODFLOWIDTemp(:,1)==0)=[];  % delete OD flow that does not exist
        ODChoice{1,i}=ODFLOWIDTemp';    % total OD covered by i-th choice
        for j=1:size(ODChoice{1,i},2)
            if ismember(ODChoice{1,i}(j),PrCorrMean(:,1)) % correlated
                ODChoiceIndCorr(PrCorrMean(:,1)==ODChoice{1,i}(j),i)=1;
            else
                ODChoiceIndUncorr(PrUncorr(:,1)==ODChoice{1,i}(j),i)=1;
            end
        end
        
        % AGGREGATE TRUTHS AND PRIORS ON FLOW LEVEL TO CHOICE LEVEL
        Mu(i,1)=TrUncorr(:,2)'*ODChoiceIndUncorr(:,i)+TrCorrMean(:,2)'*ODChoiceIndCorr(:,i);
        % There are many possible ways to assume initial prior of mean
        % flow not observed during the pilot. Here, the algorithm assumes
        % they are mean of all observed OD flows during pilots.
        AllMeans=[PrCorrMean(:,2);PrUncorr(:,2)];
        meanAll=mean(AllMeans(AllMeans>0));
        % Prior OD flow by choosing a link
        for j=1:nUncorrFlow
            if ODChoiceIndUncorr(j,i)>0 % j-th flow is covered by i-th choice
                if PrUncorr(j,2)==0 % if not observed
                    Mu0(i,1)=Mu0(i,1)+meanAll; % add mean value of the observed instead of zero
                else % if observed
                    Mu0(i,1)=Mu0(i,1)+PrUncorr(j,2); % add obseved mean flow
                end
            end
        end
        for j=1:nCorrFlow
            if ODChoiceIndCorr(j,i)>0 % j-th flow is covered by i-th choice
                if PrCorrMean(j,2)==0 % if not observed
                    Mu0(i,1)=Mu0(i,1)+meanAll; % add mean value of the observed instead of zero
                else % if observed
                    Mu0(i,1)=Mu0(i,1)+PrCorrMean(j,2); % add obseved mean flow
                end
            end
        end
        % True and prior std.dev by choosing a link
        Sigma(i,1)=(TrUncorr(:,3)').^2*ODChoiceIndUncorr(:,i)+diag(TrCorrCov)'*ODChoiceIndCorr(:,i);
        Sigma0(i,1)=(PrUncorr(:,3)').^2*ODChoiceIndUncorr(:,i)+diag(PrCorrCov)'*ODChoiceIndCorr(:,i);
    end

    % DIFFERENT LEARNING POLICIES FOR LINKE EXTENSIONS
    if mabsetting==1 % Multi-armed bandit policy (MAB)
        [Mu_est, OC, ObsOD, PrUncorr, PrCorrMean,PrecisU,PrecisC]=PolicyMAB(Mu,Mu0,M,ObsOD,ODChoice,ODSystemNet,PrUncorr,TrUncorr,PrCorrMean,TrCorrMean,TrCorrCov,PrecisU,PrecisC);
    elseif kggeneral==1 % Knowledge gradient policy (KG)
        [Mu_est, OC, ObsOD, PrUncorr, PrCorrMean,PrecisU,PrecisC]=PolicyKG(Mu,Mu0,Sigma0,M,ObsOD,ODChoice,ODSystemNet,PrUncorr,TrUncorr,PrCorrMean,TrCorrMean,TrCorrCov,PrecisU,PrecisC);
    elseif kgcb==1 % Knowledge gradient with correlated belief policy (KGCB)
        [Mu_est, OC, ObsOD, PrUncorr, PrCorrMean, PrCorrCov,PrecisU,PrecisC]=PolicyKGCB(Mu,Mu0,M,ObsOD,ODChoice,ODSystemNet,PrUncorr,TrUncorr,PrCorrMean,PrCorrCov,TrCorrMean,TrCorrCov,ODChoiceIndCorr,PrecisU,PrecisC);
    else
        [Mu_est, OC, ObsOD, PrUncorr, PrCorrMean,PrecisU,PrecisC]=PolicyGreedy(Mu,Mu0,M,ObsOD,ODChoice,ODSystemNet,PrUncorr,TrUncorr,PrCorrMean,TrCorrMean,TrCorrCov,PrecisU,PrecisC);
    end
    [~,I]=max(Mu_est);  % choice with the largest estimated flow
    [~,I2]=max(Mu);     % choice with the largest true flow
    nextLink=LinkChoice1(I);        % chosen link
    nextNode=setdiff(EdgeInfo(2:3,nextLink),NodeSystemNet);  % chosen node
    
    % for random results
    if policy<0
        I=randi(nCandLink);
        I2=0;
    end
    
    % UPDATE NODE ARCHIVE
    NodeRoute=[NodeRoute,nextNode];                 % attach new node
    Result(nSim+1).Route(RT).NodeRoute=NodeRoute;	% archive to Result
    NodeRouteList{NodeRoute(end),2}=RT;             % update incidence between node and route
    NodeSystemNet=[NodeSystemNet,nextNode];         % archive to NodeSystemNet
    
    % UPDATE OD PAIR ARCHIVE
    ODSystemNet=ODChoice{1,I};          % total covered OD
    ODSystemNetIndUncorr=ODChoiceIndUncorr(:,I);    % total covered OD (indexed)
    ODSystemNetIndCorr=ODChoiceIndCorr(:,I);    % total covered OD (indexed)
    
    % UPDATE OPPORTUNITY COST ARCHIVE
    OppCost{2,1}=RT;        % route ID
    OppCost{2,2}=nextNode;  % chosen node
    OppCost{2,3}=I;         % choice made
    OppCost{2,4}=I2;        % true best
    OppCost{2,5}=OC;        % opportunity cost
    
    % UPDATE LINK ARCHIVE
    IncludedLink=[IncludedLink,RemainingLink(:,RemainingLink(1,:)==nextLink)];	% add chosen link pair
    RemainingLink(:,RemainingLink(1,:)==nextLink)=[]; % delete chosen link pair
    IncludedLink=sortrows(IncludedLink')';  % sort included links              
    
    % VIRTUAL NETWORK SHOULD BE CREATED FOR TRAVEL TIME CALCULATION WITH TRANSFER LATER
    
    %% REPEAT LINK-LEVEL EXTENSION UNTIL CRITERIA ARE MET
    while length(NodeRoute)<maxL % # of nodes < maximum #
        % CANDIDATE LINKS FROM BOTH ROUTE ENDS
        LinkChoice2=[RemainingLink(1,RemainingLink(2,:)==NodeRoute(1)|RemainingLink(2,:)==NodeRoute(end)),RemainingLink(1,RemainingLink(3,:)==NodeRoute(1)|RemainingLink(3,:)==NodeRoute(end))]; % create next round candidates from remaining links
        nCandLink=size(LinkChoice2,2);  % # of options
        % DELETE LINKS ALREADY COVERED BY CURRENT ROUTE
        for i=1:nCandLink
            if sum(ismember(EdgeInfo(2:3,LinkChoice2(nCandLink-i+1))',NodeRoute))==2
                LinkChoice2(nCandLink-i+1)=[]; % inversely indexed to avoid error
            end
        end
        
        % EVALUATE LINKS BASED ON COLLECTED INFORMATION (ASSUMED) 
        % Only look into available links, not the entire network
        nCandLink=size(LinkChoice2,2);
        ODChoice=cell(1,nCandLink);
        ODChoiceIndUncorr=zeros(nUncorrFlow,nCandLink);
        ODChoiceIndCorr=zeros(nCorrFlow,nCandLink);
        for i=1:nCandLink
            ODChoice{1,i}=ODSystemNet;          % load previously included OD pairs
            ODChoiceIndUncorr(:,i)=ODSystemNetIndUncorr;    % load indicators of ODs that covered by choices
            ODChoiceIndCorr(:,i)=ODSystemNetIndCorr;
        end
        
        Mu=zeros(nCandLink,1);
        Mu0=zeros(nCandLink,1);
        Sigma=zeros(nCandLink,1);
        Sigma0=zeros(nCandLink,1);
        
        for i=1:nCandLink
            % IDENTIFY OD PAIRS NEWLY AVAILABLE PER CHOICE
            % Create from existing nodes
            ODTemp=NodeSystemNet';  % nodes included previously
            ODTemp(:,end+1)=setdiff(EdgeInfo(2:3,LinkChoice2(i)),NodeSystemNet); % pair with end node of new link (possible new node to the system)
            ODTemp=unique(sort(ODTemp,2,'ascend'),'rows');           % delete duplicates
            ODFLOWIDTemp=zeros(size(ODTemp,1),1);	% array for ODflowID
            for j=1:size(ODTemp,1)
                [~,ODFLOWIDTemp(j,1)]=ismember(ODTemp(j,:),PUMAinfo(:,2:3),'rows'); % find ODflowID
            end
            ODFLOWIDTemp(ODFLOWIDTemp(:,1)==0,:)=[]; % exclude OD pairs does not exist
            ODChoice{1,i}=[ODChoice{1,i},ODFLOWIDTemp']; % attach new OD pairs to ODChoice
            ODChoice{1,i}=unique(ODChoice{1,i}); % delete duplicates
            for j=1:size(ODChoice{1,i},2)
                if ismember(ODChoice{1,i}(j),PrCorrMean(:,1)) % correlated
                    ODChoiceIndCorr(PrCorrMean(:,1)==ODChoice{1,i}(j),i)=1;
                else
                    ODChoiceIndUncorr(PrUncorr(:,1)==ODChoice{1,i}(j),i)=1;
                end
            end

            % PREPARE FUNCTION INPUTS
            Mu(i,1)=TrUncorr(:,2)'*ODChoiceIndUncorr(:,i)+TrCorrMean(:,2)'*ODChoiceIndCorr(:,i);
            AllMeans=[PrCorrMean(:,2);PrUncorr(:,2)];
            meanAll=mean(AllMeans(AllMeans>0));
            for j=1:nUncorrFlow
                if ODChoiceIndUncorr(j,i)>0 % j-th flow is covered by i-th choice
                    if PrUncorr(j,2)==0 % if not observed
                        Mu0(i,1)=Mu0(i,1)+meanAll; % add mean value of the observed instead of zero
                    else % if observed
                        Mu0(i,1)=Mu0(i,1)+PrUncorr(j,2); % add obseved mean flow
                    end
                end
            end
            for j=1:nCorrFlow
                if ODChoiceIndCorr(j,i)>0 % j-th flow is covered by i-th choice
                    if PrCorrMean(j,2)==0 % if not observed
                        Mu0(i,1)=Mu0(i,1)+meanAll; % add mean value of the observed instead of zero
                    else % if observed
                        Mu0(i,1)=Mu0(i,1)+PrCorrMean(j,2); % add obseved mean flow
                    end
                end
            end
            Sigma(i,1)=(TrUncorr(:,3)').^2*ODChoiceIndUncorr(:,i)+diag(TrCorrCov)'*ODChoiceIndCorr(:,i);
            Sigma0(i,1)=(PrUncorr(:,3)').^2*ODChoiceIndUncorr(:,i)+diag(PrCorrCov)'*ODChoiceIndCorr(:,i);
        end
        
        % LEARNING PROCEDURE
        if mabsetting==1 % Multi-armed bandit
            [Mu_est2, OC2, ObsOD, PrUncorr, PrCorrMean,PrecisU,PrecisC]=PolicyMAB(Mu,Mu0,M,ObsOD,ODChoice,ODSystemNet,PrUncorr,TrUncorr,PrCorrMean,TrCorrMean,TrCorrCov,PrecisU,PrecisC);
        elseif kggeneral==1 % Knowledge gradient
            [Mu_est2, OC2, ObsOD, PrUncorr, PrCorrMean,PrecisU,PrecisC]=PolicyKG(Mu,Mu0,Sigma0,M,ObsOD,ODChoice,ODSystemNet,PrUncorr,TrUncorr,PrCorrMean,TrCorrMean,TrCorrCov,PrecisU,PrecisC);
        elseif kgcb==1 % Knowledge gradient with correlated belief
            [Mu_est2, OC2, ObsOD, PrUncorr, PrCorrMean, PrCorrCov,PrecisU,PrecisC]=PolicyKGCB(Mu,Mu0,M,ObsOD,ODChoice,ODSystemNet,PrUncorr,TrUncorr,PrCorrMean,PrCorrCov,TrCorrMean,TrCorrCov,ODChoiceIndCorr,PrecisU,PrecisC);
        else
            [Mu_est2, OC2, ObsOD, PrUncorr, PrCorrMean,PrecisU,PrecisC]=PolicyGreedy(Mu,Mu0,M,ObsOD,ODChoice,ODSystemNet,PrUncorr,TrUncorr,PrCorrMean,TrCorrMean,TrCorrCov,PrecisU,PrecisC);
        end
        
        [~,II]=max(Mu_est2);    % choose the link with the largest estimated flow
        [~,II2]=max(Mu);        % chosen link

        % for random results
        if policy<0
            II=randi(nCandLink);
            II2=0;
        end


        nextLink=LinkChoice2(II);
        nextNode=setdiff(EdgeInfo(2:3,nextLink),NodeSystemNet);   % chosen node
        
        % UPDATE NODE ARCHIVE
        ATT=EdgeInfo(2:3,nextLink)==NodeRoute(end);
        if sum(ATT)==1
            NodeRoute=[NodeRoute,nextNode];   % in a sequence of connection
        else
            NodeRoute=[nextNode,NodeRoute];  % extend in the opposite direction
        end
        NodeRouteList{nextNode,2}=RT;   % indicate current
        Result(nSim+1).Route(RT).NodeRoute=NodeRoute; % update NodeRoute
        NodeSystemNet=NodeRoute; % update NodeSystemNet (equivalent only for 1st route)
        
        % UPDATE OPPORTUNITY COST ARCHIVE
        OppCost{length(NodeRoute),1}=RT;
        OppCost{length(NodeRoute),2}=nextNode;
        OppCost{length(NodeRoute),3}=II;
        OppCost{length(NodeRoute),4}=II2;
        OppCost{length(NodeRoute),5}=OC2;
        
        % UPDATE LINK ARCHIVE
        IncludedLink=[IncludedLink,RemainingLink(:,RemainingLink(1,:)==nextLink)];	% add chosen link pair
        RemainingLink(:,RemainingLink(1,:)==nextLink)=[]; % delete chosen link pair
        IncludedLink=sortrows(IncludedLink')';  % sort included links
        
        % UPDATE OD PAIR ARCHIVE
        ODSystemNet=ODChoice{1,II};
        ODSystemNetIndUncorr=ODChoiceIndUncorr(:,II);    % total covered OD (indexed)
        ODSystemNetIndCorr=ODChoiceIndCorr(:,II);
    end
    
    EstMeanODFlow=[EstMeanODFlow,[PrUncorr(:,1:2);PrCorrMean(:,1:2)]];
    Result(nSim+1).ServingOD{RT,1}=ODSystemNet;
    OPPCOST=[OPPCOST;OppCost];
    
    %% PREPARE TRANSFER INFORMATION (NODE COVERAGE)
    % 1. Relationship between routes: transferable?
    % 2. Transfer node: intersections
    % 3. Available OD pairs: either choice- or route-level
    
    % COVERED OD PAIRS PER ROUTE PAIR
    ODwithRTpair=cell(nReqRoute); % each cell will have OD pairs covered by i-/j-th routes
    ODwithRTpair{1,1}=Result(nSim+1).ServingOD{RT,1}; % result of 1st route    
    
    %% SUBSEQUENT ROUTES (RT>1)
    for RT = 2:nReqRoute
        % NODES COVERED BY CURRENT NETWORK
        % When RT>1: nodes and links of previous routes
        
        %{
        FOR EVERY 1ST EXTENSION, THE ALGORITHM SHOULD INVESTIGATE:
        1. IS A NODE COVERED BY ANY OTHER ROUTES?
           -> YES: 2, NO: CASE I
        2. DOES A ROUTE ALREADY INTERSECT WITH CURRENT ONE?
           -> YES: CASE II, NO: CASE III
        
        CASE I: GENERAL NODE
                AVAILABLE NEW OD PAIRS: BETWEEN NODE AND CONNECTABLE ROUTES
        CASE II: COVERED AND ACCESSIBLE NODE
                AVAILABLE NEW OD PAIRS: NONE
        CASE III: COVERED BUT INACCESSIBLE NODE
                AVAILABLE NEW OD PAIRS: BETWEEN CURRENT ROUTE AND NEW ONE
                                        (BETWEEN ALL NODES OF BOTH)
        %}
        
        OppCost=cell(maxL,5);
        OppCost{1,1}=RT;
        % LARGEST IN- AND OUTFLOW TO ALL OTHER NODES & SUBTRACT COVERED : REMAINING OD FLOWS
        % Identify potential links to extend
        CandStartNode = NodeSystemNet;      % all existing nodes can be candidates
        nCandStart = size(CandStartNode,2); % # of candidates 
        
        % FIND THE EXISTING NODE WITH THE LARGEST REMAINING FLOW
        NNF1=zeros(1,nCandStart); % net flow per candidate node
        for i=1:nCandStart
            PairsU=setdiff(PairperNode{CandStartNode(i),1},ODSystemNet); % uncorrelated pairs
            PairsC=setdiff(PairperNode{CandStartNode(i),2},ODSystemNet); % correlated pairs
            for j=1:length(PairsU)
                NNF1(i)=NNF1(i)+PrUncorr(PrUncorr(:,1)==PairsU(j),2);
            end
            for j=1:length(PairsC)
                NNF1(i)=NNF1(i)+PrCorrMean(PrCorrMean(:,1)==PairsC(j),2);
            end
        end
        
        NNF2=zeros(1,4);
        k=0;
        for i=1:nCandStart
            RN=ReachableNodes{CandStartNode(i),2};
            for j=1:length(RN)
                % Exclude 2nd nodes without vacant adjacent nodes
                if sum(setdiff(ReachableNodes{RN(j),2},NodeSystemNet))>0
                    k=k+1;
                    NNF2(k,1:3)=[CandStartNode(i),RN(j),NNF1(i)];
                    PairsU=setdiff(PairperNode{RN(j),1},ODSystemNet); % uncorrelated pairs
                    PairsC=setdiff(PairperNode{RN(j),2},ODSystemNet); % correlated pairs
                    for m=1:length(PairsU)
                        NNF2(k,4)=NNF2(k,4)+PrUncorr(PrUncorr(:,1)==PairsU(m),2);
                    end
                    for m=1:length(PairsC)
                        NNF2(k,4)=NNF2(k,4)+PrCorrMean(PrCorrMean(:,1)==PairsC(m),2);
                    end
                    edgeDup=intersect(PairperNode{RN(j),1},PairperNode{CandStartNode(i),1});
                    if ~ismember(edgeDup,ODSystemNet)
                        if edgeDup>0
                            NNF2(k,4)=NNF2(k,4)-PrUncorr(PrUncorr(:,1)==edgeDup,2);
                        else
                            edgeDup=intersect(PairperNode{RN(j),2},PairperNode{CandStartNode(i),2});
                            NNF2(k,4)=NNF2(k,4)-PrCorrMean(PrCorrMean(:,1)==edgeDup,2);
                        end
                    end
                end
            end
        end
        
        NNF2(:,5)=sum(NNF2(:,3:4),2);
        
        NodeNetFlow=NNF2(:,5)';
        [~,maxSet]=max(NodeNetFlow);
        if policy<0
            maxSet=randi(length(NodeNetFlow));
        end
        StartNodeSet=NNF2(maxSet,1:2);
        start=StartNodeSet(1);
        nextNode=StartNodeSet(2);
        SNS=sort(StartNodeSet);
        [~,nextLink]=ismember(SNS,EdgeInfo(2:3,:)','rows');
        
        % CHOSEN NODE: "start"
        NodeRoute=start;	% add this node to route
        NodeSystemNet=unique([NodeSystemNet,start]); % update total nodes in the system
        OppCost{1,2}=start; % 1st choice is "start"
        Result(nSim+1).Route(RT).NodeRoute=start;   % add "start" to Rt-th route as 1st node
        NodeRouteList{start,2}=unique([NodeRouteList{start,2},RT]); % update routes visiting "start"
        RTStartA=NodeRouteList{start,2};    % update "RTStartA"
        nRT_A=length(RTStartA); % # of routes serving "start"
        for i=1:nRT_A
            RTintsct(RTStartA(i),RT)=1; % update RTintsct since routes in RTStartA intersect with RT
        end
        
        % WE NEED TO IDENTIFY WHICH OD PAIRS ARE ACTIVATED WHEN CHOOSING A
        % LINK AMONG AVAILABLE ONES
        LinkChoice1=EdgeInfo(:,nextLink);
        nCandLink=size(LinkChoice1,2);      % # of available links
        ODChoice=cell(1,nCandLink);    % [total OD covered]
        ODChoiceIndUncorr=zeros(nUncorrFlow,nCandLink);
        ODChoiceIndCorr=zeros(nCorrFlow,nCandLink);
        for i=1:nCandLink
            ODChoice{1,i}=ODSystemNet;          % load previously included OD pairs
            ODChoiceIndUncorr(:,i)=ODSystemNetIndUncorr;    % load indicators of ODs that covered by choices
            ODChoiceIndCorr(:,i)=ODSystemNetIndCorr;
        end

        Mu=zeros(nCandLink,1);      % prepare mean flow based on truth
        Mu0=zeros(nCandLink,1);     % prepare mean flow based on prior
        Sigma=zeros(nCandLink,1);	% prepare std.dev of flow based on truth
        Sigma0=zeros(nCandLink,1);	% prepare std.dev of flow based on prior
        
        % IDENTIFY OD FLOWS INVOLVED IN OPTIONS
        for i=1:nCandLink
            ODTemp=[];	% temporary space for total new OD flows from choosing i
            ODTemp2=[]; % from routes serving "start" to opposite end of link
            ODTemp3=[]; % from "start" to routes serving opposite end of link
            RTStartB=NodeRouteList{LinkChoice1(3,i),2}; % routes covering opposite end of link
            nRT_B=length(RTStartB); % # of routes covering opposite end of link
            
            % CLASSIFICATION OF LINK CASE
            % Case I: new OD pairs are generated between new node and 
            %         existing transferable routes (RTintsct(:,RT)==1)
            % Case II: infeasible
            % Case III: new OD pairs are generated between
            %           1) A and route covering B and not accessible from A
            %           2) B and route covering A and not accessible from B
            
            % Regardless of link case, routes serving "start" and opposite end must be connected
            for j=1:nRT_A
                % NEW OD PAIRS: BETWEEN NEW NODE(B) AND NODES IN "RTStartA"
                ODTemp2=Result(nSim+1).Route(RTStartA(j)).NodeRoute'; % nodes included in j-th route
                snode=setdiff(LinkChoice1(2:3,i)',ODTemp2);
                if snode>0
                    ODTemp2(:,2)=snode;
                    ODTemp=[ODTemp;ODTemp2];
                end
            end
            
            % If i-th link is Case III, it should consider connections between "start" and routes serving opposite end
            if nRT_B>0 % there exist routes covering opposite end
                for j=1:nRT_B
                    % NEW OD PAIRS: BETWEEN "start" AND NODES IN "RTStartB"
                    ODTemp3=Result(nSim+1).Route(RTStartB(j)).NodeRoute'; % nodes included in j-th route
                    snode=setdiff(LinkChoice1(2:3,i)',ODTemp3);
                    if snode>0
                        ODTemp3(:,2)=snode;
                        ODTemp=[ODTemp;ODTemp3];
                    end
                end
            end
            ODTemp=unique(ODTemp,'rows');           % All unique OD pairs (i,j) of choice
            ODTemp=sort(ODTemp,2);
            ODFLOWIDTemp=zeros(size(ODTemp,1),1);   % Convert (i,j) to ODflowID
            for j=1:size(ODTemp,1)
                [~,ODFLOWIDTemp(j,1)]=ismember(ODTemp(j,:),PUMAinfo(:,2:3),'rows');
            end
            ODFLOWIDTemp(ODFLOWIDTemp(:,1)==0,:)=[];    % delete empty rows
            ODChoice{1,i}=ODFLOWIDTemp';   % ODs covered by choice
            ODChoice{1,i}=setdiff(ODChoice{1,i},ODSystemNet);
            
            for j=1:size(ODChoice{1,i},2)
                if ismember(ODChoice{1,i}(j),PrCorrMean(:,1)) % correlated
                    ODChoiceIndCorr(PrCorrMean(:,1)==ODChoice{1,i}(j),i)=1;
                else
                    ODChoiceIndUncorr(PrUncorr(:,1)==ODChoice{1,i}(j),i)=1;
                end
            end
            % PREPARE FUNCTION INPUTS
            Mu(i,1)=TrUncorr(:,2)'*ODChoiceIndUncorr(:,i)+TrCorrMean(:,2)'*ODChoiceIndCorr(:,i);
            AllMeans=[PrCorrMean(:,2);PrUncorr(:,2)];
            meanAll=mean(AllMeans(AllMeans>0));
            for j=1:nUncorrFlow
                if ODChoiceIndUncorr(j,i)>0 % j-th flow is covered by i-th choice
                    if PrUncorr(j,2)==0 % if not observed
                        Mu0(i,1)=Mu0(i,1)+meanAll; % add mean value of the observed instead of zero
                    else % if observed
                        Mu0(i,1)=Mu0(i,1)+PrUncorr(j,2); % add obseved mean flow
                    end
                end
            end
            for j=1:nCorrFlow
                if ODChoiceIndCorr(j,i)>0 % j-th flow is covered by i-th choice
                    if PrCorrMean(j,2)==0 % if not observed
                        Mu0(i,1)=Mu0(i,1)+meanAll; % add mean value of the observed instead of zero
                    else % if observed
                        Mu0(i,1)=Mu0(i,1)+PrCorrMean(j,2); % add obseved mean flow
                    end
                end
            end
            Sigma(i,1)=(TrUncorr(:,3)').^2*ODChoiceIndUncorr(:,i)+diag(TrCorrCov)'*ODChoiceIndCorr(:,i);
            Sigma0(i,1)=(PrUncorr(:,3)').^2*ODChoiceIndUncorr(:,i)+diag(PrCorrCov)'*ODChoiceIndCorr(:,i);
        end
        
        % NEW LINE
        I=1;

        IncludedLink=[];
       
        NodeRouteList{nextNode,2}=RT;   % indicate current
        
        NodeSystemNet=unique([NodeSystemNet,nextNode]);
        NodeRoute=[NodeRoute,nextNode];   % in a sequence of connection
        Result(nSim+1).Route(RT).NodeRoute=NodeRoute;
        NodeRouteList{NodeRoute(end),2}=[NodeRouteList{NodeRoute(end),2},RT];

        ODSystemNet = unique([ODSystemNet,ODChoice{I}]);
        ODSystemNetIndUncorr=ODChoiceIndUncorr(:,I);    % total covered OD (indexed)
        ODSystemNetIndCorr=ODChoiceIndCorr(:,I);    % total covered OD (indexed)
        OppCost{2,2}=nextNode;
        
        for i=1:RT
            ODTemp=[];
            if RTintsct(i,RT)==1
                ODTemp=Result(nSim+1).Route(i).NodeRoute';
                ODTemp(:,2)= nextNode;
            end
            ODTemp=unique(ODTemp,'rows');           % All unique OD pairs (i,j) of choice
            ODTemp=sort(ODTemp,2);
            ODFLOWIDTemp=zeros(size(ODTemp,1),1);   % Convert (i,j) to ODflowID
            for j=1:size(ODTemp,1)
                [~,ODFLOWIDTemp(j,1)]=ismember(ODTemp(j,:),PUMAinfo(:,2:3),'rows');
            end

            ODFLOWIDTemp(ODFLOWIDTemp(:,1)==0,:)=[];
            ODwithRTpair{i,RT}=[ODwithRTpair{i,RT},ODFLOWIDTemp'];
        end 
        
        % UPDATE LINK ARCHIVE
        IncludedLink=[IncludedLink,RemainingLink(:,RemainingLink(1,:)==nextLink)];	% add chosen link pair
        RemainingLink(:,RemainingLink(1,:)==nextLink)=[]; % delete chosen link pair
        IncludedLink=sortrows(IncludedLink')';  % sort included links
        
        for i=1:RT
            K=[];
            for j=1:length(Result(nSim+1).Route(i).NodeRoute)
                K=[K,NodeRouteList{Result(nSim+1).Route(i).NodeRoute(j),2}];
            end
            RTGroups{i,2}=unique(K);            
        end
        
        %% 5. IDENTIFY ADJACENT LINKS TO EXTEND AND REPEAT

        while length(NodeRoute)<maxL
            LinkChoice2=[RemainingLink(:,RemainingLink(2,:)==NodeRoute(1)|RemainingLink(2,:)==NodeRoute(end)),RemainingLink(:,RemainingLink(3,:)==NodeRoute(1)|RemainingLink(3,:)==NodeRoute(end))]; % create next round candidates
            LinkChoice22=LinkChoice2;
            nCandLink=size(LinkChoice2,2);
            for i=1:nCandLink % exclude nodes already on this route
                if sum(ismember(EdgeInfo(2:3,LinkChoice2(1,nCandLink-i+1))',NodeRoute))==2
                    LinkChoice2(:,nCandLink-i+1)=[]; % inversely indexed to avoid error
                end
            end
            nCandLink=size(LinkChoice2,2);
            % if all nodes are covered, randomly pick one route to extend
            % among adjacent nodes
 
            if nCandLink==0 % all links are already used: randomly choose one
                LinkChoice2=LinkChoice22; % create next round candidates
                nCandLink=size(LinkChoice2,2);
                LinkChoice2=LinkChoice2(:,randi(nCandLink));
                nCandLink=1;
            end
            
            ODChoice=cell(1,nCandLink);    % [total OD covered]
            ODChoiceIndUncorr=zeros(nUncorrFlow,nCandLink);
            ODChoiceIndCorr=zeros(nCorrFlow,nCandLink);
            for i=1:nCandLink
                ODChoiceIndUncorr(:,i)=ODSystemNetIndUncorr;    % load indicators of ODs that covered by choices
                ODChoiceIndCorr(:,i)=ODSystemNetIndCorr;
            end
            Mu=zeros(nCandLink,1);      % prepare mean flow based on truth
            Mu0=zeros(nCandLink,1);     % prepare mean flow based on prior
            Sigma=zeros(nCandLink,1);	% prepare std.dev of flow based on truth
            Sigma0=zeros(nCandLink,1);	% prepare std.dev of flow based on prior

            for i=1:nCandLink
                ODTemp=[]; ODTemp2=[]; ODTemp3=[];
                RTStartA=RTGroups{RT,2};
                nRT_A=length(RTStartA);
                t=setdiff(LinkChoice2(2:3,i)',[NodeRoute(1),NodeRoute(end)]);
                RTStartB=NodeRouteList{t,2}; % routes covering opposite end of link
                nRT_B=length(RTStartB); % # of routes covering opposite end of link
                % CLASSIFICATION OF LINK CASE
                for j=1:nRT_A
                    ODTemp2=Result(nSim+1).Route(RTStartA(j)).NodeRoute'; % nodes included in j-th route
                    ODTemp2(:,2)=t;
                    ODTemp=[ODTemp;ODTemp2];
                end


                if nRT_B>0 % some routes covering opposite end
                    % Case III: new OD pairs are generated between
                    % 1) A and route covering B and not accessible from A
                    % 2) B and route covering A and not accessible from B
                    for j=1:nRT_B
                        % NEW OD PAIRS
                        % : BETWEEN A AND RTStartB
                        % : BETWEEN B AND RTStartA
                        for k=1:length(NodeRoute)
                            ODTemp3=Result(nSim+1).Route(RTStartB(j)).NodeRoute'; % nodes included in j-th route
                            ODTemp3(:,2)=NodeRoute(k);
                            ODTemp=[ODTemp;ODTemp3];
                        end
                    end
                end
                
                ODTemp=unique(ODTemp,'rows');           % All unique OD pairs (i,j) of choice
                ODTemp=sort(ODTemp,2);
                ODFLOWIDTemp=zeros(size(ODTemp,1),1);   % Convert (i,j) to ODflowID
                for j=1:size(ODTemp,1)
                    [~,ODFLOWIDTemp(j,1)]=ismember(ODTemp(j,:),PUMAinfo(:,2:3),'rows');
                end
 
                ODFLOWIDTemp(ODFLOWIDTemp(:,1)==0,:)=[];
                ODChoice{1,i}=ODFLOWIDTemp';   % ODs covered by choice
                ODChoice{1,i}=setdiff(ODChoice{1,i},ODSystemNet);
                
                for j=1:size(ODChoice{1,i},2)
                    if ismember(ODChoice{1,i}(j),PrCorrMean(:,1)) % correlated
                        ODChoiceIndCorr(PrCorrMean(:,1)==ODChoice{1,i}(j),i)=1;
                    else
                        ODChoiceIndUncorr(PrUncorr(:,1)==ODChoice{1,i}(j),i)=1;
                    end
                end
                % PREPARE FUNCTION INPUTS
                Mu(i,1)=TrUncorr(:,2)'*ODChoiceIndUncorr(:,i)+TrCorrMean(:,2)'*ODChoiceIndCorr(:,i);
                AllMeans=[PrCorrMean(:,2);PrUncorr(:,2)];
                meanAll=mean(AllMeans(AllMeans>0));
                for j=1:nUncorrFlow
                    if ODChoiceIndUncorr(j,i)>0 % j-th flow is covered by i-th choice
                        if PrUncorr(j,2)==0 % if not observed
                            Mu0(i,1)=Mu0(i,1)+meanAll; % add mean value of the observed instead of zero
                        else % if observed
                            Mu0(i,1)=Mu0(i,1)+PrUncorr(j,2); % add obseved mean flow
                        end
                    end
                end
                for j=1:nCorrFlow
                    if ODChoiceIndCorr(j,i)>0 % j-th flow is covered by i-th choice
                        if PrCorrMean(j,2)==0 % if not observed
                            Mu0(i,1)=Mu0(i,1)+meanAll; % add mean value of the observed instead of zero
                        else % if observed
                            Mu0(i,1)=Mu0(i,1)+PrCorrMean(j,2); % add obseved mean flow
                        end
                    end
                end
                Sigma(i,1)=(TrUncorr(:,3)').^2*ODChoiceIndUncorr(:,i)+diag(TrCorrCov)'*ODChoiceIndCorr(:,i);
                Sigma0(i,1)=(PrUncorr(:,3)').^2*ODChoiceIndUncorr(:,i)+diag(PrCorrCov)'*ODChoiceIndCorr(:,i);
            end
            
            if sum(Mu)>0
                for i=1:nCandLink
                    if Mu(nCandLink-i+1,1)==0
                        LinkChoice2(:,nCandLink-i+1)=[];
                        Mu(nCandLink-i+1,:)=[];
                        Mu0(nCandLink-i+1,:)=[];
                        ODChoice(:,nCandLink-i+1)=[];
                    end
                end
                if mabsetting==1 % Multi-armed bandit
                    [Mu_est, OC2, ObsOD, PrUncorr, PrCorrMean,PrecisU,PrecisC]=PolicyMAB(Mu,Mu0,M,ObsOD,ODChoice,ODSystemNet,PrUncorr,TrUncorr,PrCorrMean,TrCorrMean,TrCorrCov,PrecisU,PrecisC);
                elseif kggeneral==1 % Knowledge gradient
                    [Mu_est, OC2, ObsOD, PrUncorr, PrCorrMean,PrecisU,PrecisC]=PolicyKG(Mu,Mu0,Sigma0,M,ObsOD,ODChoice,ODSystemNet,PrUncorr,TrUncorr,PrCorrMean,TrCorrMean,TrCorrCov,PrecisU,PrecisC);
                elseif kgcb==1 % Knowledge gradient with correlated belief
                    [Mu_est, OC2, ObsOD, PrUncorr, PrCorrMean, PrCorrCov,PrecisU,PrecisC]=PolicyKGCB(Mu,Mu0,M,ObsOD,ODChoice,ODSystemNet,PrUncorr,TrUncorr,PrCorrMean,PrCorrCov,TrCorrMean,TrCorrCov,ODChoiceIndCorr,PrecisU,PrecisC);
                else
                    [Mu_est, OC2, ObsOD, PrUncorr, PrCorrMean,PrecisU,PrecisC]=PolicyGreedy(Mu,Mu0,M,ObsOD,ODChoice,ODSystemNet,PrUncorr,TrUncorr,PrCorrMean,TrCorrMean,TrCorrCov,PrecisU,PrecisC);
                end
                [~,II]=max(Mu_est);
                [~,II2]=max(Mu);
            else
                II=randi(nCandLink);
                II2=0;
            end

            % for random results
            if policy<0
                II=randi(nCandLink);
                II2=0;
            end

            nextLink=LinkChoice2(1,II);
            nextNode=setdiff(EdgeInfo(2:3,nextLink),NodeRoute);   % chosen node
            NodeSystemNet=unique([NodeSystemNet,nextNode]);            
            
            ATT=EdgeInfo(2:3,nextLink)==NodeRoute(end);
            if sum(ATT)==1
                NodeRoute=[NodeRoute,nextNode];   % in a sequence of connection
            else
                NodeRoute=[nextNode,NodeRoute];  % extend in the opposite direction
            end
            NodeRouteList{nextNode,2}=[NodeRouteList{nextNode,2},RT];
            
            OppCost{length(NodeRoute),1}=RT;
            OppCost{length(NodeRoute),2}=nextNode;
            OppCost{length(NodeRoute),3}=II;
            OppCost{length(NodeRoute),4}=II2;
            OppCost{length(NodeRoute),5}=OC2;
            
            RTStartB=NodeRouteList{nextNode,2}; % routes covering opposite end of link
            nRT_B=length(RTStartB);
            for i=1:nRT_B
                RTintsct(RTStartB(i),RT)=1;
            end
            
            IncludedLink=[IncludedLink,RemainingLink(:,RemainingLink(1,:)==nextLink)];	% add chosen link pair
            RemainingLink(:,RemainingLink(1,:)==nextLink)=[]; % delete chosen link pair
            IncludedLink=sortrows(IncludedLink')';  % sort included links
            
            ODSystemNet=unique([ODChoice{1,II},ODSystemNet]);
            ODSystemNetIndUncorr=ODChoiceIndUncorr(:,II);    % total covered OD (indexed)
            ODSystemNetIndCorr=ODChoiceIndCorr(:,II);
            Result(nSim+1).Route(RT).NodeRoute=NodeRoute;
            for i=1:RT
                K=[];
                for j=1:length(Result(nSim+1).Route(i).NodeRoute)
                    K=[K,NodeRouteList{Result(nSim+1).Route(i).NodeRoute(j),2}];
                end
                RTGroups{i,2}=unique(K);            
            end
            
            for i=1:RT
                ODTemp=[];
                if RTintsct(i,RT)==1
                    ODTemp=Result(nSim+1).Route(i).NodeRoute';
                    ODTemp(:,2)=nextNode;
                 end
                ODTemp=unique(ODTemp,'rows');           % All unique OD pairs (i,j) of choice
                ODFLOWIDTemp=zeros(size(ODTemp,1),1);   % Convert (i,j) to ODflowID
                for j=1:size(ODTemp,1)
                    [~,ODFLOWIDTemp(j,1)]=ismember(ODTemp(j,:),PUMAinfo(:,2:3),'rows');
                end

                ODFLOWIDTemp(ODFLOWIDTemp(:,1)==0,:)=[];
                ODwithRTpair{i,RT}=[ODwithRTpair{i,RT},ODFLOWIDTemp'];
            end 
        end

        EstMeanODFlow=[EstMeanODFlow,[PrUncorr(:,1:2);PrCorrMean(:,1:2)]];
        Result(nSim+1).ServingOD{RT,1}=ODSystemNet;
        OPPCOST=[OPPCOST;OppCost];
        Result(nSim+1).Route(RT).NodeRoute=NodeRoute;
    
    end
    TrODMean=[TrUncorr(:,1:2);TrCorrMean(:,1:2)];
    nSim=nSim+1
    Result(nSim).EstMeanODFlow=EstMeanODFlow;
    Result(nSim).OPPCOST=OPPCOST;
    Result(nSim).Err(1,:)=Result(nSim).Err(1,:)+1;
    T(nSim,1)=toc;
    catch
        Result(nSim+1).Err(1,2)=Result(nSim+1).Err(1,2)+1;
    end
end

%% 6. YIELD RESULTS
% COVERED TRUE DEMAND
TotFlow=zeros(Run,2*RT);   % column1: est, column2: true
for i=1:Run
    for j=1:RT
        for k=1:size(Result(i).ServingOD{j,1},2)
            TotFlow(i,j)=TotFlow(i,j)+Result(i).EstMeanODFlow(Result(i).EstMeanODFlow(:,1)==Result(i).ServingOD{j,1}(k),end);
            TotFlow(i,j+RT)=TotFlow(i,j+RT)+TrODMean(TrODMean(:,1)==Result(i).ServingOD{j,1}(k),2);
        end
    end
end

% OPTIMAL CHOICE RATE
OptChc=zeros(Run,1);
for j=1:Run
    Chc=zeros(1,6); %[Total iteration, # of optimal choice, # of optimal+random choice]
    for k=1:nReqRoute*maxL
        if Result(j).OPPCOST{k,3}>0
            Chc(1,1)=Chc(1,1)+1;
            if Result(j).OPPCOST{k,3}==Result(j).OPPCOST{k,4}
                Chc(1,2:3)=Chc(1,2:3)+1;
            elseif Result(j).OPPCOST{k,4}==0
                Chc(1,3)=Chc(1,3)+1;
            end
        end
    end
    OptChc(j,1)=Chc(1,2)/Chc(1,1);
    OptChc(j,2)=Chc(1,3)/Chc(1,1);
end

% ERROR RATE
ErrRate=zeros(1,2);
for j=1:Run
    ErrRate=ErrRate+Result(j).Err;
end
Result = rmfield(Result,'Err');
end