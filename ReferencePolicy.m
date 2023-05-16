function [EstDemand]=ReferencePolicy(gen,InputData)
%{
Motivated by:
Chow, J. Y., & Sayarshad, H. R. (2016). Reference policies for non-myopic 
sequential network design and timing problems. Networks and Spatial 
Economics, 16(4), 1183-1209. doi: 10.1007/s11067-015-9315-5.

[INPUT]
gen: # of draws of estimated demand from reference policy
Mat: 
[OUTPUT]
EstDemand: gen-by-1 vector with estimated demand

[Pseudocode for reference policy]
1. Recall network topology
    (not used) a. Identify node-incidence matrix
    b. Create a table relating nodes and their adjacent ones
2. Generate random route sets
    a. Randomly choose a series of node
    b. Start the next route from existing nodes
        i. Choosing nodes already included in the system is allowed
        ii. Excessively verlapping the existing route is prohibited.
    c. Repeat until the number of routes reaches the limit.
3. Estimate the demand coverage of chosen route sets

4. Estimate the maximum value by assuming extreme values follows either Weibull or Gumbel distribution

%}
    
    %% 1.RECALL NETWORK TOPOLOGY
%     load(Mat, 'ReachableNodes','PUMAinfo','TrueCorrPairMeanH','TrueUncorrPairH'); % bring adjacent node information
    EstDemand = zeros(gen,3); % prepare vector for estimated demand from different variability scenarios
    load(InputData, 'ReachableNodes','PUMAinfo'); % bring adjacent node information
    for vrblty = 1:3
        if vrblty == 1 % low variability
            load(InputData, 'TrueCorrPairMeanL','TrueUncorrPairL');
            CorrPairMean=TrueCorrPairMeanL;
            UncorrPairMean=TrueUncorrPairL;
        elseif vrblty == 2 % medium variability
            load(InputData, 'TrueCorrPairMeanM','TrueUncorrPairM');
            CorrPairMean=TrueCorrPairMeanM;
            UncorrPairMean=TrueUncorrPairM;
        else % vrblty == 3: high variability
            load(InputData, 'TrueCorrPairMeanH','TrueUncorrPairH');
            CorrPairMean=TrueCorrPairMeanH;
            UncorrPairMean=TrueUncorrPairH;
        end
   
        NodeAdj=ReachableNodes(:,1:2);  % select 1st and 2nd column only
        nNode=size(NodeAdj,1);  % # of nodes
        
        %% 2. GENERATE RANDOM ROUTE SETS
        rtlength=8; % route length limit in # of nodes
        nRoutes=5;
        Route=cell(1,nRoutes);
        p=1;
        
        while p<=gen
            % Initial route
            initnode=randi(nNode);  % randomly choose one of nodes
            Route{p,1}=initnode;
            SYSNODE=initnode;
            k=1;
            while k<rtlength
                t=rand;
                if t>0.5
                    i=Route{p,1}(k);
                else
                    i=Route{p,1}(1);
                end
                nOption=length(NodeAdj{i,2});
                j=NodeAdj{i,2}(randi(nOption));
                if sum(Route{p,1}==j)==0
                    if t>0.5
                        Route{p,1}=[Route{p,1},j];
                    else
                        Route{p,1}=[j,Route{p,1}];
                    end
                    SYSNODE=[SYSNODE,j];
                    k=k+1;
                else
                    if sum(setdiff(NodeAdj{i,2},Route{p,1}))==0
                        k=1;
                        break
                    end
                end
            end
            % Subsequent routes
            m=1;    % # of complete route
            %{
            We should determine the level of overlapping between routes
            i) Only one node: may be too dispersed
            -> Current route should dodge the existing route system
            ii) Limited to a certain number: hard to define (CHOSEN)
            -> Current route can share some nodes with the existing route system
            iii) Don't care: cannot avoid redundancy
            -> Current route randomly expands regardless of the existing route system
            %}
            maxOvrlp=4;	% # of allowed overlapping
        
            while m<=nRoutes
                SYSNODE=[];
                for a=1:m
                    SYSNODE=[SYSNODE,Route{p,a}];
                end
                SYSNODE=unique(SYSNODE);
                nAvailInit=length(SYSNODE);	% # of available nodes to start the next route
                % Initial node and link
                i=SYSNODE(randi(nAvailInit));   % randomly choose one of them
                nOvrlp=1;   % initial node is overlapped with the existing system
                nOption=length(NodeAdj{i,2});   % # of available neighbor from i
                j=NodeAdj{i,2}(randi(nOption)); % choose the next node to connect
                if any(SYSNODE==j)
                    nOvrlp=2;   % both i and j are already covered by the existing system
                end
                Route{p,m}=[i,j];   % add i and j to m-th route
                k=2;    % two nodes are determined
                % Subsequent nodes
                while k<rtlength
                    t=rand;
                    if t>0.5
                        i=Route{p,m}(k); % bring the current last node of m-th route
                    else
                        i=Route{p,m}(1);
                    end
                    nOption=length(NodeAdj{i,2});   % # of available nodes from i
                    j=NodeAdj{i,2}(randi(nOption)); % choose one of available nodes
                    if any(SYSNODE==j)
                        nOvrlp=nOvrlp+1;   % both i and j are already covered by the existing system
                    end
        
                    % CONDITIONS FOR RESTART THE M-TH ROUTE BUILDING
                    if nOvrlp>maxOvrlp
                        k=2;
                        break % stop and restart the current route
                    end
        
                    if sum(Route{p,m}==j)==0
                        if t>0.5
                            Route{p,m}=[Route{p,m},j];
                        else
                            Route{p,m}=[j,Route{p,m}];
                        end
                        SYSNODE=unique([SYSNODE,j]);
                        k=k+1;
                    else
                        if sum(setdiff(NodeAdj{i,2},Route{p,m}))==0
                            k=2;
                            break
                        end
                    end
                end
                if k==rtlength
                    m=m+1;
                end
            end
            p=p+1;
        end
        
        %% 3. ESTIMATE THE MAXIMUM VALUE
        % Estimate the demand coverage of chosen route sets
        %{
        Identify transferability among routes
        Extract available OD pairs and aggregate their flows
        Estimate the maximum value from the distribution of aggregated flow
        %}
        
        for r=1:gen
            % Create route incidence matrix and array
            AdjRoutesInd=eye(nRoutes);    % 1: intersect, 0: not intersect
            TransferNodes=cell(nRoutes);    % transfer nodes
            for i=1:nRoutes-1
                for j=i+1:nRoutes
                    A=intersect(Route{r,i},Route{r,j});
                    if isempty(A)==0
                        AdjRoutesInd(i,j)=1;
                        AdjRoutesInd(j,i)=1;
                        TransferNodes{i,j}=A;
                        TransferNodes{j,i}=A;
                    end
                end
            end
        
            % Find OD pairs available with different route pairs
            ODbyRoute=cell(nRoutes);
            % OD pairs that can be covered within a route
            for i=1:nRoutes
                for j=1:rtlength-1
                    for k=j+1:rtlength
                        if Route{r,i}(j)<Route{r,i}(k)
                            [~,u1]=ismember([Route{r,i}(j),Route{r,i}(k)],PUMAinfo(:,2:3),'rows');
                        else
                            [~,u1]=ismember([Route{r,i}(k),Route{r,i}(j)],PUMAinfo(:,2:3),'rows');
                        end
                        ODbyRoute{i,i}=[ODbyRoute{i,i},u1];
                    end
                end
            end
        
            % OD pairs that can be covered with 2 different routes
            for i=1:nRoutes-1
                for j=i+1:nRoutes
                    if AdjRoutesInd(i,j)==1
                        Ni=setdiff(Route{r,i},TransferNodes{i,j});
                        Nj=setdiff(Route{r,j},TransferNodes{i,j});
                        for k=1:length(Ni)
                            for m=1:length(Nj)
                                if Ni(k)<Nj(m)
                                    [~,u1]=ismember([Ni(k),Nj(m)],PUMAinfo(:,2:3),'rows');
                                else
                                    [~,u1]=ismember([Nj(m),Ni(k)],PUMAinfo(:,2:3),'rows');
                                end
                                if u1==0
                                    [k,m]
                                end
                                ODbyRoute{i,j}=[ODbyRoute{i,j},u1];
                                ODbyRoute{j,i}=ODbyRoute{i,j};
                            end
                        end
                    end
                end
            end
        
            % Aggregate all OD pairs in ODbyRoute
            TotalOD=[];
            for i=1:nRoutes
                for j=1:nRoutes
                    if AdjRoutesInd(i,j)==1
                        TotalOD=[TotalOD,ODbyRoute{i,j}];
                    end
                end
            end
            TotalOD=unique(TotalOD);
            estdemand=0;
            for i=1:length(TotalOD)
                if ismember(TotalOD(i),CorrPairMean(:,1))
                    estdemand=estdemand+CorrPairMean(CorrPairMean(:,1)==TotalOD(i),2);
                else
                    estdemand=estdemand+UncorrPairMean(UncorrPairMean(:,1)==TotalOD(i),2);
                end
            end
            EstDemand(r,vrblty)=estdemand;
        end
    end
%% 4. DISTRIBUTION
% [parmHat,parmCI] = wblfit(EstDemand(1:50,1))
end