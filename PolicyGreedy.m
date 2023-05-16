function [Mu_est,OC,ObsOD,PrUncorr,PrCorrMean,PrecisU,PrecisC]=PolicyGreedy(Mu,Mu_0,M,ObsOD,CoveredOD,ODSysNet,PrUncorr,TrUncorr,PrCorrMean,TrCorrMean,TrCorrCov,PrecisU,PrecisC)
%{
Greedy Policy: expand the system to the segment with the highest reward
               based on the current knowledge, ignore exploration

[INPUT]
Mu: true mean value of options
Mu_0: prior mean value of options
M: # of trials with the same option set
ObsOD: observed demand per OD pair
CoveredOD: IDs of OD pairs covered per option
ODSysNet: IDs of OD pairs covered by current system
PrUncorr: prior mean vector of uncorrelated flows
TrUncorr: true mean vector of uncorrelated flows
PrCorrMean: prior mean vector of correlated flows
TrCorrMean: true mean vector of correlated flows
TrCorrCov: true covariance matrix of correlated flows
PrecisU: observation precision of uncorrelated flows
PrecisC: observation precision of correlated flows

[OUTPUT]
Mu_est: estimated mean value of options
OC: opportunity cost
ObsOD: updated
PrUncorr: updated
PrCorrMean: updated
PrecisU: updated
PrecisC: updated
%}

    OC = [];    % space for opportunity cost
   
    % CREATE CovM, COVARIANCE MATRIX OF ACTIONS (CHOOSING A LINK), FROM CoveredOD
    nCandLink=size(CoveredOD,2);
    
    % RUN INITIAL nCandLink TRIALS TO TAKE AT LEAST ONE OBSERVATIONS FROM EACH
    ObsGRDY = [Mu_0,ones(nCandLink,1)];
    
    % GREEDILY DETERMINE WHICH ACTION TO TAKE
    for k = 1:M 
        [~,x]=max(ObsGRDY(:,1));
        
        % RECEIVE OBSERVED OD FLOWS (RANDOM NUMBER REGARDING NORMAL DISTRIBUTION)
        % SINCE OUR PROBLEM INVOLVES DIFFERENT LEVEL OF FLOW, THE FUNCTION
        % SHOULD UPDATE THE INFORMATION ABOUT LOWER LEVEL (OD). UPDATE "ObsOD"
        % AND "PriorODMean".
        nObs = size(ObsOD,2);
        TotalOD = unique([CoveredOD{x},ODSysNet]);
        nCoverOD = size(TotalOD,2);
        
        MeanU = [];   p = 0;
        StdevU = [];  
        InvolvedFlowU = [];
        InvolvedFlowC = []; q = 0;
        
        for i=1:nCoverOD
            if ismember(TotalOD(i),TrUncorr(:,1))
                p = p+1;
                InvolvedFlowU(p) = TotalOD(i);
                MeanU(p,1) = TrUncorr(TrUncorr(:,1) == TotalOD(i),2);
                StdevU(p,1) = TrUncorr(TrUncorr(:,1) == TotalOD(i),3);
            else
                q = q+1;
                InvolvedFlowC(q) = TotalOD(i);
            end
        end
        WU = normrnd(MeanU,StdevU);
        
        % Replace observation in W if correlated
        InvFl = size(InvolvedFlowC,2);
        MeanC = zeros(1,InvFl);
        CovC = zeros(InvFl);
        if InvFl > 0
            for m = 1:InvFl
                MeanC(m) = TrCorrMean(TrCorrMean(:,1) == InvolvedFlowC(m),2);
                for n = m:InvFl
                    % BRING COVARIANCES AND FORM A TEMPORARY COVARIANCE MATRIX
                    CovC(m,n) = TrCorrCov(TrCorrMean(:,1) == InvolvedFlowC(m),TrCorrMean(:,1) == InvolvedFlowC(n));
                    CovC(n,m) = CovC(m,n);
                end
            end
            WC = mvnrnd(MeanC,CovC);
        else
            WC = 0;
        end
        
        % update ObsOD
        for i = 1:p
            ObsOD{InvolvedFlowU(i),nObs+1} = WU(i);        
        end
        for i = 1:q
            ObsOD{InvolvedFlowC(i),nObs+1} = WC(i);
        end
        w_k = sum(WU)+sum(WC);
    
        ObsGRDY(x,1) = (ObsGRDY(x,1)*ObsGRDY(x,2)+w_k)/(ObsGRDY(x,2)+1);
        ObsGRDY(x,2) = ObsGRDY(x,2)+1;
      
        %pick the best one to compare OC
        [~, max_choice]=max(ObsGRDY(:,1));
    
        %calculate the opportunity cost
        o_cost = max(Mu)-Mu(max_choice);
        OC = [OC,o_cost]; %update the OC matrix
        Mu_est = ObsGRDY(:,1);
    end
    
    % UPDATE "PriorODMean" and "CovOD" FOR OD PAIRS FROM ACCUMULATED "ObsOD"
    nInvUncorr = length(InvolvedFlowU);
    nInvCorr = length(InvolvedFlowC);
    nExpUncorr = size(PrUncorr,1);
    nExpCorr = size(PrCorrMean,1);
    for i = 1:nInvUncorr
        t = find(PrUncorr(:,1) == InvolvedFlowU(i));
        if t > 0 % true uncorrelated flow is in uncorrelated prior
            if PrUncorr(t,3) > 0
                precisN = 1/PrUncorr(t,3)^2;
                PrUncorr(t,2) = (precisN*PrUncorr(t,2)+PrecisU(t,1)*WU(i))/(precisN+PrecisU(t,1));
                PrUncorr(t,3) = sqrt(1/(precisN+PrecisU(t,1)));
                PrUncorr(t,4) = PrUncorr(t,4)+1;
            else
                if WU(i) > 0
                    PrUncorr(t,2) = WU(i);
                    PrUncorr(t,3) = WU(i)*0.01;
                    PrUncorr(t,4) = PrUncorr(t,4)+1;
                    PrecisU(t,1) = 1/(PrUncorr(t,3)^2);
                end
            end
        else % true uncorrelated flow is in correlated prior
            v = find(PrCorrMean(:,1) == InvolvedFlowU(i));
            if PrCorrMean(v,3) > 0
                precisN = 1/PrCorrMean(v,3)^2;
                PrCorrMean(v,2) = (precisN*PrCorrMean(v,2)+PrecisC(v,1)*WU(i))/(precisN+PrecisC(v,1));
                PrCorrMean(v,3) = sqrt(1/(precisN+PrecisC(v,1)));
                PrCorrMean(v,4) = PrCorrMean(v,4)+1;
            else
                if WU(i) > 0
                    PrCorrMean(v,2) = WU(i);
                    PrCorrMean(v,3) = WU(i)*0.01;
                    PrCorrMean(v,4) = PrCorrMean(v,4)+1;
                    PrecisC(v,1) = 1/(PrCorrMean(v,3)^2);
                end
            end
        end
    end
    
    for i = 1:nInvCorr
        t = find(PrCorrMean(:,1) == InvolvedFlowC(i));
        if t > 0 % true correlated flow is in correlated prior
            if PrCorrMean(t,3) > 0
                precisN = 1/PrCorrMean(t,3)^2;
                PrCorrMean(t,2) = (precisN*PrCorrMean(t,2)+PrecisC(t,1)*WC(i))/(precisN+PrecisC(t,1));
                PrCorrMean(t,3) = sqrt(1/(precisN+PrecisC(t,1)));
                PrCorrMean(t,4) = PrCorrMean(t,4)+1;
            else
                if WC(i) > 0
                    PrCorrMean(t,2) = WC(i);
                    PrCorrMean(t,3) = WC(i)*0.01;
                    PrCorrMean(t,4) = PrCorrMean(t,4)+1;
                    PrecisC(t,1) = 1/(PrCorrMean(t,3)^2);
                end
            end
        else % true correlated flow is not in correlated prior
            v = find(PrUncorr(:,1) == InvolvedFlowC(i));
            if PrUncorr(v,3) > 0
                precisN = 1/PrUncorr(v,3)^2;
                PrUncorr(v,2) = (precisN*PrUncorr(v,2)+PrecisU(v,1)*WC(i))/(precisN+PrecisU(v,1));
                PrUncorr(v,3) = sqrt(1/(precisN+PrecisU(v,1)));
                PrUncorr(v,4) = PrUncorr(v,4)+1;
            else
                if WC(i)>0
                    PrUncorr(v,2) = WC(i);
                    PrUncorr(v,3) = WC(i)*0.01;
                    PrUncorr(v,4) = PrUncorr(v,4)+1;
                    PrecisU(v,1) = 1/(PrUncorr(v,3)^2);
                end
            end
        end
    end

end