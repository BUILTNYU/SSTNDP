function [Mu_est,OC,ObsOD,PrUncorr,PrCorrMean,PrecisU,PrecisC]=PolicyMAB(Mu,Mu_0,M,ObsOD,CoveredOD,ODSysNet,PrUncorr,TrUncorr,PrCorrMean,TrCorrMean,TrCorrCov,PrecisU,PrecisC)
%{
Multiarmed Bandit Algorithm (MAB)
: using Upper Confidence Bound Algorithm, inspired by Li, L., Lu, Y. and 
  Zhou, D., 2017, July. Provably optimal algorithms for generalized linear 
  contextual bandits. In International Conference on Machine Learning (pp. 
  2071-2080). PMLR.

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

    K=length(Mu_0); % number of available choices
    Mu_est=Mu_0;    % estimated mu from prior
    OC=[];
    
    % CREATE CovM, COVARIANCE MATRIX OF ACTIONS (CHOOSING A LINK), FROM CoveredOD
    nCandLink=size(CoveredOD,2);
    
    % RUN INITIAL nCandLink TRIALS TO TAKE AT LEAST ONE OBSERVATIONS FROM EACH
    ObsMAB=[Mu_est,ones(nCandLink,1)];
    
    % DETERMINE WHICH ACTION TO TAKE ACCORDING TO UCB ALGORITHM
    for k=1:M % explore after nCandLink observation
        ArmValue=zeros(K,1);
        for iter1=1:K
            ArmValue(iter1)=ObsMAB(iter1,1)+sqrt(2*log(k)/ObsMAB(iter1,2));
        end
        [~,x]=max(ArmValue);
        
        % RECEIVE OBSERVED OD FLOWS (RANDOM NUMBER REGARDING NORMAL DISTRIBUTION)
        % SINCE OUR PROBLEM INVOLVES DIFFERENT LEVEL OF FLOW, THE FUNCTION
        % SHOULD UPDATE THE INFORMATION ABOUT LOWER LEVEL (OD). UPDATE "ObsOD"
        % AND "PriorODMean".
        nObs=size(ObsOD,2);
    
        TotalOD = unique([CoveredOD{x},ODSysNet]);
        nCoverOD=size(TotalOD,2);
        
        MeanU=[];   p=0;
        StdevU=[];  
        InvolvedFlowU=[];
        InvolvedFlowC=[]; q=0;
        
        for i=1:nCoverOD
            if ismember(TotalOD(i),TrUncorr(:,1))
                p=p+1;
                InvolvedFlowU(p)=TotalOD(i);
                MeanU(p,1)=TrUncorr(TrUncorr(:,1)==TotalOD(i),2);
                StdevU(p,1)=TrUncorr(TrUncorr(:,1)==TotalOD(i),3);
            else
                q=q+1;
                InvolvedFlowC(q)=TotalOD(i);
            end
        end
        WU=normrnd(MeanU,StdevU);
        
        % Replace observation in W if correlated
        InvFl=size(InvolvedFlowC,2);
        MeanC=zeros(1,InvFl);
        CovC=zeros(InvFl);
        if InvFl>0
            for m=1:InvFl
                MeanC(m)=TrCorrMean(TrCorrMean(:,1)==InvolvedFlowC(m),2);
                for n=m:InvFl
                    % BRING COVARIANCES AND FORM A TEMPORARY COVARIANCE MATRIX
                    CovC(m,n)=TrCorrCov(TrCorrMean(:,1)==InvolvedFlowC(m),TrCorrMean(:,1)==InvolvedFlowC(n));
                    CovC(n,m)=CovC(m,n);
                end
            end
            WC=mvnrnd(MeanC,CovC);
        else
            WC=0;
        end
        
        % update ObsOD
        for i=1:p
            ObsOD{InvolvedFlowU(i),nObs+1}=WU(i);        
        end
        for i=1:q
            ObsOD{InvolvedFlowC(i),nObs+1}=WC(i);
        end
        w_k=sum(WU)+sum(WC);
    
        ObsMAB(x,1)=(ObsMAB(x,1)*ObsMAB(x,2)+w_k)/(ObsMAB(x,2)+1);
        ObsMAB(x,2)=ObsMAB(x,2)+1;
        
        %pick the best one to compare OC
        [~, max_choice]=max(ObsMAB(:,1));
    
        %calculate the opportunity cost
        o_cost=max(Mu)-Mu(max_choice);
        
        OC=[OC,o_cost]; %update the OC matrix
        Mu_est=ObsMAB(:,1);
    end
    
    % UPDATE "PriorODMean" and "CovOD" FOR OD PAIRS FROM ACCUMULATED "ObsOD"
    
    nInvUncorr=length(InvolvedFlowU);
    nInvCorr=length(InvolvedFlowC);
    for i=1:nInvUncorr
        t=find(PrUncorr(:,1)==InvolvedFlowU(i));
        if t>0 % true uncorrelated flow is in uncorrelated prior
            if PrUncorr(t,3)>0
                precisN=1/PrUncorr(t,3)^2;
                PrUncorr(t,2)=(precisN*PrUncorr(t,2)+PrecisU(t,1)*WU(i))/(precisN+PrecisU(t,1));
                PrUncorr(t,3)=sqrt(1/(precisN+PrecisU(t,1)));
                PrUncorr(t,4)=PrUncorr(t,4)+1;
            else
                if WU(i)>0
                    PrUncorr(t,2)=WU(i);
                    PrUncorr(t,3)=WU(i)*0.01;
                    PrUncorr(t,4)=PrUncorr(t,4)+1;
                    PrecisU(t,1)=1/(PrUncorr(t,3)^2);
                end
            end
        else % true uncorrelated flow is in correlated prior
            v=find(PrCorrMean(:,1)==InvolvedFlowU(i));
            if PrCorrMean(v,3)>0
                precisN=1/PrCorrMean(v,3)^2;
                PrCorrMean(v,2)=(precisN*PrCorrMean(v,2)+PrecisC(v,1)*WU(i))/(precisN+PrecisC(v,1));
                PrCorrMean(v,3)=sqrt(1/(precisN+PrecisC(v,1)));
                PrCorrMean(v,4)=PrCorrMean(v,4)+1;
            else
                if WU(i)>0
                    PrCorrMean(v,2)=WU(i);
                    PrCorrMean(v,3)=WU(i)*0.01;
                    PrCorrMean(v,4)=PrCorrMean(v,4)+1;
                    PrecisC(v,1)=1/(PrCorrMean(v,3)^2);
                end
            end
        end
    end
    
    for i=1:nInvCorr
        t=find(PrCorrMean(:,1)==InvolvedFlowC(i));
        if t>0 % true correlated flow is in correlated prior
            if PrCorrMean(t,3)>0
                precisN=1/PrCorrMean(t,3)^2;
                PrCorrMean(t,2)=(precisN*PrCorrMean(t,2)+PrecisC(t,1)*WC(i))/(precisN+PrecisC(t,1));
                PrCorrMean(t,3)=sqrt(1/(precisN+PrecisC(t,1)));
                PrCorrMean(t,4)=PrCorrMean(t,4)+1;
            else
                if WC(i)>0
                    PrCorrMean(t,2)=WC(i);
                    PrCorrMean(t,3)=WC(i)*0.01;
                    PrCorrMean(t,4)=PrCorrMean(t,4)+1;
                    PrecisC(t,1)=1/(PrCorrMean(t,3)^2);
                end
            end
        else % true correlated flow is not in correlated prior
            v=find(PrUncorr(:,1)==InvolvedFlowC(i));
            if PrUncorr(v,3)>0
                precisN=1/PrUncorr(v,3)^2;
                PrUncorr(v,2)=(precisN*PrUncorr(v,2)+PrecisU(v,1)*WC(i))/(precisN+PrecisU(v,1));
                PrUncorr(v,3)=sqrt(1/(precisN+PrecisU(v,1)));
                PrUncorr(v,4)=PrUncorr(v,4)+1;
            else
                if WC(i)>0
                    PrUncorr(v,2)=WC(i);
                    PrUncorr(v,3)=WC(i)*0.01;
                    PrUncorr(v,4)=PrUncorr(v,4)+1;
                    PrecisU(v,1)=1/(PrUncorr(v,3)^2);
                end
            end
        end
    end
end