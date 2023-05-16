function [Mu_est,OC,ObsOD,PrUncorr,PrCorrMean,PrecisU,PrecisC]=PolicyKG(Mu,Mu_0,Sigma_0,M,ObsOD,CoveredOD,ODSysNet,PrUncorr,TrUncorr,PrCorrMean,TrCorrMean,TrCorrCov,PrecisU,PrecisC)
%{
Knowledge Gradient with Independent Beliefs (KG)

[INPUT]
Mu: true mean value of options
Mu_0: prior mean value of options
sigma: true std.dev of options
Sigma_0: prior std.dev of options
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

K=length(Mu_0); %number of available choices
Mu_est=Mu_0;        % estimated mu from prior
OC=[];

% CREATE CovM, COVARIANCE MATRIX OF ACTIONS (CHOOSING A LINK), FROM CoveredOD
nCandLink=size(CoveredOD,2);

% CREATE beta_W, ARRAY OF MEASUREMENT PRECISION (CHOICE LEVEL)
beta_W=zeros(nCandLink,1);
for i=1:nCandLink
    beta_W(i,1)=1/Sigma_0(i,1);
end
beta=beta_W;

for k=1:M %try the kgcb for M number of times
    
    %Plogy is the log values of KG for alternatives
    KGAlts=zeros(K,8);
    for iter1=1:K % fill 1st four columns
        KGAlts(iter1,1)=Mu_est(iter1,1);    % mean value
        KGAlts(iter1,2)=beta(iter1,1);      % precision_n
        KGAlts(iter1,3)=KGAlts(iter1,2)+beta_W(iter1,1);    % precision_n+1 
        KGAlts(iter1,4)=sqrt(1/KGAlts(iter1,2)-1/KGAlts(iter1,3));  % change in variance
    end
    
    if K>1
        % 5th column (max value except for itself)
        Mu_est_sort=[(1:K)',Mu_est];
        [~,maxalt]=max(Mu_est);
        Mu_est_sort=sortrows(Mu_est_sort,2,'descend');
        for iter1=1:K
            if iter1==maxalt
                KGAlts(iter1,5)=Mu_est_sort(2,2);
            else
                KGAlts(iter1,5)=Mu_est_sort(1,2);
            end
        end

        % 6-8th column
        for iter1=1:K
            KGAlts(iter1,6)=-abs((KGAlts(iter1,1)-KGAlts(iter1,5))/KGAlts(iter1,4)); % normalized influence of decision
            KGAlts(iter1,7)=KGAlts(iter1,6)*normcdf(KGAlts(iter1,6))+normpdf(KGAlts(iter1,6)); % function of zeta
            KGAlts(iter1,8)=KGAlts(iter1,4)*KGAlts(iter1,7);    % knowledge gradient
        end

        [~,x]=max(KGAlts(:,8));
    else
        % only 1 option available
        x=1;
    end
    %max_value is the best estimated value of the KG 
    %x is the argument that produces max_value

    %{
    HOW CAN WE UPDATE mu AND covM REGARDING DIFFERENT OBSERVATIONS OF ODS?
    I THINK WE HAVE TO DRAW RANDOM NUMBERS FROM EACH OD FLOW DISTRIBUTION
    AND UPDATE CovOD 1ST AND CovM 
    %}
    
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
    Mu_est(x)=(Mu_est(x)*(beta(x))+w_k*(beta_W(x)))/(beta(x)+beta_W(x)); % update mean
    beta(x)=beta(x)+beta_W(x); % update precision
    
    %pick the best one to compare OC
    [~, max_choice]=max(Mu_est);

    %calculate the opportunity cost
    o_cost=max(Mu)-Mu(max_choice);
    
    OC=[OC,o_cost]; %update the OC matrix
end

% update priors 
nInvUncorr=length(InvolvedFlowU);
nInvCorr=length(InvolvedFlowC);
nExpUncorr=size(PrUncorr,1);
nExpCorr=size(PrCorrMean,1);
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
    if sum(isnan(PrUncorr(:,2)))>0
        123
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