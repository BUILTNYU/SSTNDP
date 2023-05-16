function [Mu_est,OC,ObsOD,PrUncorr,PrCorrMean,PrCorrCov,PrecisU,PrecisC]=PolicyKGCB(Mu,Mu_0,M,ObsOD,CoveredOD,ODSysNet,PrUncorr,TrUncorr,PrCorrMean,PrCorrCov,TrCorrMean,TrCorrCov,ODChoiceIndCorr,PrecisU,PrecisC)
%{
Knowledge Gradient with Correlated Beliefs (KGCB)
: inspired by Software (Castle Labs) 
  http://castlelab.princeton.edu/html/software/kg/KGCB.zip

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
    nCandLink=size(CoveredOD,2);
    
    % ADD ALL OD FLOWS
    for i=1:nCandLink
        CoveredOD{i} = unique([CoveredOD{i},ODSysNet]);
    end
    
    % CREATE CovM, COVARIANCE MATRIX OF ACTIONS (CHOOSING A LINK), FROM CoveredOD
    CovM=zeros(nCandLink);
    for i=1:nCandLink
        for j=i:nCandLink
            for m=1:size(CoveredOD{i},2)
                for n=1:size(CoveredOD{j},2)
                    odim=CoveredOD{i}(m);
                    odjn=CoveredOD{j}(n);
                    if ismember(odim,PrCorrMean(:,1))&&ismember(odjn,PrCorrMean(:,1)) % both are considered as correlated
                        CovM(i,j)=CovM(i,j)+PrCorrCov(PrCorrMean(:,1)==odim,PrCorrMean(:,1)==odjn);
                    else % both are not considered as correlated
                        if odim==odjn
                            CovM(i,j)=CovM(i,j)+PrUncorr(PrUncorr(:,1)==odim,3)^2;
                        end
                    end
                end
            end
            if i~=j
                CovM(j,i)=CovM(i,j);
            end
            if i==j && CovM(i,j)<0
                i,j
            end
        end
    end
    
    % % CREATE beta_W, ARRAY OF MEASUREMENT PRECISION (CHOICE LEVEL)
    beta_W=zeros(nCandLink,1);
    for i=1:nCandLink
        beta_W(i,1)=1/CovM(i,i);
    end
    CovCheck=struct; % for saving previous version of CovM
    for k=1:M %try the kgcb for M number of times
        % CREATE beta_W, ARRAY OF MEASUREMENT PRECISION (CHOICE LEVEL)
    
        %Plogy is the log values of KG for alternatives
        Plogy=[];
        Py=[];
        
        for iter1=1:K
            a=Mu_est';  % choice-level mean
            b=CovM(iter1,:)/sqrt(1/beta_W(iter1) + CovM(iter1,iter1));  % choice-level variation
            if isreal(b)~=1 || size(a,2)~=size(b,2)
                iter1
            end
            [KG,LogKG] = LogEmaxAffine(a,b);
            
            [Plogy]=[Plogy, LogKG];
            [Py]=[Py,KG];
        end
        
        [~,x]=max(Plogy);
        
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
        
        TotalOD = CoveredOD{x};
        nCoverOD=size(TotalOD,2);
        
        MeanU=[];   p=0;
        StdevU=[];  
        InvolvedFlowU=[];
        InvolvedFlowC=[]; q=0;
        
        for i=1:nCoverOD % number of OD pairs covered by the best
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
        
        % Extract OD pairs should be used in estimating CovM
        
        w_k=sum(WU)+sum(WC);
        CovCheck(k).cov=CovM;
        e_x=zeros(K,1);       % prepare unit column vector
        e_x(x)=1;
        addscalar = (w_k - Mu_est(x))/(1/beta_W(x) + CovM(x,x));
        Mu_est=Mu_est + addscalar*CovM*e_x;
        CovM = CovM - ((CovM*e_x) * (e_x' * CovM)) / ((1 / beta_W(x)) + CovM(x,x));
        
        
        for i=1:nCandLink
            if CovM(i,i)<0
                [k,i]
            end
        end
        
        
        %pick the best one to compare OC
        [~, max_choice]=max(Mu_est);
    
        %calculate the opportunity cost
        o_cost=max(Mu)-Mu(max_choice);
        
        OC=[OC,o_cost]; %update the OC matrix
    end
    
    % UPDATE "PriorODMean" and "CovOD" FOR OD PAIRS FROM ACCUMULATED "ObsOD"
    % Prior mean, std.dev, and precision of uncorrelated flows
    nInvUncorr=length(InvolvedFlowU);
    nInvCorr=length(InvolvedFlowC);
    nExpCorr=size(PrCorrMean,1);
    Y=zeros(nExpCorr,1);
    
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
            Y(PrCorrMean(:,1)==InvolvedFlowU(i),1)=WU(i);
        end
    end
    
    for i=1:nInvCorr
        t=find(PrCorrMean(:,1)==InvolvedFlowC(i));
        if t>0 % true correlated flow is in correlated prior
            Y(t,1)=WC(i);
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
    
    OmegaC=zeros(nExpCorr);
    for i=1:nExpCorr
        OmegaC(i,i)=ODChoiceIndCorr(i,1)*PrecisC(i,1);
    end
    PrCorrMean(:,2)=(OmegaC+PrCorrCov^-1)^-1*(Y'*OmegaC+PrCorrMean(:,2)'*PrCorrCov^-1)';
    PrCorrCov=(OmegaC+PrCorrCov^-1)^-1;
end

% logy = LogEmaxAffine(a,b)
% Calculates log(Exp[max_x a_x + b_x Z]-max_x a_x), where Z is a standard
% normal random variable and a,b are 1xM input vectors.
function [y,logy, a,b,c] = LogEmaxAffine(a,b)
    if (any(isnan(a)) || any(isnan(b)))
        warning('a or b is NaN');
    end
    assert(all(isreal(a)));
    assert(all(isreal(b)));

    a = a';
    b = b';

    % Check that a and b are column vectors of the right size
    if (any(size(a) ~= size(b)))
        error('LogEmaxAffine: a and b must be column vectors of the same size');
    end
    
    [a,b] = AffineBreakpointsPrep(a,b);
    
    [c, keep] = AffineBreakpoints(a,b); 
    a = a(keep);
    b = b(keep);
    c = c([1,keep+1]);
    M = length(keep);
    assert(all(isreal(c)));

    % I need logbdiff=log(diff(b)).  I thought that the following code would be
    % more numerically stable, able for example to distinguish cases like 
    % logb = [-25 -.3] vs. logb = [-35 -.3], but it doesn't seem to be able to.
    % Indeed, in the debugging output that I have below, the difference was 0.
    %{
    logb = log(abs(b)); % If b is 0, this is -Inf.
    sgnb = sign(b); % If b is 0, this is 0.
    logbdiff = zeros(size(c(2:M)));
    for i=1:length(b)-1
	[logbdiff(i),logbdiffsgn] = LogPlusExpSigned(logb(i),sgnb(i),logb(i+1),-sgnb(i+1));
	%assert(logbdiffsgn>=0);  % The b are distinct, so bdiff(i) can't be 0.
    end
    disp(sprintf('log(b)=%s log(diff(b))=%g logbdiff=%g difference=%g',mat2str(log(b)),log(diff(b)),logbdiff,log(diff(b))-logbdiff));
    %}
    logbdiff = log(diff(b))';  
    
    if M==1
        logy=log(a);
    elseif M>=2
        logy = LogSumExp(logbdiff+LogEI(-abs(c(2:M))));
    end

    logy=real(logy);
    y=exp(logy);
end
% Prepares vectors for passing to AffineEmaxBreakpoints, changing their
% order and removing elements with duplicate slope.

function [a,b] = AffineBreakpointsPrep(a,b)
    % Make sure a and b are column vectors.
    rows = size(a); if (rows == 1), a=a'; end
    rows = size(b); if (rows == 1), b=b'; end
    
    % 11/29/2008 PF: Experimental preprocessing step, which I hope will remove
    % a large number of the entries.
    [b1, i1] = min(b); % [a1,b1] is best at z=-infinity
    [a2, i2] = max(a); % [a2,b2] is best at z=0
    [b3, i3] = max(b); % [a3,b3] is best at z=+infinity
    a1 = a(i1);
    b2 = b(i2);
    a3 = a(i3);
    cleft = (a - a1)./(b1 - b); % intersection with leftmost line. 
    cright = (a - a3)./(b3 - b); % intersection with rightmost line.
    c2left = (a2 - a1)./(b1 - b2); % intersection with leftmost line. 
    c2right = (a2 - a3)./(b3 - b2); % intersection with rightmost line.
    keep = find(b==b1 | b==b3 | cleft <= c2left | cright >= c2right);
    %disp(sprintf('Preprocessing cut %d of %d entries', length(a)-length(keep), length(a)));
    a = a(keep);
    b = b(keep);
    clear keep cleft cright
   
        
    % Form a matrix for which ba(x,1) is the slope b(x) and ba(x,2) is the
    % y-intercept a(x).  Sort this matrix in ascending order of slope, 
    % breaking ties in slope with the y-intercept.  
    ba = [b, a];
    ba = sortrows(ba,[1,2]);
    a = ba(:,2);
    b = ba(:,1);
    
    % Then, from each pair of indices with the b component equal, remove
    % the one with smaller a component.  This code works because the sort
    % above enforced the condition: if b(i) == b(i+1), then a(i) <= a(i+1).
    keep = [find(diff(b)); length(b)];
    % This previous line is equivalent to:
    % keep = [];
    % for i=[1:length(b)-1]
    %    if b(i)~=b(i+1)
    %        keep = [keep, i];
    %    end
    %end 
    %keep = [keep, length(b)];  % We always keep the last one.
    
    % Note that the elements of keep are in ascending order.
    % This makes it so that b(keep) is still sorted in ascending order.
    a = a(keep);
    b = b(keep);
end

% Inputs are two M-vectors, a and b.
% Requires that the b vector is sorted in increasing order.
% Also requires that the elements of b all be unique.
% This function is used in AffineEmax, and the preparation of generic
% vectors a and b to satisfy the input requirements of this function are
% shown there.
%
% The output is an (M+1)-vector c and a vector A ("A" is for accept).  Think of
% A as a set which is a subset of {1,...,M}.  This output has the property
% that, for any i in {1,...,M} and any real number z,
%   i \in argmax_j a_j + b_j z
% iff
%   i \in A and z \in [c(j+1),c(i+1)],
%   where j = sup {0,1,...,i-1} \cap A.
%
% A note about indexing:
% Since Matlab does not allow indexing from 0, but instead requires
% indexing from 1, what is called c_i in the paper is written in matlab as
% c(1+i).  This is because in the paper we reference c_0.  For the vectors
% a and b, however, we don't need to reference a_0 or b_0, so we reference
% a_i and b_i by a(i) and b(i) respectively, rather than a(i+1) or b(i+1).
% 
function [c,A] = AffineBreakpoints(a,b)
    % Preallocate for speed.  Instead of resizing the array A whenever we add
    % to it or delete from it, we keep it the maximal size, and keep a length
    % indicator Alen telling us how many of its entries are good.  When the
    % function ends, we remove the unused elements from A before passing
    % it.
    M = length(a);
    c = zeros(1,M+1);
    A = zeros(1,M);
    
    % Step 0
    i=0;
    c(1+i) = -inf;
    c(1+i+1) = +inf;
    A(1) = 1;
    Alen = 1;
    
    for i=[1:M-1]
        c(1+i+1) = +inf;
        while(1)
            j = A(Alen); % jindex = Alen
            c(1+j) = (a(j) - a(i+1))/(b(i+1)-b(j));
	    % The if statement below replaces these lines from version 2 of the
	    % function.
	    %    kindex = jindex-1 = Alen-1
            %    if kindex > 0 && c(1+j)<=c(1+A(kindex))
            if Alen > 1 && c(1+j)<c(1+A(Alen-1))
		Alen = Alen-1; % Remove last element j
                % continue in while(1) loop
            else
                break % quit while(1) loop
            end
        end
	A(Alen+1) = i+1;
	Alen = Alen + 1;
    end
    A = A(1:Alen);
end

% Returns the log of Exp[(s+Z)^+], where s is a constant and Z is a standard
% normal random variable.  For large negative arguments Exp[(s+Z)^+] function
% is close to 0.  For large positive arguments, the function is close to the
% argument.  For s large enough, s>-10, we use the formula
% Exp[(s+Z)^+] = s*normcdf(s) + normpdf(s).  For smaller s we use an asymptotic
% approximation based on Mill's ratio.  EI stands for "expected improvement",
% since Exp[(s+Z)^+] would be the log of the expected improvement by measuring
% an alternative with excess predictive mean s over the best other measured
% alternative, and predictive variance 0.
function logy = LogEI(s)

% Use the asymptotic approximation for these large negative s.  The
% approximation is derived via:
%   s*normcdf(s) + normpdf(s) = normpdf(s)*[1-|s|normcdf(-|s|)/normpdf(s)]
% and noting that normcdf(-|s|)/normpdf(s) is the Mill's ratio at |s|, which is
% asymptotically approximated by |s|/(s^2+1) [Gordon 1941, also documented in
% Frazier,Powell,Dayanik 2009 on page 14].  This gives,
%   s*normcdf(s) + normpdf(s) = normpdf(s)*[1-s^2/(s^2+1)] = normpdf(s)/(s^2+1).

i=find(s<-10);
if (length(i)>0)
    logy(i) = LogNormPDF(s(i)) - log(s(i).^2 + 1);
end

% Use straightforward routines for s in the more numerically stable region.
i=find(s>=-10);
if (length(i)>0)
    logy(i) = log(s(i).*normcdf(s(i))+normpdf(s(i)));
end

assert(all(isreal(logy)));
end

% logy = LogNormPDF(z)
% Returns the log of the normal pdf evaluated at z.  z can be a vector or a scalar.
function logy = LogNormPDF(z)
	const = -.5*log(2*pi); % log of 1/sqrt(2pi).
	logy = const - z.^2/2;
end

% function y=LogSumExp(x)
% Computes log(sum(exp(x))) for a vector x, but in a numerically careful way.
function y=LogSumExp(x)
xmax = max(x);
y = xmax + log(sum(exp(x-xmax)));
end
