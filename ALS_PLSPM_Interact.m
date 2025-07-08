function [WT, LD, Path,it,Flag_Converge, Gamma] = ALS_PLSPM_Interact(z0,W0,B0,B02,modetype,scheme,nnlv_index,ind_sign,Max_iter,Min_limit,N,J,P)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ALS_PlSPM_Interact() - MATLAB function to implement the ALS algorithm for%
%                Partial Least Squares Path Modeling (PLSPM) with latent  %
%                interactions.                                            % 
% Author: Heungsun Hwang & Gyeongcheol Cho                                % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

z = zscore(z0,1)./sqrt(N); 
W = double(W0);
B = double(B0);
Gamma = z*W;
Gg_sd=sqrt(diag(Gamma'*Gamma)'./N);
W=W./repmat(Gg_sd,J,1);
Gamma = zeros(N,P);                     
for p = 1:P
    windex = find(W0(:,p));
    z_j = z(:,windex);
    Gamma(:,p) = z_j*W(windex,p);
    Gamma(:,p) = Gamma(:,p)/norm(Gamma(:,p));
end

if scheme == 1              %centroid scheme
      corLV = corrcoef(Gamma);
      for p = 1:P
          bindex = B0(:,p);   % DV 
         if sum(bindex,1)>0
            B(bindex,p) = sign(corLV(bindex,p));  
         end
      end 
      for p = 1:P
         bindex = B0(p,:);   % IV 
         if sum(bindex,2)>0
            B(bindex,p) = sign(corLV(p,bindex));
         end
      end
elseif scheme == 2         % factorial scheme
      corLV = corrcoef(Gamma);
      for p = 1:P
          bindex = B0(:,p);   % DV 
         if sum(bindex,1)>0
            B(bindex,p) = corLV(bindex,p);  
         end
      end 
      for p = 1:P
         bindex = B0(p,:);   % IV 
         if sum(bindex,2)>0
            B(bindex,p) = corLV(p,bindex);
         end
      end
elseif scheme == 3       % path weighting scheme
       for p = 1:P
           bindex = B0(:,p);   % DV
           if sum(bindex,1)>0
              gp = Gamma(:,bindex);
              B(bindex,p) = (gp'*gp)\gp'*Gamma(:,p);  
           end
       end    
       corLV = corrcoef(Gamma);
       for p = 1:P
           bindex = B0(p,:);   % IV
           if sum(bindex,2)>0
               B(bindex,p) = corLV(p,bindex);
           end
       end
end

it = 0;                   % iteration counter
imp = 100000;             % initial improvement
while it < Max_iter && imp > Min_limit 
      it = it+1;   
      W_old = W;
      %% update weighted composites of latents
      F = Gamma*B;    
      %% update weights and latents 
      for p = 1:P
          windex = W0(:,p);
          Jp=sum(windex,1);
          alpha = 2 - modetype(p);
          beta =  1 - alpha;
%           H = eye(P);               
%           H(p,p) = 0;
%           F1 = z - F*H*W';
%           y1 = reshape(F1,N*J,1);
%           X1 = kron(eye(J),F(:,p));
%           X1 = X1(:,windex);
          y1 = reshape(z(:,windex),N*Jp,1);
          X1 = kron(eye(Jp),F(:,p));
          XX1 = alpha*(X1'*X1);
          Xy1 = alpha*X1'*y1;
          XX2 = beta*(z(:,windex)'*z(:,windex));
          Xy2 = beta*z(:,windex)'*F(:,p);
          w_p = (XX1 + XX2)\(Xy1 + Xy2);
          W(windex,p) = w_p;
          Gamma(:,p) = z(:,windex)*w_p; 
          Gamma(:,p) = Gamma(:,p)/norm(Gamma(:,p));
      end      
       %% update B 
       if scheme == 1        % centroid scheme
        corLV = corrcoef(Gamma);
        for p = 1:P
            bindex = B0(:,p);   % DV 
            if sum(bindex,1)>0
               B(bindex,p) = sign(corLV(bindex,p));  
            end
            bindex = B0(p,:);   % IV 
            if sum(bindex,2)>0
               B(bindex,p) = sign(corLV(p,bindex));
            end
        end
     elseif scheme == 2    % factorial scheme
            corLV = corrcoef(Gamma);
            for p = 1:P
                bindex = B0(:,p);   % DV 
                if sum(bindex,1)>0
                    B(bindex,p) = corLV(bindex,p);  
                end
            end 
            for p = 1:P
                bindex = B0(p,:);   % IV 
                if sum(bindex,2)>0
                   B(bindex,p) = corLV(p,bindex);
                end
            end
     elseif scheme == 3     % path weighting scheme 
            for p = 1:P
                bindex = B0(:,p);   % DV
                if sum(bindex,1)>0
                   gp = Gamma(:,bindex);
                   B(bindex,p) = (gp'*gp)\gp'*Gamma(:,p);  
                end
            end    
            corLV = corrcoef(Gamma);
            for p = 1:P
                bindex = B0(p,:);   % IV
                if sum(bindex,2)>0
                   B(bindex,p) = corLV(p,bindex);
                end
             end
      end
      %% check convergence (based on V. Esposito Vinzi et al., 2010)
%       f1 = 0;
%       f2 = 0;
%       for p = 1:P
%           windex = find(W0(:,p));
%           alpha = modetype(p);
%           beta = 1 - alpha;
%           dif1 = z(:,windex) - Gamma*B(:,p)*W(windex,p)';
%           f1 = f1 + alpha*trace(dif1'*dif1);
%           dif2 = F(:,p) - Gamma*B(:,p);
%           f2 = f2 + beta*(dif2'*dif2);
%       end
%      f = f1 + f2;
      imp = max(max(W - W_old));
%      imp = f0 - f;
%      [it f f0 imp];
%      f0 = f;
end
Flag_Converge=true;
if (it== Max_iter) && (imp > Min_limit); Flag_Converge=false; end
%% Estimation of loadings, weights, and paths 
%Path = zeros(P,P);      % path coefficients
WT = zeros(J,P);       % outer weights
LD = zeros(P,J);       % outer loadings

for p=1:P
    if ind_sign(1,p)>0
        if z(:,ind_sign(1,p))'*Gamma(:,p)<0
            Gamma(:,p)=-Gamma(:,p);
        end
    end
end
for p = 1:P
%    bindex = B0(:,p);
%    if sum(bindex,1)>0
%       Gamma_p = Gamma(:,bindex);
%       Path(bindex,p) = pinv(Gamma_p'*Gamma_p)*Gamma_p'*Gamma(:,p);
%    end
    windex = W0(:,p);
    z_p = z(:,windex);
    WT(windex,p) = pinv(z_p'*z_p)*z_p'*Gamma(:,p);
    LD(p,windex) = Gamma(:,p)'*z_p;
end

%% Stage 2 (Linear regression with an interaction term)
Gamma = Gamma*sqrt(N);              % standardized component scores
NLV = lint_score(nnlv_index,Gamma);
Path = double(B02);
G = [Gamma,(NLV-mean(NLV))];
for j = 1:P
    bindex = B02(:,j);
    if sum(bindex,1)>0
       X = G(:,bindex);
       y = G(:,j);
       Path(bindex,j) = pinv(X'*X)*X'*y;
    end
end
end
