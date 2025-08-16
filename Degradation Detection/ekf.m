function [xpost,y,innov,S,traceSigmapost] = ekf(p,NT,ysim,u,CI,v,KF)

    syms SOCn CSCn Cur
    A = [1 0; p.gnh/(p.Beta*(1-p.Beta)) 1-p.gnh/(p.Beta*(1-p.Beta))];
    B = [1; 1/(1-p.Beta)]*p.aa0ne;
    [h,H] = JacobianMeasMat(p);
    
    R = v*v;
    Q = diag([1e-6 1e-6]);
    
    xpred = zeros(2,NT); xpred(1:2,1) = CI;
    xpost = zeros(2,NT); xpost(1:2,1) = CI;
    Sigma = diag([2.5e-1 2.5e-1]);
    Sigma_post = Sigma;
    traceSigmapost = zeros(1,NT); traceSigmapost(1) = trace(Sigma_post);
    y = zeros(1,NT); y(1) = double(subs(h,[SOCn CSCn Cur],[xpred(1,1),xpred(2,1),u(1)]));
    innov = zeros(1,NT); innov(1) = ysim(1)-y(1); %innovation sequence
    Hk = double(subs(H,[SOCn CSCn Cur],[xpred(1,1),xpred(2,1),u(1)]));
    S = zeros(1,NT); S(1) = R+Hk*Sigma*Hk'; %innovation covariance matrix sequence
    
    % fprintf(1,'k: %3.2f | SOC-(k)/CSC-(k): %2.2f/%2.2f | V(k): %2.2f V | KF: %2.2f',...
    %     1,xpost(1,1),xpost(2,1),y(1),KF);
    fprintf(1,'k: %2.2f | KF: %2.2f',1,KF);
    disp(".");
    
    for k = 2:(NT-1)
        % Time update
        xpred(:,k) = A*xpost(:,k-1) + B*u(k-1);
        if xpred(1,k) < 0
            xpred(1,k) = 0;
        end
        if xpred(2,k) < 0
            xpred(2,k) = 0;
        end
        if xpred(1,k) > 1
            xpred(1,k) = 1;
        end
        if xpred(2,k) > 1
            xpred(2,k) = 1;
        end
        Sigma = A*Sigma_post*A' + Q;
        Hk = double(subs(H,[SOCn CSCn Cur],[xpred(1,k),xpred(2,k),u(k)]));
        y(k) = double(subs(h,[SOCn CSCn Cur],[xpred(1,k),xpred(2,k),u(k)]));
        innov(k) = ysim(k)-y(k);
        S(k) = R+Hk*Sigma*Hk';
    
        % Measurement update
        Kf = Sigma*Hk'/S(k);
        xpost(:,k) = xpred(:,k) + Kf*innov(k);
        if xpost(1,k) < 0
            xpost(1,k) = 0;
        end
        if xpost(2,k) < 0
            xpost(2,k) = 0;
        end
        if xpost(1,k) > 1
            xpost(1,k) = 1;
        end
        if xpost(2,k) > 1
            xpost(2,k) = 1;
        end
        Sigma_post = Sigma - Kf*S(k)*Kf';
        traceSigmapost(k) = trace(Sigma_post);
    
        % fprintf(1,'k: %3.2f | SOC-(k)/CSC-(k): %2.2f/%2.2f | V(k): %2.2f V | KF: %2.2f',...
        %     k,xpost(1,k),xpost(2,k),y(k),KF);
        fprintf(1,'k: %2.2f | KF: %2.2f',k,KF);
        disp(".");
    end
end