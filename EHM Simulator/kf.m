function [xpred,xpost,deltax] = kf(NT,A,B,C,R,Q,ysim,u)

    xpred = zeros(2,NT);
    xpred(:,1) = 0.5*ones(2,1);
    xpost = zeros(2,NT);
    Sigma = [1e-3 0; 0 1e-5];
    deltax = zeros(2,NT);

    for k = 1:(NT-1)
        % Measurement update
        Kf = Sigma*C'/(R+C*Sigma*C');
        xpost(:,k) = xpred(:,k) + Kf*(ysim(k)-C*xpred(:,k));
        Sigma_post = Sigma - Sigma*C'*((R+C*Sigma*C')\(C*Sigma));
        deltax(:,k) = xpred(:,k)-xpost(:,k);

        % Time update
        xpred(:,k+1) = A*xpost(:,k) + B*u(:,k);
        Sigma = A*Sigma_post*A' + Q;
    end
end