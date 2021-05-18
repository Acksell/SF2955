%% 
tau=load("coal-mine.csv");

%%

% Gibbs sampling + MH
rng(42)



params=struct();
params.d=5;
params.vartheta=2;

t=linspace(1851,1963, params.d+1);

params.ndisasters=ndisasters(tau,t);
params.rho=5;

% t=randn(1,d);
% t=[1851 t 1963];


lambda=ones(1, params.d);
theta=(4)/(sum(lambda) + params.vartheta);

K=50000;
thetaStore=zeros(K,1);
lambdaStore=zeros(K,params.d);
tStore=zeros(K,params.d+1);
for i=1:K
    [theta,lambda,t]=GibbsStep(theta, lambda, t, params);
    thetaStore(i)=theta;
    lambdaStore(i,:)=lambda;
    tStore(i,:)=t;
end
%%
burnIn=400;
figure
loglog(thetaStore(burnIn:end))
figure
plot(lambdaStore(burnIn:end,:))
figure
loglog(tStore(burnIn:end,:))

%%
function [n]=ndisasters(tau, t)
    n=histcounts(tau,t);
end

function [theta, lambda, t]=GibbsStep(~, lambda, t, params)
    theta=thetaConditional(lambda, params);
    lambda=lambdaConditional(theta, t, params);
    t=tMetropolisHastings(lambda, t, params);
end

function [theta]=thetaConditional(lambda, params)
    theta=gamrnd(4, 1/(sum(lambda) + params.vartheta));
end

function [lambda]=lambdaConditional(theta, t, params)
    lambda=zeros(params.d, 1);
    for i=1:length(lambda)
        lambda(i)=gamrnd(params.ndisasters(i)+2, 1/(t(i+1)-t(i)+theta));
    end
end

function [t]=tMetropolisHastings(lambda, t, params)
    if length(t) > 2
        tcand=t;
        epsilon=betarnd(params.rho, params.rho, 1, params.d-1);
%         size(tcand(2:end-1))
%         size(t(1:end-2))
%         size(epsilon)
%         size(t(3:end))
%         size(t(1:end-2))
        tcand(2:end-1)=t(1:end-2) + epsilon.*(t(3:end) - t(1:end-2));

        U=rand(1,params.d);
        for i=2:(params.d)
            alpha=exp(lambda(i)*(tcand(i) - t(i))-lambda(i-1)*(tcand(i) - t(i)))*(t(i+1)-tcand(i))/(tcand(i)-t(i-1));
            if U(i-1)<=alpha
                t(i)=tcand(i);
            end
        end
        
    end
end