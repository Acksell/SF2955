%% 
tau=load("coal-mine.csv");

%%
close all
% Gibbs sampling + MH
rng(42)



params=struct();
params.d=3;
params.vartheta=2;
params.tau = tau;

t=linspace(1851,1963, params.d+1);
params.rho=1;

% t=randn(1,d);
% t=[1851 t 1963];


lambda=ones(1, params.d);
theta=(4)/(sum(lambda) + params.vartheta);

K=1000;
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
burnIn= 0.3 * K;
figure
plot(thetaStore(burnIn:end))
title("Theta")
figure
plot(lambdaStore(burnIn:end,1))
title("Lambdas")
figure
plot(tStore(burnIn:end,:))
title("Breakpoints")
figure
hold on
hist(params.tau, 30)
for i = 2:params.d
  line([t(i), t(i)], ylim, 'LineWidth', 2, 'Color', 'r');
end

%%
function [n]=ndisasters(tau, t)
%   n = zeros(1, length(t) - 1);
%   prevIdx = 1;
%   for i = 2:length(t)
%     newIdx = find(tau > t(i), 1);
%     if isempty(newIdx)
%       newIdx = length(t);
%     end
%     n(i - 1) = newIdx - prevIdx;
%     prevIdx = newIdx;
%   end
  n=histcounts(tau,t);
end

function [theta, lambda, t]=GibbsStep(~, lambda, t, params)
    theta=thetaConditional(lambda, params);
    lambda=lambdaConditional(theta, t, params);
    t=tMetropolisHastings(lambda, t, params);
end

function [theta]=thetaConditional(lambda, params)
    theta= 4/(sum(lambda) + params.vartheta);%gamrnd(4, 1/(sum(lambda) + params.vartheta));
end

function [lambda]=lambdaConditional(theta, t, params)
    lambda=zeros(params.d, 1);
    nDisasters = ndisasters(params.tau, t);
    for i=1:length(lambda)
        lambda(i)= (nDisasters(i) + 2) / (t(i + 1) - t(i) + theta);%gamrnd(nDisasters(i)+2, 1/(t(i+1)-t(i) + theta));
    end
end

function [p] = tConditional(t, lambdas, params)
    nDisasters = ndisasters(params.tau, t);
    T = diff(t);
    p = exp(-T .* lambdas') * prod(lambdas' .^ nDisasters) * prod(T);
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

        U = rand(1,params.d);
        for i=2:(params.d)
            if tcand(i) > tcand(i + 1)
              alpha = 0;
            else
              %alpha=exp(lambda(i)*(tcand(i) - t(i))-lambda(i-1)*(tcand(i) - t(i)))*(t(i+1)-tcand(i))/(tcand(i)-t(i-1));
              tmp = t;
              tmp(i) = tcand(i);
              alpha = tConditional(tmp, lambda, params) / tConditional(t, lambda, params);
            end
                         
            if U(i-1)<=alpha
                t(i)=tcand(i);
            end
        end
        
    end
end