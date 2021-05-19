%% 
tau=load("coal-mine.csv");
obs = load("mixture-observations.csv");


%%
% close all
% Gibbs sampling + MH
rng(42)



params=struct();
params.d=4;
params.vartheta=1;
params.tau = tau;

% t=linspace(1851,1963, params.d+1);
c=1:(params.d-1);
tinit=ones(1,params.d-1)*1851;
tinit=tinit+c;
tinit=[1851 tinit 1963];

params.rho=1;

% t=randn(1,d);
% t=[1851 t 1963];

%%
N = 10;
figure
theta = EM(obs, N, 0.1882);
set(gca, 'YScale','log')
plot(1:N,theta, '-*')

%%
lambda=ones(1, params.d);
theta=(4)/(sum(lambda) + params.vartheta);

K=10000;
thetaStore=zeros(K,1);  
lambdaStore=zeros(K,params.d);
tStore=zeros(K+1,params.d+1);
tStore(1,:)=tinit;
acceptedStore=zeros(K,1);
for i=1:K
    [theta,lambda,t,accepted]=GibbsStep(theta, lambda, t, params);
    thetaStore(i)=theta;
    lambdaStore(i,:)=lambda;
    tStore(i+1,:)=t;
    acceptedStore(i)=accepted;
end
mean(thetaStore)
rate=sum(acceptedStore)/K
%%
plot(tStore(1:100,:))
%%
plot(cumsum(acceptedStore)./(1:K)')
%%
burnIn= 0.3 * K;
figure
histogram(thetaStore(burnIn:end))
title("Theta")
% figure
% histogram(lambdaStore(burnIn:end,1))
% title("Lambdas")
% figure
% plot(tStore(burnIn:end,:))
% title("Breakpoints")

%%
subplot(2,2,2)
hist(params.tau, 30)
title("Location of breakpoints")
xlabel("Year")
ylabel("Number of disasters")
for i = 2:params.d
  point = mean(tStore(burnIn:end, i));
  line([point, point], ylim, 'LineWidth', 2, 'Color', 'r');
end

%% Breakpoints plotting
figure
hold on
for i=2:(params.d)
    histogram(tStore(burnIn:end,i),20)
end
title("Breakpoints")
xlabel("Year")
ylabel("Number of breakpoints in bin")

%% Intensities lambda
figure
hold on
for i=1:(params.d)
    histogram(lambdaStore(burnIn:end,i),20)
end
title("Intensities \lambda_i")
xlabel("Year")
ylabel("Number in bin")
%%
function [theta] = EM(data, N, theta0)
    theta = zeros(N, 1);
    theta(1) = theta0;
    
    for i = 2:N
      probabilities = zeros(1, length(data));
      for x = 0:1
        probabilities = probabilities + log(normpdf(data, x, x + 1))+log((1 - theta(i-1)))*(1 - x)+log(theta(i-1))*x;
      end
      theta(i) = mean(probabilities);
    end
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

function [theta, lambda, t, accepted]=GibbsStep(~, lambda, t, params)
    theta=thetaConditional(lambda, params);
    lambda=lambdaConditional(theta, t, params);
    [t,accepted]=tMetropolisHastings(lambda, t, params);
end

function [theta]=thetaConditional(lambda, params)
    theta= gamrnd(2*params.d+2, 1/(sum(lambda) + params.vartheta));
end

function [lambda]=lambdaConditional(theta, t, params)
    lambda=zeros(params.d, 1);
    nDisasters = ndisasters(params.tau, t);
    for i=1:length(lambda)
        lambda(i)= gamrnd(nDisasters(i)+2, 1/(t(i+1)-t(i) + theta));
    end
end

function [p] = tConditional(t, lambdas, params)
    nDisasters = ndisasters(params.tau, t);
    T = diff(t);
    p = exp(-T * lambdas) * prod(lambdas' .^ nDisasters) * prod(T);
end

function [t,accepted]=tMetropolisHastings(lambda, t, params)
    if length(t) > 2
        tcand=t;
        epsilon=betarnd(params.rho, params.rho, 1, params.d-1);
%         size(tcand(2:end-1))
%         size(t(1:end-2))
%         size(epsilon)
%         size(t(3:end))
%         size(t(1:end-2))
        tcand(2:end-1)=t(1:end-2) + epsilon.*(t(3:end) - t(1:end-2));
        accepted=0;
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
                accepted=accepted+1;
            end
        end
        accepted=accepted/(params.d-1);
    end
end