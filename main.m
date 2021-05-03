%%
% rng()
pi=load('stations.mat');
pi = pi.pos_vec;
alpha = 0.6;
sigmaW=0.5;
deltat=0.5; % seconds
v=90;
eta=3;
zeta=1.5;

M=500;
Zvalues = [[0 0]' [3.5 0]' [0 3.5]' [0 -3.5]' [-3.5 0]'];


X=zeros(6,M);
variances = sqrt([500 5 5 200 5 5]');
X(:,1)=randn(6,1).*variances;

s = RandStream('mlfg6331_64','Seed',1337);

P=(ones(5,5)+15*eye(5,5))/20;
Z=zeros(M,5);
Z(1,randi(5,1))=1;
for m=2:M
    Z_prob = P*Z(m-1,:)';
    R = randsample(s,1:5,1,true,Z_prob);
    Z(m,R)=1;
end

W=sigmaW*randn(2,M);

Phi = [1 deltat deltat^2/2; 0 1 deltat; 0 0 alpha];
Phi = kron(eye(2),Phi);

Psi_z = [deltat^2/2 deltat 0]';
Psi_z = kron(eye(2), Psi_z);

Psi_w = [deltat^2/2 deltat 1]';
Psi_w = kron(eye(2), Psi_w);


for m=2:M
    X(:,m)=Phi*X(:, m-1) + Psi_z*(Zvalues*Z(m-1,:)') + Psi_w*W(:,m);
end


nbase_stations = 6;
V = zeta*randn(nbase_stations, M);
Xpos = [X(1,:); X(4,:)];

Y=zeros(nbase_stations,M);
for k=1:nbase_stations
    Y(k,:) = v - 10*eta*log10( vecnorm( Xpos - pi(:,k) )) + V(k,:);
end


%% Plot trajectory
figure
plot(Xpos(1,:),Xpos(2,:),'-k')
hold on
plot(pi(1,:),pi(2,:),'*r')

%%
Y=load('RSSI-measurements.mat');
Y=Y.Y;
M=length(Y);

%% Problem

q0=getQ0(variances);
q=getQ(Phi, Psi_z, sigmaW);
p=getP(eta, pi, zeta, v);

N=10000;
% X0=randn(6,1).*variances;
% q0(X0)
[X_sis, tau_sis, w_sis] = SIS(N,M,p,Phi,Psi_z,Psi_w,Y);


%% Plotting

% Plot trajectory
figure
% plot(Xpos(1,:),Xpos(2,:),'-k')
hold on
plot(tau_sis(:,1),tau_sis(:,2),'-b')
plot(pi(1,:),pi(2,:),'*r')
% figure
% hist(w)


cmap = jet(256);
vrescaled = rescale(w_sis(1,:), 1, 256); % Nifty trick!
numValues = length(w_sis(1,:));
markerColors = zeros(numValues, 3);
% Now assign marker colors according to the value of the data.
for k = 1 : numValues
    row = round(vrescaled(k));
    markerColors(k, :) = cmap(row, :);
end
% Create the scatter plot.
figure
scatter(X_sis(1,1,:), X_sis(1,4,:), [], markerColors);
grid on;
colorbar()
%% histogram

mlist=[1, 10, 100, M];
for m=1:length(mlist)
    subplot(2,2,m);
    histogram(log(w_sis(mlist(m),:)),100);
    title(sprintf('Log-weights at m=%d',mlist(m)))
    ylabel("frequency")
    xlabel("log-weights")
end

% ESS(w(1,:))
%%

q0=getQ0(variances);
q=getQ(Phi, Psi_z, sigmaW);
p=getP(eta, pi, zeta, v);

N=10000;
% X0=randn(6,1).*variances;
% q0(X0)
[X_sisr, tau_sisr, w_sisr, Zexp] = SISR(N,M,p,Phi,Psi_z,Psi_w,Y);

%% Plotting

% Plot trajectory
figure
% plot(Xpos(1,:),Xpos(2,:),'-k')
hold on
plot(tau_sisr(:,1),tau_sisr(:,2),'-b')
plot(pi(1,:),pi(2,:),'*r')
% figure
% hist(w)

cmap = jet(256);
vrescaled = rescale(w_sisr(1,:), 1, 256); % Nifty trick!
numValues = length(w_sisr(1,:));
markerColors = zeros(numValues, 3);
% Now assign marker colors according to the value of the data.
for k = 1 : numValues
    row = round(vrescaled(k));
    markerColors(k, :) = cmap(row, :);
end
% Create the scatter plot.
figure
scatter(X_sisr(1,1,:), X_sisr(1,4,:), [], markerColors);
grid on;
colorbar()

%%
[maxZ, maxZi]=max(Zexp,[],2);
[realZ,realmaxZ]=max(Z,[],2);

correct = realmaxZ==maxZi;

sum(correct)/M

%%
Y=load('RSSI-measurements-unknown-sigma.mat');
Y=Y.Y(:,1:end);
M=length(Y);

%%

N=1000;
gridsize=100;
Omega=zeros(gridsize,1);
gridz = linspace(1.8,2.7,gridsize);

for z=1:length(gridz)
    p=getP(eta, pi, gridz(z), v);
    [X_sisr, tau_sisr, w_sisr, Zexp] = SISR(N,M,p,Phi,Psi_z,Psi_w,Y);
    
    wExp_sisr=exp(w_sisr);
    Omega(z)=sum(log(sum(wExp_sisr, 2)))-(M+1)*log(N);
end
%%
plot(gridz, Omega)
title("Log-likelihood depending on \varsigma")
xlabel("\varsigma")
ylabel("Log-likelihood")
[OmegaMax, i]=max(Omega);
zetaOpt = gridz(i)

%%
function [X, tau, w]=SIS(N,M,p,Phi,Psi_z,Psi_w,Y)
    Zvalues = [[0 0]' [3.5 0]' [0 3.5]' [0 -3.5]' [-3.5 0]'];
    P=(ones(5,5)+15*eye(5,5))/20;
    w=ones(M,N);
    g=zeros(M,N);
    Z=zeros(M,5,N);
    X=zeros(M,6,N);
    
    variances = sqrt([500 5 5 200 5 5]');
    for i=1:N
        X(1,:,i)=randn(6,1).*variances;
        Z(1,randi(5,1),i)=1;
    end
    w(1,:)=p([squeeze(X(1,1,:)), squeeze(X(1,4,:))], Y(:,1));
    
    sigmaW=0.5;
    W=sigmaW*randn(M,2,N); 
    s = RandStream('mlfg6331_64', 'Seed', 43);
    pop = 1:5;
    tau=zeros(M,2);
    for m=2:M
       X(m,:,:)=Phi*squeeze(X(m-1,:,:)) + Psi_z*(Zvalues*squeeze(Z(m-1,:,:))) + Psi_w*squeeze(W(m,:,:));
       Z_prob = P*squeeze(Z(m-1,:,:));
       for n=1:N
           R = randsample(s,pop,1,true,Z_prob(:,n));
           Z(m,R,n)=1;
       end
       w(m,:)=w(m-1,:) + p([squeeze(X(m,1,:)) squeeze(X(m,4,:))],Y(:,m));
    end
    w=normalize(w);
    for m=1:M
        tau1 = sum(w(m,:).*squeeze(X(m,1,:))');
        tau2 = sum(w(m,:).*squeeze(X(m,4,:))');
        tau(m,:)=[tau1; tau2];
    end
end

function [X, tau, w, Zexp]=SISR(N,M,p,Phi,Psi_z,Psi_w,Y)
    Zvalues = [[0 0]' [3.5 0]' [0 3.5]' [0 -3.5]' [-3.5 0]'];
    P=(ones(5,5)+15*eye(5,5))/20;
    w=ones(M,N);
    g=zeros(M,N);
    Z=zeros(M,5,N);
    X=zeros(M,6,N);
    
    variances = sqrt([500 5 5 200 5 5]');
    for i=1:N
        X(1,:,i)=randn(6,1).*variances;
        Z(1,randi(5,1),i)=1;
    end
    w(1,:)=p([squeeze(X(1,1,:)), squeeze(X(1,4,:))], Y(:,1));
    
    sigmaW=0.5;
    W=sigmaW*randn(M,2,N); 
    s = RandStream('mlfg6331_64', 'Seed', 43);
    pop = 1:5;
    tau=zeros(M,2);

    for m=2:M
       X(m,:,:)=Phi*squeeze(X(m-1,:,:)) + Psi_z*(Zvalues*squeeze(Z(m-1,:,:))) + Psi_w*squeeze(W(m,:,:));
       Z_prob = P*squeeze(Z(m-1,:,:));

       for n=1:N
           R = randsample(s,pop,1,true,Z_prob(:,n));
           Z(m,R,n)=1;
       end
       w(m,:)=p([squeeze(X(m,1,:)) squeeze(X(m,4,:))],Y(:,m));

       ind=randsample(N,N,true,normalize(w(m,:)));
       X(m,:,:)=X(m,:,ind);
       Z(m,:,:)=Z(m,:,ind);
    end
    Zexp = zeros(5,M);

    for m=1:M
        tau1 = sum(squeeze(X(m,1,:))/N);
        tau2 = sum(squeeze(X(m,4,:))/N);
        Zexp(:,m) = sum(squeeze(Z(m,:,:))/N, 2);
        tau(m,:)=[tau1; tau2];
    end
    Zexp=Zexp';
end

function [W]=normalize(wLog)
%     wLog(1,:)
    wExp=exp(wLog-max(wLog, [], 2));
    W=wExp./sum(wExp, 2);
end

function [pfunc]=getP(eta, pi, zeta, v)
    function [logprob]=p(Xpos, Y)
        logprob=zeros(1,length(Xpos));
        for k=1:6
            logprob = logprob + log(normpdf(Y(k), v - 10*eta*log10( vecnorm( Xpos' - pi(:,k) )), zeta));
        end
    end
    pfunc=@p;
end

function [ess]=ESS(w)
    ess=1/sum(w.^2);
end

function [qfunc]=getQ(Phi, Psi_z, sigmaW)
    function [prob]=q(X_new, X_old, Z_old)
        prob=normpdf(X_new, Phi*X_old + Psi_z*Z_old, sigmaW);
    end
    qfunc=@q;
end

function [q0func]=getQ0(variances)
    function [prob]=q0(X_0)
        prob=normpdf(X_0, 0, variances);
    end
    q0func=@q0;
end










