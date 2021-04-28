%%
pi=load('stations.mat');
pi = pi.pos_vec;
alpha = 0.6;
sigmaW=0.5;
deltat=0.5; % seconds
v=90;
eta=3;
zeta=1.5;

M=50;
Zvalues = [[0 0]' [3.5 0]' [0 3.5]' [0 -3.5]' [-3.5 0]'];


X=zeros(6,M);
variances = sqrt([500 5 5 200 5 5]');
X(:,1)=randn(6,1).*variances;

P=(ones(5,5)+15*eye(5,5))/20;
Z=zeros(M,5);
Z(1,randi(5,1))=1;
for m=2:M
    Z_prob = P*Z(m-1,:)';
    R = randsample(s,1:5,1,true,Z_prob);
    Z(m,R)=1;
end

W=sqrt(sigmaW)*randn(2,M);

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


Y=zeros(6,M);
for k=1:6
    Y(k,:) = v - 10*eta*log10( vecnorm( Xpos - pi(:,k) )) + V(k,:);
end

% Plot trajectory
plot(Xpos(1,:),Xpos(2,:))

%% Problem
q0=getQ0(variances);
q=getQ(Phi, Psi_z, sigmaW);
p=getP(eta, pi, zeta, v);

N=10000;
% X0=randn(6,1).*variances;
% q0(X0)
[X, tau] = SIS(N,M,q0,p,q,Phi,Psi_z,Psi_w,Y);

%% Plotting

% Plot trajectory
figure
plot(Xpos(1,:),Xpos(2,:),'-k')
hold on
plot(tau(:,1),tau(:,2),'-b')

%%
function [X, tau]=SIS(N,M,q0,p,q,Phi,Psi_z,Psi_w,Y)
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
 
    size([squeeze(X(1,1,:)), squeeze(X(1,4,:))])
    size(Y)
    size(p([squeeze(X(1,1,:)), squeeze(X(1,4,:))], Y(1,:)))
    w(i,:)=p([squeeze(X(1,1,:)), squeeze(X(1,4,:))], Y(1,:));
    sigmaW=0.5;
    W=sqrt(sigmaW)*randn(M,2,N); 
    s = RandStream('mlfg6331_64');
    pop = 1:5;
    tau=zeros(M,2);
    for m=2:M
       X(m,:,:)=Phi*squeeze(X(m-1,:,:)) + Psi_z*(Zvalues*squeeze(Z(m-1,:,:))) + Psi_w*squeeze(W(m,:,:));
       Z_prob = P*squeeze(Z(m-1,:,:));
       for n=1:N
           R = randsample(s,pop,1,true,Z_prob(:,n));
           Z(m,R,n)=1;
       end
       w(m,:)=w(m-1,:).*p([squeeze(X(m,1,:)) squeeze(X(m,4,:))],Y);
       w(m,:)=w(m,:)/sum(w(m,:));
       
       tau(m,:)=[sum(w(m,:).*squeeze(X(m,1,:))'); sum(w(m,:).*squeeze(X(m,4,:))')];
    end
end


%%
function [pfunc]=getP(eta, pi, zeta, v)
    function [prob]=p(Xpos, Y)
        prob=ones(1,length(Xpos));
        for k=1:6
            prob = prob.*normpdf(Y(k), v - 10*eta*log10( vecnorm( Xpos' - pi(:,k) )), zeta);
        end
    end
    pfunc=@p;
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










