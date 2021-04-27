pi=load('stations.mat');
pi = pi.pos_vec;
alpha = 0.6;
sigma=0.5;
deltat=0.5; % seconds
v=90;
eta=3;
zeta=1.5;

M=10000;
Zvalues = [[0 0]' [3.5 0]' [0 3.5]' [0 -3.5]' [-3.5 0]'];


X=zeros(6,M);
variances = sqrt([500 5 5 200 5 5]');
X(:,1)=randn(6,1).*variances;

Z=Zvalues(:, randi(5,M));

W=sqrt(sigma)*randn(2,M);

Phi = [1 deltat deltat^2/2; 0 1 deltat; 0 0 alpha];
Phi = kron(eye(2),Phi);

Psi_z = [deltat^2/2 deltat 0]';
Psi_z = kron(eye(2), Psi_z);

Psi_w = [deltat^2/2 deltat 1]';
Psi_w = kron(eye(2), Psi_w);


for m=2:M
    X(:,m)=Phi*X(:, m-1) + Psi_z*Z(:, m-1) + Psi_w*W(:,m);
end


nbase_stations = 6;
V = zeta*randn(nbase_stations, M);
Xpos = [X(1,:); X(4,:)];


Y=zeros(6,M);
for k=1:6
    Y(k,:) = v - 10*eta*log10( vecnorm( Xpos - pi(:,k) )) + V(k,:);
end






