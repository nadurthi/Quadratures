n=1;
mu=zeros(n,1);
P=eye(n);

%% CKF moments above 2 
[x,w]=cubature_KF_points(mu,P);
m_ckf=cal_moments_wrt_pts(x,w,4,P);
%% UT moments above 2 and tuning of kappa
global kappa
kappa=2;
[x,w]=UT_sigmapoints(mu,P,2);
m_ut=cal_moments_wrt_pts(x,w,4,P);
[m_ckf(:,end),m_ut(:,end)]
