function dxdt=Tyson1988_fun(t,stuff)
% parameters
global L1 L2 k c alpha lambda1 lambda2 s1 s2 s epsilon epsilon_prime D k1
global numofgrids h
% vector to matrix
rho=reshape(stuff(1:numofgrids^2,1),numofgrids,numofgrids);
beta=reshape(stuff((numofgrids^2+1):numofgrids.^2*2,1),numofgrids,numofgrids);
gamma=reshape(stuff((numofgrids.^2*2+1):numofgrids.^2*3,1),numofgrids,numofgrids);
% Laplacian
L=4*del2(gamma,h);
for i=1:1:numofgrids
    for j=1:1:numofgrids
        f1gamma(i,j)=(1+k.*gamma(i,j))./(1+gamma(i,j));
        f2gamma(i,j)=(L1+k.*L2.*c.*gamma(i,j))./(1+c.*gamma(i,j));
        Y1(i,j)=rho(i,j).*gamma(i,j)./(1+gamma(i,j));
        phi(i,j)=(lambda1+Y1(i,j).^2)./(lambda2+Y1(i,j).^2);
    end
end
% ODEs 
dxdt_1=k1.*(-f1gamma.*rho+f2gamma.*(ones(numofgrids,numofgrids)-rho));
dxdt_2=k1.*(s1.*phi-beta)./epsilon_prime;
dxdt_3=k1.*(epsilon.*L+(s2.*beta-gamma)./epsilon);
% matrix to vector
rho_vec=reshape(dxdt_1,[],1);
beta_vec=reshape(dxdt_2,[],1);
gamma_vec=reshape(dxdt_3,[],1);
dxdt=[rho_vec;beta_vec;gamma_vec];

end