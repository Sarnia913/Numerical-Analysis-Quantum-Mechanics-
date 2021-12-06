clc
clear;close all;


a =input('Enter the width of the well(nm):')
e=input('Enter the Electric Field(Vnm-1):')
N=input('Enter the number of eigen states:')

C=-1*e*a/2;

V1=@(e,x,C) (1*e*x+C);




B = zeros(N,N)
for k=1:N
    for n=1:N
        phin=@(a,n,x)(sqrt(2/a)*sin(n*pi*x/a));
        phik=@(a,k,x)(sqrt(2/a)*sin(k*pi*x/a));
        mult=@(x)V1(e,x,C).*phin(a,n,x).*phik(a,k,x);
        t=quad(mult,0,a);
        B(k,n)=t;
    end
end
V=B;  %electric field purturbation


for n=1:N
    A(n)=(200^2*pi^2*n^2)/(2*0.5*10^6*1*a^2);
end    

En=diag(A) % energy of unpurturbed state
T=En+V  % matrix of purturbed state
[Vec,D]=eig(T,'vector')  % energy and eigen vectors of purturbed states
[D, ind] = sort(D)
Vec = Vec(:, ind)
Vec=Vec.'

syms t 'real'

for n=1:N
     phi(n)=sqrt(2/a)*sin(n*pi*t*(1/a)) ;  
end

for n=1:N

  sai(n)= sum(phi.*Vec(n,:))
end  


for n=1:N
  subplot(round(N/2),2,n)
  fplot(sai(n),[0 a])
  leg=sprintf('E= %d En= %d',D(n,1),En(n,n));
  legend(leg);
  title([' Eigen vector =',num2str(n)], 'FontSize', 12);
  hold on
  grid on
  xlabel('x(nm)', 'FontSize',10);
  ylabel('wave function', 'FontSize', 10);
end  
sgtitle(['plot of wave function vs. x in purturbed infinite square well n=',num2str(N)], 'FontSize', 12);


