clc
clear;close all;



b =(input('Enter the width of the well(nm):'))/2
V =input('Enter potential of the well(eV):')
E =input('Enter the electron energy for finite well :')
m=input('Enter the particle mass in terms of the electron mass (m_e=1):')
n =input('Enter the quantum number:')
finite_well_wave_function(b,V,E,m,n)
infinite_well_wave_function(b,m,n)


function finite_well_wave_function(b0,V0,E0,m0,n0)
 
    miu=sqrt((2*m0*0.5*10^6*(E0+V0))/200^2)
    k=sqrt((2*m0*0.5*10^6*(-E0))/200^2)
    if mod(n0,2)==1
      B=exp(k*b0)*cos(miu*b0)*sqrt(k/(1+k*b0))
      D=sqrt(k/(1+k*b0))
      f1=@(x)B*exp(k*x)
      f3=@(x)B*exp(-k*x)
      f2=@(x)D*cos(miu*x)
    else
       
      C=sqrt(k/(1+k*b0))
      B=-exp(k*b0)*sin(miu*b0)*sqrt(k/(1+k*b0))
      f1=@(x)B*exp(k*x)
      f3=@(x)-B*exp(-k*x)
      f2=@(x)C*sin(miu*x)
    end 
    subplot(2,1,1);
    fplot(f1,[-1.5*b0,-b0]);
    title(['Wave function for finite well , n=',num2str(n0),'  E=',num2str(V0+E0)]);
    xlabel('x(nm)') ;
    grid on;hold on;fplot(f2,[-b0,b0]);hold on;fplot(f3,[b0,1.5*b0]);hold on;
    return;
    
end



function infinite_well_wave_function(b0,m0,n0)
    A=sqrt(1/b0)
    f1=@(x)0
    f3=@(x)0
    if mod(n0,2)==1
      f2=@(x)A*cos(n0*pi*x/(2*b0))
    else
      f2=@(x)A*sin(n0*pi*x/(2*b0))
    end 
    E=(200)^2*pi^2*n0^2/(2*0.5*10^6*m0*4*b0^2)
    subplot(2,1,2);
    fplot(f1,[-1.5*b0,-b0]);
    title(['Wave function for Infinite well , n=',num2str(n0),'  E=',num2str(E)]);
    xlabel('x(nm)') ;
    grid on;hold on;fplot(f2,[-b0,b0]);hold on;fplot(f3,[b0,1.5*b0])
    return;
    
end
