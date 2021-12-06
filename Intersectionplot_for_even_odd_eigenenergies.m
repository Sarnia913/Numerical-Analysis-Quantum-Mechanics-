clc
clear;close all;

V =input('Enter potential of the well(eV):')
b =(input('Enter the width of the well(nm):'))/2
m =input('Enter the particle mass in terms of the electron mass (m_e=1):')
finite_well_Energy_plot(V,b,m)



function finite_well_Energy_plot(V0,b0,m0)
    k=(2*m0*0.5*10^6*(b0^2))/(200^2)
    f1=@(x)tan(sqrt(k*(x+V0)))
    f2=@(x)sqrt(-x/(x+V0))
    f3=@(x)-sqrt((x+V0)/(-x))
    fplot(f1,[-V0,0]);
    title('Intersection Plot for even and odd eigen Energies');
    xlabel('Energy(eV)') ;
    grid on;hold on;fplot(f2,[-V0,0]);hold on;fplot(f3,[-V0,0])
    legend({'Tan function','Even Solutions','Odd Solutions'},'Location','northeast')
    return;
end


