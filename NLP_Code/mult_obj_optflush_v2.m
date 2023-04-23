function Q = mult_obj_optflush_v2(q)

%Parameter values (ALMO 532)
mu1=130.26;
mu2=305.39;
a=0.019; %cross sectional area
l=6.375; %pipe length


Tinitial=0;
Tfinal=290;
h=1; 
N = (Tfinal - Tinitial)/h;


t = [Tinitial:h:Tfinal]';

sum1 =0;


% Model equations:
for i=1:1:length(t)   

%----------- integration of equations --------------%
     
       
x1(i) = (exp(-q(i)*t(i)/(a*l)));
x2(i) = 1-(exp(-q(i)*t(i)/(a*l)));
x3(i) = [(x1(i)*(mu1^(1/3)))+(x2(i)*(mu2^(1/3)))]^3;

obj(i)=((x3(i)-(mu2))/mu2)^2; %x3 is the mixture viscosity

sum1=sum1+obj(i);  
end 

%objective function:

Q=sum1;

end

