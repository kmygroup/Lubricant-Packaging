function D1=derHdX_pipe_flushing(x,z)

load Q
mu1=130.26;
mu2=305.39;
a=0.019;
l=6.375;


dHdx(1)=z(1)*(-Q/(a*l))+z(2)*(Q/(a*l))+z(3)*((9*(mu1^(2/3))*(x(1)^2)*(Q/(a*l))*((mu2^(1/3))-(mu1^(1/3))))+(6*(2*(x(1))-(3*(x(1)^2))*(Q/(a*l))*(((mu1^(1/3))*(mu2^(2/3)))-((mu1^(2/3))*(mu2^(1/3))))))+(3*(mu2^(2/3))*((mu2^(1/3))-(mu1^(1/3)))*((Q/(a*l)))*(1-(4*x(1))+(3*(x(1)^2)))));

dHdx(2)=z(1)*(Q/(a*l))-z(2)*(Q/(a*l))+z(3)*((-9*(mu1^(2/3))*((1-x(2))^2)*(Q/(a*l))*((mu2^(1/3))-(mu1^(1/3))))+(6*(1-(4*x(2))+(3*(x(2)^2)))*(Q/(a*l))*(((mu1^(1/3))*(mu2^(2/3)))-((mu1^(2/3))*(mu2^(1/3)))))+(3*(mu2^(2/3))*((mu2^(1/3))-(mu1^(1/3)))*(Q/(a*l))*(2*x(2)-(3*(x(2)^2)))));

dHdx(3)=0;




D1=[dHdx(1) dHdx(2) dHdx(3)]';