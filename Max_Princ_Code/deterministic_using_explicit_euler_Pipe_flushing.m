% This file uses the explicit euler for integration of all the equations


clc
clear all
tic; 
% ========= initial information Test Case 1 ===================================
del_t=0.05;
mu1=130.26;
mu2=305.39;
a=0.019; 
l=6.375;      
t0=0;     tf=290;  h=del_t;
time=t0:h:tf; 
t_len=length(time);
          
%============initial conditions for state and theta==================
I0 = [1	0 130.26]';

Thi=[0 0 0]';

%============final conditions for adjoint and phi==================

Fz = [0 0 -1]';

Phif = [0 0 0]';

 %Preallocation of space:

         x1=ones(1,t_len); 
         x2=ones(1,t_len);
         x3=ones(1,t_len);
        
        
        alfa1 = ones(1,t_len);
        alfa2 = ones(1,t_len);
        alfa3 = ones(1,t_len);
        
        
% Preallocation of spaces for adjoints:

        z1= ones(1,t_len);
        z2= ones(1,t_len);
        z3= ones(1,t_len); 
        
        
        zdot1=ones(1,t_len);
        
        zdot2=ones(1,t_len);

        zdot3=ones(1,t_len);
    

dtheta1 = ones(1,t_len);

dtheta2 = ones(1,t_len);

dtheta3 = ones(1,t_len);


 theta1 = ones(1,t_len);
        
 theta2 = ones(1,t_len);
        
 theta3 = ones(1,t_len);

 
        phi1 = ones(1,t_len);
        
        phi2 = ones(1,t_len);
        
        phi3 = ones(1,t_len);
        

dphi1= ones(1,t_len);
dphi2 = ones(1,t_len);
dphi3 = ones(1,t_len);



%======================initial guess for Q(flowrate)==================

Qflow = ones(1,t_len);
Qflownew = ones(1,t_len);
Qflow = 0.00303*Qflownew;  


% Optimal flow

tolerance=1e-4*ones(1,t_len);
dHdQ = ones(1,t_len);
der1= ones(1,t_len);
der2 = ones(1,t_len);

H = ones(1,t_len);
numiter = 1; 
p=1;


% while (abs(dHdT(:)) > tolerance(:))

for j=1:10 
%================== solution of state equations========================
        x1 = 1;
        x2 = 0;
        x3 = 130.26;
        
   
  for i=1:(t_len)
        Q=Qflow(i);
        save Q

       

        alfa1(i) = -((x1(i))*(Q/(a*l))); 
        
        alfa2(i) = (x1(i))*(Q/(a*l));

        alfa3(i)= (((x1(i))*(Q/(a*l)))*((3)*x3(i)^(2/3)))*((mu2^(1/3))-(mu1^(1/3)));
        
        
        x1(i+1) = x1(i) + alfa1(i)*del_t;
        
        
        x2(i+1) = x2(i) + alfa2(i)*del_t;
        
 
        x3(i+1) = x3(i) + alfa3(i)*del_t;
         
    
      
  end
     state = [x1' x2' x3'];

    saveallstate(:,p)=x1';
    saveallstate(:,p+1)=x2';
    saveallstate(:,p+2)= x3';
    p=p+4; 

    figure(1)
    plot(time,Qflow); hold all
    ylabel ('Q_flow')
    xlabel ('Time')


%%Part II: the backward integration for adjoints:
%Final condition at time t = 290
Fz = [0 0 -1]';
%====================================================================
% Solution of adjoint equations: z equations
        zx1 = 0;
        zx2 = 0;
        zx3 = -1; 
        
    
    for i=(t_len):-1:2
        Q=Qflow(i);
        save Q
        z1(t_len)=zx1; 
        z2(t_len)= zx2;
        z3(t_len)= zx3; 

        zdot1(i)=(z1(i)*(Q/(a*l)))-(z2(i)*(Q/(a*l)))-(z3(i)*((9*(mu1^(2/3))*(x1(i)^2)*(Q/(a*l))*((mu2^(1/3))-(mu1^(1/3))))+((12*x1(i)-18*(x1(i)^2))*(Q/(a*l))*(((mu1^(1/3))*(mu2^(2/3)))-((mu1^(2/3))*(mu2^(1/3)))))+((3*(Q/(a*l)))*(1-4*(x1(i))+3*(x1(i))^2))-((mu2)-((mu1^(1/3))*(mu2^(2/3))))));

        zdot2(i)=((-z1(i))*(Q/(a*l)))+(z2(i)*(Q/(a*l)))-(z3(i)*((-9*(mu1^(2/3))*((1-x2(i))^2)*(Q/(a*l))*((mu2^(1/3))-(mu1^(1/3))))+(6*((1-(4*x2(i))+(3*x2(i))^(2)))*(Q/(a*l))*(((mu1^(1/3))*(mu2^(2/3)))-((mu1^(2/3))*(mu2^(1/3)))))+(3*(mu2^(2/3))*(2*(x2(i))-3*(x2(i)^(2)))*(Q/(a*l))*((mu2^(1/3))-(mu1^(1/3))))));

        zdot3(i)=0;

        
       
        z1(i-1) = z1(i) - zdot1(i)*del_t; 
        
        z2(i-1) = z2(i) - zdot2(i)*del_t; 
        
        z3(i-1) = z3(i) - zdot3(i)*del_t; 
        
       
     end
     adjz = [z1' z2' z3' ];
        
%=======Solving for the Hamiltonian================================

%================== solution of theta equations========================
        th1 = 0;
        th2 = 0;
        th3 = 0; 
        

  for i=1:(t_len)
        Q=Qflow(i);
        save Q
    
        theta1(1)= th1;
        theta2(1)= th2;
        theta3(1)= th3; 
        

dtheta1(i) = -x1(i)/(a*l);

dtheta2(i)=x1(i)/(a*l);

dtheta3(i) = (3*(mu1^(2/3))*(x1(i)^3)*(1/(a*l))*((mu2^(1/3))-(mu1^(1/3))))+(6*(x1(i)^2)*x2(i)*(1/(a*l))*(((mu1^(1/3))*(mu2^(2/3)))-((mu1^(2/3))*(mu2^(1/3)))))+(3*(mu2^(2/3))*x1(i)*(x2(i)^2)*(1/(a*l))*((mu2^(1/3))-(mu1^(1/3))));


 theta1(i+1) = theta1(i) + dtheta1(i)*del_t; 
        
 theta2(i+1) = theta2(i) + dtheta2(i)*del_t; 
        
 theta3(i+1) = theta3(i) + dtheta3(i)*del_t; 
 
 
 J(i)= ((x3(i)-mu2)^2);

 Objective=J(i);
 
  end
  
  sol_Theta= [theta1' theta2' theta3'];

% ====================solution of phi equations============================
   sol_Phi=zeros(t_len,3);       
   phi=Phif';
        phix1 = 0;
        phix2 = 0;
        phix3 = 0; 
        
    
    for i=(t_len):-1:2
        Q=Qflow(i);
        save Q
        phi1(t_len)=phix1;
        phi2(t_len)= phix2;
        phi3(t_len)= phix3; 
        

dphi1(i)= -z1(i)*(x1(i)/(a*l));

dphi2(i) = z2(i)*(x1(i)/(a*l));

dphi3(i) = z3(i)*((3*(mu1^(2/3))*(x1(i)^3)*(1/(a*l))*((mu2^(1/3))-(mu1^(1/3))))+(6*(x1(i)^2)*(x2(i))*((1/(a*l)))*(((mu1^(1/3))*(mu2^(2/3)))-((mu1^(2/3))*(mu2^(1/3)))))+(3*(mu2^(2/3))*(x1(i))*(x2(i)^2)*((1/(a*l)))*(((mu2^(1/3))-(mu1^(1/3))))));

       
        phi1(i-1) = phi1(i) - dphi1(i)*del_t; 
        
        phi2(i-1) = phi2(i) - dphi2(i)*del_t; 
        
        phi3(i-1) = phi3(i) - dphi3(i)*del_t; 
        
        

      
     end
     sol_Phi = [phi1' phi2' phi3'];


%=====================Hamiltonian%evaluation==========
    dHdx=zeros(t_len,3);
    dHdz=zeros(t_len,3);
    
for i=1:(t_len)
    Q=Qflow(i);
    save Q
    Xs=[state(i,1) state(i,2) state(i,3)]';
    Zs=[adjz(i,1) adjz(i,2) adjz(i,3)]';
    dHdx(i,:) =feval(@derHdX_pipe_flushing,Xs,Zs); 
end

for i=1:(t_len)
    Q=Qflow(i);
    save Q
         Xs =[state(i,1) state(i,2) state(i,3) ]';
    dHdz(i,:) =feval(@derHdZ_pipe_flushing,Xs);
end

    for i=1:(t_len)
        der1(i)=(dHdx(i,1)*theta1(i)+dHdx(i,2)*theta2(i)+dHdx(i,3)*theta3(i));
        der2(i)=(dHdz(i,1)*phi1(i)+dHdz(i,2)*phi2(i)+dHdz(i,3)*phi3(i));
        dHdQ(i)= der1(i) + der2(i);
        H(i)= z1(i)*dHdz(i,1)+z2(i)*dHdz(i,2)+z3(i)*dHdz(i,3);
    end

    figure(2)
    plot(time,(H)); hold all
    ylabel ('Hamiltonian')
    xlabel ('Time')

        
    figure(3)
    plot(time,(dHdQ)); hold all
    ylabel ('dHdQ')
    xlabel ('Time')

    for k=1:t_len
        updt = 1e-6;
        Qflownew(k)=Qflow(k) + updt*dHdQ(k);
        Q_iters(k)= Qflownew(k);
        Qflow(k)=Qflownew(k);
    
    if (Qflow(k)<0)
        Qflow(k)=0;
    else if (Qflow(k)>0)
        Qflow(k)=Qflow(k);
         end
    end
       
    end
        
    for i=1:t_len
    der_H(j,i) = dHdQ(i);
    Flowrate(j,i) = Qflow(i);
    Hamiltonian(j,i) = H(i);
    end
    numiter = numiter+1;
end
    toc;

figure(4)
plot(time, Qflow)
ylabel ('Flowrate profile')
xlabel ('Time')







