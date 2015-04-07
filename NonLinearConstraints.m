function [ nonLinConIneq, nonLinConEq ] = NonLinearConstraints( paramVars, param1 )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% param=[paramVars, param1];


g=9.8; 
a=1; 
sf=3/4; 

a0=paramVars(5); 
a1=paramVars(6);
a2=paramVars(7); 
l_a=paramVars(8);
l_p=paramVars(9); 
l_d=paramVars(10); 
t_p=paramVars(11);
t_d=paramVars(12); 

Ro0=param1(3);
Ri0=param1(4);
Ro1=param1(5);
Ri1=param1(6);
Ro2=param1(7);
Ri2=param1(8);
rho=param1(9)
failureStress= param1(10); 
Mmot1=param1(11); 
Mmot2=param1(12); 
f1=param1(13) ;
f2=param1(14); 

Mgrip=GripperMass(l_a,l_p,l_d,t_p,t_d);

m1=rho*pi*(Ro1^2-Ri1^2)*a1;
m2=rho*pi*(Ro2^2-Ri2^2)*a2;



V1max=-(m1+m2+Mmot2+Mgrip)*(g+a);
M1max=-((m1/2+Mmot2)*a1+m2*(a1+a2/2)+Mgrip*(a1+a2))*(g+a); 

V2max=-(m2+Mgrip)*(a+g);
M2max=-(m2/2+Mgrip)*(a+g)*a2; 

s0=SigmaMax(M1max,Ro0,Ri0);
s1=SigmaMax(M1max,Ro1,Ri1);
s2=SigmaMax(M2max, Ro2,Ri2);


nonLinConIneqStresses=[s0 s1 s2]; %stresses 
nonLinConIneqTorques=[M1max M2max]; %torques

% waitforbuttonpress


% failureStress;
% f1;
% f2;

nonLinConIneqStresses=abs(nonLinConIneqStresses)-sf*failureStress;
nonLinConIneqTorques=abs(nonLinConIneqTorques)-[f1 f2];

nonLinConIneq= [nonLinConIneqStresses, nonLinConIneqTorques];
nonLinConEq=[]; 

%constraint 1

% waitforbuttonpress

end

   function [sigma]= SigmaMax(M, Ro, Ri)
        
        Z=0.785*(Ro^4-Ri^4)/Ro;
       sigma=M/Z; 
        
   end
    
   
    function [tau]= TauMax(V, Ro, Ri)
        
        A=pi*(Ro^2-Ri^2);
       tau=2*V/A; 
        
    end
