function mass = GripperMass(l_a,l_p,l_d,t_p,t_d)

global thickness density beta x0 obj_sq w1 mass success_num alpha

thickness = 3;
alpha = 0;

volume = thickness*(l_d*t_d + 3.14*((t_d)^2)/4 + l_p*t_p + 0.5*l_p*tan(alpha*(3.14/180)) + (3.14/4)*((t_p + l_p*tan(alpha*(3.14/180)))^2) + l_a*((t_p + l_p*tan(alpha*(3.14/180)))^2));

mass = 2*volume*density;

end
