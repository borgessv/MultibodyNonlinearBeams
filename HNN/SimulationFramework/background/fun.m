function f = fun(q,p,M,I,DoF)

JqCM = complexstep(@(q) element_positionCM(q,DoF),q);

f = p.'*((JqCM.'*M*JqCM + I)\p);

end