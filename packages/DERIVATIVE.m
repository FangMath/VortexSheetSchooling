function [derivative]=DERIVATIVE(f3,f2,f1,f0,dt,o)
% approximate the first and second order derivatives
% Input:    f3,f2,f1,f0
%           dt
%           o: o==11 -> first derivative at time t_1, known derivative at
%           t0
%              o==10 -> first derivative at time t_k, k>=2
%              o==21 -> second derivative at time t_1, first order
%              o==22 -> second derivative at time t_2
%              o==20 -> second derivative at time t_k, k>=3
%              o==111-> first derivative at time t_1, first order
% (input f=0 when f useless, input dot{f} at position f0 when needed)

if o==11
    derivative=(2*f3-2*f2)/dt-f0;   %f1=0 useless, f0=dot{f}(t2)
else if o==10
        derivative=(1.5*f3-2*f2+0.5*f1)/dt; %f0=0 useless
    else if o==21
            derivative=(2*f3-2*f2)/(dt^2)-(2*f0)/dt; %f1=0 useless, f0=dot{f}(t2) 7.13??
        else if o==22
                derivative=(2.5*f3-8*f2+5.5*f1)/(dt^2)+(3*f0)/dt; %f0=dot{f}(t1)
            else if o==20
                    derivative=(2*f3-5*f2+4*f1-f0)/(dt^2);
                else if o==111 %f1=f0=0 useless
                        derivative=(f3-f2)/dt;
                    else
                        sprintf('wrong in DERIVATIVE.m');
                    end
                end
            end
        end
    end
end
end


