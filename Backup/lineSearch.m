function [alpha, gnew] = lineSearch( f, xk, dk, gk )

% In : f  ... objective function (handle)
%      xk ... current point
%      dk ... chosen direction of descent
%      gk ... gradient of f in xk
%
% Out: alpha ... a parameter satisfying (W1) and (W2)
%      gnew  ... gradient of f in xk + alpha*dk
%
% parameters: c1 = 1e-4, c2 = 0.99, alpha0 = 0, alpha1 = 1, 
%             zoom takes alpha_i = (alpha_lo + alpha_hi)/2
%
    alpha0=0;
    alpha1=1;
    c1 = 1e-4;
    c2 = 0.99;
    fk = f(xk);
    slope0 = dot(gk,dk);
    cond =  true;
    while cond
      alpha_i = 0.5*(alpha0 + alpha1);
      f1 = f(xk + alpha_i*dk); 
      if (f1 > fk + alpha_i*c1*slope0) || (f1 >= f(xk + alpha0*dk))
        alpha1 = alpha_i;
        continue;
      else
        gnew = apGrad(f,xk + alpha_i*dk);
        slope = dot(gnew,dk); 
        if abs(slope) <= -c2*slope0
          alpha=alpha_i;
          cond=false;
          break;
        else
          if (alpha1-alpha0)*slope >= 0
            alpha1=alpha_i;
          else
            alpha0 = alpha_i;
          end
        end
      end
          
      
      
    end

end