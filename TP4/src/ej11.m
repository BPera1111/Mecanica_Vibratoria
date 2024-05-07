function ej11 ;clc; close all;

    tf=0.1;
    zeta = 0.01;
    beta = beta_0(tf,zeta);
    disp("beta: " + beta)

    disp("T_f: "+trans(0,beta))



end


function tf =trans(zeta,beta)

    tf = sqrt((1+(2*zeta*beta)^2)/((1-beta^2)^2+(2*zeta*beta)^2));
    
end

function beta = beta_0(tf,zeta)
    
        beta = 0;
        while true
            if trans(zeta,beta) < tf
                break;
            end
            beta = beta + 0.01;
        end
end
