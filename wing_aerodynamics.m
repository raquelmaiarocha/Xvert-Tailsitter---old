function [Fa, Ma] = wing_aerodynamics(xvert,Vb, Wb, Vwb, rho, dfi)
    
    % Realocate variables
    [Si, bi, ci, cfi] = matsplit(xvert.wing,1);
    riG = xvert.wing(:,[5:7]);

    Fi = zeros(11,3);
    Mi = zeros(11,3);
    Mi_ac = zeros(11,3);
    
    % Horizontal Wing Segments
    for i=1:9
        [ui, vi, wi] = airspeed_vector(Vb, Wb, [riG(i,:)]', Vwb);
        [Vi, aai, bbi] = airspeed(Vb, Wb, [riG(i,:)]', Vwb);
        [CLwi,CDwi,Cmwi] = aerodynamics(aai,dfi);

        Fi(i,1) = 0.5*rho*bi(i)*ci(i) * (ui^2+wi^2) * (CLwi*sin(aai)-CDwi*cos(aai));
        Fi(i,3) = 0.5*rho*bi(i)*ci(i) * (ui^2+wi^2) * (-CLwi*cos(aai)-CDwi*sin(aai));

        Mi(i,2) = 0.5*rho*bi(i)*ci(i)^2 * (ui^2+wi^2) * Cmwi;
        
        Mi_ac_aux = Mi(i,:)' + cross(riG(i,:)',Fi(i,:)');
        [Mi_ac(i,1), Mi_ac(i,2), Mi_ac(i,3)] = matsplit(Mi_ac_aux);

    end
    
    % Vertical Wing Segments
    for i=10:11
        [ui, vi, wi] = airspeed_vector(Vb, Wb, riG(i,:), Vwb);
        [Vi, aai, bbi] = airspeed(Vb, Wb, riG(i,:), Vwb);
        [CLwi,CDwi,Cmwi] = aerodynamics(aai,dfi);

        Fi(i,1) = 0.5*rho*bi(i)*ci(i) * (ui^2+vi^2) * (CLwi*sin(bbi)-CDwi*cos(bbi));
        Fi(i,2) = 0.5*rho*bi(i)*ci(i) * (ui^2+vi^2) * (-CLwi*cos(bbi)-CDwi*sin(bbi));

        Mi(i,3) = 0.5*rho*bi(i)*ci(i)^2 * (ui^2+vi^2) * Cmwi;

        Mi_ac_aux = Mi(i,:)' + cross(riG(i,:)',Fi(i,:)');
        [Mi_ac(i,1), Mi_ac(i,2), Mi_ac(i,3)] = matsplit(Mi_ac_aux);

    end
    
    Fa = sum(Fi)'; % Sum of Matrix Rows
    Ma = sum(Mi_ac)'; % Sum of Matrix Rows
end