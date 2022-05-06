% AUXILIAR
clear
%close all
clc
%clf

%% WING SEGMENTS AERODYNAMICS

alpha = -180:1:180; % length(alpha) = 361
CL = zeros(1,length(alpha));
CD = zeros(1,length(alpha));
Cm = zeros(1,length(alpha));

for i=1:length(alpha)
    [CL(i),CD(i),Cm(i)] = aerodynamics(alpha(i));
end

figure(1)
plot(alpha,CL,'k',alpha,CD,'b')
axis([-180 180 -1 1.5])
grid on
grid minor
%hold on

figure(2)
plot(alpha,Cm,'k')
axis([-180 180 -0.6 0.6])
grid on
grid minor


function [CLi,CDi,Cmi] = aerodynamics(aaid)
    
    % WING PARAMETTERS
    delta = deg2rad(19.8); % sweep angle --> Tese 0 eq 2.16
    AR = 3.13; % wing aspect ratio --> Table H.1
    % ARt = 0.217; % winglet aspect ratio --> Table H.1
    
    % c = 0.2; % [m] cord?
    % cf = 0.063; % [m] flap cord
    
    CD0 = 0.02; % skin Friction Coefficient [-] --> Table H.1 [0.02 ; 0.04]
    k0=0.87; % oswald s efficiency factor [-] --> Table H.1 [0.85 ; 0.9]
    
    % Low Aspect Ratio (AR = 1.5-4) --> Tese PHD eq 3.11
    xp = 0.25;
    xe = 0.42;
    
    % Effect of Stall - Semi empirical coefficients for retangular flat plates (AR = 3)
    aLE = deg2rad(9.95);
    aTE = deg2rad(5.26);
    aaLE = deg2rad(14.28);
    aaTE = deg2rad(19.76);
    aahighS = 22; % High alpha regime starts at this angle

    % Linear Regime (low angle)
    CLalpha = 2*pi*cos(delta)/(2*cos(delta)/AR + sqrt(1+(2*cos(delta)/AR)^2));
    
    %CLalpha = 2*pi* (AR/(AR+2*(AR+4)/(AR+2))); % eq 3.9 Tese 21 PHD

    % eq 3.11 Tese 21 PHD
    Kp = CLalpha;
    Kv = pi;

    % High Alpha
    Cd90 = 1.98; % eq 3.17 Tese 21 PHD

    switch aaid
        %case -180:-171 % NEGATIVE LOW ANGLE
            
        %case -170:-91 % NEGATIVE HIGH ANGLE
        
        %case -90:-11 % NEGATIVE HIGH ANGLE
        
        %case -10:-1 % NEGATIVE LOW ANGLE
        
        case num2cell(0:10) % POSITIVE LOW ANGLE
            aai = deg2rad(aaid);

            % eq 3.8 Tese 21 PHD
            CLi = CLalpha*aai;
            CDi = CD0 + CLi^2/(pi*k0*AR);
            % Cmi = 0; % thin flat plate -> the a.c. and the center of pressure are coincident 
            
            % Low Aspect Ratio --> eq 3.10 Tese 21 PHD
            % CLi = Kp*sin(aai)*cos(aai)^2 + Kv*abs(sin(aai))*sin(aai)*cos(aai);
            % CDi = CD0 + CLi*abs(tan(aai));
            Cmi = -(xp-0.25)*Kp*sin(aai)*cos(aai) - (xe-0.25)*Kv*abs(sin(aai))*sin(aai);
         
        case num2cell(11:22) % POSITIVE STALL REGIME ANGLE
            aai = deg2rad(aaid);

            fLE = 0.5*(1-tanh(aLE*(aai-aaLE))); % ???
            fTE = 0.5*(1-tanh(aTE*(aai-aaTE))); % ???
    
            CLi = 0.25*(1+sqrt(fTE))^2*(Kp*sin(aai)*cos(aai)^2 + fLE^2*Kv*abs(sin(aai))*sin(aai)*cos(aai));
            CDi = CD0 + CLi*abs(tan(aai));
            Cmi = -0.25*(1+sqrt(fTE))^2*(0.0625*(-1+6*sqrt(fTE)-5*fTE)*Kp*sin(aai)*cos(aai) + 0.17*fLE^2*Kv*abs(sin(aai))*sin(aai));
            
        case num2cell(23:90) % POSITIVE HIGH ANGLE
            aai = deg2rad(aaid);

            Cn = Cd90 * sin(aai)/(0.56+0.44*sin(aai));
            CLi = Cn*cos(aai);  
            CDi = Cn*sin(aai);  
            
            % Low Aspect Ratio --> eq 3.18 Tese 21 PHD
            %Cn = Cd90 * sin(aai)/(0.56+0.44*sin(aai) - 0.41*(1-exp(-17/AR)));
            %Ca = 0.5*0.03*cos(aai);
    
            %CLi = Cn*cos(aai) - Ca*sin(aai); 
            %CDi = Cn*sin(aai) + Ca*cos(aai);
            Cmi = -Cn*(0.25-0.175*(1-2*aai/pi));

        %case 91:170 % POSITIVE HIGH ANGLE

        %case 171:180 % POSITIVE LOW ANGLE

        otherwise
            CLi = 0;
            CDi = 0;
            Cmi = 0;
    end
        
end

%         % EFFECT  OF CONTROL SURFACE DEFLECTION
%         aa = 0; % [deg] angle of attack
%         df = 39; % [deg] flap deflection
%         aa_0 = cf/c*df; % effective zero lift angle of attack due to flap deflection
%         aa_ = aa - aa_0; % effective angle of attack
%         ttf = acos(2*cf/c - 1);
%         tauf = 1 - (ttf - sin(ttf))/pi;
%         ettaf = correction_factor(df);
%         dCL = CLalpha * tauf * ettaf * deg2rad(df);
% 
%         fLE_ = 0.5*(1-tanh(aaLE*(aa*tLE*daa-aaLE)));
%         fTE_ = 0.5*(1-tanh(aaTE*(aa*tTE*daa-aaTE)));
%     
%         dCL2 = 0.25*(1+sqrt(fTE_)).^2 .* (Kp*sin(aa_).*cos(aa_).^2 + fLE_.^2*Kv.*abs(sin(aa_)).*sin(aa_).*cos(aa_))
%         
%         dCLmax = 0.47*dCL; % Mc Corninck pag 113 (depend on the flap cord and wing cord ratio)
%         [CLmax, index] = max(CL);
%         aamax = rad2deg(alpha(index));
%         CLmax_ = CLmax + dCLmax;
%         
% 
%         fun = @(a) 0.25*(1+sqrt(fTE)).^2 .* (Kp*sin(a).*cos(a).^2 + fLE.^2*Kv.*abs(sin(a)).*sin(a).*cos(a)) - CLmax_;
%         options = optimoptions('fsolve','Algorithm','levenberg-marquardt');
%         aa_CLmax_ = fsolve(fun, (aamax-0.1), options);
% end
%     
%     
%     Cd90 = 1.98;
%     
%     switch n
%         case 1 % HIGH ALPHA        
%     
%             Cn = Cd90 * sin(alphahigh)./(0.56+0.44*sin(alphahigh));
%             
%             CL = Cn.*cos(alphahigh);  CL1 = CL;
%             CD = Cn.*sin(alphahigh);  CD1 = CD;
%     
%             %figure(1)
%             plot(rad2deg(alphahigh),CL1,'k',rad2deg(alphahigh),CD1,'b') 
%         
%         case 2 % LOW ASPECT RATIO (AR = [1.5 ; 4])
%     
%             Cn = Cd90 * sin(alphahigh)./(0.56+0.44*sin(alphahigh) - 0.41*(1-exp(-17/AR)));
%             Ca = 0.5*0.03*cos(alphahigh);
%     
%             CL = Cn.*cos(alphahigh) - Ca.*sin(alphahigh);  CL2 = CL;
%             CD = Cn.*sin(alphahigh) + Ca.*cos(alphahigh);  CD2 = CD;
%             Cm = -Cn.*(0.25-0.175*(1-2*alphahigh/pi));
%     
%             %figure(2)
%             plot(rad2deg(alphahigh),CL2,'k',rad2deg(alphahigh),CD2,'b') 
%         
%         case 3 % EFFECT  OF CONTROL SURFACE DEFLECTION
%             
%             df = deg2rad(df);
%             Cd90_ = -4.26e-2* df^2 + 2.1e-1*df + 1.98 % eq 3.21 PHD Tese
%             c_ = sqrt((c-cf)^2 + cf^2 + 2*cf*(c-cf)*cos(df));
%             aaf = asin(cf/c_*sin(df));
%             aa_ = alphahigh + aaf;
%     
%             %figure(3)
%             plot(rad2deg(alphahigh),CL3,'k',rad2deg(alphahigh),CD3,'b') 
%         
%     end
% end