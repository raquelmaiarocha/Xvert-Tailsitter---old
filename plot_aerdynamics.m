function plot_aerdynamics()

    alpha = -180:1:180; % length(alpha) = 361
    df = [-39,0,39];
    
    CL = zeros(3,length(alpha)); % dfmin(1) df0(2) dfmax(3)
    CD = zeros(3,length(alpha)); % dfmin(1) df0(2) dfmax(3)
    Cm = zeros(3,length(alpha)); % dfmin(1) df0(2) dfmax(3)
    
    for j=1:3
        for i=1:361
            [CL(j,i),CD(j,i),Cm(j,i)] = aerodynamics(alpha(i),df(j));
        end
    end
    
    figure(1)
    subplot(2,1,1)
    plot(alpha,CL(2,:),'k.-',alpha,CL(1,:),'k-',alpha,CL(3,:),'k--')
    hold on
    plot(alpha,CD(2,:),'b.-',alpha,CD(1,:),'b-',alpha,CD(3,:),'b--')
    axis([-180 180 -1 1.5])
    grid on
    grid minor
    xlabel('$\alpha$ [ deg ]','Interpreter','latex','FontSize', 12)
    ylabel('$C_L, C_D$','Interpreter','latex','FontSize', 12)
    legend({'\delta = 0','\delta_{min}','\delta_{max}'},'FontSize', 12,'NumColumns',3)
    
    subplot(2,1,2)
    plot(alpha,Cm(2,:),'k.-',alpha,Cm(1,:),'k-',alpha,Cm(3,:),'k--')
    axis([-180 180 -0.6 0.6])
    grid on
    grid minor
    xlabel('$\alpha$ [ deg ]','Interpreter','latex','FontSize', 12)
    ylabel('$C_m, C_D$','Interpreter','latex','FontSize', 12)
    legend({'\delta = 0','\delta_{min}','\delta_{max}'},'FontSize', 12,'NumColumns',3)

end