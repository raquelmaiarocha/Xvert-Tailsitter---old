function plot_propeller(xvert)

    J = 0:0.01:0.8;

    CT = (xvert.CT2*(J.^2)) + (xvert.CT1*(J)) + (xvert.CT0);
    CP = (xvert.CP2*(J.^2)) + (xvert.CP1*(J)) + (xvert.CP0);
    
    figure(2)
    subplot(2,1,1)
    plot(J,CT,'k',J,CP,'b')
    axis([0 0.8 -0.05 0.15])
    grid on
    grid minor
    xlabel('advanced ratio J [ 1/rev ]','Interpreter','latex','FontSize', 12)
    legend({'$C_T$','$C_P$'},'Interpreter','latex','FontSize', 12)
    hold on

    dth = 0:0.01:1;
    Vbatt = [7.0 7.5 8.0 8.5]';
    omega = Vbatt.^(0.8)*(xvert.wp2*dth.^2 + xvert.wp1*dth + xvert.wp0);
    
    subplot(2,1,2)
    plot(dth,omega,'k')
    axis([0 1 0 1500])
    grid on
    grid minor
    xlabel('throttle $\tau$','Interpreter','latex','FontSize', 12)
    ylabel('angular speed $\omega_p$ [rad/s]','Interpreter','latex','FontSize', 12)
    legend({'$V_{batt}=7.0$ [V]','$V_{batt}=7.5$ [V]','$V_{batt}=8.0$ [V]','$V_{batt}=8.5$ [V]',},'Interpreter','latex','FontSize', 12,'Location','southeast')
end