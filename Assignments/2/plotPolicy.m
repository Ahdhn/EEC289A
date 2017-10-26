function plotPolicy (policy, pi_num)
    figure;
    hold on;
    %pcolor(policy);
    policy_size = length(policy);
    [X,Y]=meshgrid(1:policy_size+1);
    plot(X,Y,'k');
    plot(Y,X,'k');
    axis off;
    if pi_num == 0
        title('\pi_0');
    elseif pi_num ==1
         title('\pi_1');
     elseif pi_num ==2
         title('\pi_2');
     elseif pi_num ==3
         title('\pi_3');
     elseif pi_num == 4
         title('\pi_4');
     elseif pi_num == 5
         title('\pi_5');           
     end
     
    text(6,0,'# Cars at second location');
    blah = text(-0.5,6,'# Cars at first location');  
    set(blah,'Rotation',90)
    
    
    [policy_size_x,policy_size_y] = size(policy);
    for p=1:policy_size_y
        for n=1:policy_size_x        
            text(n+.1,p+0.5,num2str(policy(p,n)),'FontWeight','bold');              
        end
    end
end