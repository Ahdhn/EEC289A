function plotData(data, mytitle)

    figure;
    axisX = 12:21;
    axisY = 1:10;
    surf(axisY,axisX,data, 'FaceColor','b');
    title(mytitle);
    ylabel('Player sum');
    xlabel('Dealer showing');

    
end