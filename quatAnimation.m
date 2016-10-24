function quatAnimation(quaternion)
    % INPUT:
    % quaternion - quaternions array [length x 4], in [w i j k] order
    baseVector = [0.5774 0.5774 0.5774];
    basePosition = [0 0 0];
    pointerVector = [0 0.2 0];
    len = length(quaternion);
    skipper = round(len/1000);
%     skipper = 1;
    figure
    grid on
    daspect([1 1 1]);
    rotate3d on
    view(105, 40);
    xlim([-1 1]); xlabel('X lived here');
    ylim([-1 1]); ylabel('Y goes there');
    zlim([-1 1]); zlabel('Z is homeless');
    line([-1 1], [0 0], [0 0], 'LineWidth', 2, 'Color', 'b');
    line([0 0], [-1 1], [0 0], 'LineWidth', 2, 'Color', 'g');
    line([0 0], [0 0], [-1 1], 'LineWidth', 2, 'Color', 'r');
    legend X Y Z
    line([1 0.9], [0 0], [0 0.1], 'LineWidth', 2, 'Color', 'b');
    line([1 0.9], [0 0.1], [0 0], 'LineWidth', 2, 'Color', 'b');

    line([0 0], [1 0.9], [0 0.1], 'LineWidth', 2, 'Color', 'g');
    line([0 0.1], [1 0.9], [0 0], 'LineWidth', 2, 'Color', 'g');

    line([0 0.1], [0 0], [1 0.9], 'LineWidth', 2, 'Color', 'r');
    line([0 0], [0 0.1], [1 0.9], 'LineWidth', 2, 'Color', 'r');
    
    line([1 -1], [-1 1], [-1 1], 'LineWidth', 1, 'LineStyle','--', 'Color', [0.5, 0.5, 0.5]);
    line([1 -1], [1 -1], [-1 1], 'LineWidth', 1, 'LineStyle','--', 'Color', [0.5, 0.5, 0.5]);
    line([-1 1], [-1 1], [-1 1], 'LineWidth', 1, 'LineStyle','--', 'Color', [0.5, 0.5, 0.5]);
    line([-1 1], [1 -1], [-1 1], 'LineWidth', 1, 'LineStyle','--', 'Color', [0.5, 0.5, 0.5]);
    
    pic = line([basePosition(1) baseVector(1)],[basePosition(2) baseVector(2)],[basePosition(3) baseVector(3)], 'LineWidth', 2, 'Color', 'k');
    pVec = baseVector + pointerVector;
%     ptr = line([baseVector(1) pVec(1)],[baseVector(2) pVec(2)],[baseVector(3) pVec(3)], 'LineWidth', 2, 'Color', 'k');
%     pause
    prevPos = pVec;
    for k=1:skipper:len
        pause(0.01);
        delete(pic);
%         delete(ptr);
        vector = rotateVQ(baseVector, quaternion(k,:));
        pVec = rotateVQ(pointerVector, quaternion(k,:));
        pVec = vector + pVec;
        pic = line([basePosition(1) vector(1)],[basePosition(2) vector(2)],[basePosition(3) vector(3)], 'LineWidth', 2, 'Color', 'k');
%         ptr = line([vector(1) pVec(1)],[vector(2) pVec(2)],[vector(3) pVec(3)], 'LineWidth', 2, 'Color', 'k');
        pVec = vector;
        dotVec = pVec + 0.02;
%         line([pVec(1) dotVec(1)],[pVec(2) dotVec(2)],[pVec(3) dotVec(3)], 'LineWidth', 2, 'Color', 'k');
        line([prevPos(1) pVec(1)],[prevPos(2) pVec(2)],[prevPos(3) pVec(3)], 'LineWidth', 2, 'Color', [0.8,0.4,0]);
        prevPos = pVec;
    end;
    
   
end