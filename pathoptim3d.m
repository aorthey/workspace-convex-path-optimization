clf
clear all

robotMaxHeight = 1.539;
robotMinKneeHeight = 0.2;
robotMaxKneeHeight = 0.3;

robotKneeAperture = pi/8;
robotHipAperture = pi/100;
robotWaistAperture = pi/16;
robotNeckAperture = pi/16;
robotHeadAperture = pi/16;
robotKneeHipAperture = pi/4;

robotKneeHeadAperture = pi/32; %%TODO:debugging only

distanceWaypoints = 0.5;
robotRadiusSweptVolume = 0.1;

distanceKneeFoot = 0.4;
distanceHipKnee = 0.45;
distanceWaistHip = 0.2;
distanceNeckWaist = 0.3;
distanceHeadNeck = 0.20;

distanceConstants = [distanceKneeFoot distanceHipKnee distanceWaistHip distanceNeckWaist distanceHeadNeck];
slackRelaxation=0.02;

minDistanceKneeFoot = distanceKneeFoot*(cos(robotKneeAperture));
minDistanceHipFoot = minDistanceKneeFoot+distanceHipKnee*(cos(robotHipAperture));
minDistanceWaistFoot = minDistanceHipFoot+distanceWaistHip*(cos(robotWaistAperture));
minDistanceNeckFoot = minDistanceWaistFoot+distanceNeckWaist*(cos(robotNeckAperture));
minDistanceHeadFoot = minDistanceNeckFoot+distanceHeadNeck*(cos(robotHeadAperture));

maxDistanceKneeFoot = distanceKneeFoot;
maxDistanceHipFoot = maxDistanceKneeFoot+distanceHipKnee;
maxDistanceWaistFoot = maxDistanceHipFoot+distanceWaistHip;
maxDistanceNeckFoot = maxDistanceWaistFoot+distanceNeckWaist;
maxDistanceHeadFoot = maxDistanceNeckFoot+distanceHeadNeck;

az = [0,0,1]';
a_surface = az;
b_surface = 0;
A = eye(3);
A(3,3)=0;

cknee = [0;0;tan(robotKneeAperture)];
ckneehip = [0;0;tan(robotKneeHipAperture)];
chip = [0;0;tan(robotHipAperture)];
cwaist = [0;0;tan(robotWaistAperture)];
cneck = [0;0;tan(robotNeckAperture)];
chead = [0;0;tan(robotHeadAperture)];
ckneehead = [0;0;tan(robotKneeHeadAperture)];

M_ws = 23;
t = linspace(0,1,M_ws);
F_ws = [t.^0;t.^1;t.^2;t.^3;t.^4;t.^5;t.^6;t.^7;t.^8;t.^9;t.^10];
dF_ws = [0*t; 1*t.^0; 2*t.^1; 3*t.^2; 4*t.^3; 5*t.^4; 6*t.^5; 7*t.^6; 8*t.^7;9*t.^8;10*t.^9];
K_ws = size(F_ws,1);

Mk = 2;
t = linspace(0,1,Mk);
Fk= [t.^0;t.^1;];
Kk= size(Fk,1);

Loffset=200;
dimN = 3;
dualN = M_ws*5;
%cvx_precision(1e-3)
cvx_precision low
cvx_begin 
        variable X(dimN,M_ws)
        variable W_ws(K_ws,dimN)

        variable Z(dimN,3,M_ws)
        variable W_wk(Kk,dimN,M_ws)
        variable W_wk90(Kk,dimN,M_ws)

        variable Zknee(dimN,M_ws)
        variable Zhip(dimN,M_ws)
        variable Zwaist(dimN,M_ws)
        variable Zneck(dimN,M_ws)
        variable Zhead(dimN,M_ws)

        variable lambda(M_ws)
        variable lambda2(M_ws)
        variable lambda3(M_ws)
        variable lambda4(M_ws)
        variable lambda5(M_ws)
        variable slack_var(M_ws,5)
        %dual variable y{dualN}

        minimize sum(lambda)+norm(slack_var)+norm(lambda2)+sum(lambda3)+sum(lambda4)+sum(lambda5)
        %minimize norm(slack_var)
        
        subject to
                lambda >= 0
                lambda2 >= 0
                lambda3 >= 0
                lambda4 >= 0
                lambda5 >= 0
                slack_var >= 0
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%% objective value epigraph minimizations
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                for i=1:M_ws
                        norm(Z(:,1,i)-Z(:,2,i))+norm(Z(:,2,i)-Z(:,3,i)) <= lambda(i)
                        norm(Zhip(:,i)-[0;0;20000]) <= lambda2(i)
                        norm(Zwaist(:,i)-[0;0;20000]) <= lambda3(i)
                        norm(Zneck(:,i)-[0;0;20000]) <= lambda4(i)
                        norm(Zhead(:,i)-[0;0;20000]) <= lambda5(i)
                end
                %norm(W_ws) <= lambda(M_ws+1)

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%% orthogonal distance constraint
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                for i=1:M_ws
                        j=i+1;
                        distLambda=-Loffset;
                        if i==M_ws
                                j=i-1;
                                distLambda=Loffset;
                        end

                        X(:,i) == W_wk(:,:,i)'*[1;0]
                        X(:,j) == W_wk(:,:,i)'*[1;1]
                        W_wk90(:,:,i) == [W_wk(1,1,i) W_wk(1,2,i) W_wk(1,3,i); W_wk(2,2,i) -W_wk(2,1,i) 0]

                        Z(:,3,i) == W_wk90(:,:,i)'*[1;distLambda] + a_surface
                        Z(:,1,i) == W_wk90(:,:,i)'*[1;0] + a_surface

                        %%% constraints of the knee wrt the foot
                        %%% (cone+sphere+orthogonality)

                        %robotMinKneeHeight <= a_surface'*Z(:,2,i) <= robotMaxKneeHeight
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%% Constraints on hip
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                for i=1:M_ws
                        Zknee(:,i) == Z(:,2,i)

                        %% TODO: for debug:
                        quad_form(Zhead(:,i)-Zknee(:,i),A) <= ckneehead'*(Zhead(:,i)-Zknee(:,i))
                        %%%%%%%%%%%%%%%%%%%%%%%%

                        quad_form(Zknee(:,i)-X(:,i),A) <= cknee'*(Zknee(:,i)-X(:,i))
                        quad_form(Zhip(:,i)-X(:,i),A) <= chip'*(Zhip(:,i)-X(:,i))
                        quad_form(Zhip(:,i)-Zknee(:,i),A) <= ckneehip'*(Zhip(:,i)-Zknee(:,i))
                        quad_form(Zwaist(:,i)-Zhip(:,i),A) <= cwaist'*(Zwaist(:,i)-Zhip(:,i))
                        quad_form(Zneck(:,i)-Zwaist(:,i),A) <= cneck'*(Zneck(:,i)-Zwaist(:,i))
                        quad_form(Zhead(:,i)-Zneck(:,i),A) <= chead'*(Zhead(:,i)-Zneck(:,i))

                        norm(Zknee(:,i) - X(:,i)) <= slack_var(i,1)
                        norm(Zhip(:,i) - Zknee(:,i)) <= slack_var(i,2)
                        norm(Zwaist(:,i) - Zhip(:,i)) <= slack_var(i,3)
                        norm(Zneck(:,i) - Zwaist(:,i)) <= slack_var(i,4)
                        norm(Zhead(:,i) - Zneck(:,i)) <= slack_var(i,5)

                        %y{(i-1)*5+1}: norm(Zknee(:,i) - X(:,i)) <= slack_var(i,1)
                        %y{(i-1)*5+2}: norm(Zhip(:,i) - Zknee(:,i)) <= slack_var(i,2)
                        %y{(i-1)*5+3}: norm(Zwaist(:,i) - Zhip(:,i)) <= slack_var(i,3)
                        %y{(i-1)*5+4}: norm(Zneck(:,i) - Zwaist(:,i)) <= slack_var(i,4)
                        %y{(i-1)*5+5}: norm(Zhead(:,i) - Zneck(:,i)) <= slack_var(i,5)

                        maxDistanceKneeFoot >= az'*Zknee(:,i) >= minDistanceKneeFoot
                        maxDistanceHipFoot >= az'*Zhip(:,i) >= minDistanceHipFoot
                        maxDistanceWaistFoot >= az'*Zwaist(:,i) >= minDistanceWaistFoot
                        maxDistanceNeckFoot >= az'*Zneck(:,i) >= minDistanceNeckFoot
                        maxDistanceHeadFoot >= az'*Zhead(:,i) >= minDistanceHeadFoot

                        distanceKneeFoot-slackRelaxation <= slack_var(i,1) <= distanceKneeFoot
                        distanceHipKnee-slackRelaxation <= slack_var(i,2) <= distanceHipKnee
                        distanceWaistHip-slackRelaxation <= slack_var(i,3) <= distanceWaistHip
                        distanceNeckWaist-slackRelaxation <= slack_var(i,4) <= distanceNeckWaist
                        distanceHeadNeck-slackRelaxation <= slack_var(i,5) <= distanceHeadNeck

                        %%% hard constraints for safety (could be removed after
                        %%% serious testing)
                        %norm(X(:,i)-Zhead(:,i)) <= robotMaxHeight

                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%% Constraints on neck
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%% Constraints on ground points X
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                for i=1:M_ws-1
                        norm(X(:,i)-X(:,i+1)) <= distanceWaypoints
                end
                %% on surface
                for i=1:M_ws-1
                        az'*X(:,i) == b_surface
                end
                %% in functional space
                X(:,:) == W_ws'*F_ws

                %% through specific points
                X(:,1) == [0.9;1.1;b_surface]
                X(:,floor(M_ws/3)) == [3;1;b_surface]
                X(:,floor(2*M_ws/3)) == [3;3;b_surface]
                X(:,M_ws) == [5;4;b_surface]
cvx_end

co=cvx_optval;
if co < +inf
        %% evaluate the orthogonality of the knee vectors

        for i=1:M_ws-1
                if norm(X(:,i)-Z(:,1,i)) <= norm(0.01+az)
                        z = Zknee(:,i)-X(:,i);
                        x = X(:,i+1)-X(:,i);
                        d = z'*x;
                        d_orthogonality_sum(i) = d;
                        if d > 0.001
                                disp(sprintf('Warning: non-orthogonality for i=%d', i))
                        end

                        zz = Zknee(:,i) - X(:,i);
                        zh = Zhip(:,i) - Zknee(:,i);
                        zw = Zwaist(:,i) - Zhip(:,i);
                        zn = Zneck(:,i) - Zwaist(:,i);
                        zhead = Zhead(:,i) - Zneck(:,i);

                        d_kneedistance(i) = abs(norm(zz)-distanceKneeFoot);
                        d_hipdistance(i) = abs(norm(zh)-distanceHipKnee);
                        d_waistdistance(i) = abs(norm(zw)-distanceWaistHip);
                        d_neckdistance(i) = abs(norm(zn)-distanceNeckWaist);
                        d_headdistance(i) = abs(norm(zhead)-distanceHeadNeck);
                        
                else
                        disp(sprintf('Warning: Z and X values not aligned for i=%d', i))
                end
        end
        disp(sprintf('Orthogonality SCORE=%f (should be 0)', sum(d_orthogonality_sum)));
        disp(sprintf('Knee  distance SCORE=%f (should be 0)', sum(d_kneedistance)));
        disp(sprintf('Hip   distance SCORE=%f (should be 0)', sum(d_hipdistance)));
        disp(sprintf('Waist distance SCORE=%f (should be 0)', sum(d_waistdistance)));
        disp(sprintf('Neck  distance SCORE=%f (should be 0)', sum(d_neckdistance)));
        disp(sprintf('Head  distance SCORE=%f (should be 0)', sum(d_headdistance)));

        fig=figure(1);
        set(fig, 'Color', [1 1 1]);

        
        plot3(X(1,:),X(2,:),X(3,:),'*-k','LineWidth',3);
        hold on;
        [xS,yS,zS] = sphere;
        xS = robotRadiusSweptVolume*xS;
        yS = robotRadiusSweptVolume*yS;
        zS = robotRadiusSweptVolume*zS;

        Ncyl = 5;
        LW = 2;
        colorSkeleton = '*-k';
        colorSV = 'r';
        alphaSV = 0.1;
        for i=1:size(Z,3)
                P = [X(:,i) Zknee(:,i)];
                plot3(P(1,:),P(2,:),P(3,:),colorSkeleton,'LineWidth',LW)
                hold on;

                Cylinder(P(:,1),P(:,2),robotRadiusSweptVolume,Ncyl, colorSV, alphaSV, 0,0);
                SpherePlot = surf(xS+P(1,1),yS+P(2,1),zS+P(3,1));
                set(SpherePlot,'FaceColor',colorSV)
                set(SpherePlot,'FaceAlpha',alphaSV)
                set(SpherePlot,'EdgeAlpha',0)

                hold on;
                P = [Zknee(:,i) Zhip(:,i)];
                plot3(P(1,:),P(2,:),P(3,:),colorSkeleton,'LineWidth',LW)

                Cylinder(P(:,1),P(:,2),robotRadiusSweptVolume,Ncyl, colorSV, alphaSV, 0,0);
                SpherePlot = surf(xS+P(1,1),yS+P(2,1),zS+P(3,1));
                set(SpherePlot,'FaceColor',colorSV)
                set(SpherePlot,'FaceAlpha',alphaSV)
                set(SpherePlot,'EdgeAlpha',0)

                hold on;
                P = [Zhip(:,i) Zwaist(:,i)];
                plot3(P(1,:),P(2,:),P(3,:),colorSkeleton,'LineWidth',LW)

                Cylinder(P(:,1),P(:,2),robotRadiusSweptVolume,Ncyl, colorSV, alphaSV, 0,0);
                SpherePlot = surf(xS+P(1,1),yS+P(2,1),zS+P(3,1));
                set(SpherePlot,'FaceColor',colorSV)
                set(SpherePlot,'FaceAlpha',alphaSV)
                set(SpherePlot,'EdgeAlpha',0)
                hold on;
                %%%%%%%%%%%%%%%%%%%%%%%
                hold on;
                P = [Zwaist(:,i) Zneck(:,i)];
                plot3(P(1,:),P(2,:),P(3,:),colorSkeleton,'LineWidth',LW)

                Cylinder(P(:,1),P(:,2),robotRadiusSweptVolume,Ncyl, colorSV, alphaSV, 0,0);
                SpherePlot = surf(xS+P(1,1),yS+P(2,1),zS+P(3,1));
                set(SpherePlot,'FaceColor',colorSV)
                set(SpherePlot,'FaceAlpha',alphaSV)
                set(SpherePlot,'EdgeAlpha',0)
                hold on;
                %%%%%%%%%%%%%%%%%%%%%%%
                P = [Zneck(:,i) Zhead(:,i)];
                plot3(P(1,:),P(2,:),P(3,:),colorSkeleton,'LineWidth',LW)

                Cylinder(P(:,1),P(:,2),robotRadiusSweptVolume,Ncyl, colorSV, alphaSV, 0,0);
                SpherePlot = surf(xS+P(1,1),yS+P(2,1),zS+P(3,1));
                set(SpherePlot,'FaceColor',colorSV)
                set(SpherePlot,'FaceAlpha',alphaSV)
                set(SpherePlot,'EdgeAlpha',0)

                SpherePlot = surf(xS+P(1,2),yS+P(2,2),zS+P(3,2));
                set(SpherePlot,'FaceColor',colorSV)
                set(SpherePlot,'FaceAlpha',alphaSV)
                set(SpherePlot,'EdgeAlpha',0)
        end
        axis equal
end
