clf
clear all

robotMaxHeight = 1.539;
robotMinKneeHeight = 0.2;
robotMaxKneeHeight = 0.3;
robotKneeAperture = pi/8;
robotHipAperture = pi/200;
robotNeckAperture = pi/16;
robotHeadAperture = pi/16;
robotKneeHeadAperture = pi/32;
distanceWaypoints = 0.5;
robotRadiusSweptVolume = 0.1;

az = [0,0,1]';
a_surface = az;
b_surface = 0;
A = eye(3);
A(3,3)=0;
cknee = [0;0;tan(robotKneeAperture)];
chip = [0;0;tan(robotHipAperture)];
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
cvx_begin 
        variable X(dimN,M_ws)
        variable W_ws(K_ws,dimN)

        variable Z(dimN,3,M_ws)
        variable W_wk(Kk,dimN,M_ws)
        variable W_wk90(Kk,dimN,M_ws)

        variable Zknee(dimN,M_ws)
        variable Zhip(dimN,M_ws)
        variable Zneck(dimN,M_ws)
        variable Zhead(dimN,M_ws)

        variable lambda(M_ws)

        minimize sum(lambda)
        
        subject to
                lambda >= 0
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%% objective value epigraph minimizations
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                for i=1:M_ws
                        norm(Z(:,1,i)-Z(:,2,i))+norm(Z(:,2,i)-Z(:,3,i)) <= lambda(i)
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
                        quad_form(Z(:,2,i)-X(:,i),A) <= cknee'*(Z(:,2,i)-X(:,i))

                        robotMinKneeHeight <= a_surface'*Z(:,2,i) <= robotMaxKneeHeight
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%% Constraints on hip
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                for i=1:M_ws
                        Zknee(:,i) == Z(:,2,i)
                        quad_form(Zhead(:,i)-Zknee(:,i),A) <= ckneehead'*(Zhead(:,i)-Zknee(:,i))

                        quad_form(Zhip(:,i)-X(:,i),A) <= chip'*(Zhip(:,i)-X(:,i))
                        a_surface'*Z(:,2,i) <= a_surface'*Zhip(:,i) <= a_surface'*Zneck(:,i)
                        quad_form(Zneck(:,i)-Zhip(:,i),A) <= cneck'*(Zneck(:,i)-Zhip(:,i))
                        a_surface'*Zneck(:,i) <= a_surface'*Zhead(:,i)
                        quad_form(Zhead(:,i)-Zneck(:,i),A) <= chead'*(Zhead(:,i)-Zneck(:,i))
                        a_surface'*Zhead(:,i) <= robotMaxHeight
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
                        dsum(i) = d;
                        if d > 0.001
                                disp(sprintf('Warning: non-orthogonality for i=%d', i))
                        end
                else
                        disp(sprintf('Warning: Z and X values not aligned for i=%d', i))
                end
        end
        disp(sprintf('Orthogonality SCORE=%f (should be 0)', sum(dsum)));

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
                P = [Zhip(:,i) Zneck(:,i)];
                plot3(P(1,:),P(2,:),P(3,:),colorSkeleton,'LineWidth',LW)

                Cylinder(P(:,1),P(:,2),robotRadiusSweptVolume,Ncyl, colorSV, alphaSV, 0,0);
                SpherePlot = surf(xS+P(1,1),yS+P(2,1),zS+P(3,1));
                set(SpherePlot,'FaceColor',colorSV)
                set(SpherePlot,'FaceAlpha',alphaSV)
                set(SpherePlot,'EdgeAlpha',0)
                hold on;
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
