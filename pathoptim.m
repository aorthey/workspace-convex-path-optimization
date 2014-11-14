clf
clear all
theta = pi/2;
Rz90 = [cos(theta) -sin(theta) 0;
sin(theta) cos(theta) 0;
0 0 1;];

ay = [0,1,0]';
az = [0,0,1]';
b = 0.5;
A = eye(3);
A(3,3)=0;
tantheta = tan(pi/16);
c = [0;0;tantheta];
csmall = [0;0;tan(pi/200)];

Mh = 4;
t = linspace(0,1,Mh);
Fh= [t.^0;t.^1];
Kh= size(Fh,1);

Mk = 5;
t = linspace(0,1,Mk);
Fk= [t.^0;t.^1];
Kk= size(Fk,1);

M_ws = 30;
t = linspace(0,1,M_ws);
F_ws = [t.^0;t.^1;t.^2;t.^3;t.^4;t.^5;t.^6;t.^7;t.^8;t.^9;t.^10];
dF_ws = [0*t; 1*t.^0; 2*t.^1; 3*t.^2; 4*t.^3; 5*t.^4; 6*t.^5; 7*t.^6; 8*t.^7;9*t.^8;10*t.^9];
K_ws = size(F_ws,1);

F_wk = F_ws;
dF_wk = dF_ws;

cvx_begin 
        variable X(3,M_ws)
        variable dX(3,M_ws)

        variable W_ws(K_ws,3)

        variable W_wk(K_ws,3)

        variable Zh(3,Mh,M_ws)
        variable Wh(Kh,3,M_ws)

        variable Zk(3,Mk,M_ws)
        variable Wk(Kk,3,M_ws)

        variable zshift(3,1)

        variable lambda(1,M_ws)

        minimize norm(W_ws)+sum(sum(sum(Zh(:,:,:))))
        subject to
                %% specific Z constraints in each section
                for i=1:M_ws

                        Zk(:,:,i) == Wk(:,:,i)'*Fk
                        %% fix first pos of knee pos
                        Zk(:,1,i) == X(:,i)
                        %% last pos should be above certain height
                        az'*Zk(:,Mk,i) >= 3 %%TODO: how to stay on top of cone?

                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %% hip constraints
                        %% Rz90*lambda(:,i)*(Zk(:,Mk,i)-X(:,i)) == dX(:,i)
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                        Zh(:,1,i) == Zk(:,Mk,i)
                        Zh(:,:,i) == Wh(:,:,i)'*Fh
                        az'*Zh(:,Mh,i) >= 7 
                        quad_form(Zh(:,Mh,i)-X(:,i),A) <= csmall'*(Zh(:,Mh,i)-X(:,i))
                end
                %%% KNEE constraints
                for i=1:Mk
                        for j=1:M_ws
                                %% z stays in cone above X_j
                                quad_form(Zk(:,i,j)-X(:,j),A) <= c'*(Zk(:,i,j)-X(:,j))
                                %% z stays in ball around X_j
                                norm(Zk(:,i,j)-X(:,j)) <= 4
                        end
                end
                %% hip constraints
                for i=1:Mh
                        for j=1:M_ws
                                %% z stays in cone above X_j
                                quad_form(Zh(:,i,j)-Zk(:,Mk,j),A) <= c'*(Zh(:,i,j)-Zk(:,Mk,j))
                                %% z stays in ball around X_j
                                norm(Zh(:,i,j)-Zk(:,Mk,j)) <= 4
                        end
                end

                %% constraints on the position of the knee pos
                for i=1:M_ws-1
                        z = Zk(:,Mk,i);
                        x0 = X(:,i);
                        x1 = X(:,i+1);
                end

                %% additional floor constraints
                az'*X == b
                X == W_ws'*F_ws
                dX == W_ws'*dF_ws
                X(:,1) == [2;0;0.5]
                X(:,floor(M_ws/2)) == [2;3;0.5]
                X(:,M_ws) == [5;5;0.5]
                for i=1:M_ws-1
                        norm(X(:,i)-X(:,i+1)) <= 0.5
                end
cvx_end

dX = W_ws'*dF_ws

co =cvx_optval
if co < +inf
        figure(1)
        plot3(X(1,:),X(2,:),X(3,:),'*-r')
        hold on;
        for j=1:M_ws
                dX(:,j)=(dX(:,j)/norm(dX(:,j)))
                %P = [Zk(:,Mk,j)+dX(:,j) Zk(:,Mk,j)]
                P = [X(:,j)+Rz90*dX(:,j) X(:,j)]
                plot3(P(1,:),P(2,:),P(3,:),'*-m')
                hold on;
                P = [X(:,j)+dX(:,j) X(:,j)]
                plot3(P(1,:),P(2,:),P(3,:),'*-m')
                hold on;
                plot3(Zk(1,:,j),Zk(2,:,j),Zk(3,:,j),'*-b')
                hold on;
                plot3(Zh(1,:,j),Zh(2,:,j),Zh(3,:,j),'*-k')
                hold on;
        end
end
