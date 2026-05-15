function [hittning_time_sim,hittning_time_th] = run_walk_cont(lambda,T_max,start_node,end_node,sleep_time,plot_) 
    close all
    if plot_
        figure
        set(gcf,'color','white')
    end
    
    node_count = length(lambda);
    
    x = zeros(node_count,T_max); % x is one for the current particle state
    x(start_node,1) = 1; % starting state is node 1
    fun = 0;
    if fun
        cords = round(rand(node_count,2)*2,1)*2;
    else
        cords = [0 0.5
                 1 1
                 1 0
                 2 0
                 2 1];
    end
    
    % Transition probability matrix for the directed ring
    omega = lambda*ones(node_count,1);
    omega_star = max(omega);
    P_bar = lambda/omega_star;
    for i = 1:node_count
        P_bar(i,i) = 1 - P_bar(i,:)*ones(node_count,1) + P_bar(i,i);
    end
    if start_node == end_node
        [V,D] = eig(P_bar');
        diag(D);
        pi_bar = V(:,1);
        norm_pi_bar = pi_bar/sum(pi_bar);
        hittning_time_th = 1/(omega(start_node)*norm_pi_bar(start_node));
    else
        omega_invers = 1./omega;
        P = diag(omega_invers)*lambda;
        P_hat = -eye(node_count) + P;
        P_hat(end_node,:) = [];
        P_hat(:,end_node) = [];
        omega_invers(end_node) = [];
        tau = P_hat\-omega_invers;
        if end_node < start_node
            hittning_time_th = tau(start_node-1);
        else
            hittning_time_th = tau(start_node);
        end
    end
    % Plot the graph and mark the node that the particle is in with red
    if plot_
        subplot(211)
        gplot(P_bar,cords,'-k');
        hold on
        for i = 1:5
            if x(i, 1) == 1
                scatter(cords(i,1),cords(i,2),200,'markeredgecolor','k','markerfacecolor', 'r');
                %quiver(cords(1,1),cords(1,2),cords(2,1),cords(2,2),"off")
                %pause(100)
            else
                scatter(cords(i,1),cords(i,2),200,'markeredgecolor','k','markerfacecolor', 'w');
            end
        
        end
        set(gca,'xtick',[],'ytick',[],'xcolor','w','ycolor','w')
    end
    %---- Simulate the particle moving around ----%
    r = omega_star;
    time_elapsed = 0;
    for i = 2:T_max
        t_next = -log(rand(1))/r;
        time_elapsed = time_elapsed + t_next;
        distribution = P_bar(x(:,i-1) == 1,:); %move the particle
        inter = cumsum(distribution);
        idx = find(rand(1) <= inter,1);
        x(idx, i) = 1; %update the state vector
        
        % plot the new location of the node
        if plot_
            pause(sleep_time) % sleep for some time
            subplot(211)
            for k = 1:5
                if x(k, i) == 1
                    scatter(cords(k,1),cords(k,2),200,'markeredgecolor','k','markerfacecolor', 'r');
                else
                    scatter(cords(k,1),cords(k,2),200,'markeredgecolor','k','markerfacecolor', 'w');
                end
            end
            subplot(212)
            tvec = [0 1:i];
            plot(tvec(1:end-1),(x(:, 1:i)'*[1 2 3 4 5]'), '-o')
        end
        if (idx == end_node) && (any(x(start_node,1:i) == 0))
            hittning_time_sim = time_elapsed;
            break
        end
        hittning_time_sim = -1;
    end
    if hittning_time_sim == -1
        disp("To long time to run")
    end 
end

