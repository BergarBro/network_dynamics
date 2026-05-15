% Simple code-skeleton to help you visualize a single particle moving
% around in a directed ring-graph.
%
% Again, this is in no way compulsory or needed in order to solve the
% handin, but might prove useful for you. We will simulate a particle
% moving through a directed ring in a deterministic manner. The
% particle stays in each node exactly 1 time-unit, and then moves on
% to the next node.
%
% What you should (if you wish) take away from this code-skeleton is
% how you could visualize a particle moving around in a network.
clear all
close all
clc

% number of iterations or expected running time
Tmax = 20; 

figure
set(gcf,'color','white')

node_count = 5;

x = zeros(node_count,Tmax); % x is one for the current particle state
x(1,1) = 1; % starting state is node 1
if fun
    cords = round(rand(node_count,2)*2,1)*2;
else
    cords = [1 0
             2 2
             0 2
             4 0
             4 2];
end

% Transition probability matrix for the directed ring
W = [0 2 1 0 0;
     1 0 5 0 1;
     0 1 0 4 0;
     0 0 6 0 1;
     0 0 0 0 1];
W = rand(node_count,node_count);
omega = W*ones(node_count,1);
P = diag(1./omega)*W;
% Plot the graph and mark the node that the particle is in with red
subplot(211)
gplot(P,cords,'-k');
hold on
for i = 1:5
    if x(i, 1) == 1
        scatter(cords(i,1),cords(i,2),200,'markeredgecolor','k','markerfacecolor', 'r');
    else
        scatter(cords(i,1),cords(i,2),200,'markeredgecolor','k','markerfacecolor', 'w');
    end

end
set(gca,'xtick',[],'ytick',[],'xcolor','w','ycolor','w')

%---- Simulate the particle moving around ----%

zzz = 0; % time that the particle waits in the node before moving
         % on to the next one. Here we have just unit-time and no
         % randomness...   
n = 1; %node that the particle is currently in
%
for i = 2:Tmax
    delta = P(x(:,i-1) == 1,:); %move the particle
    inter = cumsum(delta);
    r = rand(1);
    idx = find(r <= inter,1);
    x(idx, i) = 1; %update the state vector
    pause(zzz) % sleep for some time
    
    % plot the new location of the node
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

