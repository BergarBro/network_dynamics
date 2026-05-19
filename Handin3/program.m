%% Hand In 3
%% Part 1
clc
clear

lambda = [ 0  2/5 1/5  0   0 ;
           0   0  3/4 1/4  0 ;
          1/2  0   0  1/2  0 ;
           0   0  1/3  0  2/3;
           0  1/3  0  1/3  0 ];

start_node = 1;
end_node = 5;
sleep_time = 0;
plot_visualizaer = 1;
times = 1;


time_acc = 0;
max_time = 0;
for i = 1:times
    [hitting_time_sim,hitting_time_th] = run_walk_cont(lambda,1000,start_node,end_node,sleep_time,plot_visualizaer);
    time_acc = time_acc + hitting_time_sim;
    max_time = max(max_time,hitting_time_sim);
end
average_hittning_time_sim = time_acc/times;
fprintf("Simulated Return/Hitting Time: %f \n", round(average_hittning_time_sim,4))
fprintf("Theretical Return/Hitting Time: %f \n", round(hitting_time_th,4))
radio = average_hittning_time_sim/hitting_time_th;
fprintf("Ratio: %f \n", round(radio,4))
fprintf("Maximum Return/Hitting Time: %f \n",max_time)


%% Part 2.a

clc
clear
close all

nbr_of_colors = 2;
nbr_of_nodes = 10;
all_colors = [1 0 0; 0 1 0; 0 0 1; 1 1 0; 1 0 1; 0 1 1; 1 1 1; 0 0 0];
colors_rgb = all_colors(1:nbr_of_colors,:);
colors = 1:nbr_of_colors;

show_plot = 1;
min_potential = 0;

nodes = ones(nbr_of_nodes,1);
eta = @(x) x/100;

W = zeros(nbr_of_nodes) + [zeros(nbr_of_nodes,1) [eye(nbr_of_nodes-1);zeros(1,nbr_of_nodes-1)]] + [zeros(nbr_of_nodes,1) [eye(nbr_of_nodes-1);zeros(1,nbr_of_nodes-1)]]';
x = 1:nbr_of_nodes;
y = ones(nbr_of_nodes,1);

plot_nodes(colors_rgb, nodes, eta, @isequal, W, x, y, nbr_of_nodes, nbr_of_colors, colors, show_plot, min_potential);

%% Part 2.b

clc
clear
close all

load('data/wifi.mat','-ascii')
load('data/coord.mat','-ascii')

x = coord(:,1);
y = coord(:,2);
W = wifi;

eta_b = @(x) x/100;
eta_1 = @(x) x/1000;
eta_2 = @(x) 50;
eta_3 = @(x) 50/(1+exp(0.5*(100-x)));
eta = {eta_b,eta_1,eta_2,eta_3};
all_colors = [1 0 0; 0 1 0; 0 0 1; 1 1 0; 1 0 1; 0 1 1; 1 1 1; 0 0 0];

nbr_of_colors = 8;
nbr_of_nodes = length(wifi);
colors_rgb = all_colors(1:nbr_of_colors,:);
colors = 1:nbr_of_colors;
show_plot = 0;
min_potential = 4;
times = 1;                    % Times to run the simulation
nodes = ones(nbr_of_nodes,1); % Initial state

nbr_of_eta = length(eta);
all_u = zeros(times,nbr_of_eta);
all_time_steps = zeros(times,nbr_of_eta);
for j = 1:nbr_of_eta
    e = eta{j};
    for i = 1:times
        [all_u(i,j), all_time_steps(i,j)] = plot_nodes(colors_rgb, nodes, e, @custom_c, W, x, y, nbr_of_nodes, nbr_of_colors, colors, show_plot, min_potential);
    end
end
res_u = mean(all_u,1);
res_time_steps = mean(all_time_steps,1);

fprintf("\nTimes Run: %d \n" ,times)
for j = 1:nbr_of_eta
    fprintf("Using %s(x) = %s \n",char(951),func2str(eta{j}))
    fprintf("Avrages Min Potential: %f \n",round(res_u(j),2))
    fprintf("Avrages Timesteps: %f \n \n",round(res_time_steps(j),2))
end

function [u_t, time_steps] = plot_nodes(colors_rgb, nodes, eta, c, W, x, y, nbr_of_nodes, nbr_of_colors, colors, show_plot,best_value)
    close all
    G = graph(W);
    scale = 0.7;
    if show_plot
        figure('Position', [100, 100, 1000, 500]);
        subplot(1,2,1)
        p = plot(G, 'XData', x, 'YData', y);
        p.MarkerSize = 15*scale;
        p.NodeColor = colors_rgb(nodes,:);
        
        p.NodeLabel = {};
    
        text_obj = gobjects(numnodes(G),1);
    
        for i = 1:numnodes(G)
        text_obj(i) = text(x(i), y(i), sprintf('%d', nodes(i)), ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', ...
            'FontSize', 12*scale, ...
            'FontWeight', 'bold');
        end
    end
    t = 1;
    u = zeros(1,1);
    while true
        the_one_node = randperm(nbr_of_nodes,1);
        prob_dis = zeros(nbr_of_colors,1);
        for i = 1:nbr_of_colors
            prob_dis(i) = prob_of_new_color(W,nodes,the_one_node,colors(i),c,eta(t),colors);
        end
        inter = cumsum(prob_dis);
        new_color = colors(find(rand(1) <= inter,1));
        nodes(the_one_node) = new_color;
        if show_plot
            text_obj(the_one_node).String = sprintf('%d', nodes(the_one_node));
            if (new_color == 8)||(new_color == 3)
                text_obj(the_one_node).Color = "w";
            else
                text_obj(the_one_node).Color = "black";
            end
            p.NodeColor = colors_rgb(nodes,:);
            if mod(t,10) == 0
                subplot(2,2,2)
                [counts,edges] = histcounts(nodes);
                centers = (edges(1:end-1)+edges(2:end))/2;
                b = bar(centers,counts);
                b.FaceColor = 'flat';
                b.CData = colors_rgb;
                xlabel("Colors")
                ylabel("Amounts")
                ylim([0 30])
                drawnow;
            end
        end
        u(t) = utility(W,nodes,c);
        span = 100;
        if (u(t) == 0)||(t>span)&&(sum(abs(diff(u(t-span:t))))==0)||(t>20000)||(u(t)==best_value)
            if show_plot
                subplot(2,2,2)
                [counts,edges] = histcounts(nodes);
                centers = (edges(1:end-1)+edges(2:end))/2;
                b = bar(centers,counts);
                b.FaceColor = 'flat';
                b.CData = colors_rgb;
                xlabel("Colors")
                ylabel("Amounts")
                ylim([0 30])
                subplot(2,2,4)
                plot(u)
                xlabel("Timesteps")
                ylabel("Potential")
                title("Potential U(t) = " + int2str(u(t)))
            end
            break
        end
        t = t + 1;
    end
    time_steps = t;
    u_t = u(t);
    fprintf("Timestep: %i  U_t: %d \n", t, u_t);
end

function c = custom_c(x,y)
    c = max(0,2-abs(x-y));
end

function u = utility(W,nodes,c)
    acc = 0;
    for i = 1:length(nodes)
        acc = acc + sum_neighbors(W,nodes,i,nodes(i),c);
    end
    u = 0.5*acc;
end

function prob = prob_of_new_color(W,nodes,chosen_node,new_color,c,eta_t,colors)
    numerator = exp(-eta_t*sum_neighbors(W,nodes,chosen_node,new_color,c));
    denomenator = 0;
    for col = colors
        denomenator = denomenator + exp(-eta_t*sum_neighbors(W,nodes,chosen_node,col,c));
    end
    prob = numerator / denomenator;
end

function value = sum_neighbors(W,nodes,chosen_node,new_color,c)
    acc = 0;
    for j = 1:length(nodes)
        acc = acc + W(chosen_node,j)*c(nodes(j),new_color);
    end
    value = acc;
end
