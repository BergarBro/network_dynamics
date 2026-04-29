% Handin2/program.m

clc
clear

load('traffic.mat','-ascii')
%traffic
load('capacities.mat', '-ascii')
%capacities
load('traveltime.mat', '-ascii')
%traveltime
load('flow.mat', '-ascii')
%flow

n = size(traffic,1);
m = length(flow);

c = capacities;

A_time = zeros(n,n);
A_capacity = zeros(n,n);
s = zeros(m,1);
t = zeros(m,1);
for i = 1:length(traffic)
    j = find(traffic(:,i) == 1);
    k = find(traffic(:,i) == -1);
    A_time(j,k) = traveltime(i);
    A_capacity(j,k) = capacities(i);
    s(i) = j;
    t(i) = k;
end

G_time = digraph(A_time);
G_capacity = digraph(A_capacity);
plot_grapth_2(A_capacity,s,t,c)
%plot(g)

path = shortestpath(G_time,1,17);

acc = 0;
for i = 1:length(path)-1
    acc = acc + A_time(path(i), path(i+1));
end
acc*60

maxflow(G_capacity,1,17)

traffix_flow = traffic*flow;
ext_flow = sum(traffix_flow,2)

%% Example

B = traffic;
nu = zeros(n,1);
nu(1) = ext_flow(1);
nu(n) = -ext_flow(1);
c = capacities;
zero_flow = zeros(m,1);
l = traveltime;

cvx_begin
    variables f(m)
    minimize(sum((c.*l).*inv_pos(1 - f.*inv_pos(c)) - c.*l))
    subject to
        B*f == nu
        zero_flow <= f <= c
cvx_end
f

plot_grapth_1(s,t,f,c)


%%

cvx_begin
    variable f(m)
    minimize sum(c.*l.*(1 - log(1 - f.*inv_pos(c))))
    subject to
        B*f == nu
        zero_flow <= f <= c
cvx_end
f

plot_grapth_1(s,t,f,c)



function plot_grapth_1(s,t,f,c)
    n = max(max(s,t));
    W = zeros(n,n);
    for i = 1:length(f)
        W(s(i),t(i)) = f(i);
    end
    plot_grapth_2(W,s,t,c);
end

function plot_grapth_2(W,s,t,c)
    G = digraph(W);
    f = abs(G.Edges.Weight);
    plot_flow(G,f,s,t,c);
end

function plot_flow(G, f,s,t,c)
    x = [ ...
        1, 2, 3, 4, 5, ...   % 1–5 (top row)
        1, 2, 3, ...         % 6–8
        4, ...               % 9 (center)
        2, 3, 4, ...         % 10–12
        5, 6, ...            % 13–14
        3, 4, 5 ];           % 15–17
    
    y = [ ...
        5, 5, 5, 5, 4.8, ... % 1–5
        4, 4, 4, ...         % 6–8
        4, ...               % 9
        3, 3, 3, ...         % 10–12
        3, 3, ...            % 13–14
        2, 2, 2 ];           % 15–17

    flow = f;

    labels = strings(numedges(G),1);
    
    for i = 1:numedges(G)
        edge_index = findedge(G, s(i), t(i));
        labels(edge_index) = sprintf('l_{%d} (%.1f/%i)', i, flow(edge_index), c(i));
    end
    
    figure('Position', [100, 100, 1000, 500]);
    
    p = plot(G, 'XData', x, 'YData', y);

    p.EdgeLabel = labels;
    p.MarkerSize = 20;
    %p.Marker = 'none';
    p.NodeColor = "w";

    % Overlay custom nodes with different edge colors
    hold on
    h = scatter(x, y);
    h.LineWidth = 5;
    h.MarkerFaceColor = "b";
    h.Marker = 'hexagram';
    hold off
    
    % Scale edge thickness
    p.LineWidth = 0.5 + 20 * flow / max(flow);
end