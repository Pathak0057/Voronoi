clc; clear; close all;

excelFile      = 'final_fiber_data_CentrelinePoints_Abaqus_static_S1_AF.xlsx';    % Load excel file with fiber coordinates

grid_min       = -70;                                           % These grids are used sphere movement 
grid_max       =  70;
grid_step      =  50;
sphereRadius   =  50;                                        % Sphere radius will determin number of voronoi cells
cube_min       = -50;                                           % These cube dimension is of RVE
cube_max       =  50;

neighborRadius = 80;

outputFolder   = 'VoronoiSTL';
if ~exist(outputFolder,'dir')
    mkdir(outputFolder);
end

outputCSV = 'voronoi_fiber_summary.csv';                      % The output csv records the fiber cooridnates for each voronoi cell with Vf.

tolDet  = 1e-10;
tolIneq = 1e-7;
snapTol = 1e-6;

%% ---------------- FIGURE OPTIONS ---------------------------------------
faceAlpha         = 0.5;
showEdges         = false;
maxPlotCells      = inf;

fiberLineWidth    = 2.0;
clipFibersToCube  = true;

fiberColorMode    = "byCell";     % "byCell" or "uniform"
fiberUniformColor = [0.85 0.1 0.1];

%% ---------------- VISUALIZATION (CYLINDERS + GIF) ----------------------
fiberVisRadius      = 1.2;     % <<< visualization-only cylinder radius (units)
fiberVisNCirc       = 18;      % cylinder smoothness (higher = smoother, slower)

enable_gif_generation = true;  % set false to skip GIF creation
gif_filename          = 'Voronoi_model_sphere_sliding.gif';
gif_delay_time        = 0.15;  % seconds per frame
gif_max_frames        = inf;   % set e.g. 50 to limit frames for faster export

sphereRGBA            = [0 0 1 0.25]; % blue + transparency
fiberRGBA             = [1 0 0 1.0];  % red

%% ---------------- REGION FIBER VOLUME PARAMS ---------------------------
r_fiber          = 4;     % fiber radius
insideTolPlane   = 1e-9;  % tolerance for convex clipping

%% ---------------- LOAD FIBER ENDPOINTS ---------------------------------
tbl     = readtable(excelFile, 'VariableNamingRule','preserve');
splines = unique(tbl.Spline_name);
nFibers = numel(splines);

P1s = zeros(nFibers,3);
P2s = zeros(nFibers,3);

for k = 1:nFibers
    data = tbl(strcmp(tbl.Spline_name, splines{k}), :);
    pts  = [data.Position_x_global, data.Position_y, data.Position_z];
    P1s(k,:) = pts(1,:);
    P2s(k,:) = pts(end,:);
end
fprintf('Loaded %d fibers.\n', nFibers);

%% ---------------- BUILD SEED GRID --------------------------------------
[gx, gy, gz] = ndgrid(grid_min:grid_step:grid_max, ...
                      grid_min:grid_step:grid_max, ...
                      grid_min:grid_step:grid_max);
sphereCenters = [gx(:), gy(:), gz(:)];
nSeeds_grid   = size(sphereCenters,1);

seeds = zeros(nSeeds_grid,3);
for i = 1:nSeeds_grid
    seeds(i,:) = computeSeed(sphereCenters(i,:), sphereRadius, P1s, P2s);
end

% REMOVE DUPLICATES
[seeds, ia] = unique(seeds, 'rows', 'stable');
sphereCenters = sphereCenters(ia,:);
nSeeds = size(seeds,1);

%% ---------------- PRECOMPUTE SEED DISTANCES ----------------------------
D = pdist2(seeds, seeds);

%% ---------------- HALF-SPACES FOR CUBE ---------------------------------
A_cube = [
     1  0  0;
    -1  0  0;
     0  1  0;
     0 -1  0;
     0  0  1;
     0  0 -1];
b_cube = [
    cube_max;
   -cube_min;
    cube_max;
   -cube_min;
    cube_max;
   -cube_min];

%% ---------------- INITIALIZE STORAGE -----------------------------------
allVertices = cell(nSeeds,1);
allFaces    = cell(nSeeds,1);

%% ---------------- BUILD VORONOI CELLS ----------------------------------
for i = 1:nSeeds
    si = seeds(i,:);
    neighIdx = find(D(i,:) > 0 & D(i,:) <= neighborRadius);

    if isempty(neighIdx) && any(si < cube_min | si > cube_max)
        continue;
    end

    A = A_cube;
    b = b_cube;

    for j = neighIdx
        sj = seeds(j,:);
        n_ij = sj - si;
        c_ij = 0.5*(dot(sj,sj) - dot(si,si));
        A = [A; n_ij];
        b = [b; c_ij];
    end

    mPlanes = size(A,1);
    verts = [];

    for p = 1:mPlanes-2
        Ap = A(p,:);
        for q = p+1:mPlanes-1
            Aq = A(q,:);
            for r = q+1:mPlanes
                Ar = A(r,:);
                M = [Ap; Aq; Ar];

                if abs(det(M)) < tolDet
                    continue;
                end

                rhs = [b(p); b(q); b(r)];
                x = M\rhs;

                if all(A*x <= b + tolIneq)
                    verts = [verts; x'];
                end
            end
        end
    end

    if isempty(verts), continue; end

    verts = round(verts / snapTol) * snapTol;
    verts = unique(verts, 'rows', 'stable');
    if size(verts,1) < 4, continue; end

    try
        F = convhulln(verts);
    catch
        continue;
    end

    allVertices{i} = verts;
    allFaces{i}    = F;

    % Export STL
    TR = triangulation(F, verts);
    stlName = fullfile(outputFolder, sprintf('region_%03d.stl', i));
    stlwrite(TR, stlName);
end

save('voronoi_generated.mat','allVertices','allFaces','seeds','P1s','P2s','splines');
fprintf('Done. Saved geometry.\n');

%% ------------------------------------------------------------------------
% CSV SUMMARY OVER ALL REGIONS (This CSV data is used for discrete
% modelling later on)
%% ------------------------------------------------------------------------
RegionID           = (1:nSeeds)';  % region_001 corresponds to i=1, etc.
VoronoiVolume      = zeros(nSeeds,1);
TotalFiberVolume   = zeros(nSeeds,1);
FiberVolumeFrac    = zeros(nSeeds,1);

FiberIndices       = strings(nSeeds,1);
FiberLengthsInside = strings(nSeeds,1);
FiberP1_inside     = strings(nSeeds,1);
FiberP2_inside     = strings(nSeeds,1);

MinCorner          = strings(nSeeds,1);
CuboidSize         = strings(nSeeds,1);

fprintf('\nComputing per-region fiber summary for CSV...\n');

for i = 1:nSeeds
    V = allVertices{i};
    F = allFaces{i};
    if isempty(V) || isempty(F)
        VoronoiVolume(i)      = 0;
        TotalFiberVolume(i)   = 0;
        FiberVolumeFrac(i)    = 0;
        FiberIndices(i)       = "";
        FiberLengthsInside(i) = "";
        FiberP1_inside(i)     = "";
        FiberP2_inside(i)     = "";
        MinCorner(i)          = "";
        CuboidSize(i)         = "";
        continue;
    end

    mn = min(V,[],1);
    mx = max(V,[],1);
    MinCorner(i) = sprintf('[%.3f %.3f %.3f]', mn);
    CuboidSize(i)= sprintf('[%.3f %.3f %.3f]', mx - mn);

    try
        [~, vol] = convhulln(V);
    catch
        vol = 0;
    end
    VoronoiVolume(i) = vol;

    centroidCell = mean(V,1);
    planes = buildConvexHalfspacesFromFaces(V, F, centroidCell);

    f_ids     = [];
    f_lengths = [];
    fP1s      = {};
    fP2s      = {};
    totalVol  = 0;

    for j = 1:nFibers
        A = P1s(j,:); B = P2s(j,:);
        [hit, Aclip, Bclip, lenInside] = clipSegmentWithConvexHalfspaces(A, B, planes, insideTolPlane);

        if hit && lenInside > 0
            f_ids(end+1)     = j;          %#ok<AGROW>
            f_lengths(end+1) = lenInside;  %#ok<AGROW>

            fP1s{end+1} = sprintf('[%.4f %.4f %.4f]', Aclip); %#ok<AGROW>
            fP2s{end+1} = sprintf('[%.4f %.4f %.4f]', Bclip); %#ok<AGROW>

            totalVol = totalVol + pi*r_fiber*r_fiber * lenInside;
        end
    end

    TotalFiberVolume(i) = totalVol;
    FiberVolumeFrac(i)  = totalVol / max(vol, 1e-12);

    FiberIndices(i)        = strjoin(string(f_ids), ', ');
    FiberLengthsInside(i)  = strjoin(string(round(f_lengths,4)), ', ');
    FiberP1_inside(i)      = strjoin(string(fP1s), ', ');
    FiberP2_inside(i)      = strjoin(string(fP2s), ', ');
end

T = table(RegionID, VoronoiVolume, TotalFiberVolume, FiberVolumeFrac, ...
          FiberIndices, FiberLengthsInside, FiberP1_inside, FiberP2_inside, ...
          MinCorner, CuboidSize);

writetable(T, outputCSV);
fprintf('\nSaved CSV summary: %s\n', outputCSV);

%% ========================= FIGURE (1): ALL FIBERS ======================
figure('Color','w'); hold on; axis equal; grid on; view(3);
xlabel('X'); ylabel('Y'); zlabel('Z');
title('Figure (1): All Fibers');
drawCubeWireframe(cube_min, cube_max);

for k = 1:nFibers
    A = P1s(k,:); B = P2s(k,:);
    if clipFibersToCube
        [ok, A, B] = clipSegmentToAABB(A, B, cube_min, cube_max);
        if ~ok, continue; end
    end
    drawCylinderSegment(A, B, fiberVisRadius, fiberVisNCirc, fiberRGBA);
end
camlight headlight; lighting gouraud; set(gca,'Box','on');

%% ===================== GIF: SPHERE SLIDING THROUGH FIBERS ===============
if enable_gif_generation

    figGif = figure('Color','w'); clf(figGif);
    hold on; view(3); grid on;

    xlabel('X'); ylabel('Y'); zlabel('Z');
    title('GIF: Sphere sliding through fibers');

    ax = gca;
    ax.Box = 'on';
    try
        ax.BoxStyle = 'full';   
    catch
        box on                 
    end
    ax.LineWidth = 1.2;
    ax.TickDir = 'out';
    ax.TickLength = [0.015 0.015];
    axis equal

    %% ---- Cuboid faces ----
    drawCubeFaces(-50, 50, [0.85 0.85 0.85], 0.08, 'k', 1.0);

    %% ---- Lock axis & camera ----
    pad = 20;
    xlim([cube_min-pad, cube_max+pad]);
    ylim([cube_min-pad, cube_max+pad]);
    zlim([cube_min-pad, cube_max+pad]);
    axis manual
    daspect([1 1 1]);

    axG = gca;
    view(3);
    camproj('perspective');
    campos0    = campos(axG);
    camtarget0 = camtarget(axG);
    camup0     = camup(axG);
    camva0     = camva(axG);

    %% ---- Draw ALL fibers once ----
    for k = 1:nFibers
        A = P1s(k,:); B = P2s(k,:);
        if clipFibersToCube
            [ok, A, B] = clipSegmentToAABB(A, B, cube_min, cube_max);
            if ~ok, continue; end
        end
        drawCylinderSegment(A, B, fiberVisRadius, fiberVisNCirc, fiberRGBA);
    end

    %% ---- Sphere animation ----
    hSphere = [];
    nFrames = min(size(sphereCenters,1), gif_max_frames);

    for iF = 1:nFrames
        c0 = sphereCenters(iF,:);

        if ~isempty(hSphere) && isvalid(hSphere)
            delete(hSphere);
        end
        hSphere = drawSimpleSphere(c0, sphereRadius, sphereRGBA);

        % keep camera fixed
        campos(axG, campos0);
        camtarget(axG, camtarget0);
        camup(axG, camup0);
        camva(axG, camva0);

        drawnow;

        frame = getframe(gcf);
        [im, cmap] = rgb2ind(frame.cdata,256);

        if iF==1
            imwrite(im,cmap,gif_filename,'gif','LoopCount',Inf,'DelayTime',gif_delay_time);
        else
            imwrite(im,cmap,gif_filename,'gif','WriteMode','append','DelayTime',gif_delay_time);
        end
    end

    fprintf('GIF saved: %s\n', gif_filename);
end

%% ===================== FIGURE (2): VORONOI ONLY ========================
figure('Color','w'); hold on; axis equal; grid on; view(3);
xlabel('X'); ylabel('Y'); zlabel('Z');
title('Figure (2): Voronoi Regions Only');
drawCubeWireframe(cube_min, cube_max);
scatter3(seeds(:,1), seeds(:,2), seeds(:,3), 18, 'k', 'filled');

plotted = 0;
cellColors = nan(nSeeds,3);

for i = 1:nSeeds
    V = allVertices{i}; F = allFaces{i};
    if isempty(V) || isempty(F), continue; end

    plotted = plotted + 1;
    if plotted > maxPlotCells, break; end

    rng(i); c = rand(1,3);
    cellColors(i,:) = c;

    if showEdges
        edgeCol = [0 0 0]; edgeA = 0.1;
    else
        edgeCol = 'none'; edgeA = 0;
    end

    patch('Vertices', V, 'Faces', F, ...
        'FaceColor', c, 'FaceAlpha', faceAlpha, ...
        'EdgeColor', edgeCol, 'EdgeAlpha', edgeA);
end
camlight headlight; lighting gouraud; set(gca,'Box','on');

%% =========== FIGURE (3): VORONOI + FIBER OVERLAY ========================
figure('Color','w'); hold on; axis equal; grid on; view(3);
xlabel('X'); ylabel('Y'); zlabel('Z');
title('Figure (3): Voronoi Regions + Fiber Overlay');
drawCubeWireframe(cube_min, cube_max);
scatter3(seeds(:,1), seeds(:,2), seeds(:,3), 18, 'k', 'filled');

plotted = 0;
for i = 1:nSeeds
    V = allVertices{i}; F = allFaces{i};
    if isempty(V) || isempty(F), continue; end

    plotted = plotted + 1;
    if plotted > maxPlotCells, break; end

    c = cellColors(i,:);
    if any(isnan(c)), rng(i); c = rand(1,3); end

    if showEdges
        edgeCol = [0 0 0]; edgeA = 0.1;
    else
        edgeCol = 'none'; edgeA = 0;
    end

    patch('Vertices', V, 'Faces', F, ...
        'FaceColor', c, 'FaceAlpha', faceAlpha, ...
        'EdgeColor', edgeCol, 'EdgeAlpha', edgeA);
end

midPts = 0.5*(P1s + P2s);
idxCell = knnsearch(seeds, midPts);

for k = 1:nFibers
    A = P1s(k,:); B = P2s(k,:);
    if clipFibersToCube
        [ok, A, B] = clipSegmentToAABB(A, B, cube_min, cube_max);
        if ~ok, continue; end
    end

    if fiberColorMode == "byCell"
        ci = idxCell(k);
        c  = cellColors(ci,:);
        if any(isnan(c)), c = [0 0 0]; end
    else
        c = fiberUniformColor;
    end

    plot3([A(1) B(1)], [A(2) B(2)], [A(3) B(3)], '-', ...
        'LineWidth', fiberLineWidth, 'Color', c);
end

camlight headlight; lighting gouraud; set(gca,'Box','on');

%% =========== FIGURE (4): SPECIFIC REGION + TABLES BESIDE ===============
regionIdx = input('Enter Voronoi region index (e.g., 1 for region_001): ');

if isempty(regionIdx) || regionIdx < 1 || regionIdx > nSeeds || isempty(allVertices{regionIdx})
    warning('Invalid region index or region is empty. Skipping Figure (4).');
else
    V = allVertices{regionIdx};
    F = allFaces{regionIdx};

    % Pull values from CSV-summary arrays (already computed)
    VorVol    = VoronoiVolume(regionIdx);
    TotFibVol = TotalFiberVolume(regionIdx);
    FVF       = FiberVolumeFrac(regionIdx);

    % Recompute detailed per-fiber clipped endpoints/lengths (to show in tables)
    centroidCell = mean(V,1);
    planes = buildConvexHalfspacesFromFaces(V, F, centroidCell);

    f_ids     = [];
    f_lengths = [];
    fP1_inside = [];
    fP2_inside = [];

    for j = 1:nFibers
        A = P1s(j,:); B = P2s(j,:);
        [hit, Aclip, Bclip, lenInside] = clipSegmentWithConvexHalfspaces(A, B, planes, insideTolPlane);
        if hit && lenInside > 0
            f_ids(end+1,1)     = j;          %#ok<AGROW>
            f_lengths(end+1,1) = lenInside;  %#ok<AGROW>
            fP1_inside(end+1,:)= Aclip;      %#ok<AGROW>
            fP2_inside(end+1,:)= Bclip;      %#ok<AGROW>
        end
    end

    fig4 = figure('Color','w', 'Name', sprintf('Region %03d detail', regionIdx));
    fig4.Units = 'normalized'; fig4.Position = [0.05 0.08 0.90 0.82];

    ax = axes(fig4, 'Units','normalized', 'Position',[0.05 0.10 0.55 0.85]);
    hold(ax,'on'); axis(ax,'equal'); grid(ax,'on'); view(ax,3);
    xlabel(ax,'X'); ylabel(ax,'Y'); zlabel(ax,'Z');
    title(ax, sprintf('Figure (4): region%03d (Cell + Fibers Inside)', regionIdx));
    drawCubeWireframe(cube_min, cube_max);

    patch('Vertices', V, 'Faces', F, ...
        'FaceColor', [0.2 0.6 0.9], 'FaceAlpha', 0.55, ...
        'EdgeColor', 'k', 'LineWidth', 0.8, 'Parent', ax);

    scatter3(ax, seeds(regionIdx,1), seeds(regionIdx,2), seeds(regionIdx,3), 70, 'k', 'filled');

    for ii = 1:numel(f_ids)
        Aclip = fP1_inside(ii,:);
        Bclip = fP2_inside(ii,:);
        if clipFibersToCube
            [ok, Aplot, Bplot] = clipSegmentToAABB(Aclip, Bclip, cube_min, cube_max);
            if ~ok, Aplot=Aclip; Bplot=Bclip; end
        else
            Aplot=Aclip; Bplot=Bclip;
        end
        plot3(ax, [Aplot(1) Bplot(1)], [Aplot(2) Bplot(2)], [Aplot(3) Bplot(3)], '-', ...
            'LineWidth', 2.6, 'Color', [0.85 0.1 0.1]);
    end

    camlight(ax,'headlight'); lighting(ax,'gouraud'); set(ax,'Box','on');

    % Right-side tables
    summaryData = {
        'RegionID',           regionIdx;
        'VoronoiVolume',      VorVol;
        'TotalFiberVolume',   TotFibVol;
        'FiberVolumeFrac',    FVF;
        'NumFibersInside',    numel(f_ids);
        'FiberRadius',        r_fiber;
        'CSV_Output',         outputCSV
        };

    uitable(fig4, 'Units','normalized', 'Position',[0.63 0.78 0.34 0.20], ...
        'Data', summaryData, 'ColumnName', {'Property','Value'}, ...
        'RowName', [], 'FontSize', 10);

    uitable(fig4, 'Units','normalized', 'Position',[0.63 0.58 0.16 0.18], ...
        'Data', num2cell(f_ids), 'ColumnName', {'FiberIndices'}, ...
        'RowName', [], 'FontSize', 10);

    uitable(fig4, 'Units','normalized', 'Position',[0.81 0.58 0.16 0.18], ...
        'Data', num2cell(f_lengths), 'ColumnName', {'FiberLengthsInside'}, ...
        'RowName', [], 'FontSize', 10);

    P1_cell = arrayfun(@(k) sprintf('[%.4f %.4f %.4f]', fP1_inside(k,1), fP1_inside(k,2), fP1_inside(k,3)), ...
                       (1:size(fP1_inside,1))', 'UniformOutput', false);
    P2_cell = arrayfun(@(k) sprintf('[%.4f %.4f %.4f]', fP2_inside(k,1), fP2_inside(k,2), fP2_inside(k,3)), ...
                       (1:size(fP2_inside,1))', 'UniformOutput', false);

    uitable(fig4, 'Units','normalized', 'Position',[0.63 0.34 0.34 0.22], ...
        'Data', P1_cell, 'ColumnName', {'FiberP1_inside'}, ...
        'RowName', [], 'FontSize', 10);

    uitable(fig4, 'Units','normalized', 'Position',[0.63 0.10 0.34 0.22], ...
        'Data', P2_cell, 'ColumnName', {'FiberP2_inside'}, ...
        'RowName', [], 'FontSize', 10);

    fprintf('Region %03d:\n  VoronoiVolume = %.6g\n  TotalFiberVolume = %.6g\n  FiberVolumeFrac = %.6g\n  FibersInside = %d\n', ...
        regionIdx, VorVol, TotFibVol, FVF, numel(f_ids));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function centroid = computeSeed(sc, R, P1s, P2s)
fiber_centroids = [];
for j = 1:size(P1s,1)
    [hit, t1, t2] = sphereSegmentIntersection(sc, R, P1s(j,:), P2s(j,:));
    if hit
        A = P1s(j,:) + t1*(P2s(j,:)-P1s(j,:));
        B = P1s(j,:) + t2*(P2s(j,:)-P1s(j,:));
        fiber_centroids(end+1,:) = (A + B)/2; %#ok<AGROW>
    end
end
if isempty(fiber_centroids)
    centroid = sc;
else
    centroid = mean([fiber_centroids; sc],1);
end
end

function [hit,t1,t2] = sphereSegmentIntersection(sc, R, A, B)
AB = B - A;
a  = dot(AB,AB);
if a < 1e-14
    hit=false; t1=0; t2=0; return;
end
AS = A - sc;
b  = 2*dot(AB,AS);
c  = dot(AS,AS) - R^2;

D = b^2 - 4*a*c;
if D < 0
    hit=false; t1=0; t2=0; return;
end

sd = sqrt(D);
u1 = (-b - sd)/(2*a);
u2 = (-b + sd)/(2*a);

t1 = max(0, min(1, min(u1,u2)));
t2 = max(0, min(1, max(u1,u2)));

hit = (t1 < t2);
end

function drawCubeWireframe(cmin, cmax)
P = [cmin cmin cmin;
     cmax cmin cmin;
     cmax cmax cmin;
     cmin cmax cmin;
     cmin cmin cmax;
     cmax cmin cmax;
     cmax cmax cmax;
     cmin cmax cmax];

E = [1 2; 2 3; 3 4; 4 1;
     5 6; 6 7; 7 8; 8 5;
     1 5; 2 6; 3 7; 4 8];

for k = 1:size(E,1)
    p1 = P(E(k,1),:);
    p2 = P(E(k,2),:);
    plot3([p1(1) p2(1)], [p1(2) p2(2)], [p1(3) p2(3)], 'k-', 'LineWidth', 1.2);
end
end

function [ok, A2, B2] = clipSegmentToAABB(A, B, bmin, bmax)
d = B - A;
t0 = 0; t1 = 1;

for axis = 1:3
    if abs(d(axis)) < 1e-14
        if A(axis) < bmin || A(axis) > bmax
            ok = false; A2=A; B2=B; return;
        end
    else
        invD = 1/d(axis);
        tNear = (bmin - A(axis))*invD;
        tFar  = (bmax - A(axis))*invD;
        if tNear > tFar
            tmp = tNear; tNear = tFar; tFar = tmp;
        end
        t0 = max(t0, tNear);
        t1 = min(t1, tFar);
        if t0 > t1
            ok = false; A2=A; B2=B; return;
        end
    end
end

A2 = A + t0*d;
B2 = A + t1*d;
ok = true;
end

function planes = buildConvexHalfspacesFromFaces(V, F, centroid)
nF = size(F,1);
planes = zeros(nF,4);

for i = 1:nF
    p1 = V(F(i,1),:);
    p2 = V(F(i,2),:);
    p3 = V(F(i,3),:);

    n = cross(p2-p1, p3-p1);
    nn = norm(n);
    if nn < 1e-14
        continue;
    end
    n = n/nn;
    d = dot(n, p1);

    if dot(n, centroid) > d
        n = -n;
        d = -d;
    end

    planes(i,:) = [n d];
end

planes = planes(any(abs(planes(:,1:3))>0,2),:);
end

function [hit, Aclip, Bclip, lenInside] = clipSegmentWithConvexHalfspaces(A, B, planes, tol)
dvec = B - A;
t0 = 0; t1 = 1;

for i = 1:size(planes,1)
    n = planes(i,1:3);
    d = planes(i,4);

    num = d - dot(n, A);
    den = dot(n, dvec);

    if abs(den) < 1e-14
        if dot(n, A) > d + tol
            hit = false; Aclip=A; Bclip=B; lenInside=0; return;
        else
            continue;
        end
    end

    t = num / den;

    if den > 0
        t1 = min(t1, t);
    else
        t0 = max(t0, t);
    end

    if t0 > t1
        hit = false; Aclip=A; Bclip=B; lenInside=0; return;
    end
end

t0 = max(0, min(1, t0));
t1 = max(0, min(1, t1));

if t0 >= t1
    hit = false; Aclip=A; Bclip=B; lenInside=0; return;
end

Aclip = A + t0*dvec;
Bclip = A + t1*dvec;
lenInside = norm(Bclip - Aclip);
hit = (lenInside > 0);
end


%% ===================== ADDED HELPERS (CYLINDER + SPHERE) ===============
function h = drawSimpleSphere(center, radius, rgba)
% Draw a semi-transparent sphere.
% rgba = [r g b a]
if numel(rgba) < 4, rgba(4) = 0.25; end
[xs, ys, zs] = sphere(24);
xs = xs * radius + center(1);
ys = ys * radius + center(2);
zs = zs * radius + center(3);

h = surf(xs, ys, zs, ...
    'FaceColor', rgba(1:3), ...
    'FaceAlpha', rgba(4), ...
    'EdgeColor', 'none');
end

function drawCylinderSegment(p1, p2, radius, nCirc, rgba)
% Draw a cylinder between p1 and p2 with given radius.
% rgba = [r g b a]
if numel(rgba) < 4, rgba(4) = 1.0; end

v = p2 - p1;
L = norm(v);
if L < 1e-12
    return;
end
vhat = v / L;

% Build an orthonormal basis with vhat as the cylinder axis
% Pick a vector not parallel to vhat
if abs(vhat(3)) < 0.9
    a = [0 0 1];
else
    a = [0 1 0];
end
u = cross(vhat, a); u = u / norm(u);
w = cross(vhat, u);

theta = linspace(0, 2*pi, nCirc);
z = [0, L];

[tt, zz] = meshgrid(theta, z);
X = radius * cos(tt);
Y = radius * sin(tt);

% Map to 3D: p = p1 + u*X + w*Y + vhat*Z
Px = p1(1) + u(1)*X + w(1)*Y + vhat(1)*zz;
Py = p1(2) + u(2)*X + w(2)*Y + vhat(2)*zz;
Pz = p1(3) + u(3)*X + w(3)*Y + vhat(3)*zz;

surf(Px, Py, Pz, ...
    'FaceColor', rgba(1:3), ...
    'FaceAlpha', rgba(4), ...
    'EdgeColor', 'none');
end

function h = drawCubeFaces(cube_min, cube_max, faceColor, faceAlpha, edgeColor, edgeLW)

xmin = cube_min; xmax = cube_max;

verts = [
xmin xmin xmin
xmax xmin xmin
xmax xmax xmin
xmin xmax xmin
xmin xmin xmax
xmax xmin xmax
xmax xmax xmax
xmin xmax xmax];

faces = [
1 2 3 4
5 6 7 8
1 2 6 5
2 3 7 6
3 4 8 7
4 1 5 8];

h = patch('Vertices',verts,'Faces',faces,...
    'FaceColor',faceColor,'FaceAlpha',faceAlpha,...
    'EdgeColor',edgeColor,'LineWidth',edgeLW);
end