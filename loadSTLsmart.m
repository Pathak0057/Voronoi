function [V, F] = loadSTLsmart(filePath)
% FINAL universal STL reader supporting:
% - ASCII STL
% - Binary STL (MATLAB STL 2018 header)
% - Struct formats

fid = fopen(filePath, 'r');
if fid < 0
    error('Cannot open STL file: %s', filePath);
end

% Read first 256 bytes to check if ASCII or binary
header = fread(fid, 256, 'uint8=>char')';
frewind(fid);

% If header contains "solid" but not "facet normal", treat as Binary
isASCII = startsWith(strtrim(header), 'solid') && contains(header, 'facet');

if isASCII
    % ---------------- ASCII STL ----------------
    [F, V] = READ_ascii_stl(fid);
else
    % ---------------- BINARY STL ----------------
    [F, V] = READ_binary_stl(fid);
end

fclose(fid);

% Force double precision
V = double(V);
F = double(F);

end

%% ===============================================================
%% =============== ASCII STL READER ==============================
%% ===============================================================
function [F, V] = READ_ascii_stl(fid)

C = textscan(fid, '%s', 'Delimiter', '\n');
C = C{1};

V = [];
F = [];

for i = 1:length(C)
    line = strtrim(C{i});
    if startsWith(line, 'vertex')
        nums = sscanf(line(7:end), '%f %f %f');
        V(end+1,:) = nums';
    end
end

% STL ASCII repeats 3 vertices per triangle
nTri = size(V,1) / 3;
F = reshape(1:(3*nTri), 3, nTri)' ;

end

%% ===============================================================
%% =============== BINARY STL READER =============================
%% ===============================================================
function [F, V] = READ_binary_stl(fid)

% Skip header (80 bytes)
fread(fid, 80, 'uint8');

% Number of triangles
nTri = fread(fid, 1, 'uint32');

V = zeros(nTri*3, 3);
F = zeros(nTri, 3);

for i = 1:nTri
    % normal vector (ignored)
    fread(fid, 3, 'float32');

    % vertices
    v1 = fread(fid, 3, 'float32')';
    v2 = fread(fid, 3, 'float32')';
    v3 = fread(fid, 3, 'float32')';

    V(3*(i-1)+1,:) = v1;
    V(3*(i-1)+2,:) = v2;
    V(3*(i-1)+3,:) = v3;

    F(i,:) = [3*(i-1)+1, 3*(i-1)+2, 3*(i-1)+3];

    % attribute byte count (ignored)
    fread(fid, 1, 'uint16');
end

end
