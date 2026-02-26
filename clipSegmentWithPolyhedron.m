function [P1c, P2c, lenInside] = clipSegmentWithPolyhedron(P1, P2, verts, faces)

dir = P2 - P1;
t0 = 0; t1 = 1;

TR = triangulation(faces, verts);
N  = faceNormal(TR);

for k = 1:size(faces,1)
    faceIdx = faces(k,:);
    pt = verts(faceIdx(1),:);
    n  = N(k,:);

    num = dot(n, pt - P1);
    den = dot(n, dir);

    if abs(den) < 1e-12
        if num < 0
            P1c=P1; P2c=P1; lenInside=0; return;
        end
        continue
    end

    t = num / den;

    if den > 0
        t0 = max(t0, t);
    else
        t1 = min(t1, t);
    end

    if t0 > t1
        P1c=P1; P2c=P1; lenInside=0;
        return;
    end
end

P1c = P1 + t0*dir;
P2c = P1 + t1*dir;
lenInside = norm(P2c - P1c);

end
