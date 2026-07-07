
function outputFinal = reconstructInterpolatedSurfaces(InterpFace, H_i, RefNUM, Corner)
% RECONSTRUCTINTERPOLATEDSURFACES Reconstruct intermediate surfaces from Beltrami data.
% Author: Joey.Chow
%
% Input:
%   InterpFace - Shared triangulation, sampled surface points, and mu values.
%   H_i        - Interpolation heights.
%   RefNUM     - Reference surface index for each interpolated layer.
%   Corner     - Boundary interpolation data.
%
% Output:
%   outputFinal - Reconstructed point coordinates for each interpolated layer.

outputFinal = zeros(length(InterpFace.Point{1}), 3, length(H_i));

for num = 1:length(H_i)
    numT = RefNUM(num);
    longestEdgeIndex = zeros(length(InterpFace.faces), 1);

    for i = 1:length(InterpFace.faces)
        length1 = norm(InterpFace.Point{numT}(InterpFace.faces(i, 2), :) - InterpFace.Point{numT}(InterpFace.faces(i, 1), :));
        length2 = norm(InterpFace.Point{numT}(InterpFace.faces(i, 3), :) - InterpFace.Point{numT}(InterpFace.faces(i, 2), :));
        length3 = norm(InterpFace.Point{numT}(InterpFace.faces(i, 1), :) - InterpFace.Point{numT}(InterpFace.faces(i, 3), :));
        [~, longestEdgeIndex(i)] = max([length1, length2, length3]);
    end

    global A_1F A_2F A_3F B_1F B_2F B_3F Bel_coeff
    v_k = zeros(length(InterpFace.faces), 2);

    for i = 1:length(InterpFace.faces)
        P1 = InterpFace.Point{numT}(InterpFace.faces(i, 1), 1:2);
        P2 = InterpFace.Point{numT}(InterpFace.faces(i, 2), 1:2);
        P3 = InterpFace.Point{numT}(InterpFace.faces(i, 3), 1:2);

        if longestEdgeIndex(i) == 2
            P_pause = P1;
            P1 = P2;
            P2 = P3;
            P3 = P_pause;
        elseif longestEdgeIndex(i) == 3
            P_pause1 = P1;
            P_pause2 = P2;
            P1 = P3;
            P2 = P_pause1;
            P3 = P_pause2;
        end

        AF = computeTriangleArea(P1, P2, P3);
        A_1F = (P2(2) - P3(2)) / AF;
        A_2F = (P3(2) - P1(2)) / AF;
        A_3F = (P1(2) - P2(2)) / AF;
        B_1F = (P3(1) - P2(1)) / AF;
        B_2F = (P1(1) - P3(1)) / AF;
        B_3F = (P2(1) - P1(1)) / AF;
        Bel_coeff = InterpFace.mu_t{num}(i, :);

        x0 = [0, 0];
        v_k(i, :) = fsolve(@evaluateBeltramiEquationResidual, x0);
    end

    V1_i = [zeros(length(InterpFace.faces), 1), zeros(length(InterpFace.faces), 1)];
    V2_j = [ones(length(InterpFace.faces), 1), zeros(length(InterpFace.faces), 1)];
    V3_k = v_k;

    for i = 1:length(InterpFace.faces)
        if longestEdgeIndex(i) == 2
            V_pause = V1_i(i, :);
            V1_i(i, :) = V3_k(i, :);
            V3_k(i, :) = V2_j(i, :);
            V2_j(i, :) = V_pause;
        elseif longestEdgeIndex(i) == 3
            V_pause = V1_i(i, :);
            V1_i(i, :) = V2_j(i, :);
            V2_j(i, :) = V3_k(i, :);
            V3_k(i, :) = V_pause;
        end
    end

    G_Fprime = zeros(6, 2, length(InterpFace.faces));
    matrixA1 = zeros(6 * length(InterpFace.faces), 2 * length(InterpFace.Point{numT}));
    matrixB = zeros(2 * length(InterpFace.Point{numT}), 1);

    H_Fprime = [-1 0 1 0 0 0;
                 0 -1 0 1 0 0;
                 0 0 -1 0 1 0;
                 0 0 0 -1 0 1;
                 1 0 0 0 -1 0;
                 0 1 0 0 0 -1];

    TransFlag2_1 = [-1 0 1 0 0 0;
                     0 -1 0 1 0 0];
    TransFlag2_2 = [0 0 -1 0 1 0;
                     0 0 0 -1 0 1];
    TransFlag2_3 = [1 0 0 0 -1 0;
                     0 1 0 0 0 -1];

    for i = 1:length(InterpFace.faces)
        T1 = [V1_i(i, 1) - V2_j(i, 1), V1_i(i, 2) - V2_j(i, 2);
              V1_i(i, 2) - V2_j(i, 2), -V1_i(i, 1) + V2_j(i, 1)];

        T2 = [V2_j(i, 1) - V3_k(i, 1), V2_j(i, 2) - V3_k(i, 2);
              V2_j(i, 2) - V3_k(i, 2), -V2_j(i, 1) + V3_k(i, 1)];

        T3 = [V3_k(i, 1) - V1_i(i, 1), V3_k(i, 2) - V1_i(i, 2);
              V3_k(i, 2) - V1_i(i, 2), -V3_k(i, 1) + V1_i(i, 1)];

        G_Fprime(:, :, i) = [V2_j(i, 1) - V1_i(i, 1), V2_j(i, 2) - V1_i(i, 2);
                             V2_j(i, 2) - V1_i(i, 2), -(V2_j(i, 1) - V1_i(i, 1));
                             V3_k(i, 1) - V2_j(i, 1), V3_k(i, 2) - V2_j(i, 2);
                             V3_k(i, 2) - V2_j(i, 2), -(V3_k(i, 1) - V2_j(i, 1));
                             V1_i(i, 1) - V3_k(i, 1), V1_i(i, 2) - V3_k(i, 2);
                             V1_i(i, 2) - V3_k(i, 2), -(V1_i(i, 1) - V3_k(i, 1))];

        Trans1 = TransFlag2_1 - T1 * inv(G_Fprime(:, :, i)' * G_Fprime(:, :, i)) * G_Fprime(:, :, i)' * H_Fprime;
        Trans2 = TransFlag2_2 - T2 * inv(G_Fprime(:, :, i)' * G_Fprime(:, :, i)) * G_Fprime(:, :, i)' * H_Fprime;
        Trans3 = TransFlag2_3 - T3 * inv(G_Fprime(:, :, i)' * G_Fprime(:, :, i)) * G_Fprime(:, :, i)' * H_Fprime;

        matrixA1(6 * i - 5, 2 * InterpFace.faces(i, 1) - 1) = Trans1(1, 1);
        matrixA1(6 * i - 5, 2 * InterpFace.faces(i, 2) - 1) = Trans1(1, 3);
        matrixA1(6 * i - 5, 2 * InterpFace.faces(i, 3) - 1) = Trans1(1, 5);
        matrixA1(6 * i - 5, 2 * InterpFace.faces(i, 1)) = Trans1(1, 2);
        matrixA1(6 * i - 5, 2 * InterpFace.faces(i, 2)) = Trans1(1, 4);
        matrixA1(6 * i - 5, 2 * InterpFace.faces(i, 3)) = Trans1(1, 6);

        matrixA1(6 * i - 4, 2 * InterpFace.faces(i, 1) - 1) = Trans1(2, 1);
        matrixA1(6 * i - 4, 2 * InterpFace.faces(i, 2) - 1) = Trans1(2, 3);
        matrixA1(6 * i - 4, 2 * InterpFace.faces(i, 3) - 1) = Trans1(2, 5);
        matrixA1(6 * i - 4, 2 * InterpFace.faces(i, 1)) = Trans1(2, 2);
        matrixA1(6 * i - 4, 2 * InterpFace.faces(i, 2)) = Trans1(2, 4);
        matrixA1(6 * i - 4, 2 * InterpFace.faces(i, 3)) = Trans1(2, 6);

        matrixA1(6 * i - 3, 2 * InterpFace.faces(i, 1) - 1) = Trans2(1, 1);
        matrixA1(6 * i - 3, 2 * InterpFace.faces(i, 2) - 1) = Trans2(1, 3);
        matrixA1(6 * i - 3, 2 * InterpFace.faces(i, 3) - 1) = Trans2(1, 5);
        matrixA1(6 * i - 3, 2 * InterpFace.faces(i, 1)) = Trans2(1, 2);
        matrixA1(6 * i - 3, 2 * InterpFace.faces(i, 2)) = Trans2(1, 4);
        matrixA1(6 * i - 3, 2 * InterpFace.faces(i, 3)) = Trans2(1, 6);

        matrixA1(6 * i - 2, 2 * InterpFace.faces(i, 1) - 1) = Trans2(2, 1);
        matrixA1(6 * i - 2, 2 * InterpFace.faces(i, 2) - 1) = Trans2(2, 3);
        matrixA1(6 * i - 2, 2 * InterpFace.faces(i, 3) - 1) = Trans2(2, 5);
        matrixA1(6 * i - 2, 2 * InterpFace.faces(i, 1)) = Trans2(2, 2);
        matrixA1(6 * i - 2, 2 * InterpFace.faces(i, 2)) = Trans2(2, 4);
        matrixA1(6 * i - 2, 2 * InterpFace.faces(i, 3)) = Trans2(2, 6);

        matrixA1(6 * i - 1, 2 * InterpFace.faces(i, 1) - 1) = Trans3(1, 1);
        matrixA1(6 * i - 1, 2 * InterpFace.faces(i, 2) - 1) = Trans3(1, 3);
        matrixA1(6 * i - 1, 2 * InterpFace.faces(i, 3) - 1) = Trans3(1, 5);
        matrixA1(6 * i - 1, 2 * InterpFace.faces(i, 1)) = Trans3(1, 2);
        matrixA1(6 * i - 1, 2 * InterpFace.faces(i, 2)) = Trans3(1, 4);
        matrixA1(6 * i - 1, 2 * InterpFace.faces(i, 3)) = Trans3(1, 6);

        matrixA1(6 * i, 2 * InterpFace.faces(i, 1) - 1) = Trans3(2, 1);
        matrixA1(6 * i, 2 * InterpFace.faces(i, 2) - 1) = Trans3(2, 3);
        matrixA1(6 * i, 2 * InterpFace.faces(i, 3) - 1) = Trans3(2, 5);
        matrixA1(6 * i, 2 * InterpFace.faces(i, 1)) = Trans3(2, 2);
        matrixA1(6 * i, 2 * InterpFace.faces(i, 2)) = Trans3(2, 4);
        matrixA1(6 * i, 2 * InterpFace.faces(i, 3)) = Trans3(2, 6);
    end

    Corner.PointIndex = Corner.BoundaryIndex;
    Corner.Corner = [Corner.boundaryX(num, :); Corner.boundaryY(num, :); Corner.boundaryZ(num, :)];
    Corner.Corner = Corner.Corner';

    conA = zeros(2 * length(Corner.PointIndex), 2 * length(InterpFace.Point{1}));
    conB = zeros(2 * length(Corner.PointIndex), 1);

    for i = 1:length(Corner.PointIndex)
        conA(2 * i - 1, 2 * Corner.PointIndex(i) - 1) = 1;
        conA(2 * i, 2 * Corner.PointIndex(i)) = 1;
        conB(2 * i - 1) = Corner.Corner(i, 1);
        conB(2 * i) = Corner.Corner(i, 2);
    end

    opts = optimset('Algorithm', 'interior-point-convex');
    tem_output = quadprog(matrixA1' * matrixA1, matrixB, [], [], conA, conB, [], [], [], opts);

    output = zeros(length(InterpFace.Point{numT}), 2);
    for i = 1:length(InterpFace.Point{numT})
        output(i, 1) = tem_output(2 * i - 1, 1);
        output(i, 2) = tem_output(2 * i, 1);
    end

    outputFinal(:, :, num) = [output, zeros(length(InterpFace.Point{1}), 1)];
end

end
