function P_2d = equation_parts(n_voxel,neighbor_id_2d,simulation_input)

% neighbor indices
delta_2d = sign(neighbor_id_2d); % the indicating variable delta, and it is a 2D variable

neighbor_id_2d_2 = neighbor_id_2d; % Matlab cannot index 0, so 0 is replaced with 1, 
neighbor_id_2d_2(neighbor_id_2d_2==0) = 1; % and this does not matter because the associated delta_2d will be 0

% equation parts
D = simulation_input.D0;
D11 = zeros(n_voxel,1);
D12 = zeros(n_voxel,1);
D13 = zeros(n_voxel,1);
D21 = zeros(n_voxel,1);
D22 = zeros(n_voxel,1);
D23 = zeros(n_voxel,1);
D31 = zeros(n_voxel,1);
D32 = zeros(n_voxel,1);
D33 = zeros(n_voxel,1);
for n = 1:n_voxel
    D11(n) = D{n}(1,1);
    D12(n) = D{n}(1,2);
    D13(n) = D{n}(1,3);
    D21(n) = D{n}(2,1);
    D22(n) = D{n}(2,2);
    D23(n) = D{n}(2,3);
    D31(n) = D{n}(3,1);
    D32(n) = D{n}(3,2);
    D33(n) = D{n}(3,3);
end

P_2d = zeros(n_voxel,21); % parts of the equation, and it is a 2D variable
P_2d(:,1) = 4*delta_2d(:,1).*D11;
P_2d(:,2) = 4*delta_2d(:,2).*D11;
P_2d(:,3) = 4*delta_2d(:,3).*D22;
P_2d(:,4) = 4*delta_2d(:,4).*D22;
P_2d(:,5) = 4*delta_2d(:,5).*D33;
P_2d(:,6) = 4*delta_2d(:,6).*D33;
P_2d(:,7) = delta_2d(:,1).*delta_2d(:,2) .* ...
                ( delta_2d(:,1).*delta_2d(:,2).*(D11(neighbor_id_2d_2(:,1))-D11(neighbor_id_2d_2(:,2))) + ...
                delta_2d(:,3).*delta_2d(:,4).*(D21(neighbor_id_2d_2(:,3))-D21(neighbor_id_2d_2(:,4))) + ...
                delta_2d(:,5).*delta_2d(:,6).*(D31(neighbor_id_2d_2(:,5))-D31(neighbor_id_2d_2(:,6))) );
P_2d(:,8) = delta_2d(:,3).*delta_2d(:,4) .* ...
                ( delta_2d(:,1).*delta_2d(:,2).*(D12(neighbor_id_2d_2(:,1))-D12(neighbor_id_2d_2(:,2))) + ...
                delta_2d(:,3).*delta_2d(:,4).*(D22(neighbor_id_2d_2(:,3))-D22(neighbor_id_2d_2(:,4))) + ...
                delta_2d(:,5).*delta_2d(:,6).*(D32(neighbor_id_2d_2(:,5))-D32(neighbor_id_2d_2(:,6))) );
P_2d(:,9) = delta_2d(:,5).*delta_2d(:,6) .* ...
                ( delta_2d(:,1).*delta_2d(:,2).*(D13(neighbor_id_2d_2(:,1))-D13(neighbor_id_2d_2(:,2))) + ...
                delta_2d(:,3).*delta_2d(:,4).*(D23(neighbor_id_2d_2(:,3))-D23(neighbor_id_2d_2(:,4))) + ...
                delta_2d(:,5).*delta_2d(:,6).*(D33(neighbor_id_2d_2(:,5))-D33(neighbor_id_2d_2(:,6))) );
P_2d(:,10) = 2*delta_2d(:,7).*delta_2d(:,9).*D12;
P_2d(:,11) = 2*delta_2d(:,8).*delta_2d(:,10).*D12;
P_2d(:,12) = 2*delta_2d(:,15).*delta_2d(:,17).*D13;
P_2d(:,13) = 2*delta_2d(:,16).*delta_2d(:,18).*D13;
P_2d(:,14) = 2*delta_2d(:,11).*delta_2d(:,13).*D23;
P_2d(:,15) = 2*delta_2d(:,12).*delta_2d(:,14).*D23;
P_2d(:,16) = simulation_input.tau_open_voxel;
P_2d(:,17) = simulation_input.tau_close_voxel;
P_2d(:,18) = simulation_input.tau_in_voxel;
P_2d(:,19) = simulation_input.tau_out_voxel;
P_2d(:,20) = simulation_input.v_gate_voxel;
P_2d(:,21) = simulation_input.c_voxel;

end
