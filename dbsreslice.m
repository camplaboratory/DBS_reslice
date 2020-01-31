function dbsreslice(data)
% function dbsreslice --- this is a rebuild of hippotaxy(data), to reslice
% dbs shank and sites into a pseudocoronal plane
%     Called by function kjm_DBS_parse. Part of the "Hippotaxy" tool.
%     This is a program for rotating and reslicing brain imaging data into
%     Hippocampal stereotactic space, as described by the manuscript:
%     "Hippocampal stereotaxy: A novel mesial temporal stereotactic
%     coordinate system", by Kai Miller and colleagues, and is currently in
%     submission. Please cite this manuscript in any setting (manuscripts,
%     talks) where this program was used. 
%     
%     Copyright (C) 2015, Kai J Miller, Stanford Neurosurgery
%     kai.miller@stanford.edu, kjmiller@gmail.com


%% ac-pc
    acpc_vec=(data.AC-data.PC)/(sum((data.AC-data.PC).^2).^.5); %unit norm ac-pc vector
    acpc_mcp=(data.AC+data.PC)/2; % mid-commisural point
    data.brain_info.mat(1:3,4)=-data.brain_info.mat(1:3,1:3)*acpc_mcp.'; % re-center image on mcp

%% plane of midline symmetry
    MLPlanepts=[data.AC; data.PC; data.ML1; data.ML2; data.ML3];
    [pc_vecs,tmp,pc_vals]=pca(MLPlanepts); 
    [tmp,n]=min(abs(pc_vals)); %smallest eigenvalue = normal to plane   
    MLPlane_vec=pc_vecs(:,n).'; % vector defining norm of plane of symmetry (midline plane)
    %     NOTE: need to insert automatic left-right check, and reflect MLPlane vec if pointing to left (because PCA can have "arbitrary" inversions.
    VD_vec=cross(MLPlane_vec,acpc_vec); % ventral-dorsal vector

%%  set up rotation tracker
    dpts=[[MLPlane_vec];...
        [acpc_vec];...
        [VD_vec];...
        [0 0 0]];

%% first need to rotate y-dim into ac-pc vec, which means rotate around Zaxis and Xaxis

    % rotate about z -- (spherical-polar thought) - theta = atan(x/y)
    theta0=atan(acpc_vec(1)/acpc_vec(2));        
    Rz=[...
        [cos(theta0) -sin(theta0) 0];...
        [sin(theta0) cos(theta0) 0];...
        [0 0 1]];    
    a=[Rz*data.brain_info.mat(1:3,1:3)*dpts.'].'; %rotation tracker
    
    % rotate about x -- calculate rotation angle incorporate contraction in acpc projection onto x after rotation about z (factor of 1/theta0)
    phi1=atan(acpc_vec(3)/(acpc_vec(2)/cos(theta0)));  
    Rx=[...
        [1 0 0];...
        [0 cos(phi1) sin(phi1)];...
        [0 -sin(phi1)  cos(phi1)]];    
    b=[Rx*a.'].'; %rotation tracker   

%% then need to rotate x-dim into vector defining norm of plane of symmetry, by rotating about Yaxis
    
    % rotate about y -- y-axis is now coincident with rotated acpc, calculate angle for rotation about y
    eta2=-atan(b(3,1)/b(3,3)); % brute force approach to get angle
    
    Ry=[...
        [cos(eta2) 0 sin(eta2)];...
        [0 1 0];...
        [-sin(eta2) 0  cos(eta2)]];    
       
%% housekeeping
    clear a b n dpts tmp acpc_mcp acpc_vec MLP* VD* pc_* dt_projection 
    clear theta0 phi1 eta2

%% perform rotations and clipping
    rot_mat_acpc=[[Ry*Rx*Rz*data.brain_info.mat(1:3,1:3) Ry*Rx*Rz*data.brain_info.mat(1:3,1:3)*data.brain_info.mat(1:3,4)];[0 0 0 1]]; % matrix for rotation
    [data_acpc] = mri_dbs_rot_clip(data, rot_mat_acpc);
    data_acpc.brain_info.mat(1:3,4)=-data_acpc.MCP.'; % recenter brain around middle commisural point
    
%% package and print acpc data
    data_acpc.brain_info.dim=size(data_acpc.brain_vol);
    data_acpc.brain_info.mat(isnan(data_acpc.brain_info.mat))=1;
    data_acpc.brain_info.dt=data.brain_info.dt;
    tmp=data.brain_info.fname;tmp([-3:0]+end)=[]; data_acpc.brain_info.fname=[tmp '_acpc.nii'];  %delete .nii from end, and add _acpc.nii
    spm_write_vol(data_acpc.brain_info,5000*data_acpc.brain_vol); % need to scale into standard ranges
    
    %% ride along volumes
    for k=1:length(data.ridealong)
        tmpfname=data_acpc.ridealong(k).name; tmpfname([-3:0]+end)=[]; tmpfname=[tmpfname '_acpc.nii'];  %delete .nii from end, and add _acpc.nii
        tmp=data_acpc.brain_info; tmp.fname=tmpfname;
        spm_write_vol(tmp,data_acpc.ridealong(k).vol); % need to scale into standard ranges
    end
        
%% Now Shank rotations and reslicing -- from ac-pc coordinates
if mean(data_acpc.dbs(:,1)-data_acpc.MCP(1))>0
    Hside='R';
elseif mean(data_acpc.dbs(:,1)-data_acpc.MCP(1))<0
    Hside='L';
else
    error('aaa','dbs not lateralized w.r.t. ac-pc')
end

    % dbs shank axis
    [pc_vecs,dt_projection,pc_vals]=pca(data_acpc.dbs);    
    dbsAxis=pc_vecs(:,find(max(abs(pc_vals)))); % Note: This is in AC-PC space

    % then rotate about X - angle contracted by cos(hz_ang0)
    hx_ang1=atan(dbsAxis(3)/(dbsAxis(2)));  
    RxH=[...
        [1  0            0];...
        [0  cos(hx_ang1) sin(hx_ang1)];...
        [0 -sin(hx_ang1) cos(hx_ang1)]];     
    
    
%     % first rotate about Z
%     hz_ang0=atan(dbsAxis(1)/dbsAxis(2));    
%     RzH=[...
%         [cos(hz_ang0) -sin(hz_ang0) 0];...
%         [sin(hz_ang0)  cos(hz_ang0) 0];...
%         [0             0            1]];
% 
%     % then rotate about "new X" - angle contracted by cos(hz_ang0)
%     hx_ang1=atan(dbsAxis(3)/(dbsAxis(2)/cos(hz_ang0)));  
%     RxH=[...
%         [1  0            0];...
%         [0  cos(hx_ang1) sin(hx_ang1)];...
%         [0 -sin(hx_ang1) cos(hx_ang1)]]; 
    
%% Housekeeping
    clear ans dt_projection pc_*
    
%% rotate image and ride-alongs (need to have one function call for all of this and clipping in actual program)
%     rot_mat_dbs=[[RxH*RzH*Ry*Rx*Rz*data.brain_info.mat(1:3,1:3) data.brain_info.mat(1:3,4)];[0 0 0 1]]; % matrix for rotation
    rot_mat_dbs=[[RxH*Ry*Rx*Rz*data.brain_info.mat(1:3,1:3) data.brain_info.mat(1:3,4)];[0 0 0 1]]; % matrix for rotation
    [data_dbs] = mri_dbs_rot_clip(data, rot_mat_dbs);
    data_dbs.brain_info.mat(1:3,4)=-data_dbs.dbs(1,:); % recenter brain in dbstactic coordinates

%% package and print dbscampal data
    data_dbs.brain_info.dt=data.brain_info.dt;
    data_dbs.brain_info.dim=size(data_dbs.brain_vol);
    tmp=data.brain_info.fname;tmp([-3:0]+end)=[]; data_dbs.brain_info.fname=[tmp '_dbs_' Hside '.nii'];  %delete .nii from end, and add _dbs.nii
    spm_write_vol(data_dbs.brain_info,5000*data_dbs.brain_vol); % need to scale into standard ranges
    
    %% ride along volumes
    for k=1:length(data.ridealong)
        tmpfname=data_dbs.ridealong(k).name; tmpfname([-3:0]+end)=[]; 
        tmpfname=[tmpfname '_dbs_' Hside '.nii'];  %delete .nii from end, and add _dbs.nii
        tmp=data_dbs.brain_info; tmp.fname=tmpfname;
        spm_write_vol(tmp,data_dbs.ridealong(k).vol); % need to scale into standard ranges
    end  
%%

save([data.brain_info.fname(1:(end-4)) '_rotated_' Hside],'data_acpc','data_dbs','dbsAxis')

