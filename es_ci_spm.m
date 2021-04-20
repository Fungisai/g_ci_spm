function [] = es_ci_spm(t_map,con,X,mask_img,confLevel,out_name)
%es_ci_spm.m Estimates effect size g and its CI from SPM t maps
%
%   Use by es_ci_spm(t_map,con,X,mask_img,confLevel,out_name) with
%   t_map:      string with the name of a SPM results nifti file with t values (usually something like 'spmT_000X.nii')
%   con:        contrast vector used to estimate the t map
%   X:          filtered and pre-whitened SPM design matrix (This should be SPM.xX.xKXs.X in the respective SPM.mat file. Please note that SPM is mean centering covariates) 
%   mask_img:   string with the name of the SPM mask file for the analysis (usually 'mask.nii')
%   confLevel:  confidence level of the estimated condfidence interval (usually something like .90 or .95)
%   out_name:   string with a prefix to add to the name of the results nifti files 
%
%   The results are saved as three nifti files:[out_name '_g.nii']
%   containing the g values map, and [out_name '_g_ci_l.nii'] and
%   [out_name '_g_ci_u.nii'] containing the lower and upper CI limit maps
%   respectively.
%
% Please note that the calculations are computationally expensive and can take up to several
% hours.
%
% This function is dependent on SPM 12 (https://www.fil.ion.ucl.ac.uk/spm/software/spm12/)
% and the Measures of Effect Size (MES) toolbox (https://github.com/hhentschke/measures-of-effect-size-toolbox).
%
% This work is is distributed under the terms of 
% the GNU General Public Licence as published by the Free Software Foundation
% (either version 3, or at your option, any later version).
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
% 
% If this file is of help for your own work, please cite
%
% Copyright (C) 2021
%
% Martin Fungisai Gerchen & Gordon Feld
% Department of Clinical Psychology
% Medical Faculty Mannheim/Heidelberg University 
% Central Institute of Mental Health
% Mannheim, Germany


%% preparations

% make contrast a column vector, if not already
if size(con,2)>size(con,1)
    con=con';
end

% pad contrast vector with zeros if not of sufficient length
if size(con,1) < size(X,2)
    con = [con; zeros(size(X,2)-size(con,1),1)];
end

% estimate degrees of freedom
DoF = size(X,1)-rank(X);

% estimate conversion factor from t to d
td_fac = sqrt(con'*pinv(X'*X)*con);

% flag for 'simple' two sample t-test without covariates
tst_flag = 0;
if size(X,2) == 2 && sum(con)==0 && sum(abs(con))==2 && (sum(X(:,1))+sum(X(:,2)))==size(X,1) && X(:,1)'*X(:,2)==0
    tst_flag = 1;
end

%% estimate effect size g
% load t map
V=spm_vol(t_map);
Y=spm_read_vols(V);

% keep t values for CI estimation 
t=Y;

% d
Y = Y.*td_fac;

% bias correction
J=(1-(3./(4*DoF-1))); % bias correction factor
Y = Y.*J;

%% estimate confidence interval for g 

if tst_flag % 'simple' two sample t-test
    
    n1 = sum(X(:,1));
    n2 = sum(X(:,2));
    
    alpha=1-confLevel;
    % critical z value corresponding to alpha
    zCrit=norminv(1-alpha/2);

    % approximate analytical CI
    se=sqrt((n1+n2)./(n1.*n2) + (Y.^2./(2*n1+2*n2-4))); 
    ci_l=Y-zCrit.*se; % CI lower limit
    ci_u=Y+zCrit.*se; % CI upper limit

else % all other tests
    
    % load mask file to extract voxel values
    V_m=spm_vol(mask_img);
    Y_m=spm_read_vols(V_m);

    ts = zeros(length(find(Y_m)),1); % vector of ts
    ts(:,1) = t(logical(Y_m));

    % estimate exact CI

        ci = zeros(2,numel(ts));

        % loop over t values and estimate CI for g
        for i=1:numel(ts)
            t_tmp = ts(i);
            ci_tmp=ncpci(t_tmp,'t',DoF,'confLevel',confLevel)'*td_fac; % uses the 'ncpci.m' function from the MES toolbox
            ci(:,i) = ci_tmp;
        end

    % put lower and upper CI limit into 3D maps
    ci_l = nan(size(Y_m));
    ci_l(logical(Y_m)) = ci(1,:);

    ci_u = nan(size(Y_m));
    ci_u(logical(Y_m)) = ci(2,:);
    
end

%% save images
V_out=V;

% g map
V_out.fname=[out_name '_g.nii'];
spm_write_vol(V_out,Y);

% g CI lower limit
V_out.fname=[out_name '_g_ci_l.nii'];
spm_write_vol(V_out,ci_l);

% g CI upper limit
V_out.fname=[out_name '_g_ci_u.nii'];
spm_write_vol(V_out,ci_u);

end
