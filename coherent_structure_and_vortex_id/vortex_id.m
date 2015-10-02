function [ centers, circL, area, cmat, circA] = vortex_id( cal_mat, thresh, U, V, pix, filt, spacing, vorticity, method)
% This code works to:
% 1) Find the center of vortices dervied from the following coherent
% structure identification measures:
%   a) Q-Criterion
%   b) D-Criterion
%   c) Lambda 2
%   d) Lambda CI
% 2) Finds the areas in which the structure occupies
% 3) Calculates the circL of each identified structure using
%   a) Line Integral
%   b) Vortex Area Integral
%
% NOTE: Code has been adopted from previous vortex ID methods written by
% Dave Hubble & Chris Weiland
%
% v1: Brett Meyers, 11/21/2014

% Initializes new cal_mat & dimensions for calculation
cmat        = cal_mat;
[s1,s2,~]   = size(cal_mat);
Vmag_min    = zeros([s1 s2]); % Magnitude Velocity Surface
W_surface   = zeros([s1 s2]); % Vorticity Surface
C_surface   = zeros([s1 s2]); % cal_mat Surface

% This part of the function computes the contour set (CS) for which there
% is an identified coherent structure.  The values in CS(1,:) represent the
% perimeter points for each structure in the x-domain.  The values in
% CS(2,:) represent the perimeter points for each structure in the y-domain

figure(100)
[CS,~] = contourf(cal_mat,[thresh  thresh]);
CS(isnan(CS))=0;
% close(100)

% At this point CS is evaluated ("checked") to ensure that structures
% exist.  If no values are found or CS is small, the loop sets check to 1.
check = isempty(CS); 

if check==0;
    if size(CS,2)<5;
        check=1;
    end
end

% initialize a counter
counter=0 ;

if check==0
    [~  ,x1  ]  = find(CS(1,:)==thresh);  % find start of each region
    num_points  = CS(2,x1);               % number of points in each vortex
    
    % coherent structure check is present to ensure that a structure is
    % present in the cell.  if not, the matrix will be filled in with zeros
    % accordingly.
    % (this only exists for conditions where the filtered area removes all
    % structures from the field due to their small size)
    
    coherent_structure_check = 0; 
    
    % loop through each vortex region
    for j=1:1:max(size(num_points))
        % Check for exeeding matrix dimensions when extra points crop up;
        % set to zero
        xtemp = CS(1,(x1(j)+1):((x1(j)+1)+num_points(j))-1); % x points of vortex perimeter
        ytemp = CS(2,(x1(j)+1):((x1(j)+1)+num_points(j))-1); % y points of vortex perimeter
        
        % Calculate the Velocity Magnitude and an empty minimum matrix
        Vmag     = (U.^2+V.^2).^0.5;% Velocity Magnitude
        
        % Generate linear spaced mesh grid for domain
        [X,Y]=meshgrid(1:s2,1:s1);
        
        % Set the temporary holders to the valid holders
        xval=xtemp;
        yval=ytemp;
        
        % Remove any 0s from the temporary hold
        temp1=0;
        for ww=1:length(xval)
            if xval(ww)~=0 && yval(ww)~=0
                temp1=temp1+1;
                xtemp(temp1)=xval(ww);
                ytemp(temp1)=yval(ww);
            end
        end
        
        % Generate Temporary vectors (unsure of use)
        tempx=xtemp;
        tempy=ytemp;
        
        % Check to see if the area of the vortex is in the filter area
        if (polyarea(tempx.*pix*spacing,tempy.*pix*spacing))>=filt                    % area
            % A check to make sure current CS is a set continuous vortex 
            % points, not a string of independent variables
            if max(diff(sort(tempx))) <= 1.001 && max(diff(sort(tempy))) <= 1.001;
                % coherent structure check is triggered here for cases
                % where structures do, indeed, exist in field
                coherent_structure_check = 1;
                
                % two counteres are initialized here to count in the
                % following section
                %
                % NOTE: this could be modified if the number of elements
                % that exist in the region are determined.  This allows for
                % the initialization of list matrices (ie W_list =
                % zeros(N,1))
                
                counter=counter+1;
                k=1;
                                
                % determine if point is within roi
                cspace = inpolygon(X,Y,xtemp,ytemp);
                
                % if points within cspace exist, generate list vectors for
                % properties of: Magnitude Velocity, Vorticity, cal_mat
                % Additionally, generate Magnitdue Velocity Vorticity & 
                % cal_mat surfaces
                %
                % if points are not within the roi, set surface points to 
                % NaN
                
                for ii=1:s1
                    for jj=1:s2
                        if cspace(ii,jj)==1
                            Vmag_min(ii,jj)     = Vmag(ii,jj);      % Magnitude Velocity Surface
                            V_list(k)           = Vmag(ii,jj);      % Magnitude Velocity List
                            W_surface(ii,jj)    = vorticity(ii,jj); % Vorticity Surface
                            W_list(k)           = vorticity(ii,jj); % Vorticity List
                            C_surface(ii,jj)    = cal_mat(ii,jj);   % cal_mat Surface
                            C_list(k)           = cal_mat(ii,jj);   % cal_mat List
                            k=k+1;
                        else  
                            Vmag_min(ii,jj)     = nan;   % Magnitude Surface
                            W_surface(ii,jj)    = nan;   % Vorticity Surface
                            C_surface(ii,jj)    = nan;   % cal_mat Surface
                        end
                    end
                end
                
                % Perform vortex center identificiation, vortex area
                % calculation & circulation by line & area
                
                % This part of the if statement runs vortex identification
                % by Vorticity surface
                if strcmp(method,'weightW')   
                    indx     = find(~isnan(W_surface));  
                    center_x = W_surface(indx).*X(indx);
                    W_mass   = W_surface(indx);
                    center_y = W_surface(indx).*Y(indx);
                    
                    if exist('center_x') && exist('center_y')
                        X_loc=sum(center_x)/sum(W_mass);
                        Y_loc=sum(center_y)/sum(W_mass);
                    else
                        X_loc=0;
                        Y_loc=0;
                    end
                    
                    clear W_mass center_x center_y
                    
                    centers(counter,1)  = X_loc;
                    centers(counter,2)  = Y_loc;
                    area(counter)       = polyarea(tempx.*pix*spacing,tempy.*pix*spacing);            
                    
                    [L_circ,A_circ]     = vortex_circulation( U, V, tempx, tempy, pix, spacing, W_surface);
                    circL(counter)      = L_circ;                                                                    % Circulation calculated using velocity
                    circA(counter)      = A_circ;                                                            % Circulation calculated using vorticity
                
                    keyboard
                    
                % This part of the if statement runs vortex identification
                % by cal_mat surface    
                elseif strcmp(method,'weightC')
                    indx        = find(~isnan(C_surface));
                    center_x    = C_surface(indx).*X(indx);
                    W_mass      = C_surface(indx);
                    center_y    = C_surface(indx).*Y(indx);
                    
                    if exist('center_x') && exist('center_y');
                        X_loc=sum(center_x)/sum(W_mass);
                        Y_loc=sum(center_y)/sum(W_mass);
                    else
                        X_loc=0;
                        Y_loc=0;
                    end
                    
                    clear W_mass center_x center_y
                    
                    centers(counter,1)  = X_loc;
                    centers(counter,2)  = Y_loc;
                    area(counter)=polyarea(tempx.*pix*spacing,tempy.*pix*spacing);
                    
                    [L_circ,A_circ]     = vortex_circulation( U, V, tempx, tempy, pix, spacing, C_surface);
                    circL(counter)      = L_circ;                                                                   
                    circA(counter)      = A_circ;
                    
                    keyboard
                end
            end
        else
            % Find indices of all points in the discardable area 
            % (vorticies that did not meet area restrictions) and set
            % them to zero inside cal_mat
            cal_mat2    = 1-roipoly(cal_mat,tempx,tempy);
            cmat        = cmat.*cal_mat2;
        end
        clear tempvecx tempvecy xval yval xtemp ytemp W_surface W_list C_surface
    end
    
    % If structures have been fully filtered from the field, this exists to
    % make zero
    %
    % NOTE: in the future, if possible, this can be initialized, removing
    % the need for this check
    if coherent_structure_check == 0
        centers(1:10,1:2)=0;
        circL(1,1:10)=0;
        circA(1,1:10)=0;
        area(1,1:10)=0;
        dEdt_core(1,1:10)=0;
    end
    % If the field is completely empty (prior to filtering) the values are
    % all zero in this case
else
    centers(1:10,1:2)=0;
    circL(1,1:10)=0;
    area(1,1:10)=0;
    circA(1,1:10)=0;
    dEdt_core(1,1:10)=0;
end

% These are random initializers in the case that values are NaNs
centers(isnan(centers)) = 5;
circA(isnan(circA))     = 0;