function [output info] = RawToPersistenceCell(raw_data, sample_size, dim_start, dim_end, num_samples, signal_image)
% This function takes a raw data set, along with parameters describing
% sampling behavior, and returns a cell array containing persistence
% diagrams for the desired sampled point clouds. At the outermost level is
% a cell for each dimension within the specified range; each "dimension
% cell" contains further cells corresponding to each sample taken at that
% dimension; each of these "sample cells" contains two n-by-2 arrays
% corresponding to 0- and 1-dimensional persistence.
%   @param data_raw: vector or matrix representing the image or signal
%   @param sample_size: samples taken from the point cloud will be of size
%           sample_size. If the sample_size given is -1, or if it's larger
%           than the size of the data itself, the full size of data_raw
%           will automatically be used.
%   @param dim_start: starting dimension of desired range of dimensions of ambient
%           point coud (inclusive)
%   @param dim_end: ending dimension of desired range of dimensions of ambient
%           point cloud (inclusive)
%   @param num_samples: for each desired dimension of ambient point cloud and for
%           both 0- and 1-dimensional homology groups, the number samples
%           taken (and, by consequence, the number of persistence diagrams
%           returned) will be equal to num_samples (this will automatically
%           be changed to 1 if there is no downsampling).
%   @param signal_image: if this is equal to the string 'image', raw_data will be
%           treated as representing an image and the size of the masking
%           matrix will be dim_start by dim_end. If this is equal to the
%           string 'signal', raw_data will be interpreted as a signal and
%           only its right-most column will be used to create the point
%           cloud.
%   return- output: 1xn cell array, where n=dim_end-dim_start+1. Each cell
%           is a 2xm cell array where m=num_samples. The top row consists
%           of 0-dimensional homology persistence diagrams, and the bottom
%           row are the 1-dimensional diagrams.
global curdir comptopodir;

cd (comptopodir);

%We must check whether it is a signal or an image so we can do the sliding
%window technique or masking
if strcmp(signal_image, 'image') 
    output = cell(1, 1); %output of function
    c = cell(2, num_samples); %cell array contained in the first and only cell of output
    point_cloud = (im2col(raw_data,[dim_start dim_end]))';%point cloud resulting from masking
    point_cloud = normalize(point_cloud); %point cloud projected onto the unit ball
    [rows cols] = size(point_cloud); 
    for i=1:num_samples
        %must check for invalid or out-of-range sample size
        if sample_size==-1||sample_size>rows
            sample_size=rows;
        end
        perm = randperm(rows); %perm is now a vector containing the numbers from 1 to rows randomly permuted
        perm = perm(1:sample_size); %we truncate perm at the desired sample size
        sampled = zeros(sample_size, cols); %sampled will be our new, down-sampled point cloud
        for j=1:sample_size
            sampled(j,:)=point_cloud(perm(j),:); %randomly sample points from the point cloud by using the permutation vector
        end
        [p q] = compTopoPC(sampled, max(pdist(sampled))); %p contains the 0-dimensional diagram, q the 1-dimensional diagram
        c{1,i}=p;
        c{2,i}=q;
    end
    output{1}=c;
%we go through the same process with a signal, but we do sliding window
%instead of masking
elseif strcmp(signal_image, 'signal')
    output = cell(1, dim_end-dim_start+1);%output now can have many columns if the user wants multiple sliding window dimensions
    for z=dim_start:dim_end %for each sliding window dimension, we go through the same process as before
        c = cell(2, num_samples);
        [rows cols] = size(raw_data);
        point_cloud = (im2col(raw_data(:,cols),[z, 1]))';
        point_cloud = normalize(point_cloud);
        for i=1:num_samples
            if sample_size==-1||sample_size>size(point_cloud,1)
                sample_size=size(point_cloud,1);
            end
            perm = randperm(rows-(z-1)); %rows -(z-1) is equivalent to the number of rows in the point_cloud because our variable rows corresponds to the original signal
            perm = perm(1:sample_size);
            sampled = zeros(sample_size, z);
            for j=1:sample_size
                sampled(j,:)=point_cloud(perm(j),:);
            end
            [p q] = compTopoPC(sampled, max(pdist(sampled)));
            c{1,i}=p;
            c{2,i}=q;
        end
        output{1,z-dim_start+1}=c;
    end
end
info = ['Input Name: ???'  ', Sample Size: ' int2str(sample_size) ', Dim Start: ' int2str(dim_start) ', Dim End: ' int2str(dim_end)];
cd (curdir);
end

%projection of the point cloud onto the unit ball
function [ normalized ] = normalize( point_cloud )
biggest_norm = 0;
for i=1:size(point_cloud,1)
    if norm(point_cloud(i,:))>biggest_norm
        biggest_norm=norm(point_cloud(i,:));
    end
end
for i=1:size(point_cloud,1)
    point_cloud(i,:)=point_cloud(i,:)/biggest_norm;
end
normalized = point_cloud;
end