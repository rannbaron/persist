function [ pcloud ] = spheresample( dim,samplesize,noisemultiplier )
%spheresample Generates point clouds from spheres
%   This function will return a size samplesize point cloud, sampled
%   uniformly from a sphere of dimension dim. Radius = 1.

    pcloud = spheresampler(dim,1,samplesize,noisemultiplier);
    
end

