function [ interval0, interval1 ] = compTopoDM( dm, distanceBoundOnEdges )

% points is a Matlab matrix (supplied in the format of a CompTopo distance matrix)
% distanceBoundOnEdges is a real number, specifying the bound on the length of the edges

global curdir comptopodir;

cd (comptopodir);

% Specify to CompTopo that the supplied data is actually in the form of a distance matrix
% (we didn't have to do this for 'pointCloud', because it was the default)
supplyDataAs = 'distanceMatrix';

% Set up the parameters for calling the Java app from Matlab
param2 = strcat( 'distanceBoundOnEdges=', num2str( distanceBoundOnEdges ) );
param3 = strcat( 'supplyDataAs=', supplyDataAs );

echo off;

% Run the Java app
import edu.duke.math.comptopo.application.*;

ct=CompTopo();

ct.executeCompTopo( { 'settingsFile=comptopo.settings.txt', param2, param3 },points );

% Set the outputs
I0 = ct.getResultsM12( 0 ).getIntervals();
I1 = ct.getResultsM12( 1 ).getIntervals();

% Get rid of the third column (edge label numbers) for backward
% compatibility
interval0=I0(:,1:2);
interval1=I1(:,1:2);

cd (curdir);
end



