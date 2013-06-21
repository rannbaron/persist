javaclasspath( strcat(pwd,'/comptopo/jars/comptopo.jar') );

global curdir comptopodir;
curdir=pwd;
comptopodir=strcat(curdir,'/comptopo');

addpath(strcat(pwd,'/comptopo'),strcat(pwd,'/generators'),strcat(pwd,'/means'),strcat(pwd,'/plotters'));

% The following (horrible hack) will ensure that the m_hungarian mex file
% points to the right libhungarian library on Darwin (OSX).
if strcmpi(mexext,'mexmaci64')
    [status,out]=system('otool -L means/m_hungarian.mexmaci64 | grep libhungarian.dylib | cut -d " " -f 1');
    if strcmpi(strtrim(out),'libhungarian.dylib')
        system('install_name_tool -change "libhungarian.dylib" "@loader_path/libhungarian.dylib" means/m_hungarian.mexmaci64');
    end
end

setenv('LD_LIBRARY_PATH',[getenv('LD_LIBRARY_PATH'),':',curdir,'/means'])