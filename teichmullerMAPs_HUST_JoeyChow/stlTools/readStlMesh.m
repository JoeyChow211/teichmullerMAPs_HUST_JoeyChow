function [v, f, n, name] = readStlMesh(fileName)
% READSTLMESH Read an STL file in either ASCII or binary format.
% Author: Joey.Chow
% READSTLMESH reads any STL file regardless of its format
%V are the vertices
%F are the faces
%N are the normals
%NAME is the name of the STL object (NOT the name of the STL file)

format = detectStlFormat(fileName);
if strcmp(format,'ascii')
  [v,f,n,name] = readAsciiStlMesh(fileName);
elseif strcmp(format,'binary')
  [v,f,n,name] = readBinaryStlMesh(fileName);
end
