%%xwavefront.m is a piece of code designed for finding the wavefront of
%%travelling wave solution. Often there are more than one value that are
%%numerically similar to the wavefront. This code identifies the
%%wavefrontand finds the midpoint of the wave to take a measurement of the
%%speed. This assumes that the population peaks at 1.
function X=xwavefront(a)
%finds the maximum value in the solution - this is just prior to the
%wavefront. The output is the value found maxNum, idx is the index of the
%value.
[~, idx]=max((a{1}(:)));
%u is a vector from the max solution value to the end of the original
%solution.
u=a{1}(idx:end);
%This finds the midpoint of the wave by finding the value closest to 0.5.
%val is the value, idx is the matrix index.
[~,idx]=min(abs(0.5-u));
%This finds the value in the original solution.
[row, col]=find(a{1}(:)==u(idx));
%This changes the linear indexing of the value to the more familiar
%(row,column) indexing which is what we need.
[row, col]=ind2sub(size(a{1}),row);
%The output is then the row where this value is which corresponds to
%spatial value pinpointing the wavefront.
X=max(row);