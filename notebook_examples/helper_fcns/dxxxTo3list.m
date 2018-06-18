function chemlist = dxxxTo3list(dxxx,dict)
%function that takes a string such as 'dfag' and converts it to a cell list
%of 3-letter chemistries, such as {'PAS','PHE','ALA','GLY'}
%strips off anything after '-' for core chemistries like 'dfag-ndi'
dxxx = strsplit(dxxx,'-');
chemlist = {};
if length(dxxx) > 1
   dxxx = dxxx(1); 
end
dxxx = char(dxxx);
for l = 1:length(dxxx)
   letter = dxxx(l);
   chemlist{l} = dict(letter);
end