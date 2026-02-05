function convFTmsh(name,f,v,bubbleNo, outerFluidDensity, ...
   density,  outerFluidViscosity, viscosity, surfaceTensionCoeff, Eotvos, ...
   Morton,diameter,convertToMeters)
% This function gets face and vertices data and write them in a format
% readable by FTFoam solver as an input mesh. 
%
% USAGE: convFTmsh(outputName,f,v,bubbleNo, outerFluidDensity, ...
%   density,  outerFluidViscosity, viscosity, surfaceTensionCoeff, Eotvos, ...
%   Morton,diameter,convertToMeters);
% 


    [path,file,fext] = fileparts(name) ;
   
 
    try
%-- try to write data to file
    
    ffid = fopen(name, 'w') ;
    
    nver = +3;
    
%-- write Header and Bubble Data 
    disp('Writing Header and Bubble Data...');

    fprintf(ffid,[ ...
    '  %u\n  %u %u\n  %g %g %g %g %g %g %g\n  %g\n'],bubbleNo,size(v,1),size(f,1),...
    outerFluidDensity,density,outerFluidViscosity,viscosity,surfaceTensionCoeff,Eotvos,Morton,diameter) ;

    save_FT_mesh(ffid,nver,f,v,convertToMeters) ; 
    
    fclose(ffid);
    disp('Conversion was Successful.');

    catch err
    
%-- ensure that we close the file
    disp('Somthing went Wrong!');

    if (ffid>-1)
    fclose(ffid) ;
    end
    rethrow(err) ;
        
    end
    
end
        
function save_FT_mesh(ffid,nver,f,v,convertToMeters)
  
    npts = +0;
%-- write "POINT" data 
    disp('Writing Point Data...');

    ndim = 3 ;
    npts = size(v,1) - 0 ;    
            
    if (isa(v,'double')) 
        vstr = sprintf('%%1.%uf ',16);
    else
        vstr = sprintf('%%1.%uf ',8);
    end
    temp = v;
    temp = temp * convertToMeters;
    fprintf(ffid, ...
    [repmat(vstr,1,ndim),'\n'],temp');
     
%-- write Face Index data
    disp('Writing Face Data...');

    index = f;
    index(:,1:3) = index(:,1:3)-1 ; % file is zero-indexed           
    fprintf( ...
       ffid,[repmat('%u ',1,3),'\n'],index'); 
    
end



