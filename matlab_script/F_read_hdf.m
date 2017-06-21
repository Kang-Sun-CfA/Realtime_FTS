function datavar = F_read_hdf(filename,varname)

file_id = H5F.open (filename);

for i = 1:length(varname)
    % Open the dataset
    data_id = H5D.open (file_id, varname{i});
    % Read attributes.
%     nattr  = H5A.get_num_attrs(dset_id);
%     for iattr = 0:nattr
        
%     try
%         ATTRIBUTE = 'Offset';
%         attr_id = H5A.open_name (data_id, ATTRIBUTE);
%         datavar.(varname{i}).(ATTRIBUTE) = H5A.read (attr_id, 'H5ML_DEFAULT');
%         
%         ATTRIBUTE = 'ScaleFactor';
%         attr_id = H5A.open_name (data_id, ATTRIBUTE);
%         datavar.(varname{i}).(ATTRIBUTE) = H5A.read (attr_id, 'H5ML_DEFAULT');
%     catch
%         warning('No attributes to read!')
%     end
    % Read the dataset.
    datavar.(varname{i}).data=H5D.read (data_id);
    %     datavar.(varname{i}).name = long_name(:)';
end

% Close and release resources.
% H5A.close (attr_id)
H5D.close (data_id);
H5F.close (file_id);