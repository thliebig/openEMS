#include <octave/oct.h>
#include <octave/ov-struct.h>
#include "hdf5.h"

DEFUN_DLD (h5readatt_octave, args, nargout, "h5readatt_octave(<File_Name>,<DataSet_Name>,<Attribute_Name>)")
{
	octave_value retval;
	int nargin = args.length();
	if (nargin != 3)
	{
		print_usage();
		return retval;
	}
	if ((args(0).is_string()==false) || (args(1).is_string()==false) || (args(2).is_string()==false))
	{
		print_usage();
		return retval;
	}

	//suppress hdf5 error output
	H5Eset_auto(NULL, NULL);

	hid_t file = H5Fopen( args(0).string_value().c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );
	if (file==-1)
	{
	  error("h5readatt_octave: opening the given File failed");
	  return retval;
	}

	hid_t ds = H5Dopen(file, args(1).string_value().c_str());
	if (ds==-1)
	{
	  error("h5readatt_octave: opening the given DataSet failed");
	  return retval;
	}

	hid_t attr = H5Aopen_name(ds, args(2).string_value().c_str());
	if (attr==-1)
	{
	  error("h5readatt_octave: opening the given Attribute failed");
	  return retval;
	}

	float value;
	if (H5Aread(attr, H5T_NATIVE_FLOAT, &value)<0)
	{
	  error("h5readatt_octave: reading the given Attribute failed");
	  return retval;
	}

	H5Aclose(attr);
	H5Dclose(ds);
	H5Fclose(file);
	retval = octave_value(value);
	return retval;
}

