#include <octave/oct.h>
#include <octave/ov-struct.h>
#include "hdf5.h"

// this special treatment is necessary because Win32-Octave ships with a very old hdf5 version (1.6.10)
void CloseH5Object(hid_t obj)
{
#if ((H5_VERS_MAJOR == 1) && (H5_VERS_MINOR == 6))
	// try group close, than Dataset close
	if (H5Gclose(obj)<0)
		H5Dclose(obj);
#else
	H5Oclose(obj);
#endif
}

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
	H5Eset_auto1(NULL, NULL);

	hid_t file = H5Fopen( args(0).string_value().c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );
	if (file==-1)
	{
	  error("h5readatt_octave: opening the given File failed");
	  return retval;
	}

#if ((H5_VERS_MAJOR == 1) && (H5_VERS_MINOR == 6))
	// this special treatment is necessary because Win32-Octave ships with a very old hdf5 version (1.6.10)
	hid_t obj = -1;
	//try opening the group
	obj = H5Gopen(file, args(1).string_value().c_str());
	//try opening the dataset if group failed
	if (obj==-1)
		obj = H5Dopen(file, args(1).string_value().c_str());
#else
	hid_t obj = H5Oopen(file, args(1).string_value().c_str(), H5P_DEFAULT);
#endif

	if (obj==-1)
	{
	  CloseH5Object(obj);
	  H5Fclose(file);
	  error("h5readatt_octave: opening the given Object failed");
	  return retval;
	}

	hid_t attr = H5Aopen_name(obj, args(2).string_value().c_str());
	if (attr==-1)
	{
	  CloseH5Object(obj);
	  H5Fclose(file);
	  error("h5readatt_octave: opening the given Attribute failed");
	  return retval;
	}

	hid_t type = H5Aget_type(attr);
	if (type<0)
	{
	  H5Aclose(attr);
	  CloseH5Object(obj);
	  H5Fclose(file);
	  error("h5readatt_octave: dataset type error");
	  return retval;
	}

	if (H5Tget_class(type)!=H5T_FLOAT)
	{
	  H5Aclose(attr);
	  CloseH5Object(obj);
	  H5Fclose(file);
	  error("h5readatt_octave: attribute type not supported");
	  return retval;
	}

	size_t numVal = H5Aget_storage_size(attr)/H5Tget_size(type);
	double value[numVal];
	if (H5Tget_size(type)==sizeof(float))
	{
	  float f_value[numVal];
	  if (H5Aread(attr, H5T_NATIVE_FLOAT, f_value)<0)
	  {
	    H5Aclose(attr);
	    CloseH5Object(obj);
	    H5Fclose(file);
	    error("h5readatt_octave: reading the given Attribute failed");
	    return retval;
	  }
	  for (size_t n=0;n<numVal;++n)
	    value[n] = f_value[n];
	}
	else if (H5Tget_size(type)==sizeof(double))
	{
	  if (H5Aread(attr, H5T_NATIVE_DOUBLE, value)<0)
	  {
	    H5Aclose(attr);
	    CloseH5Object(obj);
	    H5Fclose(file);
	    error("h5readatt_octave: reading the given Attribute failed");
	    return retval;
	  }
	}
	else
	{
	    H5Aclose(attr);
	    CloseH5Object(obj);
	    H5Fclose(file);
	    error("h5readatt_octave: reading the given Attribute failed: unknown type");
	    return retval;
	}

	H5Aclose(attr);
	CloseH5Object(obj);
	H5Fclose(file);
	Matrix mat(numVal,1);
	for (size_t n=0;n<numVal;++n)
		mat(n)=value[n];
	retval = octave_value(mat);
	return retval;
}

