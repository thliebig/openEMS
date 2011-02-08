#!/bin/bash

#clear LD_LIBRARY_PATH
export LD_LIBRARY_PATH=

#get path to openEMS
openEMS_PATH=`dirname $0`

#execute openEMS
exec mpirun -l -n 2 $openEMS_PATH/openEMS $@

