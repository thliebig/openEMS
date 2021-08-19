#!/bin/sh

#clear LD_LIBRARY_PATH
export LD_LIBRARY_PATH=

#get path to openEMS
openEMS_PATH=`dirname $0`

#execute openEMS
exec $openEMS_PATH/openEMS "$@"

