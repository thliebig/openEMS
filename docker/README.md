# openEMS / RF Miner Docker

## Image List

Output from `make list`:

```
(.venv) oberstet@intel-nuci7:~/scm/typedefint/rfminer/docker$ make list
docker images typedefint/openems
REPOSITORY           TAG                  IMAGE ID       CREATED          SIZE
typedefint/openems   master-from-source   fd60eb55a856   11 minutes ago   4.8GB
typedefint/openems   master               74f3431ca420   6 hours ago      3.55GB
typedefint/openems   v0.0.36              cfcfc6de46a7   6 hours ago      3.47GB
typedefint/openems   v0.0.35-110          637932cd8062   6 hours ago      3.47GB
```

## Versions ("from-source")

Output from `make version_openems_from_source_master`:

```
(.venv) oberstet@intel-nuci7:~/scm/typedefint/rfminer/docker$ make version_openems_from_source_master
docker run -it --rm \
    --entrypoint /openems/.local/bin/openEMS \
	typedefint/openems:master-from-source \

 ----------------------------------------------------------------------
 | openEMS 64bit -- version v0.0.36-12-ge5db9de
 | (C) 2010-2023 Thorsten Liebig <thorsten.liebig@gmx.de>  GPL license
 ----------------------------------------------------------------------
	Used external libraries:
		CSXCAD -- Version: v0.6.3-2-gc6a1587
		hdf5   -- Version: 1.10.7
		          compiled against: HDF5 library version: 1.10.7
		tinyxml -- compiled against: 2.6.2
		fparser
		boost  -- compiled against: 1_74
		vtk -- Version: 9.3.0
		       compiled against: 9.3.0

 Usage: openEMS <FDTD_XML_FILE> [<options>...]

 <options>
	--disable-dumps		Disable all field dumps for faster simulation
	--debug-material	Dump material distribution to a vtk file for debugging
	--debug-PEC		Dump metal distribution to a vtk file for debugging
	--debug-operator	Dump operator to vtk file for debugging
	--debug-boxes		Dump e.g. probe boxes to vtk file for debugging
	--debug-CSX		Write CSX geometry file to debugCSX.xml
	--engine=<type>		Choose engine type
		--engine=fastest		fastest available engine (default)
		--engine=basic			basic FDTD engine
		--engine=sse			engine using sse vector extensions
		--engine=sse-compressed		engine using compressed operator + sse vector extensions
		--engine=multithreaded		engine using compressed operator + sse vector extensions + multithreading
	--numThreads=<n>	Force use n threads for multithreaded engine (needs: --engine=multithreaded)
	--no-simulation		only run preprocessing; do not simulate
	--dump-statistics	dump simulation statistics to 'openEMS_run_stats.txt' and 'openEMS_stats.txt'

	 Additional global arguments
	--showProbeDiscretization	Show probe discretization information
	--nativeFieldDumps		Dump all fields using the native field components
	-v,-vv,-vvv			Set debug level: 1 to 3

failed to resize tty, using default size
make: *** [Makefile:47: version_openems_from_source_master] Fehler 255
(.venv) oberstet@intel-nuci7:~/scm/typedefint/rfminer/docker$
```

## Versions ("from-upstream")

Output from `make version_openems_from_upstream_all`:

```
(.venv) oberstet@intel-nuci7:~/scm/typedefint/rfminer/docker$ make version_openems_from_upstream_all
for version in v0.0.35-110 v0.0.36 master; do \
	docker run -it --rm \
		--entrypoint /openems/.local/bin/openEMS \
		typedefint/openems:$version ; \
done
 ----------------------------------------------------------------------
 | openEMS 64bit -- version v0.0.35-110-g782a738
 | (C) 2010-2023 Thorsten Liebig <thorsten.liebig@gmx.de>  GPL license
 ----------------------------------------------------------------------
	Used external libraries:
		CSXCAD -- Version: v0.6.2-125-ge571949
		hdf5   -- Version: 1.10.7
		          compiled against: HDF5 library version: 1.10.7
		tinyxml -- compiled against: 2.6.2
		fparser
		boost  -- compiled against: 1_74
		vtk -- Version: 9.1.0
		       compiled against: 9.1.0

 Usage: openEMS <FDTD_XML_FILE> [<options>...]

 <options>
	--disable-dumps		Disable all field dumps for faster simulation
	--debug-material	Dump material distribution to a vtk file for debugging
	--debug-PEC		Dump metal distribution to a vtk file for debugging
	--debug-operator	Dump operator to vtk file for debugging
	--debug-boxes		Dump e.g. probe boxes to vtk file for debugging
	--debug-CSX		Write CSX geometry file to debugCSX.xml
	--engine=<type>		Choose engine type
		--engine=fastest		fastest available engine (default)
		--engine=basic			basic FDTD engine
		--engine=sse			engine using sse vector extensions
		--engine=sse-compressed		engine using compressed operator + sse vector extensions
		--engine=multithreaded		engine using compressed operator + sse vector extensions + multithreading
	--numThreads=<n>	Force use n threads for multithreaded engine (needs: --engine=multithreaded)
	--no-simulation		only run preprocessing; do not simulate
	--dump-statistics	dump simulation statistics to 'openEMS_run_stats.txt' and 'openEMS_stats.txt'

	 Additional global arguments
	--showProbeDiscretization	Show probe discretization information
	--nativeFieldDumps		Dump all fields using the native field components
	-v,-vv,-vvv			Set debug level: 1 to 3

 ----------------------------------------------------------------------
 | openEMS 64bit -- version v0.0.36
 | (C) 2010-2023 Thorsten Liebig <thorsten.liebig@gmx.de>  GPL license
 ----------------------------------------------------------------------
	Used external libraries:
		CSXCAD -- Version: v0.6.3
		hdf5   -- Version: 1.10.7
		          compiled against: HDF5 library version: 1.10.7
		tinyxml -- compiled against: 2.6.2
		fparser
		boost  -- compiled against: 1_74
		vtk -- Version: 9.1.0
		       compiled against: 9.1.0

 Usage: openEMS <FDTD_XML_FILE> [<options>...]

 <options>
	--disable-dumps		Disable all field dumps for faster simulation
	--debug-material	Dump material distribution to a vtk file for debugging
	--debug-PEC		Dump metal distribution to a vtk file for debugging
	--debug-operator	Dump operator to vtk file for debugging
	--debug-boxes		Dump e.g. probe boxes to vtk file for debugging
	--debug-CSX		Write CSX geometry file to debugCSX.xml
	--engine=<type>		Choose engine type
		--engine=fastest		fastest available engine (default)
		--engine=basic			basic FDTD engine
		--engine=sse			engine using sse vector extensions
		--engine=sse-compressed		engine using compressed operator + sse vector extensions
		--engine=multithreaded		engine using compressed operator + sse vector extensions + multithreading
	--numThreads=<n>	Force use n threads for multithreaded engine (needs: --engine=multithreaded)
	--no-simulation		only run preprocessing; do not simulate
	--dump-statistics	dump simulation statistics to 'openEMS_run_stats.txt' and 'openEMS_stats.txt'

	 Additional global arguments
	--showProbeDiscretization	Show probe discretization information
	--nativeFieldDumps		Dump all fields using the native field components
	-v,-vv,-vvv			Set debug level: 1 to 3

 ----------------------------------------------------------------------
 | openEMS 64bit -- version v0.0.36-12-ge5db9de
 | (C) 2010-2023 Thorsten Liebig <thorsten.liebig@gmx.de>  GPL license
 ----------------------------------------------------------------------
	Used external libraries:
		CSXCAD -- Version: v0.6.3-2-gc6a1587
		hdf5   -- Version: 1.10.7
		          compiled against: HDF5 library version: 1.10.7
		tinyxml -- compiled against: 2.6.2
		fparser
		boost  -- compiled against: 1_74
		vtk -- Version: 9.1.0
		       compiled against: 9.1.0

 Usage: openEMS <FDTD_XML_FILE> [<options>...]

 <options>
	--disable-dumps		Disable all field dumps for faster simulation
	--debug-material	Dump material distribution to a vtk file for debugging
	--debug-PEC		Dump metal distribution to a vtk file for debugging
	--debug-operator	Dump operator to vtk file for debugging
	--debug-boxes		Dump e.g. probe boxes to vtk file for debugging
	--debug-CSX		Write CSX geometry file to debugCSX.xml
	--engine=<type>		Choose engine type
		--engine=fastest		fastest available engine (default)
		--engine=basic			basic FDTD engine
		--engine=sse			engine using sse vector extensions
		--engine=sse-compressed		engine using compressed operator + sse vector extensions
		--engine=multithreaded		engine using compressed operator + sse vector extensions + multithreading
	--numThreads=<n>	Force use n threads for multithreaded engine (needs: --engine=multithreaded)
	--no-simulation		only run preprocessing; do not simulate
	--dump-statistics	dump simulation statistics to 'openEMS_run_stats.txt' and 'openEMS_stats.txt'

	 Additional global arguments
	--showProbeDiscretization	Show probe discretization information
	--nativeFieldDumps		Dump all fields using the native field components
	-v,-vv,-vvv			Set debug level: 1 to 3

make: *** [Makefile:87: version_openems_from_upstream_all] Fehler 255
(.venv) oberstet@intel-nuci7:~/scm/typedefint/rfminer/docker$
```
