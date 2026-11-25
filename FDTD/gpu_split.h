/*
*	Copyright (C) 2026 Roel Budé
*
*	This program is free software: you can redistribute it and/or modify
*	it under the terms of the GNU General Public License as published by
*	the Free Software Foundation, either version 3 of the License, or
*	(at your option) any later version.
*
*	This program is distributed in the hope that it will be useful,
*	but WITHOUT ANY WARRANTY; without even the implied warranty of
*	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*	GNU General Public License for more details.
*
*	You should have received a copy of the GNU General Public License
*	along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef GPU_SPLIT_H
#define GPU_SPLIT_H

#include <string>
#include <vector>

#include "gpu_backend.h"
#include "gpu_ext.h"

class Engine_Extension;

//! A grid split along x over several backends of one kind, e.g. one per GPU (--gpu-devices)
/*!
  Each part is a backend of its own (GPU_Backend::New() creates them): it holds
  the x lines it owns and, at each cut, the line on the other side of it, a
  halo line (GPU_Slab), with the cell model of those lines. The main updates
  and the device extensions run on every part, the extensions on the points
  the part owns (GPU_CreatePartExtension()). Before a main update the parts
  copy the halo lines it reads from their neighbours (GPU_Backend::SendPlane()),
  and before a hook of an extension the halo lines it reads
  (GPU_Extension::HaloReads()): the results are those of one backend.

  Engine_GPU runs it like one backend, through GPU_CreateExtension() and the
  methods here; it has the fused step of the parts (PrepareFused()), but no
  temporal blocking, no snapshots and no device sums of the dumps (the dumps
  read the host mirror), and no cylindrical multi-grid. Any number of parts
  may share a device, for tests.
  */
class GPU_Backend_Split : public GPU_Backend
{
public:
	//! The grid split over \a parts, which it takes over: at least two, which can copy to each other (GPU_Backend::CanSendPlane())
	explicit GPU_Backend_Split(const std::vector<GPU_Backend*>& parts);
	virtual ~GPU_Backend_Split();

	virtual std::string GetName() const;

	//! Split the grid of \a model at the cuts (--gpu-domain-cuts, else equal parts) and initialize the parts with their slices of it
	virtual bool Init(const Operator* op, const GPU_CellModel& model);

	//! Not for a split grid: the device extensions run on the parts (CreateExtension())
	virtual void* DeviceAlloc(size_t bytes, const void* host=NULL);
	virtual void DeviceDownload(void* host, const void* dev, size_t bytes);
	virtual float* DeviceField(bool currents);
	virtual void RunPoints(GPU_PointOp kind, const void* op, uint64_t begin, uint64_t end, bool serial, const CM_PointFilter& filter);

	//! The main update of every part on its lines in \a box, after the copies of the halo lines it reads
	virtual void UpdateVoltages(const Engine_Box& box);
	virtual void UpdateCurrents(const Engine_Box& box);

	//! The fused step of every part, \a late in nodes of the whole grid; the currents of the last line a part owns are late too (except on the last part): they read the halo line above
	virtual bool PrepareFused(const std::vector<uint64_t>& late);
	virtual void ReleaseFused();
	//! The fused step of every part, after the copies of the halo lines of the currents below
	virtual void FusedUpdate();
	//! The late currents of every part, after the copies of the halo lines of the voltages above
	virtual void SeamUpdate();

	virtual void Download(bool currents, ArrayLib::ArrayNIJK<FDTD_FLOAT>& field);
	virtual void Synchronize();
	virtual bool CalcFastEnergy(const unsigned int numNodes[3], double& E_energy, double& H_energy);

	//! No cylindrical multi-grid on a split grid: throws
	virtual GPU_Backend* NewSubGridBackend();

	//! The device implementation of \a eng_ext on every part (GPU_CreatePartExtension()); NULL if there is none
	virtual GPU_Extension* CreateExtension(Engine_Extension* eng_ext, const unsigned int numLines[3]);
	virtual const char* NoTemporalBlocking() const {return "a grid split over several devices";}
	//! The z line from the part that owns it
	virtual void ReadLine(bool currents, unsigned int n, unsigned int x, unsigned int y, const unsigned int numLines[3], FDTD_FLOAT* line);

	size_t Parts() const {return m_Parts.size();}
	const GPU_Slab& Slab(size_t p) const {return m_Slabs.at(p);}
	//! Copy the halo lines of the voltages (\a currents false) or currents from the neighbours: 1 the line below each part, 2 the one above, 3 both
	void Exchange(bool currents, int halos);
	//! Exchange() at cut \a c only, between part \a c and part \a c+1: 1 the halo line of part \a c+1, 2 that of part \a c
	void ExchangeCut(size_t c, bool currents, int halos);
	//! The nodes of \a box (lines of the whole grid) part \a p owns, in its local lines; false if none
	bool LocalBox(size_t p, const Engine_Box& box, Engine_Box& local) const;

protected:
	std::vector<GPU_Backend*> m_Parts;
	std::vector<GPU_Slab> m_Slabs;
};

//! The x lines at which a grid of \a nx x lines is split into \a parts parts: \a parts+1 lines from 0 to \a nx
/*!
  From --gpu-domain-cuts, else the environment variable OPENEMS_GPU_DOMAIN_CUTS,
  else equal parts. Throws if it is not such a list, or a part would own fewer
  than 3 lines.
  */
std::vector<unsigned int> GPU_SplitCuts(unsigned int nx, size_t parts);

//! The comma separated numbers of the option \a option, its \a value, else of the environment variable \a env; none if neither is set
/*!
  \a what is the one used, as "name=list" for messages, and \a from_env
  whether it is the environment variable. Throws if it is not such a list.
  */
std::vector<unsigned int> GPU_NumberList(const std::string& value, const std::string& option, const char* env, std::string& what, bool& from_env);

#endif // GPU_SPLIT_H
