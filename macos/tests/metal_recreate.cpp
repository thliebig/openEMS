// Regression: an indexed Metal engine must restore operator UPML coefficient
// order before a fresh engine is constructed from the same operator.
#include "openems.h"
#include "FDTD/operator.h"
#include "FDTD/engine.h"
#include <cstring>
#include <iostream>
#include <stdexcept>
#include <vector>

class Fixture : public openEMS
{
public:
	std::vector<float> Fields(unsigned int steps)
	{
		FDTD_Eng->IterateTS(steps);
		std::vector<float> out;
		for (unsigned int x=0; x<FDTD_Op->GetNumberOfLines(0); ++x)
			for (unsigned int y=0; y<FDTD_Op->GetNumberOfLines(1); ++y)
				for (unsigned int z=0; z<FDTD_Op->GetNumberOfLines(2); ++z)
					for (unsigned int n=0; n<3; ++n) {
						out.push_back(FDTD_Eng->GetVolt(n,x,y,z));
						out.push_back(FDTD_Eng->GetCurr(n,x,y,z));
					}
		return out;
	}
	void Recreate()
	{
		delete FDTD_Eng; FDTD_Eng=nullptr;
		FDTD_Eng=FDTD_Op->CreateEngine();
	}
};

int main(int argc, char** argv)
{
	if (argc!=2) return 2;
	Fixture sim;
	sim.SetLibraryArguments({"engine=metal"});
	if (!sim.ReadFromXML(argv[1]) || sim.SetupFDTD()!=0) return 3;
	auto first=sim.Fields(400);
	sim.Recreate();
	auto got=sim.Fields(400);
	if (got.size()!=first.size() || std::memcmp(first.data(),got.data(),got.size()*sizeof(float)))
		throw std::runtime_error("Recreated engine fields differ");
	std::cout << "Recreated Metal engine: bit-identical fields" << std::endl;
}
