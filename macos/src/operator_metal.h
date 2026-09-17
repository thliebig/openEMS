/*
* Copyright (C) 2026 openEMS contributors
*
* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*/

#ifndef OPERATOR_METAL_H
#define OPERATOR_METAL_H

#include "operator_sse.h"

#include <string>
#include <vector>

//! True when a Metal GPU device can be created; otherwise fills \a reason.
//! Checked before the (expensive) CPU operator setup so a missing GPU aborts
//! early instead of after material/PEC sampling.
bool MetalDeviceAvailable(std::string& reason);

class Operator_Metal : public Operator_sse
{
public:
	static Operator_Metal* New(unsigned int threads = 0);
	virtual Engine* CreateEngine();
	virtual unsigned int GetSetupThreads() const;
	virtual const std::vector<GeometryWinner>* GetGeometryWinners(GeometryWinnerType type, bool dualMesh) const;

protected:
	Operator_Metal();
	virtual bool Calc_EC();
	virtual bool CalcPEC();
	virtual bool SetupCSXGrid(CSRectGrid* grid);
	virtual bool CanReleaseECBeforeExtensions() const;
	virtual void CalcOperatorCoefficients();
	unsigned int m_setupThreads;

	//! Filled by the Metal geometry pass, consumed by EC-consuming extensions.
	std::vector<GeometryWinner> m_geoConductingSheet;
	std::vector<GeometryWinner> m_geoDispersivePrimal;
	std::vector<GeometryWinner> m_geoDispersiveDual;
	bool m_geoWinnersValid;
	//! Set once, so a geometry-winner fallback is reported only once.
	mutable bool m_geoWinnersWarned;
};

#endif // OPERATOR_METAL_H
