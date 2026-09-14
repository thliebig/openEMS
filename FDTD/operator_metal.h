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

#include <vector>

class Operator_Metal : public Operator_sse
{
public:
	static Operator_Metal* New(unsigned int threads = 0);
	virtual Engine* CreateEngine();
	virtual unsigned int GetSetupThreads() const { return m_setupThreads; }
	virtual const std::vector<GeometryWinner>* GetGeometryWinners(GeometryWinnerType type, bool dualMesh) const;

protected:
	Operator_Metal();
	virtual bool Calc_EC();
	virtual bool CalcPEC();
	virtual bool CanReleaseECBeforeExtensions() const;
	virtual void CalcOperatorCoefficients();
	unsigned int m_setupThreads;

	//! Filled by the Metal geometry pass, consumed by EC-consuming extensions.
	std::vector<GeometryWinner> m_geoConductingSheet;
	std::vector<GeometryWinner> m_geoDispersivePrimal;
	std::vector<GeometryWinner> m_geoDispersiveDual;
	bool m_geoWinnersValid;
};

#endif // OPERATOR_METAL_H
