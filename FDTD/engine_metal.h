/*
* Copyright (C) 2026 openEMS contributors
*
* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*/

#ifndef ENGINE_METAL_H
#define ENGINE_METAL_H

#include "engine_sse.h"

class Engine_Metal : public Engine_sse
{
public:
	static Engine_Metal* New(const Operator_sse* op);
	virtual ~Engine_Metal();

	virtual void Init();
	virtual void Reset();

	virtual void DoPreVoltageUpdates();
	virtual void DoPostVoltageUpdates();
	virtual void DoPreCurrentUpdates();
	virtual void DoPostCurrentUpdates();

protected:
	Engine_Metal(const Operator_sse* op);
	virtual void UpdateVoltages(unsigned int startX, unsigned int numX);
	virtual void UpdateCurrents(unsigned int startX, unsigned int numX);

private:
	void InitUPML();
	void RunUPMLExtensions(bool voltage, bool pre);
	void FinishMetalCommands();

	struct MetalState;
	MetalState* m_Metal;
};

#endif // ENGINE_METAL_H
