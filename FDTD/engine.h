#ifndef ENGINE_H
#define ENGINE_H

class Engine
{
public:
    Engine();
	virtual ~Engine();

	//!Iterate a number of timesteps
	bool IterateTS(unsigned int numTS);

protected:

};

#endif // ENGINE_H
