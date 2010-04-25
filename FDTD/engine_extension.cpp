#include "engine_extension.h"

#include "engine.h"

Engine_Extension::Engine_Extension(Operator_Extension* op_ext, Engine* eng)
{
	m_Op_ext = op_ext;
	m_Eng = eng;
}
