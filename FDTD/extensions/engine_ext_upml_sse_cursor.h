/*
 * Exact dynamic-type gate and experimental scalar cursor for the SSE UPML access map.
 * The type gate admits only exact SSE, compressed SSE, and multithread engine types.
 * The cursor caller supplies the exact start vector offset and strides from ArrayENG;
 * this helper only advances that vector/lane cursor along physical z cells.
 */
#ifndef ENGINE_EXT_UPML_SSE_CURSOR_H
#define ENGINE_EXT_UPML_SSE_CURSOR_H

#include <cstddef>
#include <typeinfo>

template <typename EngineSSE, typename EngineSSECompressed, typename EngineMultiThread>
inline bool IsExactUpmlSSECursorEngineType(const std::type_info& dynamicType)
{
	return dynamicType == typeid(EngineSSE) ||
		dynamicType == typeid(EngineSSECompressed) ||
		dynamicType == typeid(EngineMultiThread);
}

template <typename Vector4>
class Engine_Ext_UPML_SSE_Cursor
{
public:
	void BeginRow(
		Vector4* data,
		std::size_t vectorOffset,
		std::size_t vectorCount,
		std::size_t componentStride,
		std::size_t vectorStride,
		std::size_t vectorIndex,
		std::size_t lane
	)
	{
		m_data = data;
		m_vectorOffset = vectorOffset;
		m_vectorCount = vectorCount;
		m_componentStride = componentStride;
		m_vectorStride = vectorStride;
		m_vectorIndex = vectorIndex;
		m_lane = lane;
	}

	float& Component(std::size_t component) const
	{
		return m_data[m_vectorOffset + component * m_componentStride].f[m_lane];
	}

	void Advance()
	{
		++m_vectorIndex;
		if (m_vectorIndex == m_vectorCount)
		{
			m_vectorIndex = 0;
			m_vectorOffset -= (m_vectorCount - 1) * m_vectorStride;
			++m_lane;
		}
		else
		{
			m_vectorOffset += m_vectorStride;
		}
	}

private:
	Vector4* m_data = NULL;
	std::size_t m_vectorOffset = 0;
	std::size_t m_vectorCount = 0;
	std::size_t m_componentStride = 0;
	std::size_t m_vectorStride = 0;
	std::size_t m_vectorIndex = 0;
	std::size_t m_lane = 0;
};

#endif // ENGINE_EXT_UPML_SSE_CURSOR_H
