// Temporary, opt-in diagnostics for the 2026-09-14 server overflow.
// This observes solver output; it must never clamp or modify physics state.
#pragma once

#include <Jolt/Physics/Body/Body.h>
#include <atomic>
#include <bit>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <string>
#include <type_traits>

JPH_NAMESPACE_BEGIN

namespace VelocityDiagnostics
{
inline std::atomic<bool> sRecorded { false };

inline const std::string &Path()
{
	static const std::string path = []()
	{
#ifdef JPH_COMPILER_MSVC
		char *value = nullptr;
		size_t size = 0;
		_dupenv_s(&value, &size, "VJOLT_DIAGNOSE_VELOCITY_FILE");
		std::string result = value? value : "";
		std::free(value);
		return result;
#else
		const char *value = std::getenv("VJOLT_DIAGNOSE_VELOCITY_FILE");
		return std::string(value? value : "");
#endif
	}();
	return path;
}

inline bool Enabled()
{
	return !Path().empty() && !sRecorded.load(std::memory_order_relaxed);
}

// Integer comparisons also recognize NaN/Inf without raising another FP exception.
inline bool Extreme(Vec3Arg inValue, float inThreshold = 1.0e6f)
{
	const uint32 threshold = std::bit_cast<uint32>(inThreshold);
	for (uint i = 0; i < 3; ++i)
		if ((std::bit_cast<uint32>(inValue[i]) & 0x7fffffffU) > threshold)
			return true;
	return false;
}

struct State
{
	Vec3 mV1, mW1, mV2, mW2;
	State(const Body &inBody1, const Body &inBody2)
		: mV1(inBody1.GetLinearVelocity()), mW1(inBody1.GetAngularVelocity()),
		  mV2(inBody2.GetLinearVelocity()), mW2(inBody2.GetAngularVelocity()) { }
	bool IsExtreme() const { return Extreme(mV1) || Extreme(mW1) || Extreme(mV2) || Extreme(mW2); }
};

inline void PrintVector(FILE *inFile, const char *inName, Vec3Arg inValue)
{
	// Exact float bits, including non-finite values, without float-to-double conversion.
	std::fprintf(inFile, "%s=%08x,%08x,%08x\n", inName,
		std::bit_cast<uint32>(inValue.GetX()), std::bit_cast<uint32>(inValue.GetY()),
		std::bit_cast<uint32>(inValue.GetZ()));
}

inline void Check(const char *inStage, uint32 inIndex, const Body &inBody1, const Body &inBody2,
	const State &inBefore, Vec3Arg inDisplacement = Vec3::sZero(), float inNormalLambda = 0.0f,
	const char *inEdge = "output", uint32 inSubType = ~uint32(0),
	void (*inWriteDetails)(FILE *, const void *) = nullptr, const void *inDetails = nullptr)
{
	const State after(inBody1, inBody2);
	if (!inBefore.IsExtreme() && !after.IsExtreme()
		&& !Extreme(inDisplacement, 1.0e12f)
		&& (std::bit_cast<uint32>(inNormalLambda) & 0x7fffffffU) <= std::bit_cast<uint32>(1.0e16f))
		return;
	// Only the first observed anomaly writes a file. No normal-path file IO or locks.
	if (sRecorded.exchange(true, std::memory_order_relaxed))
		return;
	FILE *file = std::fopen(Path().c_str(), "ab");
	if (file == nullptr)
		return;
	std::fprintf(file, "[DEBUG-velocity-0914] unix_time=%lld stage=%s edge=%s index=%u subtype=%u before_extreme=%d after_extreme=%d\n"
		"b1=%u user1=%llx b2=%u user2=%llx\nfloat_vectors_are_hex_bits\n",
		static_cast<long long>(std::time(nullptr)), inStage, inEdge, inIndex, inSubType,
		int(inBefore.IsExtreme()), int(after.IsExtreme()),
		inBody1.GetID().GetIndexAndSequenceNumber(), static_cast<unsigned long long>(inBody1.GetUserData()),
		inBody2.GetID().GetIndexAndSequenceNumber(), static_cast<unsigned long long>(inBody2.GetUserData()));
	PrintVector(file, "before_v1", inBefore.mV1); PrintVector(file, "before_w1", inBefore.mW1);
	PrintVector(file, "before_v2", inBefore.mV2); PrintVector(file, "before_w2", inBefore.mW2);
	PrintVector(file, "after_v1", after.mV1); PrintVector(file, "after_w1", after.mW1);
	PrintVector(file, "after_v2", after.mV2); PrintVector(file, "after_w2", after.mW2);
	PrintVector(file, "displacement", inDisplacement);
	std::fprintf(file, "normal_lambda=%08x\n", std::bit_cast<uint32>(inNormalLambda));
	if (inWriteDetails != nullptr)
		inWriteDetails(file, inDetails);
	std::fclose(file);
}

template <class Callback>
inline auto Observe(const char *inStage, uint32 inIndex, const Body &inBody1, const Body &inBody2,
	Callback &&inCallback, uint32 inSubType = ~uint32(0))
{
	if (!Enabled())
		return inCallback();
	const State before(inBody1, inBody2);
	if (before.IsExtreme())
		Check(inStage, inIndex, inBody1, inBody2, before, Vec3::sZero(), 0.0f, "input", inSubType);
	if constexpr (std::is_void_v<decltype(inCallback())>)
	{
		inCallback();
		Check(inStage, inIndex, inBody1, inBody2, before, Vec3::sZero(), 0.0f, "output", inSubType);
	}
	else
	{
		const auto result = inCallback();
		Check(inStage, inIndex, inBody1, inBody2, before, Vec3::sZero(), 0.0f, "output", inSubType);
		return result;
	}
}
}

JPH_NAMESPACE_END
