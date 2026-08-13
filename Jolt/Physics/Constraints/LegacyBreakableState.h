// SPDX-License-Identifier: MIT

#pragma once

#include <Jolt/Jolt.h>
#include <Jolt/Math/Mat44.h>

JPH_NAMESPACE_BEGIN

/// Exact scalar tests used by Havana's BreakableConstraintHelper. Limits are impulse magnitudes:
/// linear uses mass * |delta_v|^2, angular uses |I_body * delta_omega_body|^2.
inline bool LegacyBreaksLinear(float inMass, Vec3Arg inVelocityDelta, float inLinearLimit)
{
	return inMass * inVelocityDelta.LengthSq() >= Square(inLinearLimit);
}

inline Vec3 LegacyAngularMomentumDelta(
	QuatArg inBodyRotation,
	QuatArg inInertiaRotation,
	Vec3Arg inInverseInertiaDiagonal,
	Vec3Arg inAngularVelocityDelta)
{
	const Quat inertia_to_world = inBodyRotation * inInertiaRotation;
	const Vec3 inertia_space_delta = inertia_to_world.Conjugated() * inAngularVelocityDelta;
	const UVec4 zero_inverse_inertia = Vec3::sEquals(inInverseInertiaDiagonal, Vec3::sZero());
	const Vec3 denominator = Vec3::sSelect(inInverseInertiaDiagonal, Vec3::sOne(), zero_inverse_inertia);
	const Vec3 inverse_mask = Vec3::sSelect(Vec3::sOne(), Vec3::sZero(), zero_inverse_inertia);
	const Vec3 inertia_diagonal = inverse_mask / denominator;
	return inertia_to_world * (inertia_diagonal * inertia_space_delta);
}

inline bool LegacyBreaksAngular(
	QuatArg inBodyRotation,
	QuatArg inInertiaRotation,
	Vec3Arg inInverseInertiaDiagonal,
	Vec3Arg inAngularVelocityDelta,
	float inAngularLimit)
{
	return LegacyAngularMomentumDelta(inBodyRotation, inInertiaRotation, inInverseInertiaDiagonal, inAngularVelocityDelta).LengthSq()
		>= Square(inAngularLimit);
}

JPH_NAMESPACE_END
