// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-FileCopyrightText: 2021 Jorrit Rouwe
// SPDX-License-Identifier: MIT

#include <Jolt/Jolt.h>

#include <Jolt/Physics/Constraints/ContactConstraintManager.h>
#include <Jolt/Physics/Constraints/CalculateSolverSteps.h>
#include <Jolt/Physics/Body/Body.h>
#include <Jolt/Physics/PhysicsUpdateContext.h>
#include <Jolt/Physics/PhysicsSettings.h>
#include <Jolt/Physics/PhysicsSystem.h>
#include <Jolt/Physics/IslandBuilder.h>
#include <Jolt/Physics/DeterminismLog.h>
#include <Jolt/Core/TempAllocator.h>
#include <Jolt/Core/QuickSort.h>
#include <Jolt/Core/Prefetch.h>
#include <cstdio>
#include <cstdlib>
#ifdef JPH_DEBUG_RENDERER
	#include <Jolt/Renderer/DebugRenderer.h>
#endif // JPH_DEBUG_RENDERER

JPH_NAMESPACE_BEGIN

using namespace literals;

static bool IsPlayerSupportSpreadPenetrationEnabled()
{
	static const bool enabled = std::getenv("VJOLT_TEST_PLAYER_SUPPORT_SPREAD_PENETRATION") != nullptr;
	return enabled;
}

static bool IsPlayerSupportCacheTraceEnabled()
{
	static const bool enabled = std::getenv("VJOLT_TRACE_PLAYER_SUPPORT_CACHE_ACTIVE") != nullptr;
	return enabled;
}

static bool IsVPhysicsImpactTraceEnabled()
{
	static const bool enabled = std::getenv("VJOLT_TRACE_VPHYSICS_IMPACT") != nullptr;
	return enabled;
}

static bool IsLegacyFrictionTraceEnabled()
{
	static const bool enabled = std::getenv("VJOLT_TRACE_LEGACY_FRICTION") != nullptr;
	return enabled;
}

#ifdef JPH_DEBUG_RENDERER
bool ContactConstraintManager::sDrawContactPoint = false;
bool ContactConstraintManager::sDrawSupportingFaces = false;
bool ContactConstraintManager::sDrawContactPointReduction = false;
bool ContactConstraintManager::sDrawContactManifolds = false;
#endif // JPH_DEBUG_RENDERER

//#define JPH_MANIFOLD_CACHE_DEBUG

////////////////////////////////////////////////////////////////////////////////////////////////////////
// ContactConstraintManager::WorldContactPoint
////////////////////////////////////////////////////////////////////////////////////////////////////////

template <EMotionType Type1, EMotionType Type2>
JPH_INLINE void ContactConstraintManager::WorldContactPoint<Type1, Type2>::CalculateFrictionAndNonPenetrationConstraintProperties(float inDeltaTime, Vec3Arg inGravity, const Body &inBody1, const Body &inBody2, float inInvM1, float inInvM2, Mat44Arg inInvI1, Mat44Arg inInvI2, RVec3Arg inWorldSpacePosition1, RVec3Arg inWorldSpacePosition2, Vec3Arg inWorldSpaceNormal, Vec3Arg inWorldSpaceTangent1, Vec3Arg inWorldSpaceTangent2, const ContactSettings &inSettings, float inMinVelocityForRestitution, float inContactSeparationTarget, float inContactVelocityBiasFactorOutsideTarget, float inContactVelocityBiasFactorInsideTarget, float inContactRescueVelocityFactor, bool inAllowContactRescue, bool inAllowRestitution)
{
	JPH_DET_LOG("CalculateFrictionAndNonPenetrationConstraintProperties: p1: " << inWorldSpacePosition1 << " p2: " << inWorldSpacePosition2
		<< " normal: " << inWorldSpaceNormal << " tangent1: " << inWorldSpaceTangent1 << " tangent2: " << inWorldSpaceTangent2
		<< " restitution: " << inSettings.mCombinedRestitution << " friction: " << inSettings.mCombinedFriction << " minv: " << inMinVelocityForRestitution
		<< " surface_vel: " << inSettings.mRelativeLinearSurfaceVelocity << " surface_ang: " << inSettings.mRelativeAngularSurfaceVelocity);

	// Calculate collision points relative to body
	RVec3 p = 0.5_r * (inWorldSpacePosition1 + inWorldSpacePosition2);
	Vec3 r1 = Vec3(p - inBody1.GetCenterOfMassPosition());
	Vec3 r2 = Vec3(p - inBody2.GetCenterOfMassPosition());

	const MotionProperties *mp1 = inBody1.GetMotionPropertiesUnchecked();
	const MotionProperties *mp2 = inBody2.GetMotionPropertiesUnchecked();

	// Calculate velocity of collision points
	Vec3 relative_velocity;
	if constexpr (Type1 != EMotionType::Static && Type2 != EMotionType::Static)
		relative_velocity = mp2->GetPointVelocityCOM(r2) - mp1->GetPointVelocityCOM(r1);
	else if constexpr (Type1 != EMotionType::Static)
		relative_velocity = -mp1->GetPointVelocityCOM(r1);
	else if constexpr (Type2 != EMotionType::Static)
		relative_velocity = mp2->GetPointVelocityCOM(r2);
	else
	{
		JPH_ASSERT(false, "Static vs static makes no sense");
		relative_velocity = Vec3::sZero();
	}
	float normal_velocity = relative_velocity.Dot(inWorldSpaceNormal);
	const float contact_gap = Vec3(inWorldSpacePosition2 - inWorldSpacePosition1).Dot(inWorldSpaceNormal);

	// IVP distinguishes a pending mindist from an established friction contact. A point without
	// cached support may live in the broad contact range, but it only constrains motion against the
	// inner collision distance (half the persistent friction target). Once it has carried a normal
	// impulse, the ordinary friction-distance target and its gentle convergence rates take over.
	const bool disable_penetration_recovery = inSettings.mVPhysicsDisablePenetrationRecovery;
	const bool pending_legacy_contact = !disable_penetration_recovery
		&& inAllowContactRescue && inContactRescueVelocityFactor > 0.0f;
	const float solver_separation_target = pending_legacy_contact
		? 0.5f * inContactSeparationTarget
		: inContactSeparationTarget;
	// Signed distance from the active velocity target (> 0 if outside / farther apart, < 0 if inside).
	float separation = Vec3(inWorldSpacePosition2 - inWorldSpacePosition1).Dot(inWorldSpaceNormal) - solver_separation_target;

	// Outside the target, native Jolt closes the entire gap in one step. Legacy friction contacts use
	// a fixed per-second convergence rate instead. They may also request a much gentler outward rate
	// when already inside the target rather than leaving all correction to the position solver.
	// See: GDC 2013: "Physics for Game Programmers; Continuous Collision" - Erin Catto
	float speculative_contact_velocity_bias;
	if (separation > 0.0f)
	{
		float outside_factor = !pending_legacy_contact && inContactVelocityBiasFactorOutsideTarget > 0.0f
			? inContactVelocityBiasFactorOutsideTarget
			: 1.0f / inDeltaTime;
		speculative_contact_velocity_bias = separation * outside_factor;
	}
	else
	{
		speculative_contact_velocity_bias = separation
			* (pending_legacy_contact || disable_penetration_recovery
				? 0.0f : inContactVelocityBiasFactorInsideTarget);
	}

	// Determine if the velocity is big enough for restitution
	float normal_velocity_bias;
	if (inAllowRestitution && inSettings.mCombinedRestitution > 0.0f && normal_velocity < -inMinVelocityForRestitution)
	{
		// We have a velocity that is big enough for restitution. This is where speculative contacts don't work
		// great as we have to decide now if we're going to apply the restitution or not. If the relative
		// velocity is big enough for a hit, we apply the restitution (in the end, due to other constraints,
		// the objects may actually not collide and we will have applied restitution incorrectly). Another
		// artifact that occurs because of this approximation is that the object will bounce from its current
		// position rather than from a position where it is touching the other object. This causes the object
		// to appear to move faster for 1 frame (the opposite of time stealing).
		if (normal_velocity < -speculative_contact_velocity_bias)
		{
			// The gravity / constant forces are applied in the beginning of the time step.
			// If we get here, there was a collision at the beginning of the time step, so we've applied too much force.
			// This means that our calculated restitution can be too high resulting in an increase in energy.
			// So, when we apply restitution, we cancel the added velocity due to these forces.
			Vec3 relative_acceleration;

			// Calculate effect of gravity
			if constexpr (Type1 != EMotionType::Static && Type2 != EMotionType::Static)
				relative_acceleration = inGravity * (mp2->GetGravityFactor() - mp1->GetGravityFactor());
			else if constexpr (Type1 != EMotionType::Static)
				relative_acceleration = -inGravity * mp1->GetGravityFactor();
			else if constexpr (Type2 != EMotionType::Static)
				relative_acceleration = inGravity * mp2->GetGravityFactor();
			else
			{
				JPH_ASSERT(false, "Static vs static makes no sense");
				relative_acceleration = Vec3::sZero();
			}

			// Calculate effect of accumulated forces
			if constexpr (Type1 == EMotionType::Dynamic)
				relative_acceleration -= mp1->GetAccumulatedForce() * mp1->GetInverseMass();
			if constexpr (Type2 == EMotionType::Dynamic)
				relative_acceleration += mp2->GetAccumulatedForce() * mp2->GetInverseMass();

			// We only compensate forces towards the contact normal.
			float force_delta_velocity = min(0.0f, relative_acceleration.Dot(inWorldSpaceNormal) * inDeltaTime);

			normal_velocity_bias = inSettings.mCombinedRestitution * (normal_velocity - force_delta_velocity);
		}
		else
		{
			// In this case we have predicted that we don't hit the other object, but if we do (due to other constraints changing velocities)
			// the speculative contact will prevent penetration but will not apply restitution leading to another artifact.
			normal_velocity_bias = speculative_contact_velocity_bias;
		}
	}
	else
	{
		// No restitution. We can safely apply our contact velocity bias.
		normal_velocity_bias = speculative_contact_velocity_bias;
	}

	// IVP resolves a newly impacting contact in the velocity domain instead of relying solely on
	// split-position correction. A zero cached normal lambda means this point has not supported an
	// active contact yet; this also lets a speculative manifold wait until it crosses the inner
	// collision tolerance before consuming its one impact rescue.
	const bool player_support_panic_rescue = inSettings.mVPhysicsPlayerSupportStartupContact
		&& contact_gap < 0.0f && normal_velocity >= 0.0f;
	if (!disable_penetration_recovery && inAllowContactRescue
		&& (normal_velocity < 0.0f || player_support_panic_rescue)
		&& inContactRescueVelocityFactor > 0.0f && inContactSeparationTarget > 0.0f && inDeltaTime > 0.0f)
	{
		// An exact-touch player support manifold has already consumed IVP's inner
		// collision-distance lead-in. During startup, charge rescue only for real
		// penetration below geometric contact instead of adding the full tolerance again.
		float collision_tolerance = inSettings.mVPhysicsPlayerSupportStartupContact
			? 0.0f : 0.5f * inContactSeparationTarget;
		if (contact_gap <= collision_tolerance)
		{
			const float startup_rescue_penetration =
				inSettings.mVPhysicsPlayerSupportStartupRescuePenetration >= 0.0f
					? inSettings.mVPhysicsPlayerSupportStartupRescuePenetration
					: inSettings.mVPhysicsPlayerSupportStartupPenetration;
			const float raw_rescue_penetration = inSettings.mVPhysicsPlayerSupportStartupContact
				&& startup_rescue_penetration >= 0.0f
				? startup_rescue_penetration
				: max(0.0f, collision_tolerance - contact_gap);
			// IVP's ordinary impact rescue is evaluated by an event-time mindist at the
			// inner collision boundary. A Jolt manifold can instead be discovered after
			// the complete discrete step and report a much deeper overlap. Feeding that
			// post-event depth into 2 * depth / dt injects energy that IVP would route to
			// its bounded anomaly manager instead. Keep the ordinary rescue within one
			// collision-tolerance interval; deep overlap has its own explicit path.
			const float rescue_penetration = min(
				0.5f * inContactSeparationTarget, raw_rescue_penetration);
			float rescue_addon = 2.0f * rescue_penetration / inDeltaTime;
			float rescue_velocity = inContactRescueVelocityFactor * (inContactSeparationTarget + rescue_addon);
			normal_velocity_bias = min(normal_velocity_bias,
				player_support_panic_rescue
					? -(normal_velocity + rescue_velocity)
					: -rescue_velocity);
		}
	}

	mNonPenetrationConstraint.CalculateConstraintProperties(inInvM1, inInvI1, r1, inInvM2, inInvI2, r2, inWorldSpaceNormal, normal_velocity_bias);

	// Calculate friction part
	if (inSettings.mCombinedFriction > 0.0f)
	{
		// Get surface velocity relative to tangents
		Vec3 ws_surface_velocity = inSettings.mRelativeLinearSurfaceVelocity + inSettings.mRelativeAngularSurfaceVelocity.Cross(r1);
		float surface_velocity1 = inWorldSpaceTangent1.Dot(ws_surface_velocity);
		float surface_velocity2 = inWorldSpaceTangent2.Dot(ws_surface_velocity);

		// Implement friction as 2 ContactConstraintParts
		mFrictionConstraint1.CalculateConstraintProperties(inInvM1, inInvI1, r1, inInvM2, inInvI2, r2, inWorldSpaceTangent1, surface_velocity1);
		mFrictionConstraint2.CalculateConstraintProperties(inInvM1, inInvI1, r1, inInvM2, inInvI2, r2, inWorldSpaceTangent2, surface_velocity2);
	}
	else
	{
		// Turn off friction constraint
		mFrictionConstraint1.Deactivate();
		mFrictionConstraint2.Deactivate();
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
// ContactConstraintManager::ContactConstraint
////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef JPH_DEBUG_RENDERER
template <EMotionType Type1, EMotionType Type2>
void ContactConstraintManager::ContactConstraint<Type1, Type2>::Draw(DebugRenderer *inRenderer, ColorArg inManifoldColor) const
{
	if (mNumContactPoints == 0)
		return;

	// Get body transforms
	RMat44 transform_body1 = mBody1->GetCenterOfMassTransform();
	RMat44 transform_body2 = mBody2->GetCenterOfMassTransform();

	RVec3 prev_point = transform_body1 * Vec3::sLoadFloat3Unsafe(mContactPoints[mNumContactPoints - 1].mContactPoint->mPosition1);
	for (uint32 i = 0; i < mNumContactPoints; ++i)
	{
		const WorldContactPoint<Type1, Type2> &wcp = mContactPoints[i];

		// Test if any lambda from the previous frame was transferred
		float radius = wcp.mNonPenetrationConstraint.GetTotalLambda() == 0.0f
					&& wcp.mFrictionConstraint1.GetTotalLambda() == 0.0f
					&& wcp.mFrictionConstraint2.GetTotalLambda() == 0.0f? 0.1f :  0.2f;

		RVec3 next_point = transform_body1 * Vec3::sLoadFloat3Unsafe(wcp.mContactPoint->mPosition1);
		inRenderer->DrawMarker(next_point, Color::sCyan, radius);
		inRenderer->DrawMarker(transform_body2 * Vec3::sLoadFloat3Unsafe(wcp.mContactPoint->mPosition2), Color::sPurple, radius);

		// Draw edge
		inRenderer->DrawArrow(prev_point, next_point, inManifoldColor, 0.05f);
		prev_point = next_point;
	}

	// Draw normal
	RVec3 wp = transform_body1 * Vec3::sLoadFloat3Unsafe(mContactPoints[0].mContactPoint->mPosition1);
	inRenderer->DrawArrow(wp, wp + GetWorldSpaceNormal(), Color::sRed, 0.05f);

	// Get tangents
	Vec3 t1, t2;
	GetTangents(t1, t2);

	// Draw tangents
	inRenderer->DrawLine(wp, wp + t1, Color::sGreen);
	inRenderer->DrawLine(wp, wp + t2, Color::sBlue);
}
#endif // JPH_DEBUG_RENDERER

////////////////////////////////////////////////////////////////////////////////////////////////////////
// ContactConstraintManager::CachedContactPoint
////////////////////////////////////////////////////////////////////////////////////////////////////////

void ContactConstraintManager::CachedContactPoint::SaveState(StateRecorder &inStream) const
{
	inStream.Write(mPosition1);
	inStream.Write(mPosition2);
	inStream.Write(mNonPenetrationLambda);
	inStream.Write(mFrictionLambda);
	inStream.Write(mVPhysicsFrictionDisplacement);
	inStream.Write(mVPhysicsFrictionNormalLambda);
	inStream.Write(mVPhysicsContactEstablished);
}

void ContactConstraintManager::CachedContactPoint::RestoreState(StateRecorder &inStream)
{
	inStream.Read(mPosition1);
	inStream.Read(mPosition2);
	inStream.Read(mNonPenetrationLambda);
	inStream.Read(mFrictionLambda);
	inStream.Read(mVPhysicsFrictionDisplacement);
	inStream.Read(mVPhysicsFrictionNormalLambda);
	inStream.Read(mVPhysicsContactEstablished);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
// ContactConstraintManager::CachedManifold
////////////////////////////////////////////////////////////////////////////////////////////////////////

void ContactConstraintManager::CachedManifold::SaveState(StateRecorder &inStream) const
{
	inStream.Write(mContactNormal);
	inStream.Write(mVPhysicsPlayerSupportStartupPenetration);
}

void ContactConstraintManager::CachedManifold::RestoreState(StateRecorder &inStream)
{
	inStream.Read(mContactNormal);
	inStream.Read(mVPhysicsPlayerSupportStartupPenetration);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
// ContactConstraintManager::CachedBodyPair
////////////////////////////////////////////////////////////////////////////////////////////////////////

void ContactConstraintManager::CachedBodyPair::SaveState(StateRecorder &inStream) const
{
	inStream.Write(mDeltaPosition);
	inStream.Write(mDeltaRotation);
}

void ContactConstraintManager::CachedBodyPair::RestoreState(StateRecorder &inStream)
{
	inStream.Read(mDeltaPosition);
	inStream.Read(mDeltaRotation);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
// ContactConstraintManager::ManifoldCache
////////////////////////////////////////////////////////////////////////////////////////////////////////

void ContactConstraintManager::ManifoldCache::Init(uint inMaxBodyPairs, uint inMaxContactConstraints, uint inCachedManifoldsSize)
{
	JPH_ASSERT(inMaxContactConstraints <= cMaxContactConstraintsLimit); // Should have been enforced by caller

	uint max_body_pairs = min(inMaxBodyPairs, cMaxBodyPairsLimit);
	JPH_ASSERT(max_body_pairs == inMaxBodyPairs, "Cannot support this many body pairs!");
	max_body_pairs = max(max_body_pairs, 4u); // Because our hash map requires at least 4 buckets, we need to have a minimum number of body pairs

	mAllocator.Init(uint(min(uint64(max_body_pairs) * sizeof(BPKeyValue) + inCachedManifoldsSize, uint64(~uint(0)))));

	mCachedManifolds.Init(GetNextPowerOf2(inMaxContactConstraints));
	mCachedBodyPairs.Init(GetNextPowerOf2(max_body_pairs));
}

void ContactConstraintManager::ManifoldCache::Clear()
{
	JPH_PROFILE_FUNCTION();

	mCachedManifolds.Clear();
	mCachedBodyPairs.Clear();
	mAllocator.Clear();

#ifdef JPH_ENABLE_ASSERTS
	// Mark as incomplete
	mIsFinalized = false;
#endif
}

void ContactConstraintManager::ManifoldCache::Prepare(uint inExpectedNumBodyPairs, uint inExpectedNumManifolds)
{
	// Minimum amount of buckets to use in the hash map
	constexpr uint32 cMinBuckets = 1024;

	// Use the next higher power of 2 of amount of objects in the cache from last frame to determine the amount of buckets in this frame
	mCachedManifolds.SetNumBuckets(min(max(cMinBuckets, GetNextPowerOf2(inExpectedNumManifolds)), mCachedManifolds.GetMaxBuckets()));
	mCachedBodyPairs.SetNumBuckets(min(max(cMinBuckets, GetNextPowerOf2(inExpectedNumBodyPairs)), mCachedBodyPairs.GetMaxBuckets()));
}

const ContactConstraintManager::MKeyValue *ContactConstraintManager::ManifoldCache::Find(const SubShapeIDPair &inKey, uint64 inKeyHash) const
{
	JPH_ASSERT(mIsFinalized);
	return mCachedManifolds.Find(inKey, inKeyHash);
}

ContactConstraintManager::MKeyValue *ContactConstraintManager::ManifoldCache::Create(ContactAllocator &ioContactAllocator, const SubShapeIDPair &inKey, uint64 inKeyHash, int inNumContactPoints)
{
	JPH_ASSERT(!mIsFinalized);
	MKeyValue *kv = mCachedManifolds.Create(ioContactAllocator, inKey, inKeyHash, CachedManifold::sGetRequiredExtraSize(inNumContactPoints));
	if (kv == nullptr)
	{
		ioContactAllocator.mErrors |= EPhysicsUpdateError::ManifoldCacheFull;
		return nullptr;
	}
	kv->GetValue().mNumContactPoints = uint16(inNumContactPoints);
	++ioContactAllocator.mNumManifolds;
	return kv;
}

ContactConstraintManager::MKVAndCreated ContactConstraintManager::ManifoldCache::FindOrCreate(ContactAllocator &ioContactAllocator, const SubShapeIDPair &inKey, uint64 inKeyHash, int inNumContactPoints)
{
	MKeyValue *kv = const_cast<MKeyValue *>(mCachedManifolds.Find(inKey, inKeyHash));
	if (kv != nullptr)
		return { kv, false };

	return { Create(ioContactAllocator, inKey, inKeyHash, inNumContactPoints), true };
}

uint32 ContactConstraintManager::ManifoldCache::ToHandle(const MKeyValue *inKeyValue) const
{
	JPH_ASSERT(!mIsFinalized);
	return mCachedManifolds.ToHandle(inKeyValue);
}

const ContactConstraintManager::MKeyValue *ContactConstraintManager::ManifoldCache::FromHandle(uint32 inHandle) const
{
	JPH_ASSERT(mIsFinalized);
	return mCachedManifolds.FromHandle(inHandle);
}

const ContactConstraintManager::BPKeyValue *ContactConstraintManager::ManifoldCache::Find(const BodyPair &inKey, uint64 inKeyHash) const
{
	JPH_ASSERT(mIsFinalized);
	return mCachedBodyPairs.Find(inKey, inKeyHash);
}

ContactConstraintManager::BPKeyValue *ContactConstraintManager::ManifoldCache::Create(ContactAllocator &ioContactAllocator, const BodyPair &inKey, uint64 inKeyHash)
{
	JPH_ASSERT(!mIsFinalized);
	BPKeyValue *kv = mCachedBodyPairs.Create(ioContactAllocator, inKey, inKeyHash, 0);
	if (kv == nullptr)
	{
		ioContactAllocator.mErrors |= EPhysicsUpdateError::BodyPairCacheFull;
		return nullptr;
	}
	++ioContactAllocator.mNumBodyPairs;
	return kv;
}

void ContactConstraintManager::ManifoldCache::GetAllBodyPairsSorted(Array<const BPKeyValue *> &outAll) const
{
	JPH_ASSERT(mIsFinalized);
	mCachedBodyPairs.GetAllKeyValues(outAll);

	// Sort by key
	QuickSort(outAll.begin(), outAll.end(), [](const BPKeyValue *inLHS, const BPKeyValue *inRHS) {
		return inLHS->GetKey() < inRHS->GetKey();
	});
}

void ContactConstraintManager::ManifoldCache::GetAllManifoldsSorted(const CachedBodyPair &inBodyPair, Array<const MKeyValue *> &outAll) const
{
	JPH_ASSERT(mIsFinalized);

	// Iterate through the attached manifolds
	for (uint32 handle = inBodyPair.mFirstCachedManifold; handle != ManifoldMap::cInvalidHandle; handle = FromHandle(handle)->GetValue().mNextWithSameBodyPair)
	{
		const MKeyValue *kv = mCachedManifolds.FromHandle(handle);
		outAll.push_back(kv);
	}

	// Sort by key
	QuickSort(outAll.begin(), outAll.end(), [](const MKeyValue *inLHS, const MKeyValue *inRHS) {
		return inLHS->GetKey() < inRHS->GetKey();
	});
}

void ContactConstraintManager::ManifoldCache::GetAllCCDManifoldsSorted(Array<const MKeyValue *> &outAll) const
{
	mCachedManifolds.GetAllKeyValues(outAll);

	for (int i = (int)outAll.size() - 1; i >= 0; --i)
		if ((outAll[i]->GetValue().mFlags & (uint16)CachedManifold::EFlags::CCDContact) == 0)
		{
			outAll[i] = outAll.back();
			outAll.pop_back();
		}

	// Sort by key
	QuickSort(outAll.begin(), outAll.end(), [](const MKeyValue *inLHS, const MKeyValue *inRHS) {
		return inLHS->GetKey() < inRHS->GetKey();
	});
}

void ContactConstraintManager::ManifoldCache::ContactPointRemovedCallbacks(ContactListener *inListener)
{
	JPH_PROFILE_FUNCTION();

	for (MKeyValue &kv : mCachedManifolds)
		if ((kv.GetValue().mFlags & uint16(CachedManifold::EFlags::ContactPersisted)) == 0)
			inListener->OnContactRemoved(kv.GetKey());
}

#ifdef JPH_ENABLE_ASSERTS

void ContactConstraintManager::ManifoldCache::Finalize()
{
	mIsFinalized = true;

#ifdef JPH_MANIFOLD_CACHE_DEBUG
	Trace("ManifoldMap:");
	mCachedManifolds.TraceStats();
	Trace("BodyPairMap:");
	mCachedBodyPairs.TraceStats();
#endif // JPH_MANIFOLD_CACHE_DEBUG
}

#endif

void ContactConstraintManager::ManifoldCache::SaveState(StateRecorder &inStream, const StateRecorderFilter *inFilter) const
{
	JPH_ASSERT(mIsFinalized);

	// Get contents of cache
	Array<const BPKeyValue *> all_bp;
	GetAllBodyPairsSorted(all_bp);

	// Determine which ones to save
	Array<const BPKeyValue *> selected_bp;
	if (inFilter == nullptr)
		selected_bp = std::move(all_bp);
	else
	{
		selected_bp.reserve(all_bp.size());
		for (const BPKeyValue *bp_kv : all_bp)
			if (inFilter->ShouldSaveContact(bp_kv->GetKey().mBodyA, bp_kv->GetKey().mBodyB))
				selected_bp.push_back(bp_kv);
	}

	// Write body pairs
	uint32 num_body_pairs = uint32(selected_bp.size());
	inStream.Write(num_body_pairs);
	for (const BPKeyValue *bp_kv : selected_bp)
	{
		// Write body pair key
		inStream.Write(bp_kv->GetKey());

		// Write body pair
		const CachedBodyPair &bp = bp_kv->GetValue();
		bp.SaveState(inStream);

		// Get attached manifolds
		Array<const MKeyValue *> all_m;
		GetAllManifoldsSorted(bp, all_m);

		// Write num manifolds
		uint32 num_manifolds = uint32(all_m.size());
		inStream.Write(num_manifolds);

		// Write all manifolds
		for (const MKeyValue *m_kv : all_m)
		{
			// Write key
			inStream.Write(m_kv->GetKey());
			const CachedManifold &cm = m_kv->GetValue();
			JPH_ASSERT((cm.mFlags & (uint16)CachedManifold::EFlags::CCDContact) == 0);

			// Write amount of contacts
			inStream.Write(cm.mNumContactPoints);

			// Write manifold
			cm.SaveState(inStream);

			// Write contact points
			for (uint32 i = 0; i < cm.mNumContactPoints; ++i)
				cm.mContactPoints[i].SaveState(inStream);
		}
	}

	// Get CCD manifolds
	Array<const MKeyValue *> all_m;
	GetAllCCDManifoldsSorted(all_m);

	// Determine which ones to save
	Array<const MKeyValue *> selected_m;
	if (inFilter == nullptr)
		selected_m = std::move(all_m);
	else
	{
		selected_m.reserve(all_m.size());
		for (const MKeyValue *m_kv : all_m)
			if (inFilter->ShouldSaveContact(m_kv->GetKey().GetBody1ID(), m_kv->GetKey().GetBody2ID()))
				selected_m.push_back(m_kv);
	}

	// Write all CCD manifold keys
	uint32 num_manifolds = uint32(selected_m.size());
	inStream.Write(num_manifolds);
	for (const MKeyValue *m_kv : selected_m)
		inStream.Write(m_kv->GetKey());
}

bool ContactConstraintManager::ManifoldCache::RestoreState(const ManifoldCache &inReadCache, StateRecorder &inStream, const StateRecorderFilter *inFilter)
{
	JPH_ASSERT(!mIsFinalized);

	bool success = true;

	// Create a contact allocator for restoring the contact cache
	ContactAllocator contact_allocator(GetContactAllocator());

	// When validating, get all existing body pairs
	Array<const BPKeyValue *> all_bp;
	if (inStream.IsValidating())
		inReadCache.GetAllBodyPairsSorted(all_bp);

	// Read amount of body pairs
	uint32 num_body_pairs;
	if (inStream.IsValidating())
		num_body_pairs = uint32(all_bp.size());
	inStream.Read(num_body_pairs);

	// Read entire cache
	for (uint32 i = 0; i < num_body_pairs; ++i)
	{
		// Read key
		BodyPair body_pair_key;
		if (inStream.IsValidating() && i < all_bp.size())
			body_pair_key = all_bp[i]->GetKey();
		inStream.Read(body_pair_key);

		// Check if we want to restore this contact
		if (inFilter == nullptr || inFilter->ShouldRestoreContact(body_pair_key.mBodyA, body_pair_key.mBodyB))
		{
			// Create new entry for this body pair
			uint64 body_pair_hash = body_pair_key.GetHash();
			BPKeyValue *bp_kv = Create(contact_allocator, body_pair_key, body_pair_hash);
			if (bp_kv == nullptr)
			{
				// Out of cache space
				success = false;
				break;
			}
			CachedBodyPair &bp = bp_kv->GetValue();

			// Read body pair
			if (inStream.IsValidating() && i < all_bp.size())
				memcpy(&bp, &all_bp[i]->GetValue(), sizeof(CachedBodyPair));
			bp.RestoreState(inStream);

			// When validating, get all existing manifolds
			Array<const MKeyValue *> all_m;
			if (inStream.IsValidating() && i < all_bp.size())
				inReadCache.GetAllManifoldsSorted(all_bp[i]->GetValue(), all_m);

			// Read amount of manifolds
			uint32 num_manifolds = 0;
			if (inStream.IsValidating())
				num_manifolds = uint32(all_m.size());
			inStream.Read(num_manifolds);

			uint32 handle = ManifoldMap::cInvalidHandle;
			for (uint32 j = 0; j < num_manifolds; ++j)
			{
				// Read key
				SubShapeIDPair sub_shape_key;
				if (inStream.IsValidating() && j < all_m.size())
					sub_shape_key = all_m[j]->GetKey();
				inStream.Read(sub_shape_key);
				uint64 sub_shape_key_hash = sub_shape_key.GetHash();

				// Read amount of contact points
				uint16 num_contact_points = 0;
				if (inStream.IsValidating() && j < all_m.size())
					num_contact_points = all_m[j]->GetValue().mNumContactPoints;
				inStream.Read(num_contact_points);

				// Read manifold
				MKeyValue *m_kv = Create(contact_allocator, sub_shape_key, sub_shape_key_hash, num_contact_points);
				if (m_kv == nullptr)
				{
					// Out of cache space
					success = false;
					break;
				}
				CachedManifold &cm = m_kv->GetValue();
				if (inStream.IsValidating() && j < all_m.size())
				{
					memcpy(&cm, &all_m[j]->GetValue(), CachedManifold::sGetRequiredTotalSize(num_contact_points));
					cm.mNumContactPoints = uint16(num_contact_points); // Restore num contact points
				}
				cm.RestoreState(inStream);
				cm.mNextWithSameBodyPair = handle;
				handle = ToHandle(m_kv);

				// Read contact points
				for (uint32 k = 0; k < num_contact_points; ++k)
					cm.mContactPoints[k].RestoreState(inStream);
			}
			bp.mFirstCachedManifold = handle;
		}
		else
		{
			// Skip the contact
			CachedBodyPair bp;
			bp.RestoreState(inStream);
			uint32 num_manifolds = 0;
			inStream.Read(num_manifolds);
			for (uint32 j = 0; j < num_manifolds; ++j)
			{
				SubShapeIDPair sub_shape_key;
				inStream.Read(sub_shape_key);
				uint16 num_contact_points;
				inStream.Read(num_contact_points);
				CachedManifold cm;
				cm.RestoreState(inStream);
				for (uint32 k = 0; k < num_contact_points; ++k)
					cm.mContactPoints[0].RestoreState(inStream);
			}
		}
	}

	// When validating, get all existing CCD manifolds
	Array<const MKeyValue *> all_m;
	if (inStream.IsValidating())
		inReadCache.GetAllCCDManifoldsSorted(all_m);

	// Read amount of CCD manifolds
	uint32 num_manifolds;
	if (inStream.IsValidating())
		num_manifolds = uint32(all_m.size());
	inStream.Read(num_manifolds);

	for (uint32 j = 0; j < num_manifolds; ++j)
	{
		// Read key
		SubShapeIDPair sub_shape_key;
		if (inStream.IsValidating() && j < all_m.size())
			sub_shape_key = all_m[j]->GetKey();
		inStream.Read(sub_shape_key);

		// Check if we want to restore this contact
		if (inFilter == nullptr || inFilter->ShouldRestoreContact(sub_shape_key.GetBody1ID(), sub_shape_key.GetBody2ID()))
		{
			// Create CCD manifold
			uint64 sub_shape_key_hash = sub_shape_key.GetHash();
			MKeyValue *m_kv = Create(contact_allocator, sub_shape_key, sub_shape_key_hash, 0);
			if (m_kv == nullptr)
			{
				// Out of cache space
				success = false;
				break;
			}
			CachedManifold &cm = m_kv->GetValue();
			cm.mFlags |= (uint16)CachedManifold::EFlags::CCDContact;
		}
	}

#ifdef JPH_ENABLE_ASSERTS
	// We don't finalize until the last part is restored
	if (inStream.IsLastPart())
		mIsFinalized = true;
#endif

	return success;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
// ContactConstraintManager
////////////////////////////////////////////////////////////////////////////////////////////////////////

ContactConstraintManager::ContactConstraintManager(const PhysicsSettings &inPhysicsSettings) :
	mPhysicsSettings(inPhysicsSettings)
{
#ifdef JPH_ENABLE_ASSERTS
	// For the first frame mark this empty buffer as finalized
	mCache[mCacheWriteIdx ^ 1].Finalize();
#endif
}

ContactConstraintManager::~ContactConstraintManager()
{
	JPH_ASSERT(mConstraints == nullptr);
	JPH_ASSERT(mConstraintIdxToOffset == nullptr);
}

void ContactConstraintManager::Init(uint inMaxBodyPairs, uint inMaxContactConstraints)
{
	// Limit the number of constraints so that the allocation size fits in an unsigned integer
	mMaxConstraints = min(inMaxContactConstraints, cMaxContactConstraintsLimit);
	JPH_ASSERT(mMaxConstraints == inMaxContactConstraints, "Cannot support this many contact constraints!");
	mMaxConstraints = max(mMaxConstraints, 4u); // Because our hash map requires at least 4 buckets, we need to have a minimum number of constraints

	// Calculate worst case cache usage
	constexpr uint cMaxManifoldSizePerConstraint = sizeof(MKeyValue) + CachedManifold::sGetRequiredExtraSize(MaxContactPoints);
	static_assert(cMaxManifoldSizePerConstraint < cMaxConstraintSize); // If not true, then the next line can overflow
	uint cached_manifolds_size = mMaxConstraints * cMaxManifoldSizePerConstraint;

	// Init the caches
	mCache[0].Init(inMaxBodyPairs, mMaxConstraints, cached_manifolds_size);
	mCache[1].Init(inMaxBodyPairs, mMaxConstraints, cached_manifolds_size);
}

void ContactConstraintManager::PrepareConstraintBuffer(PhysicsUpdateContext *inContext)
{
	// Store context
	mUpdateContext = inContext;

	// Allocate temporary constraint buffer
	JPH_ASSERT(mConstraints == nullptr);
	mConstraints = (uint8 *)inContext->mTempAllocator->Allocate(mMaxConstraints * cMaxConstraintSize);
	JPH_ASSERT(mConstraintIdxToOffset == nullptr);
	mConstraintIdxToOffset = (uint32 *)inContext->mTempAllocator->Allocate(mMaxConstraints * sizeof(uint32));
}

template <EMotionType Type1, EMotionType Type2>
JPH_INLINE ContactConstraintManager::ContactConstraint<Type1, Type2> *ContactConstraintManager::CreateConstraint(bool &ioActivateAndLinkBodies, Body &inBody1, Body &inBody2, uint64 inSortKey, Vec3Arg inWorldSpaceNormal, const ContactSettings &inSettings, uint32 inNumContactPoints)
{
	// Calculate the size of this constraint
	uint32 constraint_size = sizeof(ContactConstraint<Type1, Type2>) + (inNumContactPoints - 1) * sizeof(WorldContactPoint<Type1, Type2>);
	JPH_ASSERT(constraint_size <= cMaxConstraintSize);

	// Reserve space for constraint
	uint64 constraint_idx_and_constraint_offset = mNumConstraintsAndNextConstraintOffset.fetch_add((uint64(constraint_size) << 32) + 1, memory_order_relaxed);
	uint32 constraint_idx = uint32(constraint_idx_and_constraint_offset);
	if (constraint_idx >= mMaxConstraints)
		return nullptr;

	bool body1_dynamic = inBody1.IsDynamic();
	bool body2_dynamic = inBody2.IsDynamic();

	if (ioActivateAndLinkBodies)
	{
		// Do this only once
		ioActivateAndLinkBodies = false;

		// Wake up sleeping bodies
		BodyID body_ids[2];
		int num_bodies = 0;
		if (body1_dynamic && !inBody1.IsActive())
			body_ids[num_bodies++] = inBody1.GetID();
		if (body2_dynamic && !inBody2.IsActive())
			body_ids[num_bodies++] = inBody2.GetID();
		if (num_bodies > 0)
			mUpdateContext->mBodyManager->ActivateBodies(body_ids, num_bodies);

		// Link the two bodies only if both are dynamic. If one of them is static or kinematic they don't need to go into
		// the same simulation island as a constraint cannot affect the velocity of a kinematic body.
		if (body1_dynamic && body2_dynamic)
			mUpdateContext->mIslandBuilder->LinkBodies(inBody1.GetIndexInActiveBodiesInternal(), inBody2.GetIndexInActiveBodiesInternal());
	}

	// Link the contact to the first dynamic body
	if (body1_dynamic)
		mUpdateContext->mIslandBuilder->LinkContact(constraint_idx, inBody1.GetIndexInActiveBodiesInternal());
	else
	{
		JPH_ASSERT(body2_dynamic);
		mUpdateContext->mIslandBuilder->LinkContact(constraint_idx, inBody2.GetIndexInActiveBodiesInternal());
	}

	// Store offset for constraint
	uint32 constraint_offset = uint32(constraint_idx_and_constraint_offset >> 32);
	JPH_ASSERT(constraint_offset + constraint_size <= mMaxConstraints * cMaxConstraintSize);
	mConstraintIdxToOffset[constraint_idx] = constraint_offset;

	// Construct constraint
	ContactConstraint<Type1, Type2> *constraint = reinterpret_cast<ContactConstraint<Type1, Type2> *>(mConstraints + constraint_offset);
	new (constraint) ContactConstraint<Type1, Type2>();
	constraint->mBody1 = &inBody1;
	constraint->mBody2 = &inBody2;
	constraint->mSortKey = inSortKey;
	inWorldSpaceNormal.StoreFloat3(&constraint->mWorldSpaceNormal);
	constraint->mCombinedFriction = inSettings.mCombinedFriction;
	constraint->mCombinedRestitution = inSettings.mCombinedRestitution;
	constraint->mMinVelocityForRestitution = mPhysicsSettings.mMinVelocityForRestitution;
	constraint->mInvInertiaScale1 = inSettings.mInvInertiaScale1;
	constraint->mInvInertiaScale2 = inSettings.mInvInertiaScale2;
	constraint->mNumContactPoints = inNumContactPoints;
	constraint->mUseVPhysicsSolver = inSettings.mVPhysicsUseLegacySolver;
	constraint->mDisablePenetrationRecovery = inSettings.mVPhysicsDisablePenetrationRecovery;
	constraint->mSevereInterpenetration = inSettings.mVPhysicsSevereInterpenetration;
	constraint->mPlayerSupportStartup = inSettings.mVPhysicsPlayerSupportStartupContact;
	constraint->mPlayerSupportStartupPenetration =
		inSettings.mVPhysicsPlayerSupportStartupRescuePenetration >= 0.0f
			? inSettings.mVPhysicsPlayerSupportStartupRescuePenetration
			: inSettings.mVPhysicsPlayerSupportStartupPenetration;
	if (constraint->mPlayerSupportStartup)
	{
		constraint->mPlayerSupportStartupInitialNormalVelocity1 = inBody1.IsDynamic()
			? inBody1.GetLinearVelocity().Dot(inWorldSpaceNormal) : 0.0f;
		constraint->mPlayerSupportStartupInitialNormalVelocity2 = inBody2.IsDynamic()
			? inBody2.GetLinearVelocity().Dot(inWorldSpaceNormal) : 0.0f;
		constraint->mPlayerSupportStartupDeltaTime = mUpdateContext->mStepDeltaTime;
	}

#ifdef JPH_TRACK_SIMULATION_STATS
	// Track new contact constraints
	if constexpr (Type1 != EMotionType::Static)
		inBody1.GetMotionPropertiesUnchecked()->GetSimulationStats().mNumContactConstraints.fetch_add(1, memory_order_relaxed);
	if constexpr (Type2 != EMotionType::Static)
		inBody2.GetMotionPropertiesUnchecked()->GetSimulationStats().mNumContactConstraints.fetch_add(1, memory_order_relaxed);
#endif

	return constraint;
}

template <EMotionType Type1, EMotionType Type2>
void ContactConstraintManager::TemplatedGetContactsFromCache(ContactAllocator &ioContactAllocator, Body &inBody1, Body &inBody2, const CachedBodyPair &inCachedBodyPair, CachedBodyPair &outCachedBodyPair)
{
	// Get body transforms
	RMat44 transform_body1 = inBody1.GetCenterOfMassTransform();
	RMat44 transform_body2 = inBody2.GetCenterOfMassTransform();

	// Get time step and gravity
	float delta_time = mUpdateContext->mStepDeltaTime;
	Vec3 gravity = mUpdateContext->mPhysicsSystem->GetGravity();

	const ManifoldCache &read_cache = mCache[mCacheWriteIdx ^ 1];
	ManifoldCache &write_cache = mCache[mCacheWriteIdx];

	// Copy manifolds
	uint32 output_handle = ManifoldMap::cInvalidHandle;
	uint32 input_handle = inCachedBodyPair.mFirstCachedManifold;
	bool link_bodies = true;
	do
	{
		JPH_PROFILE("Add Constraint From Cached Manifold");

		// Find the existing manifold
		const MKeyValue *input_kv = read_cache.FromHandle(input_handle);
		const SubShapeIDPair &input_key = input_kv->GetKey();
		const CachedManifold &input_cm = input_kv->GetValue();
		JPH_ASSERT(input_cm.mNumContactPoints > 0); // There should be contact points in this manifold!

		// Create room for manifold in write buffer and copy data
		uint64 input_hash = input_key.GetHash();
		MKeyValue *output_kv = write_cache.Create(ioContactAllocator, input_key, input_hash, input_cm.mNumContactPoints);
		if (output_kv == nullptr)
			break; // Out of cache space
		CachedManifold *output_cm = &output_kv->GetValue();
		memcpy(output_cm, &input_cm, CachedManifold::sGetRequiredTotalSize(input_cm.mNumContactPoints));
		// This flag only bridges the discrete and CCD passes of one update. A manifold copied
		// from the previous update must become ordinary support again.
		output_cm->mFlags.fetch_and(
			(uint16)~(uint16)CachedManifold::EFlags::VPhysicsDeferPlayerSupport,
			memory_order_relaxed);
		for (uint32 point_index = 0; point_index < output_cm->mNumContactPoints; ++point_index)
			if (output_cm->mContactPoints[point_index].mVPhysicsContactEstablished)
			{
				output_cm->mFlags.fetch_and(
					(uint16)~(uint16)CachedManifold::EFlags::VPhysicsPlayerSupportStartup,
					memory_order_relaxed);
				break;
			}

		// Link the object under the body pairs
		output_cm->mNextWithSameBodyPair = output_handle;
		output_handle = write_cache.ToHandle(output_kv);

		// Calculate default contact settings
		ContactSettings settings;
		settings.mCombinedFriction = mCombineFriction(inBody1, input_key.GetSubShapeID1(), inBody2, input_key.GetSubShapeID2());
		settings.mCombinedRestitution = mCombineRestitution(inBody1, input_key.GetSubShapeID1(), inBody2, input_key.GetSubShapeID2());
		settings.mIsSensor = inBody1.IsSensor() || inBody2.IsSensor();
		settings.mVPhysicsPlayerSupportStartupContact =
			(output_cm->mFlags.load(memory_order_relaxed)
				& (uint16)CachedManifold::EFlags::VPhysicsPlayerSupportStartup) != 0;
		settings.mVPhysicsPlayerSupportStartupPenetration =
			output_cm->mVPhysicsPlayerSupportStartupPenetration;

		// Calculate world space contact normal
		Vec3 world_space_normal = transform_body2.Multiply3x3(Vec3::sLoadFloat3Unsafe(output_cm->mContactNormal)).Normalized();

		// Call contact listener to update settings
		if (mContactListener != nullptr)
		{
			// Convert constraint to manifold structure for callback
			ContactManifold manifold;
			manifold.mWorldSpaceNormal = world_space_normal;
			manifold.mSubShapeID1 = input_key.GetSubShapeID1();
			manifold.mSubShapeID2 = input_key.GetSubShapeID2();
			manifold.mBaseOffset = transform_body1.GetTranslation();
			manifold.mRelativeContactPointsOn1.resize(output_cm->mNumContactPoints);
			manifold.mRelativeContactPointsOn2.resize(output_cm->mNumContactPoints);
			Mat44 local_transform_body2 = transform_body2.PostTranslated(-manifold.mBaseOffset).ToMat44();
			float penetration_depth = -FLT_MAX;
			for (uint32 i = 0; i < output_cm->mNumContactPoints; ++i)
			{
				const CachedContactPoint &ccp = output_cm->mContactPoints[i];
				manifold.mRelativeContactPointsOn1[i] = transform_body1.Multiply3x3(Vec3::sLoadFloat3Unsafe(ccp.mPosition1));
				manifold.mRelativeContactPointsOn2[i] = local_transform_body2 * Vec3::sLoadFloat3Unsafe(ccp.mPosition2);
				penetration_depth = max(penetration_depth, (manifold.mRelativeContactPointsOn1[i] - manifold.mRelativeContactPointsOn2[i]).Dot(world_space_normal));
			}
			manifold.mPenetrationDepth = penetration_depth; // We don't have the penetration depth anymore, estimate it
			if (settings.mVPhysicsPlayerSupportStartupContact
				&& settings.mVPhysicsPlayerSupportStartupPenetration < 0.0f)
			{
				settings.mVPhysicsPlayerSupportStartupPenetration = max(0.0f, penetration_depth);
				output_cm->mVPhysicsPlayerSupportStartupPenetration =
					settings.mVPhysicsPlayerSupportStartupPenetration;
			}
			else if (settings.mVPhysicsPlayerSupportStartupContact
				&& IsPlayerSupportSpreadPenetrationEnabled())
			{
				const float live_penetration = max(0.0f, penetration_depth);
				const float remaining_feature_count = max(1.0f, float(output_cm->mNumContactPoints - 1));
				settings.mVPhysicsPlayerSupportStartupRescuePenetration =
					settings.mVPhysicsPlayerSupportStartupPenetration + max(0.0f,
					live_penetration - settings.mVPhysicsPlayerSupportStartupPenetration)
					/ remaining_feature_count;
			}

			// Notify callback
			mContactListener->OnContactPersisted(inBody1, inBody2, manifold, settings);
		}

		JPH_ASSERT(settings.mIsSensor || !(inBody1.IsSensor() || inBody2.IsSensor()), "Sensors cannot be converted into regular bodies by a contact callback!");
		if (!settings.mIsSensor // If one of the bodies is a sensor, don't actually create the constraint
			&& ((Type1 == EMotionType::Dynamic && settings.mInvMassScale1 != 0.0f) // One of the bodies must have mass to be able to create a contact constraint
				|| (Type2 == EMotionType::Dynamic && settings.mInvMassScale2 != 0.0f)))
		{
			// Create a new constraint
			ContactConstraint<Type1, Type2> *constraint = CreateConstraint<Type1, Type2>(link_bodies, inBody1, inBody2, input_hash, world_space_normal, settings, output_cm->mNumContactPoints);
			if (constraint == nullptr)
			{
				ioContactAllocator.mErrors |= EPhysicsUpdateError::ContactConstraintsFull;
				break;
			}
			constraint->mContactImpacting = false;
			constraint->mSubShapePair = input_key;

			JPH_DET_LOG("GetContactsFromCache: id1: " << inBody1.GetID() << " id2: " << inBody2.GetID() << " key: " << constraint->mSortKey);

			// Calculate scaled mass and inertia
			Mat44 inv_i1;
			if constexpr (Type1 == EMotionType::Dynamic)
			{
				const MotionProperties *mp1 = inBody1.GetMotionPropertiesUnchecked();
				constraint->mInvMass1 = settings.mInvMassScale1 * mp1->GetInverseMass();
				inv_i1 = settings.mInvInertiaScale1 * mp1->GetInverseInertiaForRotation(transform_body1.GetRotation());
			}
			else
			{
				constraint->mInvMass1 = 0.0f;
				inv_i1 = Mat44::sZero();
			}

			Mat44 inv_i2;
			if constexpr (Type2 == EMotionType::Dynamic)
			{
				const MotionProperties *mp2 = inBody2.GetMotionPropertiesUnchecked();
				constraint->mInvMass2 = settings.mInvMassScale2 * mp2->GetInverseMass();
				inv_i2 = settings.mInvInertiaScale2 * mp2->GetInverseInertiaForRotation(transform_body2.GetRotation());
			}
			else
			{
				constraint->mInvMass2 = 0.0f;
				inv_i2 = Mat44::sZero();
			}

			// Calculate tangents
			Vec3 t1, t2;
			constraint->GetTangents(t1, t2);

			// Setup velocity constraint properties
			float min_velocity_for_restitution = mPhysicsSettings.mMinVelocityForRestitution;
			for (uint32 i = 0; i < constraint->mNumContactPoints; ++i)
			{
				CachedContactPoint &ccp = output_cm->mContactPoints[i];
				WorldContactPoint<Type1, Type2> &wcp = constraint->mContactPoints[i];

				RVec3 p1_ws = transform_body1 * Vec3::sLoadFloat3Unsafe(ccp.mPosition1);
				RVec3 p2_ws = transform_body2 * Vec3::sLoadFloat3Unsafe(ccp.mPosition2);

				wcp.mNonPenetrationConstraint.SetTotalLambda(ccp.mNonPenetrationLambda);
				wcp.mFrictionConstraint1.SetTotalLambda(ccp.mFrictionLambda[0]);
				wcp.mFrictionConstraint2.SetTotalLambda(ccp.mFrictionLambda[1]);
				wcp.mContactPoint = &ccp;

				bool has_cached_support = wcp.mNonPenetrationConstraint.GetTotalLambda() > 0.0f
					|| ccp.mVPhysicsFrictionNormalLambda > 0.0f;
				const bool has_established_contact = has_cached_support
					|| ccp.mVPhysicsContactEstablished;
				if (IsPlayerSupportCacheTraceEnabled())
					std::fprintf(stderr,
						"PLAYER_SUPPORT_CACHE cached b1=%u b2=%u point=%u lambda=%.9g legacy=%.9g established=%d startup=%d impacting=%d\n",
						inBody1.GetID().GetIndexAndSequenceNumber(), inBody2.GetID().GetIndexAndSequenceNumber(),
						i, wcp.mNonPenetrationConstraint.GetTotalLambda(), ccp.mVPhysicsFrictionNormalLambda,
						int(ccp.mVPhysicsContactEstablished), int(settings.mVPhysicsPlayerSupportStartupContact),
						int(settings.mVPhysicsImpactingContact));
				constraint->mContactHasCachedSupport |= has_cached_support;
				bool allow_contact_rescue = !has_established_contact;
				const bool fresh_impact = settings.mVPhysicsSevereInterpenetration
					|| (settings.mVPhysicsImpactingContact
						&& !settings.mVPhysicsDeferPlayerSupportContact
						&& !has_established_contact);
				if (fresh_impact && settings.mVPhysicsUseLegacySolverForImpactOnly)
					constraint->mUseVPhysicsSolver = true;
				constraint->mContactImpacting |= fresh_impact;
				if (fresh_impact && !constraint->mUseVPhysicsSolver
					&& !constraint->mPlayerSupportStartup)
				{
					const Vec3 relative_velocity = inBody2.GetPointVelocity(p2_ws)
						- inBody1.GetPointVelocity(p1_ws);
					const float closing_speed = -relative_velocity.Dot(world_space_normal);
					const float gap = Vec3(p2_ws - p1_ws).Dot(world_space_normal);
					const float collision_tolerance = 0.5f * mPhysicsSettings.mContactSeparationTarget;
					float impact_time = 0.0f;
					if (gap > collision_tolerance && closing_speed > 1.0e-4f)
						impact_time = Clamp((gap - collision_tolerance) / closing_speed,
							0.0f, delta_time);
					if (!constraint->mVPhysicsNativeImpactCorrection
						|| impact_time < constraint->mVPhysicsNativeImpactTime)
						constraint->mVPhysicsNativeImpactTime = impact_time;
					constraint->mVPhysicsNativeImpactCorrection = true;
				}
				// A pending speculative manifold can persist for several updates before it actually
				// crosses IVP's inner impact distance. That first supported update still gets the one
				// restitution response; an already-supported rolling contact never does.
				bool allow_restitution = !mPhysicsSettings.mContactRestitutionOnlyOnNewContact || fresh_impact;
				wcp.CalculateFrictionAndNonPenetrationConstraintProperties(delta_time, gravity, inBody1, inBody2, constraint->mInvMass1, constraint->mInvMass2, inv_i1, inv_i2, p1_ws, p2_ws, world_space_normal, t1, t2, settings, min_velocity_for_restitution, mPhysicsSettings.mContactSeparationTarget, mPhysicsSettings.mContactVelocityBiasFactorOutsideTarget, mPhysicsSettings.mContactVelocityBiasFactorInsideTarget, mPhysicsSettings.mContactRescueVelocityFactor, allow_contact_rescue, allow_restitution);
			}

		#ifdef JPH_DEBUG_RENDERER
			// Draw the manifold
			if (sDrawContactManifolds)
				constraint->Draw(DebugRenderer::sInstance, Color::sYellow);
		#endif // JPH_DEBUG_RENDERER
		}

		// Mark contact as persisted so that we won't fire OnContactRemoved callbacks
		input_cm.mFlags |= (uint16)CachedManifold::EFlags::ContactPersisted;

		// Fetch the next manifold
		input_handle = input_cm.mNextWithSameBodyPair;
	}
	while (input_handle != ManifoldMap::cInvalidHandle);
	outCachedBodyPair.mFirstCachedManifold = output_handle;
}

void ContactConstraintManager::GetContactsFromCache(ContactAllocator &ioContactAllocator, Body &inBody1, Body &inBody2, bool &outPairHandled)
{
	// Start with not handled
	outPairHandled = false;

	// Swap bodies so that body 1 id < body 2 id
	Body *body1, *body2;
	if (inBody1.GetID() < inBody2.GetID())
	{
		body1 = &inBody1;
		body2 = &inBody2;
	}
	else
	{
		body1 = &inBody2;
		body2 = &inBody1;
	}

	// Find the cached body pair
	BodyPair body_pair_key(body1->GetID(), body2->GetID());
	uint64 body_pair_hash = body_pair_key.GetHash();
	const ManifoldCache &read_cache = mCache[mCacheWriteIdx ^ 1];
	const BPKeyValue *kv = read_cache.Find(body_pair_key, body_pair_hash);
	if (kv == nullptr)
		return;
	const CachedBodyPair &input_cbp = kv->GetValue();

	// Get relative translation
	Quat inv_r1 = body1->GetRotation().Conjugated();
	Vec3 delta_position = inv_r1 * Vec3(body2->GetCenterOfMassPosition() - body1->GetCenterOfMassPosition());

	// Get old position delta
	Vec3 old_delta_position = Vec3::sLoadFloat3Unsafe(input_cbp.mDeltaPosition);

	// Check if bodies are still roughly in the same relative position
	if ((delta_position - old_delta_position).LengthSq() > mPhysicsSettings.mBodyPairCacheMaxDeltaPositionSq)
		return;

	// Determine relative orientation
	Quat delta_rotation = inv_r1 * body2->GetRotation();

	// Reconstruct old quaternion delta
	Quat old_delta_rotation = Quat::sLoadFloat3Unsafe(input_cbp.mDeltaRotation);

	// Check if bodies are still roughly in the same relative orientation
	// The delta between 2 quaternions p and q is: p q^* = [rotation_axis * sin(angle / 2), cos(angle / 2)]
	// From the W component we can extract the angle: cos(angle / 2) = px * qx + py * qy + pz * qz + pw * qw = p . q
	// Since we want to abort if the rotation is smaller than -angle or bigger than angle, we can write the comparison as |p . q| < cos(angle / 2)
	if (abs(delta_rotation.Dot(old_delta_rotation)) < mPhysicsSettings.mBodyPairCacheCosMaxDeltaRotationDiv2)
		return;

	// The cache is valid, return that we've handled this body pair
	outPairHandled = true;

	// Copy the cached body pair to this frame
	ManifoldCache &write_cache = mCache[mCacheWriteIdx];
	BPKeyValue *output_bp_kv = write_cache.Create(ioContactAllocator, body_pair_key, body_pair_hash);
	if (output_bp_kv == nullptr)
		return; // Out of cache space
	CachedBodyPair *output_cbp = &output_bp_kv->GetValue();
	memcpy(output_cbp, &input_cbp, sizeof(CachedBodyPair));

	// If there were no contacts, we have handled the contact
	if (input_cbp.mFirstCachedManifold == ManifoldMap::cInvalidHandle)
		return;

	// Build dispatch table
	// Note: Non-dynamic vs non-dynamic can happen in this case due to one body being a sensor, so we need to have an extended table here
	using DispatchFunc = void (ContactConstraintManager::*)(ContactAllocator &, Body &, Body &, const CachedBodyPair &, CachedBodyPair &);
	static const DispatchFunc table[3][3] = {
		{
			nullptr, // Static vs static doesn't exist
			&ContactConstraintManager::TemplatedGetContactsFromCache<EMotionType::Static, EMotionType::Kinematic>,
			&ContactConstraintManager::TemplatedGetContactsFromCache<EMotionType::Static, EMotionType::Dynamic>
		},
		{
			&ContactConstraintManager::TemplatedGetContactsFromCache<EMotionType::Kinematic, EMotionType::Static>,
			&ContactConstraintManager::TemplatedGetContactsFromCache<EMotionType::Kinematic, EMotionType::Kinematic>,
			&ContactConstraintManager::TemplatedGetContactsFromCache<EMotionType::Kinematic, EMotionType::Dynamic>
		},
		{
			&ContactConstraintManager::TemplatedGetContactsFromCache<EMotionType::Dynamic, EMotionType::Static>,
			&ContactConstraintManager::TemplatedGetContactsFromCache<EMotionType::Dynamic, EMotionType::Kinematic>,
			&ContactConstraintManager::TemplatedGetContactsFromCache<EMotionType::Dynamic, EMotionType::Dynamic>
		}
	};

	// Dispatch to the correct templated form
	(this->*table[(int)body1->GetMotionType()][(int)body2->GetMotionType()])(ioContactAllocator, *body1, *body2, input_cbp, *output_cbp);
}

ContactConstraintManager::BodyPairHandle ContactConstraintManager::AddBodyPair(ContactAllocator &ioContactAllocator, const Body &inBody1, const Body &inBody2)
{
	// Swap bodies so that body 1 id < body 2 id
	const Body *body1, *body2;
	if (inBody1.GetID() < inBody2.GetID())
	{
		body1 = &inBody1;
		body2 = &inBody2;
	}
	else
	{
		body1 = &inBody2;
		body2 = &inBody1;
	}

	// Add an entry
	BodyPair body_pair_key(body1->GetID(), body2->GetID());
	uint64 body_pair_hash = body_pair_key.GetHash();
	BPKeyValue *body_pair_kv = mCache[mCacheWriteIdx].Create(ioContactAllocator, body_pair_key, body_pair_hash);
	if (body_pair_kv == nullptr)
		return nullptr; // Out of cache space
	CachedBodyPair *cbp = &body_pair_kv->GetValue();
	cbp->mFirstCachedManifold = ManifoldMap::cInvalidHandle;

	// Get relative translation
	Quat inv_r1 = body1->GetRotation().Conjugated();
	Vec3 delta_position = inv_r1 * Vec3(body2->GetCenterOfMassPosition() - body1->GetCenterOfMassPosition());

	// Store it
	delta_position.StoreFloat3(&cbp->mDeltaPosition);

	// Determine relative orientation
	Quat delta_rotation = inv_r1 * body2->GetRotation();

	// Store it
	delta_rotation.StoreFloat3(&cbp->mDeltaRotation);

	return cbp;
}

template <EMotionType Type1, EMotionType Type2>
void ContactConstraintManager::TemplatedAddContactConstraint(ContactAllocator &ioContactAllocator, bool &ioActivateAndLinkBodies, BodyPairHandle inBodyPairHandle, Body &inBody1, Body &inBody2, const ContactManifold &inManifold)
{
	// Calculate hash
	SubShapeIDPair key { inBody1.GetID(), inManifold.mSubShapeID1, inBody2.GetID(), inManifold.mSubShapeID2 };
	uint64 key_hash = key.GetHash();

	// Determine number of contact points
	int num_contact_points = (int)inManifold.mRelativeContactPointsOn1.size();
	JPH_ASSERT(num_contact_points <= MaxContactPoints);
	JPH_ASSERT(num_contact_points == (int)inManifold.mRelativeContactPointsOn2.size());

	// Reserve space for new contact cache entry
	// Note that for dynamic vs dynamic we always require the first body to have a lower body id to get a consistent key
	// under which to look up the contact
	ManifoldCache &write_cache = mCache[mCacheWriteIdx];
	MKeyValue *new_manifold_kv = write_cache.Create(ioContactAllocator, key, key_hash, num_contact_points);
	if (new_manifold_kv == nullptr)
		return; // Out of cache space
	CachedManifold *new_manifold = &new_manifold_kv->GetValue();
	new_manifold->mVPhysicsPlayerSupportStartupPenetration = -1.0f;

	// Transform the world space normal to the space of body 2 (this is usually the static body)
	RMat44 inverse_transform_body2 = inBody2.GetInverseCenterOfMassTransform();
	inverse_transform_body2.Multiply3x3(inManifold.mWorldSpaceNormal).Normalized().StoreFloat3(&new_manifold->mContactNormal);

	// Settings object that gets passed to the callback
	ContactSettings settings;
	settings.mCombinedFriction = mCombineFriction(inBody1, inManifold.mSubShapeID1, inBody2, inManifold.mSubShapeID2);
	settings.mCombinedRestitution = mCombineRestitution(inBody1, inManifold.mSubShapeID1, inBody2, inManifold.mSubShapeID2);
	settings.mIsSensor = inBody1.IsSensor() || inBody2.IsSensor();

	// Get the contact points for the old cache entry
	const ManifoldCache &read_cache = mCache[mCacheWriteIdx ^ 1];
	const MKeyValue *old_manifold_kv = read_cache.Find(key, key_hash);
	const bool contact_persisted = old_manifold_kv != nullptr;
	const CachedManifold *old_manifold = nullptr;
	const CachedContactPoint *ccp_start;
	const CachedContactPoint *ccp_end;
	if (old_manifold_kv != nullptr)
	{
		const CachedManifold &old_cached_manifold = old_manifold_kv->GetValue();
		settings.mVPhysicsPlayerSupportStartupContact =
			(old_cached_manifold.mFlags.load(memory_order_relaxed)
				& (uint16)CachedManifold::EFlags::VPhysicsPlayerSupportStartup) != 0;
		settings.mVPhysicsPlayerSupportStartupPenetration =
			old_cached_manifold.mVPhysicsPlayerSupportStartupPenetration;
		for (uint32 point_index = 0;
			settings.mVPhysicsPlayerSupportStartupContact
				&& point_index < old_cached_manifold.mNumContactPoints;
			++point_index)
		{
			if (old_cached_manifold.mContactPoints[point_index].mVPhysicsContactEstablished)
				settings.mVPhysicsPlayerSupportStartupContact = false;
		}
		if (settings.mVPhysicsPlayerSupportStartupContact
			&& settings.mVPhysicsPlayerSupportStartupPenetration < 0.0f)
		{
			settings.mVPhysicsPlayerSupportStartupPenetration =
				max(0.0f, inManifold.mPenetrationDepth);
		}
		else if (settings.mVPhysicsPlayerSupportStartupContact
			&& IsPlayerSupportSpreadPenetrationEnabled())
		{
			const float live_penetration = max(0.0f, inManifold.mPenetrationDepth);
			const float remaining_feature_count = max(1.0f, float(inManifold.mRelativeContactPointsOn1.size() - 1));
			settings.mVPhysicsPlayerSupportStartupRescuePenetration =
				settings.mVPhysicsPlayerSupportStartupPenetration + max(0.0f,
				live_penetration - settings.mVPhysicsPlayerSupportStartupPenetration)
				/ remaining_feature_count;
		}
		// Call point persisted listener
		if (mContactListener != nullptr)
			mContactListener->OnContactPersisted(inBody1, inBody2, inManifold, settings);

		// Fetch the contact points from the old manifold
		old_manifold = &old_manifold_kv->GetValue();
		ccp_start = old_manifold->mContactPoints;
		ccp_end = ccp_start + old_manifold->mNumContactPoints;

		// Mark contact as persisted so that we won't fire OnContactRemoved callbacks
		old_manifold->mFlags |= (uint16)CachedManifold::EFlags::ContactPersisted;
	}
	else
	{
		// Call point added listener
		if (mContactListener != nullptr)
			mContactListener->OnContactAdded(inBody1, inBody2, inManifold, settings);

		// No contact points available from old manifold
		ccp_start = nullptr;
		ccp_end = nullptr;
	}
	if (settings.mVPhysicsDeferPlayerSupportContact)
		new_manifold->mFlags |= (uint16)CachedManifold::EFlags::VPhysicsDeferPlayerSupport;
	if (settings.mVPhysicsPlayerSupportStartupContact)
		new_manifold->mFlags |= (uint16)CachedManifold::EFlags::VPhysicsPlayerSupportStartup;
	new_manifold->mVPhysicsPlayerSupportStartupPenetration =
		settings.mVPhysicsPlayerSupportStartupPenetration;

	// Get inverse transform for body 1
	RMat44 inverse_transform_body1 = inBody1.GetInverseCenterOfMassTransform();

	// If one of the bodies is a sensor, don't actually create the constraint
	JPH_ASSERT(settings.mIsSensor || !(inBody1.IsSensor() || inBody2.IsSensor()), "Sensors cannot be converted into regular bodies by a contact callback!");
	if (!settings.mIsSensor
		&& ((Type1 == EMotionType::Dynamic && settings.mInvMassScale1 != 0.0f) // One of the bodies must have mass to be able to create a contact constraint
			|| (Type2 == EMotionType::Dynamic && settings.mInvMassScale2 != 0.0f)))
	{
		// Create a new constraint
		ContactConstraint<Type1, Type2> *constraint = CreateConstraint<Type1, Type2>(ioActivateAndLinkBodies, inBody1, inBody2, key_hash, inManifold.mWorldSpaceNormal, settings, num_contact_points);
		if (constraint == nullptr)
		{
			ioContactAllocator.mErrors |= EPhysicsUpdateError::ContactConstraintsFull;

			// Manifold has been created already, we're not filling it in, so we need to reset the contact number of points.
			// Note that we don't hook it up to the body pair cache so that it won't be used as a cache during the next simulation.
			new_manifold->mNumContactPoints = 0;
			return;
		}
		constraint->mContactImpacting = false;
		constraint->mSubShapePair = key;

		JPH_DET_LOG("AddContactConstraint: id1: " << constraint->mBody1->GetID() << " id2: " << constraint->mBody2->GetID() << " key: " << constraint->mSortKey);

		// Get time step and gravity
		float delta_time = mUpdateContext->mStepDeltaTime;
		Vec3 gravity = mUpdateContext->mPhysicsSystem->GetGravity();

		// Calculate scaled mass and inertia
		Mat44 inv_i1;
		if constexpr (Type1 == EMotionType::Dynamic)
		{
			const MotionProperties *mp1 = inBody1.GetMotionPropertiesUnchecked();
			constraint->mInvMass1 = settings.mInvMassScale1 * mp1->GetInverseMass();
			inv_i1 = settings.mInvInertiaScale1 * mp1->GetInverseInertiaForRotation(inverse_transform_body1.Transposed3x3());
		}
		else
		{
			constraint->mInvMass1 = 0.0f;
			inv_i1 = Mat44::sZero();
		}

		Mat44 inv_i2;
		if constexpr (Type2 == EMotionType::Dynamic)
		{
			const MotionProperties *mp2 = inBody2.GetMotionPropertiesUnchecked();
			constraint->mInvMass2 = settings.mInvMassScale2 * mp2->GetInverseMass();
			inv_i2 = settings.mInvInertiaScale2 * mp2->GetInverseInertiaForRotation(inverse_transform_body2.Transposed3x3());
		}
		else
		{
			constraint->mInvMass2 = 0.0f;
			inv_i2 = Mat44::sZero();
		}

		// Calculate tangents
		Vec3 t1, t2;
		constraint->GetTangents(t1, t2);

		// IVP transfers a friction contact's pretension through world space when the
		// underlying point/edge/triangle contact is replaced. Jolt's normal/friction
		// warm-start cache intentionally requires nearly identical local points, which is
		// too strict for a rolling wheel. Keep that world-space spring independently of
		// the regular Jolt lambda cache. In particular, do not reconstruct it through
		// mContactNormal: that normal is stored in body 2 space, which is often the rolling
		// wheel and would incorrectly rotate the ground friction spring with the wheel.
		bool legacy_contact_used[MaxContactPoints] = { false };

		for (int i = 0; i < num_contact_points; ++i)
		{
			// Convert to world space and set positions
			WorldContactPoint<Type1, Type2> &wcp = constraint->mContactPoints[i];
			RVec3 p1_ws = inManifold.mBaseOffset + inManifold.mRelativeContactPointsOn1[i];
			RVec3 p2_ws = inManifold.mBaseOffset + inManifold.mRelativeContactPointsOn2[i];

			// Convert to local space to the body
			Vec3 p1_ls = Vec3(inverse_transform_body1 * p1_ws);
			Vec3 p2_ls = Vec3(inverse_transform_body2 * p2_ws);

			// Check if we have a close contact point from last update
			bool lambda_set = false;
			int matching_contact_index = -1;
			Vec3 legacy_friction_displacement = Vec3::sZero();
			float legacy_friction_normal_lambda = 0.0f;
			bool legacy_contact_established = false;
			for (const CachedContactPoint *ccp = ccp_start; ccp < ccp_end; ccp++)
				if (Vec3::sLoadFloat3Unsafe(ccp->mPosition1).IsClose(p1_ls, mPhysicsSettings.mContactPointPreserveLambdaMaxDistSq)
					&& Vec3::sLoadFloat3Unsafe(ccp->mPosition2).IsClose(p2_ls, mPhysicsSettings.mContactPointPreserveLambdaMaxDistSq))
				{
					// Get lambdas from previous frame
					wcp.mNonPenetrationConstraint.SetTotalLambda(ccp->mNonPenetrationLambda);
					wcp.mFrictionConstraint1.SetTotalLambda(ccp->mFrictionLambda[0]);
					wcp.mFrictionConstraint2.SetTotalLambda(ccp->mFrictionLambda[1]);
					matching_contact_index = int(ccp - ccp_start);
					lambda_set = true;
					break;
				}
			if (!lambda_set)
			{
				wcp.mNonPenetrationConstraint.SetTotalLambda(0.0f);
				wcp.mFrictionConstraint1.SetTotalLambda(0.0f);
				wcp.mFrictionConstraint2.SetTotalLambda(0.0f);
			}

			// The regular warm-start match may legitimately fail as a sphere rolls. IVP does
			// not discard the friction spring in that case: it transfers pressure and the
			// world-space pretension to the replacement contact. Match each new point to the
			// nearest unused old point inside the already-persisted sub-shape manifold.
			int legacy_contact_index = matching_contact_index;
			if (legacy_contact_index < 0 || legacy_contact_used[legacy_contact_index])
			{
				float closest_distance_sq = FLT_MAX;
				for (const CachedContactPoint *ccp = ccp_start; ccp < ccp_end; ++ccp)
				{
					const int candidate_index = int(ccp - ccp_start);
					if (legacy_contact_used[candidate_index])
						continue;
					const float distance_sq = (Vec3::sLoadFloat3Unsafe(ccp->mPosition1) - p1_ls).LengthSq()
						+ (Vec3::sLoadFloat3Unsafe(ccp->mPosition2) - p2_ls).LengthSq();
					if (distance_sq < closest_distance_sq)
					{
						closest_distance_sq = distance_sq;
						legacy_contact_index = candidate_index;
					}
				}
			}
			if (legacy_contact_index >= 0)
			{
				legacy_contact_used[legacy_contact_index] = true;
				const CachedContactPoint &legacy_contact = ccp_start[legacy_contact_index];
				legacy_friction_displacement = Vec3::sLoadFloat3Unsafe(
					legacy_contact.mVPhysicsFrictionDisplacement);
				legacy_friction_normal_lambda = legacy_contact.mVPhysicsFrictionNormalLambda;
				legacy_contact_established = legacy_contact.mVPhysicsContactEstablished;
			}

			// Create new contact point
			CachedContactPoint &cp = new_manifold->mContactPoints[i];
			p1_ls.StoreFloat3(&cp.mPosition1);
			p2_ls.StoreFloat3(&cp.mPosition2);
			legacy_friction_displacement.StoreFloat3(&cp.mVPhysicsFrictionDisplacement);
			cp.mVPhysicsFrictionNormalLambda = legacy_friction_normal_lambda;
			cp.mVPhysicsContactEstablished = legacy_contact_established;
			wcp.mContactPoint = &cp;

			// Setup velocity constraint
			bool has_cached_support = wcp.mNonPenetrationConstraint.GetTotalLambda() > 0.0f
				|| legacy_friction_normal_lambda > 0.0f;
			const bool has_established_contact = has_cached_support
				|| legacy_contact_established;
			if (IsPlayerSupportCacheTraceEnabled())
				std::fprintf(stderr,
					"PLAYER_SUPPORT_CACHE added b1=%u b2=%u point=%d lambda=%.9g legacy=%.9g established=%d startup=%d impacting=%d\n",
					inBody1.GetID().GetIndexAndSequenceNumber(), inBody2.GetID().GetIndexAndSequenceNumber(),
					i, wcp.mNonPenetrationConstraint.GetTotalLambda(), legacy_friction_normal_lambda,
					int(legacy_contact_established), int(settings.mVPhysicsPlayerSupportStartupContact),
					int(settings.mVPhysicsImpactingContact));
			constraint->mContactHasCachedSupport |= has_cached_support;
			bool allow_contact_rescue = !has_established_contact;
			const bool fresh_impact = settings.mVPhysicsSevereInterpenetration
				|| (settings.mVPhysicsImpactingContact
					&& !settings.mVPhysicsDeferPlayerSupportContact
					&& !has_established_contact);
			if (fresh_impact && settings.mVPhysicsUseLegacySolverForImpactOnly)
				constraint->mUseVPhysicsSolver = true;
			constraint->mContactImpacting |= fresh_impact;
			if (fresh_impact && !constraint->mUseVPhysicsSolver
				&& !constraint->mPlayerSupportStartup)
			{
				const Vec3 relative_velocity = inBody2.GetPointVelocity(p2_ws)
					- inBody1.GetPointVelocity(p1_ws);
				const float closing_speed = -relative_velocity.Dot(inManifold.mWorldSpaceNormal);
				const float gap = Vec3(p2_ws - p1_ws).Dot(inManifold.mWorldSpaceNormal);
				const float collision_tolerance = 0.5f * mPhysicsSettings.mContactSeparationTarget;
				float impact_time = 0.0f;
				if (gap > collision_tolerance && closing_speed > 1.0e-4f)
					impact_time = Clamp((gap - collision_tolerance) / closing_speed,
						0.0f, delta_time);
				if (!constraint->mVPhysicsNativeImpactCorrection
					|| impact_time < constraint->mVPhysicsNativeImpactTime)
					constraint->mVPhysicsNativeImpactTime = impact_time;
				constraint->mVPhysicsNativeImpactCorrection = true;
			}
			bool allow_restitution = !mPhysicsSettings.mContactRestitutionOnlyOnNewContact
				|| !contact_persisted || fresh_impact;
			wcp.CalculateFrictionAndNonPenetrationConstraintProperties(delta_time, gravity, inBody1, inBody2, constraint->mInvMass1, constraint->mInvMass2, inv_i1, inv_i2, p1_ws, p2_ws, inManifold.mWorldSpaceNormal, t1, t2, settings, mPhysicsSettings.mMinVelocityForRestitution, mPhysicsSettings.mContactSeparationTarget, mPhysicsSettings.mContactVelocityBiasFactorOutsideTarget, mPhysicsSettings.mContactVelocityBiasFactorInsideTarget, mPhysicsSettings.mContactRescueVelocityFactor, allow_contact_rescue, allow_restitution);
		}

	#ifdef JPH_DEBUG_RENDERER
		// Draw the manifold
		if (sDrawContactManifolds)
			constraint->Draw(DebugRenderer::sInstance, Color::sOrange);
	#endif // JPH_DEBUG_RENDERER
	}
	else
	{
		// Store the contact manifold in the cache
		for (int i = 0; i < num_contact_points; ++i)
		{
			// Convert to local space to the body
			Vec3 p1 = Vec3(inverse_transform_body1 * (inManifold.mBaseOffset + inManifold.mRelativeContactPointsOn1[i]));
			Vec3 p2 = Vec3(inverse_transform_body2 * (inManifold.mBaseOffset + inManifold.mRelativeContactPointsOn2[i]));

			// Create new contact point
			CachedContactPoint &cp = new_manifold->mContactPoints[i];
			p1.StoreFloat3(&cp.mPosition1);
			p2.StoreFloat3(&cp.mPosition2);

			// Reset contact impulses, we haven't applied any
			cp.mNonPenetrationLambda = 0.0f;
			cp.mFrictionLambda[0] = 0.0f;
			cp.mFrictionLambda[1] = 0.0f;
			Vec3::sZero().StoreFloat3(&cp.mVPhysicsFrictionDisplacement);
			cp.mVPhysicsFrictionNormalLambda = 0.0f;
			cp.mVPhysicsContactEstablished = false;
		}
	}

	// Store cached contact point in body pair cache
	CachedBodyPair *cbp = reinterpret_cast<CachedBodyPair *>(inBodyPairHandle);
	new_manifold->mNextWithSameBodyPair = cbp->mFirstCachedManifold;
	cbp->mFirstCachedManifold = write_cache.ToHandle(new_manifold_kv);
}

void ContactConstraintManager::AddContactConstraint(ContactAllocator &ioContactAllocator, bool &ioActivateAndLinkBodies, BodyPairHandle inBodyPairHandle, Body &inBody1, Body &inBody2, const ContactManifold &inManifold)
{
	JPH_PROFILE_FUNCTION();

	JPH_DET_LOG("AddContactConstraint: id1: " << inBody1.GetID() << " id2: " << inBody2.GetID()
		<< " subshape1: " << inManifold.mSubShapeID1 << " subshape2: " << inManifold.mSubShapeID2
		<< " normal: " << inManifold.mWorldSpaceNormal << " pendepth: " << inManifold.mPenetrationDepth);

	JPH_ASSERT(inManifold.mWorldSpaceNormal.IsNormalized());

	// Swap bodies so that body 1 id < body 2 id
	const ContactManifold *manifold;
	Body *body1, *body2;
	ContactManifold temp;
	if (inBody2.GetID() < inBody1.GetID())
	{
		body1 = &inBody2;
		body2 = &inBody1;
		temp = inManifold.SwapShapes();
		manifold = &temp;
	}
	else
	{
		body1 = &inBody1;
		body2 = &inBody2;
		manifold = &inManifold;
	}

	// Build dispatch table
	// Note: Non-dynamic vs non-dynamic can happen in this case due to one body being a sensor, so we need to have an extended table here
	using DispatchFunc = void (ContactConstraintManager::*)(ContactAllocator &, bool &, BodyPairHandle, Body &, Body &, const ContactManifold &);
	static const DispatchFunc table[3][3] = {
		{
			nullptr, // Static vs static doesn't exist
			&ContactConstraintManager::TemplatedAddContactConstraint<EMotionType::Static, EMotionType::Kinematic>,
			&ContactConstraintManager::TemplatedAddContactConstraint<EMotionType::Static, EMotionType::Dynamic>
		},
		{
			&ContactConstraintManager::TemplatedAddContactConstraint<EMotionType::Kinematic, EMotionType::Static>,
			&ContactConstraintManager::TemplatedAddContactConstraint<EMotionType::Kinematic, EMotionType::Kinematic>,
			&ContactConstraintManager::TemplatedAddContactConstraint<EMotionType::Kinematic, EMotionType::Dynamic>
		},
		{
			&ContactConstraintManager::TemplatedAddContactConstraint<EMotionType::Dynamic, EMotionType::Static>,
			&ContactConstraintManager::TemplatedAddContactConstraint<EMotionType::Dynamic, EMotionType::Kinematic>,
			&ContactConstraintManager::TemplatedAddContactConstraint<EMotionType::Dynamic, EMotionType::Dynamic>
		}
	};

	// Dispatch to the correct templated form
	return (this->*table[(int)body1->GetMotionType()][(int)body2->GetMotionType()])(ioContactAllocator, ioActivateAndLinkBodies, inBodyPairHandle, *body1, *body2, *manifold);
}

void ContactConstraintManager::OnCCDContactAdded(ContactAllocator &ioContactAllocator, const Body &inBody1, const Body &inBody2, const ContactManifold &inManifold, ContactSettings &outSettings)
{
	JPH_ASSERT(inManifold.mWorldSpaceNormal.IsNormalized());

	// Calculate contact settings
	outSettings.mCombinedFriction = mCombineFriction(inBody1, inManifold.mSubShapeID1, inBody2, inManifold.mSubShapeID2);
	outSettings.mCombinedRestitution = mCombineRestitution(inBody1, inManifold.mSubShapeID1, inBody2, inManifold.mSubShapeID2);
	outSettings.mIsSensor = false; // For now, no sensors are supported during CCD
	outSettings.mVPhysicsCCDContact = true;
	outSettings.mVPhysicsDeferPlayerSupportContact = false;
	outSettings.mVPhysicsPlayerSupportStartupContact = false;

	// The remainder of this function only deals with calling contact callbacks, if there's no contact callback we also don't need to do this work
	if (mContactListener != nullptr)
	{
		// Swap bodies so that body 1 id < body 2 id
		const ContactManifold *manifold;
		const Body *body1, *body2;
		ContactManifold temp;
		if (inBody2.GetID() < inBody1.GetID())
		{
			body1 = &inBody2;
			body2 = &inBody1;
			temp = inManifold.SwapShapes();
			manifold = &temp;
		}
		else
		{
			body1 = &inBody1;
			body2 = &inBody2;
			manifold = &inManifold;
		}

		// Calculate hash
		SubShapeIDPair key { body1->GetID(), manifold->mSubShapeID1, body2->GetID(), manifold->mSubShapeID2 };
		uint64 key_hash = key.GetHash();

		// Check if we already created this contact this physics update
		ManifoldCache &write_cache = mCache[mCacheWriteIdx];
		MKVAndCreated new_manifold_kv = write_cache.FindOrCreate(ioContactAllocator, key, key_hash, 0);
		if (new_manifold_kv.first != nullptr)
		{
			const uint16 vphysics_flags =
				new_manifold_kv.first->GetValue().mFlags.load(memory_order_relaxed);
			if ((vphysics_flags & (uint16)CachedManifold::EFlags::VPhysicsDeferPlayerSupport) != 0)
			{
				// LinearCast can report the same pair again after the discrete pass. Carry the
				// one-PSI Source player-support deferral into CCD so it cannot reintroduce the
				// response that the discrete manifold deliberately suppressed.
				outSettings.mVPhysicsDeferPlayerSupportContact = true;
			}
			outSettings.mVPhysicsPlayerSupportStartupContact =
				(vphysics_flags & (uint16)CachedManifold::EFlags::VPhysicsPlayerSupportStartup) != 0;
			outSettings.mVPhysicsPlayerSupportStartupPenetration =
				new_manifold_kv.first->GetValue().mVPhysicsPlayerSupportStartupPenetration;
		}
		if (new_manifold_kv.second)
		{
			// This contact is new for this physics update, check if previous update we already had this contact.
			const ManifoldCache &read_cache = mCache[mCacheWriteIdx ^ 1];
			const MKeyValue *old_manifold_kv = read_cache.Find(key, key_hash);
			if (old_manifold_kv == nullptr)
			{
				// New contact
				mContactListener->OnContactAdded(*body1, *body2, *manifold, outSettings);
			}
			else
			{
				// Existing contact
				mContactListener->OnContactPersisted(*body1, *body2, *manifold, outSettings);

				// Mark contact as persisted so that we won't fire OnContactRemoved callbacks
				old_manifold_kv->GetValue().mFlags |= (uint16)CachedManifold::EFlags::ContactPersisted;
			}

			// Check if the cache is full
			if (new_manifold_kv.first != nullptr)
			{
				// We don't store any contact points in this manifold as it is not for caching impulses, we only need to know that the contact was created
				CachedManifold &new_manifold = new_manifold_kv.first->GetValue();
				new_manifold.mContactNormal = { 0, 0, 0 };
				new_manifold.mFlags |= (uint16)CachedManifold::EFlags::CCDContact;
				if (outSettings.mVPhysicsDeferPlayerSupportContact)
					new_manifold.mFlags |= (uint16)CachedManifold::EFlags::VPhysicsDeferPlayerSupport;
				if (outSettings.mVPhysicsPlayerSupportStartupContact)
					new_manifold.mFlags |= (uint16)CachedManifold::EFlags::VPhysicsPlayerSupportStartup;
				new_manifold.mVPhysicsPlayerSupportStartupPenetration =
					outSettings.mVPhysicsPlayerSupportStartupPenetration;
			}
		}
		else
		{
			// Already found this contact this physics update.
			// Note that we can trigger OnContactPersisted multiple times per physics update, but otherwise we have no way of obtaining the settings
			mContactListener->OnContactPersisted(*body1, *body2, *manifold, outSettings);
		}

		// If we swapped body1 and body2 we need to swap the mass scales back
		if (manifold == &temp)
		{
			std::swap(outSettings.mInvMassScale1, outSettings.mInvMassScale2);
			std::swap(outSettings.mInvInertiaScale1, outSettings.mInvInertiaScale2);
			// Note we do not need to negate the relative surface velocity as it is not applied by the CCD collision constraint
		}
	}

	JPH_ASSERT(outSettings.mIsSensor || !(inBody1.IsSensor() || inBody2.IsSensor()), "Sensors cannot be converted into regular bodies by a contact callback!");
}

void ContactConstraintManager::ConstraintIdxToConstraintOffset(uint32 *ioConstraintIdxBegin, const uint32 *inConstraintIdxEnd) const
{
	JPH_PROFILE_FUNCTION();

	for (uint32 *i = ioConstraintIdxBegin; i < inConstraintIdxEnd; ++i)
		*i = mConstraintIdxToOffset[*i];
}

void ContactConstraintManager::SortContacts(uint32 *ioConstraintOffsetBegin, uint32 *inConstraintOffsetEnd) const
{
	JPH_PROFILE_FUNCTION();

	QuickSort(ioConstraintOffsetBegin, inConstraintOffsetEnd, [this](uint32 inLHS, uint32 inRHS) {
		const ContactConstraintBase &lhs = *reinterpret_cast<const ContactConstraintBase *>(mConstraints + inLHS);
		const ContactConstraintBase &rhs = *reinterpret_cast<const ContactConstraintBase *>(mConstraints + inRHS);

		// Most of the time the sort key will be different so we sort on that
		if (lhs.mSortKey != rhs.mSortKey)
			return lhs.mSortKey < rhs.mSortKey;

		// If they're equal we use the IDs of body 1 to order
		if (lhs.mBody1 != rhs.mBody1)
			return lhs.mBody1->GetID() < rhs.mBody1->GetID();

		// If they're still equal we use the IDs of body 2 to order
		if (lhs.mBody2 != rhs.mBody2)
			return lhs.mBody2->GetID() < rhs.mBody2->GetID();

		JPH_ASSERT(inLHS == inRHS, "Hash collision, ordering will be inconsistent");
		return false;
	});
}

void ContactConstraintManager::FinalizeContactCacheAndCallContactPointRemovedCallbacks(uint inExpectedNumBodyPairs, uint inExpectedNumManifolds)
{
	JPH_PROFILE_FUNCTION();

#ifdef JPH_ENABLE_ASSERTS
	// Mark cache as finalized
	ManifoldCache &old_write_cache = mCache[mCacheWriteIdx];
	old_write_cache.Finalize();

	// Check that the count of body pairs and manifolds that we tracked outside of the cache (to avoid contention on an atomic) is correct
	JPH_ASSERT(old_write_cache.GetNumBodyPairs() == inExpectedNumBodyPairs);
	JPH_ASSERT(old_write_cache.GetNumManifolds() == inExpectedNumManifolds);
#endif

	// Buffers are now complete, make write buffer the read buffer
	mCacheWriteIdx ^= 1;

	// Get the old read cache / new write cache
	ManifoldCache &old_read_cache = mCache[mCacheWriteIdx];

	// Call the contact point removal callbacks
	if (mContactListener != nullptr)
		old_read_cache.ContactPointRemovedCallbacks(mContactListener);

	// We're done with the old read cache now
	old_read_cache.Clear();

	// Use the amount of contacts from the last iteration to determine the amount of buckets to use in the hash map for the next iteration
	old_read_cache.Prepare(inExpectedNumBodyPairs, inExpectedNumManifolds);
}

bool ContactConstraintManager::WereBodiesInContact(const BodyID &inBody1ID, const BodyID &inBody2ID) const
{
	// The body pair needs to be in the cache and it needs to have a manifold (otherwise it's just a record indicating that there are no collisions)
	const ManifoldCache &read_cache = mCache[mCacheWriteIdx ^ 1];
	BodyPair key;
	if (inBody1ID < inBody2ID)
		key = BodyPair(inBody1ID, inBody2ID);
	else
		key = BodyPair(inBody2ID, inBody1ID);
	uint64 key_hash = key.GetHash();
	const BPKeyValue *kv = read_cache.Find(key, key_hash);
	return kv != nullptr && kv->GetValue().mFirstCachedManifold != ManifoldMap::cInvalidHandle;
}

bool ContactConstraintManager::ClearVPhysicsFrictionState(const BodyID &inBody1ID, const SubShapeID &inSubShapeID1, const BodyID &inBody2ID, const SubShapeID &inSubShapeID2)
{
	JPH_ASSERT(mUpdateContext == nullptr, "Cannot modify the contact cache during PhysicsSystem::Update");

	ManifoldCache &read_cache = mCache[mCacheWriteIdx ^ 1];
	auto find_manifold = [&read_cache](const SubShapeIDPair &inKey) -> MKeyValue *
	{
		return const_cast<MKeyValue *>(read_cache.Find(inKey, inKey.GetHash()));
	};

	SubShapeIDPair key(inBody1ID, inSubShapeID1, inBody2ID, inSubShapeID2);
	MKeyValue *kv = find_manifold(key);
	if (kv == nullptr)
	{
		// Friction snapshots are published from both bodies' point of view, while
		// Jolt stores only the narrow phase's original body ordering.
		key = SubShapeIDPair(inBody2ID, inSubShapeID2, inBody1ID, inSubShapeID1);
		kv = find_manifold(key);
	}
	if (kv == nullptr)
		return false;

	CachedManifold &manifold = kv->GetValue();
	for (uint32 point_idx = 0; point_idx < manifold.mNumContactPoints; ++point_idx)
	{
		CachedContactPoint &point = manifold.mContactPoints[point_idx];
		point.mFrictionLambda[0] = 0.0f;
		point.mFrictionLambda[1] = 0.0f;
		Vec3::sZero().StoreFloat3(&point.mVPhysicsFrictionDisplacement);
	}
	return true;
}

template <EMotionType Type1, EMotionType Type2>
void ContactConstraintManager::sGetVelocities(const MotionProperties *inMotionProperties1, const MotionProperties *inMotionProperties2, Vec3 &outLinearVelocity1, Vec3 &outAngularVelocity1, Vec3 &outLinearVelocity2, Vec3 &outAngularVelocity2)
{
	if constexpr (Type1 != EMotionType::Static)
	{
		outLinearVelocity1 = inMotionProperties1->GetLinearVelocity();
		outAngularVelocity1 = inMotionProperties1->GetAngularVelocity();
	}
	else
	{
		JPH_IF_DEBUG(outLinearVelocity1 = Vec3::sNaN();)
		JPH_IF_DEBUG(outAngularVelocity1 = Vec3::sNaN();)
	}

	if constexpr (Type2 != EMotionType::Static)
	{
		outLinearVelocity2 = inMotionProperties2->GetLinearVelocity();
		outAngularVelocity2 = inMotionProperties2->GetAngularVelocity();
	}
	else
	{
		JPH_IF_DEBUG(outLinearVelocity2 = Vec3::sNaN();)
		JPH_IF_DEBUG(outAngularVelocity2 = Vec3::sNaN();)
	}
}

template <EMotionType Type1, EMotionType Type2>
void ContactConstraintManager::sSetVelocities(MotionProperties *ioMotionProperties1, MotionProperties *ioMotionProperties2, Vec3Arg inLinearVelocity1, Vec3Arg inAngularVelocity1, Vec3Arg inLinearVelocity2, Vec3Arg inAngularVelocity2)
{
	if constexpr (Type1 == EMotionType::Dynamic)
	{
		ioMotionProperties1->ApplyLinearVelocityStep(inLinearVelocity1);
		ioMotionProperties1->ApplyAngularVelocityStep(inAngularVelocity1);
	}

	if constexpr (Type2 == EMotionType::Dynamic)
	{
		ioMotionProperties2->ApplyLinearVelocityStep(inLinearVelocity2);
		ioMotionProperties2->ApplyAngularVelocityStep(inAngularVelocity2);
	}
}

template <EMotionType Type1, EMotionType Type2, bool WarmStartFriction, bool WarmStartNonPenetration>
void ContactConstraintManager::sWarmStartConstraint(ContactConstraintBase &ioConstraint, MotionProperties *ioMotionProperties1, MotionProperties *ioMotionProperties2, float inWarmStartImpulseRatio)
{
	ContactConstraint<Type1, Type2> &constraint = static_cast<ContactConstraint<Type1, Type2> &>(ioConstraint);

	bool any_impulse_applied = false;

	// Calculate tangents
	Vec3 t1, t2;
	constraint.GetTangents(t1, t2);

	// Get velocities
	Vec3 linear_velocity1, angular_velocity1, linear_velocity2, angular_velocity2;
	sGetVelocities<Type1, Type2>(ioMotionProperties1, ioMotionProperties2, linear_velocity1, angular_velocity1, linear_velocity2, angular_velocity2);

	Vec3 ws_normal = constraint.GetWorldSpaceNormal();

	for (uint32 i = 0; i < constraint.mNumContactPoints; ++i)
	{
		WorldContactPoint<Type1, Type2> &wcp = constraint.mContactPoints[i];

		// Warm starting: Apply impulse from last frame
		if constexpr (WarmStartFriction)
		{
			if (wcp.mFrictionConstraint1.IsActive() || wcp.mFrictionConstraint2.IsActive())
			{
				if (wcp.mFrictionConstraint1.WarmStart(linear_velocity1, angular_velocity1, linear_velocity2, angular_velocity2, constraint.mInvMass1, constraint.mInvMass2, t1, inWarmStartImpulseRatio))
					any_impulse_applied = true;
				if (wcp.mFrictionConstraint2.WarmStart(linear_velocity1, angular_velocity1, linear_velocity2, angular_velocity2, constraint.mInvMass1, constraint.mInvMass2, t2, inWarmStartImpulseRatio))
					any_impulse_applied = true;
			}
		}
		if constexpr (WarmStartNonPenetration)
			if (wcp.mNonPenetrationConstraint.WarmStart(linear_velocity1, angular_velocity1, linear_velocity2, angular_velocity2, constraint.mInvMass1, constraint.mInvMass2, ws_normal, inWarmStartImpulseRatio))
				any_impulse_applied = true;
	}

	// Apply changed velocities
	if (any_impulse_applied)
		sSetVelocities<Type1, Type2>(ioMotionProperties1, ioMotionProperties2, linear_velocity1, angular_velocity1, linear_velocity2, angular_velocity2);
}

template <bool WarmStartFriction, bool WarmStartNonPenetration, class MotionPropertiesCallback>
void ContactConstraintManager::TemplatedWarmStartVelocityConstraints(const uint32 *inConstraintOffsetBegin, const uint32 *inConstraintOffsetEnd, float inWarmStartImpulseRatio, MotionPropertiesCallback &ioCallback)
{
	JPH_PROFILE_FUNCTION();

	// Build dispatch table
	using DispatchFunc = void (*)(ContactConstraintBase &, MotionProperties *, MotionProperties *, float);
	static const DispatchFunc table[3][3] = {
		{
			nullptr, // Static vs static doesn't exist
			nullptr, // Static vs kinematic doesn't exist
			sWarmStartConstraint<EMotionType::Static, EMotionType::Dynamic, WarmStartFriction, WarmStartNonPenetration>
		},
		{
			nullptr, // Kinematic vs static doesn't exist
			nullptr, // Kinematic vs kinematic doesn't exist
			sWarmStartConstraint<EMotionType::Kinematic, EMotionType::Dynamic, WarmStartFriction, WarmStartNonPenetration>
		},
		{
			sWarmStartConstraint<EMotionType::Dynamic, EMotionType::Static, WarmStartFriction, WarmStartNonPenetration>,
			sWarmStartConstraint<EMotionType::Dynamic, EMotionType::Kinematic, WarmStartFriction, WarmStartNonPenetration>,
			sWarmStartConstraint<EMotionType::Dynamic, EMotionType::Dynamic, WarmStartFriction, WarmStartNonPenetration>
		}
	};

	if (inConstraintOffsetBegin >= inConstraintOffsetEnd)
		return;

	ContactConstraintBase *next_constraint = reinterpret_cast<ContactConstraintBase *>(mConstraints + *inConstraintOffsetBegin);
	for (const uint32 *next_constraint_offset = inConstraintOffsetBegin + 1; next_constraint != nullptr; ++next_constraint_offset)
	{
		ContactConstraintBase &constraint = *next_constraint;
		if (next_constraint_offset < inConstraintOffsetEnd)
		{
			next_constraint = reinterpret_cast<ContactConstraintBase *>(mConstraints + *next_constraint_offset);
			PrefetchL1(next_constraint);
		}
		else
			next_constraint = nullptr;

		// Dispatch to the correct templated form
		Body &body1 = *constraint.mBody1;
		Body &body2 = *constraint.mBody2;
		MotionProperties *motion_properties1 = body1.GetMotionPropertiesUnchecked();
		MotionProperties *motion_properties2 = body2.GetMotionPropertiesUnchecked();
		table[(int)body1.GetMotionType()][(int)body2.GetMotionType()](constraint, motion_properties1, motion_properties2, inWarmStartImpulseRatio);

		// Call callbacks
		if (body1.IsDynamic())
			ioCallback(motion_properties1);
		if (body2.IsDynamic())
			ioCallback(motion_properties2);
	}
}

template <class MotionPropertiesCallback>
void ContactConstraintManager::WarmStartVelocityConstraints(const uint32 *inConstraintOffsetBegin, const uint32 *inConstraintOffsetEnd, float inWarmStartImpulseRatio, MotionPropertiesCallback &ioCallback)
{
	TemplatedWarmStartVelocityConstraints<true, true>(inConstraintOffsetBegin, inConstraintOffsetEnd, inWarmStartImpulseRatio, ioCallback);
}

template <class MotionPropertiesCallback>
void ContactConstraintManager::WarmStartFrictionConstraints(const uint32 *inConstraintOffsetBegin, const uint32 *inConstraintOffsetEnd, float inWarmStartImpulseRatio, MotionPropertiesCallback &ioCallback)
{
	TemplatedWarmStartVelocityConstraints<true, false>(inConstraintOffsetBegin, inConstraintOffsetEnd, inWarmStartImpulseRatio, ioCallback);
}

// Specialize for the two body callback types
template void ContactConstraintManager::WarmStartVelocityConstraints<CalculateSolverSteps>(const uint32 *inConstraintOffsetBegin, const uint32 *inConstraintOffsetEnd, float inWarmStartImpulseRatio, CalculateSolverSteps &ioCallback);
template void ContactConstraintManager::WarmStartVelocityConstraints<DummyCalculateSolverSteps>(const uint32 *inConstraintOffsetBegin, const uint32 *inConstraintOffsetEnd, float inWarmStartImpulseRatio, DummyCalculateSolverSteps &ioCallback);
template void ContactConstraintManager::WarmStartFrictionConstraints<CalculateSolverSteps>(const uint32 *inConstraintOffsetBegin, const uint32 *inConstraintOffsetEnd, float inWarmStartImpulseRatio, CalculateSolverSteps &ioCallback);
template void ContactConstraintManager::WarmStartFrictionConstraints<DummyCalculateSolverSteps>(const uint32 *inConstraintOffsetBegin, const uint32 *inConstraintOffsetEnd, float inWarmStartImpulseRatio, DummyCalculateSolverSteps &ioCallback);

template <class MotionPropertiesCallback>
void ContactConstraintManager::RegisterVelocityConstraints(const uint32 *inConstraintOffsetBegin, const uint32 *inConstraintOffsetEnd, MotionPropertiesCallback &ioCallback)
{
	for (const uint32 *offset = inConstraintOffsetBegin; offset < inConstraintOffsetEnd; ++offset)
	{
		const ContactConstraintBase &constraint = *reinterpret_cast<const ContactConstraintBase *>(mConstraints + *offset);
		if (constraint.mBody1->IsDynamic())
			ioCallback(constraint.mBody1->GetMotionPropertiesUnchecked());
		if (constraint.mBody2->IsDynamic())
			ioCallback(constraint.mBody2->GetMotionPropertiesUnchecked());
	}
}

template void ContactConstraintManager::RegisterVelocityConstraints<CalculateSolverSteps>(const uint32 *inConstraintOffsetBegin, const uint32 *inConstraintOffsetEnd, CalculateSolverSteps &ioCallback);
template void ContactConstraintManager::RegisterVelocityConstraints<DummyCalculateSolverSteps>(const uint32 *inConstraintOffsetBegin, const uint32 *inConstraintOffsetEnd, DummyCalculateSolverSteps &ioCallback);

template <EMotionType Type1, EMotionType Type2, bool SolveFriction, bool SolveNonPenetration>
bool ContactConstraintManager::sSolveVelocityConstraint(ContactConstraintBase &ioConstraint, MotionProperties *ioMotionProperties1, MotionProperties *ioMotionProperties2)
{
	ContactConstraint<Type1, Type2> &constraint = static_cast<ContactConstraint<Type1, Type2> &>(ioConstraint);

	bool any_impulse_applied = false;

	// Calculate tangents
	Vec3 t1, t2;
	constraint.GetTangents(t1, t2);

	// Get velocities
	Vec3 linear_velocity1, angular_velocity1, linear_velocity2, angular_velocity2;
	sGetVelocities<Type1, Type2>(ioMotionProperties1, ioMotionProperties2, linear_velocity1, angular_velocity1, linear_velocity2, angular_velocity2);

	if constexpr (SolveNonPenetration)
	{
		// IVP's impact event runs after PSI controllers / Havana have updated speed_change and
		// computes energy return from that current velocity. Jolt normally freezes restitution
		// bias during contact setup and subtracts current-step forces. For a fresh legacy impact,
		// refresh the bias here, immediately before the normal solve, from the post-Havana speed.
		if (constraint.mContactImpacting && constraint.mCombinedRestitution > 0.0f)
		{
			const Vec3 ws_normal = constraint.GetWorldSpaceNormal();
			for (uint32 i = 0; i < constraint.mNumContactPoints; ++i)
			{
				WorldContactPoint<Type1, Type2> &wcp = constraint.mContactPoints[i];
				const float closing_speed = wcp.mNonPenetrationConstraint.GetVelocityConstraintValue(
					linear_velocity1, angular_velocity1, linear_velocity2, angular_velocity2, ws_normal);
				if (closing_speed > constraint.mMinVelocityForRestitution)
					wcp.mNonPenetrationConstraint.SetVelocityBias(
						min(wcp.mNonPenetrationConstraint.GetVelocityBias(),
							-constraint.mCombinedRestitution * closing_speed));
			}
		}
	}

	// First apply all friction constraints (non-penetration is more important than friction)
	if constexpr (SolveFriction)
	for (uint32 i = 0; i < constraint.mNumContactPoints; ++i)
	{
		WorldContactPoint<Type1, Type2> &wcp = constraint.mContactPoints[i];

		// Check if friction is enabled
		if (wcp.mFrictionConstraint1.IsActive() || wcp.mFrictionConstraint2.IsActive())
		{
			// Calculate impulse to stop motion in tangential direction
			float lambda1 = wcp.mFrictionConstraint1.SolveVelocityConstraintGetTotalLambda(linear_velocity1, angular_velocity1, linear_velocity2, angular_velocity2, t1);
			float lambda2 = wcp.mFrictionConstraint2.SolveVelocityConstraintGetTotalLambda(linear_velocity1, angular_velocity1, linear_velocity2, angular_velocity2, t2);
			float total_lambda_sq = Square(lambda1) + Square(lambda2);

			// Calculate max impulse that can be applied. Note that we're using the non-penetration impulse from the previous iteration here.
			// We do this because non-penetration is more important so is solved last (the last things that are solved in an iterative solver
			// contribute the most).
			float max_lambda_f = constraint.mCombinedFriction * wcp.mNonPenetrationConstraint.GetTotalLambda();

			// If the total lambda that we will apply is too large, scale it back
			if (total_lambda_sq > Square(max_lambda_f))
			{
				float scale = max_lambda_f / Sqrt(total_lambda_sq);
				lambda1 *= scale;
				lambda2 *= scale;
			}

			// Apply the friction impulse
			if (wcp.mFrictionConstraint1.SolveVelocityConstraintApplyLambda(linear_velocity1, angular_velocity1, linear_velocity2, angular_velocity2, constraint.mInvMass1, constraint.mInvMass2, t1, lambda1))
				any_impulse_applied = true;
			if (wcp.mFrictionConstraint2.SolveVelocityConstraintApplyLambda(linear_velocity1, angular_velocity1, linear_velocity2, angular_velocity2, constraint.mInvMass1, constraint.mInvMass2, t2, lambda2))
				any_impulse_applied = true;
		}
	}

	if constexpr (SolveNonPenetration)
	{
		Vec3 ws_normal = constraint.GetWorldSpaceNormal();

		// Then apply all non-penetration constraints
		for (uint32 i = 0; i < constraint.mNumContactPoints; ++i)
		{
			WorldContactPoint<Type1, Type2> &wcp = constraint.mContactPoints[i];

			// Calculate impulse
			float total_lambda = wcp.mNonPenetrationConstraint.SolveVelocityConstraintGetTotalLambda(linear_velocity1, angular_velocity1, linear_velocity2, angular_velocity2, ws_normal);

			// Contact constraints can only push and not pull
			total_lambda = max(total_lambda, 0.0f);

			// Apply impulse
			if (wcp.mNonPenetrationConstraint.SolveVelocityConstraintApplyLambda(linear_velocity1, angular_velocity1, linear_velocity2, angular_velocity2, constraint.mInvMass1, constraint.mInvMass2, ws_normal, total_lambda))
				any_impulse_applied = true;
		}
	}

	if (!any_impulse_applied)
		return false;

	sSetVelocities<Type1, Type2>(ioMotionProperties1, ioMotionProperties2, linear_velocity1, angular_velocity1, linear_velocity2, angular_velocity2);
	return true;
}

template <EMotionType Type1, EMotionType Type2>
bool ContactConstraintManager::sSolveVPhysicsImpactConstraint(ContactConstraintBase &ioConstraint, MotionProperties *ioMotionProperties1, MotionProperties *ioMotionProperties2, float inDeltaTime, float inContactSeparationTarget, float inInterpenetrationRescueSpeed, uint32 &ioRemainingSystemImpacts, VPhysicsSolverStats &ioStats)
{
	ContactConstraint<Type1, Type2> &constraint = static_cast<ContactConstraint<Type1, Type2> &>(ioConstraint);
	const bool trace_impact = IsVPhysicsImpactTraceEnabled();
	if (trace_impact)
		std::fprintf(stderr, "VPHYSICS_IMPACT_ROUTE b1=%u b2=%u active=%d points=%u\n",
			constraint.mBody1->GetID().GetIndexAndSequenceNumber(),
			constraint.mBody2->GetID().GetIndexAndSequenceNumber(),
			int(constraint.mContactImpacting), constraint.mNumContactPoints);
	if (!constraint.mContactImpacting)
		return false;

	// sGetVelocities intentionally leaves static-body outputs untouched in release builds. This
	// solver evaluates relative velocity explicitly (rather than through ConstraintPart), so seed
	// every lane with zero before filling the non-static sides. Otherwise a static floor contributes
	// uninitialized stack data to the impact speed and impulse.
	Vec3 linear_velocity1 = Vec3::sZero();
	Vec3 angular_velocity1 = Vec3::sZero();
	Vec3 linear_velocity2 = Vec3::sZero();
	Vec3 angular_velocity2 = Vec3::sZero();
	sGetVelocities<Type1, Type2>(ioMotionProperties1, ioMotionProperties2,
		linear_velocity1, angular_velocity1, linear_velocity2, angular_velocity2);

	if (constraint.mSevereInterpenetration)
	{
		// IVP_Anomaly_Manager::solve_inter_penetration_simple applies one bounded
		// asynchronous speed change per body along the line between their cores.
		// Crucially, it does not turn penetration depth into a velocity and does not
		// repeat the full rescue for every point in the contact manifold.
		const Vec3 rescue_direction = Vec3(
			constraint.mBody2->GetCenterOfMassPosition()
			- constraint.mBody1->GetCenterOfMassPosition()).NormalizedOr(Vec3::sAxisX());
		// The rescue is one pair impulse, not an equal velocity assignment. Applying
		// the same speed to both bodies makes a 1 kg prop kick an 85 kg player as hard
		// as the prop itself, while IVP conserves pair momentum and lets the light body
		// take nearly all of the separation velocity.
		const float inverse_mass_sum = constraint.mInvMass1 + constraint.mInvMass2;
		const uint32 dynamic_body_count =
			(Type1 == EMotionType::Dynamic ? 1u : 0u)
			+ (Type2 == EMotionType::Dynamic ? 1u : 0u);
		const float relative_speed_change =
			float(dynamic_body_count) * inInterpenetrationRescueSpeed;
		const float rescue_impulse = inverse_mass_sum > 0.0f
			? relative_speed_change / inverse_mass_sum : 0.0f;
		bool applied = rescue_impulse > 0.0f;
		if constexpr (Type1 == EMotionType::Dynamic)
			linear_velocity1 -= rescue_direction
				* (constraint.mInvMass1 * rescue_impulse);
		if constexpr (Type2 == EMotionType::Dynamic)
			linear_velocity2 += rescue_direction
				* (constraint.mInvMass2 * rescue_impulse);
		if (applied)
			sSetVelocities<Type1, Type2>(ioMotionProperties1, ioMotionProperties2,
				linear_velocity1, angular_velocity1, linear_velocity2, angular_velocity2);
		if (trace_impact)
			std::fprintf(stderr,
				"VPHYSICS_INTERPENETRATION_RESCUE b1=%u b2=%u dir=(%.9g,%.9g,%.9g) speed=%.9g impulse=%.9g\n",
				constraint.mBody1->GetID().GetIndexAndSequenceNumber(),
				constraint.mBody2->GetID().GetIndexAndSequenceNumber(),
				rescue_direction.GetX(), rescue_direction.GetY(), rescue_direction.GetZ(),
				relative_speed_change, rescue_impulse);
		return applied;
	}

	const RMat44 transform1 = constraint.mBody1->GetCenterOfMassTransform();
	const RMat44 transform2 = constraint.mBody2->GetCenterOfMassTransform();
	const Vec3 normal = constraint.GetWorldSpaceNormal();
	Vec3 tangent1, tangent2;
	constraint.GetTangents(tangent1, tangent2);
	// JoltPhysicsEnvironment's combine callback already converts IVP's retained-energy product
	// to the velocity restitution coefficient expected by Jolt. Do not take a second square root
	// here: that would turn a low-elasticity player contact into its fourth root.
	const float material_velocity_restitution = constraint.mCombinedRestitution;
	const float impact_friction = constraint.mCombinedFriction
		* (1.0f + material_velocity_restitution);
	const float friction_cos = 1.0f / Sqrt(1.0f + Square(impact_friction));
	const float friction_sin = friction_cos * impact_friction;
	// The adapter can deliberately disable IVP's penetration recovery for contacts
	// owned by Source game movement (for example, a player already overlapping a
	// wall). The impact response must still block an incoming velocity, but it must
	// not add the separate rescue speed that would push the body out of the wall.
	const bool allow_penetration_rescue = !constraint.mDisablePenetrationRecovery;
	constexpr float cImpactEpsilon = 1.0e-4f;
	constexpr uint32 cMaximumImpactsPerManifoldPoint = 32;
	const float collision_tolerance = 0.5f * inContactSeparationTarget;
	const Vec3 velocity_before_impact1 = linear_velocity1;
	const Vec3 angular_velocity_before_impact1 = angular_velocity1;
	const Vec3 velocity_before_impact2 = linear_velocity2;
	const Vec3 angular_velocity_before_impact2 = angular_velocity2;
	// IVP keeps a committed core velocity and a pending speed_change. Impact solving sees their
	// sum, public GetVelocity sees their sum, but impact-system mindist selection deliberately
	// uses only the committed speed. This delayed-push path is essential to the mass-dependent
	// player/prop response (light props delay the player side; heavy props delay the prop side).
	Vec3 committed_linear_velocity1 = linear_velocity1;
	Vec3 committed_angular_velocity1 = angular_velocity1;
	Vec3 committed_linear_velocity2 = linear_velocity2;
	Vec3 committed_angular_velocity2 = angular_velocity2;

	// Jolt gives us a complete manifold at PSI start. IVP does not solve that array in geometry
	// order. One mindist reaches the collision tolerance first, the cores synchronize to that
	// event time, and IVP_Impact_System::test_loop_all_pairs repeatedly recalculates every affected
	// mindist and chooses the one with the smallest predicted distance. Replaying the manifold once
	// from point 0 to N leaves the torque from the first corner unopposed after later corners change
	// its velocity; a centered player push then produces a large, entirely artificial box spin.
	float earliest_impact_time = inDeltaTime;
	uint32 initial_impact_point = MaxContactPoints;
	float initial_feature_order = -FLT_MAX;
	bool late_discrete_impact_point[MaxContactPoints] = {};
	float late_discrete_impact_correction = 0.0f;
	for (uint32 point_index = 0; point_index < constraint.mNumContactPoints; ++point_index)
	{
		WorldContactPoint<Type1, Type2> &wcp = constraint.mContactPoints[point_index];
		const RVec3 point1 = transform1 * Vec3::sLoadFloat3Unsafe(wcp.mContactPoint->mPosition1);
		const RVec3 point2 = transform2 * Vec3::sLoadFloat3Unsafe(wcp.mContactPoint->mPosition2);
		const Vec3 r1 = Vec3(point1 - constraint.mBody1->GetCenterOfMassPosition());
		const Vec3 r2 = Vec3(point2 - constraint.mBody2->GetCenterOfMassPosition());
		const Vec3 relative_velocity = linear_velocity2 + angular_velocity2.Cross(r2)
			- linear_velocity1 - angular_velocity1.Cross(r1);
		const float closing_speed = -relative_velocity.Dot(normal);
		const float raw_gap = Vec3(point2 - point1).Dot(normal);
		// A body can cross IVP's collision boundary in one update while remaining
		// below Jolt's radius-scaled CCD threshold. The resulting new manifold is a
		// late ordinary impact, not a failed interpenetrating mindist. Reconstruct
		// the boundary event so a fixed negative gap cannot trigger repeated rescue.
		const float penetration_to_boundary = collision_tolerance - raw_gap;
		const bool late_discrete_impact = !constraint.mContactHasCachedSupport
			&& closing_speed > cImpactEpsilon
			&& raw_gap < collision_tolerance
			&& penetration_to_boundary <= closing_speed * inDeltaTime
				+ collision_tolerance;
		late_discrete_impact_point[point_index] = late_discrete_impact;
		if (late_discrete_impact)
			late_discrete_impact_correction = max(
				late_discrete_impact_correction, penetration_to_boundary);
		const float gap = late_discrete_impact ? collision_tolerance : raw_gap;
		float impact_time = inDeltaTime + 1.0f;
		if (gap <= collision_tolerance)
			impact_time = 0.0f;
		else if (closing_speed > cImpactEpsilon && inDeltaTime > 0.0f)
			impact_time = (gap - collision_tolerance) / closing_speed;
		const float feature_order = r1.Dot(tangent1) + r1.Dot(tangent2);
		const bool earlier = impact_time < earliest_impact_time - 1.0e-7f;
		const bool same_time_earlier_feature = abs(impact_time - earliest_impact_time) <= 1.0e-7f
			&& feature_order > initial_feature_order;
		if (impact_time >= 0.0f && (earlier || same_time_earlier_feature))
		{
			earliest_impact_time = impact_time;
			initial_impact_point = point_index;
			initial_feature_order = feature_order;
		}
	}
	if (initial_impact_point == MaxContactPoints)
		return false;
	earliest_impact_time = Clamp(earliest_impact_time, 0.0f, inDeltaTime);

	auto advance_rotation = [](QuatArg inRotation, Vec3Arg inAngularVelocity, float inTime)
	{
		const float angular_speed = inAngularVelocity.Length();
		if (angular_speed * inTime <= 1.0e-6f)
			return Quat(inRotation);
		return (Quat::sRotation(inAngularVelocity / angular_speed, angular_speed * inTime)
			* inRotation).Normalized();
	};
	const Quat event_rotation1 = advance_rotation(
		constraint.mBody1->GetRotation(), angular_velocity_before_impact1, earliest_impact_time);
	const Quat event_rotation2 = advance_rotation(
		constraint.mBody2->GetRotation(), angular_velocity_before_impact2, earliest_impact_time);
	const RVec3 event_com1 = constraint.mBody1->GetCenterOfMassPosition()
		+ velocity_before_impact1 * earliest_impact_time;
	const RVec3 event_com2 = constraint.mBody2->GetCenterOfMassPosition()
		+ velocity_before_impact2 * earliest_impact_time;
	const RMat44 event_transform1 = RMat44::sRotationTranslation(event_rotation1, event_com1);
	const RMat44 event_transform2 = RMat44::sRotationTranslation(event_rotation2, event_com2);

	Mat44 inv_i1 = Mat44::sZero();
	if constexpr (Type1 == EMotionType::Dynamic)
		inv_i1 = constraint.mInvInertiaScale1
			* ioMotionProperties1->GetInverseInertiaForRotation(Mat44::sRotation(event_rotation1));
	Mat44 inv_i2 = Mat44::sZero();
	if constexpr (Type2 == EMotionType::Dynamic)
		inv_i2 = constraint.mInvInertiaScale2
			* ioMotionProperties2->GetInverseInertiaForRotation(Mat44::sRotation(event_rotation2));

	RVec3 event_point1[MaxContactPoints];
	RVec3 event_point2[MaxContactPoints];
	Vec3 event_r1[MaxContactPoints];
	Vec3 event_r2[MaxContactPoints];
	float event_gap[MaxContactPoints] = {};
	Vec3 total_push_on_body1[MaxContactPoints];
	uint32 impacts_while_system[MaxContactPoints] = {};
	bool point_impacted[MaxContactPoints] = {};
	for (uint32 point_index = 0; point_index < constraint.mNumContactPoints; ++point_index)
	{
		WorldContactPoint<Type1, Type2> &wcp = constraint.mContactPoints[point_index];
		event_point1[point_index] = event_transform1
			* Vec3::sLoadFloat3Unsafe(wcp.mContactPoint->mPosition1);
		event_point2[point_index] = event_transform2
			* Vec3::sLoadFloat3Unsafe(wcp.mContactPoint->mPosition2);
		event_r1[point_index] = Vec3(event_point1[point_index] - event_com1);
		event_r2[point_index] = Vec3(event_point2[point_index] - event_com2);
		const float raw_event_gap = Vec3(
			event_point2[point_index] - event_point1[point_index]).Dot(normal);
		event_gap[point_index] = late_discrete_impact_point[point_index]
			? max(collision_tolerance, raw_event_gap) : raw_event_gap;
		total_push_on_body1[point_index] = Vec3::sZero();
	}

	bool any_impulse_applied = false;
	auto get_relative_velocity_from = [&](uint32 inPointIndex,
		Vec3Arg inLinearVelocity1, Vec3Arg inAngularVelocity1,
		Vec3Arg inLinearVelocity2, Vec3Arg inAngularVelocity2)
	{
		const Vec3 velocity1 = inLinearVelocity1
			+ inAngularVelocity1.Cross(event_r1[inPointIndex]);
		const Vec3 velocity2 = inLinearVelocity2
			+ inAngularVelocity2.Cross(event_r2[inPointIndex]);
		return velocity2 - velocity1;
	};
	auto get_relative_velocity = [&](uint32 inPointIndex)
	{
		return get_relative_velocity_from(inPointIndex,
			linear_velocity1, angular_velocity1, linear_velocity2, angular_velocity2);
	};
	auto get_committed_relative_velocity = [&](uint32 inPointIndex)
	{
		return get_relative_velocity_from(inPointIndex,
			committed_linear_velocity1, committed_angular_velocity1,
			committed_linear_velocity2, committed_angular_velocity2);
	};
	auto apply_source_push = [&](uint32 inPointIndex, Vec3Arg inPushOnBody1)
	{
		if constexpr (Type1 == EMotionType::Dynamic)
		{
			linear_velocity1 += constraint.mInvMass1 * inPushOnBody1;
			angular_velocity1 += inv_i1.Multiply3x3(
				event_r1[inPointIndex].Cross(inPushOnBody1));
		}
		if constexpr (Type2 == EMotionType::Dynamic)
		{
			linear_velocity2 -= constraint.mInvMass2 * inPushOnBody1;
			angular_velocity2 -= inv_i2.Multiply3x3(
				event_r2[inPointIndex].Cross(inPushOnBody1));
		}
	};
	auto solve_impact_point = [&](uint32 inPointIndex, uint32 inImpactsWhileSystem)
	{
		++ioStats.mImpactCount;
		++ioStats.mSystemImpactCount;
		WorldContactPoint<Type1, Type2> &wcp = constraint.mContactPoints[inPointIndex];
		const Vec3 &r1 = event_r1[inPointIndex];
		const Vec3 &r2 = event_r2[inPointIndex];
		const float initial_closing_speed = -get_relative_velocity(inPointIndex).Dot(normal);
		const float gap = event_gap[inPointIndex];
		if (trace_impact)
		{
			const Vec3 relative_velocity = get_relative_velocity(inPointIndex);
			std::fprintf(stderr,
				"VPHYSICS_IMPACT_INPUT point=%u system=%u normal=(%.9g,%.9g,%.9g) rel=(%.9g,%.9g,%.9g) closing=%.9g bias=%.9g "
				"v1=(%.9g,%.9g,%.9g) w1=(%.9g,%.9g,%.9g) r1=(%.9g,%.9g,%.9g) "
				"v2=(%.9g,%.9g,%.9g) w2=(%.9g,%.9g,%.9g) r2=(%.9g,%.9g,%.9g)\n",
				inPointIndex, inImpactsWhileSystem, normal.GetX(), normal.GetY(), normal.GetZ(),
				relative_velocity.GetX(), relative_velocity.GetY(), relative_velocity.GetZ(),
				initial_closing_speed, wcp.mNonPenetrationConstraint.GetVelocityBias(),
				linear_velocity1.GetX(), linear_velocity1.GetY(), linear_velocity1.GetZ(),
				angular_velocity1.GetX(), angular_velocity1.GetY(), angular_velocity1.GetZ(),
				r1.GetX(), r1.GetY(), r1.GetZ(),
				linear_velocity2.GetX(), linear_velocity2.GetY(), linear_velocity2.GetZ(),
				angular_velocity2.GetX(), angular_velocity2.GetY(), angular_velocity2.GetZ(),
				r2.GetX(), r2.GetY(), r2.GetZ());
		}

		const Vec3 r1_cross_normal = r1.Cross(normal);
		const Vec3 r2_cross_normal = r2.Cross(normal);
		float inverse_virtual_mass1 = constraint.mInvMass1;
		float inverse_virtual_mass2 = constraint.mInvMass2;
		if constexpr (Type1 == EMotionType::Dynamic)
			inverse_virtual_mass1 += inv_i1.Multiply3x3(r1_cross_normal).Dot(r1_cross_normal);
		if constexpr (Type2 == EMotionType::Dynamic)
			inverse_virtual_mass2 += inv_i2.Multiply3x3(r2_cross_normal).Dot(r2_cross_normal);
		const float inverse_virtual_mass = inverse_virtual_mass1 + inverse_virtual_mass2;
		if (inverse_virtual_mass <= 0.0f)
			return false;
		const bool repeated_mindist = point_impacted[inPointIndex];
		point_impacted[inPointIndex] = true;

		// IVP estimates the perfectly elastic normal impulse, then advances in fixed 10% pushes
		// while recomputing the full relative-speed direction after every push. Because tangential
		// effective mass differs from normal effective mass, this intentionally overshoots and
		// changes direction; replacing it with a one-shot Coulomb projection changes wheel spin.
		const float push_step = initial_closing_speed > cImpactEpsilon
			? 0.2f * initial_closing_speed / inverse_virtual_mass : 0.0f;
		if (trace_impact)
			std::fprintf(stderr, "VPHYSICS_IMPACT_MASS point=%u inverse_mass=%.9g push_step=%.9g\n",
				inPointIndex, inverse_virtual_mass, push_step);
		Vec3 point_push = Vec3::sZero();
		float virtual_speed = initial_closing_speed;
		int push_count = 0;
		while (virtual_speed > cImpactEpsilon && push_count < 100)
		{
			++push_count;
			Vec3 relative_velocity = get_relative_velocity(inPointIndex);
			const float relative_length_sq = relative_velocity.LengthSq();
			if (relative_length_sq <= 1.0e-16f)
				break;
			Vec3 push_direction = relative_velocity / Sqrt(relative_length_sq);
			const float normal_part = push_direction.Dot(normal);
			if (normal_part > 0.0f)
				break;
			if (normal_part > -friction_cos)
			{
				Vec3 surface_part = push_direction - normal * normal_part;
				const float surface_length_sq = surface_part.LengthSq();
				if (surface_length_sq > 1.0e-16f)
					push_direction = surface_part * (friction_sin / Sqrt(surface_length_sq))
						- normal * friction_cos;
				else
					push_direction = -normal;
			}

			const Vec3 push = push_step * push_direction;
			apply_source_push(inPointIndex, push);
			point_push += push;
			any_impulse_applied = true;
			virtual_speed = -get_relative_velocity(inPointIndex).Dot(normal);
		}

		// A new Jolt manifold point maps to a new IVP contact feature and enters with the material
		// restitution. Only when the impact system returns to the same point do we apply IVP's
		// per-mindist conservation ramp and system-impact speed add-on.
		const float material_conservation = Square(constraint.mCombinedRestitution);
		const float used_conservation = repeated_mindist
			? 1.0f - (1.0f / (0.5f * float(inImpactsWhileSystem) + 1.0f))
				* (1.0f - material_conservation)
			: material_conservation;
		float give_back_speed = initial_closing_speed > cImpactEpsilon
			? initial_closing_speed * Sqrt(max(0.0f, used_conservation)) + virtual_speed
				+ (repeated_mindist ? 0.01f * Sqrt(float(inImpactsWhileSystem)) : 0.0f)
			: 0.0f;
		if (give_back_speed > cImpactEpsilon && push_count != 100)
		{
			const Vec3 push = -normal * (give_back_speed / inverse_virtual_mass);
			apply_source_push(inPointIndex, push);
			point_push += push;
			any_impulse_applied = true;
		}

		// Preserve IVP's minimum opening speed and penetration rescue when the adapter has
		// not explicitly disabled that recovery path.
		float rescue_speed_addon = 0.0f;
		float desired_opening_speed = 0.0f;
		if (allow_penetration_rescue)
		{
			// Preserve IVP's minimum opening speed and penetration rescue. Do not reuse
			// Jolt's precomputed velocity bias here: it is based on a different separation
			// target and was frozen before Havana changed the body's velocity.
			rescue_speed_addon = inDeltaTime > 0.0f
				? 2.0f * min(collision_tolerance,
					max(0.0f, collision_tolerance - gap)) / inDeltaTime
				: 0.0f;
			desired_opening_speed = 1.2f
				* (inContactSeparationTarget + rescue_speed_addon);
			const float current_opening_speed = get_relative_velocity(inPointIndex).Dot(normal);
			// If the committed-speed scan selected this mindist while the effective speed
			// (including a delayed speed_change) is already separating, IVP enters panic
			// rescue mode. It clamps the existing separating speed to zero for the
			// calculation and still adds one full rescue-speed push; without this, the same
			// mindist is selected forever with a zero impulse.
			const bool panic_rescue = initial_closing_speed <= cImpactEpsilon;
			const float opening_speed_change = panic_rescue
				? desired_opening_speed
				: max(0.0f, desired_opening_speed - current_opening_speed);
			if (opening_speed_change > 0.0f)
			{
				const Vec3 push = -normal
					* (opening_speed_change / inverse_virtual_mass);
				apply_source_push(inPointIndex, push);
				point_push += push;
				any_impulse_applied = true;
			}
		}
		total_push_on_body1[inPointIndex] += point_push;

		// Match IVP_Impact_Solver::delay_decision. The side with the larger virtual mass is the
		// candidate for delaying. If committing only the other side would still leave the contact
		// closing faster than the rescue threshold, keep this side's delta pending and expose only
		// the other side to the next impact-system distance scan.
		if constexpr (Type1 == EMotionType::Dynamic && Type2 == EMotionType::Dynamic)
		{
			// calc_correct_virt_mass uses the length of the complete contact-point velocity
			// response to a unit push, not merely its projection back onto the normal.
			const Vec3 response1 = constraint.mInvMass1 * normal
				+ inv_i1.Multiply3x3(r1.Cross(normal)).Cross(r1);
			const Vec3 response2 = -constraint.mInvMass2 * normal
				+ inv_i2.Multiply3x3(r2.Cross(-normal)).Cross(r2);
			const float response_length1 = response1.Length();
			const float response_length2 = response2.Length();
			const float virtual_mass1 = response_length1 > 0.0f
				? 1.0f / response_length1 : FLT_MAX;
			const float virtual_mass2 = response_length2 > 0.0f
				? 1.0f / response_length2 : FLT_MAX;
			const bool delay_body1 = virtual_mass1 > virtual_mass2;
			const Vec3 delayed_test_relative_velocity = delay_body1
				? get_relative_velocity_from(inPointIndex,
					committed_linear_velocity1, committed_angular_velocity1,
					linear_velocity2, angular_velocity2)
				: get_relative_velocity_from(inPointIndex,
					linear_velocity1, angular_velocity1,
					committed_linear_velocity2, committed_angular_velocity2);
			// IVP's helper returns candidate-minus-other and flips the result again for core 1.
			// Our common convention is body2-minus-body1, so the delay test uses the opposite sign.
			const float delay_closing_speed = -delayed_test_relative_velocity.Dot(normal);
			const bool delay_selected_body = delay_closing_speed
				< -(inContactSeparationTarget + rescue_speed_addon);
			if (delay_selected_body)
			{
				++ioStats.mDelayedImpactCount;
				if (delay_body1)
				{
					committed_linear_velocity2 = linear_velocity2;
					committed_angular_velocity2 = angular_velocity2;
				}
				else
				{
					committed_linear_velocity1 = linear_velocity1;
					committed_angular_velocity1 = angular_velocity1;
				}
			}
			else
			{
				committed_linear_velocity1 = linear_velocity1;
				committed_angular_velocity1 = angular_velocity1;
				committed_linear_velocity2 = linear_velocity2;
				committed_angular_velocity2 = angular_velocity2;
			}
			if (trace_impact)
				std::fprintf(stderr,
					"VPHYSICS_IMPACT_DELAY point=%u candidate=%u delayed=%d closing=%.9g vm1=%.9g vm2=%.9g\n",
					inPointIndex, delay_body1 ? 1u : 2u, int(delay_selected_body),
					delay_closing_speed, virtual_mass1, virtual_mass2);
		}
		else
		{
			committed_linear_velocity1 = linear_velocity1;
			committed_angular_velocity1 = angular_velocity1;
			committed_linear_velocity2 = linear_velocity2;
			committed_angular_velocity2 = angular_velocity2;
		}
		if (trace_impact)
		{
			const Vec3 relative_velocity = get_relative_velocity(inPointIndex);
			std::fprintf(stderr,
				"VPHYSICS_IMPACT_OUTPUT point=%u pushes=%d total=(%.9g,%.9g,%.9g) rel=(%.9g,%.9g,%.9g) gap=%.9g rescue=%.9g target=%.9g\n",
				inPointIndex, push_count, point_push.GetX(), point_push.GetY(),
				point_push.GetZ(), relative_velocity.GetX(), relative_velocity.GetY(),
				relative_velocity.GetZ(), gap, rescue_speed_addon, desired_opening_speed);
		}
		return !point_push.IsNearZero();
	};

	// The narrow-phase event itself is the first impact. Subsequent selections are the IVP impact
	// system and therefore use a per-mindist repeat count starting at one.
	solve_impact_point(initial_impact_point, 0);
	const uint32 manifold_impact_budget = min(ioRemainingSystemImpacts,
		cMaximumImpactsPerManifoldPoint * max(1u, constraint.mNumContactPoints));
	uint32 consumed_system_impacts = 0;
	while (consumed_system_impacts < manifold_impact_budget)
	{
		uint32 worst_point = MaxContactPoints;
		float smallest_predicted_distance = collision_tolerance;
		ioStats.mImpactCollisionChecks += constraint.mNumContactPoints;
		for (uint32 point_offset = constraint.mNumContactPoints; point_offset > 0; --point_offset)
		{
			const uint32 point_index = point_offset - 1;
			const float closing_speed = -get_committed_relative_velocity(point_index).Dot(normal);
			const float rescue_speed_addon = allow_penetration_rescue && inDeltaTime > 0.0f
				? 2.0f * min(collision_tolerance,
					max(0.0f, collision_tolerance - event_gap[point_index])) / inDeltaTime
				: 0.0f;
			const float predicted_distance = event_gap[point_index]
				- (closing_speed + 0.5f * rescue_speed_addon) * inDeltaTime;
			if (predicted_distance < smallest_predicted_distance)
			{
				smallest_predicted_distance = predicted_distance;
				worst_point = point_index;
			}
		}
		if (worst_point == MaxContactPoints)
			break;
		++impacts_while_system[worst_point];
		++consumed_system_impacts;
		if (trace_impact)
			std::fprintf(stderr,
				"VPHYSICS_IMPACT_SELECT iteration=%u point=%u predicted=%.9g count=%u\n",
				consumed_system_impacts - 1, worst_point, smallest_predicted_distance,
				impacts_while_system[worst_point]);
		if (!solve_impact_point(worst_point, impacts_while_system[worst_point]))
			break;
	}
	ioRemainingSystemImpacts -= consumed_system_impacts;
	if (trace_impact && consumed_system_impacts == manifold_impact_budget
		&& manifold_impact_budget > 0)
		std::fprintf(stderr,
			"VPHYSICS_IMPACT_BUDGET b1=%u b2=%u consumed=%u remaining=%u points=%u\n",
			constraint.mBody1->GetID().GetIndexAndSequenceNumber(),
			constraint.mBody2->GetID().GetIndexAndSequenceNumber(),
			consumed_system_impacts, ioRemainingSystemImpacts, constraint.mNumContactPoints);

	// Jolt's lambda sign is the impulse applied to body 2; IVP accumulated the opposite push on
	// body 1. Publish the complete per-mindist sums only after the impact system has converged.
	for (uint32 point_index = 0; point_index < constraint.mNumContactPoints; ++point_index)
	{
		WorldContactPoint<Type1, Type2> &wcp = constraint.mContactPoints[point_index];
		const Vec3 &push = total_push_on_body1[point_index];
		wcp.mNonPenetrationConstraint.SetTotalLambda(max(0.0f, -push.Dot(normal)));
		wcp.mFrictionConstraint1.SetTotalLambda(-push.Dot(tangent1));
		wcp.mFrictionConstraint2.SetTotalLambda(-push.Dot(tangent2));
	}

	if (any_impulse_applied)
	{
		// Rewind a one-frame discrete overshoot to IVP's impact boundary. This is a
		// position reconstruction of the missed event, not penetration velocity:
		// distribute it by inverse mass and preserve the solver's outgoing speeds.
		// Genuine deep teleports fail the one-frame plausibility test and remain on
		// the bounded anomaly-rescue path.
		if (late_discrete_impact_correction > 0.0f)
		{
			const float inverse_mass_sum = constraint.mInvMass1 + constraint.mInvMass2;
			if (inverse_mass_sum > 0.0f)
			{
				if constexpr (Type1 == EMotionType::Dynamic)
					ioMotionProperties1->AddVPhysicsIntegrationCorrection(
						-normal * (late_discrete_impact_correction
							* constraint.mInvMass1 / inverse_mass_sum), Vec3::sZero());
				if constexpr (Type2 == EMotionType::Dynamic)
					ioMotionProperties2->AddVPhysicsIntegrationCorrection(
						normal * (late_discrete_impact_correction
							* constraint.mInvMass2 / inverse_mass_sum), Vec3::sZero());
			}
		}

		// A friction contact created by this event only exists for the remainder of the PSI.
		// Seed its spring with exactly that post-impact time slice; subsequent persistent frames
		// accumulate their full interval after all velocity controllers have run.
		for (uint32 point_index = 0; point_index < constraint.mNumContactPoints; ++point_index)
			if (point_impacted[point_index])
			{
				WorldContactPoint<Type1, Type2> &wcp = constraint.mContactPoints[point_index];
				const float relative_velocity1 = wcp.mFrictionConstraint1.GetVelocityConstraintValue(
					linear_velocity1, angular_velocity1, linear_velocity2, angular_velocity2, tangent1)
					- wcp.mFrictionConstraint1.GetVelocityBias();
				const float relative_velocity2 = wcp.mFrictionConstraint2.GetVelocityConstraintValue(
					linear_velocity1, angular_velocity1, linear_velocity2, angular_velocity2, tangent2)
					- wcp.mFrictionConstraint2.GetVelocityBias();
				const float remaining_time = max(0.0f, inDeltaTime - earliest_impact_time);
				const Vec3 displacement = (tangent1 * relative_velocity1
					+ tangent2 * relative_velocity2) * remaining_time;
				displacement.StoreFloat3(&wcp.mContactPoint->mVPhysicsFrictionDisplacement);
			}

		// IVP advances to the impact event with the incoming velocity, applies the impulse there,
		// then advances the remainder of the PSI with the outgoing velocity. Jolt's symplectic
		// integrator would otherwise apply the outgoing velocity to the entire step (time stealing).
		// Queue the missing incoming segment as a transient correction consumed during integration.
		if constexpr (Type1 == EMotionType::Dynamic)
			ioMotionProperties1->AddVPhysicsIntegrationCorrection(
				(velocity_before_impact1 - linear_velocity1) * earliest_impact_time,
				(angular_velocity_before_impact1 - angular_velocity1) * earliest_impact_time);
		if constexpr (Type2 == EMotionType::Dynamic)
			ioMotionProperties2->AddVPhysicsIntegrationCorrection(
				(velocity_before_impact2 - linear_velocity2) * earliest_impact_time,
				(angular_velocity_before_impact2 - angular_velocity2) * earliest_impact_time);
		sSetVelocities<Type1, Type2>(ioMotionProperties1, ioMotionProperties2,
			linear_velocity1, angular_velocity1, linear_velocity2, angular_velocity2);
	}
	return any_impulse_applied;
}

template <EMotionType Type1, EMotionType Type2>
bool ContactConstraintManager::sSolveLegacyFrictionConstraint(ContactConstraintBase &ioConstraint, MotionProperties *ioMotionProperties1, MotionProperties *ioMotionProperties2, float inDeltaTime)
{
	if (inDeltaTime <= 0.0f)
		return false;

	ContactConstraint<Type1, Type2> &constraint = static_cast<ContactConstraint<Type1, Type2> &>(ioConstraint);
	if (constraint.mContactImpacting)
		return false;
	Vec3 t1, t2;
	constraint.GetTangents(t1, t2);
	const Vec3 ws_normal = constraint.GetWorldSpaceNormal();
	const RMat44 transform1 = constraint.mBody1->GetCenterOfMassTransform();
	const RMat44 transform2 = constraint.mBody2->GetCenterOfMassTransform();
	Mat44 inv_i1 = Mat44::sZero();
	if constexpr (Type1 == EMotionType::Dynamic)
		inv_i1 = constraint.mInvInertiaScale1 * constraint.mBody1->GetInverseInertia();
	Mat44 inv_i2 = Mat44::sZero();
	if constexpr (Type2 == EMotionType::Dynamic)
		inv_i2 = constraint.mInvInertiaScale2 * constraint.mBody2->GetInverseInertia();
	const bool trace_friction = IsLegacyFrictionTraceEnabled();

	Vec3 linear_velocity1, angular_velocity1, linear_velocity2, angular_velocity2;
	sGetVelocities<Type1, Type2>(ioMotionProperties1, ioMotionProperties2, linear_velocity1, angular_velocity1, linear_velocity2, angular_velocity2);

	bool any_impulse_applied = false;
	for (uint32 i = 0; i < constraint.mNumContactPoints; ++i)
	{
		WorldContactPoint<Type1, Type2> &wcp = constraint.mContactPoints[i];
		if (!wcp.mFrictionConstraint1.IsActive() && !wcp.mFrictionConstraint2.IsActive())
			continue;

		Vec3 displacement = Vec3::sLoadFloat3Unsafe(wcp.mContactPoint->mVPhysicsFrictionDisplacement);
		// Match IVP's contact replacement path: preserve the old world-space friction
		// vector, then project it into the current contact plane before integrating slip.
		displacement -= ws_normal * displacement.Dot(ws_normal);
		const float relative_velocity1 = wcp.mFrictionConstraint1.GetVelocityConstraintValue(
			linear_velocity1, angular_velocity1, linear_velocity2, angular_velocity2, t1)
			- wcp.mFrictionConstraint1.GetVelocityBias();
		const float relative_velocity2 = wcp.mFrictionConstraint2.GetVelocityConstraintValue(
			linear_velocity1, angular_velocity1, linear_velocity2, angular_velocity2, t2)
			- wcp.mFrictionConstraint2.GetVelocityBias();

		// The energy-friction controller has already accumulated the previous interval's slip.
		// IVP clips that stored tangential spring displacement to the distance that last frame's
		// normal pressure can traverse during one PSI. The cached normal lambda is impulse,
		// so pressure * dt collapses back to that lambda here.
		const float normal_lambda = max(0.0f, wcp.mContactPoint->mVPhysicsFrictionNormalLambda);
		const float max_friction_lambda = constraint.mCombinedFriction * normal_lambda;
		const float inverse_virtual_mass = max(
			wcp.mFrictionConstraint1.GetInverseEffectiveMass(),
			wcp.mFrictionConstraint2.GetInverseEffectiveMass());
		const float max_displacement = max_friction_lambda * inverse_virtual_mass * inDeltaTime;
		const float displacement_length_sq = displacement.LengthSq();
		if (displacement_length_sq > Square(max_displacement) && displacement_length_sq > 0.0f)
			displacement *= max_displacement / Sqrt(displacement_length_sq);
		displacement.StoreFloat3(&wcp.mContactPoint->mVPhysicsFrictionDisplacement);

		// IVP does not replay last frame's tangential impulse. It solves a fresh local
		// constraint that drives the integrated friction spring back toward zero.
		wcp.mFrictionConstraint1.SetTotalLambda(0.0f);
		wcp.mFrictionConstraint2.SetTotalLambda(0.0f);
		const float rhs1 = relative_velocity1 + displacement.Dot(t1) / inDeltaTime;
		const float rhs2 = relative_velocity2 + displacement.Dot(t2) / inDeltaTime;

		// IVP's local 2D friction constraint inverts the complete 2x2 response matrix. Solving
		// the two tangents as unrelated scalar constraints is only equivalent for a centered sphere;
		// an anisotropic .phy wheel has a non-zero cross term and the scalar approximation injects
		// velocity into the other wheel axis.
		const RVec3 p1 = transform1 * Vec3::sLoadFloat3Unsafe(wcp.mContactPoint->mPosition1);
		const RVec3 p2 = transform2 * Vec3::sLoadFloat3Unsafe(wcp.mContactPoint->mPosition2);
		const RVec3 p = 0.5_r * (p1 + p2);
		const Vec3 r1 = Vec3(p - constraint.mBody1->GetCenterOfMassPosition());
		const Vec3 r2 = Vec3(p - constraint.mBody2->GetCenterOfMassPosition());
		const Vec3 r1_cross_t1 = r1.Cross(t1);
		const Vec3 r1_cross_t2 = r1.Cross(t2);
		const Vec3 r2_cross_t1 = r2.Cross(t1);
		const Vec3 r2_cross_t2 = r2.Cross(t2);
		const float k00 = wcp.mFrictionConstraint1.GetInverseEffectiveMass();
		const float k11 = wcp.mFrictionConstraint2.GetInverseEffectiveMass();
		float k01 = 0.0f;
		if constexpr (Type1 == EMotionType::Dynamic)
			k01 += inv_i1.Multiply3x3(r1_cross_t1).Dot(r1_cross_t2);
		if constexpr (Type2 == EMotionType::Dynamic)
			k01 += inv_i2.Multiply3x3(r2_cross_t1).Dot(r2_cross_t2);
		const float determinant = k00 * k11 - k01 * k01;
		float lambda1;
		float lambda2;
		if (determinant > 1.0e-12f)
		{
			const float inverse_determinant = 1.0f / determinant;
			lambda1 = (k11 * rhs1 - k01 * rhs2) * inverse_determinant;
			lambda2 = (k00 * rhs2 - k01 * rhs1) * inverse_determinant;
		}
		else
		{
			lambda1 = k00 > 0.0f ? rhs1 / k00 : 0.0f;
			lambda2 = k11 > 0.0f ? rhs2 / k11 : 0.0f;
		}
		if (trace_friction)
			std::fprintf(stderr,
				"VPHYSICS_FRICTION_INPUT b1=%u b2=%u point=%u rel=(%.9g,%.9g) disp=(%.9g,%.9g,%.9g) rhs=(%.9g,%.9g) k=(%.9g,%.9g,%.9g) lambda=(%.9g,%.9g) max=%.9g\n",
				constraint.mBody1->GetID().GetIndexAndSequenceNumber(),
				constraint.mBody2->GetID().GetIndexAndSequenceNumber(), i,
				relative_velocity1, relative_velocity2,
				displacement.GetX(), displacement.GetY(), displacement.GetZ(), rhs1, rhs2,
				k00, k01, k11, lambda1, lambda2, max_friction_lambda);
		const float total_lambda_sq = Square(lambda1) + Square(lambda2);
		if (total_lambda_sq > Square(max_friction_lambda) && total_lambda_sq > 0.0f)
		{
			const float scale = max_friction_lambda / Sqrt(total_lambda_sq);
			lambda1 *= scale;
			lambda2 *= scale;
		}

		if (wcp.mFrictionConstraint1.SolveVelocityConstraintApplyLambda(
			linear_velocity1, angular_velocity1, linear_velocity2, angular_velocity2,
			constraint.mInvMass1, constraint.mInvMass2, t1, lambda1))
			any_impulse_applied = true;
		if (wcp.mFrictionConstraint2.SolveVelocityConstraintApplyLambda(
			linear_velocity1, angular_velocity1, linear_velocity2, angular_velocity2,
			constraint.mInvMass1, constraint.mInvMass2, t2, lambda2))
			any_impulse_applied = true;
		if (trace_friction)
		{
			const float output1 = wcp.mFrictionConstraint1.GetVelocityConstraintValue(
				linear_velocity1, angular_velocity1, linear_velocity2, angular_velocity2, t1)
				- wcp.mFrictionConstraint1.GetVelocityBias();
			const float output2 = wcp.mFrictionConstraint2.GetVelocityConstraintValue(
				linear_velocity1, angular_velocity1, linear_velocity2, angular_velocity2, t2)
				- wcp.mFrictionConstraint2.GetVelocityBias();
			std::fprintf(stderr, "VPHYSICS_FRICTION_OUTPUT b1=%u b2=%u point=%u rel=(%.9g,%.9g)\n",
				constraint.mBody1->GetID().GetIndexAndSequenceNumber(),
				constraint.mBody2->GetID().GetIndexAndSequenceNumber(), i, output1, output2);
		}
	}

	if (any_impulse_applied)
		sSetVelocities<Type1, Type2>(ioMotionProperties1, ioMotionProperties2, linear_velocity1, angular_velocity1, linear_velocity2, angular_velocity2);
	return any_impulse_applied;
}

template <EMotionType Type1, EMotionType Type2>
void ContactConstraintManager::sIntegrateLegacyFrictionState(ContactConstraintBase &ioConstraint, MotionProperties *ioMotionProperties1, MotionProperties *ioMotionProperties2, float inDeltaTime)
{
	if (inDeltaTime <= 0.0f)
		return;

	ContactConstraint<Type1, Type2> &constraint = static_cast<ContactConstraint<Type1, Type2> &>(ioConstraint);
	if (constraint.mContactImpacting)
		return;
	Vec3 t1, t2;
	constraint.GetTangents(t1, t2);
	const Vec3 ws_normal = constraint.GetWorldSpaceNormal();
	Vec3 linear_velocity1 = Vec3::sZero();
	Vec3 angular_velocity1 = Vec3::sZero();
	Vec3 linear_velocity2 = Vec3::sZero();
	Vec3 angular_velocity2 = Vec3::sZero();
	sGetVelocities<Type1, Type2>(ioMotionProperties1, ioMotionProperties2,
		linear_velocity1, angular_velocity1, linear_velocity2, angular_velocity2);

	for (uint32 i = 0; i < constraint.mNumContactPoints; ++i)
	{
		WorldContactPoint<Type1, Type2> &wcp = constraint.mContactPoints[i];
		Vec3 displacement = Vec3::sLoadFloat3Unsafe(wcp.mContactPoint->mVPhysicsFrictionDisplacement);
		displacement -= ws_normal * displacement.Dot(ws_normal);
		const float relative_velocity1 = wcp.mFrictionConstraint1.GetVelocityConstraintValue(
			linear_velocity1, angular_velocity1, linear_velocity2, angular_velocity2, t1)
			- wcp.mFrictionConstraint1.GetVelocityBias();
		const float relative_velocity2 = wcp.mFrictionConstraint2.GetVelocityConstraintValue(
			linear_velocity1, angular_velocity1, linear_velocity2, angular_velocity2, t2)
			- wcp.mFrictionConstraint2.GetVelocityBias();
		displacement += (t1 * relative_velocity1 + t2 * relative_velocity2) * inDeltaTime;
		displacement.StoreFloat3(&wcp.mContactPoint->mVPhysicsFrictionDisplacement);
	}
}

template <bool SolveFriction, bool SolveNonPenetration, bool PersistentOnly>
bool ContactConstraintManager::TemplatedSolveVelocityConstraints(const uint32 *inConstraintOffsetBegin, const uint32 *inConstraintOffsetEnd)
{
	JPH_PROFILE_FUNCTION();

	// Build dispatch table
	using DispatchFunc = bool (*)(ContactConstraintBase &, MotionProperties *, MotionProperties *);
	static const DispatchFunc table[3][3] = {
		{
			nullptr, // Static vs static doesn't exist
			nullptr, // Static vs kinematic doesn't exist
			sSolveVelocityConstraint<EMotionType::Static, EMotionType::Dynamic, SolveFriction, SolveNonPenetration>
		},
		{
			nullptr, // Kinematic vs static doesn't exist
			nullptr, // Kinematic vs kinematic doesn't exist
			sSolveVelocityConstraint<EMotionType::Kinematic, EMotionType::Dynamic, SolveFriction, SolveNonPenetration>
		},
		{
			sSolveVelocityConstraint<EMotionType::Dynamic, EMotionType::Static, SolveFriction, SolveNonPenetration>,
			sSolveVelocityConstraint<EMotionType::Dynamic, EMotionType::Kinematic, SolveFriction, SolveNonPenetration>,
			sSolveVelocityConstraint<EMotionType::Dynamic, EMotionType::Dynamic, SolveFriction, SolveNonPenetration>
		}
	};

	if (inConstraintOffsetBegin >= inConstraintOffsetEnd)
		return false;

	bool any_impulse_applied = false;

	ContactConstraintBase *next_constraint = reinterpret_cast<ContactConstraintBase *>(mConstraints + *inConstraintOffsetBegin);
	for (const uint32 *next_constraint_offset = inConstraintOffsetBegin + 1; next_constraint != nullptr; ++next_constraint_offset)
	{
		ContactConstraintBase &constraint = *next_constraint;
		if (next_constraint_offset < inConstraintOffsetEnd)
		{
			next_constraint = reinterpret_cast<ContactConstraintBase *>(mConstraints + *next_constraint_offset);
			PrefetchL1(next_constraint);
		}
		else
			next_constraint = nullptr;
		if constexpr (PersistentOnly)
			if (constraint.mContactImpacting)
				continue;

		// Dispatch to the correct templated form
		Body &body1 = *constraint.mBody1;
		Body &body2 = *constraint.mBody2;
		any_impulse_applied |= table[(int)body1.GetMotionType()][(int)body2.GetMotionType()](constraint, body1.GetMotionPropertiesUnchecked(), body2.GetMotionPropertiesUnchecked());
	}

	return any_impulse_applied;
}

bool ContactConstraintManager::SolveVelocityConstraints(const uint32 *inConstraintOffsetBegin, const uint32 *inConstraintOffsetEnd)
{
	return TemplatedSolveVelocityConstraints<true, true, false>(inConstraintOffsetBegin, inConstraintOffsetEnd);
}

bool ContactConstraintManager::SolveFrictionConstraints(const uint32 *inConstraintOffsetBegin, const uint32 *inConstraintOffsetEnd)
{
	return TemplatedSolveVelocityConstraints<true, false, false>(inConstraintOffsetBegin, inConstraintOffsetEnd);
}

bool ContactConstraintManager::SolveVPhysicsImpactConstraints(const uint32 *inConstraintOffsetBegin, const uint32 *inConstraintOffsetEnd, float inDeltaTime)
{
	using DispatchFunc = bool (*)(ContactConstraintBase &, MotionProperties *, MotionProperties *, float, float, float, uint32 &, VPhysicsSolverStats &);
	static const DispatchFunc table[3][3] = {
		{ nullptr, nullptr, sSolveVPhysicsImpactConstraint<EMotionType::Static, EMotionType::Dynamic> },
		{ nullptr, nullptr, sSolveVPhysicsImpactConstraint<EMotionType::Kinematic, EMotionType::Dynamic> },
		{
			sSolveVPhysicsImpactConstraint<EMotionType::Dynamic, EMotionType::Static>,
			sSolveVPhysicsImpactConstraint<EMotionType::Dynamic, EMotionType::Kinematic>,
			sSolveVPhysicsImpactConstraint<EMotionType::Dynamic, EMotionType::Dynamic>
		}
	};

	bool any_impulse_applied = false;
	VPhysicsSolverStats stats;
	// IVP's MAXIMUM_SYSTEM_PUSH_NUMBER is one budget for the complete friction impact
	// system. The adapter previously gave every Jolt manifold its own 5000-push budget,
	// turning a pasted constrained building into contacts * 5000 nested solves.
	uint32 remaining_system_impacts = 5000;
	const float interpenetration_rescue_speed = 2.0f
		* mUpdateContext->mPhysicsSystem->GetGravity().Length() * inDeltaTime;
	for (const uint32 *offset = inConstraintOffsetBegin; offset < inConstraintOffsetEnd; ++offset)
	{
		ContactConstraintBase &constraint = *reinterpret_cast<ContactConstraintBase *>(mConstraints + *offset);
		Body &body1 = *constraint.mBody1;
		Body &body2 = *constraint.mBody2;
		any_impulse_applied |= table[(int)body1.GetMotionType()][(int)body2.GetMotionType()](
			constraint, body1.GetMotionPropertiesUnchecked(), body2.GetMotionPropertiesUnchecked(),
			inDeltaTime, mPhysicsSettings.mContactSeparationTarget,
			interpenetration_rescue_speed, remaining_system_impacts, stats);
	}
	if (mContactListener != nullptr && !stats.IsEmpty())
		mContactListener->OnVPhysicsSolverStats(stats);
	return any_impulse_applied;
}

bool ContactConstraintManager::SolveLegacyFrictionConstraints(const uint32 *inConstraintOffsetBegin, const uint32 *inConstraintOffsetEnd, float inDeltaTime)
{
	using DispatchFunc = bool (*)(ContactConstraintBase &, MotionProperties *, MotionProperties *, float);
	static const DispatchFunc table[3][3] = {
		{
			nullptr,
			nullptr,
			sSolveLegacyFrictionConstraint<EMotionType::Static, EMotionType::Dynamic>
		},
		{
			nullptr,
			nullptr,
			sSolveLegacyFrictionConstraint<EMotionType::Kinematic, EMotionType::Dynamic>
		},
		{
			sSolveLegacyFrictionConstraint<EMotionType::Dynamic, EMotionType::Static>,
			sSolveLegacyFrictionConstraint<EMotionType::Dynamic, EMotionType::Kinematic>,
			sSolveLegacyFrictionConstraint<EMotionType::Dynamic, EMotionType::Dynamic>
		}
	};

	bool any_impulse_applied = false;
	for (const uint32 *offset = inConstraintOffsetBegin; offset < inConstraintOffsetEnd; ++offset)
	{
		ContactConstraintBase &constraint = *reinterpret_cast<ContactConstraintBase *>(mConstraints + *offset);
		Body &body1 = *constraint.mBody1;
		Body &body2 = *constraint.mBody2;
		any_impulse_applied |= table[(int)body1.GetMotionType()][(int)body2.GetMotionType()](
			constraint, body1.GetMotionPropertiesUnchecked(), body2.GetMotionPropertiesUnchecked(), inDeltaTime);
	}
	return any_impulse_applied;
}

void ContactConstraintManager::IntegrateLegacyFrictionState(const uint32 *inConstraintOffsetBegin, const uint32 *inConstraintOffsetEnd, float inDeltaTime)
{
	using DispatchFunc = void (*)(ContactConstraintBase &, MotionProperties *, MotionProperties *, float);
	static const DispatchFunc table[3][3] = {
		{ nullptr, nullptr, sIntegrateLegacyFrictionState<EMotionType::Static, EMotionType::Dynamic> },
		{ nullptr, nullptr, sIntegrateLegacyFrictionState<EMotionType::Kinematic, EMotionType::Dynamic> },
		{
			sIntegrateLegacyFrictionState<EMotionType::Dynamic, EMotionType::Static>,
			sIntegrateLegacyFrictionState<EMotionType::Dynamic, EMotionType::Kinematic>,
			sIntegrateLegacyFrictionState<EMotionType::Dynamic, EMotionType::Dynamic>
		}
	};
	for (const uint32 *offset = inConstraintOffsetBegin; offset < inConstraintOffsetEnd; ++offset)
	{
		ContactConstraintBase &constraint = *reinterpret_cast<ContactConstraintBase *>(mConstraints + *offset);
		Body &body1 = *constraint.mBody1;
		Body &body2 = *constraint.mBody2;
		table[(int)body1.GetMotionType()][(int)body2.GetMotionType()](
			constraint, body1.GetMotionPropertiesUnchecked(), body2.GetMotionPropertiesUnchecked(), inDeltaTime);
	}
}

bool ContactConstraintManager::SolveNonPenetrationConstraints(const uint32 *inConstraintOffsetBegin, const uint32 *inConstraintOffsetEnd)
{
	return TemplatedSolveVelocityConstraints<false, true, false>(inConstraintOffsetBegin, inConstraintOffsetEnd);
}

bool ContactConstraintManager::SolvePersistentNonPenetrationConstraints(const uint32 *inConstraintOffsetBegin, const uint32 *inConstraintOffsetEnd)
{
	return TemplatedSolveVelocityConstraints<false, true, true>(inConstraintOffsetBegin, inConstraintOffsetEnd);
}

template <EMotionType Type1, EMotionType Type2, bool ResetFriction, bool ResetNonPenetration>
void ContactConstraintManager::sResetVelocityConstraintLambdas(ContactConstraintBase &ioConstraint)
{
	ContactConstraint<Type1, Type2> &constraint = static_cast<ContactConstraint<Type1, Type2> &>(ioConstraint);
	for (uint32 i = 0; i < constraint.mNumContactPoints; ++i)
	{
		WorldContactPoint<Type1, Type2> &wcp = constraint.mContactPoints[i];
		if constexpr (ResetFriction)
		{
			wcp.mFrictionConstraint1.SetTotalLambda(0.0f);
			wcp.mFrictionConstraint2.SetTotalLambda(0.0f);
		}
		if constexpr (ResetNonPenetration)
			wcp.mNonPenetrationConstraint.SetTotalLambda(0.0f);
	}
}

template <EMotionType Type1, EMotionType Type2>
void ContactConstraintManager::sResetLegacyFrictionState(ContactConstraintBase &ioConstraint)
{
	ContactConstraint<Type1, Type2> &constraint = static_cast<ContactConstraint<Type1, Type2> &>(ioConstraint);
	if (!constraint.mContactImpacting)
		return;
	for (uint32 i = 0; i < constraint.mNumContactPoints; ++i)
	{
		WorldContactPoint<Type1, Type2> &wcp = constraint.mContactPoints[i];
		wcp.mFrictionConstraint1.SetTotalLambda(0.0f);
		wcp.mFrictionConstraint2.SetTotalLambda(0.0f);
		Vec3::sZero().StoreFloat3(&wcp.mContactPoint->mVPhysicsFrictionDisplacement);
		wcp.mContactPoint->mVPhysicsFrictionNormalLambda = 0.0f;
	}
}

template <bool ResetFriction, bool ResetNonPenetration>
void ContactConstraintManager::TemplatedResetVelocityConstraintLambdas(const uint32 *inConstraintOffsetBegin, const uint32 *inConstraintOffsetEnd)
{
	using DispatchFunc = void (*)(ContactConstraintBase &);
	static const DispatchFunc table[3][3] = {
		{
			nullptr,
			nullptr,
			sResetVelocityConstraintLambdas<EMotionType::Static, EMotionType::Dynamic, ResetFriction, ResetNonPenetration>
		},
		{
			nullptr,
			nullptr,
			sResetVelocityConstraintLambdas<EMotionType::Kinematic, EMotionType::Dynamic, ResetFriction, ResetNonPenetration>
		},
		{
			sResetVelocityConstraintLambdas<EMotionType::Dynamic, EMotionType::Static, ResetFriction, ResetNonPenetration>,
			sResetVelocityConstraintLambdas<EMotionType::Dynamic, EMotionType::Kinematic, ResetFriction, ResetNonPenetration>,
			sResetVelocityConstraintLambdas<EMotionType::Dynamic, EMotionType::Dynamic, ResetFriction, ResetNonPenetration>
		}
	};

	for (const uint32 *offset = inConstraintOffsetBegin; offset < inConstraintOffsetEnd; ++offset)
	{
		ContactConstraintBase &constraint = *reinterpret_cast<ContactConstraintBase *>(mConstraints + *offset);
		table[(int)constraint.mBody1->GetMotionType()][(int)constraint.mBody2->GetMotionType()](constraint);
	}
}

void ContactConstraintManager::ResetFrictionLambdas(const uint32 *inConstraintOffsetBegin, const uint32 *inConstraintOffsetEnd)
{
	TemplatedResetVelocityConstraintLambdas<true, false>(inConstraintOffsetBegin, inConstraintOffsetEnd);
}

void ContactConstraintManager::ResetNonPenetrationLambdas(const uint32 *inConstraintOffsetBegin, const uint32 *inConstraintOffsetEnd)
{
	TemplatedResetVelocityConstraintLambdas<false, true>(inConstraintOffsetBegin, inConstraintOffsetEnd);
}

void ContactConstraintManager::ResetImpactingLegacyFrictionState(const uint32 *inConstraintOffsetBegin, const uint32 *inConstraintOffsetEnd)
{
	using DispatchFunc = void (*)(ContactConstraintBase &);
	static const DispatchFunc table[3][3] = {
		{ nullptr, nullptr, sResetLegacyFrictionState<EMotionType::Static, EMotionType::Dynamic> },
		{ nullptr, nullptr, sResetLegacyFrictionState<EMotionType::Kinematic, EMotionType::Dynamic> },
		{
			sResetLegacyFrictionState<EMotionType::Dynamic, EMotionType::Static>,
			sResetLegacyFrictionState<EMotionType::Dynamic, EMotionType::Kinematic>,
			sResetLegacyFrictionState<EMotionType::Dynamic, EMotionType::Dynamic>
		}
	};
	for (const uint32 *offset = inConstraintOffsetBegin; offset < inConstraintOffsetEnd; ++offset)
	{
		ContactConstraintBase &constraint = *reinterpret_cast<ContactConstraintBase *>(mConstraints + *offset);
		table[(int)constraint.mBody1->GetMotionType()][(int)constraint.mBody2->GetMotionType()](constraint);
	}
}

bool ContactConstraintManager::AreAnyContactsImpacting(const uint32 *inConstraintOffsetBegin, const uint32 *inConstraintOffsetEnd) const
{
	if (inConstraintOffsetBegin >= inConstraintOffsetEnd)
		return false;

	for (const uint32 *offset = inConstraintOffsetBegin; offset < inConstraintOffsetEnd; ++offset)
	{
		const ContactConstraintBase &constraint = *reinterpret_cast<const ContactConstraintBase *>(mConstraints + *offset);
		if (constraint.mContactImpacting)
			return true;
	}
	return false;
}

bool ContactConstraintManager::AreAnyContactsUsingVPhysicsSolver(const uint32 *inConstraintOffsetBegin, const uint32 *inConstraintOffsetEnd) const
{
	if (inConstraintOffsetBegin >= inConstraintOffsetEnd)
		return false;

	for (const uint32 *offset = inConstraintOffsetBegin; offset < inConstraintOffsetEnd; ++offset)
	{
		const ContactConstraintBase &constraint = *reinterpret_cast<const ContactConstraintBase *>(mConstraints + *offset);
		if (constraint.mUseVPhysicsSolver)
			return true;
	}
	return false;
}

bool ContactConstraintManager::AreAnyContactIndicesUsingVPhysicsSolver(const uint32 *inConstraintIdxBegin, const uint32 *inConstraintIdxEnd) const
{
	if (inConstraintIdxBegin >= inConstraintIdxEnd)
		return false;

	// IslandBuilder still stores contact indices when LargeIslandSplitter::Prepare runs.
	// The per-island solve converts these entries to byte offsets later, immediately before
	// sorting and solving. Inspect through the immutable index-to-offset table here so Prepare
	// can count only islands that will actually call SplitIsland.
	for (const uint32 *idx = inConstraintIdxBegin; idx < inConstraintIdxEnd; ++idx)
	{
		const uint32 offset = mConstraintIdxToOffset[*idx];
		const ContactConstraintBase &constraint = *reinterpret_cast<const ContactConstraintBase *>(mConstraints + offset);
		if (constraint.mUseVPhysicsSolver)
			return true;
	}
	return false;
}

template <EMotionType Type1, EMotionType Type2>
float ContactConstraintManager::sStoreAppliedImpulses(ContactConstraintBase &ioConstraint)
{
	ContactConstraint<Type1, Type2> &constraint = static_cast<ContactConstraint<Type1, Type2> &>(ioConstraint);
	float normal_impulse = 0.0f;
	const bool queue_native_impact_correction = constraint.mContactImpacting
		&& constraint.mVPhysicsNativeImpactCorrection
		&& !constraint.mUseVPhysicsSolver
		&& !constraint.mSevereInterpenetration
		&& !constraint.mPlayerSupportStartup
		&& constraint.mVPhysicsNativeImpactTime > 0.0f;
	Vec3 position_correction1 = Vec3::sZero();
	Vec3 rotation_correction1 = Vec3::sZero();
	Vec3 position_correction2 = Vec3::sZero();
	Vec3 rotation_correction2 = Vec3::sZero();
	Vec3 normal = Vec3::sZero();
	Vec3 tangent1 = Vec3::sZero();
	Vec3 tangent2 = Vec3::sZero();
	RMat44 transform1;
	RMat44 transform2;
	Mat44 inv_i1 = Mat44::sZero();
	Mat44 inv_i2 = Mat44::sZero();
	if (queue_native_impact_correction)
	{
		normal = constraint.GetWorldSpaceNormal();
		constraint.GetTangents(tangent1, tangent2);
		transform1 = constraint.mBody1->GetCenterOfMassTransform();
		transform2 = constraint.mBody2->GetCenterOfMassTransform();
		if constexpr (Type1 == EMotionType::Dynamic)
			inv_i1 = constraint.mInvInertiaScale1 * constraint.mBody1->GetInverseInertia();
		if constexpr (Type2 == EMotionType::Dynamic)
			inv_i2 = constraint.mInvInertiaScale2 * constraint.mBody2->GetInverseInertia();
	}

	for (uint32 i = 0; i < constraint.mNumContactPoints; ++i)
	{
		WorldContactPoint<Type1, Type2> &wcp = constraint.mContactPoints[i];

		const float point_normal_impulse = wcp.mNonPenetrationConstraint.GetTotalLambda();
		if (queue_native_impact_correction)
		{
			const Vec3 point_impulse = normal * point_normal_impulse
				+ tangent1 * wcp.mFrictionConstraint1.GetTotalLambda()
				+ tangent2 * wcp.mFrictionConstraint2.GetTotalLambda();
			const RVec3 point1 = transform1
				* Vec3::sLoadFloat3Unsafe(wcp.mContactPoint->mPosition1);
			const RVec3 point2 = transform2
				* Vec3::sLoadFloat3Unsafe(wcp.mContactPoint->mPosition2);
			const RVec3 point = 0.5_r * (point1 + point2);
			if constexpr (Type1 == EMotionType::Dynamic)
			{
				const Vec3 r1 = Vec3(point - constraint.mBody1->GetCenterOfMassPosition());
				position_correction1 += constraint.mInvMass1 * point_impulse;
				rotation_correction1 += inv_i1.Multiply3x3(r1.Cross(point_impulse));
			}
			if constexpr (Type2 == EMotionType::Dynamic)
			{
				const Vec3 r2 = Vec3(point - constraint.mBody2->GetCenterOfMassPosition());
				position_correction2 -= constraint.mInvMass2 * point_impulse;
				rotation_correction2 -= inv_i2.Multiply3x3(r2.Cross(point_impulse));
			}
		}
		// Preserve the exact solver warm-start state for the next PSI. Whether this lambda
		// represents an impact or a distance-keeper pressure is decided by the listener/cache
		// below; changing these values here alters the actual collision solver and can lock a
		// player on slopes.
		wcp.mContactPoint->mNonPenetrationLambda = point_normal_impulse;
		wcp.mContactPoint->mFrictionLambda[0] = wcp.mFrictionConstraint1.GetTotalLambda();
		wcp.mContactPoint->mFrictionLambda[1] = wcp.mFrictionConstraint2.GetTotalLambda();
		wcp.mContactPoint->mVPhysicsFrictionNormalLambda = point_normal_impulse;
		wcp.mContactPoint->mVPhysicsContactEstablished =
			wcp.mContactPoint->mVPhysicsContactEstablished || point_normal_impulse > 0.0f;
		// The listener's normal-pressure value is the distance-keeper contribution. The
		// one-shot impact lambda is still retained in the solver cache above, but is not
		// published as continuous Source friction pressure.
		if ( !constraint.mContactImpacting )
			normal_impulse += point_normal_impulse;
	}
	if (queue_native_impact_correction)
	{
		const float impact_time = constraint.mVPhysicsNativeImpactTime;
		if constexpr (Type1 == EMotionType::Dynamic)
			constraint.mBody1->GetMotionPropertiesUnchecked()->AddVPhysicsIntegrationCorrection(
				position_correction1 * impact_time, rotation_correction1 * impact_time);
		if constexpr (Type2 == EMotionType::Dynamic)
			constraint.mBody2->GetMotionPropertiesUnchecked()->AddVPhysicsIntegrationCorrection(
				position_correction2 * impact_time, rotation_correction2 * impact_time);
	}
	if (constraint.mPlayerSupportStartup && constraint.mContactImpacting
		&& !constraint.mPlayerSupportStartupCorrectionQueued
		&& constraint.mPlayerSupportStartupDeltaTime > 0.0f)
	{
		const Vec3 normal = constraint.GetWorldSpaceNormal();
		if constexpr (Type1 == EMotionType::Dynamic)
		{
			MotionProperties *motion_properties1 = constraint.mBody1->GetMotionPropertiesUnchecked();
			const float current_normal_velocity = motion_properties1->GetLinearVelocity().Dot(normal);
			motion_properties1->AddVPhysicsIntegrationCorrection(
				normal * ((constraint.mPlayerSupportStartupInitialNormalVelocity1
					- current_normal_velocity) * constraint.mPlayerSupportStartupDeltaTime),
				Vec3::sZero());
		}
		if constexpr (Type2 == EMotionType::Dynamic)
		{
			MotionProperties *motion_properties2 = constraint.mBody2->GetMotionPropertiesUnchecked();
			const float current_normal_velocity = motion_properties2->GetLinearVelocity().Dot(normal);
			motion_properties2->AddVPhysicsIntegrationCorrection(
				normal * ((constraint.mPlayerSupportStartupInitialNormalVelocity2
					- current_normal_velocity) * constraint.mPlayerSupportStartupDeltaTime),
				Vec3::sZero());
		}
		constraint.mPlayerSupportStartupCorrectionQueued = true;
	}
	return normal_impulse;
}

void ContactConstraintManager::StoreAppliedImpulses(const uint32 *inConstraintOffsetBegin, const uint32 *inConstraintOffsetEnd) const
{
	// Build dispatch table
	using DispatchFunc = float (*)(ContactConstraintBase &);
	static const DispatchFunc table[3][3] = {
		{
			nullptr, // Static vs static doesn't exist
			nullptr, // Static vs kinematic doesn't exist
			sStoreAppliedImpulses<EMotionType::Static, EMotionType::Dynamic>
		},
		{
			nullptr, // Kinematic vs static doesn't exist
			nullptr, // Kinematic vs kinematic doesn't exist
			sStoreAppliedImpulses<EMotionType::Kinematic, EMotionType::Dynamic>
		},
		{
			sStoreAppliedImpulses<EMotionType::Dynamic, EMotionType::Static>,
			sStoreAppliedImpulses<EMotionType::Dynamic, EMotionType::Kinematic>,
			sStoreAppliedImpulses<EMotionType::Dynamic, EMotionType::Dynamic>
		}
	};

	if (inConstraintOffsetBegin >= inConstraintOffsetEnd)
		return;

	// Copy back total applied impulse to cache for the next frame
	ContactConstraintBase *next_constraint = reinterpret_cast<ContactConstraintBase *>(mConstraints + *inConstraintOffsetBegin);
	for (const uint32 *next_constraint_offset = inConstraintOffsetBegin + 1; next_constraint != nullptr; ++next_constraint_offset)
	{
		ContactConstraintBase &constraint = *next_constraint;
		if (next_constraint_offset < inConstraintOffsetEnd)
		{
			next_constraint = reinterpret_cast<ContactConstraintBase *>(mConstraints + *next_constraint_offset);
			PrefetchL1(next_constraint);
		}
		else
			next_constraint = nullptr;

		// Dispatch to the correct templated form
		const bool was_impact = constraint.mContactImpacting;
		const float normal_impulse = table[(int)constraint.mBody1->GetMotionType()][(int)constraint.mBody2->GetMotionType()](constraint);
		if (mContactListener != nullptr)
			mContactListener->OnContactSolved(constraint.mSubShapePair, normal_impulse, was_impact);
	}
}

template <EMotionType Type1, EMotionType Type2>
bool ContactConstraintManager::sSolvePositionConstraint(ContactConstraintBase &ioConstraint, Body &ioBody1, Body &ioBody2, const PhysicsSettings &inSettings)
{
	ContactConstraint<Type1, Type2> &constraint = static_cast<ContactConstraint<Type1, Type2> &>(ioConstraint);
	// IVP resolves a fresh impact in the velocity domain at the event time. Running
	// Jolt's split-position correction for the same contact afterwards moves the body
	// a second time without removing kinetic energy, which can make a dropped object
	// rebound above its release height. Persistent contacts may still use the ordinary
	// position stabilizer; only the frame that consumed the impact skips it.
	if (constraint.mDisablePenetrationRecovery || constraint.mContactImpacting)
		return false;

	// Get transforms
	RMat44 transform1 = ioBody1.GetCenterOfMassTransform();
	RMat44 transform2 = ioBody2.GetCenterOfMassTransform();

	Vec3 ws_normal = constraint.GetWorldSpaceNormal();

	bool any_impulse_applied = false;

	for (uint32 i = 0; i < constraint.mNumContactPoints; ++i)
	{
		WorldContactPoint<Type1, Type2> &wcp = constraint.mContactPoints[i];

		// Calculate new contact point positions in world space (the bodies may have moved)
		RVec3 p1 = transform1 * Vec3::sLoadFloat3Unsafe(wcp.mContactPoint->mPosition1);
		RVec3 p2 = transform2 * Vec3::sLoadFloat3Unsafe(wcp.mContactPoint->mPosition2);

		// Calculate separation from the desired contact boundary (negative if too close)
		// Allow a little penetration by default (PhysicsSettings::mPenetrationSlop) to avoid jittering between contact/no-contact which wipes out the contact cache and warm start impulses
		// Clamp penetration to a max PhysicsSettings::mMaxPenetrationDistance so that we don't apply a huge impulse if we're penetrating a lot
		// IVP's outer friction distance is maintained by its velocity-level distance keeper. Its
		// hard positional boundary remains the inner collision distance (one collision tolerance),
		// regardless of whether the contact already carries support pressure.
		const float separation_target = 0.5f * inSettings.mContactSeparationTarget;
		float separation = max(Vec3(p2 - p1).Dot(ws_normal) - separation_target + inSettings.mPenetrationSlop, -inSettings.mMaxPenetrationDistance);

		// Only enforce constraint when separation < 0 (otherwise we're apart)
		if (separation < 0.0f)
		{
			// Calculate scaled inertia
			Mat44 inv_i1;
			if constexpr (Type1 == EMotionType::Dynamic)
				inv_i1 = constraint.mInvInertiaScale1 * ioBody1.GetInverseInertia();
			else
				inv_i1 = Mat44::sZero();

			Mat44 inv_i2;
			if constexpr (Type2 == EMotionType::Dynamic)
				inv_i2 = constraint.mInvInertiaScale2 * ioBody2.GetInverseInertia();
			else
				inv_i2 = Mat44::sZero();

			// Calculate collision points relative to body
			RVec3 p = 0.5_r * (p1 + p2);
			Vec3 r1 = Vec3(p - ioBody1.GetCenterOfMassPosition());
			Vec3 r2 = Vec3(p - ioBody2.GetCenterOfMassPosition());

			// Update constraint properties (bodies may have moved)
			wcp.mNonPenetrationConstraint.CalculateConstraintProperties(constraint.mInvMass1, inv_i1, r1, constraint.mInvMass2, inv_i2, r2, ws_normal);

			// Solve position errors
			if (wcp.mNonPenetrationConstraint.SolvePositionConstraint(ioBody1, constraint.mInvMass1, ioBody2, constraint.mInvMass2, ws_normal, separation, inSettings.mBaumgarte))
				any_impulse_applied = true;
		}
	}

	return any_impulse_applied;
}

bool ContactConstraintManager::SolvePositionConstraints(const uint32 *inConstraintOffsetBegin, const uint32 *inConstraintOffsetEnd)
{
	JPH_PROFILE_FUNCTION();

	// Build dispatch table
	using DispatchFunc = bool (*)(ContactConstraintBase &, Body &, Body &, const PhysicsSettings &);
	static const DispatchFunc table[3][3] = {
		{
			nullptr, // Static vs static doesn't exist
			nullptr, // Static vs kinematic doesn't exist
			sSolvePositionConstraint<EMotionType::Static, EMotionType::Dynamic>
		},
		{
			nullptr, // Kinematic vs static doesn't exist
			nullptr, // Kinematic vs kinematic doesn't exist
			sSolvePositionConstraint<EMotionType::Kinematic, EMotionType::Dynamic>
		},
		{
			sSolvePositionConstraint<EMotionType::Dynamic, EMotionType::Static>,
			sSolvePositionConstraint<EMotionType::Dynamic, EMotionType::Kinematic>,
			sSolvePositionConstraint<EMotionType::Dynamic, EMotionType::Dynamic>
		}
	};

	if (inConstraintOffsetBegin >= inConstraintOffsetEnd)
		return false;

	bool any_impulse_applied = false;

	ContactConstraintBase *next_constraint = reinterpret_cast<ContactConstraintBase *>(mConstraints + *inConstraintOffsetBegin);
	for (const uint32 *next_constraint_offset = inConstraintOffsetBegin + 1; next_constraint != nullptr; ++next_constraint_offset)
	{
		ContactConstraintBase &constraint = *next_constraint;
		if (next_constraint_offset < inConstraintOffsetEnd)
		{
			next_constraint = reinterpret_cast<ContactConstraintBase *>(mConstraints + *next_constraint_offset);
			PrefetchL1(next_constraint);
		}
		else
			next_constraint = nullptr;

		// Fetch bodies
		Body &body1 = *constraint.mBody1;
		Body &body2 = *constraint.mBody2;

		// Dispatch to the correct templated form
		any_impulse_applied |= table[(int)body1.GetMotionType()][(int)body2.GetMotionType()](constraint, body1, body2, mPhysicsSettings);
	}

	return any_impulse_applied;
}

void ContactConstraintManager::RecycleConstraintBuffer()
{
	// Reset constraint array
	mNumConstraintsAndNextConstraintOffset = 0;
}

void ContactConstraintManager::FinishConstraintBuffer()
{
	// Free constraints buffer
	mUpdateContext->mTempAllocator->Free(mConstraintIdxToOffset, mMaxConstraints * sizeof(uint32));
	mConstraintIdxToOffset = nullptr;
	mUpdateContext->mTempAllocator->Free(mConstraints, mMaxConstraints * cMaxConstraintSize);
	mConstraints = nullptr;
	mNumConstraintsAndNextConstraintOffset = 0;

	// Reset update context
	mUpdateContext = nullptr;
}

void ContactConstraintManager::SaveState(StateRecorder &inStream, const StateRecorderFilter *inFilter) const
{
	mCache[mCacheWriteIdx ^ 1].SaveState(inStream, inFilter);
}

bool ContactConstraintManager::RestoreState(StateRecorder &inStream, const StateRecorderFilter *inFilter)
{
	bool success = mCache[mCacheWriteIdx].RestoreState(mCache[mCacheWriteIdx ^ 1], inStream, inFilter);

	// If this is the last part, the cache is finalized
	if (inStream.IsLastPart())
	{
		mCacheWriteIdx ^= 1;
		mCache[mCacheWriteIdx].Clear();
	}

	return success;
}

JPH_NAMESPACE_END
