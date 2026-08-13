// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-FileCopyrightText: 2021 Jorrit Rouwe
// SPDX-License-Identifier: MIT

#include <Jolt/Jolt.h>

#include <Jolt/Physics/Constraints/Constraint.h>
#include <Jolt/Physics/Constraints/TwoBodyConstraint.h>
#include <Jolt/Physics/Constraints/LegacyBreakableState.h>
#include <Jolt/Physics/StateRecorder.h>
#include <Jolt/ObjectStream/TypeDeclarations.h>
#include <Jolt/Core/StreamUtils.h>

#include <algorithm>

JPH_NAMESPACE_BEGIN

JPH_IMPLEMENT_SERIALIZABLE_VIRTUAL(ConstraintSettings)
{
	JPH_ADD_BASE_CLASS(ConstraintSettings, SerializableObject)

	JPH_ADD_ATTRIBUTE(ConstraintSettings, mEnabled)
	JPH_ADD_ATTRIBUTE(ConstraintSettings, mDrawConstraintSize)
	JPH_ADD_ATTRIBUTE(ConstraintSettings, mConstraintPriority)
	JPH_ADD_ATTRIBUTE(ConstraintSettings, mNumVelocityStepsOverride)
	JPH_ADD_ATTRIBUTE(ConstraintSettings, mNumPositionStepsOverride)
	JPH_ADD_ATTRIBUTE(ConstraintSettings, mUserData)
}

void ConstraintSettings::SaveBinaryState(StreamOut &inStream) const
{
	inStream.Write(GetRTTI()->GetHash());
	inStream.Write(mEnabled);
	inStream.Write(mDrawConstraintSize);
	inStream.Write(mConstraintPriority);
	inStream.Write(mNumVelocityStepsOverride);
	inStream.Write(mNumPositionStepsOverride);
}

void ConstraintSettings::RestoreBinaryState(StreamIn &inStream)
{
	// Type hash read by sRestoreFromBinaryState
	inStream.Read(mEnabled);
	inStream.Read(mDrawConstraintSize);
	inStream.Read(mConstraintPriority);
	inStream.Read(mNumVelocityStepsOverride);
	inStream.Read(mNumPositionStepsOverride);
}

ConstraintSettings::ConstraintResult ConstraintSettings::sRestoreFromBinaryState(StreamIn &inStream)
{
	return StreamUtils::RestoreObject<ConstraintSettings>(inStream, &ConstraintSettings::RestoreBinaryState);
}

void Constraint::SaveState(StateRecorder &inStream) const
{
	inStream.Write(mEnabled);
}

void Constraint::RestoreState(StateRecorder &inStream)
{
	inStream.Read(mEnabled);
	mLegacyBroken.store(false, std::memory_order_release);
}

void Constraint::ToConstraintSettings(ConstraintSettings &outSettings) const
{
	outSettings.mEnabled = mEnabled;
	outSettings.mConstraintPriority = mConstraintPriority;
	outSettings.mNumVelocityStepsOverride = mNumVelocityStepsOverride;
	outSettings.mNumPositionStepsOverride = mNumPositionStepsOverride;
	outSettings.mUserData = mUserData;
#ifdef JPH_DEBUG_RENDERER
	outSettings.mDrawConstraintSize = mDrawConstraintSize;
#endif // JPH_DEBUG_RENDERER
}

bool Constraint::CheckLegacyBreakAfterSolve(Vec3Arg inLinearVelocity1, Vec3Arg inAngularVelocity1, Vec3Arg inLinearVelocity2, Vec3Arg inAngularVelocity2, bool inRollbackOnBreak)
{
	if (!HasLegacyBreakLimits() || GetType() != EConstraintType::TwoBodyConstraint
		|| mLegacyBroken.load(std::memory_order_acquire))
		return false;

	auto body_breaks = [&](const Body &inBody, Vec3Arg inOldLinearVelocity, Vec3Arg inOldAngularVelocity)
	{
		if (!inBody.IsDynamic())
			return false;

		const MotionProperties *motion = inBody.GetMotionPropertiesUnchecked();
		if (mLegacyLinearBreakLimit > 0.0f
			&& LegacyBreaksLinear(
				1.0f / motion->GetInverseMassUnchecked(),
				inBody.GetLinearVelocity() - inOldLinearVelocity,
				mLegacyLinearBreakLimit))
			return true;

		return mLegacyAngularBreakLimit > 0.0f
			&& LegacyBreaksAngular(
				inBody.GetRotation(),
				motion->GetInertiaRotation(),
				motion->GetInverseInertiaDiagonal(),
				inBody.GetAngularVelocity() - inOldAngularVelocity,
				mLegacyAngularBreakLimit);
	};

	const TwoBodyConstraint *constraint = static_cast<const TwoBodyConstraint *>(this);
	Body *body1 = constraint->GetBody1();
	Body *body2 = constraint->GetBody2();
	if (body_breaks(*body1, inLinearVelocity1, inAngularVelocity1)
		|| body_breaks(*body2, inLinearVelocity2, inAngularVelocity2))
	{
		// A standalone Source constraint owns a one-member local system. Its synchronous game
		// callback deactivates that system, so the wrapper restores both bodies to the velocity
		// snapshot from immediately before this pass. Shared LCS breaking is explicitly marked as
		// unsupported in Valve's adapter and is kept as a separate compatibility item.
		if (inRollbackOnBreak && !body1->IsStatic())
		{
			body1->SetLinearVelocity(inLinearVelocity1);
			body1->SetAngularVelocity(inAngularVelocity1);
		}
		if (inRollbackOnBreak && !body2->IsStatic())
		{
			body2->SetLinearVelocity(inLinearVelocity2);
			body2->SetAngularVelocity(inAngularVelocity2);
		}
		mLegacyBroken.store(true, std::memory_order_release);
		return true;
	}

	return false;
}

LegacyConstraintGroup::~LegacyConstraintGroup()
{
	// Havana leaves member constraints alive when their local constraint system is deleted,
	// but clears the back-pointer so their later destruction is safe. The deleted system no
	// longer drives them, so keep the Jolt members disabled instead of letting them silently
	// fall back to the ordinary per-constraint solver.
	for (Constraint *constraint : mConstraints)
		if (constraint->GetLegacyConstraintGroup() == this)
		{
			constraint->SetEnabled(false);
			constraint->SetLegacyConstraintGroup(nullptr);
		}
}

void LegacyConstraintGroup::AddConstraint(Constraint *inConstraint)
{
	if (std::find(mConstraints.begin(), mConstraints.end(), inConstraint) == mConstraints.end())
		mConstraints.push_back(inConstraint);
	inConstraint->SetLegacyConstraintGroup(this);
}

void LegacyConstraintGroup::RemoveConstraint(Constraint *inConstraint)
{
	if (inConstraint->GetLegacyConstraintGroup() == this)
		inConstraint->SetLegacyConstraintGroup(nullptr);
	std::erase(mConstraints, inConstraint);
}

bool LegacyConstraintGroup::IsSolverRepresentative(const Constraint *inCaller) const
{
	for (Constraint *constraint : mConstraints)
		if (constraint->IsActive())
			return constraint == inCaller;

	return false;
}

bool LegacyConstraintGroup::SolveVelocityConstraints(Constraint *inCaller, float inDeltaTime)
{
	if (!IsSolverRepresentative(inCaller))
		return false;
	if (!mSolvePending.exchange(false, std::memory_order_acq_rel))
		return false;

	// ConstraintManager performs the ordinary setup before contact/friction solving. Source local
	// constraint systems are controllers that run afterwards, so their predictive limits must see
	// the velocities produced by those contacts. Rebuild every member's velocity state here; merely
	// postponing SolveVelocityConstraint leaves a one-frame hole whenever a previously inactive
	// limit is crossed by a contact impulse (notably a simfphys wheel striking a curb).
	for (Constraint *constraint : mConstraints)
		if (constraint->IsActive())
			constraint->SetupVelocityConstraint(inDeltaTime);

	bool appliedImpulse = false;
	bool groupBroken = false;
	auto solveConstraint = [&](Constraint *constraint)
	{
		if (groupBroken)
			return;
		const TwoBodyConstraint *two_body = constraint->GetType() == EConstraintType::TwoBodyConstraint?
			static_cast<const TwoBodyConstraint *>(constraint) : nullptr;
		const Vec3 linear1 = two_body? two_body->GetBody1()->GetLinearVelocity() : Vec3::sZero();
		const Vec3 angular1 = two_body? two_body->GetBody1()->GetAngularVelocity() : Vec3::sZero();
		const Vec3 linear2 = two_body? two_body->GetBody2()->GetLinearVelocity() : Vec3::sZero();
		const Vec3 angular2 = two_body? two_body->GetBody2()->GetAngularVelocity() : Vec3::sZero();
		appliedImpulse |= constraint->SolveVelocityConstraint(inDeltaTime);
		if (constraint->CheckLegacyBreakAfterSolve(linear1, angular1, linear2, angular2, mRollbackOnBreak)
			&& mRollbackOnBreak)
		{
			constraint->SetEnabled(false);
			groupBroken = true;
		}
	};
	auto solveForward = [&]()
	{
		for (Constraint *constraint : mConstraints)
			if (constraint->IsActive())
				solveConstraint(constraint);
	};
	auto solveReverse = [&]()
	{
		for (auto it = mConstraints.rbegin(); it != mConstraints.rend(); ++it)
			if ((*it)->IsActive())
				solveConstraint(*it);
	};

	// Source performs one setup solve, followed by m_n_iterations passes;
	// m_n_iterations is additionalIterations + 2 and each pass has two
	// reverse/forward stages.
	solveForward();
	// Havana stores m_n_iterations as signed additionalIterations + 2. Negative input is not
	// clamped: -1 performs one iteration, while <= -2 performs none after the setup solve.
	const int64 num_iterations = int64(mAdditionalIterations) + 2;
	for (int64 iteration = 0; iteration < num_iterations && !groupBroken; ++iteration)
	{
		solveReverse();
		solveForward();
		solveReverse();
		solveForward();
	}

	(void)inCaller;
	return appliedImpulse;
}

JPH_NAMESPACE_END
