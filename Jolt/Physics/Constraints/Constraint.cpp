// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-FileCopyrightText: 2021 Jorrit Rouwe
// SPDX-License-Identifier: MIT

#include <Jolt/Jolt.h>

#include <Jolt/Physics/Constraints/Constraint.h>
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

bool LegacyConstraintGroup::SolveVelocityConstraints(Constraint *inCaller, float inDeltaTime)
{
	if (!mSolvePending.exchange(false, std::memory_order_acq_rel))
		return false;

	bool appliedImpulse = false;
	auto solveForward = [&]()
	{
		for (Constraint *constraint : mConstraints)
			if (constraint->IsActive())
				appliedImpulse |= constraint->SolveVelocityConstraint(inDeltaTime);
	};
	auto solveReverse = [&]()
	{
		for (auto it = mConstraints.rbegin(); it != mConstraints.rend(); ++it)
			if ((*it)->IsActive())
				appliedImpulse |= (*it)->SolveVelocityConstraint(inDeltaTime);
	};

	// Source performs one setup solve, followed by m_n_iterations passes;
	// m_n_iterations is additionalIterations + 2 and each pass has two
	// reverse/forward stages.
	solveForward();
	for (uint iteration = 0; iteration < mAdditionalIterations + 2; ++iteration)
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
