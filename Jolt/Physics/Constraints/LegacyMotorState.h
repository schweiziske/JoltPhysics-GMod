// SPDX-License-Identifier: MIT

#pragma once

#include <Jolt/Jolt.h>
#include <Jolt/Math/Math.h>
#include <Jolt/Physics/StateRecorder.h>

#include <cmath>

JPH_NAMESPACE_BEGIN

struct LegacyMotorStep
{
	float mBias = 0.0f;
	float mRelaxation = 1.0f;
	float mMinLambda = 0.0f;
	float mMaxLambda = 0.0f;
	bool mEnabled = false;
};

/// Persistent reference-position controller used by Havana hinge and prismatic constraints.
/// This is intentionally separate from Jolt's normal velocity motor: Havana advances the
/// reference before applying an independently clamped impulse and rolls the unsatisfied part
/// back when the force limit saturates.
class LegacyMotorState
{
public:
	void Configure(float inDesiredVelocity, float inMaxForce, float inEpsilon = 1.0f)
	{
		mDesiredVelocity = inDesiredVelocity * inEpsilon;
		mMaxForce = std::abs(inMaxForce * inEpsilon);
		mEnabled = mMaxForce != 0.0f;
	}

	bool HasReferencePosition() const
	{
		return mReferenceInitialized;
	}

	void InitializeLinearReference(float inPosition)
	{
		mReferencePosition = inPosition;
		mReferenceInitialized = true;
	}

	void InitializeAngularReference(float inAngle)
	{
		if (!mReferenceInitialized)
		{
			// hk_Constraint_Limit constructs m_ref_position at zero. The initialization call only
			// unwraps an already existing reference; it does not adopt the current hinge angle.
			mReferencePosition = 0.0f;
			mReferenceInitialized = true;
		}

		if (inAngle < 0.0f && mReferencePosition > 0.0f && mReferencePosition - inAngle > JPH_PI)
			mReferencePosition = inAngle + 2.0f * JPH_PI;
		else if (inAngle > 0.0f && mReferencePosition < 0.0f && inAngle - mReferencePosition > JPH_PI)
			mReferencePosition = inAngle - 2.0f * JPH_PI;
	}

	LegacyMotorStep BeginLinearStep(float inPosition, float inDeltaTime, float inTauFactor = 1.0f, float inDampingFactor = 1.0f)
	{
		if (!mReferenceInitialized)
			InitializeLinearReference(inPosition);

		return BeginStep(inPosition - mReferencePosition, inDeltaTime, inTauFactor, inDampingFactor, false);
	}

	LegacyMotorStep BeginAngularStep(float inAngle, float inDeltaTime, float inTauFactor = 1.0f, float inDampingFactor = 1.0f)
	{
		if (!mReferenceInitialized)
			InitializeAngularReference(inAngle);

		float position_error = mReferencePosition - inAngle;
		if (position_error > JPH_PI)
			position_error -= 2.0f * JPH_PI;
		else if (position_error < -JPH_PI)
			position_error += 2.0f * JPH_PI;

		return BeginStep(position_error, inDeltaTime, inTauFactor, inDampingFactor, true);
	}

	void FinishStep(float inRequestedLambda, float inAppliedLambda)
	{
		if (!mStepPending)
			return;

		const float requested_magnitude = std::abs(inRequestedLambda);
		if (requested_magnitude > 0.0f && std::abs(inAppliedLambda) < requested_magnitude)
		{
			const float factor = Clamp(std::abs(inAppliedLambda) / requested_magnitude, 0.0f, 1.0f);
			if (mAngularStep)
				mReferencePosition -= (1.0f - factor) * mPositionError;
			else
				mReferencePosition += (1.0f - factor) * mPositionError;
		}

		mStepPending = false;
	}

	void SetReferencePosition(float inReferencePosition)
	{
		mReferencePosition = inReferencePosition;
		mReferenceInitialized = true;
	}

	float GetReferencePosition() const { return mReferencePosition; }
	float GetDesiredVelocity() const { return mDesiredVelocity; }
	float GetMaxForce() const { return mMaxForce; }
	bool IsEnabled() const { return mEnabled; }

	void SaveState(StateRecorder &inStream) const
	{
		inStream.Write(mReferencePosition);
		inStream.Write(mDesiredVelocity);
		inStream.Write(mMaxForce);
		inStream.Write(mReferenceInitialized);
		inStream.Write(mEnabled);
	}

	void RestoreState(StateRecorder &inStream)
	{
		inStream.Read(mReferencePosition);
		inStream.Read(mDesiredVelocity);
		inStream.Read(mMaxForce);
		inStream.Read(mReferenceInitialized);
		inStream.Read(mEnabled);
		mStepPending = false;
	}

private:
	LegacyMotorStep BeginStep(float inPositionError, float inDeltaTime, float inTauFactor, float inDampingFactor, bool inAngular)
	{
		LegacyMotorStep result;
		if (!mEnabled || inDeltaTime <= 0.0f)
			return result;

		mPositionError = inPositionError;
		mReferencePosition += mDesiredVelocity * inDeltaTime;
		mAngularStep = inAngular;
		mStepPending = true;

		result.mBias = inPositionError * 0.8f * inTauFactor / inDeltaTime;
		result.mRelaxation = inDampingFactor;
		result.mMinLambda = -mMaxForce * inDeltaTime;
		result.mMaxLambda = mMaxForce * inDeltaTime;
		result.mEnabled = true;
		return result;
	}

	float mReferencePosition = 0.0f;
	float mDesiredVelocity = 0.0f;
	float mMaxForce = 0.0f;
	float mPositionError = 0.0f;
	bool mReferenceInitialized = false;
	bool mEnabled = false;
	bool mAngularStep = false;
	bool mStepPending = false;
};

JPH_NAMESPACE_END
