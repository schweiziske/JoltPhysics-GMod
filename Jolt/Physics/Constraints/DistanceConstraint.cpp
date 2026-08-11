// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-FileCopyrightText: 2021 Jorrit Rouwe
// SPDX-License-Identifier: MIT

#include <Jolt/Jolt.h>

#include <Jolt/Physics/Constraints/DistanceConstraint.h>
#include <Jolt/Physics/Body/Body.h>
#include <Jolt/ObjectStream/TypeDeclarations.h>
#include <Jolt/Core/StreamIn.h>
#include <Jolt/Core/StreamOut.h>
#ifdef JPH_DEBUG_RENDERER
	#include <Jolt/Renderer/DebugRenderer.h>
#endif // JPH_DEBUG_RENDERER

JPH_NAMESPACE_BEGIN

using namespace literals;

JPH_IMPLEMENT_SERIALIZABLE_VIRTUAL(DistanceConstraintSettings)
{
	JPH_ADD_BASE_CLASS(DistanceConstraintSettings, TwoBodyConstraintSettings)

	JPH_ADD_ENUM_ATTRIBUTE(DistanceConstraintSettings, mSpace)
	JPH_ADD_ATTRIBUTE(DistanceConstraintSettings, mPoint1)
	JPH_ADD_ATTRIBUTE(DistanceConstraintSettings, mPoint2)
	JPH_ADD_ATTRIBUTE(DistanceConstraintSettings, mMinDistance)
	JPH_ADD_ATTRIBUTE(DistanceConstraintSettings, mMaxDistance)
	JPH_ADD_ENUM_ATTRIBUTE_WITH_ALIAS(DistanceConstraintSettings, mLimitsSpringSettings.mMode, "mSpringMode")
	JPH_ADD_ATTRIBUTE_WITH_ALIAS(DistanceConstraintSettings, mLimitsSpringSettings.mFrequency, "mFrequency") // Renaming attributes to stay compatible with old versions of the library
	JPH_ADD_ATTRIBUTE_WITH_ALIAS(DistanceConstraintSettings, mLimitsSpringSettings.mDamping, "mDamping")
	JPH_ADD_ATTRIBUTE(DistanceConstraintSettings, mUseVelocityOnlySolver)
	JPH_ADD_ATTRIBUTE(DistanceConstraintSettings, mPredictiveLimitActivation)
	JPH_ADD_ATTRIBUTE(DistanceConstraintSettings, mVelocityBiasFactor)
	JPH_ADD_ATTRIBUTE(DistanceConstraintSettings, mVelocityRelaxation)
}

void DistanceConstraintSettings::SaveBinaryState(StreamOut &inStream) const
{
	ConstraintSettings::SaveBinaryState(inStream);

	inStream.Write(mSpace);
	inStream.Write(mPoint1);
	inStream.Write(mPoint2);
	inStream.Write(mMinDistance);
	inStream.Write(mMaxDistance);
	mLimitsSpringSettings.SaveBinaryState(inStream);
	inStream.Write(mUseVelocityOnlySolver);
	inStream.Write(mPredictiveLimitActivation);
	inStream.Write(mVelocityBiasFactor);
	inStream.Write(mVelocityRelaxation);
}

void DistanceConstraintSettings::RestoreBinaryState(StreamIn &inStream)
{
	ConstraintSettings::RestoreBinaryState(inStream);

	inStream.Read(mSpace);
	inStream.Read(mPoint1);
	inStream.Read(mPoint2);
	inStream.Read(mMinDistance);
	inStream.Read(mMaxDistance);
	mLimitsSpringSettings.RestoreBinaryState(inStream);
	inStream.Read(mUseVelocityOnlySolver);
	inStream.Read(mPredictiveLimitActivation);
	inStream.Read(mVelocityBiasFactor);
	inStream.Read(mVelocityRelaxation);
}

TwoBodyConstraint *DistanceConstraintSettings::Create(Body &inBody1, Body &inBody2) const
{
	return new DistanceConstraint(inBody1, inBody2, *this);
}

DistanceConstraint::DistanceConstraint(Body &inBody1, Body &inBody2, const DistanceConstraintSettings &inSettings) :
	TwoBodyConstraint(inBody1, inBody2, inSettings),
	mMinDistance(inSettings.mMinDistance),
	mMaxDistance(inSettings.mMaxDistance),
	mUseVelocityOnlySolver(inSettings.mUseVelocityOnlySolver),
	mPredictiveLimitActivation(inSettings.mPredictiveLimitActivation),
	mVelocityBiasFactor(inSettings.mVelocityBiasFactor),
	mVelocityRelaxation(inSettings.mVelocityRelaxation)
{
	if (inSettings.mSpace == EConstraintSpace::WorldSpace)
	{
		// If all properties were specified in world space, take them to local space now
		mLocalSpacePosition1 = Vec3(inBody1.GetInverseCenterOfMassTransform() * inSettings.mPoint1);
		mLocalSpacePosition2 = Vec3(inBody2.GetInverseCenterOfMassTransform() * inSettings.mPoint2);
		mWorldSpacePosition1 = inSettings.mPoint1;
		mWorldSpacePosition2 = inSettings.mPoint2;
	}
	else
	{
		// If properties were specified in local space, we need to calculate world space positions
		mLocalSpacePosition1 = Vec3(inSettings.mPoint1);
		mLocalSpacePosition2 = Vec3(inSettings.mPoint2);
		mWorldSpacePosition1 = inBody1.GetCenterOfMassTransform() * inSettings.mPoint1;
		mWorldSpacePosition2 = inBody2.GetCenterOfMassTransform() * inSettings.mPoint2;
	}

	// Store distance we want to keep between the world space points
	float distance = Vec3(mWorldSpacePosition2 - mWorldSpacePosition1).Length();
	float min_distance, max_distance;
	if (mMinDistance < 0.0f && mMaxDistance < 0.0f)
	{
		min_distance = max_distance = distance;
	}
	else
	{
		min_distance = mMinDistance < 0.0f? min(distance, mMaxDistance) : mMinDistance;
		max_distance = mMaxDistance < 0.0f? max(distance, mMinDistance) : mMaxDistance;
	}
	SetDistance(min_distance, max_distance);

	// Most likely gravity is going to tear us apart (this is only used when the distance between the points = 0)
	mWorldSpaceNormal = Vec3::sAxisY();

	// Store spring settings
	SetLimitsSpringSettings(inSettings.mLimitsSpringSettings);
}

void DistanceConstraint::NotifyShapeChanged(const BodyID &inBodyID, Vec3Arg inDeltaCOM)
{
	if (mBody1->GetID() == inBodyID)
		mLocalSpacePosition1 -= inDeltaCOM;
	else if (mBody2->GetID() == inBodyID)
		mLocalSpacePosition2 -= inDeltaCOM;
}

void DistanceConstraint::CalculateConstraintProperties(float inDeltaTime)
{
	// Update world space positions (the bodies may have moved)
	mWorldSpacePosition1 = mBody1->GetCenterOfMassTransform() * mLocalSpacePosition1;
	mWorldSpacePosition2 = mBody2->GetCenterOfMassTransform() * mLocalSpacePosition2;

	// Calculate world space normal
	Vec3 delta = Vec3(mWorldSpacePosition2 - mWorldSpacePosition1);
	float delta_len = delta.Length();
	if (delta_len > 0.0f)
		mWorldSpaceNormal = delta / delta_len;

	// Calculate points relative to body
	// r1 + u = (p1 - x1) + (p2 - p1) = p2 - x1
	Vec3 r1_plus_u = Vec3(mWorldSpacePosition2 - mBody1->GetCenterOfMassPosition());
	Vec3 r2 = Vec3(mWorldSpacePosition2 - mBody2->GetCenterOfMassPosition());

	float position_error;
	bool is_active = true;
	if (mMinDistance == mMaxDistance)
	{
		position_error = delta_len - mMinDistance;

		// Single distance, allow constraint forces in both directions
		mMinLambda = -FLT_MAX;
		mMaxLambda = FLT_MAX;
	}
	else if (delta_len <= mMinDistance)
	{
		position_error = delta_len - mMinDistance;

		// Allow constraint forces to make distance bigger only
		mMinLambda = 0;
		mMaxLambda = FLT_MAX;
	}
	else if (delta_len >= mMaxDistance)
	{
		position_error = delta_len - mMaxDistance;

		// Allow constraint forces to make distance smaller only
		mMinLambda = -FLT_MAX;
		mMaxLambda = 0;
	}
	else
	{
		position_error = 0.0f;
		is_active = false;

		if (mUseVelocityOnlySolver && mPredictiveLimitActivation)
		{
			const RVec3 next_position1 = mWorldSpacePosition1 + mBody1->GetPointVelocity(mWorldSpacePosition1) * inDeltaTime;
			const RVec3 next_position2 = mWorldSpacePosition2 + mBody2->GetPointVelocity(mWorldSpacePosition2) * inDeltaTime;
			const float next_distance = Vec3(next_position2 - next_position1).Length();
			is_active = next_distance < mMinDistance || next_distance > mMaxDistance;
		}
	}

	if (!is_active)
	{
		mAxisConstraint.Deactivate();
		return;
	}

	if (mUseVelocityOnlySolver)
	{
		mAxisConstraint.CalculateConstraintProperties(*mBody1, r1_plus_u, *mBody2, r2, mWorldSpaceNormal);
		mTargetVelocityBias = mVelocityBiasFactor * position_error / inDeltaTime;
	}
	else
		mAxisConstraint.CalculateConstraintPropertiesWithSettingsForLimit(inDeltaTime, *mBody1, r1_plus_u, *mBody2, r2, mWorldSpaceNormal, 0.0f, position_error, mLimitsSpringSettings);
}

void DistanceConstraint::SetupVelocityConstraint(float inDeltaTime)
{
	CalculateConstraintProperties(inDeltaTime);
	if (mVelocityOnlyDebugEnabled)
		mVelocityOnlyDebugStats = {};
}

void DistanceConstraint::ResetWarmStart()
{
	mAxisConstraint.Deactivate();
}

void DistanceConstraint::WarmStartVelocityConstraint(float inWarmStartImpulseRatio)
{
	if (mUseVelocityOnlySolver)
		mAxisConstraint.ResetTotalLambda();
	else
		mAxisConstraint.WarmStart(*mBody1, *mBody2, mWorldSpaceNormal, inWarmStartImpulseRatio);
}

bool DistanceConstraint::SolveVelocityConstraint(float inDeltaTime)
{
	if (mAxisConstraint.IsActive())
	{
		if (mUseVelocityOnlySolver)
		{
			if (!mVelocityOnlyDebugEnabled)
				return mAxisConstraint.SolveVelocityConstraintWithRelaxation(*mBody1, *mBody2, mWorldSpaceNormal, mTargetVelocityBias, mVelocityRelaxation);

			const Vec3 linear_velocity_before = mBody1->GetLinearVelocity();
			const Vec3 angular_velocity_before = mBody1->GetAngularVelocity();
			const float lambda_before = mAxisConstraint.GetTotalLambda();
			const bool applied_impulse = mAxisConstraint.SolveVelocityConstraintWithRelaxation(*mBody1, *mBody2, mWorldSpaceNormal, mTargetVelocityBias, mVelocityRelaxation);
			const Vec3 linear_delta = mBody1->GetLinearVelocity() - linear_velocity_before;
			const Vec3 angular_delta = mBody1->GetAngularVelocity() - angular_velocity_before;
			const float delta_lambda = mAxisConstraint.GetTotalLambda() - lambda_before;

			VelocityOnlyDebugStats &stats = mVelocityOnlyDebugStats;
			const uint step = ++stats.mNumSteps;
			stats.mTotalLinearVelocityDeltaBody1 += linear_delta;
			stats.mTotalAngularVelocityDeltaBody1 += angular_delta;
			if (linear_delta.LengthSq() > stats.mPeakLinearVelocityDeltaBody1.LengthSq())
			{
				stats.mPeakLinearStep = step;
				stats.mPeakLinearVelocityDeltaBody1 = linear_delta;
				stats.mPeakLinearDeltaLambda = delta_lambda;
			}
			if (angular_delta.LengthSq() > stats.mPeakAngularVelocityDeltaBody1.LengthSq())
			{
				stats.mPeakAngularStep = step;
				stats.mPeakAngularVelocityDeltaBody1 = angular_delta;
			}
			return applied_impulse;
		}
		else
			return mAxisConstraint.SolveVelocityConstraint(*mBody1, *mBody2, mWorldSpaceNormal, mMinLambda, mMaxLambda);
	}
	else
		return false;
}

bool DistanceConstraint::SolvePositionConstraint(float inDeltaTime, float inBaumgarte)
{
	if (mUseVelocityOnlySolver)
		return false;

	if (mLimitsSpringSettings.mFrequency <= 0.0f) // When the spring is active, we don't need to solve the position constraint
	{
		float distance = Vec3(mWorldSpacePosition2 - mWorldSpacePosition1).Dot(mWorldSpaceNormal);

		// Calculate position error
		float position_error = 0.0f;
		if (distance < mMinDistance)
			position_error = distance - mMinDistance;
		else if (distance > mMaxDistance)
			position_error = distance - mMaxDistance;

		if (position_error != 0.0f)
		{
			// Update constraint properties (bodies may have moved)
			CalculateConstraintProperties(inDeltaTime);

			return mAxisConstraint.SolvePositionConstraint(*mBody1, *mBody2, mWorldSpaceNormal, position_error, inBaumgarte);
		}
	}

	return false;
}

#ifdef JPH_DEBUG_RENDERER
void DistanceConstraint::DrawConstraint(DebugRenderer *inRenderer) const
{
	// Draw constraint
	Vec3 delta = Vec3(mWorldSpacePosition2 - mWorldSpacePosition1);
	float len = delta.Length();
	if (len < mMinDistance)
	{
		RVec3 real_end_pos = mWorldSpacePosition1 + (len > 0.0f? delta * mMinDistance / len : Vec3(0, len, 0));
		inRenderer->DrawLine(mWorldSpacePosition1, mWorldSpacePosition2, Color::sGreen);
		inRenderer->DrawLine(mWorldSpacePosition2, real_end_pos, Color::sYellow);
	}
	else if (len > mMaxDistance)
	{
		RVec3 real_end_pos = mWorldSpacePosition1 + (len > 0.0f? delta * mMaxDistance / len : Vec3(0, len, 0));
		inRenderer->DrawLine(mWorldSpacePosition1, real_end_pos, Color::sGreen);
		inRenderer->DrawLine(real_end_pos, mWorldSpacePosition2, Color::sRed);
	}
	else
		inRenderer->DrawLine(mWorldSpacePosition1, mWorldSpacePosition2, Color::sGreen);

	// Draw constraint end points
	inRenderer->DrawMarker(mWorldSpacePosition1, Color::sWhite, 0.1f);
	inRenderer->DrawMarker(mWorldSpacePosition2, Color::sWhite, 0.1f);

	// Draw current length
	inRenderer->DrawText3D(0.5_r * (mWorldSpacePosition1 + mWorldSpacePosition2), StringFormat("%.2f", (double)len));
}
#endif // JPH_DEBUG_RENDERER

void DistanceConstraint::SaveState(StateRecorder &inStream) const
{
	TwoBodyConstraint::SaveState(inStream);

	mAxisConstraint.SaveState(inStream);
	inStream.Write(mWorldSpaceNormal); // When distance = 0, the normal is used from last frame so we need to store it
}

void DistanceConstraint::RestoreState(StateRecorder &inStream)
{
	TwoBodyConstraint::RestoreState(inStream);

	mAxisConstraint.RestoreState(inStream);
	inStream.Read(mWorldSpaceNormal);
}

Ref<ConstraintSettings> DistanceConstraint::GetConstraintSettings() const
{
	DistanceConstraintSettings *settings = new DistanceConstraintSettings;
	ToConstraintSettings(*settings);
	settings->mSpace = EConstraintSpace::LocalToBodyCOM;
	settings->mPoint1 = RVec3(mLocalSpacePosition1);
	settings->mPoint2 = RVec3(mLocalSpacePosition2);
	settings->mMinDistance = mMinDistance;
	settings->mMaxDistance = mMaxDistance;
	settings->mLimitsSpringSettings = mLimitsSpringSettings;
	settings->mUseVelocityOnlySolver = mUseVelocityOnlySolver;
	settings->mPredictiveLimitActivation = mPredictiveLimitActivation;
	settings->mVelocityBiasFactor = mVelocityBiasFactor;
	settings->mVelocityRelaxation = mVelocityRelaxation;
	return settings;
}

JPH_NAMESPACE_END
