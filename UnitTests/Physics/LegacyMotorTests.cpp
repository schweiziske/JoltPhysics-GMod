// SPDX-License-Identifier: MIT

#include "UnitTestFramework.h"
#include "PhysicsTestContext.h"
#include <Jolt/Physics/Constraints/HingeConstraint.h>
#include <Jolt/Physics/StateRecorderImpl.h>
#include <cmath>

TEST_SUITE("LegacyMotorTests")
{
	TEST_CASE("TestAngularReferenceMultiTurnAndRestore")
	{
		constexpr float dt = 1.0f / 60.0f;
		for (float direction : { -1.0f, 1.0f })
		{
			LegacyMotorState state;
			state.Configure(direction * 2.0f, 1000.0f);
			// Also accept references accumulated by older saved states.
			state.SetReferencePosition(direction * 100.0f * JPH_PI + 0.25f);
			StateRecorderImpl recorder;
			state.SaveState(recorder);
			recorder.Rewind();
			LegacyMotorState restored;
			restored.RestoreState(recorder);
			CHECK_FALSE(recorder.IsFailed());
			for (int tick = 0; tick < 2400; ++tick)
			{
				if (tick == 1200)
				{
					state.Configure(-direction * 2.0f, 1000.0f);
					restored.Configure(-direction * 2.0f, 1000.0f);
				}
				const float angle = std::remainder(state.GetReferencePosition(), 2.0f * JPH_PI);
				const LegacyMotorStep step = state.BeginAngularStep(angle, dt);
				const LegacyMotorStep restored_step = restored.BeginAngularStep(angle, dt);
				CHECK_APPROX_EQUAL(step.mBias, 0.0f, 1.0e-3f);
				CHECK_APPROX_EQUAL(step.mBias, restored_step.mBias);
				state.FinishStep(0.0f, 0.0f);
				restored.FinishStep(0.0f, 0.0f);
				CHECK(std::abs(state.GetReferencePosition()) <= JPH_PI);
			}
		}
	}

	TEST_CASE("TestAngularSaturationAcrossWrap")
	{
		for (float direction : { -1.0f, 1.0f })
		{
			LegacyMotorState state;
			state.Configure(direction * 2.0f, 0.1f);
			state.SetReferencePosition(direction * (20.0f * JPH_PI + 3.1f));
			const float angle = direction * -3.0f;
			const float error = std::remainder(state.GetReferencePosition() - angle, 2.0f * JPH_PI);
			const LegacyMotorStep step = state.BeginAngularStep(angle, 0.1f);
			CHECK_APPROX_EQUAL(step.mBias, error * 8.0f, 1.0e-4f);
			// A quarter of the requested impulse is available: roll back 3/4 of the error.
			state.FinishStep(0.04f, 0.01f);
			const float expected = angle + 0.25f * error + direction * 0.2f;
			CHECK_APPROX_EQUAL(std::remainder(state.GetReferencePosition() - expected, 2.0f * JPH_PI), 0.0f, 1.0e-5f);
			CHECK(std::abs(state.GetReferencePosition()) <= JPH_PI);
		}
	}

	TEST_CASE("TestLegacyHingeMotorMultiTurnAndReversal")
	{
		for (float direction : { -1.0f, 1.0f })
			for (float torque : { 100.0f, 10000.0f })
			{
				PhysicsTestContext context;
				context.ZeroGravity();
				Body &body = context.CreateBox(RVec3::sZero(), Quat::sIdentity(), EMotionType::Dynamic,
					EMotionQuality::Discrete, Layers::MOVING, Vec3(0.1f, 0.2f, 0.3f));
				body.GetMotionProperties()->SetAngularDamping(0.0f);
				body.SetAllowSleeping(false);
				HingeConstraintSettings settings;
				settings.mSpace = EConstraintSpace::LocalToBodyCOM;
				settings.mPoint1 = settings.mPoint2 = Vec3::sZero();
				settings.mHingeAxis1 = settings.mHingeAxis2 = Vec3::sAxisX();
				settings.mNormalAxis1 = settings.mNormalAxis2 = Vec3::sAxisY();
				HingeConstraint &hinge = context.CreateConstraint<HingeConstraint>(Body::sFixedToWorld, body, settings);
				hinge.SetLegacyVelocitySolverSettings(1.0f, 1.0f);
				hinge.SetLegacyMotor(direction * 2.0f, torque);
				LegacyConstraintGroup group;
				group.AddConstraint(&hinge);
				for (int tick = 0; tick < 2400; ++tick)
				{
					if (tick == 1200)
						hinge.SetLegacyMotor(-direction * 2.0f, torque);
					CHECK(context.SimulateSingleStep() == EPhysicsUpdateError::None);
					const float velocity = body.GetAngularVelocity().GetX();
					// Check every step, including direction reversal and torque-limited acceleration.
					CHECK(std::abs(velocity) < 2.1f);
					if ((tick > 30 && tick < 1200) || tick > 1230)
						CHECK_APPROX_EQUAL(velocity, (tick < 1200? -direction : direction) * 2.0f, 0.01f);
				}
				group.RemoveConstraint(&hinge);
			}
	}
}
