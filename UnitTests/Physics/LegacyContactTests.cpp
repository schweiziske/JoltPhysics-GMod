// SPDX-License-Identifier: MIT

#include "UnitTestFramework.h"
#include "PhysicsTestContext.h"
#include <Jolt/Core/FPException.h>
#include <Jolt/Core/TempAllocator.h>
#include <Jolt/Physics/Collision/ContactListener.h>
#include <Jolt/Physics/StateRecorder.h>
#include <cmath>
#include <cstring>

namespace
{
	class FilledTempAllocator final : public TempAllocator
	{
	public:
		explicit FilledTempAllocator(int inFill) : mFill(inFill) { }
		void *Allocate(uint inSize) override
		{
			void *result = mStorage.Allocate(inSize);
			if (result != nullptr)
				std::memset(result, mFill, inSize);
			return result;
		}
		void Free(void *inAddress, uint inSize) override { mStorage.Free(inAddress, inSize); }
	private:
		TempAllocatorMalloc mStorage;
		int mFill;
	};

	// Contacts serialize their Float3 state separately. Check all such values so that
	// a bad spring cache is detected even if the outgoing body velocity is still finite.
	class ContactStateRecorder final : public StateRecorder
	{
	public:
		void WriteBytes(const void *inData, size_t inSize) override
		{
			if (inSize == sizeof(Float3))
			{
				Float3 value;
				std::memcpy(&value, inData, sizeof(value));
				mTriples.push_back(value);
				CHECK(std::isfinite(value.x));
				CHECK(std::isfinite(value.y));
				CHECK(std::isfinite(value.z));
			}
		}
		void ReadBytes(void *, size_t) override { }
		bool IsEOF() const override { return false; }
		bool IsFailed() const override { return false; }
		Array<Float3> mTriples;
	};

	class LegacyImpactListener final : public ContactListener
	{
	public:
		void Configure(const ContactManifold &inManifold, ContactSettings &ioSettings)
		{
			ioSettings.mVPhysicsUseLegacySolver = mLegacy;
			ioSettings.mVPhysicsImpactingContact = mLegacy && mImpacting;
			ioSettings.mCombinedFriction = mLegacy? mFriction : 1.0f;
			ioSettings.mCombinedRestitution = 0.0f;
			ioSettings.mRelativeLinearSurfaceVelocity = Vec3(mLegacy? 0.0f : mPrimeVelocity, 0, 0);
			mMaxPoints = max(mMaxPoints, uint(inManifold.mRelativeContactPointsOn1.size()));
		}
		void OnContactAdded(const Body &, const Body &, const ContactManifold &inManifold, ContactSettings &ioSettings) override { Configure(inManifold, ioSettings); }
		void OnContactPersisted(const Body &, const Body &, const ContactManifold &inManifold, ContactSettings &ioSettings) override { Configure(inManifold, ioSettings); }
		void OnVPhysicsSolverStats(const VPhysicsSolverStats &inStats) override { mImpacts += inStats.mImpactCount; }
		bool mLegacy = true;
		bool mImpacting = true;
		float mFriction = 0.0f;
		float mPrimeVelocity = 0.0f;
		uint mMaxPoints = 0;
		uint mImpacts = 0;
	};
}

TEST_SUITE("LegacyContactTests")
{
	TEST_CASE("TestZeroFrictionMultiPointImpactDirtyMemory")
	{
		// Preserve a test failure instead of terminating on the NaNs produced by the old code.
		FPExceptionDisableInvalid disable_invalid;
		FPExceptionDisableDivByZero disable_divzero;
		FPExceptionDisableOverflow disable_overflow;
		Array<Float3> reference;
		for (int fill : { 0, 255 })
		{
			PhysicsTestContext context;
			context.ZeroGravity();
			LegacyImpactListener listener;
			PhysicsSystem *world = context.GetSystem();
			world->SetContactListener(&listener);
			PhysicsSettings settings = world->GetPhysicsSettings();
			settings.mContactSeparationTarget = 0.02f;
			world->SetPhysicsSettings(settings);
			context.CreateFloor();
			Body &body = context.CreateBox(RVec3(0, 0.499f, 0), Quat::sIdentity(), EMotionType::Dynamic,
				EMotionQuality::Discrete, Layers::MOVING, Vec3::sReplicate(0.5f));
			body.SetLinearVelocity(Vec3(1, -1, 0));
			FilledTempAllocator allocator(fill);
			CHECK(world->Update(context.GetDeltaTime(), 1, &allocator, context.GetJobSystem()) == EPhysicsUpdateError::None);
			CHECK(listener.mImpacts > 0);
			CHECK(listener.mMaxPoints == 4);
			ContactStateRecorder recorder;
			world->SaveState(recorder, EStateRecorderState::Contacts);
			if (fill == 0)
				reference = recorder.mTriples;
			else
			{
				CHECK(recorder.mTriples.size() == reference.size());
				for (size_t i = 0; i < min(reference.size(), recorder.mTriples.size()); ++i)
					CHECK_APPROX_EQUAL(Vec3(recorder.mTriples[i]), Vec3(reference[i]), 1.0e-6f);
			}
			// A later frictional contact must not inherit a poisoned impact spring.
			listener.mFriction = 0.5f;
			listener.mImpacting = false;
			for (int tick = 0; tick < 10; ++tick)
			{
				CHECK(context.SimulateSingleStep() == EPhysicsUpdateError::None);
				CHECK_FALSE(body.GetLinearVelocity().IsNaN());
				CHECK_FALSE(body.GetAngularVelocity().IsNaN());
			}
			world->SaveState(recorder, EStateRecorderState::Contacts);
			world->SetContactListener(nullptr);
		}
	}

	TEST_CASE("TestZeroFrictionImpactAfterNativeFriction")
	{
		Array<Float3> reference;
		for (float prime_velocity : { 0.0f, 1000.0f })
		{
			PhysicsTestContext context;
			context.ZeroGravity();
			LegacyImpactListener listener;
			listener.mLegacy = false;
			listener.mPrimeVelocity = prime_velocity;
			PhysicsSystem *world = context.GetSystem();
			world->SetContactListener(&listener);
			context.CreateFloor();
			Body &body = context.CreateBox(RVec3(0, 0.499f, 0), Quat::sIdentity(), EMotionType::Dynamic,
				EMotionQuality::Discrete, Layers::MOVING, Vec3::sReplicate(0.5f));
			CHECK(context.SimulateSingleStep() == EPhysicsUpdateError::None);
			listener.mLegacy = true;
			body.SetLinearVelocity(Vec3(1, -1, 0));
			CHECK(context.SimulateSingleStep() == EPhysicsUpdateError::None);
			CHECK(listener.mImpacts > 0);
			CHECK(listener.mMaxPoints == 4);
			ContactStateRecorder recorder;
			world->SaveState(recorder, EStateRecorderState::Contacts);
			if (prime_velocity == 0.0f)
				reference = recorder.mTriples;
			else
			{
				CHECK(recorder.mTriples.size() == reference.size());
				for (size_t i = 0; i < min(reference.size(), recorder.mTriples.size()); ++i)
					CHECK_APPROX_EQUAL(Vec3(recorder.mTriples[i]), Vec3(reference[i]), 1.0e-6f);
			}
			world->SetContactListener(nullptr);
		}
	}
}
