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
#include <cstdio>
#include <cstdlib>
#include <thread>

JPH_NAMESPACE_BEGIN

static bool IsLegacyGroupSolveTraceEnabled()
{
	static const bool enabled = []()
	{
#if defined(JPH_COMPILER_MSVC)
		char *value = nullptr;
		size_t length = 0;
		const bool found = _dupenv_s(&value, &length, "VJOLT_TRACE_LEGACY_GROUP_SOLVES") == 0
			&& value != nullptr;
		std::free(value);
		return found;
#else
		return std::getenv("VJOLT_TRACE_LEGACY_GROUP_SOLVES") != nullptr;
#endif
	}();
	return enabled;
}

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
	bool enabled = false;
	inStream.Read(enabled);
	SetEnabled(enabled);
	mLegacyBroken.store(false, std::memory_order_release);
}

void Constraint::SetEnabled(bool inEnabled)
{
	if (mEnabled == inEnabled)
		return;

	mEnabled = inEnabled;
	if (mLegacyConstraintGroup != nullptr)
		mLegacyConstraintGroup->InvalidateSolverMetadata();
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

void LegacyConstraintGroup::RebuildSolverStorageLayout()
{
	// Havana keeps the total VMQ storage size on the local constraint system and only rebuilds it
	// when members change. Do the same here instead of rediscovering sizes and alignments every PSI.
	mSolverStorageOffsets.resize(mConstraints.size());
	size_t storage_size = 0;
	for (size_t index = 0; index < mConstraints.size(); ++index)
	{
		Constraint *constraint = mConstraints[index];
		const size_t constraint_size = constraint->GetLegacySolverStorageSize();
		if (constraint_size == 0)
		{
			mSolverStorageOffsets[index] = SIZE_MAX;
			continue;
		}

		const size_t alignment = constraint->GetLegacySolverStorageAlignment();
		JPH_ASSERT(alignment > 0 && alignment <= alignof(std::max_align_t)
			&& (alignment & (alignment - 1)) == 0);
		storage_size = (storage_size + alignment - 1) & ~(alignment - 1);
		mSolverStorageOffsets[index] = storage_size;
		storage_size += constraint_size;
	}

	const size_t storage_words = (storage_size + sizeof(std::max_align_t) - 1)
		/ sizeof(std::max_align_t);
	mSolverStorage.resize(storage_words);
	mActiveSolverEntries.reserve(mConstraints.size());
}

void LegacyConstraintGroup::AddConstraint(Constraint *inConstraint)
{
	if (std::find(mConstraints.begin(), mConstraints.end(), inConstraint) == mConstraints.end())
	{
		mConstraints.push_back(inConstraint);
		RebuildSolverStorageLayout();
		InvalidateSolverMetadata();
	}
	inConstraint->SetLegacyConstraintGroup(this);
}

void LegacyConstraintGroup::RemoveConstraint(Constraint *inConstraint)
{
	if (inConstraint->GetLegacyConstraintGroup() == this)
		inConstraint->SetLegacyConstraintGroup(nullptr);
	if (std::erase(mConstraints, inConstraint) != 0)
	{
		RebuildSolverStorageLayout();
		InvalidateSolverMetadata();
	}
}

void LegacyConstraintGroup::InvalidateSolverMetadata()
{
	mSolverMetadata.store(0, std::memory_order_release);
}

void LegacyConstraintGroup::RefreshSolverMetadata() const
{
	if ((mSolverMetadata.load(std::memory_order_acquire) & cSolverMetadataValid) != 0)
		return;

	Constraint *representative = nullptr;
	bool requires_post_contact_solve = false;
	for (Constraint *constraint : mConstraints)
		if (constraint->IsActive())
		{
			if (representative == nullptr)
				representative = constraint;
			if (constraint->RequiresPostContactLegacySolve())
				requires_post_contact_solve = true;
			else
			{
				requires_post_contact_solve = false;
				break;
			}
		}

	mSolverRepresentative.store(representative, std::memory_order_relaxed);
	const uint8 metadata = cSolverMetadataValid
		| (representative != nullptr && requires_post_contact_solve? cRequiresPostContactSolve : 0);
	mSolverMetadata.store(metadata, std::memory_order_release);
}

void LegacyConstraintGroup::PrepareVelocitySolve()
{
	uint8 prepare_state = mSolverPrepareState.load(std::memory_order_acquire);
	if (prepare_state == cSolverPrepareReady)
		return;

	uint8 expected = cSolverPrepareUnprepared;
	if (mSolverPrepareState.compare_exchange_strong(expected, cSolverPrepareInitializing,
		std::memory_order_acq_rel, std::memory_order_acquire))
	{
		// Body activation and runtime constraint modes can change without changing group
		// membership. Rebuild this small metadata snapshot once per collision step, after
		// Jolt has frozen the active-constraint list for that step.
		InvalidateSolverMetadata();
		RefreshSolverMetadata();
		mSolvePending.store(true, std::memory_order_release);
		mSolverPrepareState.store(cSolverPrepareReady, std::memory_order_release);
		return;
	}

	// Setup runs in parallel batches. Only the first member scans the group; other members
	// wait for that short publication instead of repeating the scan or writing mSolvePending.
	while (mSolverPrepareState.load(std::memory_order_acquire) == cSolverPrepareInitializing)
		std::this_thread::yield();
}

bool LegacyConstraintGroup::IsSolverRepresentative(const Constraint *inCaller) const
{
	RefreshSolverMetadata();
	return mSolverRepresentative.load(std::memory_order_relaxed) == inCaller;
}

bool LegacyConstraintGroup::RequiresPostContactSolve() const
{
	RefreshSolverMetadata();
	return (mSolverMetadata.load(std::memory_order_relaxed) & cRequiresPostContactSolve) != 0;
}

bool LegacyConstraintGroup::SolveVelocityConstraints(Constraint *inCaller, float inDeltaTime)
{
	if (!IsSolverRepresentative(inCaller))
		return false;
	if (!mSolvePending.exchange(false, std::memory_order_acq_rel))
		return false;

	bool appliedImpulse = false;
	bool groupBroken = false;
	// RebuildSolverStorageLayout has already packed every Havana work area contiguously. The
	// backing vector only changes when group membership changes, so the steady-state PSI path does
	// not allocate, resize or query virtual storage sizes.
	std::byte *storage_base = reinterpret_cast<std::byte *>(mSolverStorage.data());
	mActiveSolverEntries.clear();
	for (size_t index = 0; index < mConstraints.size(); ++index)
	{
		Constraint *constraint = mConstraints[index];
		if (!constraint->IsActive())
			continue;
		const size_t storage_offset = mSolverStorageOffsets[index];
		void *solver_storage = storage_offset == SIZE_MAX? nullptr : storage_base + storage_offset;
		mActiveSolverEntries.push_back({ constraint, solver_storage, index });
		if (solver_storage != nullptr)
			constraint->InitializeLegacySolverStorage(solver_storage);
	}
	auto destroy_solver_storage = [&]()
	{
		for (const SolverEntry &entry : mActiveSolverEntries)
			if (entry.mStorage != nullptr)
				entry.mConstraint->DestroyLegacySolverStorage(entry.mStorage);
	};
	const bool traceGroupSolves = IsLegacyGroupSolveTraceEnabled();
	uint32 solveOrdinal = 0;
	if (traceGroupSolves)
		std::fprintf(stderr, "LEGACY_GROUP_BEGIN group=%p count=%zu additional=%d caller=%p\n",
			static_cast<void *>(this), mConstraints.size(), mAdditionalIterations,
			static_cast<void *>(inCaller));
	auto solveConstraint = [&](const SolverEntry &entry, float tauFactor, float dampingFactor,
		const char *phase, int64 iteration, bool setup)
	{
		if (groupBroken)
			return;
		++solveOrdinal;
		Constraint *constraint = entry.mConstraint;
		const bool check_break = constraint->HasLegacyBreakLimits();
		const bool capture_body_state = traceGroupSolves || check_break;
		const TwoBodyConstraint *two_body = capture_body_state
			&& constraint->GetType() == EConstraintType::TwoBodyConstraint?
				static_cast<const TwoBodyConstraint *>(constraint) : nullptr;
		const Vec3 linear1 = two_body? two_body->GetBody1()->GetLinearVelocity() : Vec3::sZero();
		const Vec3 angular1 = two_body? two_body->GetBody1()->GetAngularVelocity() : Vec3::sZero();
		const Vec3 linear2 = two_body? two_body->GetBody2()->GetLinearVelocity() : Vec3::sZero();
		const Vec3 angular2 = two_body? two_body->GetBody2()->GetAngularVelocity() : Vec3::sZero();
		appliedImpulse |= setup
			? constraint->SetupAndSolveLegacyVelocityConstraint(
				entry.mStorage, inDeltaTime, tauFactor, dampingFactor)
			: constraint->StepLegacyVelocityConstraint(
				entry.mStorage, inDeltaTime, tauFactor, dampingFactor);
		if (traceGroupSolves && two_body)
		{
			const Vec3 postLinear1 = two_body->GetBody1()->GetLinearVelocity();
			const Vec3 postAngular1 = two_body->GetBody1()->GetAngularVelocity();
			const Vec3 postLinear2 = two_body->GetBody2()->GetLinearVelocity();
			const Vec3 postAngular2 = two_body->GetBody2()->GetAngularVelocity();
			std::fprintf(stderr,
				"LEGACY_GROUP_SOLVE group=%p ordinal=%u phase=%s iteration=%lld index=%zu subtype=%u b1=%u pre_v=(%.9g,%.9g,%.9g) pre_w=(%.9g,%.9g,%.9g) post_v=(%.9g,%.9g,%.9g) post_w=(%.9g,%.9g,%.9g) b2=%u pre_v=(%.9g,%.9g,%.9g) pre_w=(%.9g,%.9g,%.9g) post_v=(%.9g,%.9g,%.9g) post_w=(%.9g,%.9g,%.9g)\n",
				static_cast<void *>(this), solveOrdinal, phase, static_cast<long long>(iteration),
				entry.mConstraintIndex, uint(constraint->GetSubType()),
				two_body->GetBody1()->GetID().GetIndexAndSequenceNumber(),
				linear1.GetX(), linear1.GetY(), linear1.GetZ(), angular1.GetX(), angular1.GetY(), angular1.GetZ(),
				postLinear1.GetX(), postLinear1.GetY(), postLinear1.GetZ(), postAngular1.GetX(), postAngular1.GetY(), postAngular1.GetZ(),
				two_body->GetBody2()->GetID().GetIndexAndSequenceNumber(),
				linear2.GetX(), linear2.GetY(), linear2.GetZ(), angular2.GetX(), angular2.GetY(), angular2.GetZ(),
				postLinear2.GetX(), postLinear2.GetY(), postLinear2.GetZ(), postAngular2.GetX(), postAngular2.GetY(), postAngular2.GetZ());
		}
		if (check_break && constraint->CheckLegacyBreakAfterSolve(
			linear1, angular1, linear2, angular2, mRollbackOnBreak)
			&& mRollbackOnBreak)
		{
			constraint->SetEnabled(false);
			groupBroken = true;
		}
	};
	auto solveForward = [&](float tauFactor, float dampingFactor, int64 iteration)
	{
		for (const SolverEntry &entry : mActiveSolverEntries)
			solveConstraint(entry, tauFactor, dampingFactor, "forward", iteration, false);
	};
	auto solveReverse = [&](float tauFactor, float dampingFactor, int64 iteration)
	{
		for (auto entry = mActiveSolverEntries.rbegin(); entry != mActiveSolverEntries.rend(); ++entry)
			solveConstraint(*entry, tauFactor, dampingFactor, "reverse", iteration, false);
	};
	auto setupAndSolveForward = [&](float tauFactor, float dampingFactor)
	{
		// Havana's setup_and_step_constraint is interleaved per constraint. Do not set up the
		// complete group first: a later predictive Length constraint must observe velocity impulses
		// written by an earlier Ragdoll/Fixed setup solve in the same local controller pass.
		for (const SolverEntry &entry : mActiveSolverEntries)
		{
			if (groupBroken)
				break;
			solveConstraint(entry, tauFactor, dampingFactor, "setup", -1, true);
		}
	};

	// GMod's shipped vphysics binary performs one setup solve followed by one reverse/forward scan
	// per local-system iteration. This differs from the duplicated x < 2 phase present in the
	// available Havana source snapshot. A staged Length fingerprint proves the binary behavior:
	// additionalIterations -2/-1/0 produces 1/3/5 solves respectively, while the duplicated loop
	// produces 1/5/9. Match the game binary that this adapter replaces.
	setupAndSolveForward(1.0f, 1.0f);
	// Havana stores m_n_iterations as signed additionalIterations + 2. Negative input is not
	// clamped: -1 performs one iteration, while <= -2 performs none after the setup solve.
	const int64 num_iterations = int64(mAdditionalIterations) + 2;
	static constexpr float tau_factors[] = { 1.0f, 1.0f, 0.8f, 0.6f, 0.4f, 0.4f, 0.4f, 0.4f, 0.4f, 0.0f };
	static constexpr float damping_factors[] = { 1.0f, 1.0f, 0.8f, 0.8f, 0.8f, 0.8f, 0.8f, 0.8f, 0.8f, 0.0f };
	for (int64 iteration = 0;
		iteration < num_iterations && iteration < int64(std::size(tau_factors))
			&& tau_factors[iteration] != 0.0f && !groupBroken;
		++iteration)
	{
		solveReverse(tau_factors[iteration], damping_factors[iteration], iteration);
		solveForward(tau_factors[iteration], damping_factors[iteration], iteration);
	}
	if (traceGroupSolves)
		std::fprintf(stderr, "LEGACY_GROUP_END group=%p solves=%u broken=%d\n",
			static_cast<void *>(this), solveOrdinal, groupBroken ? 1 : 0);
	destroy_solver_storage();
	mSolverPrepareState.store(cSolverPrepareUnprepared, std::memory_order_release);

	(void)inCaller;
	return appliedImpulse;
}

JPH_NAMESPACE_END
