// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-FileCopyrightText: 2021 Jorrit Rouwe
// SPDX-License-Identifier: MIT

#pragma once

#include <Jolt/Physics/Constraints/Constraint.h>
#include <Jolt/Physics/Body/BodyID.h>
#include <Jolt/Physics/PhysicsLock.h>
#include <Jolt/Core/Mutex.h>
#include <Jolt/Core/JobSystem.h>
#include <Jolt/Core/UnorderedMap.h>
#include <functional>
#include <memory>
#include <span>

JPH_NAMESPACE_BEGIN

class IslandBuilder;
class BodyManager;
class Body;
class StateRecorderFilter;
#ifdef JPH_DEBUG_RENDERER
class DebugRenderer;
#endif // JPH_DEBUG_RENDERER

/// A list of constraints
using Constraints = Array<Ref<Constraint>>;

/// A constraint manager manages all constraints of the same type
///
/// WARNING: This class is an internal part of PhysicsSystem, it has no functions that can be called by users of the library.
/// Its functionality is exposed through PhysicsSystem and BodyInterface.
class JPH_EXPORT ConstraintManager : public NonCopyable
{
public:
	JPH_OVERRIDE_NEW_DELETE

#ifdef JPH_ENABLE_ASSERTS
	/// Constructor
							ConstraintManager(PhysicsLockContext inContext) : mLockContext(inContext) { }
#endif // JPH_ENABLE_ASSERTS

	/// Add a new constraint. This is thread safe.
	void					Add(Constraint **inConstraints, int inNumber);

	/// Remove a constraint. This is thread safe.
	void					Remove(Constraint **inConstraint, int inNumber);

	/// Get a list of all constraints
	Constraints				GetConstraints() const;
	Constraints				GetConstraintsNoLock() const { return mConstraints; }

	/// Get total number of constraints
	inline uint32			GetNumConstraints() const					{ return uint32(mConstraints.size()); }
	/// Registered constraints with velocity-independent legacy preparation (including sleeping ones).
	inline uint32			GetNumLegacyPreparationConstraints() const	{ return mNumLegacyPreparationConstraints; }

	/// Close enabled two-body and legacy-group dynamic wake components before the active
	/// constraint snapshot. The caller holds all body/constraint locks and grants activation.
	/// This does not join solver islands or change legacy solver member slices.
	struct WakeStatistics
	{
		uint32 mScannedConstraints = 0;
		uint32 mEligibleConstraints = 0;
		uint32 mDynamicBodies = 0;
		uint32 mMixedComponents = 0;
		uint32 mActivationBatchBodies = 0;
		uint32 mNewlyActivatedBodies = 0;
	};
	~ConstraintManager();
	WakeStatistics			PropagateWake(BodyManager &ioBodyManager);
	bool					NeedsWakePropagation(const BodyManager &inBodyManager) const;
	void					InvalidateWakeTopology() { mWakeTopologyRevision.fetch_add(1, memory_order_release); }

	/// Determine the active constraints of a subset of the constraints
	void					GetActiveConstraints(uint32 inStartConstraintIdx, uint32 inEndConstraintIdx, Constraint **outActiveConstraints, uint32 &outNumActiveConstraints) const;

	/// Link bodies to form islands
	static void				sBuildIslands(Constraint **inActiveConstraints, uint32 inNumActiveConstraints, IslandBuilder &ioBuilder, BodyManager &inBodyManager);

	/// In order to have a deterministic simulation, we need to sort the constraints of an island before solving them
	static void				sSortConstraints(Constraint **inActiveConstraints, uint32 *inConstraintIdxBegin, uint32 *inConstraintIdxEnd);

	/// Prior to solving the velocity constraints, you must call SetupVelocityConstraints once to precalculate values that are independent of velocity
	static void				sSetupVelocityConstraints(Constraint **inActiveConstraints, uint32 inNumActiveConstraints, float inDeltaTime, bool inPrepareLegacy = false);
	/// Release preparation storage for this island before publishing contact batches, even
	/// when the velocity iteration count is zero or a group stopped early after breaking.
	static void				sReleasePreparedLegacyVelocityConstraints(Constraint **inActiveConstraints, const uint32 *inConstraintIdxBegin, const uint32 *inConstraintIdxEnd);

	/// Apply last frame's impulses, must be called prior to SolveVelocityConstraints
	template <class ConstraintCallback>
	static void				sWarmStartVelocityConstraints(Constraint **inActiveConstraints, const uint32 *inConstraintIdxBegin, const uint32 *inConstraintIdxEnd, float inWarmStartImpulseRatio, ConstraintCallback &ioCallback);

	/// This function is called multiple times to iteratively come to a solution that meets all velocity constraints
	static bool				sSolveVelocityConstraints(Constraint **inActiveConstraints, const uint32 *inConstraintIdxBegin, const uint32 *inConstraintIdxEnd, float inDeltaTime);

	static bool				sHasLegacyConstraints(Constraint **inActiveConstraints, const uint32 *inConstraintIdxBegin, const uint32 *inConstraintIdxEnd);
	static bool				sSupportsDeferredContactWarmStart(Constraint **inActiveConstraints, const uint32 *inConstraintIdxBegin, const uint32 *inConstraintIdxEnd);

	/// One island's ordered groups. The scheduler owns this storage until every native
	/// velocity job returns, including readers that already observed a completed batch.
	struct LegacyVelocityBatch
	{
		struct Group
		{
			LegacyConstraintGroup *mGroup;
			Constraint *mFirstMember = nullptr;
			std::vector<Constraint *> mMembers;
			bool mPostContact = true;
			void AddMember(Constraint *inConstraint);
			std::span<Constraint *const> GetMembers() const;
			bool Solve(float inDeltaTime);
		};
		std::vector<Group> mGroups;
		std::vector<uint32> mGroupIndices, mLayerEnds;
		std::function<void(bool)> mContinuation;
		LegacyVelocityBatch *mNext = nullptr; // Immutable after publication
		float mDeltaTime;
		atomic<uint64> mCursor { 0 }; // High 32 bits: layer; low 32 bits: next group
		atomic<uint32> mProcessed { 0 };
		atomic<bool> mAppliedImpulse { false };
		void BuildGroups(Constraint **inActiveConstraints, const uint32 *inBegin, const uint32 *inEnd);
		bool BuildLayers(); // False when no parallel width or unknown body dependencies
		bool SolveSerial();
		bool TryExecute();
		bool IsDone() const { return mCursor.load(memory_order_acquire) == UINT64_MAX; }
	};

	class LegacyVelocityScheduler
	{
	public:
		void Prepare(const IslandBuilder &inBuilder, Constraint **inActiveConstraints, bool inEnabled);
		void Clear(); // All original velocity jobs must have finished
		uint32 GetNumIslands() const { return mNumIslands; }
		bool IsEnabled() const { return mEnabled; }
		bool HasPendingIslands() const { return mRemainingIslands.load(memory_order_acquire) != 0; }
		void FinishIsland()
		{
			if (mEnabled)
			{
				JPH_IF_ENABLE_ASSERTS(uint32 previous =) mRemainingIslands.fetch_sub(1, memory_order_release);
				JPH_ASSERT(previous > 0);
			}
		}
		void Publish(uint32 inIsland, std::unique_ptr<LegacyVelocityBatch> inBatch);
		bool TryExecute();
	private:
		std::vector<std::unique_ptr<LegacyVelocityBatch>> mBatches;
		atomic<LegacyVelocityBatch *> mHead { nullptr };
		atomic<uint32> mRemainingIslands { 0 };
		uint32 mNumIslands = 0;
		bool mEnabled = false;
	};

	/// Direct callers retain a synchronous, ordered controller pass.
	static bool				sSolveLegacyVelocityConstraints(Constraint **inActiveConstraints, const uint32 *inConstraintIdxBegin, const uint32 *inConstraintIdxEnd, float inDeltaTime);
	/// Publish ordered layers to the existing native velocity jobs, or invoke the
	/// continuation synchronously for small/serial islands. Never creates a pool job.
	template <class Continuation>
	static void sDispatchLegacyVelocityConstraints(Constraint **inActiveConstraints, const uint32 *inBegin, const uint32 *inEnd, float inDeltaTime, LegacyVelocityScheduler &ioScheduler, uint32 inIsland, Continuation inContinuation)
	{
		if (!ioScheduler.IsEnabled() || inEnd - inBegin < 128)
		{
			inContinuation(sSolveLegacyVelocityConstraints(inActiveConstraints, inBegin, inEnd, inDeltaTime));
			return;
		}
		auto batch = std::make_unique<LegacyVelocityBatch>();
		batch->mDeltaTime = inDeltaTime;
		batch->BuildGroups(inActiveConstraints, inBegin, inEnd);
		if (!batch->BuildLayers())
		{
			inContinuation(batch->SolveSerial());
			return;
		}
		// Only asynchronous islands need an owning callback allocation.
		batch->mContinuation = std::move(inContinuation);
		ioScheduler.Publish(inIsland, std::move(batch));
	}

	/// Refresh predictive micro-angle members after contacts. This removes contact-generated angular
	/// velocity without moving an otherwise mixed Weld/Rope group to the post-contact controller stage.
	static bool				sSolvePostContactLegacyVelocityMembers(Constraint **inActiveConstraints, const uint32 *inConstraintIdxBegin, const uint32 *inConstraintIdxEnd, float inDeltaTime);

	/// This function is called multiple times to iteratively come to a solution that meets all position constraints
	static bool				sSolvePositionConstraints(Constraint **inActiveConstraints, const uint32 *inConstraintIdxBegin, const uint32 *inConstraintIdxEnd, float inDeltaTime, float inBaumgarte);

	/// Re-apply predictive legacy angular limits after Jolt's contact position correction, including
	/// individual micro-angle members of an otherwise mixed Source constraint system.
	static bool				sSolvePostContactLegacyPositionConstraints(Constraint **inActiveConstraints, const uint32 *inConstraintIdxBegin, const uint32 *inConstraintIdxEnd, float inDeltaTime, float inBaumgarte);

#ifdef JPH_DEBUG_RENDERER
	/// Draw all constraints
	void					DrawConstraints(DebugRenderer *inRenderer) const;

	/// Draw all constraint limits
	void					DrawConstraintLimits(DebugRenderer *inRenderer) const;

	/// Draw all constraint reference frames
	void					DrawConstraintReferenceFrame(DebugRenderer *inRenderer) const;
#endif // JPH_DEBUG_RENDERER

	/// Save state of constraints
	void					SaveState(StateRecorder &inStream, const StateRecorderFilter *inFilter) const;

	/// Restore the state of constraints. Returns false if failed.
	bool					RestoreState(StateRecorder &inStream);

	/// Lock all constraints. This should only be done during PhysicsSystem::Update().
	void					LockAllConstraints()						{ PhysicsLock::sLock(mConstraintsMutex JPH_IF_ENABLE_ASSERTS(, mLockContext, EPhysicsLockTypes::ConstraintsList)); }
	void					UnlockAllConstraints()						{ PhysicsLock::sUnlock(mConstraintsMutex JPH_IF_ENABLE_ASSERTS(, mLockContext, EPhysicsLockTypes::ConstraintsList)); }

private:
	uint32					mNumLegacyPreparationConstraints = 0;
	struct WakeNode
	{
		BodyID				mBodyID;
		uint32				mParent;
		uint32				mSize = 1;
		bool				mHasActiveBody;
		bool				mHasSleepingBody;
	};

	uint32					GetOrCreateWakeNode(const Body &inBody);
	uint32					FindWakeRoot(uint32 inNode);
	void					LinkWakeNodes(uint32 inNode1, uint32 inNode2);

	// Connectivity changes only with constraint membership/eligibility or body motion type.
	// All those mutation sites invalidate; sleep/wake changes only rescan cached bodies.
	atomic<uint64>			mWakeTopologyRevision { 1 };
	uint64					mBuiltWakeTopologyRevision = 0;
	uint64					mWakeMotionTypeSequence = ~uint64(0);
	uint64					mWakeActiveBodiesRevision = ~uint64(0);
	uint32					mWakeEligibleConstraints = 0;
	Array<WakeNode>			mWakeNodes;
	UnorderedMap<uint32, uint32> mWakeNodeByBody;
	UnorderedMap<const LegacyConstraintGroup *, uint32> mWakeFirstNodeByGroup;
	Array<BodyID>			mWakeActivationBatch;

#ifdef JPH_ENABLE_ASSERTS
	PhysicsLockContext		mLockContext;
#endif // JPH_ENABLE_ASSERTS
	Constraints				mConstraints;
	mutable Mutex			mConstraintsMutex;
};

JPH_NAMESPACE_END
