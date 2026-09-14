// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-FileCopyrightText: 2021 Jorrit Rouwe
// SPDX-License-Identifier: MIT

#include <Jolt/Jolt.h>

#include <Jolt/Physics/Constraints/ConstraintManager.h>
#include <Jolt/Physics/Constraints/CalculateSolverSteps.h>
#include <Jolt/Physics/Constraints/TwoBodyConstraint.h>
#include <Jolt/Physics/Constraints/PathConstraint.h>
#include <Jolt/Physics/Body/BodyManager.h>
#include <Jolt/Physics/IslandBuilder.h>
#include <Jolt/Physics/StateRecorder.h>
#include <Jolt/Physics/PhysicsLock.h>
#include <Jolt/Core/Profiler.h>
#include <Jolt/Core/QuickSort.h>
#include <Jolt/Core/UnorderedMap.h>
#include <span>
#include <memory>


JPH_NAMESPACE_BEGIN

ConstraintManager::~ConstraintManager()
{
	// A caller can retain constraints after destroying a PhysicsSystem.
	for (const Ref<Constraint> &constraint : mConstraints)
		constraint->mConstraintManager = nullptr;
}

void ConstraintManager::Add(Constraint **inConstraints, int inNumber)
{
	UniqueLock lock(mConstraintsMutex JPH_IF_ENABLE_ASSERTS(, mLockContext, EPhysicsLockTypes::ConstraintsList));

	// One-at-a-time Source creation must not reallocate/copy every existing Ref.
	const size_t required = mConstraints.size() + inNumber;
	if (required > mConstraints.capacity())
		mConstraints.reserve(max(required, max(size_t(8), mConstraints.capacity() + mConstraints.capacity() / 2)));
	InvalidateWakeTopology();

	for (Constraint **c = inConstraints, **c_end = inConstraints + inNumber; c < c_end; ++c)
	{
		Constraint *constraint = *c;

		// Assume this constraint has not been added yet
		JPH_ASSERT(constraint->mConstraintIndex == Constraint::cInvalidConstraintIndex);

		// Add to the list
		constraint->mConstraintIndex = uint32(mConstraints.size());
		constraint->mConstraintManager = this;
		mConstraints.push_back(constraint);
		mNumLegacyPreparationConstraints += constraint->CanPrepareLegacyVelocityConstraint()? 1 : 0;
	}
}

void ConstraintManager::Remove(Constraint **inConstraints, int inNumber)
{
	UniqueLock lock(mConstraintsMutex JPH_IF_ENABLE_ASSERTS(, mLockContext, EPhysicsLockTypes::ConstraintsList));

	for (Constraint **c = inConstraints, **c_end = inConstraints + inNumber; c < c_end; ++c)
	{
		Constraint *constraint = *c;

		// Reset constraint index for this constraint
		if (constraint->CanPrepareLegacyVelocityConstraint())
		{
			JPH_ASSERT(mNumLegacyPreparationConstraints > 0);
			--mNumLegacyPreparationConstraints;
		}
		uint32 this_constraint_idx = constraint->mConstraintIndex;
		constraint->mConstraintIndex = Constraint::cInvalidConstraintIndex;
		constraint->mConstraintManager = nullptr;
		InvalidateWakeTopology();
		JPH_ASSERT(this_constraint_idx != Constraint::cInvalidConstraintIndex);

		// Check if this constraint is somewhere in the middle of the constraints, in this case we need to move the last constraint to this position
		uint32 last_constraint_idx = uint32(mConstraints.size() - 1);
		if (this_constraint_idx < last_constraint_idx)
		{
			Constraint *last_constraint = mConstraints[last_constraint_idx];
			last_constraint->mConstraintIndex = this_constraint_idx;
			mConstraints[this_constraint_idx] = last_constraint;
		}

		// Pop last constraint
		mConstraints.pop_back();
	}
}

Constraints ConstraintManager::GetConstraints() const
{
	UniqueLock lock(mConstraintsMutex JPH_IF_ENABLE_ASSERTS(, mLockContext, EPhysicsLockTypes::ConstraintsList));

	Constraints copy = mConstraints;
	return copy;
}

namespace
{
	bool IsWakeLinkEnabled(const Constraint &inConstraint)
	{
		if (inConstraint.GetType() != EConstraintType::TwoBodyConstraint || !inConstraint.GetEnabled())
			return false;
		return inConstraint.GetSubType() != EConstraintSubType::Path
			|| static_cast<const PathConstraint &>(inConstraint).GetPath() != nullptr;
	}

	const LegacyConstraintGroup *GetWakeGroup(const Constraint &inConstraint)
	{
		const LegacyConstraintGroup *group = inConstraint.GetLegacyConstraintGroup();
		// A private single-constraint group adds nothing to its ordinary dynamic edge.
		// Read the live count so 1 -> 2 -> 1 membership changes take effect this step.
		return group != nullptr && group->GetConstraintCount() > 1? group : nullptr;
	}
}

bool ConstraintManager::NeedsWakePropagation(const BodyManager &inBodyManager) const
{
	return inBodyManager.GetNumActiveBodies(EBodyType::RigidBody) != 0
		&& (mBuiltWakeTopologyRevision != mWakeTopologyRevision.load(memory_order_acquire)
			|| mWakeMotionTypeSequence != Body::GetMotionTypeChangeSequence()
			|| mWakeActiveBodiesRevision != inBodyManager.GetActiveBodiesRevision());
}

ConstraintManager::WakeStatistics ConstraintManager::PropagateWake(BodyManager &ioBodyManager)
{
	JPH_PROFILE_FUNCTION();
	WakeStatistics statistics;
	if (!NeedsWakePropagation(ioBodyManager))
		return statistics;

	const uint64 topology_revision = mWakeTopologyRevision.load(memory_order_acquire);
	const uint64 motion_sequence = Body::GetMotionTypeChangeSequence();
	if (mBuiltWakeTopologyRevision != topology_revision || mWakeMotionTypeSequence != motion_sequence)
	{
		mWakeNodes.clear();
		mWakeNodeByBody.clear();
		mWakeFirstNodeByGroup.clear();
		mWakeEligibleConstraints = 0;
		statistics.mScannedConstraints = uint32(mConstraints.size());
		for (const Ref<Constraint> &constraint : mConstraints)
		{
			if (!IsWakeLinkEnabled(*constraint))
				continue;
			const TwoBodyConstraint &two_body = static_cast<const TwoBodyConstraint &>(*constraint);
			const Body *body1 = two_body.GetBody1(), *body2 = two_body.GetBody2();
			bool eligible = false;
			if (body1->IsDynamic() && body2->IsDynamic())
			{
				const uint32 node1 = GetOrCreateWakeNode(*body1);
				const uint32 node2 = GetOrCreateWakeNode(*body2);
				LinkWakeNodes(node1, node2);
				eligible = true;
			}
			if (const LegacyConstraintGroup *group = GetWakeGroup(*constraint))
			{
				for (const Body *body : { body1, body2 })
					if (body->IsDynamic())
					{
						const uint32 node = GetOrCreateWakeNode(*body);
						auto [entry, inserted] = mWakeFirstNodeByGroup.try_emplace(group, node);
						if (!inserted)
							LinkWakeNodes(entry->second, node);
						eligible = true;
					}
			}
			if (eligible)
				++mWakeEligibleConstraints;
		}
		// Flatten once: subsequent wake events never traverse constraints or hash maps.
		for (uint32 i = 0; i < uint32(mWakeNodes.size()); ++i)
			mWakeNodes[i].mParent = FindWakeRoot(i);
		mBuiltWakeTopologyRevision = topology_revision;
		mWakeMotionTypeSequence = motion_sequence;
	}

	for (WakeNode &node : mWakeNodes)
		node.mHasActiveBody = node.mHasSleepingBody = false;
	for (const WakeNode &node : mWakeNodes)
	{
		const Body &body = ioBodyManager.GetBody(node.mBodyID);
		WakeNode &root = mWakeNodes[node.mParent];
		root.mHasActiveBody |= body.IsActive();
		root.mHasSleepingBody |= !body.IsActive();
	}
	mWakeActivationBatch.clear();
	for (uint32 i = 0; i < uint32(mWakeNodes.size()); ++i)
	{
		const WakeNode &node = mWakeNodes[i];
		const WakeNode &root = mWakeNodes[node.mParent];
		if (root.mHasActiveBody && root.mHasSleepingBody)
		{
			mWakeActivationBatch.push_back(node.mBodyID);
			if (node.mParent == i)
				++statistics.mMixedComponents;
		}
	}
	statistics.mEligibleConstraints = mWakeEligibleConstraints;
	statistics.mDynamicBodies = uint32(mWakeNodes.size());
	statistics.mActivationBatchBodies = uint32(mWakeActivationBatch.size());
	const uint32 active_before = ioBodyManager.GetNumActiveBodies(EBodyType::RigidBody);
	// Reset the whole mixed component's sleep timers, but notify only actual wakes.
	ioBodyManager.ActivateBodies(mWakeActivationBatch.data(), int(mWakeActivationBatch.size()));
	statistics.mNewlyActivatedBodies = ioBodyManager.GetNumActiveBodies(EBodyType::RigidBody) - active_before;
	mWakeActiveBodiesRevision = ioBodyManager.GetActiveBodiesRevision();
	return statistics;
}

uint32 ConstraintManager::GetOrCreateWakeNode(const Body &inBody)
{
	const uint32 index = uint32(mWakeNodes.size());
	auto [entry, inserted] = mWakeNodeByBody.try_emplace(inBody.GetID().GetIndexAndSequenceNumber(), index);
	if (inserted)
		mWakeNodes.push_back({ inBody.GetID(), index, 1, false, false });
	return entry->second;
}

uint32 ConstraintManager::FindWakeRoot(uint32 inNode)
{
	uint32 root = inNode;
	while (mWakeNodes[root].mParent != root)
		root = mWakeNodes[root].mParent;
	while (mWakeNodes[inNode].mParent != inNode)
	{
		const uint32 parent = mWakeNodes[inNode].mParent;
		mWakeNodes[inNode].mParent = root;
		inNode = parent;
	}
	return root;
}

void ConstraintManager::LinkWakeNodes(uint32 inNode1, uint32 inNode2)
{
	uint32 root1 = FindWakeRoot(inNode1), root2 = FindWakeRoot(inNode2);
	if (root1 == root2)
		return;
	if (mWakeNodes[root1].mSize < mWakeNodes[root2].mSize)
		std::swap(root1, root2);
	WakeNode &parent = mWakeNodes[root1], &child = mWakeNodes[root2];
	child.mParent = root1;
	parent.mSize += child.mSize;
}

void ConstraintManager::GetActiveConstraints(uint32 inStartConstraintIdx, uint32 inEndConstraintIdx, Constraint **outActiveConstraints, uint32 &outNumActiveConstraints) const
{
	JPH_PROFILE_FUNCTION();

	JPH_ASSERT(inEndConstraintIdx <= mConstraints.size());

	uint32 num_active_constraints = 0;
	for (uint32 constraint_idx = inStartConstraintIdx; constraint_idx < inEndConstraintIdx; ++constraint_idx)
	{
		Constraint *c = mConstraints[constraint_idx];
		JPH_ASSERT(c->mConstraintIndex == constraint_idx);
		if (c->IsActive())
		{
			*(outActiveConstraints++) = c;
			num_active_constraints++;
		}
	}

	outNumActiveConstraints = num_active_constraints;
}

void ConstraintManager::sBuildIslands(Constraint **inActiveConstraints, uint32 inNumActiveConstraints, IslandBuilder &ioBuilder, BodyManager &inBodyManager)
{
	JPH_PROFILE_FUNCTION();

	for (uint32 constraint_idx = 0; constraint_idx < inNumActiveConstraints; ++constraint_idx)
	{
		Constraint *c = inActiveConstraints[constraint_idx];
		c->BuildIslands(constraint_idx, ioBuilder, inBodyManager);
	}
}

void ConstraintManager::sSortConstraints(Constraint **inActiveConstraints, uint32 *inConstraintIdxBegin, uint32 *inConstraintIdxEnd)
{
	JPH_PROFILE_FUNCTION();

	QuickSort(inConstraintIdxBegin, inConstraintIdxEnd, [inActiveConstraints](uint32 inLHS, uint32 inRHS) {
		const Constraint *lhs = inActiveConstraints[inLHS];
		const Constraint *rhs = inActiveConstraints[inRHS];

		if (lhs->GetConstraintPriority() != rhs->GetConstraintPriority())
			return lhs->GetConstraintPriority() < rhs->GetConstraintPriority();

		return lhs->mConstraintIndex < rhs->mConstraintIndex;
	});
}

void ConstraintManager::sSetupVelocityConstraints(Constraint **inActiveConstraints, uint32 inNumActiveConstraints, float inDeltaTime, bool inPrepareLegacy)
{
	JPH_PROFILE_FUNCTION();

	for (Constraint **c = inActiveConstraints, **c_end = inActiveConstraints + inNumActiveConstraints; c < c_end; ++c)
	{
		(*c)->SetupVelocityConstraint(inDeltaTime);
		(*c)->mLegacyPostContactSolveInIsland = false;
		if (inPrepareLegacy && (*c)->GetLegacyConstraintGroup() != nullptr && (*c)->CanPrepareLegacyVelocityConstraint())
		{
			Constraint *constraint = *c;
			LegacyConstraintGroup *group = constraint->GetLegacyConstraintGroup();
			const size_t index = constraint->mLegacyGroupIndex;
			JPH_ASSERT(!constraint->mLegacySolverStoragePrepared);
			JPH_ASSERT(index < group->mConstraints.size() && group->mConstraints[index] == constraint);
			const size_t offset = group->mSolverStorageOffsets[index];
			if (offset != SIZE_MAX)
			{
				void *storage = reinterpret_cast<std::byte *>(group->mSolverStorage.data()) + offset;
				constraint->InitializeLegacySolverStorage(storage);
				constraint->PrepareLegacyVelocityConstraint(storage, inDeltaTime);
				constraint->mLegacySolverStoragePrepared = true;
			}
		}
	}
}

void ConstraintManager::sReleasePreparedLegacyVelocityConstraints(Constraint **inActiveConstraints, const uint32 *inConstraintIdxBegin, const uint32 *inConstraintIdxEnd)
{
	for (const uint32 *index = inConstraintIdxBegin; index < inConstraintIdxEnd; ++index)
	{
		Constraint *constraint = inActiveConstraints[*index];
		if (!constraint->mLegacySolverStoragePrepared)
			continue;
		LegacyConstraintGroup *group = constraint->GetLegacyConstraintGroup();
		JPH_ASSERT(group != nullptr && group->mConstraints[constraint->mLegacyGroupIndex] == constraint);
		void *storage = reinterpret_cast<std::byte *>(group->mSolverStorage.data())
			+ group->mSolverStorageOffsets[constraint->mLegacyGroupIndex];
		constraint->DestroyLegacySolverStorage(storage);
		constraint->mLegacySolverStoragePrepared = false;
	}
}

template <class ConstraintCallback>
void ConstraintManager::sWarmStartVelocityConstraints(Constraint **inActiveConstraints, const uint32 *inConstraintIdxBegin, const uint32 *inConstraintIdxEnd, float inWarmStartImpulseRatio, ConstraintCallback &ioCallback)
{
	JPH_PROFILE_FUNCTION();

	for (const uint32 *constraint_idx = inConstraintIdxBegin; constraint_idx < inConstraintIdxEnd; ++constraint_idx)
	{
		Constraint *c = inActiveConstraints[*constraint_idx];
		ioCallback(c);
		if (c->GetLegacyConstraintGroup())
			c->WarmStartVelocityConstraint(0.0f);
		else
			c->WarmStartVelocityConstraint(inWarmStartImpulseRatio);
	}
}

// Specialize for the two constraint callback types
template void ConstraintManager::sWarmStartVelocityConstraints<CalculateSolverSteps>(Constraint **inActiveConstraints, const uint32 *inConstraintIdxBegin, const uint32 *inConstraintIdxEnd, float inWarmStartImpulseRatio, CalculateSolverSteps &ioCallback);
template void ConstraintManager::sWarmStartVelocityConstraints<DummyCalculateSolverSteps>(Constraint **inActiveConstraints, const uint32 *inConstraintIdxBegin, const uint32 *inConstraintIdxEnd, float inWarmStartImpulseRatio, DummyCalculateSolverSteps &ioCallback);

bool ConstraintManager::sSolveVelocityConstraints(Constraint **inActiveConstraints, const uint32 *inConstraintIdxBegin, const uint32 *inConstraintIdxEnd, float inDeltaTime)
{
	JPH_PROFILE_FUNCTION();

	bool any_impulse_applied = false;

	for (const uint32 *constraint_idx = inConstraintIdxBegin; constraint_idx < inConstraintIdxEnd; ++constraint_idx)
	{
		Constraint *c = inActiveConstraints[*constraint_idx];
		if (c->GetLegacyConstraintGroup() == nullptr)
			any_impulse_applied |= c->SolveVelocityConstraint(inDeltaTime);
	}

	return any_impulse_applied;
}

bool ConstraintManager::sHasLegacyConstraints(Constraint **inActiveConstraints, const uint32 *inConstraintIdxBegin, const uint32 *inConstraintIdxEnd)
{
	for (const uint32 *constraint_idx = inConstraintIdxBegin; constraint_idx < inConstraintIdxEnd; ++constraint_idx)
		if (inActiveConstraints[*constraint_idx]->GetLegacyConstraintGroup() != nullptr)
			return true;
	return false;
}

void ConstraintManager::LegacyVelocityBatch::Group::AddMember(Constraint *inConstraint)
{
	if (mFirstMember == nullptr)
		mFirstMember = inConstraint;
	else
	{
		if (mMembers.empty())
			mMembers.push_back(mFirstMember);
		mMembers.push_back(inConstraint);
	}
	mPostContact &= inConstraint->RequiresPostContactLegacySolve();
}

std::span<Constraint *const> ConstraintManager::LegacyVelocityBatch::Group::GetMembers() const
{
	return mMembers.empty()? std::span<Constraint *const>(&mFirstMember, 1)
		: std::span<Constraint *const>(mMembers.data(), mMembers.size());
}

bool ConstraintManager::LegacyVelocityBatch::Group::Solve(float inDeltaTime)
{
	const auto members = GetMembers();
	for (Constraint *c : members)
		c->mLegacyPostContactSolveInIsland = mPostContact;
	return mGroup->SolveVelocityConstraints(members.data(), members.size(), inDeltaTime);
}

void ConstraintManager::LegacyVelocityBatch::BuildGroups(Constraint **inActiveConstraints, const uint32 *inBegin, const uint32 *inEnd)
{
	UnorderedMap<LegacyConstraintGroup *, size_t> group_indices;
	for (const uint32 *index = inBegin; index < inEnd; ++index)
	{
		Constraint *constraint = inActiveConstraints[*index];
		if (LegacyConstraintGroup *group = constraint->GetLegacyConstraintGroup())
		{
			auto [entry, inserted] = group_indices.try_emplace(group, mGroups.size());
			if (inserted)
				mGroups.push_back({ group });
			mGroups[entry->second].AddMember(constraint);
		}
	}
}

bool ConstraintManager::LegacyVelocityBatch::BuildLayers()
{
	if (mGroups.size() < 2)
		return false;
	std::vector<uint32> layers(mGroups.size());
	UnorderedMap<const Body *, uint32> next_body_layer;
	uint32 num_layers = 0;
	// Preserve the effective controller order, not merely group insertion order.
	for (bool post_contact : { false, true })
		for (uint32 index = 0; index < mGroups.size(); ++index)
		{
			const Group &group = mGroups[index];
			if (group.mPostContact != post_contact)
				continue;
			uint32 layer = 0;
			for (Constraint *constraint : group.GetMembers())
			{
				if (constraint->GetType() != EConstraintType::TwoBodyConstraint)
					return false;
				const auto *two_body = static_cast<const TwoBodyConstraint *>(constraint);
				for (const Body *body : { two_body->GetBody1(), two_body->GetBody2() })
					if (body->IsDynamic())
					{
						auto previous = next_body_layer.find(body);
						if (previous != next_body_layer.end())
							layer = max(layer, previous->second);
					}
			}
			// Read the complete group's dependencies before updating any endpoint, so
			// repeated bodies within the same group cannot create a self-dependency.
			for (Constraint *constraint : group.GetMembers())
			{
				const auto *two_body = static_cast<const TwoBodyConstraint *>(constraint);
				for (const Body *body : { two_body->GetBody1(), two_body->GetBody2() })
					if (body->IsDynamic())
						next_body_layer[body] = layer + 1;
			}
			layers[index] = layer;
			num_layers = max(num_layers, layer + 1);
		}
	if (num_layers == mGroups.size())
		return false; // A genuine chain has no parallel width; keep its serial fast path.
	mLayerEnds.resize(num_layers, 0);
	for (uint32 layer : layers)
		++mLayerEnds[layer];
	for (uint32 layer = 1; layer < num_layers; ++layer)
		mLayerEnds[layer] += mLayerEnds[layer - 1];
	std::vector<uint32> next(num_layers, 0);
	for (uint32 layer = 1; layer < num_layers; ++layer)
		next[layer] = mLayerEnds[layer - 1];
	mGroupIndices.resize(mGroups.size());
	for (bool post_contact : { false, true })
		for (uint32 index = 0; index < mGroups.size(); ++index)
			if (mGroups[index].mPostContact == post_contact)
				mGroupIndices[next[layers[index]]++] = index;
	return true;
}

bool ConstraintManager::LegacyVelocityBatch::SolveSerial()
{
	bool applied = false;
	for (bool post_contact : { false, true })
		for (Group &group : mGroups)
			if (group.mPostContact == post_contact)
				applied |= group.Solve(mDeltaTime);
	return applied;
}

bool ConstraintManager::LegacyVelocityBatch::TryExecute()
{
	uint64 cursor;
	uint32 layer, index;
	for (;;)
	{
		cursor = mCursor.load(memory_order_acquire);
		if (cursor == UINT64_MAX)
			return false;
		layer = uint32(cursor >> 32);
		index = uint32(cursor);
		if (index >= mLayerEnds[layer])
			return false; // Claimed groups are still running; try another island.
		if (mCursor.compare_exchange_weak(cursor, cursor + 1, memory_order_acquire, memory_order_relaxed))
			break;
	}
	JPH_PROFILE("HavanaGroup");
	if (mGroups[mGroupIndices[index]].Solve(mDeltaTime))
		mAppliedImpulse.store(true, memory_order_relaxed);
	const uint32 layer_begin = layer == 0? 0 : mLayerEnds[layer - 1];
	if (mProcessed.fetch_add(1, memory_order_acq_rel) + 1 == mLayerEnds[layer] - layer_begin)
	{
		if (layer + 1 < mLayerEnds.size())
		{
			mProcessed.store(0, memory_order_relaxed);
			// Publish all this layer's velocity writes before any dependent group runs.
			mCursor.store((uint64(layer + 1) << 32) | mLayerEnds[layer], memory_order_release);
		}
		else
		{
			// Only this last worker resumes contacts and releases prepared storage.
			// Keep the consumed cursor until it returns: no second worker may resume it.
			mContinuation(mAppliedImpulse.load(memory_order_relaxed));
			mCursor.store(UINT64_MAX, memory_order_release);
		}
	}
	return true;
}

void ConstraintManager::LegacyVelocityScheduler::Prepare(const IslandBuilder &inBuilder, Constraint **inActiveConstraints, bool inEnabled)
{
	Clear();
	mEnabled = false;
	mNumIslands = 0;
	// Empty islands form a suffix in Jolt's sorted island list. They cannot publish work.
	for (; mNumIslands < inBuilder.GetNumIslands(); ++mNumIslands)
	{
		uint32 *cb, *ce, *pb, *pe;
		const bool constraints = inBuilder.GetConstraintsInIsland(mNumIslands, cb, ce);
		const bool contacts = inBuilder.GetContactsInIsland(mNumIslands, pb, pe);
		if (!constraints && !contacts)
			break;
		// The input only sizes publication slots here. Groups and their dependencies
		// are built by the owning solver after deterministic constraint sorting.
		if (!mEnabled && inEnabled && constraints && ce - cb >= 128)
			mEnabled = sHasLegacyConstraints(inActiveConstraints, cb, ce);
	}
	if (mEnabled)
		mBatches.resize(mNumIslands);
	mRemainingIslands.store(mEnabled? mNumIslands : 0, memory_order_relaxed);
}

void ConstraintManager::LegacyVelocityScheduler::Clear()
{
	JPH_ASSERT(!HasPendingIslands());
	mHead.store(nullptr, memory_order_relaxed);
	mBatches.clear();
}

void ConstraintManager::LegacyVelocityScheduler::Publish(uint32 inIsland, std::unique_ptr<LegacyVelocityBatch> inBatch)
{
	JPH_ASSERT(mEnabled && inIsland < mBatches.size() && mBatches[inIsland] == nullptr);
	LegacyVelocityBatch *batch = inBatch.get();
	mBatches[inIsland] = std::move(inBatch);
	LegacyVelocityBatch *head = mHead.load(memory_order_relaxed);
	do { batch->mNext = head; }
	while (!mHead.compare_exchange_weak(head, batch, memory_order_release, memory_order_relaxed));
}

bool ConstraintManager::LegacyVelocityScheduler::TryExecute()
{
	for (LegacyVelocityBatch *batch = mHead.load(memory_order_acquire); batch != nullptr; batch = batch->mNext)
	{
		if (batch->IsDone())
		{
			// Retire completed prefixes, but retain their storage until every solver job
			// returns. A worker that already read this pointer can safely finish its scan.
			LegacyVelocityBatch *expected = batch;
			mHead.compare_exchange_strong(expected, batch->mNext, memory_order_acq_rel, memory_order_relaxed);
		}
		else if (batch->TryExecute())
			return true;
	}
	return false;
}

bool ConstraintManager::sSolveLegacyVelocityConstraints(Constraint **inActiveConstraints, const uint32 *inBegin, const uint32 *inEnd, float inDeltaTime)
{
	JPH_PROFILE_FUNCTION();
	if (inEnd - inBegin == 1)
	{
		Constraint *constraint = inActiveConstraints[*inBegin];
		if (LegacyConstraintGroup *group = constraint->GetLegacyConstraintGroup())
		{
			constraint->mLegacyPostContactSolveInIsland = constraint->RequiresPostContactLegacySolve();
			return group->SolveVelocityConstraints(&constraint, 1, inDeltaTime);
		}
		return false;
	}
	LegacyVelocityBatch batch;
	batch.mDeltaTime = inDeltaTime;
	batch.BuildGroups(inActiveConstraints, inBegin, inEnd);
	return batch.SolveSerial();
}

bool ConstraintManager::sSolvePostContactLegacyVelocityMembers(Constraint **inActiveConstraints, const uint32 *inConstraintIdxBegin, const uint32 *inConstraintIdxEnd, float inDeltaTime)
{
	JPH_PROFILE_FUNCTION();

	bool any_impulse_applied = false;
	for (const uint32 *constraint_idx = inConstraintIdxBegin; constraint_idx < inConstraintIdxEnd; ++constraint_idx)
	{
		Constraint *c = inActiveConstraints[*constraint_idx];
		LegacyConstraintGroup *group = c->GetLegacyConstraintGroup();
		if (group == nullptr || c->RequiresPostContactSolveInIsland() || !c->RequiresPostContactLegacySolve())
			continue;

		// The LCS already completed its normal Havana pass. Correct only this predictive member
		// against the final contact velocity; running a complete mixed group here regressed Dukes.
		any_impulse_applied |= c->SolvePostContactVelocityConstraintInLegacyGroup(inDeltaTime);
	}
	return any_impulse_applied;
}

bool ConstraintManager::sSolvePositionConstraints(Constraint **inActiveConstraints, const uint32 *inConstraintIdxBegin, const uint32 *inConstraintIdxEnd, float inDeltaTime, float inBaumgarte)
{
	JPH_PROFILE_FUNCTION();

	bool any_impulse_applied = false;

	for (const uint32 *constraint_idx = inConstraintIdxBegin; constraint_idx < inConstraintIdxEnd; ++constraint_idx)
	{
		Constraint *c = inActiveConstraints[*constraint_idx];
		if (c->GetLegacyConstraintGroup())
		{
			if (!c->RequiresPostContactSolveInIsland())
				any_impulse_applied |= c->SolvePositionConstraintInLegacyGroup(inDeltaTime, inBaumgarte);
		}
		else
			any_impulse_applied |= c->SolvePositionConstraint(inDeltaTime, inBaumgarte);
	}

	return any_impulse_applied;
}

bool ConstraintManager::sSolvePostContactLegacyPositionConstraints(Constraint **inActiveConstraints, const uint32 *inConstraintIdxBegin, const uint32 *inConstraintIdxEnd, float inDeltaTime, float inBaumgarte)
{
	JPH_PROFILE_FUNCTION();

	bool any_impulse_applied = false;
	for (const uint32 *constraint_idx = inConstraintIdxBegin; constraint_idx < inConstraintIdxEnd; ++constraint_idx)
	{
		Constraint *c = inActiveConstraints[*constraint_idx];
		// Keep a mixed Source LCS in its ordinary velocity-controller stage, but still let an
		// individual predictive micro-angle member repair the pose change made by Jolt's contact
		// position pass. simfphys rear wheels live in the same Weld/Rope group as the chassis, so
		// classifying this correction by the whole group leaves their +/-0.01 degree limits open.
		if (c->GetLegacyConstraintGroup() != nullptr && c->RequiresPostContactLegacySolve())
			any_impulse_applied |= c->SolvePostContactPositionConstraintInLegacyGroup(inDeltaTime, inBaumgarte);
	}
	return any_impulse_applied;
}

#ifdef JPH_DEBUG_RENDERER
void ConstraintManager::DrawConstraints(DebugRenderer *inRenderer) const
{
	JPH_PROFILE_FUNCTION();

	UniqueLock lock(mConstraintsMutex JPH_IF_ENABLE_ASSERTS(, mLockContext, EPhysicsLockTypes::ConstraintsList));

	for (const Ref<Constraint> &c : mConstraints)
		c->DrawConstraint(inRenderer);
}

void ConstraintManager::DrawConstraintLimits(DebugRenderer *inRenderer) const
{
	JPH_PROFILE_FUNCTION();

	UniqueLock lock(mConstraintsMutex JPH_IF_ENABLE_ASSERTS(, mLockContext, EPhysicsLockTypes::ConstraintsList));

	for (const Ref<Constraint> &c : mConstraints)
		c->DrawConstraintLimits(inRenderer);
}

void ConstraintManager::DrawConstraintReferenceFrame(DebugRenderer *inRenderer) const
{
	JPH_PROFILE_FUNCTION();

	UniqueLock lock(mConstraintsMutex JPH_IF_ENABLE_ASSERTS(, mLockContext, EPhysicsLockTypes::ConstraintsList));

	for (const Ref<Constraint> &c : mConstraints)
		c->DrawConstraintReferenceFrame(inRenderer);
}
#endif // JPH_DEBUG_RENDERER

void ConstraintManager::SaveState(StateRecorder &inStream, const StateRecorderFilter *inFilter) const
{
	UniqueLock lock(mConstraintsMutex JPH_IF_ENABLE_ASSERTS(, mLockContext, EPhysicsLockTypes::ConstraintsList));

	// Write state of constraints
	if (inFilter != nullptr)
	{
		// Determine which constraints to save
		Array<Constraint *> constraints;
		constraints.reserve(mConstraints.size());
		for (const Ref<Constraint> &c : mConstraints)
			if (inFilter->ShouldSaveConstraint(*c))
				constraints.push_back(c);

		// Save them
		uint32 num_constraints = (uint32)constraints.size();
		inStream.Write(num_constraints);
		for (const Constraint *c : constraints)
		{
			inStream.Write(c->mConstraintIndex);
			c->SaveState(inStream);
		}
	}
	else
	{
		// Save all constraints
		uint32 num_constraints = (uint32)mConstraints.size();
		inStream.Write(num_constraints);
		for (const Ref<Constraint> &c : mConstraints)
		{
			inStream.Write(c->mConstraintIndex);
			c->SaveState(inStream);
		}
	}
}

bool ConstraintManager::RestoreState(StateRecorder &inStream)
{
	UniqueLock lock(mConstraintsMutex JPH_IF_ENABLE_ASSERTS(, mLockContext, EPhysicsLockTypes::ConstraintsList));

	if (inStream.IsValidating())
	{
		// Read state of constraints
		uint32 num_constraints = (uint32)mConstraints.size(); // Initialize to current value for validation
		inStream.Read(num_constraints);
		if (num_constraints != mConstraints.size())
		{
			JPH_ASSERT(false, "Cannot handle adding/removing constraints");
			return false;
		}
		for (const Ref<Constraint> &c : mConstraints)
		{
			uint32 constraint_index = c->mConstraintIndex;
			inStream.Read(constraint_index);
			if (constraint_index != c->mConstraintIndex)
			{
				JPH_ASSERT(false, "Unexpected constraint index");
				return false;
			}
			c->RestoreState(inStream);
		}
	}
	else
	{
		// Not validating, use more flexible reading, read number of constraints
		uint32 num_constraints = 0;
		inStream.Read(num_constraints);

		for (uint32 idx = 0; idx < num_constraints; ++idx)
		{
			uint32 constraint_index;
			inStream.Read(constraint_index);
			if (mConstraints.size() <= constraint_index)
			{
				JPH_ASSERT(false, "Restoring state for non-existing constraint");
				return false;
			}
			mConstraints[constraint_index]->RestoreState(inStream);
		}
	}

	return true;
}

JPH_NAMESPACE_END
