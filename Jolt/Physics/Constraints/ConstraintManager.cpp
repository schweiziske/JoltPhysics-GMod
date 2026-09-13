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
#include <unordered_map>
#include <memory>
#include <thread>

JPH_NAMESPACE_BEGIN

void ConstraintManager::Add(Constraint **inConstraints, int inNumber)
{
	UniqueLock lock(mConstraintsMutex JPH_IF_ENABLE_ASSERTS(, mLockContext, EPhysicsLockTypes::ConstraintsList));

	mConstraints.reserve(mConstraints.size() + inNumber);

	for (Constraint **c = inConstraints, **c_end = inConstraints + inNumber; c < c_end; ++c)
	{
		Constraint *constraint = *c;

		// Assume this constraint has not been added yet
		JPH_ASSERT(constraint->mConstraintIndex == Constraint::cInvalidConstraintIndex);

		// Add to the list
		constraint->mConstraintIndex = uint32(mConstraints.size());
		mConstraints.push_back(constraint);
	}
}

void ConstraintManager::Remove(Constraint **inConstraints, int inNumber)
{
	UniqueLock lock(mConstraintsMutex JPH_IF_ENABLE_ASSERTS(, mLockContext, EPhysicsLockTypes::ConstraintsList));

	for (Constraint **c = inConstraints, **c_end = inConstraints + inNumber; c < c_end; ++c)
	{
		Constraint *constraint = *c;

		// Reset constraint index for this constraint
		uint32 this_constraint_idx = constraint->mConstraintIndex;
		constraint->mConstraintIndex = Constraint::cInvalidConstraintIndex;
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

void ConstraintManager::PropagateWake(BodyManager &ioBodyManager)
{
	JPH_PROFILE_FUNCTION();

	// No active rigid body means no dynamic wake seed. Do not copy or scan the list.
	if (ioBodyManager.GetNumActiveBodies(EBodyType::RigidBody) == 0 || mConstraints.empty())
		return;

	constexpr uint8 cActive = 1, cSleeping = 2;
	mWakeGroupState.clear();
	bool has_wake_frontier = false;
	for (const Ref<Constraint> &constraint : mConstraints)
	{
		if (!IsWakeLinkEnabled(*constraint))
			continue;
		const TwoBodyConstraint &two_body = static_cast<const TwoBodyConstraint &>(*constraint);
		const Body *body1 = two_body.GetBody1(), *body2 = two_body.GetBody2();
		if (body1->IsDynamic() && body2->IsDynamic() && body1->IsActive() != body2->IsActive())
		{
			has_wake_frontier = true;
			break;
		}
		if (const LegacyConstraintGroup *group = GetWakeGroup(*constraint))
		{
			uint8 state = 0;
			if (body1->IsDynamic()) state |= body1->IsActive()? cActive : cSleeping;
			if (body2->IsDynamic()) state |= body2->IsActive()? cActive : cSleeping;
			if (state == 0)
				continue;
			uint8 &group_state = mWakeGroupState[group];
			group_state |= state;
			if (group_state == (cActive | cSleeping))
			{
				has_wake_frontier = true;
				break;
			}
		}
	}
	if (!has_wake_frontier)
		return;

	mWakeNodes.clear();
	mWakeNodeByBody.clear();
	mWakeFirstNodeByGroup.clear();
	mWakeActivationBatch.clear();
	for (const Ref<Constraint> &constraint : mConstraints)
	{
		// IsActive() would omit wholly sleeping links and only wake one wave per step.
		if (!IsWakeLinkEnabled(*constraint))
			continue;
		const TwoBodyConstraint &two_body = static_cast<const TwoBodyConstraint &>(*constraint);
		const Body *body1 = two_body.GetBody1(), *body2 = two_body.GetBody2();
		if (body1->IsDynamic() && body2->IsDynamic())
		{
			const uint32 node1 = GetOrCreateWakeNode(*body1);
			const uint32 node2 = GetOrCreateWakeNode(*body2);
			LinkWakeNodes(node1, node2);
		}
		if (const LegacyConstraintGroup *group = GetWakeGroup(*constraint))
		{
			// Group dynamic members share a wake state even when separately anchored.
			// Static/kinematic endpoints are neither seeds nor bridges between groups.
			auto link_group_body = [&](const Body &inBody)
			{
				if (!inBody.IsDynamic())
					return;
				const uint32 node = GetOrCreateWakeNode(inBody);
				auto [entry, inserted] = mWakeFirstNodeByGroup.try_emplace(group, node);
				if (!inserted)
					LinkWakeNodes(entry->second, node);
			};
			link_group_body(*body1);
			link_group_body(*body2);
		}
	}

	for (uint32 index = 0; index < uint32(mWakeNodes.size()); ++index)
	{
		const WakeNode &root = mWakeNodes[FindWakeRoot(index)];
		if (root.mHasActiveBody && root.mHasSleepingBody)
			mWakeActivationBatch.push_back(mWakeNodes[index].mBodyID);
	}
	// Include awake members to reset the component's sleep timers together. Activation
	// callbacks are emitted only for bodies that were actually sleeping.
	ioBodyManager.ActivateBodies(mWakeActivationBatch.data(), int(mWakeActivationBatch.size()));
}

uint32 ConstraintManager::GetOrCreateWakeNode(const Body &inBody)
{
	const uint32 index = uint32(mWakeNodes.size());
	auto [entry, inserted] = mWakeNodeByBody.try_emplace(inBody.GetID().GetIndexAndSequenceNumber(), index);
	if (inserted)
		mWakeNodes.push_back({ inBody.GetID(), index, 1, inBody.IsActive(), !inBody.IsActive() });
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
	parent.mHasActiveBody |= child.mHasActiveBody;
	parent.mHasSleepingBody |= child.mHasSleepingBody;
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

void ConstraintManager::sSetupVelocityConstraints(Constraint **inActiveConstraints, uint32 inNumActiveConstraints, float inDeltaTime)
{
	JPH_PROFILE_FUNCTION();

	for (Constraint **c = inActiveConstraints, **c_end = inActiveConstraints + inNumActiveConstraints; c < c_end; ++c)
	{
		(*c)->SetupVelocityConstraint(inDeltaTime);
		(*c)->mLegacyPostContactSolveInIsland = false;
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

bool ConstraintManager::sSolveLegacyVelocityConstraints(Constraint **inActiveConstraints, const uint32 *inConstraintIdxBegin, const uint32 *inConstraintIdxEnd, float inDeltaTime, JobSystem *inJobSystem, JobSystem::Barrier *inBarrier, atomic<uint> *ioJobBudget)
{
	JPH_PROFILE_FUNCTION();

	struct IslandGroup
	{
		LegacyConstraintGroup *mGroup;
		std::vector<Constraint *> mMembers;
		bool mPostContact = true;

		bool Solve(float inStepDeltaTime)
		{
			for (Constraint *c : mMembers)
				c->mLegacyPostContactSolveInIsland = mPostContact;
			return mGroup->SolveVelocityConstraints(mMembers.data(), mMembers.size(), inStepDeltaTime);
		}
	};
	std::vector<IslandGroup> groups;
	std::unordered_map<LegacyConstraintGroup *, size_t> group_indices;
	for (const uint32 *constraint_idx = inConstraintIdxBegin; constraint_idx < inConstraintIdxEnd; ++constraint_idx)
	{
		Constraint *c = inActiveConstraints[*constraint_idx];
		if (LegacyConstraintGroup *group = c->GetLegacyConstraintGroup())
		{
			auto [entry, inserted] = group_indices.try_emplace(group, groups.size());
			if (inserted)
				groups.push_back({ group, {}, true });
			IslandGroup &batch = groups[entry->second];
			batch.mMembers.push_back(c);
			batch.mPostContact &= c->RequiresPostContactLegacySolve();
		}
	}

	// A contact can join many vehicles into one Jolt island without creating a dependency
	// between their Havana controller passes. Keep each Source group whole (including its
	// insertion order, scratch storage and rollback-on-break boundary), and merge groups
	// that share a writable body. Native constraints and contacts still bracket this pass.
	// Small islands / direct callers retain the original serial path without scheduling cost.
	constexpr size_t cMinParallelConstraints = 128;
	if (inJobSystem != nullptr && inBarrier != nullptr && ioJobBudget != nullptr
		&& ioJobBudget->load(memory_order_relaxed) != 0 && inJobSystem->GetMaxConcurrency() > 1
		&& groups.size() > 1 && size_t(inConstraintIdxEnd - inConstraintIdxBegin) >= cMinParallelConstraints)
	{
		std::vector<size_t> parents(groups.size());
		for (size_t i = 0; i < parents.size(); ++i)
			parents[i] = i;
		auto root = [&](size_t index)
		{
			while (parents[index] != index)
			{
				parents[index] = parents[parents[index]];
				index = parents[index];
			}
			return index;
		};
		std::unordered_map<const Body *, size_t> body_owners;
		body_owners.reserve(size_t(inConstraintIdxEnd - inConstraintIdxBegin));
		bool known_dependencies = true;
		for (size_t i = 0; i < groups.size() && known_dependencies; ++i)
			for (Constraint *c : groups[i].mMembers)
			{
				if (c->GetType() != EConstraintType::TwoBodyConstraint)
				{
					known_dependencies = false;
					break;
				}
				const auto *two_body = static_cast<const TwoBodyConstraint *>(c);
				for (const Body *body : { two_body->GetBody1(), two_body->GetBody2() })
					if (!body->IsStatic())
					{
						// Include kinematic bodies: the existing break rollback restores their
						// velocities too. Sharing a read-only static/world anchor is safe.
						auto [entry, inserted] = body_owners.try_emplace(body, i);
						if (!inserted)
						{
							const size_t a = root(i), b = root(entry->second);
							parents[max(a, b)] = min(a, b);
						}
					}
			}
		if (known_dependencies)
		{
			std::vector<std::vector<size_t>> components;
			std::vector<size_t> component_index(groups.size(), SIZE_MAX);
			for (size_t i = 0; i < groups.size(); ++i)
			{
				size_t &index = component_index[root(i)];
				if (index == SIZE_MAX)
				{
					index = components.size();
					components.emplace_back();
				}
				components[index].push_back(i);
			}
			if (components.size() > 1)
			{
				struct ParallelBatch
				{
					std::vector<IslandGroup> mGroups;
					std::vector<std::vector<size_t>> mComponents;
					atomic<size_t> mNext { 0 }, mRemaining { 0 };
					atomic<bool> mAppliedImpulse { false };

					void Execute(float inStepDeltaTime, JobSystem &inJobs)
					{
						for (;;)
						{
							const size_t index = mNext.fetch_add(1, memory_order_relaxed);
							if (index >= mComponents.size())
								return;
							inJobs.CooperativeYield();
							bool applied = false;
							// Preserve both controller phases and first-encounter order for every
							// pair of groups that can observe each other's velocity writes.
							for (bool post_contact : { false, true })
								for (size_t group_index : mComponents[index])
									if (mGroups[group_index].mPostContact == post_contact)
										applied |= mGroups[group_index].Solve(inStepDeltaTime);
							if (applied)
								mAppliedImpulse.store(true, memory_order_relaxed);
							// Publish completion only after all body and solver-storage accesses.
							mRemaining.fetch_sub(1, memory_order_release);
						}
					}
				};
				auto batch = std::make_shared<ParallelBatch>();
				batch->mGroups = std::move(groups);
				batch->mComponents = std::move(components);
				batch->mRemaining.store(batch->mComponents.size(), memory_order_relaxed);
				// Avoid one tiny pool job per constraint/group. The owner also drains the
				// batch, so progress never depends on finding a free worker for a child job.
				const size_t num_helpers = min(size_t(7), min(size_t(inJobSystem->GetMaxConcurrency() - 1),
					min(batch->mComponents.size() - 1, size_t(inConstraintIdxEnd - inConstraintIdxBegin) / 64)));
				for (size_t i = 0; i < num_helpers; ++i)
				{
					uint remaining_budget = ioJobBudget->load(memory_order_relaxed);
					while (remaining_budget != 0 && !ioJobBudget->compare_exchange_weak(
						remaining_budget, remaining_budget - 1, memory_order_relaxed)) { }
					if (remaining_budget == 0)
						break;
					JobHandle helper = inJobSystem->CreateJob("SolveLegacyConstraintComponents", Color::sYellow,
						[batch, inDeltaTime, inJobSystem]()
						{
						#ifdef JPH_ENABLE_ASSERTS
							BodyAccess::Grant grant(BodyAccess::EAccess::ReadWrite, BodyAccess::EAccess::Read);
						#endif
							batch->Execute(inDeltaTime, *inJobSystem);
						}, 1);
					// Reuse Update's barrier, never nest WaitForJobs / enter a second owner
					// callback scheduler. A late helper owns its batch and only finds no work.
					inBarrier->AddJob(helper);
					helper.RemoveDependency();
				}
				batch->Execute(inDeltaTime, *inJobSystem);
				while (batch->mRemaining.load(memory_order_acquire) != 0)
				{
					inJobSystem->CooperativeYield();
					std::this_thread::yield();
				}
				return batch->mAppliedImpulse.load(memory_order_relaxed);
			}
		}
	}

	// A connected Gauss-Seidel chain must keep its velocity-write order. Geometry and
	// effective-mass preparation can share bodies, so batch that read-only work across
	// all Source groups in this island. Already-parallel independent components above
	// keep their complete group work and do not pay for a second scheduling layer.
	struct PreparedEntry
	{
		Constraint *mConstraint;
		void *mStorage;
	};
	struct PreparationBatch
	{
		std::vector<PreparedEntry> mEntries;
		atomic<size_t> mNext { 0 }, mRemaining { 0 };

		void Execute(float inStepDeltaTime)
		{
		#ifdef JPH_ENABLE_ASSERTS
			// The caller is already inside an island BodyAccess grant; temporarily narrow
			// its velocity rights instead of nesting Grant (which is forbidden).
			const BodyAccess::EAccess old_velocity_access = BodyAccess::sVelocityAccess();
			BodyAccess::sVelocityAccess() = BodyAccess::EAccess::None;
		#endif
			for (;;)
			{
				constexpr size_t cBatchSize = 16;
				const size_t begin = mNext.fetch_add(cBatchSize, memory_order_relaxed);
				if (begin >= mEntries.size())
					break;
				const size_t end = min(begin + cBatchSize, mEntries.size());
				for (size_t i = begin; i < end; ++i)
				{
					const PreparedEntry &entry = mEntries[i];
					entry.mConstraint->InitializeLegacySolverStorage(entry.mStorage);
					entry.mConstraint->PrepareLegacyVelocityConstraint(entry.mStorage, inStepDeltaTime);
					entry.mConstraint->mLegacySolverStoragePrepared = true;
				}
				// Publish all writes, including object initialization, before the original
				// group solve is allowed to borrow these storage ranges.
				mRemaining.fetch_sub(end - begin, memory_order_release);
			}
		#ifdef JPH_ENABLE_ASSERTS
			BodyAccess::sVelocityAccess() = old_velocity_access;
		#endif
		}
	};
	struct PreparationScope
	{
		std::shared_ptr<PreparationBatch> mBatch;
		~PreparationScope()
		{
			if (mBatch != nullptr)
				for (const PreparedEntry &entry : mBatch->mEntries)
				{
					// Cleanup belongs to this scope, not to the last (possibly late) helper.
					// It also covers members skipped by rollback-on-break or deactivation.
					JPH_ASSERT(entry.mConstraint->mLegacySolverStoragePrepared);
					entry.mConstraint->DestroyLegacySolverStorage(entry.mStorage);
					entry.mConstraint->mLegacySolverStoragePrepared = false;
				}
		}
	} prepared;
	constexpr size_t cMinParallelPreparationConstraints = 64;
	if (inJobSystem != nullptr && inBarrier != nullptr && ioJobBudget != nullptr
		&& inJobSystem->GetMaxConcurrency() > 1 && ioJobBudget->load(memory_order_relaxed) != 0
		&& size_t(inConstraintIdxEnd - inConstraintIdxBegin) >= cMinParallelPreparationConstraints)
	{
		auto batch = std::make_shared<PreparationBatch>();
		batch->mEntries.reserve(size_t(inConstraintIdxEnd - inConstraintIdxBegin));
		for (const IslandGroup &group : groups)
			for (Constraint *c : group.mMembers)
				if (c->IsActive() && c->CanPrepareLegacyVelocityConstraint())
				{
					const size_t index = c->mLegacyGroupIndex;
					JPH_ASSERT(!c->mLegacySolverStoragePrepared);
					JPH_ASSERT(index < group.mGroup->mConstraints.size() && group.mGroup->mConstraints[index] == c);
					const size_t offset = group.mGroup->mSolverStorageOffsets[index];
					if (offset != SIZE_MAX)
						batch->mEntries.push_back({ c, reinterpret_cast<std::byte *>(group.mGroup->mSolverStorage.data()) + offset });
				}
		if (batch->mEntries.size() >= cMinParallelPreparationConstraints)
		{
			batch->mRemaining.store(batch->mEntries.size(), memory_order_relaxed);
			const size_t num_helpers = min(size_t(7), min(size_t(inJobSystem->GetMaxConcurrency() - 1), batch->mEntries.size() / 32 - 1));
			for (size_t i = 0; i < num_helpers; ++i)
			{
				uint remaining_budget = ioJobBudget->load(memory_order_relaxed);
				while (remaining_budget != 0 && !ioJobBudget->compare_exchange_weak(
					remaining_budget, remaining_budget - 1, memory_order_relaxed)) { }
				if (remaining_budget == 0)
					break;
				JobHandle helper = inJobSystem->CreateJob("PrepareLegacyConstraintStorage", Color::sYellow,
					[batch, inDeltaTime]()
					{
					#ifdef JPH_ENABLE_ASSERTS
						BodyAccess::Grant grant(BodyAccess::EAccess::None, BodyAccess::EAccess::Read);
					#endif
						batch->Execute(inDeltaTime);
					}, 1);
				// Share Update's existing barrier and cumulative helper budget. A late
				// helper owns its batch and cannot access storage after the caller drains it.
				inBarrier->AddJob(helper);
				helper.RemoveDependency();
			}
			batch->Execute(inDeltaTime);
			while (batch->mRemaining.load(memory_order_acquire) != 0)
			{
				// All unclaimed work was drained by the caller. Only running, callback-free
				// preparation can remain; it cannot depend on servicing an owner callback.
				// Do not add a callback safe point that could invalidate prepared poses.
				std::this_thread::yield();
			}
			prepared.mBatch = std::move(batch);
		}
	}
	// Preserve the existing two controller phases and first-encounter group order. The hash
	// table is only an index; neither addresses nor hash iteration determine solve order.
	bool any_impulse_applied = false;
	for (bool post_contact : { false, true })
		for (IslandGroup &batch : groups)
			if (batch.mPostContact == post_contact)
				any_impulse_applied |= batch.Solve(inDeltaTime);
	return any_impulse_applied;
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
