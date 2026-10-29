// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-License-Identifier: MIT

#pragma once

#include <Jolt/Physics/Collision/CollidePointResult.h>
#include <Jolt/Physics/Collision/CollideShape.h>
#include <Jolt/Physics/Collision/Shape/BoxShape.h>
#include <Jolt/Physics/Collision/Shape/ConvexHullShape.h>
#include <Jolt/Physics/Collision/TransformedShape.h>

JPH_NAMESPACE_BEGIN

namespace CompoundSurfaceContactDetail
{
	inline TransformedShape GetLeaf(const Shape *inShape, SubShapeID inID)
	{
		TransformedShape leaf(RVec3::sZero(), Quat::sIdentity(), inShape, BodyID());
		while (leaf.mShape->GetType() == EShapeType::Compound || leaf.mShape->GetType() == EShapeType::Decorated)
		{
			SubShapeID remainder;
			leaf = leaf.mShape->GetSubShapeTransformedShape(inID, Vec3(leaf.mShapePositionCOM),
				leaf.mShapeRotation, leaf.GetShapeScale(), remainder);
			inID = remainder;
		}
		return leaf;
	}
}

/// Detect a contact on an internal convex boundary of a static compound.
/// The caller must opt in with its geometry rounding tolerance. Unsupported
/// leaves and exposed surfaces retain their original collision response.
inline bool IsContactBuriedInCompound(const Shape *inMoving, const Shape *inTarget,
	Mat44Arg inTargetTransform, const CollideShapeResult &inHit, float inTolerance, const ShapeFilter &inFilter)
{
	if (inTolerance <= 0.0f || inTarget->GetType() != EShapeType::Compound)
		return false;
	Vec3 normal = inHit.mPenetrationAxis.NormalizedOr(Vec3::sZero());
	if (normal.IsNearZero())
		return false;
	// Probe beyond the rounded surface, then require a full radius of clearance
	// in a different leaf. A mere touch at a floor/wall corner is insufficient.
	Vec3 point = inTargetTransform.InversedRotationTranslation() * (inHit.mContactPointOn2 - 2.0f * inTolerance * normal);
	TransformedShape moving = CompoundSurfaceContactDetail::GetLeaf(inMoving, inHit.mSubShapeID1);
	class LeafFilter final : public ShapeFilter
	{
	public:
		LeafFilter(const Shape *inMovingLeaf, const SubShapeID &inID, const ShapeFilter &inFilter) :
			mMoving(inMovingLeaf), mID(inID), mFilter(inFilter) { }
		virtual bool ShouldCollide(const Shape *inShape, const SubShapeID &inID) const override
		{
			return mFilter.ShouldCollide(mMoving, mID, inShape, inID);
		}
		const Shape *mMoving;
		SubShapeID mID;
		const ShapeFilter &mFilter;
	} filter(moving.mShape, inHit.mSubShapeID1, inFilter);
	class InteriorCollector final : public CollidePointCollector
	{
	public:
		InteriorCollector(const Shape *inTarget, const SubShapeID &inID, Vec3Arg inPoint, float inMargin) :
			mTarget(inTarget), mID(inID), mPoint(inPoint), mMargin(inMargin) { }
		virtual void AddHit(const CollidePointResult &inPointHit) override
		{
			if (inPointHit.mSubShapeID2 == mID)
				return;
			TransformedShape leaf = CompoundSurfaceContactDetail::GetLeaf(mTarget, inPointHit.mSubShapeID2);
			const Shape *shape = leaf.mShape;
			Vec3 scale = leaf.GetShapeScale();
			Vec3 local = Vec3(leaf.GetInverseCenterOfMassTransform() * RVec3(mPoint)) / scale;
			if (shape->GetSubType() == EShapeSubType::ConvexHull)
			{
				for (const Plane &plane : static_cast<const ConvexHullShape *>(shape)->GetPlanes())
					if (plane.SignedDistance(local) >= -mMargin * (plane.GetNormal() / scale).Length())
						return;
			}
			else if (shape->GetSubType() == EShapeSubType::Box)
			{
				Vec3 clearance = static_cast<const BoxShape *>(shape)->GetHalfExtent() - local.Abs();
				if ((clearance * scale.Abs()).ReduceMin() <= mMargin)
					return;
			}
			else
				return;
			mBuried = true;
			ForceEarlyOut();
		}
		const Shape *mTarget;
		SubShapeID mID;
		Vec3 mPoint;
		float mMargin;
		bool mBuried = false;
	} interior(inTarget, inHit.mSubShapeID2, point, inTolerance);
	inTarget->CollidePoint(point, SubShapeIDCreator(), interior, filter);
	return interior.mBuried;
}

JPH_NAMESPACE_END
