// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-License-Identifier: MIT

#pragma once

#include <Jolt/Physics/Body/Body.h>
#include <Jolt/Physics/Collision/ManifoldBetweenTwoFaces.h>
#include <Jolt/Physics/Collision/TransformedShape.h>
#include <Jolt/Physics/Collision/Shape/SphereShape.h>

JPH_NAMESPACE_BEGIN

/// Local updates of an existing contact feature. Failure means that the caller must
/// rediscover this body pair, not keep using the old points. No history is discarded.
namespace PersistentContactManifold
{
struct Leaf
{
	const Shape *mShape;
	Mat44 mTransform;
	Vec3 mScale = Vec3::sReplicate(1.0f);
	SubShapeID mID;
};

inline Leaf GetLeaf(const Body &inBody, const SubShapeID &inID, RVec3Arg inBaseOffset)
{
	Leaf leaf { inBody.GetShape(), inBody.GetCenterOfMassTransform().PostTranslated(-inBaseOffset).ToMat44(), Vec3::sReplicate(1.0f), inID };
	while (leaf.mShape->GetType() == EShapeType::Compound || leaf.mShape->GetType() == EShapeType::Decorated)
	{
		SubShapeID remainder;
		TransformedShape child = leaf.mShape->GetSubShapeTransformedShape(leaf.mID,
			leaf.mTransform.GetTranslation(), leaf.mTransform.GetQuaternion(), leaf.mScale, remainder);
		leaf = { child.mShape, child.GetCenterOfMassTransform().ToMat44(), child.GetShapeScale(), remainder };
	}
	return leaf;
}

inline void GetFace(const Leaf &inLeaf, Vec3Arg inOutwardNormal, Shape::SupportingFace &outFace)
{
	// Jolt's supporting-face API takes the direction towards the shape.
	inLeaf.mShape->GetSupportingFace(inLeaf.mID, -inLeaf.mTransform.Multiply3x3Transposed(inOutwardNormal),
		inLeaf.mScale, inLeaf.mTransform, outFace);
}

inline bool IsInsideFace(const Shape::SupportingFace &inFace, Vec3Arg inNormal, Vec3Arg inPoint)
{
	// Keep boundary contacts on the narrow phase path, where active-edge rules and
	// adjacent triangles can select the correct feature and contact normal.
	constexpr float cBoundaryMargin = 1.0e-5f;
	for (uint i = 0; i < inFace.size(); ++i)
	{
		Vec3 edge = inFace[(i + 1) % inFace.size()] - inFace[i];
		if (edge.Cross(inPoint - inFace[i]).Dot(inNormal) <= cBoundaryMargin * edge.Length())
			return false;
	}
	return true;
}

inline bool RefreshFace(const Leaf &inReference, const Leaf &inIncident, Vec3Arg inExpectedNormal,
	bool inReferenceIsBody1, float inMaxSeparation, ContactManifold &ioManifold)
{
	if (inIncident.mShape->GetType() != EShapeType::Convex)
		return false;
	Shape::SupportingFace face;
	GetFace(inReference, inExpectedNormal, face);
	if (face.size() < 3)
		return false;
	Vec3 normal = (face[1] - face[0]).Cross(face[2] - face[0]).NormalizedOr(Vec3::sZero());
	if (normal.Dot(inExpectedNormal) < 1.0f - 1.0e-5f)
		return false; // Old contact was not with this face (e.g. an edge or corner).
	const ContactPoints &old_reference = inReferenceIsBody1? ioManifold.mRelativeContactPointsOn1 : ioManifold.mRelativeContactPointsOn2;
	for (Vec3 p : old_reference)
		if (!IsInsideFace(face, normal, p - (p - face[0]).Dot(normal) * normal))
			return false; // A reduced manifold may also contain points from adjacent faces.

	ConvexShape::SupportBuffer buffer;
	const auto *convex = static_cast<const ConvexShape *>(inIncident.mShape);
	const auto *support = convex->GetSupportFunction(ConvexShape::ESupportMode::IncludeConvexRadius, buffer, inIncident.mScale);
	Vec3 point = inIncident.mTransform * support->GetSupport(-inIncident.mTransform.Multiply3x3Transposed(normal));
	float gap = (point - face[0]).Dot(normal);
	if (abs(gap) > inMaxSeparation || !IsInsideFace(face, normal, point - gap * normal))
		return false;

	Shape::SupportingFace incident_face;
	GetFace(inIncident, -normal, incident_face);
	if (incident_face.empty())
		incident_face.push_back(point); // Sphere, capsule cap, or another smooth support point.
	ContactPoints on_reference, on_incident;
	for (Vec3 p : incident_face)
	{
		float separation = (p - face[0]).Dot(normal);
		Vec3 projected = p - separation * normal;
		if (!IsInsideFace(face, normal, projected))
			return false; // Clipping could change the feature; rediscover the pair.
		if (separation <= inMaxSeparation)
		{
			on_reference.push_back(projected);
			on_incident.push_back(p);
		}
	}
	if (on_reference.empty())
		return false;
	if (on_reference.size() > 4)
		PruneContactPoints(normal, on_reference, on_incident JPH_IF_DEBUG_RENDERER(, ioManifold.mBaseOffset));
	ioManifold.mWorldSpaceNormal = inReferenceIsBody1? normal : -normal;
	ioManifold.mRelativeContactPointsOn1 = inReferenceIsBody1? on_reference : on_incident;
	ioManifold.mRelativeContactPointsOn2 = inReferenceIsBody1? on_incident : on_reference;
	ioManifold.mPenetrationDepth = -FLT_MAX;
	for (uint i = 0; i < on_reference.size(); ++i)
		ioManifold.mPenetrationDepth = max(ioManifold.mPenetrationDepth, (on_reference[i] - on_incident[i]).Dot(normal));
	return true;
}

inline bool Refresh(const Body &inBody1, const Body &inBody2, Vec3Arg inNormalFollowing1,
	Vec3Arg inNormalFollowing2, float inMaxSeparation, ContactManifold &ioManifold)
{
	ioManifold.mBaseOffset = inBody1.GetCenterOfMassPosition();
	Leaf leaf1 = GetLeaf(inBody1, ioManifold.mSubShapeID1, ioManifold.mBaseOffset);
	Leaf leaf2 = GetLeaf(inBody2, ioManifold.mSubShapeID2, ioManifold.mBaseOffset);
	if (leaf1.mShape->GetSubType() == EShapeSubType::Sphere && leaf2.mShape->GetSubType() == EShapeSubType::Sphere)
	{
		Vec3 delta = leaf2.mTransform.GetTranslation() - leaf1.mTransform.GetTranslation();
		float distance = delta.Length();
		if (distance < 1.0e-6f)
			return false;
		float radius1 = static_cast<const SphereShape *>(leaf1.mShape)->GetRadius() * leaf1.mScale.Abs().GetX();
		float radius2 = static_cast<const SphereShape *>(leaf2.mShape)->GetRadius() * leaf2.mScale.Abs().GetX();
		if (abs(distance - radius1 - radius2) > inMaxSeparation)
			return false;
		ioManifold.mWorldSpaceNormal = delta / distance;
		ioManifold.mPenetrationDepth = radius1 + radius2 - distance;
		ioManifold.mRelativeContactPointsOn1.clear();
		ioManifold.mRelativeContactPointsOn2.clear();
		ioManifold.mRelativeContactPointsOn1.push_back(leaf1.mTransform.GetTranslation() + radius1 * ioManifold.mWorldSpaceNormal);
		ioManifold.mRelativeContactPointsOn2.push_back(leaf2.mTransform.GetTranslation() - radius2 * ioManifold.mWorldSpaceNormal);
		return true;
	}
	// Prefer a static supporting surface, independent of BodyID / creation order.
	if (inBody2.IsStatic())
		return RefreshFace(leaf2, leaf1, -inNormalFollowing2, false, inMaxSeparation, ioManifold)
			|| RefreshFace(leaf1, leaf2, inNormalFollowing1, true, inMaxSeparation, ioManifold);
	return RefreshFace(leaf1, leaf2, inNormalFollowing1, true, inMaxSeparation, ioManifold)
		|| RefreshFace(leaf2, leaf1, -inNormalFollowing2, false, inMaxSeparation, ioManifold);
}
}

JPH_NAMESPACE_END
